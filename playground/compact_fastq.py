#!/usr/bin/env python3
"""Compact FastQ encoder/decoder.

Binary format (little-endian):
  magic      : 4 bytes  b"CQFQ"
  version    : 1 byte   1
  n_records  : 4 bytes  uint32

Per record:
  name_len   : 2 bytes  uint16
  name       : name_len bytes  ASCII (no trailing \\0)
  seq_len    : 4 bytes  uint32   (number of bases)
  seq_packed : ceil(seq_len / 4) bytes
               2 bits per base, MSB first.
               A=00  C=01  G=10  T=11  N=00 (fallback)
               Last byte is right-padded with 00 if seq_len not a multiple of 4.
  qual_packed: ceil(seq_len / 2) bytes
               4 bits per qual, MSB first (high nibble = first qual).
               qual value is clamped to [0, 15].
               Last byte low nibble is 0 if seq_len is odd.
"""

from __future__ import annotations

import argparse
import gzip
import os
import struct
import sys
from pathlib import Path

MAGIC = b"CQFQ"
VERSION = 1
HEADER_FMT = "<4sBI"          # magic, version, n_records
RECORD_HEAD_FMT = "<HBI"       # name_len (uint16), seq_len (uint32) — note: we pack name_len + seq_len together

BASE_TO_2BIT = {b"A": 0, b"a": 0, b"C": 1, b"c": 1,
                b"G": 2, b"g": 2, b"T": 3, b"t": 3}
BIT2_TO_BASE = [b"A", b"C", b"G", b"T"]

# ---------------------------------------------------------------------------
# FASTQ parser (minimal, no external deps)
# ---------------------------------------------------------------------------

def iter_fastq(path: str | Path):
    """Yield (name: bytes, seq: bytes, qual: bytes) from a plain or gzipped FASTQ."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            if not header.startswith(b"@"):
                raise ValueError(f"Expected '@' at start of header, got: {header!r}")
            seq = fh.readline().rstrip(b"\r\n")
            fh.readline()                # '+' line
            qual = fh.readline().rstrip(b"\r\n")
            if len(seq) != len(qual):
                raise ValueError(f"seq/qual length mismatch in record: {header!r}")
            name = header[1:].rstrip(b"\r\n")
            yield name, seq, qual


# ---------------------------------------------------------------------------
# Packing helpers
# ---------------------------------------------------------------------------

def pack_bases(seq: bytes) -> bytes:
    """Pack bases into 2-bit-per-base bytes, MSB first."""
    out = bytearray()
    acc, shift = 0, 6
    for b in seq:
        code = BASE_TO_2BIT.get(bytes([b]), 0)
        acc |= code << shift
        if shift == 0:
            out.append(acc)
            acc, shift = 0, 6
        else:
            shift -= 2
    if shift != 6:                  # leftover bases
        out.append(acc)
    return bytes(out)


def pack_quals(qual: bytes) -> bytes:
    """Pack Phred+33 qual values into 4-bit-per-value bytes, high-nibble first.

    Values are clamped to [0, 15] — anything above Phred 15 maps to 15.
    """
    data, _ = _pack_quals_counted(qual)
    return data


def _pack_quals_counted(qual: bytes) -> tuple[bytes, int]:
    """Like pack_quals but also returns the number of values clamped to 15."""
    out = bytearray()
    hi = True
    clamped = 0
    for q in qual:
        v = q - 33                   # strip Phred+33 offset
        v= v // 3
        if v < 0:
            v = 0
        elif v > 15:
            v = 15
            clamped += 1
        if hi:
            out.append(v << 4)
            hi = False
        else:
            out[-1] |= v
            hi = True
    return bytes(out), clamped


def unpack_bases(data: bytes, seq_len: int) -> bytes:
    out = bytearray()
    for i in range(seq_len):
        byte_idx = i // 4
        shift = 6 - 2 * (i % 4)
        code = (data[byte_idx] >> shift) & 0x03
        out.extend(BIT2_TO_BASE[code])
    return bytes(out)


def unpack_quals(data: bytes, seq_len: int) -> bytes:
    out = bytearray()
    for i in range(seq_len):
        byte_idx = i // 2
        if i % 2 == 0:
            v = (data[byte_idx] >> 4) & 0x0F
        else:
            v = data[byte_idx] & 0x0F
        out.append(v + 33)
    return bytes(out)


# ---------------------------------------------------------------------------
# Encode: FASTQ -> compact binary
# ---------------------------------------------------------------------------

def encode(fastq_path: str | Path, out_path: str | Path, *, compress: bool = True) -> None:
    out_path = Path(out_path)
    records = list(iter_fastq(fastq_path))

    clamped = 0
    total_quals = 0

    opener = gzip.open if compress else open
    with opener(out_path, "wb") as fh:
        fh.write(struct.pack(HEADER_FMT, MAGIC, VERSION, len(records)))
        for name, seq, qual in records:
            seq_packed = pack_bases(seq)
            qual_packed, n_clamped = _pack_quals_counted(qual)
            clamped += n_clamped
            total_quals += len(qual)
            # header per record: name_len (2), seq_len (4), seq_packed_len (4), qual_packed_len (4)
            fh.write(struct.pack("<HIII",
                                 len(name), len(seq),
                                 len(seq_packed), len(qual_packed)))
            fh.write(name)
            fh.write(seq_packed)
            fh.write(qual_packed)

    print(f"Records     : {len(records)}")
    print(f"Output      : {out_path}  ({os.path.getsize(out_path):,} bytes)")
    if clamped:
        print(f"WARNING     : {clamped}/{total_quals} qual values clamped to Phred 15 "
              f"(4-bit max). Data above Phred 15 is lost.")


# ---------------------------------------------------------------------------
# Decode: compact binary -> FASTQ
# ---------------------------------------------------------------------------

def _is_gzip(path: Path) -> bool:
    """Check by magic bytes so we don't rely on file extension."""
    try:
        with open(path, "rb") as f:
            return f.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def decode(in_path: str | Path, out_path: str | Path) -> None:
    in_path = Path(in_path)
    opener = gzip.open if _is_gzip(in_path) else open
    out_opener = gzip.open if str(out_path).endswith(".gz") else open
    with opener(in_path, "rb") as fh, out_opener(out_path, "wb") as out:
        magic, ver, n_records = struct.unpack(HEADER_FMT, fh.read(struct.calcsize(HEADER_FMT)))
        if magic != MAGIC:
            sys.exit(f"Bad magic: {magic!r}")
        if ver != VERSION:
            sys.exit(f"Unsupported version: {ver}")

        for _ in range(n_records):
            name_len, seq_len, seq_packed_len, qual_packed_len = struct.unpack(
                "<HIII", fh.read(struct.calcsize("<HIII")))
            name = fh.read(name_len)
            seq_packed = fh.read(seq_packed_len)
            qual_packed = fh.read(qual_packed_len)

            seq = unpack_bases(seq_packed, seq_len)
            qual = unpack_quals(qual_packed, seq_len)

            out.write(b"@")
            out.write(name)
            out.write(b"\n")
            out.write(seq)
            out.write(b"\n+\n")
            out.write(qual)
            out.write(b"\n")

    print(f"Decoded {n_records} records -> {out_path}")


# ---------------------------------------------------------------------------
# Stats
# ---------------------------------------------------------------------------

def stats(in_path: str | Path, original_path: str | Path | None = None) -> None:
    size = os.path.getsize(in_path)
    print(f"File           : {in_path}")
    print(f"Size           : {size:,} bytes")
    if original_path:
        orig = os.path.getsize(original_path)
        print(f"Original FASTQ : {orig:,} bytes")
        print(f"Ratio          : {size / orig:.2%}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    p = argparse.ArgumentParser(
        description="Encode FASTQ into a compact 2-bit base / 4-bit qual binary format, "
                    "optionally gzip-compressed.")
    sub = p.add_subparsers(dest="cmd", required=True)

    enc = sub.add_parser("encode", help="FASTQ -> compact binary")
    enc.add_argument("input", help="Input FASTQ (plain or .gz)")
    enc.add_argument("-o", "--output", default=None,
                     help="Output path (default: <input>.cqfq.gz)")
    enc.add_argument("--no-gzip", action="store_true",
                     help="Do not gzip-compress the output")

    dec = sub.add_parser("decode", help="compact binary -> FASTQ")
    dec.add_argument("input", help="Input .cqfq or .cqfq.gz")
    dec.add_argument("-o", "--output", default=None,
                     help="Output FASTQ path (default: stdout)")

    st = sub.add_parser("stats", help="Print file size info")
    st.add_argument("input")
    st.add_argument("--original", default=None, help="Original FASTQ for ratio")

    args = p.parse_args()

    if args.cmd == "encode":
        if args.output:
            out = args.output
            if args.no_gzip:
                compress = False
            else:
                # infer from extension: .gz -> compress, .cqfq -> no, anything else -> compress
                compress = not out.endswith(".cqfq")
        else:
            base = Path(args.input).name
            if base.endswith(".gz"):           # strip input .gz
                base = base[:-3]
            stem = base.rsplit(".", 1)[0] if "." in base else base
            suffix = ".cqfq" if args.no_gzip else ".cqfq.gz"
            out = stem + suffix
            compress = not args.no_gzip
        encode(args.input, out, compress=compress)

    elif args.cmd == "decode":
        out = args.output or "/dev/stdout"
        decode(args.input, out)

    elif args.cmd == "stats":
        stats(args.input, args.original)


if __name__ == "__main__":
    main()
