#!/usr/bin/env python3
"""Compact FastQ encoder/decoder - Version 2.

Binary format (little-endian):
  magic      : 4 bytes  b"CQFQ"
  version    : 1 byte   2
  n_records  : 4 bytes  uint32

Per record:
  name_len   : 2 bytes  uint16
  name       : name_len bytes  ASCII (no trailing \\0)
  seq_len    : 4 bytes  uint32   (number of bases)
  data_packed: seq_len bytes
              1 byte per base+qual pair:
              - bits 7-2 (6 bits, high): qual value (clamped to [0, 63])
              - bits 1-0 (2 bits, low) : base (A=00, C=01, G=10, T=11, N=00)
"""

from __future__ import annotations

import argparse
import gzip
import os
import struct
import sys
from pathlib import Path

MAGIC = b"CQFQ"
VERSION = 2
HEADER_FMT = "<4sBI"          # magic, version, n_records

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
# Packing helpers - v2: 6-bit qual + 2-bit base in single byte
# ---------------------------------------------------------------------------

def pack_bases_and_quals(seq: bytes, qual: bytes) -> tuple[bytes, int]:
    """Pack bases and quals into 1-byte-per-position format.

    Each byte:
      - bits 7-2 (6 bits): qual value (clamped to [0, 63])
      - bits 1-0 (2 bits): base (A=00, C=01, G=10, T=11, N=00)

    Returns:
        (packed_data, clamped_count)
    """
    out = bytearray(len(seq))
    clamped = 0

    for i, (b, q) in enumerate(zip(seq, qual)):
        # Encode base to 2-bit
        base_code = BASE_TO_2BIT.get(bytes([b]), 0)  # N or unknown -> 0

        # Encode qual to 6-bit (Phred+33 -> raw, clamp to 0-63)
        v = q - 33  # strip Phred+33 offset
        if v < 0:
            v = 0
        elif v > 63:
            v = 63
            clamped += 1

        # Pack: qual in high 6 bits, base in low 2 bits
        out[i] = (v << 2) | base_code

    return bytes(out), clamped


def unpack_bases_and_quals(data: bytes, seq_len: int) -> tuple[bytes, bytes]:
    """Unpack 1-byte-per-position format to bases and quals.

    Returns:
        (seq_bytes, qual_bytes)
    """
    seq = bytearray(seq_len)
    qual = bytearray(seq_len)

    for i in range(seq_len):
        byte = data[i]
        # Extract base (low 2 bits)
        base_code = byte & 0x03
        seq[i] = BIT2_TO_BASE[base_code][0]

        # Extract qual (high 6 bits)
        q = (byte >> 2) & 0x3F
        qual[i] = q + 33  # add back Phred+33 offset

    return bytes(seq), bytes(qual)


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
            data_packed, n_clamped = pack_bases_and_quals(seq, qual)
            clamped += n_clamped
            total_quals += len(qual)

            # header per record: name_len (2), seq_len (4), data_packed_len (4)
            fh.write(struct.pack("<HII",
                                 len(name), len(seq), len(data_packed)))
            fh.write(name)
            fh.write(data_packed)

    print(f"Records     : {len(records)}")
    print(f"Output      : {out_path}  ({os.path.getsize(out_path):,} bytes)")
    if clamped:
        print(f"WARNING     : {clamped}/{total_quals} qual values clamped to Phred 63 "
              f"(6-bit max). Data above Phred 63 is lost.")


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
            sys.exit(f"Unsupported version: {ver} (expected {VERSION})")

        for _ in range(n_records):
            name_len, seq_len, data_packed_len = struct.unpack(
                "<HII", fh.read(struct.calcsize("<HII")))
            name = fh.read(name_len)
            data_packed = fh.read(data_packed_len)

            seq, qual = unpack_bases_and_quals(data_packed, seq_len)

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
        description="Encode FASTQ into a compact 6-bit qual + 2-bit base binary format, "
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