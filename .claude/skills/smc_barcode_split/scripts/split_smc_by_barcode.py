#!/usr/bin/env python3
"""Split an SMC (consensus) all-reads BAM -- and, optionally, the adapter
(subread) BAM it was made from -- into per-barcode BAMs.

SMC pass
--------
The grouping key is the subread id. The barcode-split FASTQs (one per
barcode, in ``--barcode-dir``) hold header lines that identify which
subreads belong to which barcode. The SMC BAM's queryname and the FASTQ
header are not guaranteed to be byte-identical, so BOTH sides are
normalized with :func:`map_qname` before matching: the mapping keeps only
the first two ``/``-separated fields (``"/".join(q.split("/")[:2])``).
This is applied unconditionally to the FASTQ headers and the BAM
querynames so that a header and a queryname that differ only after the
second ``/`` still join correctly.

Adapter pass (``--adapter-bam``)
--------------------------------
``<Run>_adapter.bam`` is the subread BAM that smicing consumed to produce
``<Run>.smc_all_reads.bam``; it holds every subread of every channel, so
it is 10-50x larger than the SMC BAM. The two files join on the channel
id: the ``ch`` tag (present on both sides; equal to the second ``/``-field
of either queryname -- ``run/123456`` for SMC, ``read_123456/123456/
subread/7`` for adapter records). So the adapter split is driven by the SMC
result, never by the FASTQs directly: ``ch -> barcode`` is collected while
splitting the SMC BAM, or, with ``--smc-dir``, read back out of per-barcode
SMC BAMs that a previous run already wrote. Every adapter record whose
channel belongs to a barcode is written to that barcode's adapter BAM --
all subreads of the channel, in source order.

Output names are the *source* BAM's stem with the barcode name inserted as
an infix before the extension, so the two passes never collide:
``Run.smc_all_reads.bam`` + ``Barcode01`` -> ``Run.smc_all_reads.Barcode01.bam``,
``Run_adapter.bam`` + ``Barcode01`` -> ``Run_adapter.Barcode01.bam``.
"""
import argparse
import collections
import glob
import os
import subprocess
import sys
import time

try:
    import pysam
except ImportError:
    sys.exit("pysam is required: pip install pysam")

DEFAULT_GLOB = "Barcode*.fastq"
DEFAULT_SMC_GLOB = "*.bam"
CHANNEL_TAG = "ch"
UNASSIGNED = "unassigned"


def map_qname(qname):
    """Normalize a queryname / FASTQ header to the matching key.

    Keeps only the first two ``/``-separated fields, e.g.
    ``20260831_250302Y0001_Run0001/1076757`` -> itself (already 2 fields)
    or ``A/B/1076757`` -> ``A/B``. Applied to both the FASTQ header and the
    BAM queryname so they always join on the same key.
    """
    return "/".join(qname.split("/")[:2])


def norm_channel(value):
    """Normalize a channel id to a string key.

    Both sides store the channel as an integer (``ch`` tag) or as a decimal
    field inside the queryname; canonicalizing to ``str(int(v))`` makes the
    tag value and the qname-derived value interchangeable.
    """
    s = str(value).strip()
    try:
        return str(int(s))
    except ValueError:
        return s


def channel_of(rec, tag=CHANNEL_TAG, qname_field=1):
    """Channel id of an alignment record: the ``ch`` tag, else qname field.

    Falls back to the ``qname_field``-th ``/``-separated field of the
    queryname when the tag is absent (0 = whole name, 1 = second field,
    which is the channel in both ``run/123`` and ``read_123/123/subread/0``).
    Returns None if neither is usable.
    """
    try:
        return norm_channel(rec.get_tag(tag))
    except KeyError:
        pass
    fields = rec.query_name.split("/")
    if len(fields) > qname_field:
        return norm_channel(fields[qname_field])
    return None


def header_run_name(header):
    """Run name recorded in a BAM header's read group, or None.

    This pipeline stores it under the non-standard ``rn`` key (``SM``/``LB``
    are checked as fallbacks). An adapter BAM and the SMC BAM derived from it
    carry the same value, so it is the only run identity available on the
    adapter side -- its querynames are ``read_<ch>/<ch>/subread/<i>`` and
    channel ids repeat between runs, which is exactly the mispairing this
    catches.
    """
    for rg in (header.to_dict().get("RG") or []):
        for key in ("rn", "SM", "LB"):
            val = rg.get(key)
            if val:
                return str(val)
    return None


def open_bam(path):
    """Open an alignment file for reading, tolerating a header with no @SQ.

    Some adapter BAMs in this pipeline carry no reference lines at all (every
    record is unmapped), which pysam refuses to open without ``check_sq``.
    """
    return pysam.AlignmentFile(path, "rb", check_sq=False, threads=os.cpu_count() // 2)


class BamWriters:
    """Lazily opened ``<stem>.<label>.bam`` writers, one per barcode.

    Writers are created on first use, so barcodes with no reads produce no
    files. Each output carries the source BAM's header verbatim (no PG
    record is injected -- provenance lives in this tool's stdout summary).
    """

    def __init__(self, outdir, stem, header, threads=1):
        self.outdir = outdir
        self.stem = stem
        self.header = header
        self.threads = threads
        self.writers = {}

    def path_for(self, label):
        return os.path.join(self.outdir, f"{self.stem}.{label}.bam")

    def write(self, label, rec):
        w = self.writers.get(label)
        if w is None:
            w = pysam.AlignmentFile(self.path_for(label), "wb", threads=self.threads,
                                    header=self.header, check_sq=False)
            self.writers[label] = w
        w.write(rec)

    def close(self):
        for w in self.writers.values():
            w.close()
        self.writers = {}


def add_channel(ch2label, ch, label, stats):
    """Record ch -> label, keeping the first barcode on a channel collision."""
    prev = ch2label.get(ch)
    if prev is None:
        ch2label[ch] = label
    elif prev != label:
        stats["collisions"] += 1


def load_id_map(barcode_dir, glob_pat):
    """Return ({qname: barcode label}, per-file-id-counts).

    Reads each FASTQ matching ``glob_pat``, taking every 4th line (the
    header, 0-indexed %4==0, starting with ``@``), strips the ``@`` and
    normalizes with :func:`map_qname`. The barcode label is the FASTQ base
    name without ``.fastq``; it becomes the infix of the output BAM name.
    """
    id2label = {}
    per_file = {}
    for fq in sorted(glob.glob(os.path.join(barcode_dir, glob_pat))):
        base = os.path.splitext(os.path.basename(fq))[0]
        n = 0
        with open(fq) as fh:
            for i, line in enumerate(fh):
                if i % 4 == 0 and line.startswith("@"):
                    id2label[map_qname(line[1:].strip())] = base
                    n += 1
        per_file[base] = n
    return id2label, per_file


def channel_map_from_smc_dir(smc_dir, smc_glob, tag, qname_field):
    """Build (ch -> label, per-label counts, stats, run names) from split SMC BAMs.

    Lets the adapter pass run on its own against the ``barcoded_smc/`` output
    of an earlier run -- the channel set is then exactly what the SMC split
    produced, with no re-read of the big all-reads SMC BAM. Each file name is
    split at the *last* dot: the label is the trailing component
    (``Run.smc_all_reads.Barcode07.bam`` -> ``Barcode07``; a plain
    ``Barcode07.bam`` -> ``Barcode07`` too). ``.unassigned.bam`` files are
    skipped: an unassigned channel has no barcode to propagate.
    """
    ch2label = {}
    counts = collections.Counter()
    stats = collections.Counter()
    runs = set()
    # Name first, read second: a channel id is only unique *within* a run, so a
    # directory holding the splits of several source BAMs would silently merge
    # their channels. Detect that from the file names before streaming anything.
    wanted = []
    for path in sorted(glob.glob(os.path.join(smc_dir, smc_glob))):
        name = os.path.splitext(os.path.basename(path))[0]
        label = name.split(".")[-1]
        if label == UNASSIGNED:
            print(f"  skipping {os.path.basename(path)} (not a barcode)")
            continue
        wanted.append((path, ".".join(name.split(".")[:-1]), label))
    if not wanted:
        sys.exit(f"no BAMs matching {smc_glob!r} in --smc-dir {smc_dir}")
    label2stems = collections.defaultdict(set)
    for _, stem, label in wanted:
        label2stems[label].add(stem)
    clash = sorted(lbl for lbl, st in label2stems.items() if len(st) > 1)
    if clash:
        stems = {st for _, st, _ in wanted}
        sys.exit(f"--smc-dir {smc_dir} holds splits of {len(stems)} "
                 f"different source BAMs that repeat barcode names "
                 f"({clash[:4]}{'...' if len(clash) > 4 else ''}); channel ids only "
                 "collide across runs, so their channels would be merged. Point "
                 "--smc-dir at one run's output, or narrow --smc-glob to that source "
                 "BAM's stem (e.g. 'Run0001.smc_all_reads.*.bam').")
    for path, _, label in wanted:
        n = 0
        with open_bam(path) as fh:
            run = header_run_name(fh.header)
            if run:
                runs.add(run)
            for rec in fh.fetch(until_eof=True):
                n += 1
                ch = channel_of(rec, tag, qname_field)
                if ch is None:
                    stats["no_channel"] += 1
                    continue
                add_channel(ch2label, ch, label, stats)
        counts[label] += n
    return ch2label, counts, stats, runs


def parse_args():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--barcode-dir",
                    help="directory with per-barcode FASTQs (e.g. barcode_assign); "
                         "required with --smc-bam")
    ap.add_argument("--smc-bam",
                    help="the SMC all-reads BAM to split")
    ap.add_argument("--outdir", required=True,
                    help="output directory for the per-barcode SMC BAMs")
    ap.add_argument("--adapter-bam",
                    help="also split this subread/adapter BAM (e.g. Run_adapter.bam) "
                         "using the channels found in the SMC reads")
    ap.add_argument("--adapter-outdir",
                    help="output directory for the per-barcode adapter BAMs "
                         "(default: --outdir)")
    ap.add_argument("--smc-dir",
                    help="directory of per-barcode SMC BAMs from a previous run, to "
                         "derive ch -> barcode from instead of re-reading --smc-bam "
                         "(adapter-only mode)")
    ap.add_argument("--smc-glob", default=DEFAULT_SMC_GLOB,
                    help=f"glob for --smc-dir BAMs (default {DEFAULT_SMC_GLOB})")
    ap.add_argument("--glob", default=DEFAULT_GLOB,
                    help=f"glob for barcode FASTQs in --barcode-dir (default {DEFAULT_GLOB})")
    ap.add_argument("--channel-tag", default=CHANNEL_TAG,
                    help=f"BAM tag holding the channel id (default {CHANNEL_TAG}); "
                         "falls back to a queryname field when absent")
    ap.add_argument("--channel-qname-field", type=int, default=1,
                    help="0-based /-separated queryname field holding the channel "
                         "when the tag is missing (default 1)")
    ap.add_argument("--unassigned", action="store_true",
                    help=f"dump SMC reads matching no barcode into "
                         f"outdir/<smc-bam-stem>.{UNASSIGNED}.bam")
    ap.add_argument("--unassigned-adapters", action="store_true",
                    help=f"dump adapter records whose channel matches no barcode into "
                         f"adapter-outdir/<adapter-bam-stem>.{UNASSIGNED}.bam (can be huge)")
    ap.add_argument("--allow-run-mismatch", action="store_true",
                    help="proceed with the adapter pass even when the adapter BAM's "
                         "read-group run name differs from the SMC reads' (the join "
                         "is then almost certainly wrong; default: abort)")
    ap.add_argument("--index", action="store_true",
                    help="samtools index each output BAM after writing")
    ap.add_argument("--threads", type=int, default=4,
                    help="bgzf writer threads per output BAM (default 4); writing "
                         "is the bottleneck here, ~10x faster than threads=1. Peak "
                         "thread count is threads x (barcodes touched), so lower it "
                         "if the box is small")
    args = ap.parse_args()

    if not args.smc_bam and not args.adapter_bam:
        ap.error("nothing to do: pass --smc-bam and/or --adapter-bam")
    if args.smc_bam and not args.barcode_dir:
        ap.error("--smc-bam requires --barcode-dir")
    if args.adapter_bam and not (args.smc_bam or args.smc_dir):
        ap.error("--adapter-bam needs a channel source: pass --smc-bam (split the "
                 "SMC BAM in the same run) or --smc-dir (pre-split per-barcode SMC BAMs)")
    if args.smc_dir and args.smc_bam:
        ap.error("--smc-dir and --smc-bam are mutually exclusive; drop one "
                 "(both write the same channel map)")
    return args


def main():
    args = parse_args()
    t_start = time.time()
    os.makedirs(args.outdir, exist_ok=True)
    adapter_outdir = args.adapter_outdir or args.outdir
    if args.adapter_bam and args.adapter_outdir:
        os.makedirs(args.adapter_outdir, exist_ok=True)

    smc_stem = (os.path.splitext(os.path.basename(args.smc_bam))[0]
                if args.smc_bam else None)
    ch2label = {}
    smc_counts = collections.Counter()
    smc_stats = collections.Counter()
    smc_runs = set()

    # ---------------- SMC pass ----------------
    if args.smc_bam:
        id2label, per_file = load_id_map(args.barcode_dir, args.glob)
        print("== SMC split ==")
        print(f"qname mapping: first 2 /-fields (forced)")
        print(f"smc bam stem: {smc_stem}")
        print(f"barcode fastqs: {len(per_file)}")
        for base, n in sorted(per_file.items()):
            print(f"  {base}: {n} ids")
        print(f"total mapped ids: {len(id2label)}")

        total = 0
        unmatched = 0
        with open_bam(args.smc_bam) as src:
            writers = BamWriters(args.outdir, smc_stem, src.header, args.threads)
            run = header_run_name(src.header)
            if run:
                smc_runs.add(run)
            for rec in src.fetch(until_eof=True):
                total += 1
                label = id2label.get(map_qname(rec.query_name))
                if label is None:
                    if not args.unassigned:
                        unmatched += 1
                        continue
                    label = UNASSIGNED
                writers.write(label, rec)
                smc_counts[label] += 1
                ch = channel_of(rec, args.channel_tag, args.channel_qname_field)
                if ch is None:
                    smc_stats["no_channel"] += 1
                else:
                    add_channel(ch2label, ch, label, smc_stats)
        writers.close()
        smc_written = total - unmatched
        print(f"\ntotal SMC records:   {total}")
        print(f"unmatched:           {unmatched}")
        print(f"written records:     {smc_written}")
        print(f"channels mapped:     {len(ch2label)}")
        if smc_stats["collisions"]:
            print(f"WARN channel->barcode collisions (kept first): "
                  f"{smc_stats['collisions']}")
        if smc_stats["no_channel"]:
            print(f"WARN SMC records with no usable channel: {smc_stats['no_channel']}")
        for label in sorted(smc_counts):
            print(f"  {os.path.basename(writers.path_for(label))}: {smc_counts[label]}")
        print(f"outdir: {args.outdir}")
        if total and not smc_written:
            sys.exit("no SMC record matched any barcode FASTQ header "
                     f"({total} read, all unmatched) -- the qname convention is not "
                     "'same first two /-fields' for this pair; see the qname notes "
                     "in SKILL.md. No output BAMs were written to "
                     f"{args.outdir}.")

    # ---------------- channel map from pre-split SMC ----------------
    elif args.smc_dir:
        print("== channel map from pre-split SMC BAMs ==")
        ch2label, smc_counts, smc_stats, smc_runs = channel_map_from_smc_dir(
            args.smc_dir, args.smc_glob, args.channel_tag,
            args.channel_qname_field)
        print(f"source dir: {args.smc_dir}")
        print(f"barcodes:   {len(smc_counts)}")
        print(f"channels mapped: {len(ch2label)}")
        if len(smc_runs) > 1:
            print(f"WARN --smc-dir BAMs name {len(smc_runs)} different runs: "
                  f"{sorted(smc_runs)} -- channel ids collide across runs, so this "
                  "map is probably not what you want")
        if smc_stats["collisions"]:
            print(f"WARN channel->barcode collisions (kept first): "
                  f"{smc_stats['collisions']}")
        if smc_stats["no_channel"]:
            print(f"WARN SMC records with no usable channel: {smc_stats['no_channel']}")
        for label in sorted(smc_counts):
            print(f"  {label}: {smc_counts[label]} smc reads")

    # ---------------- Adapter pass ----------------
    adapter_counts = collections.Counter()
    adapter_ch = collections.Counter()
    if args.adapter_bam:
        if not ch2label:
            sys.exit("empty channel map -- refusing to write an adapter split of "
                     "nothing. Check that the SMC querynames/`ch` tags line up with "
                     "the barcode FASTQ headers (or that --smc-dir BAMs carry ch).")
        print(f"\n== adapter split ==")
        adapter_stem = os.path.splitext(os.path.basename(args.adapter_bam))[0]
        print(f"adapter bam stem: {adapter_stem}")
        print(f"adapter outdir:   {adapter_outdir}")
        print(f"streaming (source order is not channel order, so one full pass): "
              f"{args.adapter_bam}")

        total = 0
        written = 0
        unmatched = 0
        adapter_stats = collections.Counter()
        seen = set()
        t0 = time.time()
        with open_bam(args.adapter_bam) as src:
            adapter_run = header_run_name(src.header)
            if adapter_run and smc_runs and adapter_run not in smc_runs:
                msg = (f"run mismatch: {args.adapter_bam} declares RG run "
                       f"{adapter_run!r} but the channel map came from "
                       f"{sorted(smc_runs)}. Channel ids are only unique within a "
                       "run, so this pairing would mix reads across runs.")
                if not args.allow_run_mismatch:
                    sys.exit(msg + " Aborted before writing anything (pass "
                                    "--allow-run-mismatch to do it anyway).")
                print(f"WARN {msg}")
            writers = BamWriters(adapter_outdir, adapter_stem, src.header,
                                args.threads)
            for rec in src.fetch(until_eof=True):
                total += 1
                ch = channel_of(rec, args.channel_tag, args.channel_qname_field)
                if ch is None:
                    adapter_stats["no_channel"] += 1
                label = ch2label.get(ch) if ch is not None else None
                if label is None:
                    if not args.unassigned_adapters:
                        unmatched += 1
                        continue
                    label = UNASSIGNED
                writers.write(label, rec)
                adapter_counts[label] += 1
                written += 1
                if ch is not None and ch not in seen:
                    seen.add(ch)
                    adapter_ch[label] += 1
                if total % 1000000 == 0:
                    rate = total / (time.time() - t0)
                    print(f"  ...{total} records read "
                          f"({time.time() - t0:.0f}s, {rate:.0f} rec/s)", flush=True)
        writers.close()
        covered = sum(adapter_ch.values())
        print(f"\ntotal adapter records:   {total}")
        print(f"written records:         {written}")
        print(f"unmatched records:       {unmatched}")
        print(f"channels in map:         {len(ch2label)}")
        print(f"channels covered:        {covered} "
              f"({100 * covered / len(ch2label):.1f}% of the SMC channels)")
        if adapter_stats["no_channel"]:
            print(f"WARN adapter records with no usable channel: "
                  f"{adapter_stats['no_channel']}")
        if not written:
            sys.exit(f"no adapter record matched any SMC channel -- "
                     f"{args.adapter_bam} is not the BAM these SMC reads came "
                     "from (or its channel tag/field differs). No output files "
                     "were written; see the qname/channel notes in SKILL.md.")

    # ---------------- reconciliation ----------------
    if adapter_counts:
        print("\n== per-barcode reconciliation ==")
        print(f"{'barcode':<16}{'smc_reads':>11}{'adapter_ch':>12}"
              f"{'adapter_recs':>14}{'recs/ch':>9}")
        for label in sorted(set(smc_counts) | set(adapter_counts)):
            s = smc_counts.get(label, 0)
            a = adapter_counts.get(label, 0)
            c = adapter_ch.get(label, 0)
            ratio = f"{a / c:.1f}" if c else "-"
            print(f"{label:<16}{s:>11}{c:>12}{a:>14}{ratio:>9}")
        print("(adapter_ch should track smc_reads: 1 SMC consensus read per "
              "channel. A gap = that channel has no record in this adapter BAM; "
              "recs/ch is the mean subread depth, ~12-17 for a real run.)")

    if args.index:
        paths = []
        for d, stem in ((args.outdir, smc_stem),
                        (adapter_outdir, args.adapter_bam
                         and os.path.splitext(os.path.basename(args.adapter_bam))[0])):
            if not stem:
                continue
            paths += sorted(glob.glob(os.path.join(d, f"{stem}.*.bam")))
        for path in paths:
            subprocess.run(["samtools", "index", path], check=True)
        print(f"\nindexed {len(paths)} BAMs")

    print(f"\nelapsed: {time.time() - t_start:.0f}s")


if __name__ == "__main__":
    sys.exit(main())
