import subprocess
from collections import defaultdict
from tqdm import tqdm


# -------------------------
# Union-Find (path compression)
# -------------------------
class UnionFind:
    def __init__(self):
        self.parent = {}

    def find(self, x):
        if x not in self.parent:
            self.parent[x] = x
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra


# -------------------------
# run minimap2 as STREAM
# -------------------------
def run_minimap2_stream(fasta, threads=16):
    cmd = [
        "minimap2",
        "-x", "ava-pb",
        "-w", "10", "-e500", "-m40", 
        "-t", str(threads),
        fasta,
        fasta
    ]

    print(" ".join(cmd))

    # IMPORTANT: stdout PIPE, no file
    return subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
        bufsize=1
    )


# -------------------------
# streaming parser + clustering
# -------------------------
def streaming_cluster(fasta, min_cov=0.99, threads=16):
    uf = UnionFind()

    proc = run_minimap2_stream(fasta, threads)
    
    for line in tqdm(proc.stdout, desc="processing ..."):
        print(line)
        cols = line.strip().split("\t")

        if len(cols) < 12:
            continue

        qname = cols[0]
        qlen = int(cols[1])
        tname = cols[5]
        tlen = int(cols[6])

        matches = int(cols[9])
        aln_len = int(cols[10])

        if aln_len == 0:
            continue

        # -------------------------
        # strict B-type rule
        # -------------------------
        identity = matches / aln_len
        if identity < 0.999:
            continue
        # if identity != 1.0:
            # continue

        cov = aln_len / min(qlen, tlen)
        if cov < min_cov:
            continue

        uf.union(qname, tname)

    proc.stdout.close()
    proc.wait()

    return uf


# -------------------------
# build clusters
# -------------------------
def build_clusters(uf, fasta):
    from Bio import SeqIO

    clusters = defaultdict(list)

    for rec in SeqIO.parse(fasta, "fasta"):
        root = uf.find(rec.id)
        clusters[root].append(rec.id)

    return clusters


# -------------------------
# write output
# -------------------------
def write_clusters(clusters, out):
    with open(out, "w") as f:
        for cid, members in enumerate(clusters.values()):
            rep = members[0]
            for m in members:
                f.write(f"{m}\tcluster_{cid}\t{rep}\n")


# -------------------------
# main
# -------------------------
if __name__ == "__main__":
    fasta = "/data1/ccs_data/lingen_16S/deref-uniques/Barcode10_uniques.fasta"
    fasta = "/data1/ccs_data/lingen_16S/deref-uniques/Barcode10.partial.fasta"

    uf = streaming_cluster(
        fasta,
        min_cov=0.99,
        threads=128
    )

    clusters = build_clusters(uf, fasta)

    write_clusters(clusters, "clusters.tsv")
