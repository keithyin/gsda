import os
import glob
import ctypes
import pysam
import multiprocessing as mp
import time
from tqdm import tqdm

libc = ctypes.CDLL("libc.so.6", use_errno=True)
POSIX_FADV_DONTNEED = 4


def posix_fadvise(fd, offset, length, advice):
    ret = libc.posix_fadvise(
        ctypes.c_int(fd),
        ctypes.c_long(offset),
        ctypes.c_long(length),
        ctypes.c_int(advice),
    )
    if ret != 0:
        errno = ctypes.get_errno()
        raise OSError(errno, os.strerror(errno))


def process_bam(path):
    fd = os.open(path, os.O_RDONLY)

    try:
        bam = pysam.AlignmentFile(
            path, "rb", check_sq=False, threads=mp.cpu_count())
        cnt = 0
        for i, read in tqdm(enumerate(bam), desc=f"reading {path}"):
            cnt = i
            pass
        print(f"{path}: cnt={cnt}")
        bam.close()

    finally:
        # posix_fadvise(fd, 0, 0, POSIX_FADV_DONTNEED)
        os.close(fd)

    time.sleep(30)


if __name__ == "__main__":
    pattern = "/data1/ccs_data/chicken-data/*_adapter.bam"

    for bam_path in glob.glob(pattern):
        print(f"Processing: {bam_path}")
        process_bam(bam_path)
