#!/usr/bin/env python3
"""Release one completed SPARTA dump from Linux page cache, without root."""
import os
from pathlib import Path
import stat
import sys


def evict(path: Path) -> int:
    if not hasattr(os, "posix_fadvise"):
        raise OSError("POSIX_FADV_DONTNEED is unavailable on this platform")
    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    try:
        info = os.fstat(fd)
        if not stat.S_ISREG(info.st_mode):
            raise OSError(f"{path} is not a regular file")
        # DONTNEED only discards clean pages; flush this completed dump first.
        os.fsync(fd)
        os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
        return info.st_size
    finally:
        os.close(fd)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: evict_dump_cache.py DUMP_FILE")
    try:
        size = evict(Path(sys.argv[1]))
    except OSError as exc:
        raise SystemExit(f"cache release failed: {exc}") from exc
    print(f"requested page-cache release for {sys.argv[1]} ({size} bytes)")
