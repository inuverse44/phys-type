#!/usr/bin/env python3
"""Fetch the BICEP/Keck r-n_s contour data used by the n_s-r figures.

BICEP/Keck XIII (arXiv:2110.00483) released ``rns_code.tgz``, the gridded
posterior densities and contour levels behind its Figure 5: Planck 2018
TT,TE,EE+lowE+lensing, and the same combined with BK18 and BAO, in the plane
(n_s, r_0.002). The data are not redistributed here; this script downloads
the archive from the BICEP/Keck release page, checks its SHA-256, and
extracts the six text files into ``data/bk18_rns/`` (ignored by git).

    python3 scripts/fetch_bk18.py

Only the Python standard library is used.
"""

import hashlib
import io
import sys
import tarfile
import urllib.request
from pathlib import Path

URL = "http://bicepkeck.org/BK18_datarelease/rns_code.tgz"
SHA256 = "0612f8b5c630bbb5ce0e5586f244334152f236b4605ac0934131a5cf9f142641"
FILES = [
    f"{prefix}density_{part}.txt"
    for prefix in ("P18", "BK18")
    for part in ("xy", "P", "levels")
]
TARGET = Path(__file__).resolve().parents[1] / "data" / "bk18_rns"


def main() -> int:
    if all((TARGET / name).is_file() for name in FILES):
        print(f"already present in {TARGET}")
        return 0
    with urllib.request.urlopen(URL, timeout=120) as response:
        archive = response.read()
    digest = hashlib.sha256(archive).hexdigest()
    if digest != SHA256:
        print(f"checksum mismatch: {digest}", file=sys.stderr)
        return 1
    TARGET.mkdir(parents=True, exist_ok=True)
    with tarfile.open(fileobj=io.BytesIO(archive), mode="r:gz") as tar:
        for name in FILES:
            member = tar.getmember(f"rns_code/{name}")
            (TARGET / name).write_bytes(tar.extractfile(member).read())
    print(f"extracted {len(FILES)} files into {TARGET}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
