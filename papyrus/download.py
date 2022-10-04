"""
download

Created by: Martin Sicho
On: 04.10.22, 14:48
"""

from papyrus_scripts.download import download_papyrus
from papyrus import OUTDIR

download_papyrus(
    outdir=OUTDIR,
    version='latest',
    descriptors=[],
    stereo=True,
    disk_margin=0.05,
)
