import setuptools
import re
import os

with open("README.md", "r") as fh:
    long_description = fh.read()

PKG = "cpdseqer"
VERSIONFILE = os.path.join(PKG, "__version__.py")
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    print ("unable to find version in %s" % VERSIONFILE)
    raise RuntimeError("if %s exists, it is required to be well-formed" % VERSIONFILE)

setuptools.setup(
    name=PKG,
    version=verstr,
    author="Quanhu Sheng",
    author_email="quanhu.sheng.1@vumc.org",
    description="CPDseq data analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/shengqh/cpdseqer",
    entry_points = {
        'console_scripts': ['cpdseqer=cpdseqer.__main__:main'],
    },
    packages=['cpdseqer'],
    package_dir={'cpdseqer': 'src/cpdseqer'},
    install_requires=['argparse', 'pysam', 'pytabix', 'biopython' ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    zip_safe=False
)

