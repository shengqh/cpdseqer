import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="cpdseqer",
    version="0.0.3",
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

