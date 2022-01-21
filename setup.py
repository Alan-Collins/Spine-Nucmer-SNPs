from setuptools import setup

setup(
    name="spine-nucmer-snps",
    version="0.1.0",
    author="Alan Collins",
    autor_email="alan.collins@bath.edu",
    url="https://github.com/Alan-Collins/Spine-Nucmer-SNPs",
    license="GPL-3.0",
    install_requires=[
        "pip",
        "matplotlib",
        "setuptools",
    ],
    scripts=[
        "get_snps_support_MP.py",
        "fasta2diffmat.py",
        "snps2fasta.py"

    ],
    packages = ["snpclasses"],
    package_dir={"": "./"}
)

