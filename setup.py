from setuptools import setup, find_packages

setup(
    name="ebg",
    version="0.11",
    packages=find_packages(),
    entry_points={
        'console_scripts': ['ebg = EBG.__main__:main']
    },
    author='Julius Wiegert',
    author_email='julius-wiegert@web.de',
    description='Tool for estimating the Felsenstein Bootstrap support of phylogenetic trees',
    long_description='EBG is a fast predictor for the standard non-parametric Felsenstein Bootstrap support (SBS) of phylogenetic trees in python. It uses RAxML-NG phylogenies as input and predicts the SBS values of all non-trivial branches',
    include_package_data=True,
    install_requires=[
        "pandas",
        "numpy",
        "ete3",
        "biopython",
        "networkx",
        "scipy",
        "lightgbm"

    ],
    package_data={
        'ebg': ['EBG/Models/*.pkl'],
    },
)
