from setuptools import setup, find_packages
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="ebg",
    version="0.13.1",
    packages=find_packages(),
    entry_points={
        'console_scripts': ['ebg = EBG.__main__:main']
    },
    author='Julius Wiegert',
    author_email='julius-wiegert@web.de',
    description='Tool for estimating the Felsenstein Bootstrap support of phylogenetic trees',
    long_description=long_description,
    long_description_content_type='text/markdown',
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
