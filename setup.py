# coding: utf-8

"""Setup file for PyPI"""

from setuptools import setup, find_packages, Extension
from codecs import open
from os import path
import glob

here = path.abspath(path.dirname("__file__"))

with open(path.join(here, "DESCRIPTION.md"), encoding="utf-8") as description:
    long_description = description.read()

setup(

    name="Oracle",
    version="0.1",

    description="A Python3 annotation program to select the best gene model in each locus",
    long_description=long_description,

    url="https://github.com/lucventurini/oracle_trainer.git",

    author="Luca Venturini",
    author_email="luca.venturini@tgac.ac.uk",

    license="GPL3",

    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Gene Annotation",
        "License :: OSI Approved :: GPL3",
        'Programming Language :: Python :: 3.4',
        "Operating System :: Linux"
    ],

    zip_safe=False,

    keywords="rna-seq annotation genomics transcriptomics",

    packages=find_packages(),

    scripts=glob.glob("bin/*.py") + glob.glob("util/*.py"),

    install_requires=[line.rstrip() for line in open("requirements.txt", "rt")],

    extras_require={
        "postgresql": ["psycopg2"],
        "mysql": ["mysqlclient>=1.3.6"],
    },

    test_suite="test",

    data_files=[("Mikado/configuration",
                 glob.glob("Mikado/configuration/*json") + glob.glob("Mikado/configuration/*yaml") )],

)