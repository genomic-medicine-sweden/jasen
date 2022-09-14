import os

from setuptools import find_packages, setup

with open("README.md") as f:
    long_description = f.read()

setup(
    name="prp",
    version="1.0.0",
    description="Pipeline result processor for JASEN pipeline",
    long_description_markdown_filename=long_description,
    long_description_content_type="text/markdown",
    author="Markus Johansson",
    license="GPLv3",
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3.10",
    ],
    install_requires=["Click", "pydantic", "pandas"],
    entry_points={"console_scripts": ["prp=prp.cli:cli"]},
    packages=find_packages(exclude=("tests")),
)
