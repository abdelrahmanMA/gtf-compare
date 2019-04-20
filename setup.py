#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from setuptools import setup, find_packages

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "gtfCompare",
    version = "0.0.1",
    author = "Abdelrahman Muhammed",
    author_email = "bme.abdelrahmanmuhammed@gmail.com",
    description = ("A GTF/GFF Compare Tool"),
    license = "BSD",
    keywords = "example documentation tutorial",
    url = "https://github.com/abdelrahmanMA/gtf-compare",
    packages= find_packages(),
    long_description=read('README'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: BSD License",
    ],
)