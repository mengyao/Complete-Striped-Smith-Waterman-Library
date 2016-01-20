#!/usr/bin/env python

from setuptools import setup, Extension

libssw_ext = {"sources": ["src/ssw.c"], "include_dirs": ["src"]}

config = {
    "name": "ssw", 
    "version": "0.1",
    "description": "Complete Striped Smith Waterman Library",
    "author": "Mengyao Zhao et al.",
    "author_email": "zhaomengyao@gmail.com",
    "package_dir": {"ssw": "src"},
    "py_modules": ["ssw.__init__", "ssw.ssw_wrap"],
    "ext_modules": [Extension("libssw", **libssw_ext)],
    "scripts":  ["src/pyssw.py"],
}

if __name__ == "__main__":
    setup(**config)
