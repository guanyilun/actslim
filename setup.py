"""Build script for the actslim C extension.

All package metadata lives in pyproject.toml; this file exists only to declare
the compiled extension (the vendored slim decoder + the CPython binding).
"""
import glob
import os

from setuptools import setup, Extension

here = os.path.dirname(os.path.abspath(__file__))
vendor = sorted(glob.glob(os.path.join("actslim", "_vendor_slim", "*.cpp")))

ext = Extension(
    "actslim._actslim",
    sources=[os.path.join("actslim", "_actslim.cpp")] + vendor,
    include_dirs=[os.path.join("actslim", "_vendor_slim")],
    language="c++",
    extra_compile_args=["-O2", "-std=c++11", "-fpermissive",
                        "-Wno-write-strings", "-Wno-unused-result"],
)

setup(ext_modules=[ext])
