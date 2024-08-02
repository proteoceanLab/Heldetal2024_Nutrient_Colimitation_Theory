# Compile using python3 setup_fast_helper_functions.py build_ext --inplace

from setuptools import setup
from Cython.Build import cythonize

setup(
    ext_modules=cythonize("fast_helper_functions.pyx"),
    zip_safe=False,
)
