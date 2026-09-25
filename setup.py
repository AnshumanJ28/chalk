"""
Cross-platform build script for the C++ engine, using setuptools + pybind11.
This is the easiest way to build on native Windows (MSVC), since it handles
all the compiler/linker flags automatically -- no manual cl.exe commands needed.

Usage (Windows, from a "Developer PowerShell for VS" so cl.exe is on PATH):
    pip install pybind11 setuptools
    python setup.py build_ext --inplace

Usage (Linux/Mac -- build.sh does the same thing more simply, but this works too):
    python3 setup.py build_ext --inplace
"""
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules = [
    Pybind11Extension(
        "calcengine",
        [
            "cpp/ast.cpp",
            "cpp/parser.cpp",
            "cpp/printing.cpp",
            "cpp/differentiation.cpp",
            "cpp/integration.cpp",
            "cpp/vector_calc.cpp",
            "cpp/ode_solver.cpp",
            "cpp/limit.cpp",
            "cpp/solver.cpp",
            "cpp/series.cpp",
            "cpp/laplace.cpp",
            "cpp/complex.cpp",
            "cpp/fourier.cpp",
            "cpp/bindings.cpp",
        ],
        cxx_std=17,
    ),
]

setup(
    name="calcengine",
    version="2.0",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)
