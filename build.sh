#!/usr/bin/env bash

set -e

cd "$(dirname "$0")"

echo "==> Installing pybind11 (if needed)"
pip3 install --break-system-packages -q pybind11 flask

PYBIND_INCLUDE=$(python3 -c "import pybind11; print(pybind11.get_include())")
PY_INCLUDE=$(python3 -c "import sysconfig; print(sysconfig.get_path('include'))")
EXT_SUFFIX=$(python3 -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))")

echo "==> Compiling C++ engine -> calcengine${EXT_SUFFIX}"
g++ -O3 -Wall -shared -std=c++17 -fPIC \
    -I"$PYBIND_INCLUDE" -I"$PY_INCLUDE" \
    cpp/ast.cpp \
    cpp/parser.cpp \
    cpp/printing.cpp \
    cpp/differentiation.cpp \
    cpp/integration.cpp \
    cpp/vector_calc.cpp \
    cpp/ode_solver.cpp \
    cpp/limit.cpp \
    cpp/solver.cpp \
    cpp/series.cpp \
    cpp/laplace.cpp \
    cpp/complex.cpp \
    cpp/fourier.cpp \
    cpp/bindings.cpp \
    -o "python_app/calcengine${EXT_SUFFIX}"

echo "==> Build complete: python_app/calcengine${EXT_SUFFIX}"
echo "    Run the app with:  cd python_app && python3 app.py"
