A first example of using pybind11

Create a new subdirectory - e.g. example1 and create the following 5 files in it:

    funcs.hpp
    funcs.cpp
    wrap.cpp
    setup.py
    test_funcs.py

First write the C++ header and implementation files

In funcs.hpp

int add(int i, int j);

In funcs.cpp

int add(int i, int j) {
    return i + j;
};

Next write the C++ wrapper code using pybind11 in wrap.cpp. The arguments "i"_a=1, "j"_a=2 in the exported function definition tells pybind11 to generate variables named i with default value 1 and j with default value 2 for the add function.

#include <pybind11/pybind11.h>
#include "funcs.hpp"

namespace py = pybind11;

using namespace pybind11::literals;

PYBIND11_PLUGIN(wrap) {
    py::module m("wrap", "pybind11 example plugin");
    m.def("add", &add, "A function which adds two numbers",
          "i"_a=1, "j"_a=2);
    return m.ptr();
}

Finally, write the setup.py file to compile the extension module. This is mostly boilerplate.

import os, sys

from distutils.core import setup, Extension
from distutils import sysconfig

cpp_args = ['-std=c++11', '-stdlib=libc++', '-mmacosx-version-min=10.7']

ext_modules = [
    Extension(
    'wrap',
        ['funcs.cpp', 'wrap.cpp'],
        include_dirs=['pybind11/include'],
    language='c++',
    extra_compile_args = cpp_args,
    ),
]

setup(
    name='wrap',
    version='0.0.1',
    author='Cliburn Chan',
    author_email='cliburn.chan@duke.edu',
    description='Example',
    ext_modules=ext_modules,
)

Now build the extension module in the subdirectory with these files

python setup.py build_ext -i

