# Order of commands:

1. ```c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) assembly.cpp -I/usr/include/eigen3 -I/usr/lib/python3/dist-packages/ffc/backends/ufc/  -o pyeafe$(python3-config --extension-suffix) -L/usr/lib/x86_64-linux-gnu/ -l:libdolfin.so```
2. ```python test.py```