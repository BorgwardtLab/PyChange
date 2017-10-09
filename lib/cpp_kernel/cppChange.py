
"""
To compile it and run the test, enter:
g++ -std=c++14 -c -fPIC _main.cpp -o MC.o
g++ -shared -Wl,-install_name,MC.so -o MC.so MC.o
python cppChange.py
"""


import time
import numpy as np
import gc
from matplotlib import pyplot as plt
from ctypes import *

# A class for a C++ vector with a changepoint function in Python


def create_lib(n):
    class Vector(object):
        lib = cdll.LoadLibrary('./lib/cpp_kernel/MC.so')  # class level loading lib
        lib.new_vector.restype = c_void_p
        lib.new_vector.argtypes = []
        lib.delete_vector.restype = None
        lib.delete_vector.argtypes = [c_void_p]
        lib.vector_size.restype = c_int
        lib.vector_size.argtypes = [c_void_p]
        lib.vector_get.restype = c_double
        lib.vector_get.argtypes = [c_void_p, c_int]
        lib.vector_push_back.restype = None
        lib.vector_push_back.argtypes = [c_void_p, c_double]
        lib.change.argtypes = [c_void_p, c_char_p, c_double * n, c_double * n, c_int]
        lib.change.restype = None

        def __init__(self):
            self.vector = Vector.lib.new_vector()  # pointer to new vector

        def __del__(self):  # when reference count hits 0 in Python,
            if Vector:
                Vector.lib.delete_vector(self.vector)  # call C++ vector destructor

        def __len__(self):
            return Vector.lib.vector_size(self.vector)

        def __getitem__(self, i):  # access elements in vector at index
            if 0 <= i < len(self):
                return Vector.lib.vector_get(self.vector, c_int(i))
            raise IndexError('Vector index out of range')

        def __repr__(self):
            return '[{}]'.format(', '.join(str(self[i]) for i in range(len(self))))

        def push(self, j):  # push calls vector's push_back
            Vector.lib.vector_push_back(self.vector, c_double(j))

        def change(self, filename, seq, time):  # foo in Python calls foo in C++
            arr = (c_double * len(seq))(*seq)
            timearr = (c_double * len(time))(*time)
            Vector.lib.change(self.vector, c_char_p(filename), arr, timearr, c_int(len(seq)))

    return Vector


def decompose(seq):
    if len(seq) == 2 and len(seq[0]) > 0:  # The case of seq[0]=vals, seq[1]=timepoints
        return seq[0], seq[1]
    if len(seq) > 2 and len(seq[0]) > 0:  # The case of seq[i]=vals_i
        seqarr = []
        timearr = []
        for subseq in seq:
            seqarr.extend(subseq)
            timearr.extend(range(len(subseq)))
        return seqarr, timearr
    else:
        return seq, range(len(seq))


class CppChange:

    def __init__(self, seq, method):
        gc.collect()

        seqarr, timearr = decompose(seq)  # Create or extract timepoints

        Vector = create_lib(len(seq))
        test = Vector()

        start = time.time()
        loc = test.change(method, seqarr, timearr)
        self.duration = time.time() - start

        number = int(test[0])

        self.changepoints = []
        for n in range(number):
            self.changepoints.append(int(test[n + 1]))
        self.p = test[number + 1]
        self.window = int(test[number + 2])
        self.threshold = test[number + 3]

        del test
