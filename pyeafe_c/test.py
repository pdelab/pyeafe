import pyeafe
from dolfin import *


A = Matrix()
print(A.array())
mesh = UnitSquareMesh(8, 8)
pyeafe.poisson_matrix(A,mesh)
print(A.array())