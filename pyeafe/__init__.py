from pyeafe.assembly import Coefficient, eafe_assemble

__version__ = "0.1.1"

try:
    import dolfin
    if int(dolfin.__version__.split('.')[0])<2019:
        print("WARNING Pyeafe needs Dolfin > 2019 and the version install is {}".format(dolfin.__version__))
        print("Update your version of fenics to >= 2020".format(dolfin.__version__))
        print("See https://fenicsproject.org/download/")
        print("or download our Docker https://github.com/thepnpsolver/pdelab")
except ModuleNotFoundError:
    print('Dolfin (FENiCS) is not install')
    exit()

__all__ = [Coefficient, eafe_assemble]