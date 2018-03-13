import sys
from setuptools import setup

try:
    import dolfin
except ImportError:
    print"PyEAFE requires dolfin to be installed,"
    print "please install dolfin before installing PyEAFE."
    sys.exit()

setup(name='pyeafe',
      version='0.1',
      description='Edge-Averaged Finite Elements (EAFE) for FENiCS',
      url='https://github.com/thepnpsolver/pyeafe',
      author='Maximilian Metti, Shuonan Wu, Arthur Bousquet',
      author_email='',
      license='GNU General Public License v3.0',
      packages=['pyeafe'],
      zip_safe=False)
