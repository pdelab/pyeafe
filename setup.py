from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="pyeafe",
    version="0.1.1",
    description="Edge-Averaged Finite Elements (EAFE) for FENiCS",
    url="https://thepnpsolver.github.io",
    author="The PNP Solver",
    author_email="thepnpsolver@gmail.com",
    license="MIT License",
    packages=["pyeafe"],
    zip_safe=False,
    install_requires=[
        "fenics>=2019.1.0",
        "fenics-dolfin>=2019.1.0",
        "numpy>=1.19.0",
        "petsc4py>=3.12.0",
        "typing>=3.7.4.3",
    ],
    long_description=long_description,
    long_description_content_type="text/markdown",
    tests_require=["pytest", "flake8", "black"],
    setup_requires=["pytest-runner"],
    project_urls={"Source Code": "https://github.com/thepnpsolver/pyeafe"},
)
