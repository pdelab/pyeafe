# Contributing
The PyEAFE development environment is managed by Anaconda.
If Anaconda is not installed,
please [install Anaconda](https://docs.continuum.io/anaconda/install/)
before proceeding.

To create and activate the virtual environment, run:
```
  conda env create --file environment.yml  # create the pyeafe virtual environment
  conda activate pyeafe
```
To exit the virtual environment,
run `conda deactivate` and the `(pyeafe)` prefix should no longer be visible in the command prompt.


## Updating the virtual environment
When packages are updated, run:
```
  conda env update --file environment.yml
```
to get the latest environment installed.
If new dependencies are added to the virtual environment,
they must be tracked in the `environment.yml` file.
Newly introduced dependencies should be added to the `environment.yml`,
as well as tracked separately in `setup.py`,
where only required dependencies (not dev-dependencies) belong in the `install_requires`
field.

To see the package version installed by conda, run `conda list ${package-name}`
inside the conda environment.

If pre-commit hooks are updated (e.g. new tests require updating dependencies),
run:
```
  pre-commit install && pre-commit run --all-files
```


## Code Styling
Code smells and linting automatically takes place in the pre-commit hook;
however, these operations can always be executed manually with:
```
  black . && flake8
```


## Testing
A suite of unit and functional tests verify the functionality of the repo.
To run the tests, run `pytest` at the top-level directory.

Ensure that all tests pass before pushing to remote.
Failing tests prevent pull requests from being reviewed and merged.


## Releasing
Before building a release,
update the package version in `pyeafe/__init__.py` following
[semantic versioning standards](https://semver.org/).
Tests and code smells should be run before bundling for release:
```
  python setup.py test && python setup.py sdist bdist_wheel
```

If all goes well, push to the git remote and merge into the `master` branch.
Create a release on GitHub with the matching version, and upload to pypi:
```
  twine upload dist/*
```


## Uninstalling
To uninstall, simply remove the `pyeafe` virtual environment with
`conda env remove -n pyeafe`.
