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
Writing the list of dependencies to this file should be done by
```
  conda env export --name pyeafe > environment.yml
```

If pre-commit hooks are updated (e.g. new tests require updating dependencies),
run:
```
  pre-commit install && pre-commit run --all-files
```


## Testing
A suite of unit and functional tests verify the functionality of the repo.
To run the tests, run `python -m unittest discover` at the top-level directory.

All tests are required to pass before each commit.
This is automatically verified by use of pre-commit hooks.


## Uninstalling
To uninstall, simply remove the `pyeafe` virtual environment with
`conda env remove -n pyeafe`.
