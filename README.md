# Amuse_Project
This is the github of the AMUSE project *Merging of quenched dwarf galaxies* by Kristian Knudsgaard, Lea Krarup and Carol Dong.

The group report can be found in the reports directory, slides of both presentations can be found in the presentations directory.

The code used to produce the data is in the src directory, while the scripts used to analyze the data are located in the scripts directory. The data used in the report is located in the data directory.

## How to run our code
In order to run a merger simulation, one needs to run the **merger.py** code in the src directory. The parameters of this can be changed using the option parser options.  Most notably, one can add an initial velocity (-v) to the two galaxies, which is the parameter we varied in our project. Additionally, one can make the code run faster by reducing the number of particles (--n_bulge for stars, --n_cloud for gas).
This produces a folder of hdf5 files, that can be read and turned into density plots using the **read_and_plot_data.py** script in the scripts directory. When running the script, one has to specify the path to the data directory to read it using --dir.

## Group rubric
The research question that we agreed on is:
"Do mergers of two quenched dwarf galaxies trigger star formation?"

The agreed upon minimum requirement was to be able to model galaxies made up of DM/stellar particles and gas, and merge them.

The expected figures were a snapshot of the gas before and after the merger, as well as a histogram showing the density distribution of the gas.

Additionally, to improve our grade we could vary the impact of the galaxies, to see how angular momentum and velocity influence the SFR after the merger.

## Project Organization

```
├── AUTHORS.md              <- List of developers and maintainers.
├── CHANGELOG.md            <- Changelog to keep track of new features and fixes.
├── CONTRIBUTING.md         <- Guidelines for contributing to this project.
├── Dockerfile              <- Build a docker container with `docker build .`.
├── LICENSE.txt             <- License as chosen on the command-line.
├── README.md               <- The top-level README for developers.
├── configs                 <- Directory for configurations of model & application.
├── data
│   ├── external            <- Data from third party sources.
│   ├── interim             <- Intermediate data that has been transformed.
│   ├── processed           <- The final, canonical data sets for modeling.
│   └── raw                 <- The original, immutable data dump.
├── docs                    <- Directory for Sphinx documentation in rst or md.
├── environment.yml         <- The conda environment file for reproducibility.
├── models                  <- Trained and serialized models, model predictions,
│                              or model summaries.
├── notebooks               <- Jupyter notebooks. Naming convention is a number (for
│                              ordering), the creator's initials and a description,
│                              e.g. `1.0-fw-initial-data-exploration`.
├── pyproject.toml          <- Build configuration. Don't change! Use `pip install -e .`
│                              to install for development or to build `tox -e build`.
├── references              <- Data dictionaries, manuals, and all other materials.
├── reports                 <- Generated analysis as HTML, PDF, LaTeX, etc.
│   └── figures             <- Generated plots and figures for reports.
├── scripts                 <- Analysis and production scripts which import the
│                              actual PYTHON_PKG, e.g. train_model.
├── setup.cfg               <- Declarative configuration of your project.
├── setup.py                <- [DEPRECATED] Use `python setup.py develop` to install for
│                              development or `python setup.py bdist_wheel` to build.
├── src
│   └── Amuse_Project       <- Actual Python package where the main functionality goes.
├── tests                   <- Unit tests which can be run with `pytest`.
├── .coveragerc             <- Configuration for coverage reports of unit tests.
├── .isort.cfg              <- Configuration for git hook that sorts imports.
└── .pre-commit-config.yaml <- Configuration of pre-commit git hooks.
```

<!-- pyscaffold-notes -->

## Note

This project has been set up using [PyScaffold] 4.5 and the [dsproject extension] 0.7.2.

[conda]: https://docs.conda.io/
[pre-commit]: https://pre-commit.com/
[Jupyter]: https://jupyter.org/
[nbstripout]: https://github.com/kynan/nbstripout
[Google style]: http://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings
[PyScaffold]: https://pyscaffold.org/
[dsproject extension]: https://github.com/pyscaffold/pyscaffoldext-dsproject
