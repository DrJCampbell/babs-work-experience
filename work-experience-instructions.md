# Work experience
Using scanpy do this: scRNA-seq pbmc3k

## Table of Contents
1. [Setting up the environment](#setting-up-the-environment)
2. [Test](#test)


### Setting up the environment
Firstly install micromamba with this command
`"${SHELL}" <(curl -L micro.mamba.pm/install.sh)`

Once installed navigate to the directory of your choice and create your project using:
`micromamba create -n [project name] --yes`

To enter your new environment, enter:
`micromamba activate [project name]`

Install the latest version of python and use it to install the required packages.
```micromamba install python

pip install scanpy
pip install leidenalg
pip install jupyter
pip install pandas
pip install numpy```

# test
