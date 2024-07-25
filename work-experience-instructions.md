
# Work experience

Using scanpy to analyse genes found in fetal blood

## Table of Contents

1. [Setting up the environment](#env-setup)
   1. [Installing dependencies](#dependencies)
   2. [Starting the jupyter notebook](#notebook)
3. [Programming the software](#programming-software)
   1. [Obtaining the data](#obtaining-data)


### Setting up the environment <a name="env-setup"></a>

#### Installing dependencies <a name="dependencies"></a>

Firstly install micromamba with this command  
`"${SHELL}" <(curl -L micro.mamba.pm/install.sh)`

Once installed navigate to the directory of your choice and create your project using  
`micromamba create -n [project name] --yes`

To enter your new environment, enter  
`micromamba activate [project name]`

Install the latest version of python and use it to install the required packages  
```sh
micromamba install python

pip install scanpy
pip install leidenalg
pip install jupyter
pip install pandas
pip install numpy
```

#### Starting the jupyter notebook <a name="notebook"></a>

Once you are in your environment, you can start the jupyter notebook by running  
`jupyter notebook`  
in the terminal.

In order to save your command history for later refrence for a write up just like the one I am making right now, put  
`history > [file name]`  
into the terminal.

Once the command has finished executing, look for a log that is purple, it should look something like  
`[C 2024-07-25 10:07:14.500 ServerApp]`  
There should be some urls that are written out, with some instructions to open them.  
__Ctrl/Cmd click__ on the link starting with http[]()://localhost: and a tab should open in your browser.  
This is your main environment gui.

To start programming, click on __New__ then __Python 3 (ipykernel)__

> Some useful shortcuts to know are:
> - 0 0 restarts the kernel
> - d d deletes a codeblock
> - a/b adds a codeblock above/below the currently selected codeblock
> - arrow keys navigate up and down codeblocks
> - shift enter runs the current codeblock and advances forward
> - ctrl/cmd enter runs the current codeblock and does not advance forward
> - alt enter runs the current codeblock and inserts a blank codeblock below

### Programming the software <a name="programming-software"></a>

#### Obtaining the data <a name="obtaining-data"></a>

Once inside your jupyter notebook, download your data of choice by using  
`adata = sc.datasets.ebi_expression_atlas("name of your dataset")`  

In order to find your dataset head to the [single cell expression atlas](https://www.ebi.ac.uk/gxa/sc/experiments) and find a dataset that interests you. Click on the dataset and in the url there should be a unique identifier at the end of the url e.g. sc/experiments/[unique-identifier]

Once you run the code inside your jupyter notebook, the dataset should start downloading. You can find the downloaded files in project/data/[unique-identifyer]/

