
# Work experience

Using scanpy to analyse genes found in fetal blood

## Table of Contents

1. [Setting up the environment](#env-setup)
   1. [Installing dependencies](#dependencies)
   2. [Starting the jupyter notebook](#notebook)
3. [Programming the software](#programming-software)
   1. [Obtaining the data](#obtaining-data)
   2. [Processing data that uses gene IDs](#gene-ids)
      1. [Identifying mitochondrial genes using gene IDs](#mitochondrial-gene-ids)
      2. [Plotting using colors with gene IDs](#gene-ids-colors)
      3. [Getting marker genes with gene IDs](#marker-gene-ids)


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

#### Processing data which uses gene IDs <a name="gene-ids"></a>

##### Identifying mitochondrial genes using gene IDs <a name="mitochondrial-gene-ids"></a>
You might encounter some errors when trying to process data that uses gene IDs instead of gene names. In this case, the method for identifying gene symbols that have mitochondria changes as the gene IDs are similar, regardless if the gene symbol has mitochondria e.g.  
when using gene names, you can identify the mitochondrial genes by checking the prefix 'MT-'  
However when using gene IDs, all genes start with ENSG0000...  
This can be fixed by finding a list of all the mitochondrial genes and their gene IDs in a file, reading in the file as a csv and iterating through both the original data and list of mitochondrial genes and identifying any matching gene IDs. This can be achieved with the following code:
```
# this is the line of code that would normally be used if gene names were included
# adata.var["mt"] = adata.var_names.str.startswith("MT-")
# instead since we don't know which genes are mitochondrial (yet) we just init this with
adata.var["mt"] = False

bdata = pd.read_csv("mt_genes_table.txt", sep="\t", header=0).get("ensembl")

for a in adata.var_names:
    for b in bdata:
        if a == b:
            adata.var.loc[a, "mt"] = True
```
where `adata` is the original data and `"mt_genes_table.txt"` is the list of all mitochondrial genes. If a different deliminator is used, the `sep` value can be changed. The parameter `header=0` means a header is included.

##### Plotting using colors with gene IDs <a name="gene-ids-colors"></a>

When it comes to plotting using the gene names as colors such as in `sc.pl.umap(adata, color=["CST3", "NKG7", "PPBP"])`, you can simply replace the gene names with gene IDs, finding out the respective genes IDs from their names- or you could pick random gene IDs that give cool looking results.

##### Getting marker genes with gene IDs <a name="marker-gene-ids"></a>

There will come an unfortunate moment in writing the program, when you will encounter a certain line of code containing 17 different gene names in a list called marker genes:
```
marker_genes = [
   *["IL7R", "CD79A", "MS4A1", "CD8A", "LYZ", "CD14"],
   *["LGALS3", "S100A8", "GNLY", "NKG7", "KLRB1"],
   *["FCGR3A", "MS4A7", "FCER1A", "CST3", "PPBP"],
]
```
While you could also pick out random gene IDs and replace these as marker genes, your results will look better if you take a more reasoned approach. This consists of grabbing the genes from a previously made table: `pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(2)` and exporting it as a csv as shown below:
```
marker_genes = pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(2)

marker_genes.to_csv("marker_genes.csv", header=True) # export table to csv file for use in biomart ensembl
```
Once the file has been exported, you should then head to [biomart](https://www.ensembl.org/info/data/biomart/index.html) and click on __BioMart__ on the blue website header which should take you to the tool for identifying the gene names from gene IDs in bulk.  
In the __Dataset__ tab set choose the correct database and dataset for the data being analysed. In my case the database was Ensembl genes 112 and the dataset was Human Genes (GRCh38.p14)  
In the __Filters__ tab select the file you exported using the above code for the external refrences ID list.  
In the __Attributes__ tab make sure __Features__ is selected and under __GENE__ for __Ensembl__ deselect everything only keeping Gene stable ID and Gene name since this is the information we are interested in.  
Finally on the left side under the website header click results and export the data to a CSV file.

> if you encounter any warnings when running through your code, the best thing to do is ignore them! (don't actually follow this advice) The easiest way to hide the error messages is by including the following line at the top of your code:
> ```
> from warnings import simplefilter
> simplefilter(action="ignore", category=[warning type])
> ```
> where warning type is the type of warning that was displayed
