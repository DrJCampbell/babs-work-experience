# Re-processing & annotating a published endometrium dataset

**Dataset.** Non-pregnant human uterus (endometrium, epithelial), Garcia-Alonso *et al.*,
*Nature Genetics* 2021 — [doi.org/10.1038/s41588-021-00972-2](https://doi.org/10.1038/s41588-021-00972-2),
downloaded from [CELLxGENE](https://cellxgene.cziscience.com/).

## Setup

Prerequisite: [install uv](https://docs.astral.sh/uv/getting-started/installation/).


```bash
git clone https://github.com/DrJCampbell/babs-work-experience.git
cd babs-work-experience/2026/Fabian
uv sync
```

`uv sync` reads `pyproject.toml` + `uv.lock` and installs the exact package versions this
project needs into a local `.venv/`.

Then open `fabian_endometrium_project.ipynb` and select `.venv/bin/python` as the kernel
(VS Code: **Select Kernel → Python Environments**). The dataset itself downloads
automatically on first run via `pooch`.

![Their UMAP vs my attempt, same CellTypist labels](figures/start_vs_finish.png)

![my UMAP coloured by different labels](figures/What-drives-the-UMAP-islands.png)

![Their UMAP coloured by different labels](figures/on-the-authors-original-embedding.png)