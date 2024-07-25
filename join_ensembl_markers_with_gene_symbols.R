
setwd("~/git/babs-work-experience/")

library(tidyverse)


mg = read.csv("marker_genes.csv")

bm = read.delim("mart_export.txt")

mg |> 
  pivot_longer(cols=-1, names_to = "cluster", values_to = "ensembl_gene") |> 
  arrange(ensembl_gene) -> mg_long

mg_long |> 
  full_join(y=bm, by=c("ensembl_gene" = "Gene.stable.ID")) -> mg_long_with_bm


write_csv(mg_long_with_bm, "mg_long_with_bm.csv")


