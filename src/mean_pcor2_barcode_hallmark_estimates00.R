# Update Gene Set Annotation:
#
# Keep only genes present in the gene expression background


rm(list=ls())
options(stringsAsFactors = F)

cmd_args <- commandArgs(trailingOnly = T)

# ==== Hallmark Gene Sets ===
# Gene set annotation from MSigDB
fname = paste0("../data/", cmd_args[1], sep = "")
h_gs = readRDS(fname)

# ==== Background Genes ====
# Gene present in the microarrays from the gene expression background
background_genes = readRDS("../data/barcode_genes.RDS")

# for each gene set annotation, keep only the genes present in the gene expression background
my_gs = lapply(h_gs,function(x){x[x %in% background_genes]})

# save filtered gene set annotations
saveRDS(my_gs,fname)
