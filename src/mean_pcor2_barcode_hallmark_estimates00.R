# Update Gene Set Annotation:
#
# Keep only genes present in the gene expression background


rm(list=ls())
options(stringsAsFactors = F)

# ==== Hallmark Gene Sets ===
# Gene set annotation from MSigDB
h_gs = readRDS("pcxn2017/Data/Gene Sets/h_gs_v5.1.RDS")

# ==== Background Genes ====
# Gene present in the microarrays from the gene expression background
background_genes = readRDS("pcxn/data/barcode_genes.RDS")

# for each gene set annotation, keep only the genes present in the gene expression background
my_gs = lapply(h_gs,function(x){x[x %in% background_genes]})

# save filtered gene set annotations
saveRDS(my_gs,"pcxn/data/h_v51.RDS")
