# Aggregate results into a single data frame

rm(list = ls(all=TRUE)) 
options(stringsAsFactors = F)

cmd_args <- commandArgs(trailingOnly = T)

# ==== INPUTS ====
geneset_file <- cmd_args[1]
output_folder <- cmd_args[2]

# ==== Pathway Annotation ====
# Filtered gene set annotation
gs_lst = readRDS(paste("../data/",geneset_file, sep=""))


# indices for pathway pairs 
number_of_pathways = choose(length(gs_lst),2)
# split pathway pairs in chunks
pairs_chunks <- split(1:number_of_pathways, ceiling(1:number_of_pathways/1000))

# empty variable to store results
pcxn = c()
pb = txtProgressBar(min=0,max=2,initial=0,style=3)
for(k in 1:length(pairs_chunks)){
  pcxn = rbind(pcxn, readRDS(paste0("../",output_folder,"/mean_pcor2_barcode/res/pcxn_mean_pcor2_barcode_part",k,".RDS")))
  setTxtProgressBar(pb,k)
}
close(pb)

# adjust p-values for multiple comparison
pcxn$p.Adjust = p.adjust(p = pcxn$p.value, method = "fdr")

# save results
saveRDS(pcxn, paste0("../",output_folder,"/improved_PCxN_",geneset_file))
