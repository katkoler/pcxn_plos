# Aggregate results into a single data frame

rm(list = ls(all=TRUE)) 
options(stringsAsFactors = F)

cmd_args <- commandArgs(trailingOnly = T)

# ==== INPUTS ====
geneset_file <- cmd_args[1]
output_folder <- cmd_args[2]


# empty variable to store results
pcxn = c()
pb = txtProgressBar(min=0,max=2,initial=0,style=3)
for(k in 1:1){
  pcxn = rbind(pcxn, readRDS(paste0("../",output_folder,"/mean_pcor2_barcode/res/pcxn_mean_pcor2_barcode_part",k,".RDS")))
  setTxtProgressBar(pb,k)
}
close(pb)

# adjust p-values for multiple comparison
pcxn$p.Adjust = p.adjust(p = pcxn$p.value, method = "fdr")

# save results
saveRDS(pcxn, paste0("../",output_folder,"/improved_PCxN_",geneset_file))
