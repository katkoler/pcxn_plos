# Aggregate results into a single data frame

rm(list = ls(all=TRUE)) 
options(stringsAsFactors = F)


# empty variable to store results
pcxn = c()
pb = txtProgressBar(min=0,max=2,initial=0,style=3)
for(k in 1:2){
  pcxn = rbind(pcxn, readRDS(paste0("pcxn/output/mean_pcor2_barcode/res/pcxn_mean_pcor2_barcode_part",k,".RDS")))
  setTxtProgressBar(pb,k)
}
close(pb)

# adjust p-values for multiple comparison
pcxn$p.Adjust = p.adjust(p = pcxn$p.value, method = "fdr")

# save results
saveRDS(pcxn, "/pcxn/output/pcxn_mean_pcor2_barcode_hallmark.RDS")
