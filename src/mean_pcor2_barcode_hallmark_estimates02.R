# Combine experiment-level estimates
#
# Aggregate the experiment-level correlation estimates using
# a weighted average, and the experiment-level p-values
# using Liptak's method.

# = Submission Script (Odyssey - slurm) =
#
# mean_pcor2_barcode_hallmark_estimates02.sh [1-2]
# 
# Submit jobs in Odyssey using a jobArray. 
# The pathway pairs were split in two blocks. The first index (SLURM_ARRAY_TASK_ID) 
# corresponds to the pathway pairs block while the second index (SLURM_CPUS_PER_TASK) 
# is the number of cores



rm(list = ls(all=TRUE)) 
gc()
options(stringsAsFactors = F)

library(parallel)
library(metap)

# directory with gene expression background
# availalble in Google Drive)
barcode_dir <- "/shared/hidelab2/shared/Sokratis/PCxN_Plos/data/HGU133plus2/"

# ==== Arguments ====
# read command line arguments
# cmd_args <- as.numeric(commandArgs(trailingOnly = T))
cmd_args <- commandArgs(trailingOnly = T)

# ==== INPUTS ====
geneset_file <- cmd_args[2]
output_folder <- cmd_args[4]
id <- cmd_args[1]

# cmd_args <- as.numeric(cmd_args[c(1,3)])

# ==== Functions ====
# adjust p-values and correlation estimates,
# otherwise the functions to combine the p-values
# cannot handle values near 0 and values near 1
AdjustPmat = function(p_mat,eps=1E-16){
  res = t(apply(p_mat,1,function(pval){
    pval[pval <= eps] = eps
    pval[pval >= 1-eps] = 1 - eps
    return(pval)
  }))
  return(res)
}

# ==== Barcode Annotation ====
# Sample annotation for the gene expression background
tissue_annot <- readRDS( "../data/Barcode3.tissue.RDS" )

# ==== GSE Series ====
# GSE series per tissue
tmp <- subset(tissue_annot,select=c(tissue,series))
gse_lst <- lapply(split(tmp, tmp$tissue),function(x){table(x$series)})
# order by number of samples
gse_lst <- gse_lst[order(sapply(gse_lst,sum),decreasing=T)]
rm(tmp)
# filter series with at least 5 samples
res <- lapply(gse_lst,function(x){x[which(x >= 10)]})
res <- res[lapply(res,length)>0]
# order by number of samples
res <- res[order(sapply(res,sum),decreasing=T)]


fname = c()
for(tissue_select in names(res)){
    # tissue meta-data
    tissue_meta <- subset(tissue_annot[tissue_annot$tissue == tissue_select,],select=c(sample,series))
    # tissue GSE series
    tissue_series <- names(res[[ tissue_select ]])
    tissue_series <- gsub(";","_",tissue_series)
    # RDS file name
    fname = c(fname, paste0(make.names(tissue_select),"_",tissue_series))
}





# ==== Pathway Annotation ====
# Filtered gene set annotation
gs_lst = readRDS(paste("../data/",geneset_file, sep=""))


# indices for pathway pairs 
number_of_pathways = choose(length(gs_lst),2)
# split pathway pairs in chunks
pairs_chunks <- split(1:number_of_pathways, ceiling(1:number_of_pathways/1000))

# ==== Pathway Names =====
# load a set of experiment-level estimates
my_rds = paste0("../",output_folder,"/mean_pcor2_barcode/",fname[13],"_cpad_pathcor.RDS")

myLst = readRDS(my_rds)[ pairs_chunks[[ as.numeric(cmd_args[1]) ]] ]
# get names for pathway pairs
res = as.data.frame( t(sapply(myLst,function(x){ c(x[["Pathway.A"]],x[["Pathway.B"]]) })) )
colnames(res) = c("Pathway.A","Pathway.B")

# ==== Overlap Coefficient =====
# get overlap coefficients for pathway pairs
res$Overlap.Coefficient = unlist( sapply(myLst,function(x){x[[ "Overlap.Coeff" ]]}), use.names = F )

# ===== Correlation Estimates ====
GetCorEstimate = function(ic){
    # load experiment-level estimates
    my_rds = paste0("../",output_folder,"/mean_pcor2_barcode/",fname[ic],"_cpad_pathcor.RDS")
    myLst = readRDS(my_rds)[ pairs_chunks[[ as.numeric(cmd_args[1]) ]] ]
    # extract correlation estimates
    tmp = unlist( sapply(myLst, function(x){x[[ "estimate" ]]}), use.names = F )
    setTxtProgressBar(pb,ic)
    return( tmp )
}

pb = txtProgressBar(min=0,max=length(fname),style=3,initial=0,char="-")
cat("\n")
cor_estimates = mclapply( seq_along(fname), GetCorEstimate, mc.cores = as.numeric(cmd_args[3]) ) 
close(pb)


cor_estimates = Reduce(f = cbind, x = cor_estimates)
colnames(cor_estimates) = fname

# number of samples per experiment
n_vec = readRDS( paste0( "../",output_folder,"/mean_pcor2_barcode/res/n_vec.RDS" ) )

# weighted average for the correlation estimates
n_mult = n_vec/sum(n_vec)
res$PathCor = c( cor_estimates%*%n_mult )

# ==== P-Values ====
GetPvals = function(ic){
    # load experiment-level estimates
    my_rds = paste0("../",output_folder,"/mean_pcor2_barcode/",fname[ic],"_cpad_pathcor.RDS")
    myLst = readRDS(my_rds)[ pairs_chunks[[ as.numeric(cmd_args[1]) ]] ]
    # extract p-value
    tmp = unlist( sapply(myLst, function(x){x[[ "p.value" ]]}), use.names = F )
    setTxtProgressBar(pb,ic)
    return( tmp )
}


pb = txtProgressBar(min=0,max=length(fname),style=3,initial=0,char="-")
cat("\n")
pvals = mclapply( seq_along(fname), GetPvals, mc.cores = as.numeric(cmd_args[3]) ) 
close(pb)

pvals = Reduce(f = cbind, x = pvals)
colnames(pvals) = fname


# adjust p-values, otherwise the functions to combine the p-values
# cannot handle values near 0 and values near 1
pvals = AdjustPmat(pvals)

# Liptak's Method to combine p-values
CombinePval = function(ic){
    if( any( is.na( pvals[ic,] ) ) ){
        tmp = NA
    }else{
        tmp = c(sumz(p=pvals[ic,], weights=n_vec)$p)
    }
    setTxtProgressBar(pb,ic)
    return(tmp)
}

pb = txtProgressBar(min=0,max=length(fname),style=3,initial=0,char="-")
cat("\n")
combined_pvals = mclapply( 1:nrow(pvals), CombinePval, mc.cores = as.numeric(cmd_args[3]) ) 
close(pb)

res$p.value = unlist( combined_pvals )

# save results in data frame
saveRDS(res, paste0("../",output_folder,"/mean_pcor2_barcode/res/pcxn_mean_pcor2_barcode_part",as.numeric(cmd_args[1]),".RDS"))

rm(list = ls())
gc()