# PCxN
Version submitted to PLOS Computational Biology

The scripts to get the pathway correlation estimates are in the `src` folder. Below is a brief overview of the scripts.

1. **Update gene set annotations**: Keep only genes present in the gene expression background. (`mean_pcor2_barcode_hallmark_estimates00.R`)

2. **Get experiment-level estimates**: For each experiment, estimate all pairwise pathway correlation coefficients along with the corresponding p-values. (`mean_pcor2_barcode_hallmark_estimates01.R`, `mean_pcor2_barcode_hallmark_estimates01.sh`)

3. **Combine experiment-level estimates**: Aggregate the experiment-level correlation estimates using a weighted average, and the experiment-level p-values using Liptak's method. (`mean_pcor2_barcode_hallmark_estimates02.R`, `mean_pcor2_barcode_hallmark_estimates02.sh`)

4. **Aggregate results into a single data frame** (`mean_pcor2_barcode_hallmark_estimates03.R`)
