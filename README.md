# ChIP-seq Downstream Analysis Tools
Ashley Woodfin

Collection of my downstream analysis scripts to automate ChIP-seq and transcriptional pausing analysis

# Overview

1. createROtracks.R              - Make bigWig files from a bam file for either stranded or unstranded data

2. extendReads.R                 - Extend reads in a bam file and save a bed for either stranded or unstranded data

3. ChIP_pickProteinCodingGenes.R - Pick the most highly occupied protein coding transcripts based on ChIPseq signal in the TSS

4. ChIP_CalcPausingIndex.R       - Calculate the average coverage in a defined promoter and gene body region in order to calculate a pausing index

5. ChIP_plotEcdf.R               - Plot the pausing index or body promoter occupancy ratio from average coverages reported in 4

6. ChIP_metaGeneMatrix.R         - Approximate genes coverages from a bigWig to the same length and save the matrix

7. ChIP_heatmapMatrix.R          - For a window around a Tss, Tes or Peak center extract the coverage and bin into a smaller window size

8. ChIP_makeMetaGenePlot.R       - Plot meta gene average profiles for matrices from 6

9. ChIP_kmeansHeatMap.R          - Cluster heatmap matrices from 6 & 7 for multiple samples with a defined control using kmeans clustering and save into cdt format

10. ChIP_makeMetaPlot.R          - Plot average profiles for feature window matrices from 7

11. ChIP_CalcCovGeneSections.R   - Split a gene into a set number of windows and calculate the total coverage in each window and for the entire gene










