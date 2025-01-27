#QDNAseq example

# may need to increase the max size of objects for future package to work properly
# options(future.globals.maxSize = size_in_mb*1024^2)

library(QDNAseq)
library(QDNAseq.hg38)
library(future)

args <- commandArgs(trailingOnly = TRUE)
outdir <- args[2]

print(outdir)
#set up futures for parallel processing
plan("multisession", workers=10)

# generate bin annotations
bins <- getBinAnnotations(binSize=15, genome="hg38")

# process bamfiles
bam_name <- basename(args[1])

# give chunkSize a non-NULL non-int value to chunk by chromosome
readCounts <- binReadCounts(bins, bamfiles=args[1])

# Plot raw read counts, and highlight bins to filter
png(filename=sprintf("%s/%s.raw_reads.png", outdir, bam_name))
plot(readCounts, logTransform=FALSE, ylim=c(-50, 200))
dev.off()

png(filename=sprintf("%s/%s.raw_reads_filter.png", outdir, bam_name))
highlightFilters(readCounts, logTransform=FALSE,
                 residual=TRUE, blacklist=TRUE)
dev.off()

# Apply filters and plot median read counts as a function of GC content and mappability
readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)
png(filename=sprintf("%s/%s.isobar_plot.png", outdir, bam_name))
isobarPlot(readCountsFiltered)
dev.off()

# Estimate the correction for GC content and mappability, and make a plot for 
# the relationship between the observed standard deviation in the data and its read depth
readCountsFiltered <- estimateCorrection(readCountsFiltered)
png(filename=sprintf("%s/%s.noise_plot.png", outdir, bam_name))
noisePlot(readCountsFiltered)
dev.off()

# apply corrections for GC content and mappability and plot copy number profile
copyNumbers <- correctBins(readCountsFiltered)
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
png(filename=sprintf("%s/%s.CN_plot.png", outdir, bam_name))
plot(copyNumbersSmooth)
dev.off()

copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
png(filename=sprintf("%s/%s.CN_segmented_plot.png", outdir, bam_name))
plot(copyNumbersSegmented)
dev.off()

copyNumbersCalled <- callBins(copyNumbersSegmented)
png(filename=sprintf("%s/%s.CN_called.png", outdir, bam_name))
plot(copyNumbersCalled)
dev.off()

exportBins(copyNumbersCalled, 
           file=sprintf("%s/%s_QDNAseq.vcf", outdir, bam_name), 
           format="vcf")
