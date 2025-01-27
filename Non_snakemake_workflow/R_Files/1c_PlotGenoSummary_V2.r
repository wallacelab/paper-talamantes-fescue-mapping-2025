#! /usr/bin/Rscript

# Test for distortion in the kmers made from each sequence

# Libraries
library(argparse)
library(ggplot2)
library(gridExtra)
library(dplyr)

# Command-line arguments
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-s", "--sitefile", help="TASSEL genotype summary file of sites")
parser$add_argument("-t", "--taxafile", help="TASSEL genotype summary file of taxa")
parser$add_argument("-d", "--depthfile", help="VCFTools depth output file (from --site-mean-depth)")
parser$add_argument("-g", "--genodepth", help="VCFTools output of genotype depth (from --geno-depth, stripped of header & CHROM and POS columns)")
parser$add_argument("-r", "--restriction-frags", help="File with lengths of restriction fragments")
parser$add_argument("-o", "--outfile", help="Output graphic")
parser$add_argument("-m", "--max-datapoints", type="integer", help="Randomly subset to this many datapoints if dataset has more")
parser$add_argument("--min-fragsize", type="integer", default=50, help="Randomly subset to this many datapoints if dataset has more")
parser$add_argument("--max-fragsize", type="integer", default=500, help="Randomly subset to this many datapoints if dataset has more")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Misc/AddissuAyele_QuinoaGBS/2_CallSnps/')
# args=parser$parse_args(c("-s","2c_sitesummary.txt", "-t", "2c_taxasummary.txt",
#                          "-d", "2c_site_depth.txt", "-g", "2c_geno_depth.txt.gz",
#                           "-o", "99_tmp.png"))


plots = list()


###
# Mean site depth summary 
###

if(!is.null(args$depthfile)){
  
  cat("Processing site depth data from",args$depthfile, "\n")
  depthvals = read.delim(args$depthfile, check.names=F)$SUM_DEPTH #read in 
  if(!is.null(args$max_datapoints) && (length(depthvals) > args$max_datapoints)){ # Subset if needed
    cat("\tRandomly subsetting to",args$max_datapoints,"data points\n")
    depthvals = sample(depthvals, args$max_datapoints, replace=F)
  }
  depthvals.sort = sort(depthvals) # Sort for plotting
  
  
  # Depth plots in various formats
  depth_raw = qplot(1:length(depthvals.sort), depthvals.sort, color=I('darkgoldenrod')) +
    labs(title = "Sites - Mean Read Depth", x="Site # (sorted)", y="Read Depth") 
  
  depth_log = qplot(1:length(depthvals.sort), log10(depthvals.sort), color=I('darkgoldenrod4')) +
    labs(title = "Sites - Mean Read Depth (Log scale)", x="Site # (sorted)", y="Log10 Read Depth") 
  
  depth_hist = qplot(log10(depthvals.sort), fill=I('khaki4'), geom="histogram", bins=50) +
    labs(title = "Sites - Mean Read Depth (Log scale)", x="Log10 Read Depth", y="# Sites") 
  
  # Safety check in case all depths are 0 (like VCF converted from Hapmap)
  if(all(depthvals==0)){
    depth_log=NULL
    depth_hist=NULL
  }
  
  plots = c(plots, list(depth_raw, depth_log, depth_hist))
  
  
}else{
  cat("Unable to create site depth plot because missing input file\n")
}


###
# Genotype depths
###

if(!is.null(args$genodepth)){
  
  cat("Processing genotype depth data from",args$genodepth, "\n")
  
  # Read the file while skipping the header (skip = 1 for header) and then removing the first two columns
  genodepths = read.table(args$genodepth, header = TRUE, skip = 1)
  
  # Remove the first two columns (non-numeric)
  genodepths = genodepths[, -c(1, 2)]
  
  # Convert to a numeric vector
  genodepths = as.vector(unlist(genodepths))
  
  # Subset (if needed, which it probably is for this)
  if(!is.null(args$max_datapoints) && (length(genodepths) > args$max_datapoints)){ # Subset if needed
    cat("\tRandomly subsetting to",args$max_datapoints,"data points\n")
    genodepths = sample(genodepths, args$max_datapoints, replace=F)
  }
  
  # Genotype depth as a histogram (probably most useful)
  geno_depth = qplot(genodepths, fill=I('darkorchid1'), geom="histogram", bins=100) +
    labs(title = "Genotypes - Read Depth", x="Read Depth", y="# Genotype Calls") 
  geno_depthlog = qplot(log10(genodepths), fill=I('darkorchid4'), geom="histogram", bins=50) +
    labs(title = "Genotypes - Read Depth (Log scale)", x="Log10 Read Depth", y="# Genotype Calls") 
  
  # Safety check in case all depths are 0 (like VCF converted from Hapmap)
  if(all(genodepths==0)){
    geno_depthlog=NULL
  }
  
  plots = c(plots, list(geno_depth, geno_depthlog))
  
}else{
  cat("Unable to create genotype depth plot because missing input file\n")
}



###
# Site summary
###

if(!is.null(args$sitefile)){
  
  cat("Processing site data from",args$sitefile, "\n")
  
  sites = read.delim(args$sitefile, check.names=F)
  
  if(!is.null(args$max_datapoints) && nrow(sites) > args$max_datapoints){ # Subset if needed
    cat("\tRandomly subsetting to",args$max_datapoints,"data points\n")
    sites = sites[sample(1:nrow(sites), args$max_datapoints, replace=F), ]
  }
  
  
  # Percent Heterozygous; for some reason qplot() was not handling this properly so did full ggplot()
  hetvals = sort(sites[,'Proportion Heterozygous'])
  hetdata=data.frame(xval=1:length(hetvals), yval=hetvals)
  hets = ggplot(hetdata, mapping=aes(x=xval, y=yval)) +
    geom_point(color="darkred") +
    labs(title="Sites - Percent Heterozygous", x="Sites (sorted)", y="Fraction Het")
  
  # Percent Missing
  missing = qplot(sites[,'Proportion Missing'], bins=50, fill=I('darkgreen')) +
    labs(title="Sites - Percent Missing", x="Fraction Missing", y="# Sites")
  
  # Minor allele frequency
  maf = qplot(x=sites[,'Minor Allele Frequency'], bins=50, fill=I('darkblue')) +
    labs(title="Sites - Minor Allele Frequency", x="MAF", y="# Sites")
  
  
  plots = c(plots, list(hets, missing, maf))
  
}else{
  cat("Unable to create site summary plot because missing input file\n")
}



###
# Expected size/count of fragments
###


if(!is.null(args$restriction_frags)){
  cat("Plotting expected restriction fragment length based on sizes in",args$restriction_frags,"\n") 
  frags = scan(args$restriction_frags)
  
  num_captured = sum(frags >= args$min_fragsize & frags <= args$max_fragsize)
  highlight = data.frame(xmin=args$min_fragsize, xmax=args$max_fragsize, ymin=-Inf, ymax=Inf)
  
  fragplot = qplot(x=log10(frags), geom='histogram') +
    labs(title="Restriction fragment distribution", x="Log10 fragment size", y="Count") +
    geom_rect(highlight, mapping=aes(xmin=log10(xmin), xmax=log10(xmax), ymin=ymin, ymax=ymax, 
                                     fill='captured'), inherit.aes=F, alpha=0.25) +
    labs(subtitle=paste("Predicted",num_captured,"fragments captured of",sum(frags),"bp")) +
    labs(fill=paste("Size of ", args$min_fragsize, "-", args$max_fragsize, "bp", sep="")) +
    theme(legend.position = c(0.8, 0.8))
  
  
  plots = c(plots, list(fragplot))
}else{
  cat("Unable to create graph of expected restriction fragments because missing length file\n")
}




###
# Binplot/heatmap of mean depth vs % present; arguably most informative
###

cat( length(depthvals), "\n" )              # Check the length of depthvals
cat(nrow(sites),"\n" )                    # Check how many rows are in the sites data frame
cat(length(sites$'Number of Taxa'),"\n" )  # Check the length of Number of Taxa
cat(length(sites$'Proportion Missing'),"\n" ) # Check the length of Proportion Missing

if(!is.null(args$depthfile) && !is.null(args$sitefile)){
    cat("Combining depth and site info for 2D bin plot\n") 
    depthvals = read.delim(args$depthfile, check.names=F)$SUM_DEPTH #read in 
    sites = read.delim(args$sitefile, check.names=F)
    
    # Binplot of depths
    depthstats = data.frame(mean_depth = depthvals / sites$'Number of Taxa',
                            fraction_present = 1 - sites$'Proportion Missing')

    ## Binplot - Everything
    bindepth = ggplot(depthstats, mapping=aes(x=log10(mean_depth), y=fraction_present)) +
        geom_bin2d() +
        scale_fill_gradient(low='blue', high='red') +
        labs(title = "Geno Calls - Read Depth vs Presence", x="log10(Mean Depth)", y="Fraction Samples Present") +
        ylim(-0.05, 1.05)
    

    # Binplot - Only things with at least 1 mean depth
    substats = subset(depthstats, depthstats$mean_depth >= 1)
    bindepth_cutoff = ggplot(substats, mapping=aes(x=log10(mean_depth), y=fraction_present)) +
        geom_bin2d() +
        scale_fill_gradient(low='blue', high='red') +
        labs(title = "Geno Calls - Read Depth vs Presence", subtitle="(Mean read >=1)", x="log10(Mean Depth)", y="Fraction Samples Present") +
        ylim(-0.05, 1.05)

    plots = c(plots, list(bindepth, bindepth_cutoff))
}else{
    cat("Unable to create 2D bin plot because missing either site depth or site summary files\n")
}


###
# Taxa summary 
###

if(!is.null(args$taxafile)){
  
  cat("Processing taxa data from",args$sitefile, "\n")
  
  taxa = read.delim(args$taxafile, check.names=F)
  
  # Percent Heterozygous
  hetvals = sort(taxa[,'Proportion Heterozygous'])
  hets = qplot(hetvals, bins=50, fill=I('coral')) +
    labs(title = "Taxa - Percent Heterozygous", x="Fraction Het", y="# Taxa") 
  
  # Percent Missing
  missing = qplot(taxa[,'Proportion Missing'], bins=50, fill=I('aquamarine3')) +
    labs(title="Taxa - Percent Missing", x="Fraction Missing", y="# Taxa")
  
  plots = c(plots, list(hets, missing))
  
}else{
  cat("Unable to create taxa summary plots because missing input file\n")
}


###
# Output
###

nplots = length(plots)
ncol=3
nrow = ceiling(nplots/ncol)

cat("Writing graphic to",args$outfile, "\n")
png(args$outfile, height=3*nrow, width=4*ncol, units="in", res=150)
grid.arrange(grobs=plots, ncol=ncol, nrow=nrow, as.table=TRUE)
dev.off()