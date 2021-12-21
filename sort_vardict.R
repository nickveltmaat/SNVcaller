#!/usr/bin/env Rscript

# This script automatically sorts Vardict output file by Chr and Pos
# Author: Nick Veltmaat
# Date: 17-11-2021

args = commandArgs(trailingOnly=TRUE)

# test if there is two arguments passed: if not, return an error
if (length(args)!=2) {
  stop("At least two arguments must be supplied (input file, output file).n", call.=FALSE)
} 

input = args[1]
output = args[2]

df <- read.table(file = toString(input), sep = '\t', header = F)

#Adding 'chr'
df$V2 <- sub("^", "chr", df$V2 )
df$V3 <- sub("^", "chr", df$V3 )

#Order
chrOrder <- c(paste("chr",1:22,sep=""),"chrX","chrY","chrM")
df$V2 <- factor(df$V2, levels=chrOrder)
df <- df[order(df$V2, df$V4, df$V5),]


df$V2 <- sub("chr", "", df$V2 )
df$V3 <- sub("chr", "", df$V3 )
write.table(df, toString(output), sep="\t",row.names=FALSE, col.names = FALSE, quote = FALSE)
