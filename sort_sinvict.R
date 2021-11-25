# This script automatically sorts Sinvict output file by Chr and Pos
# Author: Nick Veltmaat
# Date: 25-11-2021

df <- read.table(file = './temp/output-sinvict/calls_level1.sinvict', sep = '\t', header = F)

#Adding 'chr'
df$V1 <- sub("^", "chr", df$V1 )

#Order
chrOrder <- c(paste("chr",1:22,sep=""),"chrX","chrY","chrM")
df$V1 <- factor(df$V1, levels=chrOrder)
df <- df[order(df$V1, df$V2),]


df$V1 <- sub("chr", "", df$V1 )
write.table(df,'./temp/output-sinvict/calls_level1_sorted.sinvict',sep="\t",row.names=FALSE, col.names = FALSE, quote = FALSE)
