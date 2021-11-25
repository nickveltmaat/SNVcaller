# This script automatically sorts Vardict output file by Chr and Pos
# Author: Nick Veltmaat
# Date: 17-11-2021

df <- read.table(file = './temp/VD/vardict_raw.vcf', sep = '\t', header = F)

#Adding 'chr'
df$V2 <- sub("^", "chr", df$V2 )
df$V3 <- sub("^", "chr", df$V3 )

#Order
chrOrder <- c(paste("chr",1:22,sep=""),"chrX","chrY","chrM")
df$V2 <- factor(df$V2, levels=chrOrder)
df <- df[order(df$V2, df$V4, df$V5),]


df$V2 <- sub("chr", "", df$V2 )
df$V3 <- sub("chr", "", df$V3 )
write.table(df,'./temp/VD/vardict_output_sorted.vcf',sep="\t",row.names=FALSE, col.names = FALSE, quote = FALSE)
