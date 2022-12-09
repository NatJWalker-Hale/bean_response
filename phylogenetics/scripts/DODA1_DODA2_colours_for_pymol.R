require(ggplot2)

setwd("~/Dropbox/cary_projects/bean_rebuttal/")

# read divergences
df <- read.csv("bean_response/phylogenetics/data/DODA1_DODA2_jsd.csv", header=T)
str(df)

# read site correspondences
corres <- read.table("bean_response/phylogenetics/data/alignment_BvDODAa2_corres.tsv", header=T)
str(corres)

# plot to get colour palette 
p <- ggplot(data=df, aes(x=pos, y=dist)) + 
  geom_point(aes(color=dist)) +
  scale_color_gradient(low="#E6E6E6", high="#FFCC00")

p

g <- ggplot_build(p)

# output df
outdf <- data.frame(pos=df$pos,
                    hex=g$data[[1]]$colour)

# add rgb from hex for pymol
tmp <- t(col2rgb(outdf$hex))

outdf$r <- tmp[,1]
outdf$g <- tmp[,2]
outdf$b <- tmp[,3]

# order according to position
outdf <- outdf[order(outdf$pos),]

# allows us to add BvDODAa2 site correspondence
outdf$BvDODAa2_pos <- corres$sequence

# go back to divergence ordering
outdf <- outdf[match(df$pos, outdf$pos),]

# make resulting df in terms of BvDODAa2
outdf$pos <- outdf$BvDODAa2_pos
outdf <- outdf[,-6]

write.table(outdf, "BvDODAa2_col_out.tsv", sep = "\t", row.names = F, quote = F)
write.table(outdf[c(1:41),], "BvDODAa2_col_out_top41.tsv", row.names = F, quote = F)
