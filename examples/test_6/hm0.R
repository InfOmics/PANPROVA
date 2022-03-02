

mfile = file.path("./","dmm")

d = read.csv(mfile, header=TRUE, sep=' ')
d
d = data.matrix(d)
d = t(d)
mean(d)
sd(d)

d.colnames = c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19")
d.rownames = c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19")
d

library("RColorBrewer")
col <- colorRampPalette(brewer.pal(11, "RdYlBu"))(50)
#colors = c(seq(-3,-2,length=100),seq(-2,0.5,length=100),seq(0.5,6,length=100))
#colors = c(seq(0.0,0.499999,length.out=25),seq(0.5,1.0,length.out=26))

colors = unique(c(seq(0.0,0.5,length.out=26),seq(0.5,1.0,length.out=26)))
colors = unique(c(seq(0.0,0.1,length.out=26),seq(0.1,1.0,length.out=26)))


library("gplots")
heatmap.2(d,
          col=col, 
          breaks = colors,
          labCol = FALSE,
          trace="none", 
          #dendrogram = c("col"),
            dendrogram='none', 
          Rowv=FALSE,
          Colv=FALSE,
          denscol="black",
          #          keysize=1,
          key.par=list(mar=c(4,4,4,4)),
          key=T, 
          margins = c(0.5,7))

