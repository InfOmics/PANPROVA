library("treeio")
library("ggtree")

nwk = file.path("./tree/","Hao.ffn.cv5.nwk")

tree <- read.tree(nwk)
#ggtree(tree,layout="circular")

#ggtree::ggtree(tree, layout = "circular")+ ggtree::geom_tiplab2(offset=0.1, align = F, size=4)

ggtree::ggtree(tree, layout = "circular")+  geom_text2(aes(subset=!isTip, label=label), hjust=-.3, vjust=-.5) + ggtree::geom_tiplab2(offset=0, align = F, size=3)

