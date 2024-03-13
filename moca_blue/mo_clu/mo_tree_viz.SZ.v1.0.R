# Load the required packages
library(ggtree)
library(treeio)

# Set the working directory to where your file is located
setwd("/home/ibg-4/Desktop/Rhome/moca_blue/mo_clu/out")

# Read the Newick file
tree <- read.tree("rdf5_ArthD0D6_cwm-motifs.jaspar-Sandelin-Wassermann.nwk")

# Filter the tree
filtered_tree <- drop.tip(tree, tree$tip.label[!grepl("F_", tree$tip.label)])
#filtered_tree$tip.label <- substr(filtered_tree$tip.label, 1, 17)


p<- ggtree(filtered_tree) +
  geom_tiplab(align = TRUE, linesize = 0.5, size = 2) +
  geom_tippoint(aes(
    color = ifelse(grepl("D6", label), "red", "black")
  )) +
  theme(legend.position = "none") +
  theme(plot.margin = margin(1, 5, 1, 3, "cm")) +
  coord_cartesian(clip="off") 
ggsave("tree.svg", p, width =7, height = 10, units = "in", dpi = 300)





#nodes_to_label <- grep("p0", filtered_tree$tip.label)
#nodelabels(nodes_to_label, pch = 0.2, col = "red")
