pacman::p_load(ape, ggdendro, dendextend)
setwd("/data/DP_project/04_phylogeny")
mytree <- read.tree('PromoterConsensusSequences.phy')
# get pairwise taxa-taxa distance matrix
d=cophenetic(mytree)
plot(mytree, edge.width = 2)

## hclust euclidean
dist_mat <- dist(d, method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')
plot(hclust_avg)
cut_avg <- cutree(hclust_avg, k = 3)
plot(hclust_avg)
rect.hclust(hclust_avg , k = 3, border = 2:6)

## hclust manhattan
dist_mat.manhat <- dist(d, method = "manhattan")
hclust_avg.manhat <- hclust(dist_mat.manhat, method = 'complete')
plot(hclust_avg.manhat)
cut_avg.manhat <- cutree(hclust_avg.manhat, k = 3)
plot(hclust_avg.manhat)
rect.hclust(hclust_avg.manhat , k = 3, border = 2:6)

## Prune
pruned_tree = drop.tip(mytree, c("ProDP-11"))
d.prune=cophenetic(pruned_tree)
plot(pruned_tree, edge.width = 2)

dist_mat.prune <- dist(d.prune, method = 'euclidean')
hclust_avg.prune <- hclust(dist_mat.prune, method = 'average')
plot(hclust_avg.prune)
cut_avg.prune <- cutree(hclust_avg.prune, k = 3)
plot(hclust_avg.prune, horiz = T)
rect.hclust(hclust_avg.prune , k = 3, border = 2:6)

## hclust manhattan
dist_mat.manhat.prune <- dist(d.prune, method = "manhattan")
hclust_avg.manhat.prune <- hclust(dist_mat.manhat.prune, method = 'complete')
plot(hclust_avg.manhat.prune) 
cut_avg.manhat.prune <- cutree(hclust_avg.manhat.prune, k = 3)
plot(hclust_avg.manhat.prune) 
rect.hclust(hclust_avg.manhat.prune , k = 3, border = 2:6)

## try rotation the plots
hcd <- as.dendrogram(hclust_avg.manhat.prune)
plot(hcd, type = "rectangle", ylab = "Height", horiz = T, main = "Manhattan clustering")

# or 
# library ggdendro
ggdendrogram(hclust_avg.manhat.prune)


## Final plots
# png(filename="figure.png", width=900, bg="white")
svg(filename="figure.svg", width = 14, height = 7, bg="white")
par(mfrow = c(1, 2), mar=c(2,2,2,5)+.1)
## Euclidean
dend <- dist_mat.prune %>%  
  hclust(method = "ward.D2") %>% # Hierarchical clustering 
  as.dendrogram # Turn the object into a dendrogram.

# plot(dend, main = "Euclidean clustering")

dend %>% set("branches_k_color", 
             k = 4, 
             value = c("red2", "blue", "orange", "purple")
             ) %>% 
  plot(horiz = TRUE, main = "Euclidean clustering")
dend %>% 
  rect.dendrogram(k = 4, horiz = TRUE, border = 8, lty = 5, lwd = 2)

## Manhattan
dend <- dist_mat.manhat.prune %>%  
  hclust(method = "ward.D2") %>% # Hierarchical clustering 
  as.dendrogram # Turn the object into a dendrogram.
# plot(dend, main = "Manhattan clustering")

dend %>% set("branches_k_color", 
             k = 4, 
             value = c("red2", "blue", "orange", "purple")
             ) %>% 
  plot(horiz = TRUE, main = "Manhattan clustering")
dend %>% 
  rect.dendrogram(k = 4, horiz = TRUE, border = 8, lty = 5, lwd = 2)

dev.off()

# -------------------------------------------------------------------------------------------------
## Include outlier
# png(filename="figure_withOutlier.png", width=900, bg="white")
svg(filename="figure_with_Outlier.svg", width = 14, height = 7, bg="white")
par(mfrow = c(1, 2), mar=c(2,2,2,5)+.1)
## Euclidean
dend <- dist_mat %>%  
  hclust(method = "ward.D2") %>% # Hierarchical clustering 
  as.dendrogram # Turn the object into a dendrogram.
# plot(dend, main = "Euclidean clustering")

dend %>% set("branches_k_color", 
             k = 5, 
             value = c("red2", "blue", "orange", "purple", "green")
             ) %>% 
  plot(horiz = TRUE, main = "Euclidean clustering")
dend %>% 
  rect.dendrogram(k = 5, horiz = TRUE, border = 8, lty = 5, lwd = 2)

## Manhattan
dend <- dist_mat.manhat %>%  
  hclust(method = "ward.D2") %>% # Hierarchical clustering 
  as.dendrogram # Turn the object into a dendrogram.
# plot(dend, main = "Manhattan clustering")

dend %>% set("branches_k_color", 
             k = 5, 
             value = c("red2", "blue", "orange", "purple", "green")
             ) %>% 
  plot(horiz = TRUE, main = "Manhattan clustering")
dend %>% 
  rect.dendrogram(k = 5, horiz = TRUE, border = 8, lty = 5, lwd = 2)

dev.off()

