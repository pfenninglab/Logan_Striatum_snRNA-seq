## packages for data table processing 
library(tidyverse)
library(igraph)
library(network)
library(sna)
library(stringr)
library(data.table)
library(leiden)

## this script was adapted from 
## https://github.com/kowaae22/ClarkLabDocumentation/blob/main/MakeigraphFigures.R
## http://kateto.net/network-visualization

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

##############################################################################
#functions to extract genes from fgsea enrich + count matching genes between two lists
numberMatchingGenes=function(genelist1, genelist2){
  m=match(genelist1, genelist2)
  nummatch=length(m[!is.na(m)])
  nummatch
}

plotNetSimple=function(net, ...){
  set.seed(123) ## so the layout is the same
  V(net)$label <- NA
  plot(net, edge.arrow.size=0, edge.color="grey30",  vertex.label.color="#000000", ...)
}

legend.col <- function(col, lev){
  opar <- par
  n <- length(col)
  bx <- par("usr")
  box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
              bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
  box.cy <- c(bx[3], bx[3])
  box.sy <- (bx[4] - bx[3]) / n
  
  xx <- rep(box.cx, each = 2)
  
  par(xpd = TRUE)
  for(i in 1:n){
    yy <- c(box.cy[1] + (box.sy * (i - 1)), box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i)), box.cy[1] + (box.sy * (i - 1)))
    polygon(xx, yy, col = col[i], border = col[i])
  }
  par(new = TRUE)
  plot(0, 0, type = "n",
       ylim = c(min(lev), max(lev)),
       yaxt = "n", ylab = "",
       xaxt = "n", xlab = "",
       frame.plot = FALSE)
  axis(side = 4, las = 2, tick = FALSE, line = .25)
  par <- opar
}

##############################################################################
#function to make pathway network from list of enriched pathways and their genes
# good size functions are sqrt, log10, or some other way to compress a big positive number
make_igraph_from_pathways=function(curenrich,size_function = sqrt, 
                                   type_var = 'celltype', split_char = '#'){
  
  ## the unique list of pathways, merge the genes
  nodes = curenrich %>% group_by(pathway, celltype) %>% 
    ungroup() %>% arrange(padj) %>% 
    ## create the bipartite graph variable
    dplyr::rename('type' = 'celltype') %>% 
    ## create the nodes defined by both pathway and type
    mutate(group = paste0(type, split_char, pathway)) %>%
    ## create size mapping from node values
    distinct(group, leadingEdge,.keep_all = T) %>%
    mutate( size = size_function(size)) %>%
    relocate(group, .before = everything())
  
  ## list of genes per pathway
  genes = nodes %>% dplyr::select(group, leadingEdge) %>% deframe() %>% 
    lapply(str_split,',') %>% lapply(unlist) %>% lapply(unique)
  
  # create the edgelist by the maximum ratio of overlapping genes
  edges = expand.grid(nodes$group,nodes$group) %>%
    ## drop duplicate parallel and self pointing edges
    filter(Pathway1 != Pathway2) %>% 
    mutate(tmp = map2_chr(Pathway1, Pathway2, function(x, y) 
      c(x, y) %>% sort() %>% paste(collapse = ', '))) %>% 
    filter(!duplicated(tmp)) %>% dplyr::select(-tmp) %>%
    ## add the edge weights
    rename_with(~gsub('Var', 'Pathway', .)) %>% 
    mutate(simscore = map2_dbl(genes[Pathway1], genes[Pathway2], numberMatchingGenes), 
           denom = pmin(lengths(genes[Pathway1]), lengths(genes[Pathway2])), 
           overlap_ratio = simscore/denom, 
           ## CHANGE THIS IF NECESSARY: scaling factors for edge values 
           ## (pick one so you can see all edges without them being too bulky)
           width = overlap_ratio^3) %>% 
    dplyr::select(-c(simscore, denom)) %>% 
    filter(overlap_ratio > 0)
  
  ## create the custom igraph object from the pathway data
  net=graph_from_data_frame(d = edges, vertices = nodes, directed = F)
  return(net)
}

## use quantile cutoffs to trim network until target degree is reached
trim_edges = function(net, target_degree = 0.15){
  all_weights = E(net)$overlap_ratio
  quant = 0
  net2 = net
  degree = igraph::degree(net, normalized = T) %>% mean()
  while( degree > target_degree ){
    cut.off <- quantile(all_weights, quant) ## keep top NN% of edges
    net2 <- delete_edges(net2, E(net2)[E(net2)$overlap_ratio<cut.off]) #remove edges below cutoff
    degree = igraph::degree(net2, normalized = T) %>% mean()
    message(paste('Avg Degree is', signif(degree, 3)))
    quant = quant + 0.01
  }
  return(net2)
}

trim_clustered_nodes = function(net, clp, min_nodes = 2){
  #which nodes are in clusters?
  tt=table(clp$membership)
  iimodules=which(tt>=min_nodes)
  inmodules=clp$membership %in% iimodules
  
  #make new network that only contains clusters (no singletons)
  net.clust = delete_vertices(net, V(net)[which(!inmodules)])
  
  ## return list of vertices by cluster
  vl=list()
  mem=clp$membership[inmodules]
  mem=as.numeric(as.factor(mem))
  V(net.clust)$mem = mem
  for(i in 1:max(mem)){
    vl[[i]]=V(net.clust)[mem==i]
  }
  return(list(net = net.clust, vl = vl))
}

## clean up the pathway names
## other split_char could be '^' for start of sentence
clean_pathways = function(curenrich, split_char = '#', width = 35){
  to_delete = curenrich %>% filter(!is.na(MSigDb_Group)) %>% 
    mutate(MSigDb_Group = paste0(split_char, MSigDb_Group)) %>% 
    pull(MSigDb_Group) %>% unique() %>% str_replace_all(':', '|') %>% 
    paste(collapse = '|')
  
  to_label = curenrich %>%
    group_by(mem) %>% top_n(1, jitter(hub_score)) %>% ungroup() %>% 
    mutate(name2 = name %>% 
             str_replace_all('_',' ') %>% 
             str_replace_all(to_delete, '') %>% 
             str_replace_all('#', '\n'),
           name2 = name2 %>% str_trim() %>% str_wrap( width = width)) %>% 
    dplyr::select(name, name2) %>% deframe()
  
  return(to_label)
}


plot_clusterings_3set = function(net_list, prefix, to_label, to_label_num, height = 25, width = 25, ...){
  V(net_list$net)$label <- NA
  n_vert = length(V(net_list$net))
  layout = layout_with_fr(net_list$net)
  layout = norm_coords(layout, xmin = -width/2*.95,  xmax = width/2*.95,
    ymin = -height/2*.95, ymax = height/2*.95)
  
  set.seed(1234)
  pdf(paste0(prefix,'.lab.pdf'), height = height, width = width)
  par(mar = rep(0, 4))
  plot(net_list$net, mark.group=net_list$vl, 
       edge.arrow.size=0, edge.color="grey30", vertex.label.color="#000000", 
       vertex.label=to_label[V(net_list$net)$name],
       xlim = c(-width/2*.95, width/2*.95), 
       ylim = c(-height/2*.95, height/2*.95))
  dev.off()
  
  set.seed(1234)
  pdf(paste0(prefix,'.num.pdf'), height = height, width = width)
  par(mar = rep(0, 4))
  plot(net_list$net, mark.group=net_list$vl, 
       edge.arrow.size=0, edge.color="grey30", vertex.label.color="#000000", 
       vertex.label=to_label_num[V(net_list$net)$name], 
       xlim = c(-width/2*.95, width/2*.95), 
       ylim = c(-height/2*.95, height/2*.95))
  dev.off()
  
  set.seed(1234)
  pdf(paste0(prefix,'.blank.pdf'), height = height, width = width)
  par(mar = rep(0, 4))
  plot(net_list$net, mark.group=net_list$vl, 
       edge.arrow.size=0, edge.color="grey30", vertex.label.color="#000000", 
       xlim = c(-width/2*.95, width/2*.95), 
       ylim = c(-height/2*.95, height/2*.95))
  dev.off()
}



plot_pt_legend = function(vertex.size, prefix, height = 5, width = 5, ...){
  size_range = range(vertex.size)
  sizeCut<- round(seq(size_range[1],size_range[2], diff(size_range)/4 ))
  sizeCutLab <- round(sizeCut**2)
  pdf(paste0(prefix,'.pt.lgd.pdf'), height = height, width = width)
  par(mar = rep(0, 4))
  plot.new()
  legend(x = 0 , y = 0,legend=sizeCutLab, pt.cex=sizeCut/20, 
         pch=21, col='black', pt.bg='black', cex = .5, title = '# Genes',
         bty="n", ncol=1, xjust = 0, yjust = 0)
  dev.off()
}


plot_col_legend = function(prefix, col_fun, height = 5, width = 5, fontsize = 5, ...){
  # if(max_val<min_val) a = max_val; max_val = min_val; min_val = a
  pdf(paste0(prefix,'.col.lgd.pdf'), height = height, width = width)
  par(mar = rep(0, 4))
  lgd = ComplexHeatmap::Legend(
    col_fun = col_fun, grid_width = width*.5, direction = "horizontal", 
    labels_gp = grid::gpar(fontsize = fontsize, fontface = 'bold'))
  ComplexHeatmap::draw(lgd)
  dev.off()
}

