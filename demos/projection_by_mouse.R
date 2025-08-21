# Plot of line graph to show projection strengths of neurons in each cluster, color-coded by the mouse

projection_by_mouse <- function(Y,
                             mouse_label,
                             Z,
                             region_name,
                             motifs = NULL, ncol=NULL){
  
  
  M <- length(Y)
  
  if(is.null(motifs)){
    motifs = sort(unique(unlist(Z)))
  }
  
  # Projection strength of each neuron
  Y_prop <- lapply(1:M,
                   function(m){
                     
                     matrix(Y[[m]]/rep(colSums(Y[[m]]), each = nrow(Y[[m]])),
                            
                            nrow = nrow(Y[[m]]),
                            ncol = ncol(Y[[m]]))
                   })
  
  # Data frame
  list0 <- lapply(1:M,
                  function(m){
                    
                    data.frame(projection_strength = as.vector(Y_prop[[m]]),
                               region_name = rep(region_name, ncol(Y[[m]])),
                               cluster = rep(Z[[m]], each = nrow(Y[[m]])),
                               #cluster = part_motif_names$pp.regions[part_motif_names$cluster==1]
                               mouse_label = rep(mouse_label[[m]], each = nrow(Y[[m]])))
                  })
  
  df <- do.call(rbind, list0)
  
  # Change neuron into factor
  df$neuron <- factor(rep(1:length(unlist(Z)), each = nrow(Y[[1]])),
                      levels = rep(1:length(unlist(Z))))
  
  # Change cluster into factors
  df$cluster <- factor(paste('cluster', df$cluster),
                       levels = paste('cluster', 1:max(unlist(Z))))
  
  # Change region names into factors
  df$region_name <- factor(df$region_name,
                           levels = region_name)
  
  # Filter to clusters with at least 20 neurons 
  filt = apply(matrix(df$cluster,ncol=1),1, function(x){sum(x== paste('cluster',motifs))>0})
  df = df[filt==1,]
  df$cluster =  factor(df$cluster, levels = unique(df$cluster))
  
  # Line graph
  ggplot(df)+
    geom_line(mapping = aes(x = region_name,
                            y = projection_strength,
                            colour = mouse_label,
                            group = interaction(neuron, mouse_label)),alpha=0.5)+
    facet_wrap(~cluster, ncol = ncol)+
    theme_bw()+
    theme(axis.text.x = element_text(size = 14,angle = 90, vjust = 0.5, hjust = 1))+
    xlab('region')+
    ylab('projection strengths') +
    labs(color='brain')
  
}

projection_names <- function(Y,
                              mouse_label,
                              Z,
                              region_name,
                             thrshold = 0.05){
  
  
  M <- length(Y)
  R <- length(region_name)
  
  # Projection strength of each neuron
  Y_prop <- lapply(1:M,
                   function(m){
                     
                     ps = matrix(Y[[m]]/rep(colSums(Y[[m]]), each = nrow(Y[[m]])),
                            
                            nrow = nrow(Y[[m]]),
                            ncol = ncol(Y[[m]]))
                     df = data.frame(ps_avg=apply(ps, 1, mean))
                     df$cluster_regions = rep(paste(region_name[df$ps_avg>thrshold]), nrow(Y[[m]]))
                     df$cluster_number = rep(paste(m), nrow(Y[[m]]))
                   })
  return( do.call(rbind, Y_prop) )
}

projection_filter <- function(Y,
                                mouse_label,
                                Z,
                                region_name){
  
  
  M <- length(Y)
  
  # Projection strength of each neuron
  Y_prop <- lapply(1:M,
                   function(m){
                     
                     matrix(Y[[m]]/rep(colSums(Y[[m]]), each = nrow(Y[[m]])),
                            
                            nrow = nrow(Y[[m]]),
                            ncol = ncol(Y[[m]]))
                   })
  
  # Data frame
  list0 <- lapply(1:M,
                  function(m){
                    
                    data.frame(projection_strength = as.vector(Y_prop[[m]]),
                               region_name = rep(region_name, ncol(Y[[m]])),
                               cluster = rep(Z[[m]], each = nrow(Y[[m]])),
                               mouse_label = rep(mouse_label[[m]], each = nrow(Y[[m]])))
                  })
  
  df <- do.call(rbind, list0)
  
  # Change neuron into factor
  df$neuron <- factor(rep(1:length(unlist(Z)), each = nrow(Y[[1]])),
                      levels = rep(1:length(unlist(Z))))
  
  # Change cluster into factors
  df$cluster <- factor(paste('cluster', df$cluster),
                       levels = paste('cluster', 1:max(unlist(Z))))
  
  # Change region names into factors
  df$region_name <- factor(df$region_name,
                           levels = region_name)
  
  # Filter to clusters with at least 20 neurons 
  filt = apply(matrix(df$cluster,ncol=1),1, function(x){sum(x== names(table(df$cluster))[table(df$cluster)>10*6])})
  df = df[filt==1,]
  df$cluster =  factor(df$cluster, levels = unique(df$cluster))
  
  df$projection_strength_avg = apply(df,1, function(x){mean(df[df$region_name==x[2]&df$cluster==x[3],1])})
  
  # Line graph
  ggplot(df)+
    geom_line(mapping = aes(x = region_name,
                            y = projection_strength,
                            group = interaction(neuron, mouse_label)))+
    geom_line(mapping = aes(x = region_name,
                             y = projection_strength_avg,
                            group = cluster),
                col="red")+
    facet_wrap(~cluster)+
    theme_bw()+
    xlab('region')+
    ylab('projection strengths')
  
  print(df[duplicated(df[,c(2,3)])==FALSE, c(2,3,6)])
  
}

projection_vic <- function(Y,
                                mouse_label,
                                Z,
                                region_name,
                                motifs = NULL,
                           vic,ncol=3, limts = NULL){
  
  
  M <- length(Y)
  
  if(is.null(motifs)){
    motifs = sort(unique(unlist(Z)))
  }
  
  # Projection strength of each neuron
  Y_prop <- lapply(1:M,
                   function(m){
                     
                     matrix(Y[[m]]/rep(colSums(Y[[m]]), each = nrow(Y[[m]])),
                            
                            nrow = nrow(Y[[m]]),
                            ncol = ncol(Y[[m]]))
                   })
  
  # Data frame
  list0 <- lapply(1:M,
                  function(m){
                    
                    data.frame(projection_strength = as.vector(Y_prop[[m]]),
                               region_name = rep(region_name, ncol(Y[[m]])),
                               cluster = rep(Z[[m]], each = nrow(Y[[m]])),
                               #cluster = part_motif_names$pp.regions[part_motif_names$cluster==1]
                               mouse_label = rep(mouse_label[[m]], each = nrow(Y[[m]])),
                               vic = rep(vic[[m]], each = nrow(Y[[m]])))
                  })
  
  df <- do.call(rbind, list0)
  
  # Change neuron into factor
  df$neuron <- factor(rep(1:length(unlist(Z)), each = nrow(Y[[1]])),
                      levels = rep(1:length(unlist(Z))))
  
  # Change cluster into factors
  df$cluster <- factor(paste('cluster', df$cluster),
                       levels = paste('cluster', 1:max(unlist(Z))))
  
  # Change region names into factors
  df$region_name <- factor(df$region_name,
                           levels = region_name)
  
  df$projection_strength_avg = apply(df,1, function(x){mean(df[df$region_name==x[2]&df$cluster==x[3],1])})
  
  # Filter to clusters with at least 20 neurons 
  filt = apply(matrix(df$cluster,ncol=1),1, function(x){sum(x== paste('cluster',motifs))>0})
  df = df[filt==1,]
  df$cluster =  factor(df$cluster, levels = unique(df$cluster))
  
  # Line graph
  colors <- rev(sequential_hcl(5, palette = "Purple-Yellow")[1:4])
  ggplot(df)+
    geom_line(mapping = aes(x = region_name,
                            y = projection_strength,
                            colour = vic,
                            group = interaction(neuron, mouse_label)),alpha=0.75)+
    #scale_color_gradient2(mid='black',high='gray')+
    scale_color_gradientn(colours = colors, transform = "sqrt", labels = function(x) sprintf("%.4f", x), limits = limts)+
    #geom_line(mapping = aes(x = region_name,
    #                        y = projection_strength_avg,
    #                        group = cluster),
    #          col="red")+
    facet_wrap(~cluster,ncol=ncol)+
    theme_bw()+
    theme(axis.text.x = element_text(size = 14,angle = 90, vjust = 0.5, hjust = 1))+
    xlab('region')+
    ylab('projection strengths') +
    labs(color='VIC')
  
}

projection_filter_vic <- function(Y,
                              mouse_label,
                              Z,
                              region_name, 
                              vic){
  
  
  M <- length(Y)
  
  # Projection strength of each neuron
  Y_prop <- lapply(1:M,
                   function(m){
                     
                     matrix(Y[[m]]/rep(colSums(Y[[m]]), each = nrow(Y[[m]])),
                            
                            nrow = nrow(Y[[m]]),
                            ncol = ncol(Y[[m]]))
                   })
  
  # Data frame
  list0 <- lapply(1:M,
                  function(m){
                    
                    data.frame(projection_strength = as.vector(Y_prop[[m]]),
                               region_name = rep(region_name, ncol(Y[[m]])),
                               cluster = rep(Z[[m]], each = nrow(Y[[m]])),
                               mouse_label = rep(mouse_label[[m]], each = nrow(Y[[m]])),
                               vic = rep(vic[[m]], each = nrow(Y[[m]])))
                  })
  
  df <- do.call(rbind, list0)
  
  # Change neuron into factor
  df$neuron <- factor(rep(1:length(unlist(Z)), each = nrow(Y[[1]])),
                      levels = rep(1:length(unlist(Z))))
  
  # Change cluster into factors
  df$cluster <- factor(paste('cluster', df$cluster),
                       levels = paste('cluster', 1:max(unlist(Z))))
  
  # Change region names into factors
  df$region_name <- factor(df$region_name,
                           levels = region_name)
  
  # Filter to clusters with at least 20 neurons 
  filt = apply(matrix(df$cluster,ncol=1),1, function(x){sum(x== names(table(df$cluster))[table(df$cluster)>10*6])})
  df = df[filt==1,]
  df$cluster =  factor(df$cluster, levels = unique(df$cluster))
  
  df$projection_strength_avg = apply(df,1, function(x){mean(df[df$region_name==x[2]&df$cluster==x[3],1])})
  
  # Line graph
  ggplot(df)+
    geom_line(mapping = aes(x = region_name,
                            y = projection_strength,
                            col = vic,
                            group = interaction(neuron, mouse_label)))+
    scale_color_gradientn(colours = colors, transform = "sqrt", labels = function(x) sprintf("%.4f", x))+
    scale_color_gradient2(mid='black',high='gray')+
    geom_line(mapping = aes(x = region_name,
                            y = projection_strength_avg,
                            group = cluster),
              col="red")+
    facet_wrap(~cluster)+
    theme_bw()+
    xlab('region')+
    ylab('projection strengths')
  
}


opt.clustering.frequency <- function(clustering,
                                     main = ''){
  
  cluster.names <- paste(1:length(unique(unlist(clustering))))
  
  loop.result <- lapply(1:length(clustering), function(m){
    
    data.frame(cluster = factor(cluster.names,
                                levels = cluster.names),
               mouse = paste(m),
               counts = as.vector(table(factor(clustering[[m]],
                                               levels = 1:length(cluster.names)))))
  })
  
  # Frequency table
  z.frequency <- do.call(rbind, loop.result)
  
  # Bar chart
  z.frequency %>%
    ggplot(mapping = aes(x = cluster, y = counts, fill = mouse))+
    geom_bar(stat="identity")+
    theme_bw()+
    ylab('Number of neurons')+
    labs(fill="brain") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle(main)
}



opt.clustering.frequency.filter <- function(clustering,
                                     main = ''){
  
  cluster.names <- paste(1:length(unique(unlist(clustering))))
  
  loop.result <- lapply(1:length(clustering), function(m){
    
    data.frame(cluster = factor(cluster.names,
                                levels = cluster.names),
               mouse = paste(m),
               counts = as.vector(table(factor(clustering[[m]],
                                               levels = 1:length(cluster.names)))))
  })
  
  # Frequency table
  z.frequency <- do.call(rbind, loop.result)
  filt = apply(z.frequency,1, function(x){sum(z.frequency$counts[z.frequency$cluster==x[1]])>10})
  
  z.frequency = z.frequency[filt==1,]
  
  # Bar chart
  z.frequency %>%
    ggplot(mapping = aes(x = cluster, y = counts, fill = mouse))+
    geom_bar(stat="identity")+
    theme_bw()+
    ylab('Number of neurons')+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle(main)
}


projection_by_evi <- function(Y,
                                evi,
                                Z,
                                region_name){
  
  
  M <- length(Y)
  
  # Projection strength of each neuron
  Y_prop <- lapply(1:M,
                   function(m){
                     
                     matrix(Y[[m]]/rep(colSums(Y[[m]]), each = nrow(Y[[m]])),
                            
                            nrow = nrow(Y[[m]]),
                            ncol = ncol(Y[[m]]))
                   })
  
  # Data frame
  list0 <- lapply(1:M,
                  function(m){
                    
                    data.frame(projection_strength = as.vector(Y_prop[[m]]),
                               region_name = rep(region_name, ncol(Y[[m]])),
                               cluster = rep(Z[[m]], each = nrow(Y[[m]])),
                               evi_label = rep(evi[[m]], each = nrow(Y[[m]])))
                  })
  
  df <- do.call(rbind, list0)
  
  # Change neuron into factor
  df$neuron <- factor(rep(1:length(unlist(Z)), each = nrow(Y[[1]])),
                      levels = rep(1:length(unlist(Z))))
  
  # Change cluster into factors
  df$cluster <- factor(paste('cluster', df$cluster),
                       levels = paste('cluster', 1:max(unlist(Z))))
  
  # Change region names into factors
  df$region_name <- factor(df$region_name,
                           levels = region_name)
  
  # Line graph
  ggplot(df)+
    geom_line(mapping = aes(x = region_name,
                            y = projection_strength,
                            colour = evi_label,
                            group = interaction(neuron, evi_label)))+
    facet_wrap(~cluster)+
    scale_color_gradient2(mid='blue',high='red') +
    theme_bw()+
    xlab('region')+
    ylab('projection strengths') +
    labs(colour = "EVI")
  
}

projection_scatter_vic <- function(Y,
                              vi,
                              Z1,
                              Z2,
                              region_name,
                              region1,
                              region2){
  
  
  M <- length(Y)
  
  # Projection strength of each neuron
  Y_prop <- lapply(1:M,
                   function(m){
                     
                     matrix(Y[[m]]/rep(colSums(Y[[m]]), each = nrow(Y[[m]])),
                            
                            nrow = nrow(Y[[m]]),
                            ncol = ncol(Y[[m]]))
                   })
  
  # Data frame
  list0 <- lapply(1:M,
                  function(m){
                    
                    data.frame(region1 = Y_prop[[m]][region1,],
                               region2 = Y_prop[[m]][region2,],
                               cluster1 = Z1[[m]],
                               cluster2 = Z2[[m]],
                               vi = vi[[m]])
                  })
  
  df <- do.call(rbind, list0)
  
  df$projection_strength_avg_region1 = apply(df,1, function(x){mean(df[df$cluster1==x[3],1])})
  df$projection_strength_avg_region2 = apply(df,1, function(x){mean(df[df$cluster1==x[3],2])})
  
  # Line graph
  colors <- rev(sequential_hcl(5, palette = "Purple-Yellow")[1:4])
  ggplot(df)+
    geom_point(mapping = aes(x = region1,
                            y = region2,
                            colour = vi))+
    #scale_color_gradient2(mid='blue',high='red') +
    scale_color_gradientn(colours = colors, transform = "sqrt", labels = function(x) sprintf("%.4f", x))+
    theme_bw()+
    xlab(region_name[region1])+
    ylab(region_name[region2]) +
    labs(colour = "VIC")
  
}

projection_scatter <- function(Y,
                                   Z,
                                   region_name,
                                   region1,
                                   region2){
  
  
  M <- length(Y)
  
  # Projection strength of each neuron
  Y_prop <- lapply(1:M,
                   function(m){
                     
                     matrix(Y[[m]]/rep(colSums(Y[[m]]), each = nrow(Y[[m]])),
                            
                            nrow = nrow(Y[[m]]),
                            ncol = ncol(Y[[m]]))
                   })
  
  # Data frame
  list0 <- lapply(1:M,
                  function(m){
                    
                    data.frame(region1 = Y_prop[[m]][region1,],
                               region2 = Y_prop[[m]][region2,],
                               cluster = Z[[m]])
                  })
  
  df <- do.call(rbind, list0)
  
  df$projection_strength_avg_region1 = apply(df,1, function(x){mean(df[df$cluster==x[3],1])})
  df$projection_strength_avg_region2 = apply(df,1, function(x){mean(df[df$cluster==x[3],2])})
  
  df$cluster = as.factor(df$cluster)
  
  # Line graph
  ggplot(df)+
    geom_point(mapping = aes(x = region1,
                             y = region2,
                             colour = cluster))+
    geom_point(mapping = aes(x = projection_strength_avg_region1,
                             y = projection_strength_avg_region2,
                             colour = cluster), shape = 8, size=5) +
    theme_bw()+
    xlab(region_name[region1])+
    ylab(region_name[region2]) +
    labs(colour = "cluster")
  
}

projection_filter_vic2 <- function(Y,
                                  mouse_label,
                                  Z,
                                  region_name, 
                                  vic,
                                  cluster_inds){
  
  
  M <- length(Y)
  
  # Projection strength of each neuron
  Y_prop <- lapply(1:M,
                   function(m){
                     
                     matrix(Y[[m]]/rep(colSums(Y[[m]]), each = nrow(Y[[m]])),
                            
                            nrow = nrow(Y[[m]]),
                            ncol = ncol(Y[[m]]))
                   })
  
  # Data frame
  list0 <- lapply(1:M,
                  function(m){
                    
                    data.frame(projection_strength = as.vector(Y_prop[[m]]),
                               region_name = rep(region_name, ncol(Y[[m]])),
                               cluster = rep(Z[[m]], each = nrow(Y[[m]])),
                               mouse_label = rep(mouse_label[[m]], each = nrow(Y[[m]])),
                               EVI = rep(vic[[m]], each = nrow(Y[[m]])))
                  })
  
  df <- do.call(rbind, list0)
  
  # Change neuron into factor
  df$neuron <- factor(rep(1:length(unlist(Z)), each = nrow(Y[[1]])),
                      levels = rep(1:length(unlist(Z))))
  
  # Change cluster into factors
  df$cluster <- factor(paste('cluster', df$cluster),
                       levels = paste('cluster', 1:max(unlist(Z))))
  
  # Change region names into factors
  df$region_name <- factor(df$region_name,
                           levels = region_name)
  
  # Filter to clusters with at least 20 neurons 
  filt = apply(matrix(df$cluster,ncol=1),1, function(x){sum(x==paste('cluster', cluster_inds))})
  df = df[filt==1,]
  df$cluster =  factor(df$cluster, levels = unique(df$cluster))
  
  df$projection_strength_avg = apply(df,1, function(x){mean(df[df$region_name==x[2]&df$cluster==x[3],1])})
  
  # Line graph
  ggplot(df)+
    geom_line(mapping = aes(x = region_name,
                            y = projection_strength,
                            col = EVI,
                            group = interaction(neuron, mouse_label)))+
    scale_color_gradient(low='blue',high='red')+
#   geom_line(mapping = aes(x = region_name,
#                            y = projection_strength_avg,
#                            group = cluster),
#              col="red")+
    facet_wrap(~cluster)+
    theme_bw()+
    xlab('region')+
    ylab('projection strengths')
  
}

heatmap_VIC <- function(Y,
                       Z,
                       regions.name = NULL,
                       vic = NULL,
                       cluster.index = NULL,
                       col = NULL,
                       limts = NULL, 
                       title = ''){
  
  
  J <- max(unlist(Z))
  R <- nrow(Y[[1]])
  M <- length(Y)
  C <- sapply(1:M, function(m) ncol(Y[[m]]))
  C_cumsum <- c(0, cumsum(C))
  
  if(is.null(regions.name)){
    regions.name <- paste('region', 1:R)
  }
  if(is.null(vic)){
    vic <- rep(0,sum(unlist(C)))
  }
  if(is.null(cluster.index)){
    cluster.index <- 1:J
  }
  
  # turn vic into a list for each dataset
  vic.list <- lapply(1:M, function(m) vic[(1+C_cumsum[m]):C_cumsum[m+1]])
  
  Y.cbind <- do.call(cbind, Y)
  Y.prop.cbind <- matrix(as.vector(Y.cbind)/rep(colSums(Y.cbind), each = R),
                         nrow = R)
  
  Y.prop <- lapply(1:M, function(m) matrix(Y[[m]]/rep(colSums(Y[[m]]), each = R),
                                           nrow = R))
  
  
  # Re-ordered projection strength
  df0 <- lapply(cluster.index,
                function(j){
                  
                  
                  mx0 <- matrix(Y.prop.cbind[,which(unlist(Z)==j)],
                                nrow = R)
                  
                  strong.proj.j <- which(rowSums(mx0) == max(rowSums(mx0)))[1]
                  
                  df0_list <- lapply(1:M, function(m){
                    
                    if(length(which(Z[[m]]==j)) != 0){
                      
                      mx0_m <- matrix(Y.prop[[m]][,which(Z[[m]]==j)],
                                      nrow = R)
                      
                      mx0_m[,order(mx0_m[strong.proj.j,])]
                    }
                    
                  })
                  
                  do.call(cbind, df0_list)
                  
                })
  
  # group index
  vic.df <- lapply(cluster.index,
                           function(j){
                             
                             df0_list <- lapply(1:M, function(m){
                               
                               
                               if(length(which(Z[[m]]==j)) != 0){
                                 
                                 vic.list[[m]][which(Z[[m]]==j)]
                                 # rep(m, length(which(Z[[m]]==j)))
                               }
                               
                             })
                             
                             unlist(df0_list)
                             
                           })
  
  vic.df <- unlist(vic.df)
  
  
  df0 <- do.call(cbind, df0)
  
  N.j <- sapply(cluster.index,
                function(j) length(which(unlist(Z)==j)))
  
  # Convert to data frame for plotting
  # subset for cluster.index
  index <- unlist(lapply(1:M,function(m) {
    sapply(Z[[m]],function(x) x %in% cluster.index)
  }))
  Y.prop.df <- data.frame(region = rep(regions.name,
                                       ncol(Y.prop.cbind[,index])),
                          
                          neuron = rep(1:ncol(Y.prop.cbind[,index]),
                                       each = length(regions.name)),
                          
                          projection.strength = as.vector(df0),
                          
                          VIC = rep(vic.df, each = R))
  
  colors <- rev(sequential_hcl(5, palette = "Purple-Yellow")[1:4])
  if(!is.null(col)){
    colors = col
  }
  
  gg <- Y.prop.df %>%
    ggplot(mapping = aes(x = factor(region, levels = regions.name),
                         y = neuron,
                         fill = VIC))+
    
    geom_tile(mapping = aes(alpha = projection.strength))+
    scale_alpha_identity()+
    theme_bw()+
    xlab('region')+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          plot.title = element_text(size=12),
          plot.background = element_blank() ,
          panel.grid.major = element_blank() ,
          panel.grid.minor = element_blank() ,
          panel.border = element_blank() ,
          panel.background = element_blank()) +
    #draws x and y axis line
    theme(axis.line = element_line(color = 'black'))+
    # rotate region names labels
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ylab('neurons')+
    
    geom_hline(yintercept = cumsum(N.j)[-length(cumsum(N.j))],
               color = 'grey',
               linetype = 'dashed',
               linewidth = 0.5)+
    
    scale_fill_gradientn(colours = colors, transform = "sqrt", labels = function(x) sprintf("%.4f", x), limits = limts) +
    
    ggtitle(title)
  
  print(gg)
}

