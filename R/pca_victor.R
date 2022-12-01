
#-----------------------------------------

# This function is built on pca output of vegan rda function
# metadata comprises the info on the samples and must be in the same order than the PCA coordinates
# goodness.axis : the number of axis used to calculate the goodness of fit of the species
# main.group define the grouping variable used to calculate the densities and to define the colors, it cannot be empty
# second.group is used to define a second level of goruping that, if present, will be highligthed with different shapes
# nudge.x and nudge.y define the gap between the labels and the group centroid, it must have the same length as the grouping variable defined in main.group
# labels are use to labels the main groups while different color are given to the labels according to the second grouping variable if there is one
# scale.fill : colour pallette for fill aesthetics
# scale.colour : colour pallette for colour aesthetics
# scale.shape : vector of shape types for geom_point aesthetics
# print.sp.thresh : print the goodness of fit threshold for selecting species' vectors
# ext.plot.scale : scale for density plots
# add.stat1 and 2: logical. add a convex hull or a 95% CI ellipse in each of the plots
# stat1 and 2: either "chull" for a convex hull or "ellipse" for an ellipse
# add dens : logical, either to draw density plots or not
# add.centroids :  logical, either to draw centroids of each main group or not
# c.size : centroid point size
# p.alpha : point oppacity 
# axis.size : font size for axes titles
# conf.level : confidence level for ellipses

pca_victor <- function(pca,axes=c(1,2), metadata= NULL, site.scaling = 1, species.scaling = 2, goodness.axis = 2, goodness.thresh = 0.3, main.group = NULL, second.group = NULL, nudge.x = NULL, nudge.y = NULL, scale.fill = NULL, scale.shape = NULL, scale.colour = NULL, labels = NULL, print.sp.thresh = FALSE, ext.plot.scale = 5, add.stat1 = TRUE, stat1 = "chull", add.stat2 = FALSE, stat2 = NULL, add.dens = TRUE, add.centroids = FALSE, point.size = 2.5, c.size = 5, p.alpha = 1, font.size = 18/.pt, axis.size = 22, axis.text = 20, conf.level = .95, ysites = NULL, xsites = NULL, ysp = NULL, xsp = NULL, add.labels = TRUE){

  require(ggplot2)
  require(dplyr)
  require(cowplot)
  require(data.table)
  require(ggrepel)
  require(gridExtra)

  metadata <- as.data.frame(metadata)

  if(length(axes) > 2){
    stop("This function can only plot two axes simultenaously, please select only 2 axes")
  }
  if(!is.null(labels)){
    if(length(labels) != length(unique(metadata[,main.group]))){
      stop("Number of labels provided do not match with the number of groups provided in the `main.group` variable")
    }
  }

#-----------------------------

  # Retrieve sites' scores (SCALING 1)
  site_scores <- scores(pca,scaling=site.scaling,display="sites",choices=axes)

  axes_name <- colnames(site_scores)

  if(nrow(metadata) != nrow(site_scores)){
    stop("There are no metadata available or the metadata do not match the number of sites in the pca output. This function needs metadata to customize the pca plot")
  }

  # Add the metadata to sites' scores (Scaling 1)
  site_scores <- cbind(metadata,site_scores)
  
  # Retrieve the site's scores in scaling 2 for the left plot and bind it to the metadata
  site_scores2 <- cbind(metadata, scores(pca, scaling= species.scaling, display="sites",choices=axes))

  # Retrieve species' scores (Scaling 2)
  sp_scores <- scores(pca,scaling = species.scaling, display="species", choices=axes)

  # Retrieve the goodness of fit of the species
  sp_fit <- goodness(pca,model="CA",statistic="explained")

  # Assemble species' scores and the goodness of fit for the selected axis
  sp_scores <- data.frame(species=rownames(sp_scores),sp_scores, fit=sp_fit[,goodness.axis])

  # Calculate the variance represented by each axis
  var_axes <- round(pca$CA$eig/sum(pca$CA$eig)*100,2)

  if(is.null(main.group)){
    stop("main.group cannot be empty as it defines the group on which densities and colours are defined")
  }

  # Density of sites along the first axis selected in "axes"
  PC1_dens <- site_scores %>%
    group_by_(main.group) %>%
    do(ggplot2:::compute_density(select_(.,axes_name[1]) %>% pull(), NULL)) %>%
    setnames("x", "PC1")

  # Density along the second axis selected in "axes"
  PC2_dens <- site_scores %>%
    group_by_(main.group) %>%
    do(ggplot2:::compute_density(select_(.,axes_name[2]) %>% pull(), NULL)) %>%
    setnames("x", "PC2")

  # Upper limit of the density curves
  dens_limit <- max(PC1_dens$density, PC2_dens$density) * 1.2

  # Define labels position
  if(is.null(nudge.x) | is.null(nudge.y)){
    # If there is no second.group
    if(is.null(second.group)){
      label_pos <- site_scores %>%
        group_by_(main.group) %>%
        summarise_at(vars(contains("PC")),mean)  %>%
        ungroup() %>%
        mutate(nudge_x = rep(0,nrow(.)),
               nudge_y = rep(0,nrow(.)))
    }else{ # If there is a second group
      label_pos <- site_scores %>%
        group_by_(main.group,second.group) %>%
        summarise_at(vars(contains("PC")),mean)  %>%
        ungroup() %>%
        mutate(nudge_x = rep(0,nrow(.)),
               nudge_y = rep(0,nrow(.)))
    }
  }else{
    # If there is no second.group
    if(is.null(second.group)){
      label_pos <- site_scores %>%
        group_by_(main.group) %>%
        summarise_at(vars(contains("PC")),mean)  %>%
        ungroup() %>%
        mutate(nudge_x = nudge.x,
               nudge_y = nudge.y)
    }else{ # If there is a second group
      label_pos <- site_scores %>%
        group_by_(main.group,second.group) %>%
        summarise_at(vars(contains("PC")),mean)  %>%
        ungroup() %>%
        mutate(nudge_x = nudge.x,
               nudge_y = nudge.y)
    }
  }

# Renaming the two axis to make their selection easier for labelling
colnames(label_pos)[which(colnames(label_pos)%in%axes_name)] <- c("PC1","PC2")

# Calculate the coordinates of the labels
label_pos <- label_pos %>%
  mutate(x = PC1 + nudge_x, y = PC2 + nudge_y)

#-----------------------------

  # Main plot of sites
  #-------------------
  p_sites <- ggplot()
  if(is.null(second.group)){
    p_sites <- p_sites + geom_point(data=site_scores, aes_string(x = axes_name[1], y = axes_name[2], fill = main.group, colour = main.group), alpha= p.alpha, size = point.size)
    }else{
    p_sites <- p_sites + geom_point(data=site_scores, aes_string(x = axes_name[1], y = axes_name[2], fill = main.group, shape=second.group), alpha= p.alpha, size = point.size)
    }
  # Add a convex hull or ellipses around the main group
  if(add.stat1 == TRUE){
    if(!is.null(stat1)){
      if(stat1 == "chull"){
        p_sites <- p_sites + stat_chull(data = site_scores, aes_string(x = axes_name[1], y = axes_name[2], fill = main.group, col = main.group), geom="polygon", alpha=0.3, linetype = 2)+
      scale_colour_manual(values = scale.colour)
      }else{
        p_sites <- p_sites + stat_ellipse(data = site_scores, aes_string(x = axes_name[1], y = axes_name[2], fill = main.group, col = main.group), geom="polygon", alpha=0.3, linetype = 2, level = conf.level)+
          scale_colour_manual(values = scale.colour)
      }
    }
  }
  # Add centroids
  if(add.centroids == TRUE){
    centroids <-  site_scores %>%
      group_by_(main.group) %>%
      summarise_at(vars(contains("PC")),mean)  %>%
      ungroup()
    p_sites <- p_sites + geom_point(data = centroids, aes_string(x = axes_name[1], y = axes_name[2], fill = main.group), shape = 24, col = "black", size = c.size, alpha = p.alpha)
  }

  # Custmize colour of fill if provided
  if(!is.null(scale.fill)){
  p_sites <- p_sites + scale_fill_manual(values=scale.fill)
  }
  # Customize shape if needed and values are provided
  if(!is.null(second.group) & !is.null(scale.shape)){
  p_sites <- p_sites + scale_shape_manual(values=scale.shape)
  }
  # Add variance of axes
  p_sites <- p_sites + theme(legend.position="none") + xlab(paste(axes_name[1],": ",var_axes[axes_name[1]],"%")) + ylab(paste(axes_name[2],": ",var_axes[axes_name[2]],"%"))
  
  # Add same scale
  p_sites <- p_sites + coord_fixed(ratio = 1, ylim = ysites, xlim = xsites)
  
  # Add labels of main group as provided in the variable if there is no labels provided to override this. Color of labels are unique if there is not a second grouping variable
  if(add.labels == TRUE){
  if(is.null(labels) & is.null(second.group)){
    p_sites <- p_sites + geom_label(data = label_pos, aes_string(x = "x", y = "y", label = main.group, fill = main.group), alpha= 1, fontface="bold", size = font.size, colour = "black", inherit.aes = FALSE)
  }
  # If there are labels provided and no second group
  if(!is.null(labels) & is.null(second.group)){
    p_sites <- p_sites + geom_label(data = label_pos, aes_string(x = "x", y = "y", label = "labels", fill = main.group), alpha= 1, fontface="bold", size = font.size, colour = "black", inherit.aes = FALSE)
  }
  # If there are no labels provided and a second grouping variable
  if(is.null(labels) & !is.null(second.group)){
    p_sites <- p_sites + geom_label(data = label_pos, aes_string(x = "x", y = "y", label = main.group, fill = main.group, colour = main.group), alpha= 1, fontface="bold", size = font.size, colour = "black", inherit.aes = FALSE)
  }
  # If there are labels provided and a second grouping variable
  if(!is.null(labels) & !is.null(second.group)){
    p_sites <- p_sites + geom_label(data = label_pos, aes_string(x = "x", y = "y", label = "labels", fill = main.group, colour = main.group), alpha= 1, fontface="bold", size = font.size, colour = "black", inherit.aes = FALSE)
  }
}
  p_sites <- p_sites + theme(axis.title = element_text(size = axis.size, face = "bold"), axis.text = element_text(size = axis.text))
  if(add.dens == TRUE){# Density plot of sites along the first axis selected
    #----------------------------------------------------
    x_dens <- axis_canvas(p_sites, axis = "x")
    x_dens <- x_dens + geom_density(data = PC1_dens, aes_string(x = "PC1", y = "density", fill = main.group), alpha = 0.4,stat = "identity")
    if(!is.null(scale.fill)){
    x_dens <- x_dens + scale_fill_manual(values=scale.fill)
    }
    x_dens <- x_dens + scale_y_continuous(limits = c(0, dens_limit), expand = c(0, 0)) + theme(legend.position="FALSE")
  
    # Density plot of sites along the second axis selected
    #-----------------------------------------------------
    y_dens <- axis_canvas(p_sites, axis = "y", coord_flip = TRUE)
    y_dens <- y_dens + geom_density(data = PC2_dens, aes_string(x = "PC2", y = "density", fill = main.group), alpha = 0.4,stat = "identity")
    if(!is.null(scale.fill)){
    y_dens <- y_dens + scale_fill_manual(values=scale.fill)
    }
    y_dens <- y_dens + scale_y_continuous(limits = c(0, dens_limit), expand = c(0, 0))
    y_dens <- y_dens + coord_flip() + theme(legend.position="FALSE")
  
    # Assembly the 3 plots for the site
    #----------------------------------
    p1 <- insert_xaxis_grob(p_sites +
                              theme(legend.position = "none"),
                            x_dens, grid::unit(ext.plot.scale*14, "pt"), position = "top")
  
    p2 <- insert_yaxis_grob(p1, y_dens, grid::unit(ext.plot.scale*14, "pt"), position = "right")
  
    p_sites <- ggdraw(p2)
  }
#-----------------------------

  # Main plot of species
  #-------------------

  p_sp <- ggplot()
  if(is.null(second.group)){
     p_sp <- p_sp + geom_point(data = site_scores2, aes_string(x = axes_name[1], y = axes_name[2], fill = main.group, colour = main.group), alpha= p.alpha, size = point.size)
  }else{
	  p_sp <- p_sp + geom_point(data=site_scores2, aes_string(x = axes_name[1], y = axes_name[2], fill = main.group, shape=second.group), alpha= p.alpha)
  }
  #Add a convex hull or ellipses
  if(add.stat2 == TRUE){
    if(!is.null(stat2)){
      if(stat2 == "chull"){
        p_sp <- p_sp + stat_chull(data = site_scores2, aes_string(x = axes_name[1], y = axes_name[2], fill = main.group, col = main.group), geom="polygon", alpha=0.3, linetype = 2)+
          scale_colour_manual(values = scale.colour)
      }else{
        p_sp <- p_sp + stat_ellipse(data = site_scores2, aes_string(x = axes_name[1], y = axes_name[2], fill = main.group, col = main.group), geom="polygon", alpha=0.3, linetype = 2, level = conf.level)+
          scale_colour_manual(values = scale.colour)
      }
    }
  }

  # Customize shape if needed and values are provided
  if(!is.null(second.group) & !is.null(scale.shape)){
  	p_sp <- p_sp + scale_shape_manual(values=scale.shape)
  }else{
  	p_sp <- p_sp + scale_shape_manual(values=21)
  }
  if(!is.null(scale.fill)){
  p_sp <- p_sp + scale_fill_manual(values=scale.fill)
  p_sp <- p_sp + scale_colour_manual(values=scale.fill)
  }
  p_sp <- p_sp + geom_segment(data=sp_scores[sp_scores$fit > goodness.thresh,], aes_string(x=0, y=0, xend=axes_name[1], yend=axes_name[2]), colour="black")
  p_sp <- p_sp + geom_point(data=sp_scores[sp_scores$fit > goodness.thresh,], aes_string(x = axes_name[1], y = axes_name[2]), fill="black",shape=21,size= point.size)
  p_sp <- p_sp + geom_text_repel(data=sp_scores[sp_scores$fit > goodness.thresh,], aes_string(x = axes_name[1], y = axes_name[2], label = "species"), fontface="bold", colour="black", segment.colour="black", max.overlaps = Inf, size = font.size)
  p_sp <- p_sp + theme(legend.position="none") + xlab(paste(axes_name[1],": ",var_axes[axes_name[1]],"%")) + ylab(paste(axes_name[2],": ",var_axes[axes_name[2]],"%"))
  
  p_sp <- p_sp + theme(axis.title = element_text(size = axis.size, face = "bold"), axis.text = element_text(size = axis.text))
  p_sp <- p_sp + coord_fixed(ratio = 1, ylim = ysp, xlim = xsp)
  
  
  # Trick to scale the species plot as the site plot
  if(add.dens == TRUE){
  blank_plot <- ggplot()+
    theme_void()

    if(print.sp.thresh){
      p1 <- insert_xaxis_grob(p_sp +
      theme(legend.position = "none"), blank_plot + annotate(geom="text",label=paste("Threshold for modalities representation :", goodness.thresh * 100, "% on the first",goodness.axis,"axes"), x = 0, y =0,family="Comic Sans MS"), grid::unit(ext.plot.scale*14, "pt"), position = "top")
    }else{
      p1 <- insert_xaxis_grob(p_sp +
                                theme(legend.position = "none"),
                              blank_plot, grid::unit(ext.plot.scale*14, "pt"), position = "top")
    }
  

  p2 <- insert_yaxis_grob(p1, blank_plot, grid::unit(ext.plot.scale*14, "pt"), position = "right")

  p_sp <- ggdraw(p2)
  }
#-----------------------------

  # Final plots
  #------------
  grid.arrange(p_sites,p_sp, ncol=2, nrow = 1, widths =  c(.5,.5))
}


