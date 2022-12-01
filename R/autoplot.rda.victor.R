##' @title ggplot-based plot for objects of class \code{'rda'}
##'
##' @description
##' Produces a multi-layer ggplot object representing the output of objects produced by \code{\link[vegan]{rda}}.
##'
##' @details
##' TODO
##'
##' @param object an object of class \code{"rda"}, the result of a call to \code{\link[vegan]{rda}}
##' @param axes numeric; which axes to plot, given as a vector of length 2.
##' @param geom character; which geoms to use for the layers. Can be a vector of
##'   up to length 2, in which case, the first element of \code{geom} will be
##'   used for any site scores (both weighted sum or linear combination scores),
##'   and the second element will be used for species scores. The latter will be
##'   ignored if \code{arrows = TRUE}.
##' @param layers character; which scores to plot as layers
##' @param arrows logical; represent species (variables) using vectors?
##' @param legend.position character or two-element numeric vector; where to position the legend. See \code{\link[ggplot2]{theme}} for details. Use \code{"none"} to not draw the legend.
##' @param xlab character; label for the x-axis
##' @param ylab character; label for the y-axis
##' @param title character; subtitle for the plot
##' @param subtitle character; subtitle for the plot
##' @param caption character; caption for the plot
##' @param const General scaling constant to \code{rda} scores. See
##'   \code{\link[vegan]{scores.rda}} for details.
##' @param stat character: type of geom "chull" (stat_chull) for convex hull and "ellipse" for confidence ellipse (stat_ellipse). Default is NULL.
##' @param thresh numeric: model fit threshold for selecting species (variables)
##' @param lvl numeric: confidence level for ellipse
##' @param title.size numeric: size of axes and legend titles
##' @param font.size character: font size using the ggrepel syntax
##' @param scale.colour character: hex codes for colour aesthetics
##' @param scale.fill character: hex codes for fill aesthetics
##' @param ... Additional arguments passed to \code{\link{fortify.cca}}.
##'
##' @return Returns a ggplot object.
##'
##' @author Gavin L. Simpson (changed by Victor L. Jardim)
##'
##' @export
##'
##' @importFrom grid arrow unit
##' @importFrom ggplot2 autoplot ggplot geom_point geom_text geom_segment labs coord_fixed aes_string
##'
##' @examples
##'
##' data(dune)
##'
##' pca <- rda(dune)
##' autoplot(pca)
##'
##' ## Just the species scores
##' autoplot(pca, layers = "species")
`autoplot.rda.victor` <- function(object, axes = c(1,2), geom = c("point", "text"),
                           layers = c("species", "sites", "biplot", "centroids"),
                           arrows = TRUE, legend.position = "right",
                           title = NULL, subtitle = NULL, caption = NULL, metadata = NULL, scale.fill = NULL, scale.colour = NULL, stat = NULL, thresh = 0.25, lvl = .95, ylab, xlab, const,title.size = 16, font.size = 18/.pt, ylim = NULL, xlim = NULL, ...) {
    ## determine which layers to plot
    valid <- valid_layers(object)       # vector of valid layers
    ok_layers <- check_user_layers(layers, valid, message = TRUE)
    layers <- layers[ok_layers]         # subset user-supplied layers
    draw_list <- layer_draw_list(valid, layers) # what are we drawing
    
    ##select species to plot
    sp_fit <- goodness(object,model="CCA",statistic="explained", display = "species", choices = axes) 
    
    sp_fit <- as.data.frame(sp_fit) 
    colnames(sp_fit) <- c("RDA1.fit", "RDA2.fit")
    
    sp_fit <- as.data.frame(sp_fit) %>% 
        mutate(Label = rownames(.)) %>% 
        select(Label, everything())

    ## fix-up axes needed to plot
    laxes <- length(axes)
    if (laxes != 2L) {
        if (laxes > 2L) {
            axes <- rep(axes, length.out = 2L)  # shrink to required length
        } else {
            stop("Need 2 ordination axes to plot; only 1 was given.",
                 call. = FALSE)
        }
    }

    obj <- fortify(object, axes = axes, const = const, ...) # grab some scores
    available <- levels(obj[["Score"]])
    draw_list <- layer_draw_list(valid, layers, available) # what are we drawing
    layer_names <- names(draw_list)[draw_list]
    

    ## sort out x, y aesthetics
    vars <- getDimensionNames(obj)

    ## process geom arg
    geom <- match.arg(geom, several.ok = TRUE)
    geom <- unique(geom)    # simplify geom if elements are the same

    ## subset out the layers wanted
    obj <- obj[obj[["Score"]] %in% layer_names, , drop = FALSE]

    ## skeleton layer
    plt <- ggplot()

    ## draw sites, species, constraints == lc site scores
    if (any(draw_list[c("species","sites","constraints")])) {
        plt <- add_spp_site_scores(obj, plt, vars, geom, draw_list, arrows, metadata, scale.fill = scale.fill, scale.colour = scale.colour, thresh = thresh, sp_fit = sp_fit, stat = stat, lvl = lvl, font.size = font.size, title.size = title.size)
    }

    ## remove biplot arrows for centroids if present
    if(all(draw_list[c("biplot","centroids")])) {
        want <- obj[["Score"]] == "biplot"
        tmp <- obj[want, ]
        obj <- obj[!want, ]
        bnam <- tmp[, "Label"]
        cnam <- obj[obj[["Score"]] == "centroids", "Label"]
        obj <- rbind(obj, tmp[!bnam %in% cnam, , drop = FALSE])
    }

    if(isTRUE(draw_list["biplot"])) {
        want <- obj[["Score"]] == "biplot"
        if (length(layer_names) > 1) {
            mul <- arrowMul(obj[want, vars, drop = FALSE],
                            obj[!want, vars, drop = FALSE])
            obj[want, vars] <- mul * obj[want, vars]
        }
        col <- "gray35"
        plt <- plt +
            geom_segment(data = obj[want, , drop = FALSE ],
                         aes_string(x = 0, y = 0,
                                    xend = vars[1], yend = vars[2]),
                         arrow = arrow(length = unit(0.4, "cm")),
                         colour = col, size = .75)
        obj[want, vars] <- 1.1 * obj[want, vars]
        plt <- plt + geom_text_repel(data = obj[want, , drop = FALSE ],
                               aes_string(x = vars[1], y = vars[2],
                                          label = 'Label'), size = font.size, colour = col, point.size = 3, max.overlaps = Inf)
    }

    if(isTRUE(draw_list["centroids"])) {
        want <- obj[["Score"]] == "centroids"
        plt <- plt +
            geom_point(data = obj[want, , drop = FALSE],
                      aes_string(x = vars[1], y = vars[2], label = 'Label'), shape = 24, fill = scale.fill, size = 5)
    }
    var_axes <- as.numeric(round(object$CCA$eig/object$tot.chi*100,2))
    if(missing(xlab)) {
        xlab <- paste(vars[1],": ",var_axes[axes[1]],"%")
    }
    if(missing(ylab)) {
        ylab <- paste(vars[2],": ",var_axes[axes[2]],"%")
    }
    plt <- plt + labs(x = xlab, y = ylab, title = title, subtitle = subtitle,
                      caption = caption)
    ## add equal scaling
    plt <- plt + coord_fixed(ratio = 1, ylim = ylim, xlim = xlim)
    ## do we want a legend
    plt <- plt + theme(legend.position = legend.position, axis.title = element_text(size = title.size, face = "bold")) +
        guides(shape = "none")
    plt
}
