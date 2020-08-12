
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# perform meta-analysis given fits from limma
run_meta_analysis = function(tabList, method="FE") {

  # extract logFC and se for each model fit
  resList = lapply( tabList, function(tab){
    tab$se = tab$logFC /  tab$t
    tab
  })

  # extract features names
  featureNames = lapply(resList, function(res){
    rownames(res)
    })
  featureNames = unique(unlist(featureNames))

  # check that each fit has all features
  # testIdentical = lapply(resList, function(res){
  #   identical(sort(rownames(res)), sort(featureNames))
  #   })
  # if(  any(!unlist(testIdentical)) ){
  #   stop("each entry in fitList must have all features with the same names")
  # }

  # apply meta analysis for each feature
  resRMA = lapply(featureNames, function(feature){
    # extract summary statistics for this feature from all studies
    res = lapply(resList, function(res){
      res[feature,]
      })
    df = do.call("rbind", res)
    suppressWarnings(rma(yi=logFC, sei=se, data = df, method=method))
    })
  names(resRMA) = featureNames

  # extract results
  # Currently exatract simple results
  # other results may be relevant if meta-analysing more than two studies
  resTable = lapply(resRMA, function(x){
    data.frame( logFC     = x$beta,
                se        = x$se,
                P.Value   = x$pval,
                adj.P.Val = p.adjust(x$pval, "fdr"),
                Q         = x$QE,
                N         = x$k,
                I2        = x$I2)
    })
  do.call("rbind", resTable)
}

plotVolcano = function(df, nshow=10){
  # indicate significant genes
  df$isSignif = c("no","yes")[(df$adj.P.Val < 0.05)+1]

  # show gene names
  df$EnsID = gsub('\\..*$', '', rownames(df))

  df_hugo = select(org.Hs.eg.db, df$EnsID, "SYMBOL", "ENSEMBL")
  df = merge(df, df_hugo, by.x='EnsID', by.y="ENSEMBL")

  # df = df[!duplicated(df$gene_id),]

  # top significant genes
  df2 = df[df$adj.P.Val <= sort(df$adj.P.Val )[nshow], ]

  ggplot(df, aes( logFC, -log10(P.Value), color=isSignif)) + geom_point() + theme_bw(15) + theme(aspect.ratio=1, legend.position="none", plot.title = element_text(hjust = 0.5)) + xlab(bquote(log[2]~fold~change)) + ylab(bquote(-log[10]~P)) + scale_color_manual(values=c("grey", "red")) + geom_text_repel(data=df2, aes(logFC, -log10(P.Value), label=SYMBOL), segment.size=.5,  segment.color="black", color="black", force=1, nudge_x=.005, nudge_y=.5)
}


plot_forrest = function(resList, ensID) { 
  # create a single data.frame with logFC and se
  df_forrest = lapply(resList, function(x) if(ensID %in% x[ensID,]) {
    x[ensID,]
  } 
  )
  df_forrest[sapply(df_forrest, is.null)] <- NULL
  df_forrest = do.call(rbind, df_forrest)
  df_forrest$Dataset = rownames(df_forrest)
  df_forrest$se = df_forrest$logFC /  df_forrest$t

  # bind only shared columns
  df = gtools::smartbind(df_forrest,  res_dlpfc[ensID,])
  df$Dataset = factor(df$Dataset, rev(df$Dataset))

  # set range to include zero
  rng = range(c(with(df, logFC-se), with(df, logFC+se)))
  if( rng[1] > 0 ){
    rng[1] = 0
  }
  if( rng[2] < 0 ){
    rng[2] = 0
  }

  # forrest plot
  fig1 = ggplot(df, aes(Dataset, logFC)) + geom_point(color="navy") + geom_errorbar(aes(ymin = logFC - se, ymax=logFC + se), width=0.1) + theme_bw(12) + theme(aspect.ratio= .2, plot.title = element_text(hjust = 0.5)) + coord_flip() + ylim(rng) + xlab(bquote(log[2]~fold~change)) + ggtitle('RCSD1')

  # Adjusted p-value
  ymax = max(-log10(df$adj.P.Val))*1.05
  fig2 = ggplot(df, aes(Dataset, -log10(adj.P.Val))) + geom_bar(stat="identity", fill="navy") + theme_bw(12) + theme( aspect.ratio= .5, plot.margin = unit(c(0,0,0,0), "cm"),
    axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) + coord_flip() + ylab(bquote(-log[10]~adjusted~P)) + scale_y_continuous(expand=c(0,0), limits=c(0, ymax))

  # plot_grid(fig1, fig2, align='vh')
  fig1 + fig2
}
