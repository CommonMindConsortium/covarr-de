
library(qvalue)
# library(EnsDb.Hsapiens.v86)
library(AnnotationHub)
library(ensembldb)
library(MASS)
library(ggplot2)
library(viridis)
library(data.table)

# For ENSEMBL id ENSG00000279457.4, return ENSG00000279457
trim_ensembl_ids = function(x){
  gsub("(.*)\\.(.*)", "\\1", x) 
}

# Add column Symbol using column gene_id storing ENSEMBL id
getGeneSymbol = function( df, column="row.names"){

  # load load since data base can fail with rmarkdown run multiple times
  # detach("package:EnsDb.Hsapiens.v86", unload=TRUE)
  # library(EnsDb.Hsapiens.v86)

  # Load ENSEMBL v96 database
  # see https://www.bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/ensembldb.html#10_getting_or_building_ensdb_databasespackages
  ah <- AnnotationHub()
  ensdb = ah[["AH69187"]] # ENSEMBL v96

  if( column == "row.names"){    
    geneID = rownames(df)
  }else{
    geneID = df[[column]]
  }

  ensIDs = trim_ensembl_ids( geneID )

  geneSymbol = mapIds(ensdb, keys=ensIDs, keytype="GENEID", column="GENENAME")
  df$Symbol = geneSymbol

  chroms = ensembldb::select(ensdb, keys=ensIDs, keytype="GENEID", column="SEQNAME")

  df2 = merge( data.frame(GENEID = ensIDs, stringsAsFactors=FALSE), chroms, by='GENEID')

  idx = match(ensIDs, chroms$GENEID)
  df$Chrom = c()
  df$Chrom = chroms$SEQNAME[idx]

  df
}

get_density <- function(x, y, n = 250) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

plot_concordance = function( resList, showGenes = NULL, idx=c(1,2), size=15){

  df = merge( resList[[idx[1]]], resList[[idx[2]]], by="row.names")

  lab = "log2 fold change"
  lim = with(df, max(abs(c(logFC.x, logFC.y))))
  df$density <- get_density(df$logFC.x, df$logFC.y, n = 100)
  fig1 = ggplot(df, aes(logFC.x, logFC.y, color=density)) + geom_point(size=.4) + theme_bw(size) + theme(aspect.ratio=1, legend.position="bottom") + geom_abline(color="red") + xlab(paste(names(resList)[1], lab)) + ylab(paste(names(resList)[2], lab)) + geom_vline(xintercept=0, col="grey40", linetype="dashed") + geom_hline(yintercept=0, col="grey40", linetype="dashed") + xlim(-lim, lim) + ylim(-lim, lim) + scale_color_viridis() + geom_smooth(method="rlm", se=FALSE, color="darkorange")

  if( !is.null(showGenes) ){
    df2 = df[df$Symbol.x %in% showGenes,]
    fig1 = fig1 + geom_text_repel(data=df2, aes(logFC.x, logFC.y, label=Symbol.x), segment.size=.5,  segment.color="black", color="black", force=1, nudge_x=.005, nudge_y=.1)
  }


  lab = "t-statistic"
  lim = with(df, max(abs(c(t.x, t.y))))
  df$density <- get_density(df$t.x, df$t.y, n = 100)
  fig2 = ggplot(df, aes(t.x, t.y, color=density)) + geom_point(size=.4) + theme_bw(size) + theme(aspect.ratio=1, legend.position="bottom") + geom_abline(color="red") + xlab(paste(names(resList)[1], lab)) + ylab(paste(names(resList)[2], lab)) + geom_vline(xintercept=0, col="grey40", linetype="dashed") + geom_hline(yintercept=0, col="grey40", linetype="dashed") + xlim(-lim, lim) + ylim(-lim, lim) + scale_color_viridis() + geom_smooth(method="lm", se=FALSE, color="darkorange")

  if( !is.null(showGenes) ){
    df2 = df[df$Symbol.x %in% showGenes,]
    fig2 = fig2 + geom_text_repel(data=df2, aes(t.x, t.y, label=Symbol.x), segment.size=.5,  segment.color="black", color="black", force=1, nudge_x=.005, nudge_y=.3)
  }

  res = with(df, cor.test(logFC.x, logFC.y, method="spearman"))
  
  p = with(df, P.Value.y[adj.P.Val.x < 0.05])
  pi_discovery_x = tryCatch( 
    1 - qvalue(p)$pi0, 
    error = function(e){
      NA
      })

  p = with(df, P.Value.x[adj.P.Val.y < 0.05])
  pi_discovery_y = tryCatch( 
    1 - qvalue(p)$pi0, 
    error = function(e){
      NA
      })

  tab_pi = data.frame(Discovery = names(resList), pi1 = c(pi_discovery_x, pi_discovery_y))

  tab_corr = data.frame(R_spearman = res$estimate, P = res$p.value)

  list(fig_logFC  = fig1,
      fig_tstat   = fig2, 
      tab_pi      = tab_pi,
      tab_corr    = tab_corr)
}



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
# if number of observation is less then minObs, return NA for that gene
run_meta_analysis = function(tabList, method="FE", minObs = 2){

  # extract logFC and se for each model fit
  resList = lapply( tabList, function(tab){
    tab$se = tab$logFC / tab$t
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
    df = df[!is.na(df$t),]

    # if gene is only observed one, meta-analysis should be NA
    if( nrow(df) >= minObs){
      tmp = suppressWarnings(rma(yi=logFC, sei=se, data = df, method=method))
    }else{
        tmp = NA
    }
    tmp
  })
  names(resRMA) = featureNames 
  keep = !sapply(resRMA, is.logical)
  resRMA = resRMA[keep]

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


plotVolcano = function(df, showGenes = NULL, size=15, minp=1e-310 ){

  # indicate significant genes
  df$isSignif = c("no","yes")[(df$adj.P.Val < 0.05)+1]
  df$P.Value = pmax(minp, df$P.Value )

  # get gene Symbol from rownames
  df = getGeneSymbol(df)

  lim = -log10(min(df$P.Value))
  limx = c(-max(abs(df$logFC)), max(abs(df$logFC)))
  fig = ggplot(df, aes( logFC, -log10(P.Value), color=isSignif)) + geom_point() + theme_bw(size) + theme(aspect.ratio=1, legend.position="none", plot.title = element_text(hjust = 0.5)) + xlab(bquote(log[2]~fold~change)) + ylab(bquote(-log[10]~P)) + scale_color_manual(values=c("grey", "red")) + scale_y_continuous(limits=c(0, lim*1.02), expand=c(0,0)) + xlim(limx)

  if( !is.null(showGenes) ){
    # top significant genes
    df2 = df[df$Symbol %in% showGenes,]

    fig = fig + geom_text_repel(data=df2, aes(logFC, -log10(P.Value), label=Symbol), segment.size=.5,  segment.color="black", color="black", force=1, nudge_x=.005, nudge_y=.5)
  }

  fig
}



plot_forrest = function(resList, geneSymbol){

  df_forrest = lapply(resList, function(x){
    idx = (x$Symbol == geneSymbol) & ! is.na(x$Symbol)
    if( length(idx) == 0){
      stop("Gene not found: ", geneSymbol)
    }
    x[idx,] 
  })

  df_forrest = do.call(gtools::smartbind, df_forrest)
  df_forrest$Dataset = rownames(df_forrest)
  idx = which(is.na(df_forrest$se))
  df_forrest$se[idx] = df_forrest$logFC[idx] / df_forrest$t[idx]

  df_forrest$Dataset = factor(df_forrest$Dataset, rev(df_forrest$Dataset))

  # set range to include zero
  rng = range(c(with(df_forrest, logFC-se), with(df_forrest, logFC+se)))
  if( rng[1] > 0 ){
    rng[1] = 0
  }
  if( rng[2] < 0 ){
    rng[2] = 0
  }

  # forrest plot
  fig1 = ggplot(df_forrest, aes(Dataset, logFC)) + geom_point(color="navy") + geom_errorbar(aes(ymin = logFC - se, ymax=logFC + se), width=0.1) + theme_bw(15) + theme(aspect.ratio= .2, plot.title = element_text(hjust = 0.5)) + coord_flip() + ylim(rng) + xlab(bquote(log[2]~fold~change)) + ggtitle(geneSymbol) + geom_hline(yintercept=0, color="red", linetype="dashed")

  # Adjusted p-value
  ymax = max(-log10(df_forrest$adj.P.Val))*1.05
  fig2 = ggplot(df_forrest, aes(Dataset, -log10(adj.P.Val))) + geom_bar(stat="identity", fill="navy") + theme_bw(12) + theme( aspect.ratio= .5, plot.margin = unit(c(0,0,0,0), "cm"),
    axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) + coord_flip() + ylab(bquote(-log[10]~adjusted~P)) + scale_y_continuous(expand=c(0,0), limits=c(0, ymax))

  # plot_grid(fig1, fig2, align='vh')
  # fig1 + fig2

  cbind(ggplotGrob(fig1), ggplotGrob(fig2))  
}


# Functions for gene set enrichment
###################################

# Read GMT file to list
#   reads either .gmt or gmt.gz
read.gmt = function(file){
    if( ! grepl("(\\.gmt$)|(\\.gmt\\.gz$)", file)[1] ){
        stop("Pathway information must be a .gmt file")
    }
    geneSetDB = readLines(file)
    geneSetDB = strsplit(geneSetDB, "\t")
    names(geneSetDB) = sapply(geneSetDB, "[", 1)
    geneSetDB = lapply(geneSetDB, "[", -1:-2)
    geneSetDB = lapply(geneSetDB, function(x) {
        x[which(x != "")]
    })
    return(geneSetDB)
}

# Map from HGNC to ENSEMBL gene ids
convert_gene_list_to_ensembl = function( x, geneInfo){
  geneInfo[ x,]$ensembl_gene_id
}


plot_sex_compare = function( tab, minp=1e-300 ){

  tab$Chrom = factor(tab$Chrom, c(1:22, "X", "Y"))
  tab$adj.P.Val = pmax(minp, tab$adj.P.Val)

  figList = lapply( levels(tab$Chrom), function(Chrom){

    tab2 = tab[tab$Chrom==Chrom,]
    tab2$logFC.adj = with(tab2, (adj.P.Val < 0.05) * logFC)
    tab2$Color = pmax(-1, pmin(1, with(tab2, (adj.P.Val < 0.05) * logFC)))
    tab2 = tab2[order(tab2$adj.P.Val, decreasing=TRUE),]
    limits = with(tab2, range( c(AveExpr+logFC, AveExpr)))

    ggplot(tab2, aes(AveExpr, AveExpr+logFC, label=Symbol)) + geom_point(aes(color=sign(logFC.adj), size=abs(logFC.adj))) + theme_bw(15) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + ggtitle(paste0("chr",Chrom)) + geom_abline(color="black", linetype="dashed") + geom_abline(intercept=1, slope=, color="red", linetype="dashed") + geom_abline(intercept=-1, slope=, color="blue", linetype="dashed") + geom_text_repel(data=subset(tab2, adj.P.Val<1e-10),  segment.size  = 0.2, seed=1,  point.padding = unit(.1, 'lines'),  box.padding = unit(.3, 'lines') , nudge_x=-1.5, nudge_y=1.5) + xlim(limits) + ylim(limits) + xlab(bquote(Male~log[2]~expression)) + ylab(bquote(Female~log[2]~expression)) + scale_color_gradientn(name = "logFC", colors = colorRampPalette(c("blue", "grey60", "red"))(50), limits=c(-1, 1)) + scale_size_continuous("|logFC|", limits=c(0,8))

    # tab2$Color[tab2$Symbol == "XIST"] = max(abs(tab2$Color[tab2$Symbol != "XIST"]))

    # ggplot(tab2, aes(AveExpr, AveExpr+logFC, label=Symbol)) + geom_point(aes(color=Color, shape=(Color!=0))) + theme_bw(15) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + ggtitle(paste0("chr",Chrom)) + geom_abline(color="black", linetype="dashed") + geom_abline(intercept=1, slope=, color="red", linetype="dashed") + geom_abline(intercept=-1, slope=, color="blue", linetype="dashed") + geom_text_repel(data=subset(tab2, adj.P.Val<0.01),  segment.size  = 0.2, seed=1,  point.padding = unit(.1, 'lines'),  box.padding = unit(.3, 'lines') , nudge_x=-1.5, nudge_y=1.5) + xlim(limits) + ylim(limits) + xlab(bquote(Male~log[2]~expression)) + ylab(bquote(Female~log[2]~expression)) + scale_color_gradientn(name = "logFC", colors = colorRampPalette(c("blue", "grey60", "red"))(50), limits=c(-1, 1)) + scale_shape_discrete("Significant")

  })
  names(figList) = levels(tab$Chrom)

  figList
}



















