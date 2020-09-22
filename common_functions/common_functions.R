suppressPackageStartupMessages({
library(qvalue)
library(AnnotationHub)
library(ensembldb)
library(MASS)
library(ggplot2)
library(ggbio)
library(Rfast)
library(viridis)
library(data.table)
library(ggraph)
library(ggplot2)
library(tidygraph)
library(irlba)
library(EnvStats)
library(Rfast)
library(GSEABase)
})

# For ENSEMBL id ENSG00000279457.4, return ENSG00000279457
trim_ensembl_ids = function(x){
  gsub("(.*)\\.(.*)", "\\1", x) 
}

# Add column Symbol using column gene_id storing ENSEMBL id
getGeneSymbol = function( df, column="row.names"){

  if( ! is(df, "data.frame") ){
    df = as.data.frame(df)
  }

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

  entrez = mapIds(ensdb, keys=ensIDs, keytype="GENEID", column="ENTREZID")
  df$Entrez = entrez

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
  fig1 = ggplot(df, aes(logFC.x, logFC.y, color=density)) + geom_point(size=.4) + theme_bw(size) + theme(aspect.ratio=1, legend.position="bottom", plot.title = element_text(hjust = 0.5)) + geom_abline(color="red") + xlab(paste(names(resList)[1], lab)) + ylab(paste(names(resList)[2], lab)) + geom_vline(xintercept=0, col="grey40", linetype="dashed") + geom_hline(yintercept=0, col="grey40", linetype="dashed") + xlim(-lim, lim) + ylim(-lim, lim) + scale_color_viridis() + geom_smooth(method="lm", se=FALSE, color="darkorange")

  if( !is.null(showGenes) ){
    df2 = df[df$Symbol.x %in% showGenes,]
    fig1 = fig1 + geom_text_repel(data=df2, aes(logFC.x, logFC.y, label=Symbol.x), segment.size=.5,  segment.color="black", color="black", force=1, nudge_x=.005, nudge_y=.1)
  }

  lab = "t-statistic"
  lim = with(df, max(abs(c(t.x, t.y))))
  df$density <- get_density(df$t.x, df$t.y, n = 100)
  fig2 = ggplot(df, aes(t.x, t.y, color=density)) + geom_point(size=.4) + theme_bw(size) + theme(aspect.ratio=1, legend.position="bottom", plot.title = element_text(hjust = 0.5)) + geom_abline(color="red") + xlab(paste(names(resList)[1], lab)) + ylab(paste(names(resList)[2], lab)) + geom_vline(xintercept=0, col="grey40", linetype="dashed") + geom_hline(yintercept=0, col="grey40", linetype="dashed") + xlim(-lim, lim) + ylim(-lim, lim) + scale_color_viridis() + geom_smooth(method="lm", se=FALSE, color="darkorange")

  if( !is.null(showGenes) ){
    df2 = df[df$Symbol.x %in% showGenes,]
    fig2 = fig2 + geom_text_repel(data=df2, aes(t.x, t.y, label=Symbol.x), segment.size=.5,  segment.color="black", color="black", force=1, nudge_x=.005, nudge_y=.3)
  }

  res = with(df, cor.test(logFC.x, logFC.y, method="spearman"))
  tab_corr = data.frame(R_spearman = res$estimate, P = res$p.value)
  
  pi1_x = 1 - qvalue(df$P.Value.x)$pi0
  pi1_y = 1 - qvalue(df$P.Value.y)$pi0

  pi_single = c(pi1_x, pi1_y)
  names(pi_single) = names(resList)

  p = with(df, P.Value.y[adj.P.Val.x < 0.05])
  n.de.x = length(p)
  if( n.de.x > 0){
    pi_discovery_x = tryCatch( 
      1 - qvalue(p)$pi0, 
      error = function(e){
        # must have a p-value > 0.95 to work
        1 - qvalue(c(1,p))$pi0
        })
  }else{
      pi_discovery_x = 0
  }

  p = with(df, P.Value.x[adj.P.Val.y < 0.05])
  n.de.y = length(p)
  if( n.de.y > 0){
    pi_discovery_y = tryCatch( 
      1 - qvalue(p)$pi0, 
      error = function(e){
        NA
        })
  }else{
    pi_discovery_y = 0
  }

  tab_pi = data.frame(Discovery = names(resList), pi1 = c(pi_discovery_x, pi_discovery_y))

  list(fig_logFC  = fig1,
      fig_tstat   = fig2, 
      pi_single   = pi_single,
      n.de        = c(n.de.x, n.de.y),
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
  # Currently extract simple results
  # other results may be relevant if meta-analysing more than two studies
  resTable = lapply(resRMA, function(x){
    data.frame( logFC     = x$beta,
                se        = x$se,
                P.Value   = x$pval,
                # adj.P.Val = p.adjust(x$pval, "fdr"),
                Q         = x$QE,
                N         = x$k,
                I2        = x$I2)
    })
  resTable = do.call("rbind", resTable)
  resTable$adj.P.Val = p.adjust(resTable$P.Value, "fdr")

  resTable[order(resTable$P.Value, decreasing=FALSE),]
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



plot_forrest = function(resList, geneSymbol, plot2="p"){

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
  fig1 = ggplot(df_forrest, aes(Dataset, logFC)) + geom_point(color="navy") + geom_errorbar(aes(ymin = logFC - se, ymax=logFC + se), width=0.1) + theme_bw() + theme(aspect.ratio= .2, plot.title = element_text(hjust = 0.5)) + coord_flip() + ylim(rng) + ggtitle(geneSymbol) + geom_hline(yintercept=0, color="red", linetype="dashed") + ylab(bquote(log[2]~fold~change)) + xlab('')

  # Adjusted p-value
  if( plot2 == "adjusted" ){
    df_forrest$adj.P.Val = pmax(1e-300, df_forrest$adj.P.Val)
    ymax = max(-log10(df_forrest$adj.P.Val))*1.05
    fig2 = ggplot(df_forrest, aes(Dataset, -log10(adj.P.Val))) + ylab(bquote(-log[10]~adjusted~P))  
  }else{

    df_forrest$P.Value = pmax(1e-300, df_forrest$P.Value)

    ymax = max(-log10(df_forrest$P.Value))*1.05
    fig2 = ggplot(df_forrest, aes(Dataset, -log10(P.Value))) + ylab(bquote(-log[10]~P)) 
  }

  fig2 = fig2 + geom_bar(stat="identity", fill="navy") + theme_bw() + theme( aspect.ratio= .5, plot.margin = unit(c(0,0,0,0), "cm"),
    axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) + coord_flip() + scale_y_continuous(expand=c(0,0), limits=c(0, ymax)) + geom_hline(yintercept=-log10(0.05), color="red", linetype="dashed")

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


plot_sex_compare = function( tab, minp=1e-300, ngenes=10, maxSize=10 ){

  tab$Chrom = factor(tab$Chrom, c(1:22, "X", "Y"))
  tab = tab[!is.na(tab$Chrom),] # remove non-standard chromosomes
  tab$adj.P.Val = pmax(minp, tab$adj.P.Val)
  tab$logFC = pmin(maxSize, tab$logFC)
  tab$logFC = pmax(-maxSize, tab$logFC)

  figList = lapply( levels(tab$Chrom), function(Chrom){

    tab2 = tab[tab$Chrom==Chrom,]
    tab2$logFC.adj = with(tab2, (adj.P.Val < 0.05) * logFC)
    tab2$Color = pmax(-1, pmin(1, with(tab2, (adj.P.Val < 0.05) * logFC)))
    tab2 = tab2[order(tab2$adj.P.Val, decreasing=TRUE),]
    limits = with(tab2, range( c(AveExpr+logFC, AveExpr)))

    tab3 = subset(tab2, adj.P.Val<0.05)
    ngenes = min(c(ngenes, nrow(tab3)))
    tab3 = tab3[order(tab3$adj.P.Val),]
    tab3 = tab3[1:ngenes,]

    ggplot(tab2, aes(AveExpr, AveExpr+logFC, label=Symbol)) + geom_point(aes(color=sign(logFC.adj), size=abs(logFC.adj))) + theme_bw(15) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + ggtitle(paste0("chr",Chrom)) + geom_abline(color="black", linetype="dashed") + geom_abline(intercept=1, slope=, color="red", linetype="dashed") + geom_abline(intercept=-1, slope=, color="blue", linetype="dashed") + geom_text_repel(data=tab3, segment.size  = 0.2, seed=1,  point.padding = unit(.1, 'lines'),  box.padding = unit(.3, 'lines') , nudge_x=-1.5, nudge_y=1.5) + xlim(limits) + ylim(limits) + xlab(bquote(Male~log[2]~expression)) + ylab(bquote(Female~log[2]~expression)) + scale_color_gradientn(name = "logFC", colors = colorRampPalette(c("blue", "grey60", "red"))(50), limits=c(-1, 1), guide=FALSE) + scale_size_continuous("|logFC|", limits=c(0,maxSize)) 
  })
  names(figList) = levels(tab$Chrom)

  figList
}

plot_circos_logFC = function( tab, maxLogFC=2){

  tab$Chrom = factor(tab$Chrom, c(1:22, "X", "Y"))
  tab = tab[!is.na(tab$Chrom),] # remove non-standard chromosomes

  tab$logFC.adj = with(tab, (adj.P.Val < 0.05) * logFC)
  tab = tab[order(tab$adj.P.Val, decreasing=FALSE),]

  # Get TSS locations
  #-------------------
  ah <- AnnotationHub()
  ensdb = ah[["AH69187"]] # ENSEMBL v96

  tab$gene_id = trim_ensembl_ids(rownames(tab))

  prom = promoters(ensdb, 0, 0, filter=GeneIdFilter(tab$gene_id))

  df_prom = data.table(as.data.frame(prom))
  df_prom = df_prom[,data.frame(TSS=min(start)),by=gene_id]

  tab2 = merge(tab, df_prom, by="gene_id")

  # Get genome info
  hg38 = seqinfo(ensdb)
  gr_hg38 = GRanges(hg38)
  gr_hg38 = gr_hg38[seqnames(gr_hg38) %in% tab$Chrom]
  gr_hg38 = keepStandardChromosomes(gr_hg38)
  gr_hg38 = sortSeqlevels(dropSeqlevels(gr_hg38, "MT"))

  # great GRanges for plotting
  gr = with(tab2, GRanges(Chrom, IRanges(TSS, TSS+1), score =logFC.adj, seqinfo=hg38)) 
  gr = keepStandardChromosomes(gr)
  gr = sortSeqlevels(dropSeqlevels(gr, "MT"))
  gr$score2 = pmax(-maxLogFC, pmin(maxLogFC, gr$score))

  # create circos plot
  fig = ggbio() + circle(gr_hg38, geom = "text", aes(label = seqnames), vjust = 0, size = 3)

  fig + circle(gr, geom = "point", aes(y = score2, color=factor(sign(score2), -1:1)),size=1, grid = TRUE, radius = 15, trackWidth=23, grid.background="white", grid.line="grey70", space.skip=.002) + ylab(bquote(log[2]~fold~change)) + scale_color_manual(name='logFC', values = c("red", "grey60", "blue")) + theme(legend.position="right", plot.title = element_text(hjust = 0.5))
}


downloadFile <- function(id){
  fread(synGet(id)$path, data.table = F)
}
downloadFile_version <- function(id , version){
  fread(synGet(id, version = version)$path, data.table = F)
}



  # tab$Chrom = factor(tab$Chrom, c(1:22, "X", "Y"))
  # tab = tab[!is.na(tab$Chrom),] # remove non-standard chromosomes

  # tab$logFC.adj = with(tab, (adj.P.Val < 0.05) * logFC)
  # tab = tab[order(tab$adj.P.Val, decreasing=FALSE),]
  # limits = c(-max(abs(tab$logFC)), max(abs(tab$logFC)))


  # ggplot(tab, aes(Chrom, logFC, color=factor(sign(logFC.adj), -1:1), size=abs(logFC.adj))) + geom_point() + theme_bw(15) + theme(aspect.ratio=1/2) + scale_size_continuous("|logFC|", limits=c(0,max(limits))) + scale_color_manual(values = c("red", "grey60", "blue"))


  # ah <- AnnotationHub()
  # ensdb = ah[["AH69187"]] # ENSEMBL v96

  # tab$gene_id = trim_ensembl_ids(rownames(tab))

  # prom = promoters(ensdb, 0, 0, filter=GeneIdFilter(tab$gene_id))

  # df_prom = data.table(as.data.frame(prom))
  # df_prom = df_prom[,data.frame(TSS=min(start)),by=gene_id]

  # tab2 = merge(tab, df_prom, by="gene_id")



  # ggplot(gr, aes(TSS, logFC, color=factor(sign(score), -1:1), size=abs(score))) + geom_point() + theme_bw(15) + theme(aspect.ratio=1/2) + scale_size_continuous("|logFC|", limits=c(0,max(limits))) + scale_color_manual(values = c("red", "grey60", "blue"))


  # fig = ggbio() + circle(hg19sub, geom = "ideo", fill = "gray70") #+ circle(hg19sub, geom = "scale", size = 2) +
  # fig = fig + circle(hg19sub, geom = "text", aes(label = seqnames), vjust = 0, size = 3)



  # hg38 = seqinfo(ensdb)
  # # gr_hg38 = GRanges(seqnames(hg38), IRanges(1, seqlengths(hg38)), seqinfo=hg38)
  # gr_hg38 = GRanges(hg38)
  # gr_hg38 = gr_hg38[seqnames(gr_hg38) %in% tab$Chrom]
  # gr_hg38 = keepStandardChromosomes(gr_hg38)
  # gr_hg38 = sortSeqlevels(dropSeqlevels(gr_hg38, "MT"))


  # gr = with(tab2, GRanges(Chrom, IRanges(TSS, TSS+1), score =logFC.adj, seqinfo=hg38)) 
  # gr = keepStandardChromosomes(gr)
  # gr = sortSeqlevels(dropSeqlevels(gr, "MT"))
  # gr$score2 = pmax(-2, pmin(2, gr$score))

  # fig = ggbio() + circle(gr_hg38, geom = "text", aes(label = seqnames), vjust = 0, size = 3)

  # fig + circle(gr, geom = "point", aes(y = score2, color=factor(sign(score2), -1:1)),size=1, grid = TRUE, radius = 15, trackWidth=23, grid.background="white", grid.line="grey70", space.skip=.002) + scale_size(range = c(1, 2.5)) + scale_size_continuous("|logFC|", limits=c(0,max(limits))) + ylab(bquote(log[2]~fold~change)) + scale_color_manual(values = c("red", "grey60", "blue")) + theme(legend.position="none")





#' Plot mean variance trend
#' 
#' Plot mean variance trend
#' 
#' @param x object returned by voom(), voomWithQualityWeights(), voomWithDreamWeights() or rationale()
#' @param y missing
#' @param ... other arguments
#' 
#' @import ggplot2
#' @export
plot_voom = function(x,  y, ...){

    # pass R check
    statusAB = NA

    if( is.null(x$voom.xy) || is.null(x$voom.line)){
      stop("x does not contain the trend information.\nvoom() must be run with save.plot=TRUE")
    }

    # points
    if( is.null(x$statusAB) ){
      main = "Voom: Mean-variance trend"
      isRationale = FALSE
    }else{
      main = "Rationale: Mean-variance trend"
      isRationale = TRUE
    }

  df = data.frame(x = x$voom.xy$x,
          y = x$voom.xy$y)

  # trend line
  df.line = data.frame(x = x$voom.line$x,
            y = x$voom.line$y)

  if( isRationale ){
    df$statusAB = x$statusAB
    fig = ggplot(df, aes(x,y, color=statusAB))
  }else{
    fig = ggplot(df, aes(x,y))
  }
  
  fig + geom_point(size=0.01) + theme_bw(15) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + geom_line(data=df.line, aes(x,y), color="red", size=1.2) + xlab(bquote(log[2](count + 0.5))) + ylab(expression( sqrt("standard deviation"))) + ggtitle(main) + ylim(0, max(df$y)) + scale_color_manual(name="Status", values=c("#3758ba", "#33e65c")) + guides(color = guide_legend(override.aes = list(size = 1))) 
  }




run_zenith = function(fit, coefs, gs.collection, n_genes_min=10, n_genes_max=5000 ){
    
  # Map from Ensembl genes in geneSets_GO to 
  # from trimmed Ensembl names from RNA-seq data 
  geneSets.index = ids2indices( gs.collection, trim_ensembl_ids(rownames(fit)))
        
  # filter by size of gene set
  geneSets.index = geneSets.index[sapply(geneSets.index, function(x) (length(x) >= n_genes_min) & (length(x) <= n_genes_max)) ]

  # run zenith for each coefficient
  res = lapply( coefs, function(coef){

    df = zenith(fit, coef, geneSets.index, squaredStats=TRUE)
    df$Set = rownames(df)
    rownames(df) = c()
    df$Coef = coef
    df
    })
  do.call(rbind, res)
}
 
# , only_autosomes=rep(FALSE, length(coefs))
    # if( only_autosomes[which(coefs == coef)] ){
    #   tab = coef(fit) %>% getGeneSymbol()
    #   geneids = rownames(tab)[tab$Chrom %in% 1:22]
    #   fit = fit[geneids,]
    # }

zenith_meta_analysis = function(res_zenith, coefTest){

  coefTest = unique(unlist(lapply(res_zenith, function(x) unique(x$Coef))))
  rsrc = unique(unlist(lapply(res_zenith, function(x) unique(x$Resource))))

  res_meta = mclapply( rsrc, function(Resource){
    res = mclapply( coefTest, function(coef){

      df = lapply(res_zenith, function(df){
        df[(df$Resource == Resource) & (df$Coef == coef),]
      })

      df_merge = merge(df[[1]], df[[2]], by="Set")

      df_meta = data.frame(Set      = df_merge$Set,
                           NGenes   = df_merge$NGenes.x,
                           Resource = df_merge$Resource.x,
                           Coef     = coef)
      rownames(df_meta) = df_merge$Row.names

      for(i in 1:nrow(df_merge)){
        effect = c(df_merge$delta.x[i], df_merge$delta.y[i])
        se = c(df_merge$se.x[i], df_merge$se.y[i])

        ma = suppressWarnings(rma(yi=effect, sei=se, method="FE"))

        df_meta$effect[i] = ma$b[1]
        df_meta$se[i] = ma$se
        df_meta$PValue[i] = ma$pval
      }

      df_meta$FDR = p.adjust(df_meta$PValue, "BH")
      df_meta[order(df_meta$PValue),]
    }, mc.cores=4)
    do.call(rbind, res)
  })

  do.call(rbind, res_meta)
}



get_difference_network = function(resid.lst, METADATA, variable){

  # For each cohort
  C.lst = lapply( resid.lst, function(resid){

    i = match(colnames(resid), rownames(METADATA))
    info = METADATA[i,]

    # get correlation matrix for each category
    C.lst = lapply( levels(info[[variable]]), function(key){
      j = (info[[variable]] == key)
      cora( t(resid[,j]) )  
    })
    names(C.lst) = levels(info[[variable]])

    C.lst
  })
  names(C.lst) = names(resid.lst)

  C.lst
}

# like get_difference_network, but don't compute difference
get_networks = function(resid.lst, METADATA, variable){

  # For each cohort
  nets = lapply( resid.lst, function(resid){

    i = match(colnames(resid), rownames(METADATA))
    info = METADATA[i,]

    # get correlation matrix for each category
    C.lst = lapply( levels(info[[variable]]), function(key){
      j = (info[[variable]] == key)
      cora( t(resid[,j]) )  
    })
    names(C.lst) = levels(info[[variable]])

    C.lst
    })
  names(nets) = names(resid.lst)

  nets
}



sLED_adapt = function(Y1, Y2, npermute=c(1000,1e6), BPPARAM=SerialParam()){

  # perform permutations until p-values is precise enough
  # if not precise enough
  a = log10(npermute[1])
  b = log10(npermute[2])
  permArray = 10^seq(a,b, length.out=b-a +1)

  for( nperm in round(permArray) ){

      bp = BPPARAM

      if(nperm <= 10000) bp = SerialParam()

      # compare correlation structure with sLED
      res = decorate:::.sLED(X=Y1, Y=Y2, npermute=nperm, BPPARAM=bp, rho=1, sumabs.seq=.2)

      if( res$pVal * nperm > 10){
        break
      }
  }

  # stats is always positive, so multiply by -1 if negative
  res$stats = (-1)^(res$sign=='neg') * res$stats

  res
}

sLED_perm = function(Y, info, variable, nperm, rho=1, sumabs.seq=.2){

  eval_stat = function(idx, sumabs.seq){
    lvls = levels(info[[variable]])

    C.lst = lapply( levels(info[[variable]]), function(lvl){
      j = (info[[variable]] == lvl)
      cora( Y[idx,][j,]) 
    })
    names(C.lst) = levels(info[[variable]])

    D = C.lst[[lvls[2]]] - C.lst[[lvls[1]]]
    # svd(C.diff, nu=0, nv=0)$d[1]
    # partial_eigen(C.diff, n=1)$values

    res = sLED:::sLEDTestStat(D, rho=rho, sumabs.seq = sumabs.seq)

    data.frame(sumabs.seq, 
              n.active=rowSums(res$leverage > 0), 
              stats=pmax(0, res$stats))
  }

  # evaluate on real data
  res = eval_stat(1:nrow(info), sumabs.seq)#=seq(.2, 1, length.out=5))
  val = res$stats

  # permutation
  val_perm = sapply( seq_len(nperm), function(k){
    idx = sample.int(nrow(info), nrow(info))
    eval_stat( idx, sumabs.seq )$stats
  })

   # estimate null distribution as a gamma
  fit = egamma( val_perm )
  up = pgamma(val, shape=fit$parameters[1], scale=fit$parameters[2], lower.tail=FALSE)
  down = pgamma(val, shape=fit$parameters[1], scale=fit$parameters[2])
  P.Value = 2*min(c(up, down))

  data.frame(p = P.Value, n.active = res$n.active)
}





test_differential_correlation = function(resid.lst, C.diff.discovery, dynamicColors, METADATA, variable, useSLED=FALSE, sumabs.seq, nperm){

  col.array = unique(dynamicColors)

  df_test = lapply( col.array, function(col){

    # get genes in this cluster
    geneid = rownames(C.diff.discovery)[which(dynamicColors==col)]

    res = lapply( names(resid.lst), function(key){

      if(runif(1) < .1) message(key, ' ', col, ' ', match(col, col.array), '/', length(col.array))

      # extact expression residuals for genes in this cluster
      Y = t(resid.lst[[key]][geneid,])
     
      # extract metadata
      i = match(rownames(Y), rownames(METADATA))
      info = METADATA[i,]

      if( ! useSLED ){
        # eval statistical hypothesis
        res = boxM_permute( Y, info[[variable]])

        res.sled = list(p = res$p.value, n.active=ncol(Y))

      }else{
        lvl = levels(info[[variable]])
        Y1 = Y[info[[variable]] == lvl[1],]
        Y2 = Y[info[[variable]] == lvl[2],]

        # res_sLED = sLED_adapt( scale(Y1), scale(Y2), npermute=c(500,50000), BPPARAM=BPPARAM)

        # D = cora(Y2) - cora(Y1)
        # p.sled = max(0, sLED:::sLEDTestStat(D, rho=1, sumabs.seq = .5)$stats)

        res.sled = sLED_perm(Y, info, variable, nperm=nperm, sumabs.seq=sumabs.seq)
      }

      data.frame( Module  = col, 
                  # P.Value       = res$p.value, 
                  P.Value  = res.sled$p,
                  n.active  = res.sled$n.active,
                  n.genes  = length(geneid))
    })
    res = do.call(rbind, res)
    res$Cohort = names(resid.lst)
    rownames(res) = c()
    res
  })#, mc.cores=4)

  # Format results data.frame
  df_test = do.call(rbind, df_test)
  df_test = df_test[!is.na(df_test$P.Value),]
  df_test = df_test[order(df_test$P.Value),]
  df_test$Module = factor(df_test$Module, unique(df_test$Module))

  # Compute FDR separately for each cohort
  df_test$FDR = rep(NA, nrow(df_test))
  i = which(df_test$Cohort == 'MSSM-Penn-Pitt')
  df_test$FDR[i] = p.adjust(df_test$P.Value[i], "fdr")

  i = which(df_test$Cohort == 'NIMH-HBCC')
  df_test$FDR[i] = p.adjust(df_test$P.Value[i], "fdr")
  # df_test$FDR[i] = qvalue(df_test$P.Value[i])$qvalues

  df_test
}



test_differential_correlation_interaction = function(resid.lst, C.diff.discovery, dynamicColors, METADATA, nperm=1000){

  col.array = unique(dynamicColors)

  df_test = lapply( col.array, function(col){

    # get genes in this cluster
    geneid = rownames(C.diff.discovery)[which(dynamicColors==col)]

    res = lapply( names(resid.lst), function(key){

      message(key, ' ', col, ' ', match(col, col.array), '/', length(col.array))

      # extact expression residuals for genes in this cluster
      Y = t(resid.lst[[key]][geneid,])
      # W = t(vobj.lst[[key]][geneid,]$weights)

      # extract metadata
      i = match(rownames(Y), rownames(METADATA))
      info = METADATA[i,]

      eval_stat = function(idx){
        lvls = levels(info[[variable]])

        C.lst = lapply( levels(info[[variable]]), function(lvl){
          j = (info[[variable]] == lvl)
          cora( Y[idx,][j,geneid]) 
        })
        names(C.lst) = levels(info[[variable]])

        # (C[[lvls[4]]] - C[[lvls[3]]]) - (C[[lvls[2]]] - C[[lvls[1]]])
        C_alt = (C.lst[[lvls[4]]] - C.lst[[lvls[3]]])
        C_baseline = (C.lst[[lvls[2]]] - C.lst[[lvls[1]]])
        rm(C.lst)

        C.diff = C_alt - C_baseline
        # svd(C.diff, nu=0, nv=0)$d[1]
        # partial_eigen(C.diff, n=1)$values
        max(0, sLED:::sLEDTestStat(C.diff, rho=1, sumabs.seq = .5)$stats)
      }

      # evaluate on real data
      val = eval_stat(1:nrow(info))

      # permutation
      val_perm = sapply( seq_len(nperm), function(k){
        idx = sample.int(nrow(info), nrow(info))
        eval_stat( idx )
      })

      # estimate null distribution as a gamma
      fit = egamma( val_perm )
      up = pgamma(val, shape=fit$parameters[1], scale=fit$parameters[2], lower.tail=FALSE)
      down = pgamma(val, shape=fit$parameters[1], scale=fit$parameters[2])
      P.Value = 2*min(c(up, down))
  
      data.frame( Module        = col, 
                  P.Value       = P.Value, 
                  n.genes       = length(geneid))
    })
    res = do.call(rbind, res)
    res$Cohort = names(resid.lst)
    rownames(res) = c()
    res
  })

  # Format results data.frame
  df_test = do.call(rbind, df_test)
  df_test = df_test[!is.na(df_test$P.Value),]
  df_test = df_test[order(df_test$P.Value),]
  df_test$Module = factor(df_test$Module, unique(df_test$Module))

  # Compute FDR separately for each cohort
  df_test$FDR = rep(NA, nrow(df_test))
  i = which(df_test$Cohort == 'MSSM-Penn-Pitt')
  df_test$FDR[i] = p.adjust(df_test$P.Value[i], "fdr")

  i = which(df_test$Cohort == 'NIMH-HBCC')
  df_test$FDR[i] = p.adjust(df_test$P.Value[i], "fdr")
  # df_test$FDR[i] = qvalue(df_test$P.Value[i])$qvalues

  df_test
}




# Plots for each module
plot_module = function( clusterID, df_test, METADATA, dynamicColors, C.diff.discovery, resid.lst, variable, key, base_size=11, lvlidx = 1:2, zcutoff=1.3){

  clusterID = as.character(clusterID)

  # get genes in this clusters
  geneid = rownames(C.diff.discovery)[which(dynamicColors==clusterID)]

  # extract residuals for these genes
  Y = t(resid.lst[[key]][geneid,])

  # convert from ENSEMBL to SYMBOL
  df_tmp = data.frame(ENSEMBL = colnames(Y)) %>% getGeneSymbol('ENSEMBL')
  i = match(df_tmp$ENSEMBL, colnames(Y))  
  colnames(Y)[i] = df_tmp$Symbol

  # extract metadata
  i = match(rownames(Y), rownames(METADATA))
  info = METADATA[i,]

  # perform hypothesis test
  res = boxM_permute( Y, info[[variable]])

  FDR = df_test[with(df_test, Cohort==key & Module == clusterID),'FDR']

  main = paste(key, clusterID, 'FDR =', format(FDR, digits=3))

  # Evaluate correlation matrix for each category
  lvl = levels(info[[variable]])
  if( length(lvlidx) == 2 ){
    C1 = cor(Y[info[[variable]] == lvl[lvlidx[1]],])
    C2 = cor(Y[info[[variable]] == lvl[lvlidx[2]],])
  }else if( length(lvlidx) == 4 ){

    C.lst = lapply( levels(info[[variable]]), function(l){
      j = (info[[variable]] == l)
      cora( Y[j,] ) 
    })
    names(C.lst) = levels(info[[variable]])

    # (C[[lvls[4]]] - C[[lvls[3]]]) - (C[[lvls[2]]] - C[[lvls[1]]])
    C2 = (C.lst[[lvl[4]]] - C.lst[[lvl[3]]])
    C1 = (C.lst[[lvl[2]]] - C.lst[[lvl[1]]])

    lvl = c(paste0('(',lvl[4], ' - ',lvl[3], ')'),
             paste0('(',lvl[2], ' - ',lvl[1], ')'))
    lvlidx = 1:2

  }else{
    stop(length(lvlidx))
  }

  res = sLED:::sLEDTestStat(C2-C1, rho=1, sumabs.seq = .5)
  df_leverage = data.frame(Gene = colnames(C2), leverage = t(res$leverage), stringsAsFactors=FALSE)
  genes.has.leverage = df_leverage$Gene[df_leverage$leverage > 0]

  # reorder for clustering
  hcl = hclust(as.dist(1 - C1), method="ward.D2")
  C1 = C1[hcl$order,hcl$order]
  C2 = C2[hcl$order,hcl$order]

  # reorder based on clustering
  df_leverage$Gene = factor(df_leverage$Gene, rownames(C2))

  ylim = max(df_leverage$leverage)*1.05
  fig.leverage = ggplot(df_leverage, aes(Gene, leverage)) + geom_bar(stat="identity", fill="navy") + theme_bw() + theme(aspect.ratio=5) + coord_flip() + theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlab("") + ylab("Leverage") + scale_y_continuous(limits=c(0, ylim), expand=c(0,0))

  # create data.frame to compare correlation values
  ind <- which( upper.tri(C1,diag=FALSE) , arr.ind = TRUE )

  df2 = data.frame( col = dimnames(C1)[[2]][ind[,2]] ,
              row = dimnames(C2)[[1]][ind[,1]] ,
              cor1 = C1[ ind ],
              cor2 = C2[ ind ],
              pair = '' ) 
  i = order(with(df2, abs(cor1 - cor2)), decreasing=TRUE)[1:3]

  df2$pair[i] = with(df2[i,], paste(col, '/',row))

  lim = range(c(df2$cor1, df2$cor2))

  fig1 = ggplot(df2, aes(cor1, cor2, label=pair)) + geom_abline(color="grey", linetype="dashed") + geom_point(aes(color=ifelse(pair=='', "grey", "red"))) + theme_bw(base_size) + theme(aspect.ratio = 1, legend.position='none', plot.title = element_text(hjust = 0.5))  + geom_text_repel(direction='x', nudge_x=1, hjust=1, segment.size = 0.2, box.padding=.2, size=base_size/2) + scale_color_manual(values=c("grey50", "red")) + xlim(lim) + ylim(lim) + ggtitle(main) + xlab(paste('Correlation in', lvl[lvlidx[1]])) + ylab(paste('Correlation in', lvl[lvlidx[2]])) 

  # upper is C2 and lower is C1, since C1 is baseline, is lvl[1]
  C = C2
  C[lower.tri(C)] = C1[lower.tri(C1)]
  diag(C) = NA
  df = reshape2::melt(C)

  colsFun = colorRampPalette( c("blue", "white","red"))
  color = colsFun(1000)

  fig2 = ggplot(df, aes(Var1, Var2)) + geom_tile(aes(color = value, 
        fill = value)) + scale_color_gradientn(name = "Correlation", 
        colours = color, limits = c(-1, 1), na.value = "grey") + 
        scale_fill_gradientn(name = "Correlation", colours = color, 
            limits = c(-1, 1), na.value = "grey") + theme_bw(base_size) + 
        theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5), 
            legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text.x=element_text(angle=45, hjust=1)) + ggtitle(main) + ylab(lvl[lvlidx[2]]) + xlab(lvl[lvlidx[1]])

  # fig2 = plot_grid(fig2 + theme(legend.position="left", aspect.ratio=1), fig.leverage, ncol=2, align='v', rel_widths=c(4, 1), axis='t')      

  # plot network   
  D = C2 - C1
  main = paste(key, clusterID, paste0('(', lvl[lvlidx[2]], ' - ', lvl[lvlidx[1]], ')'), 'FDR =', format(FDR, digits=3))

  D = D[genes.has.leverage,genes.has.leverage]
  fig3 = plot_corr_network( D , base_size=base_size, zcutoff=zcutoff) + ggtitle(main)

  # fig_merge = plot_grid( fig1, fig2, nrow=1)
  # plot_grid( fig_merge, fig3, nrow=2, rel_heights=c(1,1.5) )  
  list(fig1, fig2, fig3, fig.leverage)  
}

# plot_module( mod, df_test, METADATA, dynamicColors, C.diff.discovery, resid.add[[resVersion]], variable, "MSSM-Penn-Pitt", 5) 




plot_corr_network = function( C, zcutoff = 1.3, seed=1, base_size=11){

  df_net = data.frame(t(combn(colnames(C), 2)))
  colnames(df_net) = c('name.from','name.to')

  df_net$weight = apply(df_net, 1, function(x){
    C[x['name.from'], x['name.to']]
  })
  # threshold weight at 1
  df_net$weight = pmax(-1.2, pmin(1.2, df_net$weight))

  df_net$z = scale(df_net$weight)
  df_net = df_net[abs(df_net$z) >= zcutoff,]

  node_names = unique(c(df_net$name.from, df_net$name.to))

  C_sub = C[node_names,node_names]

  df_net$from = match(df_net$name.from, node_names)
  df_net$to = match(df_net$name.to, node_names)

  net = tbl_graph(nodes = data.frame(name=node_names), edges = df_net[,c('from','to','weight')], directed=FALSE)
       
  # igraph_layouts <- c('star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 
  #                     'randomly', 'fr', 'kk', 'drl', 'lgl')

  # set seed for reproducability
  set.seed(seed)

  # plot
  ggraph(net, layout = "stress") + #igraph", algorithm='graphopt') + 
    geom_edge_link(aes(width = abs(weight), color=weight), alpha = 1) + 
    scale_edge_width(name='|Corr Diff|',range = c(0.2, 2)) +
    geom_node_point(size=base_size, color="grey80") +
    geom_node_text(aes(label = name), repel = FALSE,size=base_size/3) +
    scale_edge_color_gradient2(name='Corr Diff', low="blue", mid="white", high="red", lim=c(-1.2,1.2)) +
    labs(edge_width = "Correlation")  + 
    theme_graph(base_size=base_size, title_size=base_size) + 
    theme(aspect.ratio=1,plot.title = element_text(hjust = 0.5))
}
# plot_corr_network(C,base_size=4)






# df_net = data.table(df_net)

# df_a = df_net[,data.frame(minAbsZ=max(abs(z))), by="name.from"]
# df_b = df_net[,data.frame(minAbsZ=max(abs(z))), by="name.to"]

# df = merge(df_a, df_b, by.x="name.from", by.y="name.to")

# df2 = df[,data.frame(z=max(c(minAbsZ.x, minAbsZ.y))),by="name.from"]

# min(df2$z)

# plot_corr_network( C1-C2, zcutoff=1.568965)




enrich.test = function(df_module, gs_set, testModules = unique(df_module$Module), minSize=50 ){

  # only keep genes that are in at least one gene set
  df_module = df_module[df_module$ENSEMBL %in% unique(unlist(geneIds(gs_set))),]

  # Map from Ensembl genes in geneSets_GO to 
  # from trimmed Ensembl names from RNA-seq data 
  gs.list = limma::ids2indices( recodeToList(gs_set), df_module$ENSEMBL)
     
  # filter by size of gene set
  n_genes_in = minSize
  gs.list = gs.list[sapply(gs.list, length) >= n_genes_in]

  # for each module
  res = lapply( testModules, function(mod){

    # Indicator for Module
    isInModule = rep(0, nrow(df_module))
    isInModule[df_module$Module == mod] = 1 

    # for each gene set
    df = lapply( gs.list, function(gsidx){

      # Indicator for Gene
      isInGeneset = rep(0, nrow(df_module))
      isInGeneset[gsidx] = 1

      tab = table(isInGeneset, isInModule)
        
      fit = fisher.test(tab)

      data.frame( Module = mod,
            n.genes = sum(tab[,2]),
            # Geneset = gsid,
            n.overlap = tab[2,2],
            OR = fit$estimate,
            p.value = fit$p.value)
    })
    df = do.call(rbind, df)
    df$Geneset = names(gs.list) 
    df = df[,c("Module", "n.genes", "Geneset", "n.overlap", "OR", "p.value")]
    df = df[!is.na(df$p.value) & (df$p.value < 0.99999) & (df$OR !=0),]
    df$FDR = p.adjust(df$p.value, 'fdr')
    df = df[order(df$p.value),]

    df
  })
  do.call(rbind, res)
}



plot_enrich = function( df, module, base_size=11 ){

  df2 = df[df$Module == module,]
  df2$Geneset = factor(df2$Geneset, df2$Geneset)

  ylim = max(-log10(df$p.value))
  ggplot(df2[1:10,], aes(Geneset, -log10(p.value))) + geom_bar(stat='identity', fill=module) + theme_bw(base_size) + theme(aspect.ratio=1,plot.title = element_text(hjust = 0.5)) + coord_flip() + ylab(bquote(-log[10]~P)) + xlab('') + scale_y_continuous(expand=c(0,0), lim=c(0, ylim*1.05)) 
}





enrich_module_DE = function(df_module){ 
  # create gene set lists from modules
  mods = unique(df_module$Module)
  modules.gs = lapply( mods, function(mod){
    df_module$ENSEMBL[df_module$Module == mod]
    })
  names(modules.gs) = unique(df_module$Module)

  keys = c(names(df_meta), names(df_meta.inter_sex.disease), "SCHEMA", "ASD")

  res_enrich_DE = lapply( keys, function(key){

    if( key == "SCHEMA"){
      df = df_schema2
      df$stat = df$chisq
    }else if( key == "ASD"){
      df = df_asd2
      df$stat = df$chisq
    }else{
      df = df_meta[[key]]
      if( is.null(df) ) df = df_meta.inter_sex.disease[[key]]
      df$stat = with(df, (logFC/se)^2)
    }   
    
    idx = match( trim_ensembl_ids(rownames(df)), df_module$ENSEMBL)
    modules.list = limma::ids2indices( modules.gs, df_module$ENSEMBL[idx])
         
    res = cameraPR( df$stat, modules.list )
    res$PValue = res$PValue / 2
    res$FDR = c()
    res$Geneset = rownames(res)
    rownames(res) = c()
    res$Test = key
    res
  })
  res_enrich_DE = do.call(rbind, res_enrich_DE)

  res_enrich_DE
}



apply_permutation_null = function(df_test){
  df2 = lapply( unique(df_test$Cohort), function(chrt){
    df = df_test[df_test$Cohort == chrt,]
    if( chrt != 'NIMH-HBCC' ){
      df$P.Value = calc_p(df$P.Value, beta_perm[[chrt]])
    }
    df$FDR = p.adjust(df$P.Value, "fdr")
    df
  })
  df2 = do.call(rbind, df2)
  df2
}

module_meta_analysis = function(df2){

  df3 = lapply( unique(df2$Module), function(mod){
    # combine p-values with Fisher's method
    p.meta = sumlog(df2$P.Value[df2$Module == mod])$p

    data.frame(Module = mod, P.Value = p.meta)
  } )
  df3 = do.call(rbind, df3)
  df3 = df3[!is.na(df3$P.Value),]
  df3$FDR = qvalue(c(1,df3$P.Value))$qvalues[-1]

  df3$FDR = p.adjust(df3$P.Value, "fdr")
  df3
}

pickSFT = function( diff.net ){

  df = lapply( names(diff.net), function(key){
    df = lapply( names(diff.net[[key]]), function(key2){
      res = pickSoftThreshold( diff.net[[key]][[key2]] )
      df = res$fitIndices
      df$Cohort = key
      df$Dataset = key2
      df
      })
    do.call(rbind, df)
  })
  do.call(rbind, df)
}








