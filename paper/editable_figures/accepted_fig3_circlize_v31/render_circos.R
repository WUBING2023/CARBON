# Adapted from the user's brother's plot_circos.R. Original is preserved verbatim
# in reference_bundle. Frozen CARBON linkage, leaf order, PPMs and annotations
# replace the reference study's selected motifs and biological lists.
suppressPackageStartupMessages({
  .libPaths(c(normalizePath('Rlibrary'),.libPaths()))
  library(circlize); library(grid); library(ggplot2); library(ggseqlogo)
})
args <- commandArgs(trailingOnly=TRUE)
width_mm <- if(length(args)) as.numeric(args[1]) else 285
stem <- if(length(args)>1) args[2] else 'circos_carbon_285mm'
meta <- read.csv('inputs/leaf_metadata.csv',stringsAsFactors=FALSE,check.names=FALSE)
inp <- read.csv('inputs/tree_input_ids.csv',stringsAsFactors=FALSE)
m <- read.csv('inputs/tree_merge.csv'); n <- nrow(meta)
hc <- structure(list(merge=as.matrix(m[,1:2]),height=m$height,
                     order=order(inp$order),labels=inp$task_id,
                     method='average',call=quote(frozen_CARBON_average_linkage),
                     dist.method='euclidean'),class='hclust')
storage.mode(hc$merge) <- 'integer'
dend <- as.dendrogram(hc)
stopifnot(n==219,identical(unname(labels(dend)),meta$task_id))
ppms <- read.csv('inputs/motif_ppms.csv',stringsAsFactors=FALSE)
function_cols <- c('Splicing'='#E6211A','Stability/decay'='#19499C',
                   'Translation'='#6BB72D','Stability + translation'='#F4A016',
                   'No category match'='#BDBDBD')
category <- ifelse(meta$splicing=='Yes','Splicing',
             ifelse(meta$stability_decay=='Yes' & meta$translation=='Yes','Stability + translation',
             ifelse(meta$stability_decay=='Yes','Stability/decay',
             ifelse(meta$translation=='Yes','Translation','No category match'))))
names(category) <- meta$task_id
colorize_branches <- function(node) {
  if(is.leaf(node)) {
    attr(node,'edgePar') <- list(col=unname(function_cols[category[attr(node,'label')]]),lwd=.65*width_mm/285)
  } else {
    attr(node,'edgePar') <- list(col='#222222',lwd=.52*width_mm/285)
    for(i in seq_along(node)) node[[i]] <- colorize_branches(node[[i]])
  }
  node
}
dend <- colorize_branches(dend)
groups <- cutree(hc,k=6)[meta$task_id]
cluster_colors <- c('#F9F9D9','#E0ECD3','#DBF0F9','#F7BBA9','#FFE0CC','#EAEAFF')
ring_fields <- c('splicing','stability_decay','translation','RRM','ZnF','KH')
ring_positive <- c('Yes','Yes','Yes','RRM','ZnF','KH')
ring_colors <- c('#E6211A','#19499C','#6BB72D','#52C9EC','#E5D74A','#8BC01C')
base_colors <- c(A='#6BB72D',C='#19499C',G='#F4A016',U='#E6211A')
custom_col_scheme_seq <- make_col_scheme(chars=names(base_colors),groups=names(base_colors),cols=unname(base_colors))

pdf(paste0(stem,'.pdf'),width=width_mm/25.4,height=width_mm/25.4,
    pointsize=12,paper='special',family='Helvetica',useDingbats=FALSE,compress=TRUE)
par(mar=c(0,0,0,0),xaxs='i',yaxs='i')
circos.clear()
circos.par(cell.padding=c(0,0,0,0),start.degree=90,gap.degree=0,
           clock.wise=TRUE,points.overflow.warning=FALSE,
           canvas.xlim=c(-1.025,1.025),canvas.ylim=c(-1.025,1.025))
circos.initialize(sectors='all',xlim=c(0,n))
# Outside-in allocation: same six frozen annotation tracks, separate no-match
# (grey) from missing/unknown (white). No new brother-study annotations imported.
for(k in rev(seq_along(ring_fields))) {
  circos.trackPlotRegion(ylim=c(0,1),bg.border=NA,track.height=.015,
                        track.margin=c(.0035,.0035),panel.fun=function(x,y){
    vals <- meta[[ring_fields[k]]]
    for(i in seq_len(n)) {
      col <- if(vals[i]==ring_positive[k])ring_colors[k] else if(vals[i]=='Not') '#EEEEEE' else '#FFFFFF'
      circos.rect(i-1,0,i,1,col=col,border=NA)
    }
  })
}
circos.trackPlotRegion(ylim=c(0,1),bg.border=NA,track.height=.20,
                      track.margin=c(0,.007),panel.fun=function(x,y){
  circos.rect(0,0,n,1,col='#FAF8F1',border=NA)
})
logo_track <- get.current.track.index()
label_pt <- numeric(n)
circos.trackPlotRegion(ylim=c(0,1),bg.border=NA,track.height=.16,
                      track.margin=c(0,.007),panel.fun=function(x,y){
  # Fit the longest full RBP-cell names within the radial label track.
  track_width_in <- .155*width_mm/2.05/25.4
  for(i in seq_len(n)) {
    font_pt <- 7.5*width_mm/285
    available <- strwidth(meta$label[i],units='inches',cex=font_pt/12)
    if(available>track_width_in)font_pt <- font_pt*track_width_in/available
    label_pt[i] <<- font_pt
    circos.text(i-.5,0,meta$label[i],facing='clockwise',niceFacing=TRUE,
                cex=font_pt/12,adj=c(0,.5),col='#222222')
  }
})
max_height <- attr(dend,'height')
circos.trackPlotRegion(ylim=c(0,max_height),bg.border=NA,track.height=.40,
                      track.margin=c(.002,0),panel.fun=function(x,y){
  for(i in seq_len(n))
    circos.rect(i-1,0,i,max_height,col=cluster_colors[groups[i]],border=NA)
  circos.dendrogram(dend,max_height=max_height,facing='outside')
})
geometry <- vector('list',n)
# Original code's ggseqlogo/grid rotated-viewport layer, now using four sequence
# channels only. No P/U structure logo/band is fabricated for CARBON.
for(i in seq_len(n)) {
  ppm <- as.matrix(ppms[ppms$task_id==meta$task_id[i],c('A','C','G','U')])
  ppm <- ppm/rowSums(ppm);seq <- t(ppm);rownames(seq) <- c('A','C','G','U')
  plot_seq <- ggseqlogo(seq,method='bits',seq_type='rna',col_scheme=custom_col_scheme_seq)+
    scale_y_continuous(limits=c(0,2),expand=c(0,0))+
    scale_x_continuous(limits=c(.5,9.5),expand=c(0,0))+
    theme_void()+theme(plot.margin=margin(0,0,0,0))
  polar <- circlize(i-.5,.5,track.index=logo_track)
  theta <- as.numeric(polar[1,1]);radius <- as.numeric(polar[1,2])
  cart <- radius*c(cos(theta*pi/180),sin(theta*pi/180))
  x_npc <- grconvertX(cart[1],from='user',to='ndc')
  y_npc <- grconvertY(cart[2],from='user',to='ndc')
  degree <- theta%%360
  if(degree>90 && degree<270)degree <- degree-180
  radial_mm <- .196*width_mm/2.05
  tangent_mm <- .84*radius*2*pi/n*width_mm/2.05
  pushViewport(viewport(x=x_npc,y=y_npc,width=unit(radial_mm,'mm'),
                        height=unit(tangent_mm,'mm'),angle=degree))
  grid.draw(ggplotGrob(plot_seq));popViewport()
  geometry[[i]] <- data.frame(task_id=meta$task_id[i],label=meta$label[i],leaf=i,
                              theta=theta,radius=radius,x_npc=x_npc,y_npc=y_npc,
                              rotation=degree,label_pt=label_pt[i],cluster=groups[i])
}
dev.off();circos.clear()
write.csv(do.call(rbind,geometry),paste0(stem,'_geometry.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),'R_SESSION_RENDER.txt')
cat('Rendered',n,'frozen CARBON motifs with adapted circlize/ggseqlogo code.\n')
