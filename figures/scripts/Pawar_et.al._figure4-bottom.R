# data
#---- Notes ----#
# pryr is used because it is a way of wrapping a call to plot as a function.
# helps with complicated plot layouts.
# -log10(p) < 1.

rm(list=ls())

library(data.table)
library(pryr)

setwd("~/OneDrive - University College London/Projects/Internal_branch_enrichment/Pawar_et.al._code-HO/")

cd4regionGenes <- readRDS("input/cd4regionGenes.rds")
cd4.pclr <- readRDS("input/cd4.windows.rds")
windows <- readRDS("input/windows.rds")
cd4.dafs <- readRDS("input/cd4.dafs.rds")
cd4.pbsnj  <- readRDS("input/cd4.pbsnj.rds")
nonsyn <- readRDS("input/cd4.nonsyn.rds")


setkey(cd4regionGenes,chrom,start,stop)
setkey(cd4.pclr,Chr,Physpos)
setkey(cd4.dafs,chr,start,end)
setkey(cd4.pbsnj,chr,start,end)


# updating gene tracks...
ens <- fread("ens_genes_chr12_6799607-7332026.txt",select = c("name","chrom","strand","txStart","txEnd","name2"))
ens_2_name <- fread("ens_2_gene-name.txt")
ens <- ens_2_name[ens,on=.(name)]
ens <- ens[!value %in% NA]
ens[,strand:=ifelse(strand=="+",1,-1)]
ens <- ens[,.(chrom,start=txStart,stop=txEnd,gene=value,score=0,strand)]
#ens <- ens[gene %in% genes$gene] # This doesn't work, I don't know what genes is supposed to be

xmin <- 6950000
xmax <- 6987500

ens[,new.start:=ifelse(start < xmin, xmin, start)]
ens[,new.stop:=ifelse(stop > xmax, xmax, stop)]
# remove some gene names for plotting clarity - list them in the figure legend or methods?
#thin <- c("ING4","ZNF384","PIANP","PTMS","LAG3","GPR162","P3H3","GNB3","TPI1","C12orf57","LRRC23")
thin <- c("MLF2")
# plot without p-value cutoffs; annotated AA
# 3p-clr plot

# make the gene layout plot
# xmin <- cd4.dafs[,min(start)]
# xmax <- cd4.dafs[,max(end)]

#### HO ####
ens$colour <- "grey"
ens[gene=="CD4"]$colour <- "deepskyblue"

## Remove everything which lies outside the plot to prevent plots from hanging off the edge of the axis
ens=ens[ens$start<=xmax & ens$stop>=xmin,]
cd4.pclr=cd4.pclr[cd4.pclr$PhysStart<=xmax & cd4.pclr$PhysEnd>=xmin,]
cd4.dafs=cd4.dafs[cd4.dafs$start<=xmax & cd4.dafs$end>=xmin,]
cd4.pbsnj=cd4.pbsnj[cd4.pbsnj$start<=xmax & cd4.pbsnj$end>=xmin,]

# Plot genes
genes.pryr %<a-% {
  plot(0, type="n", xlab="", ylab="", xlim=c(xmin, xmax), ylim=c(-5, 5),axes=FALSE)
  abline(h = 0)
  for (g in ens[,gene]){
    if(ens[gene==g]$strand==1){
      rect(ens[gene==g]$new.start, 0, ens[gene==g]$new.stop, 3,
           col=ens[gene==g]$colour, border=par("fg"), lty=NULL, lwd=par("lwd"), xpd=FALSE)
      if(!g %in% thin){
        text(x=ens[gene==g]$new.start + ((ens[gene==g]$new.stop-ens[gene==g]$new.start)/2),
             y=4.25,
             labels="", cex=2)
      }
    }
    if(ens[gene==g]$strand== -1){
      rect(ens[gene==g]$new.start, 0, ens[gene==g]$new.stop, -3,
           col="grey", border=par("fg"), lty=NULL, lwd=par("lwd"), xpd=FALSE)
      if(!g %in% thin){
        text(x=ens[gene==g]$new.start + ((ens[gene==g]$new.stop-ens[gene==g]$new.start)/2),
           y= -4.25,
           labels = "", cex=2)
      }
    }
  }
  # #patch in two labels
  # text(x=cd4.pclr[,min(Physpos)]-3e3,
  #      y= -4.5,
  #      labels="ACRBP")
  # text(x=cd4.pclr[,max(Physpos)]+4e3,
  #      y= -4.5,
  #      labels="C1R")
}
# make the 3p-clr scores plot
pclr.pryr %<a-% {
  plot(windows[Chr==12 & Physpos >= xmin & Physpos <= xmax][order(Physpos)]$Physpos,
       -log10(windows[Chr==12 & Physpos >= xmin & Physpos <= xmax][order(Physpos)]$anc.p),
       type='l',
       pch=19,
       col='brown',
       lwd=2,
       xaxt = "n",
       yaxt = "n",
       ylim=c(0,3),
       xlim=c(xmin,xmax),
       ylab = "",
       cex.lab=2,
       bty="n"
  )
  # add axis on y
  ytick<-seq(0, 5, by=2.5)
  axis(2, at = ytick,las=1,cex=2, cex.axis=1.5)
  legend("topright",
         title=NULL,
         legend="3P-CLR",
         col="brown",
         pt.bg = 'white',
         bty="n",
         cex=1.5)
  
}

clade.DAF.pryr %<a-% {
  cd4.dafs[order(start)][cladeDiffP < 0.3][,plot(start,-log10(cladeDiffP),type='l',col='brown',ylim=c(0,3),xaxt = "n",
                                                 yaxt = "n",
                                                 ylab = "",
                                                 cex.lab=2,
                                                 lwd=2,
                                                 bty="n",
                                                 pch=19,
                                                 xlim=c(xmin,xmax))]
  legend("topright",
         title=NULL,
         legend="Clade DAF",
         col="brown",
         bty="n",
         cex=1.5)
  ytick<-seq(0, 5, by=2.5)
  axis(2, at = ytick,las=1,cex=1.5, cex.axis=1.5)
  xtick<-seq(round(xmin/1e6,1), round(xmax/1e6,1), by=0.05)
  #axis(1, at = xtick,cex=1)
  axis(side=1,at=xtick,tcl=0.4,lwd.ticks=3,mgp=c(0,0.5,0),outer=F,pos = 2, cex.axis=1.5)
  #c
}

# central pbsnj
centralpbsnj.pryr %<a-% {
  cd4.pbsnj[order(start)][c.p < 0.3][,plot(start,-log10(c.p),type='l',col='green3',ylim=c(0,5),xaxt = "n",
                                           yaxt = "n",
                                           ylab = "",
                                           cex.lab=2,
                                           lwd=2,
                                           bty="n",
                                           pch=19,
                                           xlim=c(xmin,xmax))]
  legend("topright",
         title=NULL,
         legend="\nCentral PBSnj",
         col="green3",
         bty="n",
         cex=1.5)
  ytick<-seq(0, 5, by=2.5)
  axis(2, at = ytick,las=1,cex=1.5, cex.axis=1.5)
}
# eastern pbsnj
easternpbsnj.pryr %<a-% {
  cd4.pbsnj[order(start)][e.p < 0.3][,plot(start,-log10(e.p),type='l',col='darkorange',ylim=c(0,3), xlab="Position on Chr12",
                                           yaxt = "n",
                                           ylab = "",
                                           cex.lab=2,
                                           lwd=2,
                                           bty="n",
                                           pch=19,
                                           cex.axis=1.5,    ###### This controls the size of the numbers on the x axis 
                                           xlim=c(xmin,xmax))]
  legend("topright",
         title=NULL,
         legend="Eastern PBSnj",
         col="darkorange",
         bty="n",
         cex=1.5)
  ytick<-seq(0, 5, by=2.5)
  axis(2, at = ytick,las=1,cex=1.5, cex.axis=1.5)
  axis(1, at = c(0,10000000),labels=c("",""),lwd.ticks=0, cex.axis=1.5, cex=1.5)
}

# internal pbsnj
intpbsnj.pryr %<a-% {
  cd4.pbsnj[order(start)][i.p < 0.3][,plot(start,-log10(i.p),type='l',col='brown',ylim=c(0,3),xaxt = "n",
                                           yaxt = "n",
                                           ylab = "",
                                           cex.lab=2,
                                           lwd=2,
                                           pch=19,
                                           bty="n",
                                           xlim=c(xmin,xmax))]
  legend("topright",
         title=NULL,
         legend="Internal branch\n              PBSnj",
         col="brown",
         bty="n",
         cex=1.5)
  ytick<-seq(0, 5, by=2.5)
  axis(2, at = ytick,las=1,cex=1.5, cex.axis=1.5)
}
pdf(file="output/figure4_subplot-bottom.pdf", pointsize=5, width=6, height=3.3)
par(oma=c(3,1,1,1),mar = c(1.5,4.5,1,0))
layout(matrix(c(1,1,2,2,3,3,4,4,5,5,5,6,6), 13, 1, byrow = TRUE))
genes.pryr
pclr.pryr
intpbsnj.pryr
clade.DAF.pryr
centralpbsnj.pryr
easternpbsnj.pryr
mtext("-log10(p-value)", outer = TRUE, cex = 1.5, side = 2,line = -1)
par(xpd = NA)
# Add lines
#abline(v=nonsyn$start[1], lwd=3, col='blue', lty=3) # central SNP - decided to not include as it only has a p-value of 0.014
abline(v=nonsyn$start[2], lwd=3, col='black', lty=3) # central SNP
abline(v=nonsyn$start[3], lwd=3, col='black', lty=3) # T93P in ancestor

abline(v=6963043, lwd=3, col='magenta', lty=3) # splice region variant in ancestor 

mtext("  Chromosome 12", outer = TRUE, cex = 1.3,side = 1,line = 1.5)
mtext("CD4", outer = TRUE, cex = 1.5, side = 3,line = -1.5)
dev.off()

