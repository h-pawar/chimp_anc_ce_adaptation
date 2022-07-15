# for SIV resp/SIV coexp/VIPs results  - of the ancestor
#-----------------------------------------------------------------------------------------------------------------------
# Sun 27 Mar 2022 19:06:15 CEST
# plotting on local machine
require(data.table)
# full data in - https://docs.google.com/spreadsheets/d/1kpV6HTrqqN4WR_ZEjTPYVv4hbvATeHHIUP8HZ1fETPQ/edit#gid=0
# downloaded as a tsv file
siv_tab<-fread("~/Downloads/anc_25mar22 - anc_siv.tsv")

head(siv_tab)
> head(siv_tab)
         Gene_set           Category p-value       FDR p-value
1:                                   0.50000 0.5000000 0.10000 # remove first line = of 0.5,0.5,0.1,0.1,0.05,0.05 : tails categories
2: siv_responsive ['siv_responsive'] 0.16036 0.1747000 0.09252
3:    siv_modules          ['green'] 0.00120 0.0202400 0.12806
4:    siv_modules        ['magenta'] 0.10182 0.3205067 0.00086
5:            vip         merged_IAV 0.00002 0.0003600 0.00010
6:            vip         merged_HIV 0.00046 0.0054100 0.00014
       FDR p-value          FDR
1: 0.10000    0.05         0.05
2: 0.10690 0.04346      0.05754
3: 0.39548 0.27452     0.952372
4: 0.01054 0.02052      0.16857
5: 0.00099 0.00014      0.00112
6: 0.00099 0.00398 0.0185466667

#-----------------------------------------------------------------------------------------------------------------------
# 1) concatenate all enrich_fdr values in R and p.adjust()
# take Fdr_enrich (not BHfdr)
# ie need to correct for number of lineages (i,e. columns)

#ncol(siv_tab)
#[1] 20
#> nrow(siv_tab)
#[1] 12


# extract fdr enrich cols
fdr_enrich<-siv_tab[-c(1),c(seq(4,8,2))]


#-----------------------------------------------------------------------------------------------------------------------
# 1) concatenate all enrich_fdr values in R and p.adjust()
# read in Fdr_enrich values & correct here for number of lineages
# ie need to correct for number of lineages (i,e. columns)

proc_padj<-function(x){
m = as.matrix(fdr_enrich[x,])
q = c()
for (i in seq(1:nrow(m))){
  q = c(q, m[i,])
}
p.adjust(q, method = "BH", n = length(q))
}

# loop over rows
hold.p<-list()
for (i in 1:nrow(fdr_enrich)){
hold.p[[i]]<-proc_padj(i)
}

BH_df<-data.frame(matrix(unlist(hold.p), nrow=length(hold.p), byrow=TRUE),stringsAsFactors=FALSE)
#-----------------------------------------------------------------------------------------------------------------------

# plot this 
 BH_df
 BH_df
          X1        X2         X3
#1  0.1747000 0.1603500 0.16035000
#2  0.0607200 0.5932200 0.95237200
3  0.3205067 0.0316200 0.25285500
4  0.0010800 0.0011200 0.00112000
5  0.0081150 0.0029700 0.01854667
#6  0.1336320 0.1981380 0.27564143
#7  0.1336320 0.2697000         NA
#8  0.1336320 0.2019038 0.19840364
9  0.0859180 0.0058500 0.00084000
10 0.3244472 0.0117450 0.09633500
#11 0.1288770 0.1064160 0.12955000
#12 0.3572475 0.1139900 0.35724750
#13 0.1514267 0.0675600 0.06756000

rownames(BH_df)<-x
colnames(BH_df)<-c(0.5,0.1,0.05)

# should filter again by values < 0.05
# ie remove filter out rows 3,7:11

BH_df<-BH_df[-c(1:2,6:8,11:13),]

# & compare to BH output by JS method
#siv_tab[,c(seq(4,19,3))]

#-----------------------------------------------------------------------------------------------------------------------

# can plot accordingly **
library(pheatmap)

# convert df to matrix & plot
mat_siv_fdr_enrich <- data.matrix(BH_df)

x<-c("SIV responsive",
"SIV coexpression:green",
"SIV coexpression:magenta",
"VIP: merged_IAV",
"VIP: merged_HIV",
"VIP: merged_HCMV",
"VIP: ev71",
"VIP: merged_DENV",
"VIP: merged_EBOV",
"VIP: merged_ADV",
"VIP: merged_HSV",
"VIP: hiv1",
"VIP: asfv")

x1<-x[-c(1:2,6:8,11:13)]

rownames(mat_siv_fdr_enrich)<- x1
colnames(mat_siv_fdr_enrich)<-c("0.5","0.1","0.05")

# plot as heatmap
library(RColorBrewer)
library(grid)
library(gridExtra)

setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(mat_siv_fdr_enrich,
         cluster_cols=FALSE, cluster_rows=FALSE, color = rev(colorRampPalette(rev(brewer.pal(n = 7, name =
  "RdYlBu")))(100)), show_rownames = T, show_colnames=T,legend=T,cellwidth=12, cellheight=12)
setHook("grid.newpage", NULL, "replace")
leg_label <- textGrob("BH FDR",x=0.64,y=0.62,hjust=0,vjust=0,gp=gpar(fontsize=10,fontface="bold"))
grid.draw(leg_label)

x_label <- textGrob("Population",x=0.32,y=0.4,hjust=0,vjust=0,gp=gpar(fontsize=10,fontface="bold"))
grid.draw(x_label)


#-----------------------------------------------------------------------------------------------------------------------
# for GOs  - of the ancestor
#-----------------------------------------------------------------------------------------------------------------------
# Sun 27 Mar 2022 19:06:15 CEST
# plotting on local machine
require(data.table)
# full data in - https://docs.google.com/spreadsheets/d/1kpV6HTrqqN4WR_ZEjTPYVv4hbvATeHHIUP8HZ1fETPQ/edit#gid=0
# downloaded as a tsv file
siv_tab<-fread("~/Downloads/anc_25mar22 - anc_GOs.tsv")

head(siv_tab)

# remove first line = of 0.5,0.5,0.1,0.1,0.05,0.05 : tails categories

#colnames(siv_tab)
#[1] "Type"     "Category" "p-value"  "FDR"      "p-value"  "FDR"     
#[7] "p-value"  "FDR" 

#-----------------------------------------------------------------------------------------------------------------------
# 1) concatenate all enrich_fdr values in R and p.adjust()
# take Fdr_enrich (not BHfdr)
# ie need to correct for number of lineages (i,e. columns)

#ncol(siv_tab)
#[1] 20
#> nrow(siv_tab)
#[1] 12


# extract fdr enrich cols
fdr_enrich<-siv_tab[-c(1),c(seq(4,8,2))]


#-----------------------------------------------------------------------------------------------------------------------
# 1) concatenate all enrich_fdr values in R and p.adjust()
# read in Fdr_enrich
# ie need to correct for number of lineages (i,e. columns)

proc_padj<-function(x){
m = as.matrix(fdr_enrich[x,])
q = c()
for (i in seq(1:nrow(m))){
  q = c(q, m[i,])
}
p.adjust(q, method = "BH", n = length(q))
}

# loop over rows
hold.p<-list()
for (i in 1:nrow(fdr_enrich)){
hold.p[[i]]<-proc_padj(i)
}

BH_df<-data.frame(matrix(unlist(hold.p), nrow=length(hold.p), byrow=TRUE),stringsAsFactors=FALSE)
#-----------------------------------------------------------------------------------------------------------------------

# plot this 
 BH_df
 BH_df
       X1         X2         X3
#1  0.05491667 0.05491667 0.05491667
2  0.01530000 1.00000000         NA
3  0.01530000 0.08711400 0.33307524
4  0.01530000 0.04704750 0.59815171
#5  0.05132727 0.56952164 0.56952164
#6  0.07056000 0.38668080 0.38668080
#7  0.10578000 1.00000000 1.00000000
#8  0.08237501 0.07389600         NA
9  0.01530000 1.00000000         NA
10 0.01530000 0.81388613 1.00000000
#11 0.05948001 1.00000000 1.00000000
#12 0.06255858 0.91892978 0.91892978
13 0.03694800 0.03694800 0.69836277
#14 0.07056000 0.59815171 0.59815171
15 0.01530000 0.04704750 0.37965854
16 0.01530000 0.78240803         NA
17 0.01530000 0.60479308 0.70443158
18 0.01530000 0.04704750 0.55721679
19 0.04172400 0.61185794 0.61185794
#20 0.07056000 0.38668080 0.38668080
#21 0.07520001 1.00000000 1.00000000
#22 0.10314633 1.00000000 1.00000000
#23 0.14529144         NA         NA
#24 0.11361918 0.07389600 0.36528324
#25 0.16840815 0.07389600 0.74798765

# 1,5:8,11:12,14,20:25


x<-c("retromer, cargo-selective complex",
"MCM complex",
"nucleoplasm",
"nucleus",
"I-kappaB/NF-kappaB complex",
"RNA polymerase III complex localization to nucleus",
"Cul3-RING ubiquitin ligase complex",
"protein complex involved in cell adhesion",
"positive regulation by host of viral transcription",
"translation",
"positive regulation of G2/M transition of mitotic cell",
"covalent chromatin modification",
"DNA replication",
"protein transport",
"aldehyde dehydrogenase [NAD(P)+] activity",
"pre-miRNA binding",
"protein binding",
"RNA binding",
"DNA binding",
"RNA polymerase II complex import to nucleus",
"protein complex binding",
"GTPase activity",
"protein tag",
"morphogen activity",
"poly(A) binding")

rownames(BH_df)<-x
colnames(BH_df)<-c(0.5,0.1,0.05)

# should filter again by values < 0.05
# ie remove filter out rows 3,7:11

BH_df<-BH_df[-c(1,5:8,11:12,14,20:25),]

# & compare to BH output by JS method
#siv_tab[,c(seq(4,19,3))]

#-----------------------------------------------------------------------------------------------------------------------

# can plot accordingly **
library(pheatmap)

# convert df to matrix & plot
mat_siv_fdr_enrich <- data.matrix(BH_df)


x1<-x[-c(1,5:8,11:12,14,20:25]

rownames(mat_siv_fdr_enrich)<- x1
colnames(mat_siv_fdr_enrich)<-c("0.5","0.1","0.05")

# plot as heatmap
library(RColorBrewer)
library(grid)
library(gridExtra)

setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(mat_siv_fdr_enrich,
         cluster_cols=FALSE, cluster_rows=FALSE, color = rev(colorRampPalette(rev(brewer.pal(n = 7, name =
  "RdYlBu")))(100)), show_rownames = T, show_colnames=T,legend=T,cellwidth=12, cellheight=12)
setHook("grid.newpage", NULL, "replace")
leg_label <- textGrob("BH FDR",x=0.64,y=0.62,hjust=0,vjust=0,gp=gpar(fontsize=10,fontface="bold"))
grid.draw(leg_label)

x_label <- textGrob("Population",x=0.32,y=0.4,hjust=0,vjust=0,gp=gpar(fontsize=10,fontface="bold"))
grid.draw(x_label)
