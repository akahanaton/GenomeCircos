library(RCircos)
library(plotrix)
library(stringr)
library(data.table)
#--------------------------------------------------
# library(GenomicFeatures)
#--------------------------------------------------

# load human information
data(UCSC.HG19.Human.CytoBandIdeogram);

chr.length = read.table("p_ctg.rename.fa.fai",header=F,sep="\t",stringsAsFactors=F)
gene.info = read.table("./maker.gff3.MAKER.rename.gene",header=F,sep="\t",stringsAsFactors=F)

curChr="000000F"
curPdf = paste(curChr, ".pdf", sep="")
curChr.len = chr.length[chr.length[,1]== curChr,2]
curChr.gene = gene.info[gene.info[,1]==curChr,]
curChr.gene = cbind(curChr.gene,gsub("[(ID=)(:)]", "", str_match(curChr.gene[,9],"ID=.*?:")))
curChr.gene = curChr.gene[!is.na(curChr.gene[,10]),]
curChr.gene = curChr.gene[order(curChr.gene$V4),]
curChr.gene = curChr.gene[ ! curChr.gene[,4] > curChr.gene[,5],]
curChr.gene[,1] = sub("F","",paste('chr',curChr.gene[,1],sep=""))

curChr.cyto.pos <- curChr.gene[curChr.gene[,7] == '+',c(1,4,5,6,7,10)]
colnames(curChr.cyto.pos) = colnames(UCSC.HG19.Human.CytoBandIdeogram)
colnames(curChr.cyto.pos)[6] = "Gene"
curChr.cyto.pos[,"Band"] = UCSC.HG19.Human.CytoBandIdeogram[1:dim(curChr.cyto.pos)[1],"Band"]
curChr.cyto.pos[,"Stain"] = UCSC.HG19.Human.CytoBandIdeogram[1:dim(curChr.cyto.pos)[1],"Stain"]
# for cyto plot:
# add a fake gene start from the chromosome base 0
curChr.cyto.pos = rbind(curChr.cyto.pos[1,], curChr.cyto.pos)
curChr.cyto.pos[1,2] = 0
curChr.cyto.pos[1,3] = curChr.cyto.pos[2,2] - 1
curChr.cyto.pos[1,6] = NA
# add a fake gene end at the last chromosome base
if(curChr.cyto.pos[nrow(curChr.cyto.pos), 3] < curChr.len){
    curChr.cyto.pos = rbind(curChr.cyto.pos, curChr.cyto.pos[1,])
    curChr.cyto.pos[nrow(curChr.cyto.pos),2] = curChr.cyto.pos[nrow(curChr.cyto.pos)-1,3] + 1
    curChr.cyto.pos[nrow(curChr.cyto.pos),3] = curChr.len
}
wrongGene = c()
# no overlap gene
for(i in 2:(dim(curChr.cyto.pos)[1])){
    if(curChr.cyto.pos[i,2] <= curChr.cyto.pos[i-1,3]){
        wrongGene = c(wrongGene,i)
    }
}
#--------------------------------------------------
# print(curChr.cyto.pos[wrongGene,])
#--------------------------------------------------
curChr.cyto.pos = curChr.cyto.pos[-wrongGene,]
missing.row = curChr.cyto.pos[1,]
for(i in 2:(nrow(curChr.cyto.pos)-1)){
    if(curChr.cyto.pos[i,2] > (curChr.cyto.pos[i-1,3]+1)){
        newRow = curChr.cyto.pos[i,]
        newRow$ChromStart = curChr.cyto.pos[i-1,3] + 1
        newRow$ChromEnd = curChr.cyto.pos[i,2] - 1
        newRow$Gene = i
        missing.row = rbind(missing.row,newRow)
    }
}
curChr.cyto.pos = rbind(curChr.cyto.pos, missing.row[-1,])
curChr.cyto.pos = curChr.cyto.pos[order(curChr.cyto.pos$ChromStart),]
curChr.cyto.pos[is.na(curChr.cyto.pos$Gene),"Stain"] = "gneg"
curChr.cyto.pos[!is.na(curChr.cyto.pos$Gene),"Stain"] = "gpos25"

curChr.cyto.neg <- curChr.gene[curChr.gene[,7] == '-',c(1,4,5,6,7)]
colnames(curChr.cyto.neg) = colnames(UCSC.HG19.Human.CytoBandIdeogram)
curChr.cyto.neg[,"Band"] = UCSC.HG19.Human.CytoBandIdeogram[1:dim(curChr.cyto.neg)[1],"Band"]
curChr.cyto.neg[,"Stain"] = UCSC.HG19.Human.CytoBandIdeogram[1:dim(curChr.cyto.neg)[1],"Stain"]
# for cyto plot:
# add a fake gene start from the chromosome base 0
curChr.cyto.neg = rbind(curChr.cyto.neg[1,], curChr.cyto.neg)
curChr.cyto.neg[1,2] = 0
curChr.cyto.neg[1,3] = curChr.cyto.neg[2,2] - 1
# add a fake gene end at the last chromosome base
if(curChr.cyto.neg[nrow(curChr.cyto.neg), 3] < curChr.len){
    curChr.cyto.neg = rbind(curChr.cyto.neg, curChr.cyto.neg[1,])
    curChr.cyto.neg[nrow(curChr.cyto.neg),2] = curChr.cyto.neg[nrow(curChr.cyto.neg)-1,3] + 1
    curChr.cyto.neg[nrow(curChr.cyto.neg),3] = curChr.len
}
wrongGene = c()
# no overlap gene
for(i in 2:(dim(curChr.cyto.neg)[1])){
    if(curChr.cyto.neg[i,2] <= curChr.cyto.neg[i-1,3]){
        wrongGene = c(wrongGene,i)
    }
}
#--------------------------------------------------
# print(curChr.cyto.neg[wrongGene,])
#--------------------------------------------------
curChr.cyto.neg = curChr.cyto.neg[-wrongGene,]

all.gene <- curChr.gene[,c(1,4,5,6,7)]
colnames(all.gene) = colnames(UCSC.HG19.Human.CytoBandIdeogram)
all.gene[,"Band"] = UCSC.HG19.Human.CytoBandIdeogram[1:dim(all.gene)[1],"Band"]
all.gene[,"Stain"] = UCSC.HG19.Human.CytoBandIdeogram[1:dim(all.gene)[1],"Stain"]
# for cyto plot:
wrongGene = c()
# no overlap gene
for(i in 2:(dim(all.gene)[1])){
    if(all.gene[i,2] <= all.gene[i-1,3]){
        wrongGene = c(wrongGene,i)
    }
}
#--------------------------------------------------
all.gene = all.gene[-wrongGene,]

#-------------------------------------------------- 

pdf(file=curPdf, height=15, width=15, compress=TRUE);
#--------------------------------------------------
# par(mai=c(0.25, 0.25, 0.25, 0.25));
# plot.new();
# plot.window(c(-2.5,750), c(-2.5, 100));
# rect(0,0,512,20)
# for(i in 1:512){
#     abline(v=i,col=ColorRamp[i])
# }
#--------------------------------------------------

#--------------------------------------------------
#--------------------------------------------------
chr.exclude <- NULL
tracks.inside <- 2;
tracks.outside <- 40;
RCircos.Set.Core.Components(curChr.cyto.pos, chr.exclude, tracks.inside, tracks.outside);
rcircos.params <- RCircos.Get.Plot.Parameters();
rcircos.params$radius.len <- 3;
rcircos.params$track.height <- 0.2;
rcircos.params$base.per.unit <- 50000;
rcircos.params$heatmap.width <- 1;
rcircos.params$hist.width <- 1;
rcircos.params$chrom.paddings <- 10;
rcircos.params$track.paddings <- 1;
RCircos.Reset.Plot.Parameters(rcircos.params);
RCircos.List.Parameters()
#--------------------------------------------------
# source("./RCircos.Set.Cytoband.data.R")
# my.RCircos.Set.Cytoband.data(curChr.cyto.pos);
#--------------------------------------------------
#--------------------------------------------------
# source("./RCircos.Set.Base.Plot.Positions.R")
# my.RCircos.Set.Base.Plot.Positions();
#--------------------------------------------------
print("Starting plot")
RCircos.Set.Plot.Area();

RCircos.Cyto <- RCircos.Get.Plot.Ideogram();
RCircos.Par <- RCircos.Get.Plot.Parameters();
RCircos.Pos <- RCircos.Get.Plot.Positions();

#--------------------------------------------------
# source("./RCircos.Chromosome.Ideogram.Plot.R")
#--------------------------------------------------
RCircos.Chromosome.Ideogram.Plot();
print("Done plot chromosome ideogram")

#--------------------------------------------------
# RCircos.Set.Cytoband.data(curChr.cyto.neg);
# RCircos.Set.Base.Plot.Positions();
# RCircos.Chromosome.Ideogram.Plot();
#--------------------------------------------------
#--------------------------------------------------
# RCircos.Set.Core.Components(curChr.cyto.neg, chr.exclude, tracks.inside, tracks.outside);
#--------------------------------------------------

curChr.cyto.pos.name = curChr.cyto.pos[!is.na(curChr.cyto.pos$Gene),]
side <- "in";
name.col <- 6;
track.num <- 1;
RCircos.Gene.Connector.Plot(curChr.cyto.pos.name, track.num, side);
track.num <- 2;
RCircos.Gene.Name.Plot(curChr.cyto.pos.name, name.col,track.num, side);
print("Done plot gene names and linkes")

#--------------------------------------------------
# source("./RCircos.Coverage.Plot.R")
# source("./RCircos.Variants.Plot.R")
#--------------------------------------------------
# track.begin <- 3;
# track.end <- 4;
# data.col = 10
# side="out"
# y.lim = 10
# norm.index = 20
# RCircos.Coverage.Plot(ipd.lung.m6A, data.col, track.begin, track.end, side, norm.index, y.lim,  "Genome coverage of depth, enriched");
# track.num <- 4;
# data.col = 5
# side="in"
# RCircos.Variants.Plot(cov.info1, data.col, track.num, side);
#--------------------------------------------------

#--------------------------------------------------
# cov.info2 = read.csv("./Lib.tsv",header=T,sep="\t",stringsAsFactors=F)
# cov.info2 = cov.info2[-1, ]
# track.begin <- 0;
# track.end <- 4;
# data.col = 5
# side="out"
# y.lim = 30000
# norm.index = 147848 +  143687
# RCircos.Coverage.Plot(cov.info2, data.col, track.begin, track.end, side, norm.index, y.lim,  "Genome coverage of depth, unenriched");
# track.num <- 4;
# data.col = 5
# side="in"
# RCircos.Variants.Plot(cov.info2, data.col, track.num, side);
#-------------------------------------------------- 

#--------------------------------------------------
# data(RCircos.Heatmap.Data);
# data.col <- 6;
# track.num <- 5;
# side <- "in";
# RCircos.Heatmap.Plot(RCircos.Heatmap.Data, data.col, track.num, side);
#-------------------------------------------------- 

#--------------------------------------------------
# data(RCircos.Scatter.Data);
# data.col <- 5;
# track.num <- 6;
# side <- "in";
# by.fold <- 1;
# RCircos.Scatter.Plot(RCircos.Scatter.Data, data.col, track.num, side, by.fold);
#-------------------------------------------------- 

#--------------------------------------------------
# data(RCircos.Line.Data);
#--------------------------------------------------
#--------------------------------------------------
# data.col <- 10;
# track.num <- 7;
# side <- "out";
# RCircos.Line.Plot(ipd.kidney.m6A, data.col, track.num, side);
#-------------------------------------------------- 

regionPad = 1
chromRangeNum = round(curChr.len / rcircos.params$base.per.unit / regionPad) + 1
chromRangeTable = data.table(chr=rep("chr000000", chromRangeNum), start=seq(0, chromRangeNum -1 ) * rcircos.params$base.per.unit * regionPad, end=seq(1, chromRangeNum ) * rcircos.params$base.per.unit * regionPad -1, Name=paste("region",1:chromRangeNum,sep=""), key = c("chr","start","end"))
chromRangeFrame = data.frame(chromRangeTable)

ipd.kidney = read.table("./kidney.000000F.gff.m6c",header=F,sep="\t",stringsAsFactors=F)
ipd.kidney = cbind(ipd.kidney, gsub("[(IPDRatio=)(;)]", "", str_match(ipd.kidney[,9],"IPDRatio=.*?;")),gsub("[(identificationQv=)]", "", str_match(ipd.kidney[,9],"identificationQv=.*")))
colnames(ipd.kidney) = c("Chromosome","source","type","chromStart","chromEnd","score","strand","note","des","Data","Qv")
ipd.kidney[,1] = sub("F","",paste('chr',ipd.kidney[,1],sep=""))
ipd.kidney = ipd.kidney[order(ipd.kidney$chromStart),]
ipd.kidney.m6A = ipd.kidney[ as.character(ipd.kidney[,3]) == "m6A" & as.numeric(ipd.kidney[,11]) >=10, c(1,4,5,10)]
ipd.kidney.m4C = ipd.kidney[ as.character(ipd.kidney[,3]) == "m4C" & as.numeric(ipd.kidney[,11]) >=10, c(1,4,5,10)]
ipd.kidney.m6A.table = data.table(chr=ipd.kidney.m6A[,1], start=ipd.kidney.m6A[,2],end=ipd.kidney.m6A[,3],IPDRatio=ipd.kidney.m6A[,4],key = c("chr","start","end"))
ipd.kidney.m6A.res = foverlaps(ipd.kidney.m6A.table, chromRangeTable, nomatch=0)
count = data.frame(table(ipd.kidney.m6A.res$start))
colnames(count) = c("start","count")
ipd.kidney.m6A.res = merge(chromRangeFrame, count, by="start")
ipd.kidney.m6A.plot = ipd.kidney.m6A.res[,c(2,1,3,5)]
ipd.kidney.m6A.plot[,4] = ipd.kidney.m6A.plot[,4] / max(ipd.kidney.m6A.plot[,4])
ipd.kidney.m6A.plot = cbind(ipd.kidney.m6A.plot, rep("red", nrow(ipd.kidney.m6A.plot)))
ipd.kidney.m4C.table = data.table(chr=ipd.kidney.m4C[,1], start=ipd.kidney.m4C[,2],end=ipd.kidney.m4C[,3],IPDRatio=ipd.kidney.m4C[,4],key = c("chr","start","end"))
ipd.kidney.m4C.res = foverlaps(ipd.kidney.m4C.table, chromRangeTable, nomatch=0)
count = data.frame(table(ipd.kidney.m4C.res$start))
colnames(count) = c("start","count")
ipd.kidney.m4C.res = merge(chromRangeFrame, count, by="start")
ipd.kidney.m4C.plot = ipd.kidney.m4C.res[,c(2,1,3,5)]
ipd.kidney.m4C.plot[,4] = ipd.kidney.m4C.plot[,4] / max(ipd.kidney.m4C.plot[,4])
ipd.kidney.m4C.plot = cbind(ipd.kidney.m4C.plot, rep("red", nrow(ipd.kidney.m4C.plot)))

ipd.lung = read.table("./lung.000000F.gff.m6c",header=F,sep="\t",stringsAsFactors=F)
ipd.lung = cbind(ipd.lung, gsub("[(IPDRatio=)(;)]", "", str_match(ipd.lung[,9],"IPDRatio=.*?;")),gsub("[(identificationQv=)]", "", str_match(ipd.lung[,9],"identificationQv=.*")))
colnames(ipd.lung) = c("Chromosome","source","type","chromStart","chromEnd","score","strand","note","des","Data","Qv")
ipd.lung[,1] = sub("F","",paste('chr',ipd.lung[,1],sep=""))
ipd.lung = ipd.lung[order(ipd.lung$chromStart),]
ipd.lung.m6A = ipd.lung[ as.character(ipd.lung[,3]) == "m6A" & as.numeric(ipd.lung[,11])>=10, c(1,4,5,10)]
ipd.lung.m4C = ipd.lung[ as.character(ipd.lung[,3]) == "m4C" & as.numeric(ipd.lung[,11]) >=10, c(1,4,5,10)]

ipd.lung.m6A.table = data.table(chr=ipd.lung.m6A[,1], start=ipd.lung.m6A[,2],end=ipd.lung.m6A[,3],IPDRatio=ipd.lung.m6A[,4],key = c("chr","start","end"))
ipd.lung.m6A.res = foverlaps(ipd.lung.m6A.table, chromRangeTable, nomatch=0)
count = data.frame(table(ipd.lung.m6A.res$start))
colnames(count) = c("start","count")
ipd.lung.m6A.res = merge(chromRangeFrame, count, by="start")
ipd.lung.m6A.plot = ipd.lung.m6A.res[,c(2,1,3,5)]
ipd.lung.m6A.plot[,4] = ipd.lung.m6A.plot[,4] / max(ipd.lung.m6A.plot[,4])
ipd.lung.m6A.plot = cbind(ipd.lung.m6A.plot, rep("red", nrow(ipd.lung.m6A.plot)))

ipd.lung.m4C.table = data.table(chr=ipd.lung.m4C[,1], start=ipd.lung.m4C[,2],end=ipd.lung.m4C[,3],IPDRatio=ipd.lung.m4C[,4],key = c("chr","start","end"))
ipd.lung.m4C.res = foverlaps(ipd.lung.m4C.table, chromRangeTable, nomatch=0)
count = data.frame(table(ipd.lung.m4C.res$start))
colnames(count) = c("start","count")
ipd.lung.m4C.res = merge(chromRangeFrame, count, by="start")
ipd.lung.m4C.plot = ipd.lung.m4C.res[,c(2,1,3,5)]
ipd.lung.m4C.plot[,4] = ipd.lung.m4C.plot[,4] / max(ipd.lung.m4C.plot[,4])
ipd.lung.m4C.plot = cbind(ipd.lung.m4C.plot, rep("red", nrow(ipd.lung.m4C.plot)))

# plot m6A and m4C density
data(RCircos.Histogram.Data);
source("RCircos.Histogram.Plot.R")
data.col <- 4;
side <- "out";
track.num <- 1;
colnames(ipd.kidney.m6A.plot) = c(colnames(RCircos.Histogram.Data),"PlotColor")
my.RCircos.Histogram.Plot(ipd.kidney.m6A.plot, data.col, track.num, side);
track.num <- 2;
colnames(ipd.lung.m6A.plot) = c(colnames(RCircos.Histogram.Data),"PlotColor")
my.RCircos.Histogram.Plot(ipd.lung.m6A.plot, data.col, track.num, side);
#--------------------------------------------------
# track.num <- 3;
# colnames(ipd.kidney.m4C.plot) = c(colnames(RCircos.Histogram.Data),"PlotColor")
# my.RCircos.Histogram.Plot(ipd.kidney.m4C.plot, data.col, track.num, side);
# track.num <- 4;
# colnames(ipd.lung.m4C.plot) = c(colnames(RCircos.Histogram.Data),"PlotColor")
# my.RCircos.Histogram.Plot(ipd.lung.m4C.plot, data.col, track.num, side);
#--------------------------------------------------
print("Done plot m6A m4C density")

source("./RCircos.Chromosome.Region.Plot.R")
data.col <- 4;
side <- "out";
track.num <- 5;
m5c.kidney = read.table("./kidney.000000F.gff.m5c",header=F,sep="\t",stringsAsFactors=F)[,c(1,4,5,9)]
m5c.kidney[,1] = sub("F","",paste('chr',m5c.kidney[,1],sep=""))
colnames(m5c.kidney) = colnames(RCircos.Histogram.Data)
RCircos.Chromosome.Region.Plot(m5c.kidney, track.num, side);
track.num =  track.num + 1
m5c.lung = read.table("./lung.000000F.gff.m5c",header=F,sep="\t",stringsAsFactors=F)[,c(1,4,5,9)]
m5c.lung[,1] = sub("F","",paste('chr',m5c.lung[,1],sep=""))
colnames(m5c.lung) = colnames(RCircos.Histogram.Data)
RCircos.Chromosome.Region.Plot(m5c.lung, track.num, side);
print("Done plot hypomethylated region")

exp.files = list.files(".","*.clean")
i = 1
heatmap.exp = chromRangeTable
for(i in 1:length(exp.files)){
    sample.name = str_split(exp.files[i],"\\.")[[1]][1]
    assign(sample.name, read.table(exp.files[i],header=F,sep="\t",stringsAsFactors=F))
    tmp.table = data.table(chr=rep("chr000000", nrow(get(sample.name))), start=as.numeric(get(sample.name)[,2]), end=as.numeric(get(sample.name)[,2]), key = c("chr","start","end"))
    tmp.res = foverlaps(tmp.table, chromRangeTable, nomatch=0)
    count = data.frame(table(tmp.res$start))
    colnames(count) = c("start","count")
    tmp.res = merge(chromRangeFrame, count, by="start",all.x=T)
    tmp.res[is.na(tmp.res)] = 0
    tmp.res[,5] = log(tmp.res[,5] + 1)
    heatmap.exp = merge(heatmap.exp, tmp.res[,c(1,5)], by = "start",all.x=T)
    colnames(heatmap.exp)[ncol(heatmap.exp)] = sample.name
}
heatmap.exp = heatmap.exp[-nrow(heatmap.exp),]
heatmap.exp[,c(1,2)] = heatmap.exp[,c(2,1)]


source("./RCircos.Heatmap.Plot.R")
total.track <- ncol(heatmap.exp) - 4;
for(a.track in 1:total.track)
#--------------------------------------------------
# for(a.track in 3:4)
#--------------------------------------------------
{
    data.col <- a.track + 4;
    track.num =  track.num + 1
    #--------------------------------------------------
    # RCircos.Heatmap.Plot(heatmap.exp, data.col, track.num, "out");
    #--------------------------------------------------
    exp =as.data.frame( heatmap.exp[,c(1:4,data.col),with=FALSE])
    my.RCircos.Heatmap.Plot(exp, 5, track.num, "out");
}
print("Done plot gene expression heatmap")


gc.content = read.table("./000000F.fa.gc_content",header=F,sep=" ",stringsAsFactors=F)
colnames(gc.content) = c("chr","start","end","gc")
gc.content[,2] = gc.content[,2] - 1
gc.content[,3] = gc.content[,3] - 1
gc.content[,4] = gc.content[,4] * rcircos.params$base.per.unit

m5c.kidney.site = read.table("./kidney.000000F_class.wig",header=T,sep=" ",stringsAsFactors=F)
m5c.kidney.site = m5c.kidney.site[m5c.kidney.site[,2] == 1,]
m5c.kidney.site.table = data.table(chr=rep("chr000000", nrow(m5c.kidney.site)), start=as.numeric(unlist(m5c.kidney.site[,1])), end=as.numeric(unlist(m5c.kidney.site[,1])), count=m5c.kidney.site[,2], key = c("chr","start","end"))
m5c.kidney.site.table = m5c.kidney.site.table[order(as.numeric(m5c.kidney.site.table$start)),]
m5c.kidney.site.res = foverlaps(m5c.kidney.site.table, chromRangeTable, nomatch=0)
count = data.frame(table(m5c.kidney.site.res$start))
colnames(count) = c("start","count")
m5c.kidney.site.res = merge(chromRangeFrame, count, by="start",all=T)
m5c.kidney.site.res[is.na(m5c.kidney.site.res)] = 0
m5c.kidney.site.res[,5] = m5c.kidney.site.res[,5] / gc.content[,4] * rcircos.params$base.per.unit
m5c.lung.site = read.table("./lung.000000F_class.wig",header=T,sep=" ",stringsAsFactors=F)
m5c.lung.site = m5c.lung.site[m5c.lung.site[,2] == 1,]
m5c.lung.site.table = data.table(chr=rep("chr000000", nrow(m5c.lung.site)), start=as.numeric(unlist(m5c.lung.site[,1])), end=as.numeric(unlist(m5c.lung.site[,1])), count=m5c.lung.site[,2], key = c("chr","start","end"))
m5c.lung.site.table = m5c.lung.site.table[order(as.numeric(m5c.lung.site.table$start)),]
m5c.lung.site.res = foverlaps(m5c.lung.site.table, chromRangeTable, nomatch=0)
count = data.frame(table(m5c.lung.site.res$start))
colnames(count) = c("start","count")
m5c.lung.site.res = merge(chromRangeFrame, count, by="start",all=T)
m5c.lung.site.res[is.na(m5c.lung.site.res)] = 0
m5c.lung.site.res[,5] = m5c.lung.site.res[,5] / gc.content[,4] * rcircos.params$base.per.unit

all.site.res = merge(m5c.kidney.site.res, m5c.lung.site.res, by = "start",all=T)
heatmap.m5c = data.table(Chromosome=rep("chr000000", nrow(all.site.res)),
                         chromStart=as.numeric(all.site.res$start),
                         chromEnd=as.numeric(all.site.res$start)+rcircos.params$base.per.unit-1,
                         Name=paste("region",1:nrow(all.site.res),sep=""),
                         kidney.site=as.numeric(all.site.res$count.x),
                         lung.site=as.numeric(all.site.res$count.y))
heatmap.m5c[is.na(heatmap.m5c)] = 0
heatmap.m5c = heatmap.m5c[-nrow(heatmap.m5c),]
total.track <- 2;
for(a.track in 1:total.track)
{
    data.col <- a.track + 4;
    track.num =  track.num + 1
    RCircos.Heatmap.Plot(heatmap.m5c, data.col, track.num, "out");
}
print("Done plot m5c methylation heatmap")

region.name = heatmap.exp[seq(1,nrow(heatmap.exp),4),]
name.col <- 4;
track.num <- track.num + 1;
RCircos.Gene.Connector.Plot(region.name, track.num, side);
track.num <- track.num + 1;
RCircos.Gene.Name.Plot(region.name, name.col,track.num, side);

gc.islands = read.table("./000000F.fa.island",header=T,sep="\t",stringsAsFactors=F)[,c(1,3,4,5)]
gc.islands[,1] = sub("F","",paste('chr',gc.islands[,1],sep=""))
track.num =  track.num + 2
RCircos.Chromosome.Region.Plot(gc.islands, track.num, side);
print("Done plot GC islands region")
track.num =  track.num + 1
RCircos.Chromosome.Region.Plot(all.gene, track.num, side);
print("Done plot all gene region")

reps = read.table("./rep.data.LINE",header=F,sep=" ",stringsAsFactors=F)
reps[,1] = sub("F","",paste('chr',reps[,1],sep=""))
colnames(reps) = colnames(RCircos.Histogram.Data)
reps.table = data.table(chr=reps[,1], start=reps[,2],end=reps[,3],Name=reps[,4],key = c("chr","start","end"))
reps.res = foverlaps(reps.table, chromRangeTable, nomatch=0)
count = data.frame(table(reps.res$start))
colnames(count) = c("start","count")
reps.res = merge(chromRangeFrame, count, by="start")
reps.res[,c(1,2)] = reps.res[,c(2,1)]
data.col <-  5;
track.num =  track.num + 1
RCircos.Heatmap.Plot(reps.res, data.col, track.num, "out");
print("Done plot repeat density heatmap")

dev.off();
