RCircos.Chromosome.Ideogram.Plot = function () 
{
    params <- RCircos.Get.Plot.Parameters(); 
    #--------------------------------------------------
    # params$base.per.unit <- 20; 
    # params$chrom.paddings = 10
    # params$radius.len <- 1; 
    # RCircos.Reset.Plot.Parameters(params); 
    #--------------------------------------------------
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
    RCircos.Pos <- RCircos.Get.Plot.Positions()
    RCircos.Par <- RCircos.Get.Plot.Parameters()
    outer.location <- RCircos.Par$chr.ideog.pos + RCircos.Par$chrom.width
    inner.location <- RCircos.Par$chr.ideog.pos
    chroms <- unique(RCircos.Cyto$Chromosome)
    for (a.chr in 1:length(chroms)) {
        the.chr <- RCircos.Cyto[RCircos.Cyto$Chromosome == chroms[a.chr],]
        start <- the.chr$Location[1] - the.chr$Unit[1] + 1
        end <- the.chr$Location[nrow(the.chr)]
        mid <- round((end - start + 1)/2, digits = 0) + start
        chr.color <- the.chr$ChrColor[nrow(the.chr)]
        pos.x <- c(RCircos.Pos[start:end, 1] * outer.location, RCircos.Pos[end:start, 1] * inner.location)
        pos.y <- c(RCircos.Pos[start:end, 2] * outer.location, RCircos.Pos[end:start, 2] * inner.location)
        polygon(pos.x, pos.y)
        chr.name <- sub(pattern = "chr", replacement = "", chroms[a.chr])
        #--------------------------------------------------
        # text(RCircos.Pos[mid, 1] * RCircos.Par$chr.name.pos, RCircos.Pos[mid, 2] * RCircos.Par$chr.name.pos, label = chr.name, srt = RCircos.Pos$degree[mid])
        #-------------------------------------------------- 
        chr.color="black"
        lines(RCircos.Pos[start:end, ] * RCircos.Par$highlight.pos, col = chr.color, lwd = RCircos.Par$highlight.width)

        pos.tick.in.x <- RCircos.Pos[start:end, 1] * RCircos.Par$highlight.pos
        pos.tick.in.y <- RCircos.Pos[start:end, 2] * RCircos.Par$highlight.pos
        pos.tick.out.x.big <- RCircos.Pos[start:end, 1] * ((outer.location-inner.location)/2 +RCircos.Par$highlight.pos)
        pos.tick.out.y.big <- RCircos.Pos[start:end, 2] * ((outer.location-inner.location)/2 +RCircos.Par$highlight.pos)
        pos.tick.out.x.small <- RCircos.Pos[start:end, 1] * ((outer.location-inner.location)/4 +RCircos.Par$highlight.pos)
        pos.tick.out.y.small <- RCircos.Pos[start:end, 2] * ((outer.location-inner.location)/4 +RCircos.Par$highlight.pos)
        pos.tick.text.x <- RCircos.Pos[start:end, 1] * ((outer.location-inner.location)/1.8 +RCircos.Par$highlight.pos)
        pos.tick.text.y <- RCircos.Pos[start:end, 2] * ((outer.location-inner.location)/1.8 +RCircos.Par$highlight.pos)

        for (pos in start : end){
            if((pos-start)%%10000 == 0 & pos != start){
                segments(pos.tick.in.x[pos], pos.tick.in.y[pos], pos.tick.out.x.big[pos], pos.tick.out.y.big[pos], col="black")
            }else if ( ((pos-start) * params$base.per.unit) %% 1000 == 0) {
                segments(pos.tick.in.x[pos], pos.tick.in.y[pos], pos.tick.out.x.small[pos], pos.tick.out.y.small[pos], col="black")
                text(pos.tick.text.x[pos],pos.tick.text.y[pos], paste( ((pos-start) * params$base.per.unit) / 1000, "k",sep=""),col="black",cex=0.5,font=4)
            }
        }
    }
    for (a.band in 1:nrow(RCircos.Cyto)) {
        a.color <- RCircos.Cyto$BandColor[a.band]
        if (a.color == "white") {
            next
        }
        start <- RCircos.Cyto$Location[a.band] - RCircos.Cyto$Unit[a.band] + 1
        end <- RCircos.Cyto$Location[a.band]
        pos.x <- c(RCircos.Pos[start:end, 1] * outer.location, RCircos.Pos[end:start, 1] * inner.location)
        pos.y <- c(RCircos.Pos[start:end, 2] * outer.location, RCircos.Pos[end:start, 2] * inner.location)
        polygon(pos.x, pos.y, col = a.color, border = NA)
    }
    arctext(x = "Gene", center = c(0, 0), radius = inner.location - RCircos.Par$track.height/2, middle = pi/4, cex=0.6)
}
