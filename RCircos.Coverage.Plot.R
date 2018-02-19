RCircos.Coverage.Plot= function (cov.info, data.col, track.begin, track.end, side, norm.index, y.lim, legend.text) 
{
    if (track.num < 1) {
        stop("Track number cannot be smaller than 1.\n")
    }
    side <- tolower(side)
    if (side != "in" && side != "out") {
        stop("side must be either in or out.\n")
    }
    if (data.col < 4) {
        stop("heatmap data column must be 4 or bigger.\n")
    }
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
    RCircos.Pos <- RCircos.Get.Plot.Positions()
    RCircos.Par <- RCircos.Get.Plot.Parameters()
    locations <- RCircos.Track.Positions(side, track.begin)
    inner.location <- locations[2]
    locations <- RCircos.Track.Positions(side, track.end)
    outer.location <- locations[1]
    chroms <- unique(RCircos.Cyto$Chromosome)
    a.chr=1
    for (a.chr in 1:length(chroms)) {
        the.chr <- RCircos.Cyto[RCircos.Cyto$Chromosome == chroms[a.chr],]
        start <- the.chr$Location[1] - the.chr$Unit[1] + 1
        end <- the.chr$Location[nrow(the.chr)]


        rownames(cov.info) = cov.info[,1]
        cov.top = max(cov.info$Total.Depth) / norm.index  * 1000000
        cov.length = abs(outer.location - inner.location)
        cov.pos = (start:end -start)*RCircos.Par$base.per.unit
        cov.pos[1] = cov.info[1,1] 
        cov.height = cov.info[as.character(cov.pos), "Total.Depth"] / norm.index  * 1000000

        #--------------------------------------------------
        # label.index = as.integer(cov.top/c(10^(1:5)))
        # label.top.index = which(label.index>0 & label.index<10)
        # label.top = (label.index[label.top.index] + 1) * 10^label.top.index
        #-------------------------------------------------- 
        label.top = y.lim
        label.pos = which(cov.height==min(cov.height))

        degree.sign = sign(RCircos.Pos[start:end,"degree"])
        negative.begin = min(which(degree.sign == -1))
        postive.begin = min(which(degree.sign[negative.begin:length(degree.sign)] == 1)) + negative.begin - 1
        negative.end = postive.begin - 1
        postive.end = min(which(degree.sign[postive.begin:length(degree.sign)] == -1)) + postive.begin - 1
        degree.sign[negative.begin:negative.end] = 1
        degree.sign[postive.begin:postive.end] = -1
        pos.x <- c(RCircos.Pos[start:end, 1] * inner.location + cov.length * cov.height / label.top *
                   degree.sign * cos(RCircos.Pos[start:end,"degree"]/180*pi), 
                   RCircos.Pos[end:start, 1] * inner.location)
        pos.y <- c(RCircos.Pos[start:end, 2] * inner.location + cov.length * cov.height / label.top *
                   degree.sign * sin(RCircos.Pos[start:end,"degree"]/180*pi), 
                   RCircos.Pos[end:start, 2] * inner.location)

        polygon(pos.x, pos.y,col="brown")
        for (i in track.begin:track.end){
            locations <- RCircos.Track.Positions(side, i)
            out.pos <- locations[1]
            in.pos <- locations[2]
            lines(RCircos.Pos[start:end, ] * in.pos, col = RCircos.Par$grid.line.color, lty=3)
            arctext(x = as.character(label.top/(track.end-track.begin)*(i-track.begin)), center = c(0, 0), radius = in.pos + RCircos.Par$track.height/2, 
                    middle = 4.9*pi/10,
                    cex=0.6)
        }

    }
    arctext(x = legend.text , center = c(0, 0), radius = inner.location - RCircos.Par$track.height/2, middle = pi/4, cex=0.6)
}
