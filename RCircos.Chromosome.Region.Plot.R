RCircos.Get.Plot.Region<-function(genomic.data, plot.type)
{
	#	Check chromosome names, chromStart, and chromEnd positions
	#	_______________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	genomic.data <- RCircos.Validate.Genomic.Data(genomic.data, plot.type);


	#	Calculate the point index for each chromosome location
	#	_______________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	data.points <- rep(0, nrow(genomic.data));
	start.points <- rep(0, nrow(genomic.data));
	end.points <- rep(0, nrow(genomic.data));
	for(a.row in 1:nrow(genomic.data))
	{
		chromosome <- as.character(genomic.data[a.row, 1]);
		location <- round((genomic.data[a.row, 2] + genomic.data[a.row, 3])/2, digits=0);
		start.location <- genomic.data[a.row,2]
		end.location <- genomic.data[a.row,3]
		data.points[a.row] <- RCircos.Data.Point(chromosome, location);
		start.points[a.row] <- RCircos.Data.Point(chromosome, start.location);
		end.points[a.row] <- RCircos.Data.Point(chromosome, end.location);
	}
	genomic.data["Location"] <- data.points;
	genomic.data["Start"] <- start.points;
	genomic.data["End"] <- end.points;


	#	Sort the data by chromosome then start position
	#	____________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	genomic.data <- genomic.data[order(genomic.data$Location),];


	#	The data needs to be held for the RCircos session
	#	____________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	return (genomic.data);
}

RCircos.Chromosome.Region.Plot<-function(region, track.num, side)
{
	RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
	RCircos.Pos <- RCircos.Get.Plot.Positions()
	RCircos.Par <- RCircos.Get.Plot.Parameters()


	#	Plot chromosome outlines, chromosome names, and 
	#	chromosome highlights
	#	_________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    region.data <- RCircos.Get.Plot.Region(region, "plot");

	start <-as.numeric(region.data[, "Start"]);
	end <-as.numeric(region.data[, "End"]);

	right.side <- nrow(RCircos.Pos)/2;
	outer.location <- RCircos.Par$chr.ideog.pos + RCircos.Par$chrom.width;
	inner.location <-RCircos.Par$chr.ideog.pos

	data.chroms <- as.character(region[,1]);
	chromosomes <- unique(data.chroms);
	cyto.chroms <- as.character(RCircos.Cyto$Chromosome);

    #--------------------------------------------------
    # for(a.chr in 1:length(chromosomes))
    # {
    #     cyto.rows <- which(cyto.chroms==chromosomes[a.chr]);
    #     locations <- as.numeric(RCircos.Cyto$Location[cyto.rows]);
    #     chr.start <- min(locations) - RCircos.Cyto$Unit[cyto.rows[1]];
    #     chr.end   <- max(locations);
    #--------------------------------------------------

    #--------------------------------------------------
    #     data.rows <- which(data.chroms==chromosomes[a.chr]);
    #     start[data.rows[start[data.rows]<chr.start]] <- chr.start;
    #     end[data.rows[end[data.rows]>chr.end]] <- chr.end;
    # }
    #--------------------------------------------------

	locations <- RCircos.Track.Positions(side, track.num);
	out.pos <- locations[1];
	in.pos <- locations[2];

    source("./RCircos.Track.Outline.Simple.R")
    RCircos.Track.Outline.Simple(out.pos, in.pos);
    a.band =1 
    for (a.band in 1:nrow(region.data)) {
        a.color <- "gray"
        if (a.color == "white") {
            next
        }
        a.start <- start[a.band]
        a.end <- end[a.band]

        #--------------------------------------------------
        # pos.x <- c(RCircos.Pos[start:end, 1] * outer.location, RCircos.Pos[end:start, 1] * inner.location)
        # pos.y <- c(RCircos.Pos[start:end, 2] * outer.location, RCircos.Pos[end:start, 2] * inner.location)
        # polygon(pos.x, pos.y, col = a.color, border = NA)
        # the.chr  <- RCircos.Cyto[RCircos.Cyto$Chromosome==chroms[a.chr],];
        # start <- the.chr$Location[1]- the.chr$Unit[1] + 1;
        # end   <- the.chr$Location[nrow(the.chr)];
        #--------------------------------------------------

		polygon.x<- c(RCircos.Pos[a.start:a.end,1]*out.pos, RCircos.Pos[a.end:a.start,1]*in.pos);
		polygon.y<- c(RCircos.Pos[a.start:a.end,2]*out.pos, RCircos.Pos[a.end:a.start,2]*in.pos);
		polygon(polygon.x, polygon.y, col=a.color);
    }

    #--------------------------------------------------
    # for(a.chr in 1:length(chroms))
    # {
    #     the.chr  <- RCircos.Cyto[RCircos.Cyto$Chromosome==chroms[a.chr],];
    #     start <- the.chr$Location[1]- the.chr$Unit[1] + 1;
    #     end   <- the.chr$Location[nrow(the.chr)];
    #     mid <- round((end-start+1)/2, digits=0)+start;
    #     chr.color <- the.chr$ChrColor[nrow(the.chr)];
    #     #	Draw chromosome outlines
    #     #	_________________________________________________________
    #     #	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    #     pos.x<- c(RCircos.Pos[start:end,1]*outer.location, RCircos.Pos[end:start,1]*inner.location);
    #     pos.y<- c(RCircos.Pos[start:end,2]*outer.location, RCircos.Pos[end:start,2]*inner.location);
    #     polygon(pos.x, pos.y);
    #     #	Add chromosome names
    #     #	_________________________________________________________
    #     #	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    #     chr.name <- sub(pattern="chr", replacement="", chroms[a.chr]);
    #--------------------------------------------------
        #--------------------------------------------------
        # text(RCircos.Pos[mid,1]*RCircos.Par$chr.name.pos, RCircos.Pos[mid,2]*RCircos.Par$chr.name.pos,
        #      label=chr.name, srt=RCircos.Pos$degree[mid]);
        #--------------------------------------------------
    #--------------------------------------------------
    #     #	Add chromosome highlights
    #     #	_________________________________________________________
    #     #	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    #     lines(RCircos.Pos[start:end,]*RCircos.Par$highlight.pos, col=chr.color, lwd=RCircos.Par$highlight.width);
    # }
    #--------------------------------------------------
}

