
my.RCircos.Heatmap.Plot<-function(heatmap.data, data.col, track.num, side)
{
	RCircos.Cyto <- RCircos.Get.Plot.Ideogram();
	RCircos.Pos <- RCircos.Get.Plot.Positions();
	RCircos.Par <- RCircos.Get.Plot.Parameters();


	#	Convert raw data to plot data. The raw data will be validated
	#	first during the convertion
	#	____________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	heatmap.data <- RCircos.Get.Plot.Data(heatmap.data, "plot");
	heatmap.data <- RCircos.Get.Plot.Data(exp, "plot");


	#	Colors for different data values
	#	___________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	RedRamp <- rgb( seq(1, 1, length=256),  seq(0, 1, length=256),  seq(0, 1, length=256)); 
	BlueRamp <- rgb(seq(0, 1, length=256),  seq(0, 1, length=256),  seq(1, 1, length=256));		
	ColorRamp   <- cbind(BlueRamp, rev(RedRamp));


	#	Color level. Heatmap data has to have four leading columns
	#	for genomic position and gene(lable) names. Also color 
	#	level must be calculated with all data columns.
	#	___________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	columns <- 5:(ncol(heatmap.data)-1);
	min.value <- min(as.matrix(heatmap.data[, columns]));
	max.value <- max(as.matrix(heatmap.data[, columns]));
	ColorLevel  <- seq(min.value, max.value, length=length(ColorRamp));


	#	Each heatmap cell is centered on data point location. Make
	#	sure each one will be in the range of it chromosome
	#	____________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	heatmap.locations <- as.numeric(heatmap.data[, ncol(heatmap.data)]) + 1;
	start <- heatmap.locations - RCircos.Par$heatmap.width/2;
	end   <- heatmap.locations + RCircos.Par$heatmap.width/2;

	data.chroms <- as.character(heatmap.data[,1]);
	chromosomes <- unique(data.chroms);
	cyto.chroms <- as.character(RCircos.Cyto$Chromosome);

	for(a.chr in 1:length(chromosomes))
	{
		cyto.rows <- which(cyto.chroms==chromosomes[a.chr]);
		locations <- as.numeric(RCircos.Cyto$Location[cyto.rows]);
		chr.start <- min(locations) - RCircos.Cyto$Unit[cyto.rows[1]];
		chr.end   <- max(locations);

		data.rows <- which(data.chroms==chromosomes[a.chr]);
		start[data.rows[start[data.rows]<chr.start]] <- chr.start;
		end[data.rows[end[data.rows]>chr.end]] <- chr.end;
	}


	#	Plot position for current track. 
	#	___________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	locations <- RCircos.Track.Positions(side, track.num);
	out.pos <- locations[1];
	in.pos <- locations[2];


	#	outline of chromosomes. No lines inside.
	#	___________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	chroms <- unique(RCircos.Cyto$Chromosome);
	for(a.chr in 1:length(chroms))
	{
		the.chr  <- RCircos.Cyto[RCircos.Cyto$Chromosome==chroms[a.chr],];
		the.start <- the.chr$Location[1]- the.chr$Unit[1] + 1;
		the.end   <- the.chr$Location[nrow(the.chr)];

		polygon.x<- c(RCircos.Pos[the.start:the.end,1]*out.pos, 
				RCircos.Pos[the.end:the.start,1]*in.pos);
		polygon.y<- c(RCircos.Pos[the.start:the.end,2]*out.pos, 
				RCircos.Pos[the.end:the.start,2]*in.pos);
		polygon(polygon.x, polygon.y, col="white");
	}


	#	Plot heatmap for each gene.
	#	_______________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	heatmap.value <- as.numeric(heatmap.data[, data.col]);
    a.point=1
	for(a.point in 1:length(heatmap.value))
	{
		the.level <- which(ColorLevel>=heatmap.value[a.point]);
		cell.color <- ColorRamp[min(the.level)];
		
		the.start <- start[a.point];
		the.end <- end[a.point];

        polygon.x<- c(RCircos.Pos[the.start:the.end,1]*out.pos, RCircos.Pos[the.end:the.start,1]*in.pos);
        polygon.y<- c(RCircos.Pos[the.start:the.end,2]*out.pos, RCircos.Pos[the.end:the.start,2]*in.pos);
		polygon(polygon.x, polygon.y, col=cell.color, border=NA);
	}

}
