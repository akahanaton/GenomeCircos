my.RCircos.Set.Cytoband.data<-function(cyto.band.info)
{
	
	#       Reset colors for chromosome bands. Use yellow color for unknow 
	#	______________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	stain2color <- as.character(cyto.band.info$Stain);
	band.color <- rep(colors()[652], length(stain2color));
	
	stains <- c("gneg", "acen", "stalk", "gvar", "gpos", "gpos100", 
		"gpos75", "gpos66", "gpos50", "gpos33", "gpos25");
	color.index <- c(1, 552, 615, 418, 24, 24, 193, 203, 213, 223, 233);

	for(a.stain in 1:length(stains))
	{
		bands <- which(stain2color==stains[a.stain]);
		if(length(bands)>0) 
		{ band.color[bands] <- colors()[color.index[a.stain]]; }
	}
	cyto.band.info["BandColor"] <- band.color;


	#	Assign colors to chromosome highlight. There are total 50
	#	colors and the last 26 colors are reserved for future.
	#	___________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	
	chrom.color <- c(552, 574, 645, 498, 450, 81, 26, 584, 524, 472,
			32, 57, 615, 635, 547, 254, 100, 72, 630, 589,
			8, 95, 568, 52);

	chrom2color <- as.character(cyto.band.info$Chromosome);
	chromosomes <- unique(chrom2color);


	#	In case of multiple ideogram plot, recycle the colors
	#	__________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	num.chrom <- length(chromosomes);
	num.color <- length(chrom.color);
	if(num.chrom>num.color)
	{
		recycle.time <- floor(num.chrom/num.color);
		if(recycle.time>1) 
		{ chrom.color <- rep(chrom.color, recycle.time); }

		remains <- num.chrom%%num.color
		if(remains > 0)
		{  chrom.color <- c(chrom.color, chrom.color[1:remains]); }
	}

	for(a.chr in 1:length(chromosomes))
	{
		rows <- which(chrom2color==chromosomes[a.chr]);
		if(length(rows)>0)
		{ chrom2color[rows] <- colors()[chrom.color[a.chr]]; }
	}
	cyto.band.info["ChrColor"] <- chrom2color;


	#	Total base pairs and relative length of each band
	#	__________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	plot.par <- RCircos.Get.Plot.Parameters();

	cyto.band.info$ChromStart <- as.numeric(cyto.band.info$ChromStart);
	cyto.band.info$ChromEnd <- as.numeric(cyto.band.info$ChromEnd);

	band.len <- cyto.band.info$ChromEnd - cyto.band.info$ChromStart;
	cyto.band.info["Length"] <- band.len;
	cyto.band.info["Unit"]<- round(band.len/plot.par$base.per.unit, digits=0);
	

	#	Relative locations of each band in clockwise
	#	__________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	Relative.Loc <- cyto.band.info$Unit;
	for(i in 2:length(Relative.Loc))
    { Relative.Loc[i] <- Relative.Loc[i] + Relative.Loc[i-1]; }
	cyto.band.info["Location"] <- Relative.Loc;

	if( plot.par$chrom.paddings>0) 
	{  
		chroms <- unique(cyto.band.info$Chromosome);
		chroms <- chroms[(chroms==chroms[1])==F];
		num.pad <-  plot.par$chrom.paddings;

		for(a.chr in 1:length(chroms))
		{
			index <- grep(paste(chroms[a.chr], "$", sep=""), cyto.band.info$Chromosome);
			cyto.band.info$Location[index] <- num.pad + cyto.band.info$Location[index];
			num.pad <- num.pad +  plot.par$chrom.paddings;
		}
	}
    print(cyto.band.info["Location"])


	#	Put the cyto.band.info data in RCircos environment
	# 	______________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	RCircosEnvironment <- NULL;
	RCircosEnvironment <- get("RCircos.Env", envir=globalenv());
	RCircosEnvironment[["RCircos.Cytoband"]] <- cyto.band.info;
}
