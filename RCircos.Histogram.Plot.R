my.RCircos.Data.Point<-function(chromosome, start)
{
	the.point <- 0;
	RCircos.Cyto <- RCircos.Get.Plot.Ideogram();


	#	Which band the start position is in
	#	_________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	chrom.rows <- grep(paste("^", chromosome, "$", sep=""), RCircos.Cyto$Chromosome);
	the.row <- which(RCircos.Cyto$ChromStart[chrom.rows] <= start &  RCircos.Cyto$ChromEnd[chrom.rows] >= start)[1];


	#	total length, units, and location of the band
	#	_________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	band.length <- RCircos.Cyto$Length[chrom.rows[the.row]];
	band.units <- RCircos.Cyto$Unit[chrom.rows[the.row]];
	band.location <- RCircos.Cyto$Location[chrom.rows[the.row]];


	#	How far from the chromosome start the point is (by units)
	#	_________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	the.bases <- start - RCircos.Cyto$ChromStart[chrom.rows[the.row]] ;
	the.units  <- the.bases/band.length*band.units;


	#	The exact point index of points for the circular line
	#	_________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	the.point <- band.location - band.units + the.units;


	#	Return the index. Arguments are valudated outside of this 
	#	function so that there is no need to catch exception.
	#	_________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	return (round(the.point, digits=0));
}

RCircos.Validate.Genomic.Data<-function(genomic.data, plot.type=c("plot", "link"))
{
	RCircos.Cyto <- RCircos.Get.Plot.Ideogram();
	

	#	Plot data has only one chromosome column and link data have two
	#	_______________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	plot.type <- tolower(plot.type);
	if(plot.type=="plot") { 
		chrom.col <- 1; 
	} else if(plot.type=="link") { 
		chrom.col <- c(1,4); 
	} else { stop("Plot type must be \"plot\" or \"line\""); }
	
	for (a.col in 1:length(chrom.col))
	{
		the.col <- chrom.col[a.col];


		#	Make sure chromosome names in genomic data have prefix
		#	______________________________________________________
		#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

		genomic.data[,the.col] <- as.character(genomic.data[,the.col]);
		for(a.row in 1:nrow(genomic.data)) {
			if(length(grep("chr", genomic.data[a.row,the.col]))==0) 
			{ genomic.data[a.row,the.col] <- paste("chr", 
				genomic.data[a.row,the.col], sep=""); }
		}


		#	Make sure chromosomes in input data are all included 
		#	in chromosome ideogram data
		#	______________________________________________________
		#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

		cyto.chroms <- unique(as.character(RCircos.Cyto$Chromosome));
		data.chroms <- unique(as.character(genomic.data[,the.col]));
		if(sum(data.chroms %in% cyto.chroms) < length(data.chroms)) 
		{ 
			cat(paste("Some chromosomes are in genomic data only",
			 "and have been removed.\n\n"));

            all.chroms <- as.character(genomic.data[,the.col]);
            genomic.data <- genomic.data[all.chroms %in% cyto.chroms,];
		}
		data.chroms <- unique(as.character(genomic.data[,the.col]));


		#	Make sure chromosome start and end postions in genomic
		#	data are not negative.
		#	______________________________________________________
		#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

		if(min(genomic.data[,the.col+1])<0) 
		{ stop("Error! chromStart position less than 0."); }
		if(min(genomic.data[,the.col+2])<0) 
		{ stop("Error! chromEnd position less than 0.");  }	


		#	Make sure chromosome start and end locations in genomic
		#	data are not out of chromosome length
		#	_______________________________________________________
		#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

		for(a.chr in 1:length(data.chroms))
		{
			the.chr      <- data.chroms[a.chr]; 
			in.data      <- genomic.data[genomic.data[,the.col]==the.chr,];
			cyto.data <- RCircos.Cyto[grep(the.chr, RCircos.Cyto$Chromosome),]

			if(max(in.data[,the.col+1])>max(cyto.data[,3]) | 
				max(in.data[,the.col+2])>max(cyto.data[,3]))
			{  
				cat(paste(the.chr, max(in.data[,2]), max(in.data[,3]), "\n"));
 				stop("Error! Location is outside of chromosome length.");  
			}
		}

		
		#	Make sure in genomic data all chromosome start positions 
		#	are smaller than their paired chromosome end positions
		#	_______________________________________________________________
		#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

		for(a.row in 1:nrow(genomic.data))
		{
			if(genomic.data[a.row, the.col+1]>genomic.data[a.row, the.col+2]) 
			{ 
				cat("chromStart greater than chromEnd.\n"); 
				stop(paste("Row:", a.row, genomic.data[a.row, 2],  
					genomic.data[a.row, 3]));
			}
		}
	}


	#	The validated data needs to be held for the RCircos session
	#	___________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	return (genomic.data);
}

my.RCircos.Get.Plot.Data<-function(genomic.data, plot.type)
{
	#	Check chromosome names, chromStart, and chromEnd positions
	#	_______________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	genomic.data <- RCircos.Validate.Genomic.Data(genomic.data, plot.type);


	#	Calculate the point index for each chromosome location
	#	_______________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	data.points <- rep(0, nrow(genomic.data));
	for(a.row in 1:nrow(genomic.data))
	{
		chromosome <- as.character(genomic.data[a.row, 1]);
		location <- round((genomic.data[a.row, 2] + genomic.data[a.row, 3])/2, digits=0);
		data.points[a.row] <- my.RCircos.Data.Point(chromosome, location);
	}
	genomic.data["Location"] <- data.points;


	#	Sort the data by chromosome then start position
	#	____________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	genomic.data <- genomic.data[order(genomic.data$Location),];


	#	The data needs to be held for the RCircos session
	#	____________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	return (genomic.data);
}

my.RCircos.Histogram.Plot<-function(hist.data, data.col, track.num, side)
{
	RCircos.Pos <- RCircos.Get.Plot.Positions();
	RCircos.Par <- RCircos.Get.Plot.Parameters();
	RCircos.Cyto <- RCircos.Get.Plot.Ideogram();


	#	Convert raw data to plot data. The raw data will be validated
	#	first during the convertion
	#	____________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	hist.data <- RCircos.Get.Plot.Data(hist.data, "plot");
    #--------------------------------------------------
    # hist.data <- my.RCircos.Get.Plot.Data(ipd.lung.m6A, "plot");
    #--------------------------------------------------
    #--------------------------------------------------
    # hist.data = hist.data[ !is.na(hist.data[,5]), ]
    #--------------------------------------------------
    #--------------------------------------------------
    # print(head(hist.data))
    #--------------------------------------------------


	#	Each heatmap cell is centered on data point location. Make
	#	sure each one will be in the range of it chromosome
	#	____________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	hist.locations <- as.numeric(hist.data[, ncol(hist.data)]);
	start <- hist.locations - RCircos.Par$hist.width;
	end   <- hist.locations + RCircos.Par$hist.width;

	data.chroms <- as.character(hist.data[,1]);
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


	#	Draw histogram
	#	___________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	num.subtrack <- RCircos.Par$sub.tracks;
	RCircos.Track.Outline(out.pos, in.pos, RCircos.Par$sub.tracks);

	for(a.point in 1:nrow(hist.data))
	{
		hist.height <- hist.data[a.point, data.col];

		the.start <- start[a.point];
		the.end <- end[a.point];

		#	Plot rectangle with specific height for each 
		#	data point
		#	_______________________________________________
		#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

		height <- in.pos + RCircos.Par$track.height*hist.height;

		polygon.x<- c(RCircos.Pos[the.start:the.end,1]*height,
				RCircos.Pos[the.end:the.start,1]*in.pos);
		polygon.y<- c(RCircos.Pos[the.start:the.end,2]*height,
				RCircos.Pos[the.end:the.start,2]*in.pos);
        polygon(polygon.x, polygon.y, col="red", border=NA);
	}
}
