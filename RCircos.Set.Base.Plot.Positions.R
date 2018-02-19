my.RCircos.Set.Base.Plot.Positions<-function()
{
	RCircos.Cyto <- RCircos.Get.Plot.Ideogram();
	RCircos.Par <- RCircos.Get.Plot.Parameters();


	#	Add one padding length to the last chromosome. Others
	#	are already included in RCircos.Cyto$Location
	# 	_________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	total.points <- RCircos.Cyto$Location[nrow(RCircos.Cyto)] +  RCircos.Par$chrom.paddings;


	#	x and y coordinates for a circlar line with radius of 1,
	#	circumferance of 2*PI, and interval of 2PI/total.points
	# 	_________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	interval <- 2*pi/total.points;
    print(interval)
	base.val <- seq(0, 2*pi, interval);

	cor.x <- sin(base.val);
	cor.y <- cos(base.val);


	#	Degrees for text rotating at each posint
	# 	_________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	degree <- rep(0, length(base.val));
	mid <- round((length(base.val)-1)/2, digits=0) + 1;

	for(pt in 1:mid)
	{ degree[pt] <- 90 - (base.val[pt]*180/pi);}
	
	for(pt in (mid+1):length(base.val))
	{ degree[pt] <- 270 - (base.val[pt]*180/pi); }


	#	Put the plot postions data in RCircos environment
	# 	_________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	plot.postions <- data.frame(cor.x, cor.y, degree);

	RCircosEnvironment <- NULL;
	RCircosEnvironment <- get("RCircos.Env", envir=globalenv());
	RCircosEnvironment[["RCircos.Base.Position"]] <- plot.postions;
}
