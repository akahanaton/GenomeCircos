\name{RCircos.Tile.Plot}
\alias{RCircos.Tile.Plot}
\title{
Circos-like Tile Plot
}
\description{
Draw one track of Tiles. RCircos core components and graphics device must be initialized before drawing.
}
\usage{
RCircos.Tile.Plot(tile.data, track.num, side)
}
\arguments{
  \item{tile.data}{
A data frame with three columns for chromosome name, start and end positions of each data point. Columns for data values are optional.
}
  \item{track.num}{
Integer, representing the ordinal number of the plot track where the data will be plotted.
}
  \item{side}{
Character vector, either "in" or "out", representing the position related to chromosome ideogram.
}
}
\author{
Hongen Zhang
}
\examples{
# Load RCircos libaray
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

library(RCircos);


# Load human chromosome ideogram
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

data(UCSC.HG19.Human.CytoBandIdeogram);
cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;


# Set RCircso core components
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;
chr.exclude <- c("chrX", "chrY");
tracks.inside <- 5;
tracks.outside <- 0;
\dontrun{RCircos.Set.Core.Components(cyto.info,  
            chr.exclude, num.inside, num.outside); }



# Initialize graphic device and plot chromosome  
# ideogram
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

\dontrun{RCircos.Set.Plot.Area(); }
\dontrun{RCircos.Chromosome.Ideogram.Plot(); }


# Load tile data and plot one track of tiles  
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

data(RCircos.Tile.Data);
track.num <- 4;
side <- "in";
\dontrun{RCircos.Tile.Plot(RCircos.Tile.Data,  
                         track.num, side); }
}
\keyword{ methods }

