\name{RCircos.Line.Plot}
\alias{RCircos.Line.Plot}
\title{
Circos-like Line Plot
}
\description{
Draw one track of line plot. RCircos core components and graphics device must bae initialized before drawing.
}
\usage{
RCircos.Line.Plot(line.data, data.col, track.num, side)
}
\arguments{
  \item{line.data}{
A data frame with the first three columns as chromosome names, start position and end position of each data point followed by data values for one or more samples.
}
  \item{data.col}{
Integer, representing the ordinal number of the column in input dataset that contains the data to be plotted.
}
  \item{track.num}{
Integer, representing the ordinal number of the plot track where the lines will be plotted.
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


# Load human chromosome ideogram and sample data
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

data(UCSC.HG19.Human.CytoBandIdeogram);
data(RCircos.Line.Data);


# Set RCircso core components
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;
chr.exclude <- c("chrX", "chrY");
num.inside <- 5;
num.outside <- 0;
\dontrun{RCircos.Set.Core.Components(cyto.info,  
             chr.exclude, num.inside, num.outside); }


# Initialize graphic device, plot chromosome ideogram  
# and one track of lines
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

\dontrun{RCircos.Set.Plot.Area();}
\dontrun{RCircos.Chromosome.Ideogram.Plot(); }

data.col <- 5; 
track.num <- 4;
side <- "in";
\dontrun{RCircos.Line.Plot(RCircos.Line.Data, data.col,  
                           track.num, side); }
}
\keyword{ methods }
