\name{PersomicsArray-package}
\Rdversion{1.0}
\alias{PersomicsArray-package}
%\alias{PersomicsArray}

\docType{package}

\title{PersomicsArray: an R package for automated extraction of annotated images from high-content microscopy images of Persomics plates}

\description{The PersomicsArray package contains functions to sequentially read multi-channel jpeg or tiff high content microscopy image files, identify coordinates of printed siRNAs via one of the image channels, and then export individual multi-channel tiff images for each identified array position.  The names of the exported images are supplied by an input csv file containing the array annotations.  Images are exported as tiff files to automatically-created sub-directories named after the input image files.
}

\details{
\tabular{ll}{
Package: \tab PersomicsArray\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-09-23\cr
License: \tab GPL-3\cr
LazyLoad: \tab yes\cr
}
This is version 1.0 of the PersomicsArray package.
}

\author{
John Smestad
\cr
\email{smestad.john@gmail.com}
}

\keyword{ package }
\seealso{
<\code{\link{spot_id}}>
<\code{\link{plot_low_res}}>
<\code{\link{example_annotation}}>
<\code{\link{example_plate}}>
}
