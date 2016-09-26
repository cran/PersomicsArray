\name{spot_id}
\Rdversion{1.0}
\alias{spot_id}
%\alias{PersomicsArray}

\title{spot_id}

\description{This function sequentially reads multi-channel jpeg or tiff high content microscopy image files, identifies coordinates of printed siRNAs via one of the image channels, and then exports individual multi-channel tiff images for each identified array position.  The names of the exported images are supplied by an input csv file containing the array annotations.  Images are exported to automatically-created sub-directories named after the input image files.}

\usage{spot_id(files, annotation, channel.num=3, spot.channel=1, smooth.cycle=4,
binary.cut= 0.3, channel.scaling=TRUE, scale.percentiles=c(0.01,0.99))}

\arguments{
    \item{files}{Vector containing the names of desired files to be read from the current working directory.  All array images must have invariable numbers of rows and columns.}
    \item{annotation}{Name of the csv file in the current working directly that contains array annotations.  File must have invariable number of rows and columns.}
    \item{channel.num}{Numeric value specifying the number of color channels contained by the input image.  jpeg images contain a max of 3.  tiff files can contain more than this.  The maximum number of channels used by this function is 6.}
    \item{spot.channel}{Numeric value specifying the channel number used for identification of grid array spots.}
    \item{smooth.cycle}{Numeric value specifying the number of iterations of smoothing operation applied to the extracted binary image prior to identification of pixel clusters.  Designed for elimination of rogue pixels from binary image.}
    \item{binary.cut}{Numeric value in the range of 0-1 defining the cut-off pixel intensity for defining the binary image.}
    \item{channel.scaling}{If TRUE, all color channels are scaled according to the values supplied to scale.percentiles.  This parameter scales the input data to optimize coverage of the dynamic range of displayed pixel intensity values.  If FALSE, input data is used "as is" with no scaling operation applied.}
    \item{scale.percentiles}{Vector of two numeric values in the range of 0-1 that specify the percentiles of pixel intensity values used to define min and max for individual channel signal.}
}

\author{
John Smestad
\cr
\email{smestad.john@gmail.com}
}

\seealso{
<\code{\link{PersomicsArray-package}}>
<\code{\link{plot_low_res}}>
<\code{\link{example_annotation}}>
<\code{\link{example_plate}}>
}

\examples{
# get names of all plate image files in current working directory
tif.files <- list.files(path = getwd(), pattern = ".tif", all.files = FALSE,
                    full.names = FALSE, recursive = FALSE, ignore.case = TRUE,
                    include.dirs = FALSE)

# read images, process, and extract individual image files for each array position
spot_id(files=tif.files,ann="annotation.csv",channel.num=3, spot.channel=1,
smooth.cycle=4, binary.cut= 0.3, channel.scaling=TRUE, scale.percentiles=c(0.01,0.99))

}


