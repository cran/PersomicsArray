\name{plot_low_res}
\Rdversion{1.0}
\alias{plot_low_res}
%\alias{PersomicsArray}

\title{plot_low_res}

\description{This function generates low resolution plots of pixel values stored in arrays, and
is called by function spot_id.}

\usage{plot_low_res(plot.image=NULL, rescale.factor=1000, add=FALSE,
pallete=gray.colors(20),main=NULL)}

\arguments{
    \item{plot.image}{Array-type object containing image data.}
    \item{rescale.factor}{Numeric value specifying the fold reduction of resolution from original image.  Larger values make the function run faster, but at the cost of image resolution.}
    \item{add}{If TRUE, adds the plot to the current active device.  If FALSE, a new plot is generated.}
    \item{pallete}{Specifies the color pallete to be used for generating the plot.}
    \item{main}{Character string containing the plot title.}
}

\author{
John Smestad
\cr
\email{smestad.john@gmail.com}
}

\seealso{
<\code{\link{PersomicsArray-package}}>
<\code{\link{spot_id}}>
<\code{\link{example_annotation}}>
<\code{\link{example_plate}}>
}

\examples{
# load example image data from package
data("PersomicsPlate")

# isolate single color channel
red <- example_plate[,,1]

# plot low resolution image of
plot_low_res(red, rescale.factor=10)
}
