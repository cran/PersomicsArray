# Copyright (C) 2016  John A. Smestad
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

plot_low_res <- function(plot.image=NULL, rescale.factor=1000, add=FALSE, pallete=gray.colors(20),
                         main=NULL){
  col.seq <- seq(from=1, to=ncol(plot.image), by=round(sqrt(rescale.factor)))
  row.seq <- seq(from=1, to=nrow(plot.image), by=round(sqrt(rescale.factor)))
  temp <- plot.image[rev(row.seq),]; temp <- temp[,col.seq]
  temp <- t(temp)
  # plot image
  old.par <- par()
  par(mar=c(1,1,1.5,0),omi=c(0,0,0,0))
  image(temp, col=pallete,asp=ncol(temp)/nrow(temp), add=add,xaxt="n",yaxt="n",main=main)
  par(mar=c(5.1,4.1,4.1,2.1),omi=c(0,0,0,0))
}

