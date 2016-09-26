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


spot_id <- function(files,annotation,channel.num=3, spot.channel=1, smooth.cycle=4, binary.cut= 0.3, channel.scaling=TRUE, scale.percentiles=c(0.01,0.99)){
  if(is.character(files)==FALSE | is.character(annotation)==FALSE){
    print("error: files and annotation entries must be character strings or vectors of character strings")
  }
  if(is.character(files)==TRUE & is.character(annotation)==TRUE){
    wd <- getwd() # directory where original image files are located
    if(length(files)>0){
      for(i in 1:length(files)){
        file <- files[i]
        line <- str_c("begin analysis. ","file= ",file); print(line)
        # read annotation file
        ann <- read.table(annotation, sep=",",stringsAsFactors = TRUE)
        grid.rows <- nrow(ann); grid.columns <- ncol(ann) # extract grid dimensions from annotation file
        # image type from image name
        str <- strsplit(file,"[.]")
        # create name for output directory
        wd2 <- str_c(wd,"/",str[[1]][1])
        dir.create(wd2)
        wd3 <- str_c(wd2,"/tiff out")
        dir.create(wd3)
        type <- str[[1]][length(str[[1]])]
        # read image
        if(type %in% c("jpeg","jpg","jpe","jif","jfif","jfi")){
          img <- readJPEG(file, native=FALSE)
          print("read jpeg image...")
        }
        if(type %in% c("tif","tiff")){
          img <- readTIFF(file, native=FALSE, all=TRUE,as.is=FALSE)
          print("read tiff image...")
        }
        # calculate dimensions of image and create vectors of x and y coordinates
        if(type %in% c("jpeg","jpg","jpe","jif","jfif","jfi")){
          dimy1 <- seq(1:dim(img)[2])
          dimx1 <- seq(1:dim(img)[1])
        }
        if(type %in% c("tiff","tif")){
          dimy1 <- seq(1:dim(img[[1]])[2])
          dimx1 <- seq(1:dim(img[[1]])[1])
        }
        x <- rep(dimx1,length(dimy1))
        y <- rep(dimy1,each=max(dimx1))
        dimx <- length(dimx1)
        dimy <- length(dimy1)
        print("calculated image dimensions...")
        rm(dimx1); rm(dimy1); gc()

        # define list structures (1=file, 2=ch1, 3=ch2, 4=ch3, 5=ch4, 6=ch5, 7=ch6, 8=spot.ch.scaled-->binary
        ch <- list(matrix(1:dimx*dimy, ncol = dimx))
        if(channel.num==6){
          data <- list(file,ch,ch,ch,ch,ch,ch,ch)
          print(str_c("generated list structure for 6 channels..."))
        }
        if(channel.num==5){
          data <- list(file,ch,ch,ch,ch,ch,NA,ch)
          print(str_c("generated list structure for 5 channels..."))
        }
        if(channel.num==4){
          data <- list(file,ch,ch,ch,ch,NA,NA,ch)
          print(str_c("generated list structure for 4 channels..."))
        }
        if(channel.num==3){
          data <- list(file,ch,ch,ch,NA,NA,NA,ch)
          print(str_c("generated list structure for 3 channels..."))
        }
        if(channel.num==2){
          data <- list(file,ch,ch,NA,NA,NA,NA,ch)
          print(str_c("generated list structure for 2 channels..."))
        }
        if(channel.num==1){
          data <- list(file,ch,NA,NA,NA,NA,NA,ch)
          print(str_c("generated list structure for 1 channel..."))
        }
        rm(ch); gc()

        # store separate image channels as data frames in "data" list structure
        if(type %in% c("jpeg","jpg","jpe","jif","jfif","jfi")){
          for(j in 1:channel.num){
            data[[j+1]][[1]] <- img[,,j]
            print(str_c("channel split ",j," of ",channel.num,"..."))
          }
        }
        if(type %in% c("tiff","tif")){
          for(j in 1:channel.num){
            data[[j+1]][[1]] <- img[[j]]
            print(str_c("channel split ",j," of ",channel.num,"..."))
          }
        }
        # remove "img" object (input image file) and release the system memory
        rm(img); gc()

        # scale imported color channels.  use 1% and 99% pixel intensity percentiles to define min and max
        if(channel.scaling==TRUE){
          for(j in 1:channel.num){
            quar <- quantile(data[[j+1]][[1]], scale.percentiles)
            data[[j+1]][[1]] <- 1*(data[[j+1]][[1]]-as.numeric(quar[1]))/as.numeric(quar[2]-quar[1])
            data[[j+1]][[1]][data[[j+1]][[1]] < 0] <- 0; data[[j+1]][[1]][data[[j+1]][[1]] > 1] <- 1
            print(str_c("finished scaling channel ",j," of ",channel.num,"..."))
          }
        }

        # generate scaled (0-1) image of spot identification channel.  uses 5th and 95th pixel intensity percentiles to define min and max.
        quar <- quantile(data[[spot.channel+1]][[1]], c(0.05,0.95))
        data[[8]][[1]] <- 1*(data[[spot.channel+1]][[1]]-as.numeric(quar[1]))/as.numeric(quar[2]-quar[1])
        # redefine the 5% of data outside of the 5-95% percentiles as either 0 or 1
        data[[8]][[1]][data[[8]][[1]] < 0] <- 0; data[[8]][[1]][data[[8]][[1]] > 1] <- 1
        # convert to binary image for spot channel
        data[[8]][[1]][data[[8]][[1]] <= binary.cut] <- 0; data[[8]][[1]][data[[8]][[1]] > binary.cut] <- 1

        ##### smooth image by adjacent pixel number calcularion method using low-res image
        # first create a low res image
        rescale.factor <- 1000 # defines fold loss of resolution
        col.seq <- seq(from=1, to=ncol(data[[8]][[1]]), by=round(sqrt(rescale.factor)))
        row.seq <- seq(from=1, to=nrow(data[[8]][[1]]), by=round(sqrt(rescale.factor)))
        low.res.binary <- data[[8]][[1]][row.seq,]; low.res.binary <- low.res.binary[,col.seq]

        # calculate adjacent pixel numbers
        adj.px.low.res <- low.res.binary # copy data to new channel
        adj.px.low.res[adj.px.low.res >= 0] <- 0 # assign zero value to matrix to start with
        dimy1 <- seq(1:dim(adj.px.low.res)[2])
        dimx1 <- seq(1:dim(adj.px.low.res)[1])
        dimx <- length(dimx1)
        dimy <- length(dimy1)
        x <- rep(dimx1,length(dimy1))
        y <- rep(dimy1,each=max(dimx1))
        update <- seq(from=0, to=length(x), by=100000/rescale.factor)
        # erode and smooth image based on adjacent pixel numbers --> export images along the way to show progress
        setwd(wd2)
        tiff("1-scaled_binary.tiff",width = 1500, height = 1500, compression="none")
        plot_low_res(low.res.binary,rescale.factor = 1, main="Image 1: Scaled Binary")
        dev.off()
        if(smooth.cycle >= 1){
          # repeat for each smoothing cycle
          for(k in 1:smooth.cycle){
            # for each pixel in binary image, calculate number of adjacent pixels with signal
            for(j in 1:length(x)){
              if(j %in% update){
                print(str_c("smoothing round ",k," of ",smooth.cycle," ",100*round(j/length(x),digits=3),"% complete..."))
              }
              # for pixels not on an edge
              if(x[j]>1 & x[j]<dimx & y[j]>1 &y[j]<dimy){
                m <- low.res.binary[(x[j]-1):(x[j]+1),(y[j]-1):(y[j]+1)]
                m[2,2] <- NA
                adj.px.low.res[x[j],y[j]] <- sum(m == 1, na.rm=TRUE)
              }
              # for pixels on left edge
              if(x[j]==1 & y[j]!=1 & y[j]!=dimy){
                m <- low.res.binary[(x[j]):(x[j]+1),(y[j]-1):(y[j]+1)]
                m[1,2] <- NA
                adj.px.low.res[x[j],y[j]] <- sum(m == 1, na.rm=TRUE)
              }
              # for pixels on right edge
              if(x[j]==dimx & y[j]!=1 & y[j]!=dimy){
                m <- low.res.binary[(x[j]-1):(x[j]),(y[j]-1):(y[j]+1)]
                m[2,2] <- NA
                adj.px.low.res[x[j],y[j]] <- sum(m == 1, na.rm=TRUE)
              }
              # for pixels on top edge
              if(y[j]==1 & x[j]!=1 & x[j]!=dimx){
                m <- low.res.binary[(x[j]-1):(x[j]+1),(y[j]):(y[j]+1)]
                m[2,1] <- NA
                adj.px.low.res[x[j],y[j]] <- sum(m == 1, na.rm=TRUE)
              }
              # for pixels on bottom edge
              if(y[j]==dimy & x[j]!=1 & x[j]!=dimx){
                m <- low.res.binary[(x[j]-1):(x[j]+1),(y[j]-1):(y[j])]
                m[2,2] <- NA
                adj.px.low.res[x[j],y[j]] <- sum(m == 1, na.rm=TRUE)
              }
              # for pixels in top left
              if(x[j]==1 & y[j]==1){
                m <- low.res.binary[(x[j]):(x[j]+1),(y[j]):(y[j]+1)]
                m[1,1] <- NA
                adj.px.low.res[x[j],y[j]] <- sum(m == 1, na.rm=TRUE)
              }
              # for pixels in top right
              if(x[j]==dimx & y[j]==1){
                m <- low.res.binary[(x[j]-1):(x[j]),(y[j]):(y[j]+1)]
                m[2,1] <- NA
                adj.px.low.res[x[j],y[j]] <- sum(m == 1, na.rm=TRUE)
              }
              # for pixels in bottom left
              if(x[j]==1 & y[j]==dimy){
                m <- low.res.binary[(x[j]):(x[j]+1),(y[j]-1):(y[j])]
                m[1,2] <- NA
                adj.px.low.res[x[j],y[j]] <- sum(m == 1, na.rm=TRUE)
              }
              # for pixels in bottom right
              if(x[j]==dimx & y[j]==dimy){
                m <- low.res.binary[(x[j]-1):(x[j]),(y[j]-1):(y[j])]
                m[2,2] <- NA
                adj.px.low.res[x[j],y[j]] <- sum(m == 1, na.rm=TRUE)
              }
            }
            print("finished calculating adjacent pixel numbers ...")
            # smooth clusters based on adjacent pixels
            if(k %in% c(1)){
              low.res.binary[adj.px.low.res >= 7] <- 1; low.res.binary[adj.px.low.res <= 6] <- 0
            }
            if(k > 1){
              low.res.binary[adj.px.low.res >= 4] <- 1; low.res.binary[adj.px.low.res <= 3] <- 0
            }
            tiff(str_c(k+1,"-smooth_cycle_",k,".tiff"),width = 1500, height = 1500, compression="none")
            plot_low_res(low.res.binary,rescale.factor = 1, main=str_c("Image ",k+1,": Smooth Cycle ",k))
            dev.off()
            line2 <- str_c("finished smoothing operation ",k)
            print(line2)
          }
        }
        rm(adj.px.low.res); gc()

        # identify pixel clusters
        # assign unique numberings to pixel clusters
        clus.numer <- low.res.binary # copy data to new channel
        clus.numer[clus.numer >= 0] <- NA
        num <- 0
        for(j in 1:length(x)){
          if(j %in% update){
            print(str_c("pixel cluster calculation ",100*round(j/length(x),digits=3),"% complete..."))
          }
          # do this only for pixels that have signal in the binary image
          if(low.res.binary[x[j],y[j]] == 1){
            ## select adjacent pixels in cluster numbering matrix
            # for pixels not on an edge
            if(x[j]>1 & x[j]<dimx & y[j]>1 & y[j]<dimy){
              m <- clus.numer[(x[j]-1):(x[j]+1),(y[j]-1):(y[j]+1)]
            }
            # for pixels on left edge
            if(x[j]==1 & y[j]!=1 & y[j]!=dimy){
              m <- clus.numer[(x[j]):(x[j]+1),(y[j]-1):(y[j]+1)]
            }
            # for pixels on right edge
            if(x[j]==dimx & y[j]!=1 & y[j]!=dimy){
              m <- clus.numer[(x[j]-1):(x[j]),(y[j]-1):(y[j]+1)]
            }
            # for pixels on top edge
            if(y[j]==1 & x[j]!=1 & x[j]!=dimx){
              m <- clus.numer[(x[j]-1):(x[j]+1),(y[j]):(y[j]+1)]
            }
            # for pixels on bottom edge
            if(y[j]==dimy & x[j]!=1 & x[j]!=dimx){
              m <- clus.numer[(x[j]-1):(x[j]+1),(y[j]-1):(y[j])]
            }
            # for pixels in top left
            if(x[j]==1 & y[j]==1){
              m <- clus.numer[(x[j]):(x[j]+1),(y[j]):(y[j]+1)]
            }
            # for pixels in top right
            if(x[j]==dimx & y[j]==1){
              m <- clus.numer[(x[j]-1):(x[j]),(y[j]):(y[j]+1)]
            }
            # for pixels in bottom left
            if(x[j]==1 & y[j]==dimy){
              m <- clus.numer[(x[j]):(x[j]+1),(y[j]-1):(y[j])]
            }
            # for pixels in bottom right
            if(x[j]==dimx & y[j]==dimy){
              m <- clus.numer[(x[j]-1):(x[j]),(y[j]-1):(y[j])]
            }
            ## determine cluster number to assign to pixel, if in a cluster
            m2 <- m[is.na(m) == FALSE]
            # do this if a pixel belongs to a new cluster
            if(length(unique(m2)) == 0){
              clus.numer[x[j],y[j]] <- num+1
              num <- num+1
            }
            # do this if a pixel belongs to an existing cluster
            if(length(unique(m2)) == 1){
              clus.numer[x[j],y[j]] <- m2[1]
            }
            # do this if a pixel borders 2 or more existing numbered clusters --> merge all into single cluster
            if(length(unique(m2)) >= 2){
              clus.numer[x[j],y[j]] <- min(m2)
              for(k in unique(m2)){
                clus.numer[clus.numer == k] <- min(m2)
              }
            }
          }
        }
        print(str_c((length(unique(as.vector(clus.numer)))-1)," spots identified..."))
        tiff(str_c(smooth.cycle+2,"-identified_pixel_clusters.tiff"),width = 1500, height = 1500, compression="none")
        plot_low_res(clus.numer, rescale.factor=1, pallete=sample(rainbow(num+1),replace=FALSE), main=str_c("Image ",smooth.cycle+2,": Identified Pixel Clusters"))
        dev.off()
        clus.numer[is.na(clus.numer)] <- 0 # replace NA with 0

        #### use pixel cluster positions as starting point to calculate midpoints for serial image export
        # identify x-y positions of pixel cluster midpoints
        cnum <- unique(as.vector(clus.numer)); cnum <- cnum[cnum !=0]
        rpos <- cpos <- NULL
        for(j in cnum){
          mat1 <- which(clus.numer==j, arr.ind=TRUE)
          rpos <- c(rpos, round(mean(mat1[,1]))); cpos <- c(cpos, round(mean(mat1[,2])))
        }

        # scale rpos and cpos to match coordinates of input images
        rpos2 <- rpos*nrow(data[[8]][[1]])/nrow(clus.numer); cpos2 <- cpos*ncol(data[[8]][[1]])/ncol(clus.numer)
        rm(clus.numer); gc() # remove clus.num and release system memory

        # generate guesses for spot locations (based on max/min of rpos2 and cpos2)
        xint <- (max(cpos2)-min(cpos2))/(grid.columns-1)
        yint <- (max(rpos2)-min(rpos2))/(grid.rows-1)

        seqx <- round(seq(from=min(cpos2), to=max(cpos2), by=xint))
        seqy <- round(seq(from=min(rpos2), to=max(rpos2), by=yint))
        xpts <- rep(seqx,each=grid.rows)
        ypts <- rep(seqy, grid.columns)

        # search rpos2 and cpos2 against ypts and xpts, respectively
        # keep only points nearest to grid guesses (remove extra features that may have been identified)
        df1 <- data.frame(rpos2,cpos2,"rpos.sq"=rep(NA,length(rpos2)),"cpos.sq"=rep(NA,length(rpos2)),"dist"=rep(NA,length(rpos2)))
        rpos2 <- NULL; cpos2 <- NULL
        for(j in 1:length(xpts)){
          x1 <- xpts[j]; y1 <- ypts[j]
          rpos.sq <- (df1$rpos2-y1)^2; cpos.sq <- (df1$cpos2-x1)^2
          df1$rpos.sq <- rpos.sq; df1$cpos.sq <- cpos.sq
          dist <- sqrt(rpos.sq+cpos.sq); df1$dist <- dist
          df2 <- df1[order(df1$dist),]
          rpos2 <- c(rpos2, df2[1,]$rpos2); cpos2 <- c(cpos2, df2[1,]$cpos2)
        }
        print("eliminated rogue pixel clusters from image data...")

        # run optimization algarithm to find optimal scalilng/offset parameters to make grid match
        # positions of detected spots
        # first optimize for x-position (corresponds to cpos2)
        plotFUN_A <- function(xpts,s){
          A <- s[1]
          a <- s[2]

          f <- A+a*xpts
          f
        }
        fitFUN_A <- function(s){
          xpts <- xpts
          A <- s[1]
          a <- s[2]

          f <- A+a*xpts

          lsse <- sum((f-cpos2)^2)
          lsse
        }
        RsqrFUN_A <- function(s){
          xpts <- xpts
          A <- s[1]
          a <- s[2]

          f <- A+a*xpts

          lsse <- sum((f-cpos2)^2)
          ### also called SSerr
          SSerr <- lsse
          SStot <- sum((mean(cpos2)-cpos2)^2)
          R2 <- 1-SSerr/SStot
          R2
        }
        #linear regression for native fold slope and intercept
        npar <- 2
        start <- c(1,1)
        out <- NULL
        out2 <- NULL
        ####
        nlminbFit <- nlminb(start,fitFUN_A,lower=c(-10000,0),upper=c(30000,20))
        out$nlminb=nlminbFit$par
        out2$nlminb=nlminbFit$objective
        ##print results
        best_A <- nlminbFit$par
        print("x dimension fit:")
        print(str_c("offset= ",best_A[1],", scaling= ",best_A[2]))
        print(str_c("R^2= ",RsqrFUN_A(best_A)))
        xpts2 <- best_A[1]+best_A[2]*xpts

        # second, optimize for y-position (corresponds to rpos2)
        plotFUN_B <- function(ypts,s){
          A <- s[1]
          a <- s[2]

          f <- A+a*sort(ypts)
          f
        }
        fitFUN_B <- function(s){
          ypts <- ypts
          A <- s[1]
          a <- s[2]

          f <- A+a*sort(ypts)

          lsse <- sum((f-sort(rpos2))^2)
          lsse
        }
        RsqrFUN_B <- function(s){
          ypts <- ypts
          A <- s[1]
          a <- s[2]

          f <- A+a*sort(ypts)

          lsse <- sum((f-sort(rpos2))^2)
          ### also called SSerr
          SSerr <- lsse
          SStot <- sum((mean(sort(rpos2))-sort(rpos2))^2)
          R2 <- 1-SSerr/SStot
          R2
        }
        #linear regression for native fold slope and intercept
        npar <- 2
        start <- c(0,1)
        out <- NULL
        out2 <- NULL
        ####
        nlminbFit <- nlminb(start,fitFUN_B,lower=c(-3000,0.7),upper=c(3000,3))
        out$nlminb=nlminbFit$par
        out2$nlminb=nlminbFit$objective
        ##print results
        best_B <- nlminbFit$par
        print("y dimension fit:")
        print(str_c("offset= ",best_B[1],", scaling= ",best_B[2]))
        print(str_c("R^2= ",RsqrFUN_B(best_B)))
        ypts2 <- best_B[1]+best_B[2]*ypts

        # fine tune individual column (y) position fits (adjust offset)
        for(j in 1:grid.columns){
          ypts <- ypts2[(grid.rows*(j-1)+1):(grid.rows*j)] # select single column from grid fit
          yspots <- rpos2[(grid.rows*(j-1)+1):(grid.rows*j)]

          plotFUN_B <- function(ypts,s){
            A <- s[1]

            f <- A+sort(ypts)
            f
          }
          fitFUN_B <- function(s){
            ypts <- ypts
            A <- s[1]

            f <- A+sort(ypts)
            f

            lsse <- sum((f-sort(yspots))^2)
            lsse
          }
          RsqrFUN_B <- function(s){
            ypts <- ypts
            A <- s[1]

            f <- A+sort(ypts)
            f

            lsse <- sum((f-sort(yspots))^2)
            ### also called SSerr
            SSerr <- lsse
            SStot <- sum((mean(sort(yspots))-sort(yspots))^2)
            R2 <- 1-SSerr/SStot
            R2
          }
          #linear regression for native fold slope and intercept
          npar <- 1
          start <- c(0)
          out <- NULL
          out2 <- NULL
          ####
          nlminbFit <- nlminb(start,fitFUN_B,lower=c(-3000),upper=c(3000))
          out$nlminb=nlminbFit$par
          out2$nlminb=nlminbFit$objective
          ##print results
          best_B <- nlminbFit$par
          print("y dimension fit:")
          print(str_c("offset= ",best_B[1]))
          print(str_c("R^2= ",RsqrFUN_B(best_B)))
          # store new values in ypts2
          ypts2[(grid.rows*(j-1)+1):(grid.rows*j)] <- best_B[1]+ypts
        }

        # format spot positions into grid that corresponds to annotation file
        for(j in 1:grid.columns){
          xp <- xpts2[(grid.rows*(j-1)+1):(grid.rows*j)]
          yp <- ypts2[(grid.rows*(j-1)+1):(grid.rows*j)]
          if(j==1){
            x.grid <- data.frame(xp)
            y.grid <- data.frame(yp)
          }
          if(j>1){
            x.grid <- cbind(x.grid,data.frame(xp))
            y.grid <- cbind(y.grid,data.frame(yp))
          }
        }

        # create vector of annotations from "ann" (data frame)
        for(j in 1:grid.columns){
          if(j==1){
            spot.ids <- as.vector(ann[,j])
          }
          if(j>1){
            spot.ids <- c(spot.ids,as.vector(ann[,j]))
          }
        }
        # export tiff images of single positions
        cut <- 450 # define how many pixels from the center to include in each image
        setwd(wd3)
        len <- nrow(ann)*ncol(ann)
        for(j in 1:len){
          if(xpts2[j] > cut & xpts2[j] < ncol(data[[8]][[1]])-cut & ypts2[j] > cut & ypts2[j] < nrow(data[[8]][[1]])-cut){
            mat1 <- data[[2]][[1]][(ypts2[j]-cut):(ypts2[j]+cut),(xpts2[j]-cut):(xpts2[j]+cut)]
            if(channel.num>1){
              mat2 <- data[[3]][[1]][(ypts2[j]-cut):(ypts2[j]+cut),(xpts2[j]-cut):(xpts2[j]+cut)]
            }
            if(channel.num>2){
              mat3 <- data[[4]][[1]][(ypts2[j]-cut):(ypts2[j]+cut),(xpts2[j]-cut):(xpts2[j]+cut)]
            }
            if(channel.num>3){
              mat4 <- data[[5]][[1]][(ypts2[j]-cut):(ypts2[j]+cut),(xpts2[j]-cut):(xpts2[j]+cut)]
            }
            if(channel.num>4){
              mat5 <- data[[6]][[1]][(ypts2[j]-cut):(ypts2[j]+cut),(xpts2[j]-cut):(xpts2[j]+cut)]
            }
            if(channel.num>5){
              mat6 <- data[[7]][[1]][(ypts2[j]-cut):(ypts2[j]+cut),(xpts2[j]-cut):(xpts2[j]+cut)]
            }
          }
          if(xpts2[j] < cut & xpts2[j] < ncol(data[[8]][[1]])-cut & ypts2[j] > cut & ypts2[j] < nrow(data[[8]][[1]])-cut){
            mat1 <- data[[2]][[1]][(ypts2[j]-cut):(ypts2[j]+cut),0:(xpts2[j]+cut)]
            if(channel.num>1){
              mat2 <- data[[3]][[1]][(ypts2[j]-cut):(ypts2[j]+cut),0:(xpts2[j]+cut)]
            }
            if(channel.num>2){
              mat3 <- data[[4]][[1]][(ypts2[j]-cut):(ypts2[j]+cut),0:(xpts2[j]+cut)]
            }
            if(channel.num>3){
              mat4 <- data[[5]][[1]][(ypts2[j]-cut):(ypts2[j]+cut),0:(xpts2[j]+cut)]
            }
            if(channel.num>4){
              mat5 <- data[[6]][[1]][(ypts2[j]-cut):(ypts2[j]+cut),0:(xpts2[j]+cut)]
            }
            if(channel.num>5){
              mat6 <- data[[7]][[1]][(ypts2[j]-cut):(ypts2[j]+cut),0:(xpts2[j]+cut)]
            }
          }
          if(xpts2[j] > cut & xpts2[j] < ncol(data[[8]][[1]])-cut & ypts2[j] < cut & ypts2[j] < nrow(data[[8]][[1]])-cut){
            mat1 <- data[[2]][[1]][0:(ypts2[j]+cut),(xpts2[j]-cut):(xpts2[j]+cut)]
            if(channel.num>1){
              mat2 <- data[[3]][[1]][0:(ypts2[j]+cut),(xpts2[j]-cut):(xpts2[j]+cut)]
            }
            if(channel.num>2){
              mat3 <- data[[4]][[1]][0:(ypts2[j]+cut),(xpts2[j]-cut):(xpts2[j]+cut)]
            }
            if(channel.num>3){
              mat4 <- data[[5]][[1]][0:(ypts2[j]+cut),(xpts2[j]-cut):(xpts2[j]+cut)]
            }
            if(channel.num>4){
              mat5 <- data[[6]][[1]][0:(ypts2[j]+cut),(xpts2[j]-cut):(xpts2[j]+cut)]
            }
            if(channel.num>5){
              mat6 <- data[[7]][[1]][0:(ypts2[j]+cut),(xpts2[j]-cut):(xpts2[j]+cut)]
            }
          }
          if(xpts2[j] < cut & xpts2[j] < ncol(data[[8]][[1]])-cut & ypts2[j] < cut & ypts2[j] < nrow(data[[8]][[1]])-cut){
            mat1 <- data[[2]][[1]][0:(ypts2[j]+cut),0:(xpts2[j]+cut)]
            if(channel.num>1){
              mat2 <- data[[3]][[1]][0:(ypts2[j]+cut),0:(xpts2[j]+cut)]
            }
            if(channel.num>2){
              mat3 <- data[[4]][[1]][0:(ypts2[j]+cut),0:(xpts2[j]+cut)]
            }
            if(channel.num>3){
              mat4 <- data[[5]][[1]][0:(ypts2[j]+cut),0:(xpts2[j]+cut)]
            }
            if(channel.num>4){
              mat5 <- data[[6]][[1]][0:(ypts2[j]+cut),0:(xpts2[j]+cut)]
            }
            if(channel.num>5){
              mat6 <- data[[7]][[1]][0:(ypts2[j]+cut),0:(xpts2[j]+cut)]
            }
          }
          if(xpts2[j] < cut & xpts2[j] < ncol(data[[8]][[1]])-cut & ypts2[j] > cut & ypts2[j] > nrow(data[[8]][[1]])-cut){
            mat1 <- data[[2]][[1]][(ypts2[j]-cut):nrow(data[[8]][[1]]),0:(xpts2[j]+cut)]
            if(channel.num>1){
              mat2 <- data[[3]][[1]][(ypts2[j]-cut):nrow(data[[8]][[1]]),0:(xpts2[j]+cut)]
            }
            if(channel.num>2){
              mat3 <- data[[4]][[1]][(ypts2[j]-cut):nrow(data[[8]][[1]]),0:(xpts2[j]+cut)]
            }
            if(channel.num>3){
              mat4 <- data[[5]][[1]][(ypts2[j]-cut):nrow(data[[8]][[1]]),0:(xpts2[j]+cut)]
            }
            if(channel.num>4){
              mat5 <- data[[6]][[1]][(ypts2[j]-cut):nrow(data[[8]][[1]]),0:(xpts2[j]+cut)]
            }
            if(channel.num>5){
              mat6 <- data[[7]][[1]][(ypts2[j]-cut):nrow(data[[8]][[1]]),0:(xpts2[j]+cut)]
            }
          }
          if(xpts2[j] > cut & xpts2[j] < ncol(data[[8]][[1]])-cut & ypts2[j] > cut & ypts2[j] > nrow(data[[8]][[1]])-cut){
            mat1 <- data[[2]][[1]][(ypts2[j]-cut):nrow(data[[8]][[1]]),(xpts2[j]-cut):(xpts2[j]+cut)]
            if(channel.num>1){
              mat2 <- data[[3]][[1]][(ypts2[j]-cut):nrow(data[[8]][[1]]),(xpts2[j]-cut):(xpts2[j]+cut)]
            }
            if(channel.num>2){
              mat3 <- data[[4]][[1]][(ypts2[j]-cut):nrow(data[[8]][[1]]),(xpts2[j]-cut):(xpts2[j]+cut)]
            }
            if(channel.num>3){
              mat4 <- data[[5]][[1]][(ypts2[j]-cut):nrow(data[[8]][[1]]),(xpts2[j]-cut):(xpts2[j]+cut)]
            }
            if(channel.num>4){
              mat5 <- data[[6]][[1]][(ypts2[j]-cut):nrow(data[[8]][[1]]),(xpts2[j]-cut):(xpts2[j]+cut)]
            }
            if(channel.num>5){
              mat6 <- data[[7]][[1]][(ypts2[j]-cut):nrow(data[[8]][[1]]),(xpts2[j]-cut):(xpts2[j]+cut)]
            }
          }
          if(xpts2[j] > cut & xpts2[j] > ncol(data[[8]][[1]])-cut & ypts2[j] > cut & ypts2[j] < nrow(data[[8]][[1]])-cut){
            mat1 <- data[[2]][[1]][(ypts2[j]-cut):(ypts2[j]+cut),(xpts2[j]-cut):ncol(data[[8]][[1]])]
            if(channel.num>1){
              mat2 <- data[[3]][[1]][(ypts2[j]-cut):(ypts2[j]+cut),(xpts2[j]-cut):ncol(data[[8]][[1]])]
            }
            if(channel.num>2){
              mat3 <- data[[4]][[1]][(ypts2[j]-cut):(ypts2[j]+cut),(xpts2[j]-cut):ncol(data[[8]][[1]])]
            }
            if(channel.num>3){
              mat4 <- data[[5]][[1]][(ypts2[j]-cut):(ypts2[j]+cut),(xpts2[j]-cut):ncol(data[[8]][[1]])]
            }
            if(channel.num>4){
              mat5 <- data[[6]][[1]][(ypts2[j]-cut):(ypts2[j]+cut),(xpts2[j]-cut):ncol(data[[8]][[1]])]
            }
            if(channel.num>5){
              mat6 <- data[[7]][[1]][(ypts2[j]-cut):(ypts2[j]+cut),(xpts2[j]-cut):ncol(data[[8]][[1]])]
            }
          }
          if(xpts2[j] > cut & xpts2[j] > ncol(data[[8]][[1]])-cut & ypts2[j] > cut & ypts2[j] > nrow(data[[8]][[1]])-cut){
            mat1 <- data[[2]][[1]][(ypts2[j]-cut):nrow(data[[8]][[1]]),(xpts2[j]-cut):ncol(data[[8]][[1]])]
            if(channel.num>1){
              mat2 <- data[[3]][[1]][(ypts2[j]-cut):nrow(data[[8]][[1]]),(xpts2[j]-cut):ncol(data[[8]][[1]])]
            }
            if(channel.num>2){
              mat3 <- data[[4]][[1]][(ypts2[j]-cut):nrow(data[[8]][[1]]),(xpts2[j]-cut):ncol(data[[8]][[1]])]
            }
            if(channel.num>3){
              mat4 <- data[[5]][[1]][(ypts2[j]-cut):nrow(data[[8]][[1]]),(xpts2[j]-cut):ncol(data[[8]][[1]])]
            }
            if(channel.num>4){
              mat5 <- data[[6]][[1]][(ypts2[j]-cut):nrow(data[[8]][[1]]),(xpts2[j]-cut):ncol(data[[8]][[1]])]
            }
            if(channel.num>5){
              mat6 <- data[[7]][[1]][(ypts2[j]-cut):nrow(data[[8]][[1]]),(xpts2[j]-cut):ncol(data[[8]][[1]])]
            }
          }
          if(xpts2[j] > cut & xpts2[j] > ncol(data[[8]][[1]])-cut & ypts2[j] < cut & ypts2[j] < nrow(data[[8]][[1]])-cut){
            mat1 <- data[[2]][[1]][0:(ypts2[j]+cut),(xpts2[j]-cut):ncol(data[[8]][[1]])]
            if(channel.num>1){
              mat2 <- data[[3]][[1]][0:(ypts2[j]+cut),(xpts2[j]-cut):ncol(data[[8]][[1]])]
            }
            if(channel.num>2){
              mat3 <- data[[4]][[1]][0:(ypts2[j]+cut),(xpts2[j]-cut):ncol(data[[8]][[1]])]
            }
            if(channel.num>3){
              mat4 <- data[[5]][[1]][0:(ypts2[j]+cut),(xpts2[j]-cut):ncol(data[[8]][[1]])]
            }
            if(channel.num>4){
              mat5 <- data[[6]][[1]][0:(ypts2[j]+cut),(xpts2[j]-cut):ncol(data[[8]][[1]])]
            }
            if(channel.num>5){
              mat6 <- data[[7]][[1]][0:(ypts2[j]+cut),(xpts2[j]-cut):ncol(data[[8]][[1]])]
            }
          }
          # scale image channels
          mult <- 32767
          mat1 <- mat1*mult
          if(channel.num>1){mat2 <- mat2*mult}
          if(channel.num>2){mat3 <- mat3*mult}
          if(channel.num>3){mat4 <- mat4*mult}
          if(channel.num>4){mat5 <- mat5*mult}
          if(channel.num>5){mat6 <- mat6*mult}

          # create raster images from individual channels, if more than one channel is used
          img1 <- raster(mat1, xmn=0, xmx=1, ymn=0, ymx=1, crs=NA, template=NULL)
          if(channel.num>1){img2 <- raster(mat2, xmn=0, xmx=1, ymn=0, ymx=1, crs=NA, template=NULL)}
          if(channel.num>2){img3 <- raster(mat3, xmn=0, xmx=1, ymn=0, ymx=1, crs=NA, template=NULL)}
          if(channel.num>3){img4 <- raster(mat4, xmn=0, xmx=1, ymn=0, ymx=1, crs=NA, template=NULL)}
          if(channel.num>4){img5 <- raster(mat5, xmn=0, xmx=1, ymn=0, ymx=1, crs=NA, template=NULL)}
          if(channel.num>5){img6 <- raster(mat6, xmn=0, xmx=1, ymn=0, ymx=1, crs=NA, template=NULL)}
          # create RasterStack objects
          if(channel.num==1){img <- stack(img1,layers=NULL)}
          if(channel.num==2){img <- stack(img1,img2,layers=NULL)}
          if(channel.num==3){img <- stack(img1,img2,img3,layers=NULL)}
          if(channel.num==4){img <- stack(img1,img2,img3,img4,layers=NULL)}
          if(channel.num==5){img <- stack(img1,img2,img3,img4,img5,layers=NULL)}
          if(channel.num==6){img <- stack(img1,img2,img3,img4,img5,img6,layers=NULL)}

          # write tiff image to file
          writeRaster(img, str_c("image_",j,"_",spot.ids[j],".tif"), overwrite = TRUE,datatype='INT2S')
          print(str_c("image ",j," ",spot.ids[j]," exported..."))
        }

        # plot low-res image showing detected spots and annotations
        setwd(wd2)
        rescale.factor <- 1000 # defines fold loss of resolution
        col.seq <- seq(from=1, to=ncol(data[[8]][[1]]), by=round(sqrt(rescale.factor)))
        row.seq <- seq(from=1, to=nrow(data[[8]][[1]]), by=round(sqrt(rescale.factor)))
        low.res.binary <- data[[8]][[1]][row.seq,]; low.res.binary <- low.res.binary[,col.seq] # redefine low.res.binary (to get rid of smoothing)
        col.seq <- seq(from=1, to=ncol(low.res.binary), by=round(sqrt(1)))
        row.seq <- seq(from=1, to=nrow(low.res.binary), by=round(sqrt(1)))
        temp <- low.res.binary[rev(row.seq),]; temp <- temp[,col.seq]
        temp <- t(temp)
        tiff(str_c(smooth.cycle+3,"-identified_grid_positions.tiff"),width = 1500, height = 1500, compression="none")
        par(mar=c(0,0,1.5,0))
        image(temp, col=c("#E6E6E6","#ABABAB"),asp=ncol(temp)/nrow(temp), add=FALSE,xaxt="n",yaxt="n",
              main=str_c("Image ",smooth.cycle+3,": Identified Grid Positions"),bty="n")
        rpos.scaled <- 1-ypts2/nrow(data[[8]][[1]]); cpos.scaled <- xpts2/ncol(data[[8]][[1]]) # plot fitted grid
        points(cpos.scaled,rpos.scaled, pch=18, cex=1, col="red") # plot grid
        dev.off()
        par(mar=c(5.1,4.1,4.1,2.1))
        setwd(wd)
      }
    }
  }
}

