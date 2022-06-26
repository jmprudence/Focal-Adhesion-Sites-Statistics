### Mechnobiology and Biomedicine Lab ### 
### Department of Biomedical Engineering ### 
### Wichita State University ### 
### GitHub: @mechanobiology ### 
### Dr. David Long & Jess Prudence ###

### Purpose:
#Tests if point patterns (focal adhesion sites) from two different (cell) populations are 
#statistically different using a studentized permutation test of Hahn (2012). 
#This specific iteration simulates a machine learning estimate having location errors


### Key Reference(s):
# Hahn, U., (2012) Journal of the American Statistical Association. Vol. 107, No. 498, 
# "A Studentized Permutation Test for the Comparison of Spatial Point Patterns", pp 754-764.

### Versions:
#May 2020       - original version
#June 2021      - will read in dimensional data and normalize, will check if membrane 
#                 data is order clockwise or anti-clockwise. 
#14 August 2021 - Added functionality to determine the concave hull of the membrane
#                 points, which will be the outside, only keeps outside points.
#14 August 2020 - Added functionality to interpolate between membrane points to add any
#                 number of points (to increase the overall number of membrane points)
#27 August 2021 - Added Traveling Salesman Problem to reorder the membrane (x,y) points, then
#                 then will check if those ordered points are clockwise or anticlockwise before 
#                 making window owin function in spatstat
#27 Feb    2022 - Refined comments, eliminated unused functions, included population as well 
#                 as ml and actual in the comments to reflect testing
#15 June   2022 - Added simulation to locations being incorrect using rjitter 

### Basic Overview-Steps:
#1.  Reads in FA and membrane (x,y) data and group (actual or ml) from a .csv file
#2.  The .csv header format is the following: xfa yfa xmem ymem group
#3.  Each cell has its own .csv file, a specific filename structure is not required because
###  the program looks just for all .csv files in folder
#4.  then, a list of all .csv files in folder are created by reading in each .csv file in file_list and 
###  creates a data frame with the same name as the .csv file
#5.  Normalize all membrane and focal adhesion data (if desired). We may not use in future
#6.  Find the concave hull of the membrane data, so only outside data points will be included in the window
#7.  (optional) Interpolate between membrane points and add data point to smooth membrane outline
#8.  Use the Traveling Salesman Problem to order membrane points so that points go in order
###  around membrane
#9.  Use function 'clockwise' to determine if membrane points are now ordered
###  clockwise or anti-clockwise, if clockwise will reverse before building window
#10. A ppp list is built for each cell, 
#11. Next, a hyperframe is built
#12. Next, Hahn's studentized permutation is used to compare two populations

### Notes:
#If needed, we can use 'ondemand' and run RStudio on BeoShock, the program has been checked on BeoShock (just need to change 'pwd' and 'folder')
#We can quantify if the FA distributions are random
#We could also use marks to include, in addition to FA centroid things like FA area, FA orientation

### Install libraries: 
#if needed, install libraries by uncommenting lines below. Only need to install once. 
# install.packages("spatstat")
# install.packages("readr")
# install.packages("contoureR")
# install.packages("sf")
# install.packages("smoothr")
# install.packages("concaveman")
# install.packages("TSP")
# install.packages("tidyverse")

### Start libraries:  
library(spatstat)
library(readr)
library(contoureR)
library(sf)
library(smoothr)
library(concaveman)
library(TSP)
library(tidyverse)

### Start Function Section:

#' clockwise function
#' @title Checks whether points for an owin are clockwise
#' @param x a dataframe with x coordinates in the first column and y coordinates in the second. 
#' @details Similarly to owin, the polygon should not be closed
#' @return A logical telling whether the polygon is arranged clockwise.
#' @author The idea has been scavenged from https://stackoverflow.com/a/1165943/1082004

clockwise <- function(x) {
  
  x.coords <- c(x[[1]], x[[1]][1])
  y.coords <- c(x[[2]], x[[2]][1])
  
  double.area <- sum(sapply(2:length(x.coords), function(i) {
    (x.coords[i] - x.coords[i-1])*(y.coords[i] + y.coords[i-1])
  }))
  
  double.area > 0
} 


#### Testing Random Sampling 
# from Bruno Silva: https://www.r-bloggers.com/2019/04/random-sampling-of-files/  

# path = path to folder with files to select                                                                                                            #
# percent_number = percentage or number of recordings to select. 
  # If value is between 0 and 1 percentage of files is assumed, if value 
  # greater than 1, number of files is assumed 

results_folder <- # Enter path to results folder
reps <- 30 #readline(prompt="Enter Repetitions: ")
reps <- as.numeric(reps)
results <- replicate(reps, NULL)
percent_number_big <- 390 #readline(prompt="Enter Number of Large Group Files:")
percent_number_small <- 3 #readline(prompt="Enter Number of Small Group Files:")
Radius <- readline(prompt ="Enter Radius of Permutation (image is normalized):")
Radius <- as.numeric(Radius)


# Progress Bar to see % of iterations completed  
pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = reps,   # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "*")   # Character used to create the bar


#Starts the working section
for(r in 1:reps){
  
  #Clears out folder to start random sample with a clear folder
  unlink# (Enter path to random sampling folder, recursive=TRUE)
 
#Function to Randomly Sample Population 
  random_files_big <- function(path, percent_number_big, pattern = "*.csv"){ 
                                                
    # Get file list with full path and file names
    files <- list.files(path, full.names = TRUE, pattern = pattern)
    file_names <- list.files(path, pattern = pattern)
    
    # Select the desired % or number of file by simple random sampling
    randomize <- sample(seq(files))
    files2analyse <- files[randomize]
    names2analyse <- file_names[randomize]
    if(percent_number_big <= 1){
      size <- floor(percent_number_big * length(files))
    }else{
      size <- percent_number_big
    }
    files2analyse <- files2analyse[(1:size)]
    names2analyse <- names2analyse[(1:size)]
    
    # Move files
    for(r in seq(files2analyse)){
      file.copy(from = files2analyse[r], to = paste0(results_folder, "/", names2analyse[r]) )
    }
  }

#Randomly Sample Large Group Cells
path <- # Enter path to large group files
percent_number_big <- as.numeric(percent_number_big)
random_files_big(path, percent_number_big, pattern = "*.csv")

#Function to Randomly Sample Machine Learning Group
random_files_small <- function(path, percent_number_small, pattern = "*.csv"){
  
  # Get file list with full path and file names
  files <- list.files(path, full.names = TRUE, pattern = pattern)
  file_names <- list.files(path, pattern = pattern)
  
  # Select the desired % or number of file by simple random sampling
  randomize <- sample(seq(files))
  files2analyse <- files[randomize]
  names2analyse <- file_names[randomize]
  if(percent_number_small <= 1){
    size <- floor(percent_number_small * length(files))
  }else{
    size <- percent_number_small
  }
  files2analyse <- files2analyse[(1:size)]
  names2analyse <- names2analyse[(1:size)]
  
  # Move files
  for(r in seq(files2analyse)){
    file.copy(from = files2analyse[r], to = paste0(results_folder, "/", names2analyse[r]) )
  }
}

#Randomly Sample Small Group Cells
path <- # Enter path to small group files
percent_number_small <- as.numeric(percent_number_small)
random_files_small(path, percent_number_small, pattern = "*.csv")

### End Function Section

# set current working directory to read/write data
pwd=setwd#(Enter path to working directory)

# set path to directory that with multiple .csv files of input data
folder <- (results_folder)

# create list of all .csv files in folder, file name prefix does not matter, will read
# in all .csv files in 'folder'
file_list <- list.files (path=folder, pattern="*.csv")

# the total number of cells (the total in both populations)
ncell <- length(file_list)

# create lists to store focal adhesion points and group
falist <- vector(mode = "list", length = ncell) #list of focal adhesion points
glist <- vector(mode = "list", length = ncell) #list of group names for cells

# read in each .csv file in file_list and create a data frame with the same name as the .csv file
# note depending on how membrane data is generated, you may need to use the rev() function if 
# if you get an error about negative area for the owin. owin is expecting the points to be ordered anti-clockwise

### Primary Variables:
#mx,my         = x and y points of membrane
#fx,fy         = x and y points of focal adhesion centroids
#xmax, ymax    = maximum x and y values on membrane outline, respectively
#xmin, ymin    = minimum x and y values on membrane outline, respectively
#w             = window point within
#falist        = 2D point pattern dataset of focal adhesion sites
#glist         = list of groups i.e., does pattern belong to population cell, actual cell, machine learning cell
#op            = a list of indexes of the supplied membrane (x,y) data arrange in 'anti-clockwise' order
#mdf           = membrane data frame of (x,y) data in anti-clockwise order
#num           = number of interpolation points to add
#vi, wi        = empty vector to append the interpolated membrane
#re_ordered_xy = reordered membrane points by Traveling Salesman Problem
#result        = summary statistics for studentized permutation test

### Start Program:
  
for (i in 1:length(file_list)){
  assign(file_list[i], read.csv(paste(folder, file_list[i], sep=''), colClasses=c('numeric', 'numeric', 'numeric', 'numeric', 'character')))
  mx <- as.list(eval(parse(text = file_list[i]))[3])
  my <- as.list(eval(parse(text = file_list[i]))[4])
  #Find min and max x and y for membrane outline
  xmax <- max(mx[[1]]) 
  ymax <- max(my[[1]]) 
  xmin <- min(mx[[1]]) 
  ymin <- min(my[[1]]) 
  #normalize membrane x-y points
  mx <- lapply(mx, function(x){(x-xmin)/(xmax-xmin)})
  my <- lapply(my, function(y){(y-ymin)/(ymax-ymin)})
  mdf <- data.frame(mx,my)
  
  #determine concave hull to find outer border of membrane points
  #the concave membrane points will be saved to the variable mdf
  mdf <-data.frame(concaveman(data.matrix(mdf)))
  
  #increase the number of membrane points by linear interpolation using approx,
  #so point pattern window, w, will be smooth. This may not be necessary. 
  ###
  #number of interpolation points between to membrane points
  num <- 1
  
  #create empty vector to append the interpolated membrane (x,y) values
  vi <- vector()
  wi <- vector()
  
  v<- mdf[[1]] 
  w<- mdf[[2]]
  
  for (j in 1:(length(v)-1)){
    #two x and y value of membrane data points
    vsub<-v[j:(j+1)]
    wsub<-w[j:(j+1)]
    
    if(vsub[[1]]!=vsub[[2]]){
      #if membrane x-values are not equal, linearly interpolate between to points, 
      #then append
      vsp<-approx(vsub, wsub, method="linear", n = num)
      vi<-c(vi, vsp[[1]])
      wi<-c(wi, vsp[[2]])
    } else {
      #if values are equal, no interpolation, just append and go to next pair
      vi<-c(vi, vsub[[1]],vsub[[2]])
      wi<-c(wi, wsub[[1]],wsub[[2]])
    }
  }
  
  df <- data.frame(vi, wi)
  
  #Traveling Salesman Problem (TSP) - to order membrane points so points are in order as they
  #would be if the membrane is traced, which is shortest path.  
  xytsp <- ETSP(df)
  colnames(xytsp) <- c("x", "y")
  xytour <- solve_TSP(xytsp)
  # plot(xytsp, xytour, pch=20, tour_col="red", tour_lty="solid")
  re_ordered_xy <- df[xytour, ]
  
  #Now that we have the points ordered using the TSP, check to see if clockwise, or anticlockwise. 
  #If clockwise, reverse using 'rev'
  if(clockwise(re_ordered_xy)){
    w <- owin(poly=list(x=rev(re_ordered_xy[[1]]), y=rev(re_ordered_xy[[2]])))
  }else{
    w <- owin(poly=list(x=re_ordered_xy[[1]], y=re_ordered_xy[[2]]))
  }
  
  #Remove NaN in focal adhesion data, then normalize, then build a ppp
  #Note: NaN will occur since the number of focal adhesion points is less
  #      than the number of membrane points i.e., .csv columns different lengths
  fx <- na.omit(eval(parse(text = file_list[i]))[1])
  fy <- na.omit(eval(parse(text = file_list[i]))[2])
  fx <- lapply(fx, function(x){(x-xmin)/(xmax-xmin)}) #normalize fa centroids
  fy <- lapply(fy, function(y){(y-ymin)/(ymax-ymin)})
  
  group <- as.list(eval(parse(text= file_list[i]))[5])
  
  falist[[i]] <- ppp(fx[[1]], fy[[1]], window = w)
  
  # if the group is not population, FA sites are moved a maximum raidus based 
  # on user input. If point falls outside of window, shift will be retired
  # up to 100 times
   if(group != "population"){
    falist[[i]] <- rjitter(falist[[i]], radius = Radius, retry=TRUE, 
                           giveup=100,nsim=1, drop = TRUE)
  }else{    }
  
  
  #build group list of which group each data set belongs to
  if(i==1){glist <- na.omit(eval(parse(text = file_list[1]))[1,5])
  }else{
    glist <- c(glist, na.omit(eval(parse(text = file_list[i]))[1,5]))
  }
}


# build a hyperframe and assign each cell its group
focal_adhesion = hyperframe(FA_centroids=falist, group=factor(glist))


#run studentized permutation test from Hahn (2012) for differences in spatial point patterns

results[[r]] <- studpermu.test(focal_adhesion, FA_centroids ~ group, nperm=99, minpoints = 1)

setTxtProgressBar(pb, r)

}

close(pb)
print(Radius)

for(i in 1:reps){
  print(results[[i]])
}