
################## fun ############

cat("Function arrangeOccu installed\n")
cat("Use: arrangeOccu(occData,effData,interval = 'week')\n")
cat("* occData must have 4 columns:\n")
cat("1. Site ID (camera, plot, etc.)\n")
cat("2. Date in date format\n")
cat("3. Abundance (if not availabe, repeat 1)\n")
cat("4. Species name (even if only one species)\n")
cat("* effData must have 3 columns:\n")
cat("1. Site ID (camera, plot, etc.)\n")
cat("2. Initial date in date format\n")
cat("3. Final date in date format\n")

arrangeOccu<-function(occData,effData=NULL,interval="days",plot=FALSE,occFun=min,simplify=TRUE){

  occData<-na.omit(occData)
  
  if(ncol(occData)<3){
    occData[,3]<-1
  }
  
  if(ncol(occData)<4){
    occData[,4]<-1
  }
  
  if(is.null(effData)){
    warning("No effort data provided. Assuming continuous sampling.")
    effData <- data.frame(Site=levels(as.factor(occData[,1])),initDate=min(occData[,2]),endDate=max(occData[,2]))
  }
  
  
  occData[,1]<-as.factor(occData[,1])
  effData[,1]<-factor(effData[,1],levels = levels(occData[,1]))
## Create effort matrix
date_range<-range(effData[,2],effData[,3],na.rm=TRUE)

# Create date breaks (including one additional time step to avoid cutting off last occurrences)
breaks<-seq(date_range[1],date_range[2],by = interval)

if(length(breaks)<2){
  stop("Effort interval does not contain more than one ocasion. Try to reduce the interval.")
}

addTime<-diff.Date(breaks[1:2])
breaks2<-c(breaks,breaks[length(breaks)]+addTime)

# Check if first and last date are included
cut(min(effData[,2],effData[,3],na.rm = TRUE),breaks2)
cut(max(effData[,2],effData[,3],na.rm = TRUE),breaks2)

# Cut sampling data into intervals
effData_int<-effData
effData_int[,2]<-cut(effData[,2],breaks2)
effData_int[,3]<-cut(effData[,3],breaks2)
effData_int[,4]<-as.integer(effData_int[,2])
effData_int[,5]<-as.integer(effData_int[,3])

# Create sampling effort matrix (by record on/off)
effort_mat<-matrix(NA,nrow(effData),ncol=nlevels(effData_int[,2]))

# Save levels as category numbers and initial dates
levels<-1:nlevels(effData_int[,2])
levels_date<-as.Date(levels(effData_int[,2]))

# Check which dates cameras were on (per on/off record)
for(i in 1:nrow(effData)){
  effort_mat[i,]<-ifelse(levels>=effData_int[i,4]&levels<=effData_int[i,5],1,NA)
}

# Aggregate records on/off by sampling site
effData_mat2<-aggregate(effort_mat,list(effData[,1]),sum,na.rm=TRUE,drop=FALSE)

## Finish sampling effort matrix:
# 0 --> NA = Camera off
# 1 --> = Camera on
effData_mat3<-ifelse(effData_mat2[,-1]>0,1,NA)

####
dim(effData_mat3)
length(levels_date)


### Plot intervals by records on/off ####
# plot(range(effData[,2],effData[,3],na.rm = TRUE),c(1,nrow(effData)),type="n",xlab="Date",ylab="Camera",yaxt="n")
# 
# abline(v=as.Date(levels(effData_int[,2])))
# 
# segments(effData[,2],1:nrow(effData),effData[,3],1:nrow(effData))
# 
### Plot records by sampling site ####
if(plot){
plot(range(effData[,2]-addTime,effData[,3]+addTime,na.rm = TRUE),c(1,nrow(effData_mat2)),type="n",xlab="Date",ylab="Camera",yaxt="n")

axis(2,at = 1:nrow(effData_mat2),effData_mat2[,1],las=1,cex.axis=0.4)

abline(v=levels_date,lwd=0.2,lty=2)
abline(v=levels_date+addTime/2,lwd=0.2)

for(i in 1:nrow(effData_mat2)){
  segments(levels_date[as.logical(effData_mat3[i,])]-addTime/4,i,levels_date[as.logical(effData_mat3[i,])]+addTime/4,i,lwd=2,col=2)
}

}

#### Prepare occurrence data ####

# array with 3 dimensions:
# 1. Cameras
# 2. Days
# 3. Species

# Cut dates into intervals (days,years, etc.)
occData_cut<-cut(occData[,2],breaks2)


# Create occurrence array
occ<-tapply(occData[,3],
            list(occData[,1],
                 occData_cut,
                 occData[,4]),occFun)



dim(occ)
range(occ)
range(occ,na.rm = TRUE)
sum(occ,na.rm = TRUE)

# Check if effort and occurrence data are compatible
sum(dimnames(occ)[[1]]!=effData_mat2[,1])==0

# Convert NAs to zeroes
occ[is.na(occ)]<-0
range(occ)


effData_arr3<-array(effData_mat3,dim(occ))
# Remove zeroes where true NAs
occ2<-occ*effData_arr3

if(dim(occ2)[3]==1&simplify){
   occ2<-occ2[,,1]
 }

return(occ2)

}




