
################## fun ############

cat("Function arrangeOccu installed\n")
cat("Use: arrangeOccu(occData2,effData,interval = 'week')\n")
cat("* occData2 must have 4 columns:\n")
cat("1. Site ID (camera, plot, etc.)\n")
cat("2. Date in date format\n")
cat("3. Abundance (if not availabe, repeat 1)\n")
cat("4. Species name (even if only one species)\n")
cat("* effData must have 3 columns:\n")
cat("1. Site ID (camera, plot, etc.)\n")
cat("2. Initial date in date format\n")
cat("3. Final date in date format\n")

arrangeOccu<-function(occData,effData=NULL,interval="days",plot=FALSE,occFun=length,simplify=TRUE,error='hide'){
  
  effData<-na.omit(effData)
  
  occData<-na.omit(occData)
  occData2<-occData

  error<-c("show","hide")[pmatch(error,c("show","hide"))]
  
  if(ncol(occData2)<3){
    occData2[,3]<-1
  }
  
  if(ncol(occData2)<4){
    occData2[,4]<-1
  }
  
  if(is.null(effData)){
    warning("No effort data provided. Assuming continuous sampling.")
    effData <- data.frame(Site=levels(as.factor(occData2[,1])),initDate=min(occData2[,2]),endDate=max(occData2[,2]))
  }
  
  
  occData2[,1]<-as.factor(occData2[,1])
  effData[,1]<-factor(effData[,1],levels = levels(occData2[,1]))
  
  occData2[,2]<-as.POSIXct(occData2[,2])
  effData[,2]<-as.POSIXct(effData[,2])
  effData[,3]<-as.POSIXct(effData[,3])
  
  sites_out<-!levels(occData2[,1])%in%effData[,1]
  if(sum(sites_out)>0){
    stop("Site in occData2 missing from effort data")
  }
  
  
  ## Create effort matrix
  date_range<-range(effData[,2],effData[,3],na.rm=TRUE)

  
  # Create date breaks (including one additional time step to avoid cutting off last occurrences)
  breaks<-seq(date_range[1],date_range[2],by = interval)
  
  if(length(breaks)<2){
    stop("Effort interval does not contain more than one ocasion. Try to reduce the interval.")
  }
  
  addTime<-diff.Date(breaks[1:2])
  #breaks2<-c(breaks,breaks[length(breaks)]+addTime)
  breaks2<-breaks
  breaks2[length(breaks2)+1]<-"3100-01-01"
  
  # Check if first and last date are included
#  cut(min(effData[,2],effData[,3],na.rm = TRUE),breaks2)
#  cut(max(effData[,2],effData[,3],na.rm = TRUE),breaks2)
  
  # Cut sampling data into intervals
  effData_int<-effData
  effData_int[,2]<-cut(effData[,2],breaks2)
  effData_int[,3]<-cut(effData[,3],breaks2)
  effData_int[,4]<-as.integer(effData_int[,2])
  effData_int[,5]<-as.integer(effData_int[,3])
  
  ver<-list()
  for(i in 1:nrow(effData)){
    st<-seq(as.POSIXct(effData_int[i,2]),as.POSIXct(effData_int[i,3]),by = interval)
    ver[[i]]<-data.frame(site=effData[i,1],date=st)
  }

  ver2<-do.call(rbind,ver)
  effData_cut<-cut(ver2[,2],breaks2)

  #create effort matrix
  eff<-table(ver2[,1],effData_cut)

  #### Prepare occurrence data ####
  
  # Cut dates into intervals (days,years, etc.)
  occData2_cut<-cut(occData2[,2],breaks2)
  
    if(identical(occFun,length)){
  occ<-table(occData2[,1],occData2_cut,occData2[,4])
  }else{
  # Create occurrence array
  occ<-tapply(occData2[,3],
              list(occData2[,1],
                   occData2_cut,
                   occData2[,4]),occFun)

  # Convert NAs to zeroes
  occ[is.na(occ)]<-0
  
    }

  eff_arr3<-array(eff,dim(occ))
  occ<-as.array(occ)

  effort_mismatch<-occ>0&eff_arr3==0

  errors<-{}
  
  if(sum(effort_mismatch)>0){
    cat(sum(effort_mismatch),"(final) occurrences do not overlap with temporal region of effort and will be replaced by NA\n")
    
      resh_avail<-require("reshape2")
      
      if(resh_avail){
        effort_mismatch2<-melt(effort_mismatch,varnames = NULL)
        occData22<-data.frame(Var1=occData2[,1],Var2=occData2_cut,Var3=occData2[,4],occData)
   
        errors<-merge(occData22,effort_mismatch2[effort_mismatch,1:3],by.x = c("Var1","Var2","Var3"),by.y = c("Var1","Var2","Var3"))[,colnames(occData)]
        
      }else{
        cat("Install pacakage 'reshape2' to save output missing data\n")
      }
    
  }
  # Remove zeroes where true NAs
  occ2<-occ
  occ2[eff_arr3==0]<-NA
  occ2[occ2>0]<-1
  
  if(dim(occ2)[3]==1&simplify){
    occ2<-occ2[,,1]
    eff_arr3<-eff_arr3[,,1] 
  }
  
  resu<-list(occ01=occ2,counts=occ,effort=eff_arr3)
  resu$errors<-errors
  
  class(resu)<-c("arrOcc","list")
  
  return(resu)
}

plot.arrOcc<-function(occ,xlab="Date",ylab="Sites",las=1,cex.axis=0.4,col=heat.colors(100)){
  
  #eff<-occ$effort
  counts<-occ$occ01
  countsT<-apply(counts,c(1,2),sum)
  
  dates<-as.POSIXct(colnames(counts))
  
  addTime<-difftime(dates[2],dates[1])
  
  plot(range(dates-addTime,dates+addTime,na.rm = TRUE),c(1,nrow(countsT)),type="n",xlab=xlab,ylab=ylab,yaxt="n",las=las)
    
    axis(2,at = 1:nrow(countsT),rownames(counts),las=las,cex.axis=cex.axis)
    
    abline(v=dates,lwd=0.2,lty=3,col="grey10")
    abline(v=dates+addTime/2,lwd=0.2,col="grey20")
    
    
    for(i in 1:nrow(countsT)){
      colBreaks<-seq(min(countsT,na.rm = TRUE)-1,max(countsT,na.rm = TRUE),length.out = length(col)+1)
      vals<-countsT[i,!is.na(countsT[i,])]
      colors<-col[as.numeric(cut(vals,c(-Inf,colBreaks)))-1]
      segments(dates[!is.na(countsT[i,])]-addTime/4,i,dates[!is.na(countsT[i,])]+addTime/4,i,lwd=2,col=colors)
      }
    
}


