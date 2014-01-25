##### Get command arguments
get_command_args = function(args = commandArgs(),pos=".GlobalEnv"){
  if (any(args=="--args")){
    try(args[-(1:which(args=="--args"))])->raw.args
    var_names <- sub("=.*","",raw.args)
    vars <<- sub(".*=","",raw.args) # if there is no '=' sign, the  variable name will also be its value
    names(vars) <- var_names
  
    # Default values:
  
    for(i in 1:length(vars)){
      assign(names(vars)[i],vars[i],pos=pos)               
      cat(sprintf("Assigned %s to %s\n",names(vars)[i],vars[i])) 
    }
  }
}


## Write a tab-delimited file:
write.delim = function(dat,File,sep='\t',row.names=FALSE,quote=FALSE,col.names=TRUE, ...){
  write.table(dat,File,sep=sep,row.names=row.names,quote=quote,col.names=col.names, ...)
}

## From a list, get indicies of all things that are duplicates
getAllDups = function(myList){
  which.dup = duplicated(myList)
  which.alldup = myList %in% myList[which.dup]
  which.alldup
}

## Get the minimum of a row:
rowMins = function(dat,dim=1,na.rm=T)
  apply(dat,dim,min,na.rm=na.rm)
  
## Get the maximums of a row:  
rowMaxes = function(dat,dim=1,na.rm=T)
  apply(dat,dim,max,na.rm=na.rm)

## Get the minimum of a column:
colMins = function(dat,dim=2,na.rm=T)
  apply(dat,dim,min,na.rm=na.rm)
  
## Get the maximums of a column:  
colMaxes = function(dat,dim=2,na.rm=T)
  apply(dat,dim,max,na.rm=na.rm)
  

### Add 'arrow's to a plot to indicate standard error
addErrorBars = function(X,Y,SEM,len=1,...) { # Other options are passed to arrows
  arrows(X,Y,X,Y+SEM,angle=90,length=len,...);
  arrows(X,Y,X,Y-SEM,angle=90,length=len,...);
}


## as.numeric(as.character())
as.nc = function(x) as.numeric(as.character(x))

## Take the log (base 10), but return NA if the value started at 0
log_or_NA = function(dat){
  dat2 = log10(dat)
  dat2[is.infinite(dat2)] = NA
  dat2
}

plot2Hists = function(X,main,xlab,xlim=range(X,na.rm=T)+c(-.01*diff(range(X,na.rm=T)),.01*diff(range(X,na.rm=T))),N=100,breaks=seq(xlim[1],xlim[2],length.out=N),freqMax=NULL,cols=c(rgb(0,0,1,alpha=0.4),rgb(1,0,0,alpha=0.4)),leg=NULL,leg.at='topright',...){
  par(new=F)
  mains = c(main,'');
  xaxt=c('n','s');ylab=c('Frequency','');xlab=c(xlab,'');
  cols2 = cols
  print(wilcox.test(X[[1]],X[[2]]))
     
  # Just get the frequencies once first
   a=lapply(X,function(x,main){
    hist(x,freq=F,breaks=breaks,plot=F)
   })
  #freqMax seq(0,0.02,0.14)
  RANGE = range(c(unlist(a[[1]]$intensities),unlist(a[[2]]$intensities)))
    freqMax = seq(0,max(RANGE),length.out=10)     
 
  ylim=c(0,max(freqMax))
 
  a=lapply(X,function(x,main){
    a=hist(x,freq=F,breaks=breaks,xlim=xlim,ylim=ylim,col=cols[1],border=cols[1],main=mains[1],xaxt=xaxt[1],yaxt='n',ylab=ylab[1],xlab=xlab[1],cex.lab=1.5,cex.axis=1.2,cex.main=1.5,...);
    axis(2,at=freqMax,lab=sprintf("%.3f",freqMax*diff(breaks)[1]),cex.axis=1.2)
    #abline(v=median(x,na.rm=T),col=cols[1],lwd=2);
    cols<<-cols[2];
    mains<<-mains[2];
    xaxt<<-xaxt[2];
    ylab<<-ylab[2];
    xlab<<-xlab[2];
    par(new=T)
    a
  }) 
  if(!is.null(leg)) legend(leg.at,fill=cols2,leg=leg)
  par(new=F)
  a
}  

plot2Densities = function(X,main,xlab,xlim=range(X,na.rm=T)+c(-.01*diff(range(X,na.rm=T)),.01*diff(range(X,na.rm=T))),cols=c(rgb(0,0,1),rgb(1,0,0)),leg=NULL,leg.at='topright',...){
  par(new=F)
  mains = c(main,'');
  cols2 = cols
  #print(wilcox.test(X[[1]],X[[2]]))
     
  # Just get the frequencies once first
   a=lapply(X,function(x,main){
    a=density(x[!is.na(x)]);
   })

  RANGE = range(c(unlist(a[[1]]$y),unlist(a[[2]]$y)))   
  ylim=c(0,RANGE[2])
 
  for (i in 1:length(X)){
    if (i == 1){
      plot(density(X[[i]][!is.na(X[[i]])], na.rm=T),col=cols[i],xlim=xlim,ylim=ylim,main=main,xlab=xlab,cex.lab=1.5,cex.axis=1.2,cex.main=1.5, ...);
    }else{
      lines(density(X[[i]][!is.na(X[[i]])], na.rm=T),col=cols[2],xlim=xlim,ylim=ylim, ...)
    }
  }    
  if(!is.null(leg)) legend(leg.at,fill=cols2,leg=leg)
}

heatmap1 = function(dat,main='',col=greenred,length.cols=20, maxVal=NULL,Rowv=FALSE,Colv=FALSE,margins=c(10,5),...){#samples=dimnames(dat)[[2]],geneNames=rownames(dat)){#
  require(gplots)

  if(is.null(maxVal)){
    RANGE = max(abs(range(dat,na.rm=T))) ;#*1.01, breaks=seq(-color_range,color_range))        
    breaks = seq(-RANGE,RANGE,length.out=length.cols)
  }else
    breaks = seq(-maxVal,maxVal,length.out=length.cols)
    h1=heatmap.2(dat,col=col,breaks=breaks,trace='n',Rowv=Rowv,margins=margins,Colv=Colv,main=main,...)
}

# This Function normalizes a vector so that the max is 1 and the min is 0.
UnitSize = function(dat){
  RANGE = range(dat,na.rm=T)
  dat = dat - RANGE[1]
  dat/max(dat)
}






## Here is a way to redirect plots to either go to pdf or to dev.new
options(dev_option = 'pdf',dev_override=F) #'png', 'console' #dev_override Basically is debug mode
# Computer-specific options to make plots start in the top left corner of my 2nd screen ...
options(dev_xpos=1210, dev_ypos=290)

# 
plot.dev = function(arg1='',plotOption=getOption('dev_option'), ... ){
    outStyle = ifelse(getOption('dev_override'), 'console',plotOption)
    if(outStyle=='console'){
      eval(parse(text=sprintf("dev.new(title='%s',xpos=%s,ypos=%d,...)",arg1, getOption('dev_xpos'),getOption('dev_ypos'))))
    }else{
      eval(parse(text=sprintf("%s('%s',...)",outStyle, arg1)))
    }        
}

plot.off = function()if(!getOption('dev_override'))dev.off()

debugMode= function()options(dev_override=TRUE)
fileMode= function()options(dev_override=FALSE)