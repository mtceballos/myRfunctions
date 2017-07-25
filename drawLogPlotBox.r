drawLogPlotBox <-function(xlimits=c(1E-10,1E-11), ylimits=c(1E2,1E3),
                          x2limits=NULL, y2limits=NULL, logxy="xy", 
                          xlabel="x axis", ylabel="y axis",
                          x2label="", y2label="",
                          naxes=c(T,T,F,F)){
      # draw a plot box with logarithmic axes
      # xlimits: vector with X axis limits (axis=1)
      # ylimits: vector with Y axis limits (axis=2)
      # x2limits: vector with X2 axis limits (axis=3)
      # y2limits: vector with Y2 axis limits (axis=4)
      # logxy: string specifying which axes should be plotted in log scale:
      #        logxy="x" -> X axis
      #        logxy="y" -> Y axis
      #        logxy="xy" -> X and Y axes
      # xlabel: label for X axis
      # ylabel: label for Y axis
      # x2label: label for X2 axis
      # y2label: label for Y2 axis
      # naxes: logical vector with "T" for the axis number to be plotted 
      #       (i.e.:  c(TRUE,TRUE,TRUE,TRUE) will draw the four axes, while
      #               c(TRUE,TRUE,FALSE,FALSE) will draw only X (bottom) and
      #               Y (left) axes). By default, all four axes will be plotted.
    
    getPowerTenLab <- function(expo){
        if(expo == 0){
            powLab <- "1"
        }else if(expo == 1){
            powLab <- "10"
        }else{
            powLab <- as.expression(bquote(10^.(expo)))
        }
        return(powLab)
    }    
    
  # draw plot box with no axes
  # ==========================
  plot(xlimits,ylimits, type="n",log=logxy, xlab=xlabel, cex=0.8, axes=FALSE,
         ylab=ylabel)  
  box()
  
  # check which axes require log scale and plot them
  
  # X axes
  #========
  if(grepl("x",logxy,ignore.case=TRUE)){ #if log is for X axes
    # set major ticks, major ticks labels and minor ticks for X axis
    expseq<-seq(from=round(log10(xlimits[1])),to=round(log10(xlimits[2])))
    xmjt<-10**expseq
    xmjtlabs<-c()
    xmnt<-c()
    for (expo in expseq){
        xmjtlabs<-append(xmjtlabs,getPowerTenLab(expo))
        xmnt<-append(xmnt,seq(2,9)*10**expo)
    }
    axis(1,at=xmjt,labels=xmjtlabs,las=1,col.axis="black")
    axis(1,at=xmnt,tcl=-0.2,labels=FALSE,col.axis="white")
    if (naxes[3]){
        if(!is.null(x2limits)){
            expseq<-seq(from=round(log10(x2limits[1])),to=round(log10(x2limits[2])))
            xmjt<-(10**expseq)*xlimits[1]/x2limits[1]
            xmjtlabs<-c()
            xmnt<-c()
            for (expo in expseq){
                xmjtlabs<-append(xmjtlabs,getPowerTenLab(expo))
                xmnt<-append(xmnt,seq(2,9)*10**expo)
            }  
            xmnt <- xmnt *xlimits[1]/x2limits[1]
        }
        #axis(3,at=xmjt,col.ticks="grey",cex.axis=0.9,labels=xmjtlabs,padj=0.7)
        axis(3,at=xmjt,col.ticks="grey",cex.axis=0.9,labels=FALSE,padj=0.7)
        axis(3,at=xmnt,col.ticks="grey",tcl=-0.2,labels=FALSE,col.axis="white",padj=0.7)
        title(main=x2label,cex.main=0.8)
    }
  }else{ # if X axis not to be log scaled
    axis(1)
    if (naxes[3]){axis(3,col.ticks="grey",cex.axis=0.9,padj=0.7)}
  }
  
  # Y axes
  #========
  if(grepl("y",logxy,ignore.case=TRUE)){ #if log is for Y axes
    # set major ticks, major ticks labels and minor ticks for Y axis  
    expseq<-seq(round(log10(ylimits[1])),round(log10(ylimits[2])))
    ymjt<-10**expseq
    ymjtlabs<-c()
    ymnt<-c()
    for (expo in expseq){
      ymjtlabs<-append(ymjtlabs,getPowerTenLab(expo))
      ymnt<-append(ymnt,seq(2,9)*10**expo)
    }
    
    axis(2,at=ymjt,labels=ymjtlabs,las=1)
    axis(2,at=ymnt,tcl=-0.2,labels=FALSE,col.axis="white")
    if (naxes[4]){
        if(!is.null(y2limits)){
            expseq<-seq(round(log10(ylimits[1])),round(log10(ylimits[2])))
            ymjt<-10**expseq*ylimits[1]/y2limits[1]
            ymjtlabs<-c()
            ymnt<-c()
            for (expo in expseq){
                ymjtlabs<-append(ymjtlabs,getPowerTenLab(expo))
                ymnt<-append(ymnt,seq(2,9)*10**expo)
            }
            ymnt <- ymnt *ylimits[1]/y2limits[1]
        }
        #axis(4,at=ymjt,col.ticks="grey",cex.axis=0.9,labels=ymjtlabs,las=1,hadj=0.5)
        axis(4,at=ymjt,col.ticks="grey",cex.axis=0.9,labels=FALSE,las=1,hadj=0.5)
        axis(4,at=ymnt,col.ticks="grey",tcl=-0.2,labels=FALSE,col.axis="white",hadj=0.5)
    }
  }else{# if Y axis not to be log scaled
    axis(2,las=1)
    if (naxes[4]){axis(4,col.ticks="grey",las=1,cex.axis=0.9,hadj=0.5)}
  }
}
