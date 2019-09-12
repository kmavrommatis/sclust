plot.cluster=function(x,y,cluster,curve,name="")
    {
        par(mar=c(5,5,2,2))
        plot(NULL,NULL,xlim=c(0,max(x)),ylim=c(0,max(cluster[,2])*1.3),xlab="cancer cell fraction",ylab="density",axes=F,cex.lab=1.5,main=name)
        axis(1,cex=1.3,lwd=2)
        axis(2,cex=1.3,lwd=2)
        lines(x,y,type="h",lwd=3)
        lines(curve[,1],curve[,2],lwd=2,col="red")
        nc=dim(cluster)[1]
        for(i in 1:nc)
            {
                lines(c(cluster[i,1],cluster[i,1]),c(0,cluster[i,2]*1.1),lwd=2,col="red",lty=3)
                text(cluster[i,1],cluster[i,2]*1.105,paste("cluster ",i-1,": ",cluster[i,1],sep=""),pos=3)
            }
    }
