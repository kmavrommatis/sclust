plot.profile=function(file,p.est,name="")
    {
        data=read.table(paste(file,"_cn_profile.txt",sep=""),header=T)
        options(warn=-1)
        par(mar=c(0,5,5,0))
        layout(matrix(c(1,2,3,4),4,1),1,c(4,4,6,6))

        x=data[,1]
        plot(NULL,NULL,xlim=c(min(x),max(x)),ylim=c(0,1),axes=F,main=name,cex.lab=1.4,font.lab=2,xlab="",ylab="purity",cex.main=1.6)
        y=data[,2]
        lines(x,y,lwd=3)
        points(x,y,col="grey",pch=21,bg="grey",cex=1.2)
        axis(2,at=c(0,0.25,0.5,0.75,1),label=c(0,0.25,0.5,0.75,1),cex=1.2,lwd=2)

        par(mar=c(0,5,2,0))
        
        plot(NULL,NULL,xlim=c(min(x),max(x)),ylim=c(1,6),axes=F,main="",cex.lab=1.4,font.lab=2,xlab="",ylab="ploidy",cex.main=1.6)
        y=data[,3]
        y[y>=6]=6
        lines(x,y,col="black",lwd=3)
        points(x,y,col="grey",pch=21,bg="grey",cex=1.2)
        axis(2,at=c(1,2,3,4,5,6),label=c("1","2","3","4","5",">6"),cex=1.2,lwd=2)
        lines(c(min(x),max(x)),2*c(1,1),col="grey",lwd=1,lty=2)
        lines(c(min(x),max(x)),3*c(1,1),col="grey",lwd=1,lty=2)
        lines(c(min(x),max(x)),4*c(1,1),col="grey",lwd=1,lty=2)
        lines(c(min(x),max(x)),5*c(1,1),col="grey",lwd=1,lty=2)

        y=data[,6]
        ym=y
        ym[ym > 1e9]=0
        plot(NULL,xlim=c(min(x),max(x)),ylim=c(min(c(y,data[,6]))*0.5,max(ym)*1.5),log="y",axes=F,main="",cex.lab=1.4,font.lab=2,xlab="",ylab="L-cn",cex.main=1.6)
        lines(x,y,col="black",lwd=4)
        axis(2,cex=1.2,lwd=2)
        points(x,y,col="grey",pch=21,bg="grey",cex=1.2)
        lines(c(p.est,p.est),c(min(c(y,data[,5]))*0.5,max(y)*1.5),lwd=1,lty=1,col="red")

        par(mar=c(5,5,2,0))
        
        y=data[,5]
        ym=y
        ym[ym > 1e9]=0
        plot(NULL,xlim=c(min(x),max(x)),ylim=c(min(c(y,data[,5]))*0.5,max(ym)*1.5),log="y",axes=F,main="",cex.lab=1.4,font.lab=2,xlab="expected ploidy",ylab="L-bi",cex.main=1.6)
        lines(x,y,col="black",lwd=4)
        points(x,y,col="grey",pch=21,bg="grey",cex=1.2)
        lines(c(p.est,p.est),c(min(c(y,data[,5]))*0.5,max(y)*1.5),lwd=1,lty=1,col="red")
        
        axis(1,cex=1.2,lwd=2)
        axis(2,cex=1.2,lwd=2)
    }

cmp.dist=function(x,theta,sig,p)
{
  n = length(theta)
  out = numeric(length(x))
  for(i in 1:n)
    {
      out = out + p[i]/(sqrt(2*p[i])*sig[i])*exp(-(x-theta[i])^2/(2*sig[i]*sig[i]))
    }
      
  return(out)
}
  
snps.alleles=function(d,purity,ploidy,name="")
  {
    par(mar=c(5,5,2,0))
    layout(matrix(c(1,2),1,2),c(10,3),1)

    am = d[,8]+10*d[,9]
    uam = unique(am)
    ccol = numeric(length(am))
    leg = character(0)
    exp.theta = numeric(0)
    exp.cn = numeric(0)
    nu = length(uam)
    for(i in 1:nu)
      {
        ind = which(am == uam[i])
        ccol[ind] = i
        leg = c(leg,paste(d[ind[1],8],d[ind[1],9]))
        exp.theta = c(exp.theta,d[ind[1],10])
        exp.cn = c(exp.cn,d[ind[1],1])
      }
    ps=numeric(0)
    for(i in 1:length(d[,11]))
      ps=c(ps,max(d[i,11],0.5))
        
    plot(d[,1],d[,2],col=ccol,axes=F,xlab="copy number",ylab="theta",cex.lab=2,font.lab=2,xlim=c(0,6.5),ylim=c(0,1.1),main=name,sub=paste(" purity=",purity,"  ploidy=",ploidy,sep=""),cex.main=1.5,cex=ps,pch=21,bg=ccol)
    axis(1,cex=1.6,font=2,lwd=2,at=seq(0,5,0.5),label=seq(0,5,0.5))
    axis(2,cex=1.6,font=2,lwd=2)

    points(rep(5.4,nu),(1:nu)/nu,col=1:nu,pch=21,bg=1:nu,cex=1.5)
    text(rep(5.45,nu),(1:nu)/nu,leg,pos=4,font=2,cex=1.5)
    
    points(5.4,0,col=2,pch=3,bg=2,cex=2,lwd=2)
    #points(5.4,0,col="black",pch=20,bg="black",cex=1)
    text(5.45,0,"model",pos=4,font=2,cex=1.3)
    
    points(exp.cn,exp.theta,col=1:nu,pch=3,bg=1:nu,cex=2,lwd=2)
    #points(exp.cn,exp.theta,col="black",pch=20,bg="black",cex=1)
    
    theta.t = d[,2]
    theta.sig.t = d[,3]
    theta.n = d[,4]
    theta.sig.n = d[,5]
    p = d[,7]/sum(d[,7])

    x=seq(0,1,0.001)
    dist.t=cmp.dist(x,theta.t,theta.sig.t,p)
    dist.n=cmp.dist(x,theta.n,theta.sig.n,p)
    
    
    #for(i in 1:dim(cluster)[1])
    #  {
    #    lines(c(0,5),c(cluster[i,3],cluster[i,3]),col=i,lwd=2,lty=2)
    #    text(0,cluster[i,3]+0.03,paste("theta = ",round(cluster[i,3],3),sep=""),col=i,pos=4)
    #    text(2.51,cluster[i,3]+0.03,paste("m = ",round(sum(p[d[,6]==i-1]),2),sep=""),col=i,pos=4)
    #  }
    
    #text(0,purity-0.04,paste("purity = ",round((2*purity/(1+purity))*100,1)," %",sep=""),col=which(cluster[,3]==purity),pos=4)
    par(mar=c(5,0,2,2))
    plot(dist.n,x,axes=F,xlab="dist.",ylab="",type="l",lwd=2,col="blue",cex.lab=2,font.lab=2,ylim=c(0,1.1))
    lines(dist.t,x,lwd=3,col="red")
    axis(1,cex=1.8,font=2,lwd=2)
  }

plot.cn.profile=function(fname,purity,ploidy)
    {
        layout(1,1,1)
        d=read.table(paste(fname,"_allelic_states.txt",sep=""),header=T)
        s=read.table(paste(fname,"_subclonal_cn.txt",sep=""),header=T)
        d=d[which(d$Chromosome!="chrX"),]
        chrs=as.character(unique(d$Chromosome))
        chrs.single=sub("chr","",chrs)
        chr.size=numeric(0)
        for(i in 1:length(chrs))
            {
                chr.size=c(chr.size,max(d[d$Chromosome==chrs[i],4]))
            }
        chr.start=c(0)
        for(i in 2:length(chr.size))
            {
                chr.start=c(chr.start,chr.start[i-1]+chr.size[i-1])
            }

        max.cn=5
        par(mar=c(1,2,1,1))
        max.l=sum(chr.size)
        plot(NULL,xlim=c(0,max.l),ylim=c(-2,max.cn+2),axes=F,xlab="",ylab="")

        text(0,max.cn+1.5,paste(fname,": purity=",purity," ploidy=",ploidy,sep=""),pos=4,cex=1,font=2)
        for(i in 0:max.cn)
            lines(c(0,max.l),c(i,i),col="grey",lwd=1,lty=2)
        for(i in 2:length(chr.start))
            {
                lines(chr.start[i]*c(1,1),c(-0.2,6.2),col="grey",lwd=1,lty=2)
                lines(chr.start[i]*c(1,1),c(-0.3,-0.2),col="grey",lwd=2,lty=1)
            }

        rect(0,-0.2,max.l,6.2,lwd=2)

        axis(2,at=c(0,1,2,3,4,5,6),label=c("0","1","2","3","4","5",">5"),lwd=2,pos=0)
        
        dy=0.2
        for(i in 1:length(chrs.single))
            {
                dy=dy*(-1)
                text(chr.start[i]+0.5*chr.size[i],-0.7+dy,chrs.single[i],cex=1)
            }
        dy=0.1
        for(i in 1:dim(d)[1])
            {
                start=chr.start[chrs==d$Chromosome[i]]

                if(d$Is_Subclonal_CN[i] == 0 || d$Is_Inconsistent_State[i] == 1)
                    {
                        cn=d$CopyNr[i]
                        if(cn > max.cn)
                            cn=6
                        rect(start+d$Start[i],cn-dy,start+d$End[i],cn+dy,border=NA,lwd=0,col="orange4")
                        cn=min(d$A[i],d$B[i])
                        if(cn > max.cn)
                            cn=6
                        rect(start+d$Start[i],cn-dy,start+d$End[i],cn+dy,border=NA,lwd=0,col="darkgreen")
                    }
            }
        if(dim(s)[1] != 0)
            {
                for(i in 1:dim(s)[1])
                    {
                        start=chr.start[chrs==s$Chromosome[i]]
                        cn=s$Subclonal_CN[i]
                        if(cn > max.cn)
                            cn=6
                        rect(start+s$Start[i],cn-dy,start+s$End[i],cn+dy,border=NA,lwd=0,col="orange")

                        
                        cn=s$Clone1_Fraction[i]*min(s$Clone1_A[i],s$Clone1_B[i])+s$Clone2_Fraction[i]*min(s$Clone2_A[i],s$Clone2_B[i])
                        if(cn > max.cn)
                            cn=6
                        rect(start+s$Start[i],cn-dy,start+s$End[i],cn+dy,border=NA,lwd=0,col="green3")
                    }
            }

        y=-1.5
        rect(0,y-dy,max.l*0.04,y+dy,border=NA,lwd=0,col="orange4")
        text(max.l*0.04,y,"copy number",cex=1,pos=4)

        rect(max.l*0.5,y-dy,max.l*0.54,y+dy,border=NA,lwd=0,col="orange")
        text(max.l*0.54,y,"subclonal copy number",cex=1,pos=4)

        y=-1.5-4*dy
        rect(0,y-dy,max.l*0.04,y+dy,border=NA,lwd=0,col="darkgreen")
        text(max.l*0.04,y,"minor allele copy number",cex=1,pos=4)

        rect(max.l*0.5,y-dy,max.l*0.54,y+dy,border=NA,lwd=0,col="green3")
        text(max.l*0.54,y,"minor allele subclonal copy number",cex=1,pos=4)
        
    }

