
# SHOULD ONLY OCCUR IN ONE FILE
#common.dir = "/Users/jeremyliu1/galaxy-dist/tools/my_tools"

# commonDir from region_motif_intersect.r file
dyn.load(paste(commonDir, "/region_motif_lib/regions.so",sep=""))

##reg = matrix(cbind(from,to)) from<to
##region[[chr]] = reg
##pos = unique(integer())
##poslist = list(chr,pos, optional(strand=c(-1,0,+1)))

# USED
merge.reg <- function(...,sep=1) {
  ##This function returns union of regs.
  reg = rbind(...)
  x=.C("merge_regions",as.integer(t(reg)),as.integer(nrow(reg)),as.integer(sep))
  reg=matrix(x[[1]][1:(x[[2]]*2)],ncol=2,byrow=TRUE)
  reg = matrix(reg[which(reg[,2]>reg[,1]),],ncol=2)
  reg[which(reg==0)]=1
  return(reg)
}

merge.regions<-function(...,sep=1) {
  ##This function returns union of regions.
    regions=list(...)
    chrs = unique(unlist(lapply(regions,names),use.names=F))
    region = list()
    for(chr in chrs) {
        region[[chr]] = do.call("merge.reg",c(lapply(regions,function(i) i[[chr]]),sep=sep))
  }
  return(region)
}

plot.reg<-function(reg,xlim=NULL,y=NULL,vertical=FALSE,...) {
  ##This function does not stack if reg is overlapping.
  ##new plot is made unless y is specified.
  if(nrow(reg)==0) return()
  if(is.null(xlim)) xlim=range(reg)
  if(is.null(y)) {
    plot(xlim,c(0,1),type="n",axes=FALSE,xlab=" ",ylab=" ")
    y=0.5
  }
  segments(reg[,1],y,reg[,2],...)
  if(vertical) abline(v=reg)
}

distance.to.closest.reg.of.reg <- function(reg,reg2) {
  ##for each element of reg, what is the closest distance to any element of reg2?
  reg2 = merge.reg(reg2)
  reg2 = c(-Inf,t(reg2),Inf)
  s=reg[,1]
  e=reg[,2]
  sbin = as.integer(cut(s,reg2))
  ebin = as.integer(cut(e,reg2))
  d = pmin(s-reg2[sbin], reg2[sbin+1]-s, e-reg2[ebin], reg2[ebin+1]-e)
  d[which(sbin!=ebin | sbin%%2==0)] = 0
  return(d)
}

# USED
distance.to.closest.reg.of.pos <- function(pos,reg) {
  ##for each element of pos, what is the closest distance to any element of reg?
  reg = merge.reg(reg)
  reg = c(-Inf,t(reg),Inf)
  pbin = as.integer(cut(pos,reg))
  d = pmin(pos-reg[pbin], reg[pbin+1]-pos)
  d[which(pbin%%2==0)] = 0
  return(d)
}

distance.to.closest.pos.of.reg <- function(reg,pos,pos.strand=NULL,index.return=FALSE) {
  ##for each element of reg, what is the closest distance to any element of pos?
  ##if strand is provided, distance is along strand
  o = order(pos)
  pos =  c(-Inf,pos[o],Inf)
  o = c(o[1],o,o[length(o)])
    
  s=reg[,1]
  e=reg[,2]
  sbin = as.integer(cut(s,pos))
  ebin = as.integer(cut(e,pos))
  
  d=integer(nrow(reg))
  s.is.closer = s-pos[sbin] < pos[sbin+1]-e
  if(index.return) {
    return(ifelse(s.is.closer,o[sbin],o[sbin+1]))
  }
  d = ifelse(s.is.closer, s-pos[sbin], e-pos[sbin+1])
  d[which(sbin!=ebin)] = 0
  if(!is.null(pos.strand)) {
    reg.strand = ifelse(s.is.closer,pos.strand[o][sbin],pos.strand[o][sbin+1])
    d = d * reg.strand
  }
  return(d)
}

if(F) {
    pos = sample(seq(0,1000,200))
    pos2 = sample(seq(10,1010,100))
    pos.strand = sample(c(1,-1),6,replace = T)
    pos2.strand = sample(c(1,-1),11,replace = T)
}

distance.to.closest.pos.of.pos <- function(pos,pos2,pos.strand=NULL,pos2.strand=NULL, ignore.pos.strand=TRUE,index.return=FALSE) {
  ##for each element of pos, what is the closest distance to any element of pos2?
  ##if index.return==TRUE, index of pos2 closest to pos is returned
  ##else if strand2 is provided, distance is along strand2
  ##if strand and strand2 are both provided and !ignore.pos.strand
  ##  then output is a list giving plus.up, plus.down, minus.up, minus.down
  ##    plus.up: distance to closest upstream on the same same strand etc. etc. 
  o = order(pos2)
  pos2 =  c(-Inf,pos2[o],Inf)
  if(!is.null(pos2.strand))   pos2.strand = c(-Inf,pos2.strand[o],Inf)

  if(is.null(pos2.strand) | is.null(pos.strand) | ignore.pos.strand) {
    pbin = as.integer(cut(pos,pos2))
    
    pbin = ifelse(pos-pos2[pbin] < pos2[pbin+1]-pos,pbin,pbin+1)
    d = pos-pos2[pbin]
    if(!is.null(pos2.strand)) d = d * pos2.strand[pbin]
    
    if(index.return) return(o[pbin-1])
    return(d)
  }
  strands = list(plus=1,minus=-1)
  relcoords = list(up=0,down=1)
  ind = lapply(strands,function(strand) {
    ind.p = c(1,which(pos2.strand==strand),length(pos2))
    pbin.p = cut(pos,pos2[ind.p],labels=FALSE)
    as.data.frame(lapply(relcoords,function(i) ind.p[pbin.p+i]))
  })
  ind.temp = ind
  ind.minus = which(pos.strand==-1)
  if(length(ind.minus)>0) {
      ind[[1]][ind.minus,]=ind.temp[[2]][ind.minus,2:1]
      ind[[2]][ind.minus,]=ind.temp[[1]][ind.minus,2:1]
  }
  ind = unlist(ind,recursive=FALSE)
  if(index.return) {
    return( lapply(ind,function(i) {
      i[which(i==1)]=NA
      i[which(i==length(pos2))]=NA
      o[i-1]
    }) )
  }
  return(lapply(ind,function(i) pos.strand*(pos2[i]-pos)))
}

distance.to.closest.region.of.region <- function(region,region2) {
  ##for each element of region[[chr]], what is the closest distance to any element of region2[[chr]]?
  ##returns d[[chr]]
  chrs = names(region)
  d=list()
  for(chr in chrs) {
    if(is.null(region2[[chr]])) {
      d[[chr]] = rep(Inf,nrow(region[[chr]]))
    } else {
      d[[chr]] = distance.to.closest.reg.of.reg(region[[chr]],region2[[chr]])
    }
  }
  return(d)
}

# USED
distance.to.closest.region.of.poslist <- function(poslist,region) {
  ##for each element of poslist, what is the closest distance to any element of region?
  chrs = names(table(poslist$chr))
  d=integer()
  for(chr in chrs) {
    ind = which(poslist$chr==chr)
    pos=poslist$pos[ind]
    if(is.null(region[[chr]])) {
      d[ind] = Inf
    } else {
      d[ind] = distance.to.closest.reg.of.pos(pos,region[[chr]])
    }
  }
  return(d)
}
distance.to.closest.poslist.of.region <- function(region,poslist,index.return=FALSE) {
  ##for each element of region, what is the closest distance to any element of poslist?
  chrs = names(region)
  d=list()
  for(chr in chrs) {
    ind = which(poslist$chr==chr)
    pos=poslist$pos[ind]
    pos.strand=poslist$strand[ind]
    d[[chr]] = distance.to.closest.pos.of.reg(region[[chr]],pos,pos.strand,index.return=index.return)
    if(index.return) d[[chr]] = ind[d[[chr]]]
  }
  return(d)
}

distance.to.closest.poslist.of.poslist <- function(poslist,poslist2,ignore.poslist.strand=TRUE,index.return=FALSE) {
  ##for each element of poslist, what is the closest distance to any element of poslist2?
  ##if poslist2$strand is provided, distance is along strand2
  ##if strand and strand2 are provided and no ignore.poslist.strand
  ##  then output is a list giving plus.up, plus.down, minus.up, minus.down
  ##    plus.up: distance to closest upstream on the same same strand etc. etc. 
  ##if index.return==TRUE, index of pos2 closest to pos is returned
  
  chrs = names(table(poslist$chr))
  
  d=integer()
  stranded = !(is.null(poslist2$strand) | is.null(poslist$strand) | ignore.poslist.strand)
  if(stranded) {
    brs = c("plus.up","plus.down","minus.up","minus.down")
    d=list()
    for(br in brs) d[[br]]=integer()
  }
  
  for(chr in chrs) {
    ind = which(poslist$chr==chr)
    ind2 = which(poslist2$chr==chr)
    pos=poslist$pos[ind]
    pos2=poslist2$pos[ind2]
    pos.strand=poslist$strand[ind]
    pos2.strand=poslist2$strand[ind2]
    if(!stranded) {
      d[ind] = distance.to.closest.pos.of.pos(pos,pos2,pos.strand,pos2.strand,ignore.poslist.strand,index.return=index.return)
      if(index.return) d[ind] = ind2[d[ind]]
    } else {
      x =  distance.to.closest.pos.of.pos(pos,pos2,pos.strand,pos2.strand,ignore.poslist.strand)
      for(br in brs) {
        d[[br]][ind] = x[[br]]
        if(index.return) d[[br]][ind] = ind2[d[[br]][ind]]
      }
    }
  }
  return(d)
}


reg.minus.reg <- function(reg,reg2) {
  x = .C("region_minus_region",as.integer(t(reg)),as.integer(nrow(reg)),as.integer(t(reg2)),as.integer(nrow(reg2)),integer((nrow(reg)+nrow(reg2))*2))[[5]]
  x=x[which(x>=0)]
  return(matrix(x,ncol=2,byrow=TRUE))
}

intersection.of.regs <- function(reg,reg2) {
  x = .C("intersection_of_regions",as.integer(t(reg)),as.integer(nrow(reg)),as.integer(t(reg2)),as.integer(nrow(reg2)),integer((nrow(reg)+nrow(reg2))*2))[[5]]
  x=x[which(x>=0)]
  return(matrix(x,ncol=2,byrow=TRUE))
}

region.minus.region<-function(region,region2) {
  chrs = names(region)
  for(chr in chrs) {
    if(is.null(region[[chr]])) next
    if(!is.null(region2[[chr]])) {
      region[[chr]] = reg.minus.reg(region[[chr]],region2[[chr]])
    }
  }
  return(region)
}

intersection.of.regions<-function(region,region2) {
  chrs = names(region)
  for(chr in chrs) {
    if(is.null(region2[[chr]])) {
      region[[chr]]<-NULL
    } else {
      region[[chr]] = intersection.of.regs(region[[chr]],region2[[chr]])
    }
  }
  return(region)
}

reg.around.pos <-function(pos,range=500,strand=NULL) {
  if(length(range)==1) range=c(range,range)
  if(is.null(strand)) strand = 1;
  reg = cbind(pos-range[1]*strand,pos+range[2]*strand);
  ind = which(reg[,2]<reg[,1])
  reg[ind,] =  reg[ind,2:1]
  ind = which(reg<=0)
  reg[ind] = 1
  return(reg)
}


region.around.poslist <-function(poslist,range=500) {
  chrs = names(table(poslist$chr))
  region=list()
  for(chr in chrs) {
    ind = which(poslist$chr==chr)
    pos=poslist$pos[ind]
    strand = 1
    if(!is.null(poslist$strand)) {
      strand = poslist$strand[ind]
    }
    region[[chr]] =  reg.around.pos(pos,range,strand)
  }
  return(region)
}


poslist.of.region.centers <-function(region) {
  chrs = names(region)
  n=sapply(region,nrow)
  return(data.frame(chr=rep(chrs,n),pos=unlist(lapply(region,function(chr)(chr[,1]+chr[,2])/2),use.names = FALSE)))
}

write.gff.region<-function(region,outfname) {
  region = lapply(region,function(chr) list(s=chr[,1],e=chr[,2]))
  out=unlist.chr(region)
  out$chr=rep(names(region),sapply(region,function(i) length(i$s)))
  empty=rep(".",length(out$chr))
  write.table(data.frame(out$chr,empty,empty,out$s,out$e,empty,empty,empty,empty),quote=FALSE,sep="\t",file=outfname,col.names=FALSE,row.names=FALSE)
}

number.of.regions = function(region)sum(sapply(region,nrow))
size.of.regions = function(region) sum(sapply(merge.regions(region),function(reg) sum(reg[,2]-reg[,1])))
