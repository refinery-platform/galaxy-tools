library(graphics, quietly=TRUE)

plot.verbose=F
name.cleaner<-function(...,sep="",replace="_") {
  plot.name=gsub(" ",replace,paste(...,sep=sep))
  plot.name=gsub("/",replace,plot.name)
  plot.name=gsub(",",replace,plot.name)
  plot.name=gsub("'",replace,plot.name)
  plot.name=gsub("\\+","plus",plot.name)
  plot.name=gsub("\\(","",plot.name)
  plot.name=gsub("\\)","",plot.name)
  return(plot.name)
}
plot.namer <- function(..., date=0, fig.dir=0, format="png",sep="") {
  plot.name=name.cleaner(...,sep=sep)
  if(date==0) date=gsub("-","",as.character(Sys.Date()))
  if(fig.dir==0) fig.dir="/Users/alver/allplots"
  plot.name=paste(fig.dir,"/",date,plot.name,".",format,sep="")
  if(plot.verbose) cat("  saving figure: ",plot.name,"\n")
  return(plot.name)
}

plot.scatter <- function(x,y=NULL,f=0.1,same=FALSE,n.points=-1,draw.lowess=FALSE,write.r=TRUE,cex.r=1,col=NULL,col.line=NULL,lwd.line=1,
                         draw.loess=FALSE,span=0.5,bandwidth=bandwidth,draw.prof=FALSE,xlog=FALSE,ylog=FALSE,cor.method="pearson",log="",ind=NULL,
                         draw.spread=FALSE,...) {

  ## if col is the same length as x, use col for each point matching x.
  ## if col is the same length as ind, use col for each point matching x[ind].
  ## else use densCols function based on col.
  ## if col is null, densCols is used with bluetone for first plot and redtone for same=T.
    
    #print(length(x))
    #print(length(y))

    xy <- xy.coords(x, y)
    x=xy$x
    y=xy$y
    
    output=list()
    col.use = col

    if(!is.null(ind)) {
        if(length(col.use)==length(x)) {
            col.use=col.use[ind]
        }
        x=x[ind]
        y=y[ind]
    }
    
    if(length(col.use)!=length(x)) {
        col.use=rep(NA,length(x))
    }
  
  
  take=which(is.finite(x) & is.finite(y))
  x=x[take]
  y=y[take]
  col.use=col.use[take]
  
  if(grepl("x",log)) xlog=TRUE
  if(grepl("y",log)) ylog=TRUE
  if(xlog) log="x"
  if(ylog) log=paste(log,"y",sep="")
  
  if(draw.lowess) {
    lo = lowess(x,y,f)
    output$lowess=lo
  }
  if(draw.loess | draw.spread) {
      px=x;py=y
      if(xlog)  px=log(x)
      if(ylog)  py=log(y)
      ind = which(is.finite(px+py))
      px=px[ind]
      py=py[ind]
      lo = loess(py ~ px,span=span,iterations=5)
      lo.y=as.numeric(lo$fitted)
      lo.x=as.numeric(lo$x)
      if(draw.spread) lo.sd = loess((lo.y-py)^2 ~ lo.x,span=span*1.5,iterations=5)
      if(xlog) lo.x=exp(lo.x)
      if(ylog) lo.y=exp(lo.y)
      lo =data.frame(x=lo.x,y=lo.y)
      if(draw.spread) {
          lo.sd=lo.sd$fitted
          if(ylog) lo.sd=lo.sd*lo.y*lo.y
          lo$sd=sqrt(pmax(0,lo.sd))
      }
      lo=unique(lo)
      lo = lo[order(lo$x),]
      output$loess=lo
  }
  
  if(draw.prof) {
    px=x;py=y
    if(xlog)  px=log(x)
    p=prof(px,py,50)
    if(xlog)  p$x=exp(p$x)
    output$prof=p
  }     

  r=cor(x,y,method=cor.method)
  output$cor=r
  output$cor.method=cor.method
  
  len=length(x)
  if(n.points>0 & n.points<len) {
    take=sample(1:len,n.points)
    x=x[take]
    y=y[take]
    col.use=col.use[take]
  }

  if(xlog) {
    ind = which(x>0)
    x=x[ind]
    y=y[ind]
    col.use=col.use[ind]    
  }
  xcol=x
  if(xlog) xcol=log(xcol)
  if(ylog) {
    ind = which(y>0)
    x=x[ind]
    xcol=xcol[ind]
    y=y[ind]
    col.use=col.use[ind]
  }
  ycol=y
  if(ylog) ycol=log(ycol)

  if(is.null(col)) {
    if(!same) {
      col=colorRampPalette(blues9[-(1:3)])
    } else {
      col=colorRampPalette(c("lightpink","red","darkred"))
    }
  }
  if(!is.na(col.use[1])) {
    col=col.use
  } else {
    col= suppressPackageStartupMessages(densCols(xcol,ycol,col =col,bandwidth=bandwidth,nbin=500))
  }
  if(!same) {
    plot(x,y,col=col,log=log,...)
  } else {
    points(x,y,col=col,...)
  }

  if(is.null(col.line)) {
    col.line="darkblue"
    if(same) col.line="darkred"
  }
  if(draw.lowess | draw.loess) lines(lo,col=col.line,lwd=lwd.line)
  if(draw.spread) {
      lines(lo$x,lo$y+lo$sd,col=col.line,lwd=lwd.line)
      lines(lo$x,lo$y-lo$sd,col=col.line,lwd=lwd.line)
  }
  if(draw.prof) {
    points(p)
    plot.prof(p)
  }
  if(write.r & !same) mtext(paste("r=",round(r,3),sep=""),cex=cex.r)
  return(invisible(output))
}

#color.int=c(144,586,465,257,490,100,74,24)
#coli=1
#cols = integer()
colramp.bwr = vector()
colramp.byr = vector()
colramp.bw = vector()
colramp.bw2 = vector()

plot.save=F

setup.plotting <- function() {
  pdf.options(useDingbats = FALSE)
#  cols<<-colors()[color.int]
#  cols<<-rep(cols,100)
  colramp.bwr <<- colorRampPalette(c("blue","white","red"),space="Lab")(100);
  colramp.byr <<- colorRampPalette(c("blue","yellow","red"),space="Lab")(100);
  colramp.bw  <<- colorRampPalette(c("white","black"),space="Lab")(100)
  colramp.bw2  <<- colorRampPalette(c("grey92","grey5"),space="Lab")(100)
}


plot.cluster <- function(x,k, max.points.cl=-1, image.sep=-1,col=NULL, reorder=FALSE) {
    x[which(is.na(x))]=0
    if(reorder) {
        o=hclust(dist(t(x)))$order
        x=x[,o]
    }
    if(image.sep<0) {
        if(max.points.cl>0) {
            image.sep=ceiling(0.2*max.points.cl)
        }    else {
            image.sep=ceiling(0.2 * nrow(x) / nrow(k$centers))
        }
    }
    
    distances<-dist(k$centers)
    hcl=hclust(distances)

    adjust.branch.sep <-function(ddr,lengths) {
        assign.branch.sep <- function(d,i.leaf) {
            if(is.leaf(d)) {
                attr(d,"members")<-lengths[i.leaf]
                i.leaf=i.leaf+1
                output=list(d=d,i.leaf=i.leaf)
                return(output)
            }
            else{
                input=assign.branch.sep(d[[1]],i.leaf)
                i.leaf=input$i.leaf
                d[[1]]=input$d
                
                input=assign.branch.sep(d[[2]],i.leaf)
                i.leaf=input$i.leaf
                d[[2]]=input$d
                
                attr(d,"members")<-attr(d[[1]],"members")+attr(d[[2]],"members")
                output=list(d=d,i.leaf=i.leaf)
                return(output)      
            }
        }
        ddr<-as.dendrogram(ddr)
        ddr=assign.branch.sep(ddr,1)$d
        return(ddr)
    }
    
    n.points.actual=k$size
    if(max.points.cl>0) {
    k$size[which(k$size>max.points.cl)] = max.points.cl
}

    ddr<-adjust.branch.sep(hcl,k$size[hcl$order]+image.sep)
    centers=length(hcl$order)
    
    n.points=sum(k$size)
    n.dims=ncol(x)
    z=matrix(numeric((n.points+(centers-1)*image.sep)*n.dims),ncol=n.dims)


  last.row=0
  cluster.y.pos=numeric(centers)
  for(i.c in hcl$order) {
    n.p=k$size[i.c]
    z[last.row+1:n.p,] = x[which(k$cluster==i.c)[1:n.p],]
    cluster.y.pos[i.c]=last.row+n.p/2
    last.row=last.row+n.p+image.sep
  }
  
  zlim=c(0,max(z))
  if(min(z)<0) {
    m=max(c(z,-z))
    zlim=c(-m,m)
  }
  if(is.null(col)) {
    if(min(z)>=0) {
      col= colramp.bw
    } else {
      col= colorRampPalette(c("blue","yellow","red"),space="Lab")(100);
    }
  }
  x.pl=seq1(n.dims+1)-0.5
  y.pl=seq1(nrow(z)+1)-0.5
  l <- layout(matrix(1:2,ncol=2),widths=c(1,5))
  par(mar = c(6,0.5,6,0))
  my.plot.dendrogram(ddr,horiz=T,axes=F,yaxs="i",xaxs="i",leaflab="none",center=T,lwd=10)
  par(mar = c(6,0.1,6,2.1))
  image(x=x.pl,y=y.pl,z=t(z),zlim=zlim,axes=FALSE,xlab="",col=col)   
  mtext("cluster",side=4,adj=1.1)
  mtext("points",side=4,adj=1.1,line=1)
  mtext(seq1(centers),side=4,at=cluster.y.pos)
  mtext(n.points.actual,side=4,at=cluster.y.pos,line=1)

  if(!is.null(dimnames(x)[[2]])) {
    mtext(dimnames(x)[[2]],side=1,at=seq1(n.dims),las=2)
  }
}

plot.cluster2 <- function(k, n.clusters=-1, n.clusters.per.panel=4, cols=c("black","red","blue","darkgreen","orange"),f=0,xshift=0,...) {
  if(n.clusters<=0) n.clusters=nrow(k$centers)

  n.elements=as.numeric(unlist(lapply(seq1(n.clusters), function(cl) length(which(abs(k$cluster)==cl)))))
  
  distances<-dist(k$centers)
  n.panels = ceiling(n.clusters/n.clusters.per.panel)
  n.rows=ceiling(sqrt(n.panels))
  n.cols=ceiling(n.panels/n.rows)
  n.panels.layout=n.rows*n.cols

  layout(matrix(seq1(n.panels.layout),nrow=n.rows,byrow=TRUE))
  
  min=min(k$centers)
  max=max(k$centers)

  if(f>0) {
    for(i.cluster in seq1(n.clusters)) {
      k$centers[i.cluster,]=lowess(k$centers[i.cluster,],f=f)$y
    }
  }
  
  ##  hcl=hclust(distances)
  hcl=list()
  hcl$order=1:n.clusters
  
  for(i.cluster in seq1(n.clusters)) {
    if(i.cluster %% n.clusters.per.panel == 1 ) {
      clusters.of.panel=i.cluster:(i.cluster+n.clusters.per.panel-1)
      clusters.of.panel=clusters.of.panel[which(clusters.of.panel<=n.clusters)]
      clusters.of.panel=hcl$order[clusters.of.panel]
      plot(c(0,length(k$centers[1,]))+xshift,c(min,max),type="n",...)
      mtext(paste(clusters.of.panel," (",n.elements[clusters.of.panel],")",sep=""),line=length(clusters.of.panel)-seq1(length(clusters.of.panel)),col=cols[seq1(length(clusters.of.panel)) %% n.clusters.per.panel+1] )
    }
   # lines(k$centers[hcl$order[i.cluster],],col=cols[i.cluster %% n.clusters.per.panel+1])
     lines(seq1(length(k$centers[1,]))+xshift,k$centers[hcl$order[i.cluster],],col=cols[i.cluster %% n.clusters.per.panel+1])
  }
}

my.colors <- function(n) {
  few.colors=c("black","red","blue","green3","mediumorchid3","gold2","darkcyan","sienna2")
  if(n<=length(few.colors)) return(few.colors [seq1(n)])
  col=integer(n)
  n.families=7
  n.members=ceiling(n/n.families)
  for(i in seq1(n)) {
    member=ceiling(i/n.families)
    ratio=(member-1)/(n.members-1)
    c2=0+0.8*ratio
    if(member %% 2 == 1) ratio=-ratio
    c1=0.8-0.2*ratio
    c3=0.75-0.2*ratio
    if(i %% n.families == 1) {col[i]=rgb(c2,c2,c2)}
    if(i %% n.families == 2) {col[i]=rgb(c1,c1/2,c1/2)}
    if(i %% n.families == 3) {col[i]=rgb(c1/2,0.9*c1,c1/2)}
    if(i %% n.families == 4) {col[i]=rgb(c1/2,c1/2,c1)}
    if(i %% n.families == 5) {col[i]=rgb(c3,c3,c3/2)}
    if(i %% n.families == 6) {col[i]=rgb(c3,c3/2,c3)}
    if(i %% n.families == 0) {col[i]=rgb(c3/2,c3,c3)}
  }
  return(col)
}

plot.my.colors <-function(n) {
  x11()
  col=my.colors(n)
  plot(x=c(0,n),y=c(0,1),type="n")
  segments(seq1(n)-1,runif(n),seq1(n),runif(n),col=col)
}


plot.colors <-function() {
  x11(width=10,height=10)
  plot(c(0,26),c(0,26),type="n")
  c=colors()
  n=length(c)
  i=1:n
  x=i%%26
  y=floor(i/26)
  rect(x,y,x+1,y+1,col=c[i],border=c[i])
  text(x+0.5,y+0.5,i)
}


adjust.branch.sep <-function(ddr,lengths) {
  assign.branch.sep <- function(d,i.leaf) {
    if(is.leaf(d)) {
      attr(d,"members")<-lengths[i.leaf]
      i.leaf=i.leaf+1
      output=list(d=d,i.leaf=i.leaf)
      return(output)
    }
    else{
      input=assign.branch.sep(d[[1]],i.leaf)
      i.leaf=input$i.leaf
      d[[1]]=input$d
      
      input=assign.branch.sep(d[[2]],i.leaf)
      i.leaf=input$i.leaf
      d[[2]]=input$d
      
      attr(d,"members")<-attr(d[[1]],"members")+attr(d[[2]],"members")
      output=list(d=d,i.leaf=i.leaf)
      return(output)      
    }
  }
  ddr<-as.dendrogram(ddr)
  ddr=assign.branch.sep(ddr,1)$d
  return(ddr)
}
t.dhcol <- function(dr,h,cols=c(1)) {
                                        # check child heights
  if(attr(dr[[1]],"height")<h) {
                                        # color
    ecol <- cols[coli];
    coli <<- coli+1;
    dr[[1]] <- dendrapply(dr[[1]],function(e) { attr(e,"edgePar") <- list(col=ecol); e});
    attr(dr[[1]],"edgePar") <- list(col=ecol,p.border=NA,p.col=NA,t.col=1,t.cex=1.3);
  } else {
    dr[[1]] <- t.dhcol(dr[[1]],h,cols);
  }
  
  if(attr(dr[[2]],"height")<h) {
                                        # color
    ecol <- cols[coli];
    coli <<- coli+1;
    dr[[2]] <- dendrapply(dr[[2]],function(e) { attr(e,"edgePar") <- list(col=ecol); e});
    attr(dr[[2]],"edgePar") <- list(col=ecol,p.border=NA,p.col=NA,t.col=1,t.cex=1.3);
  } else {
    dr[[2]] <- t.dhcol(dr[[2]],h,cols);
  }
  return(dr);
}



### The rest is PeterK's my.plot.dendogram

## FIXME: need larger par("mar")[1] or [4] for longish labels !
## {probably don't change, just print a warning ..}
my.plot.dendrogram <-
    function (x, type = c("rectangle", "triangle"), center = FALSE,
          edge.root = is.leaf(x) || !is.null(attr(x, "edgetext")),
          nodePar = NULL, edgePar = list(),
          leaflab = c("perpendicular", "textlike", "none"), dLeaf = NULL,
          xlab = "", ylab = "", xaxt="n", yaxt="s",
          horiz = FALSE, frame.plot = FALSE, ...)
{
    type <- match.arg(type)
    leaflab <- match.arg(leaflab)
    hgt <- attr(x, "height")
    if (edge.root && is.logical(edge.root))
    edge.root <- 0.0625 * if(is.leaf(x)) 1 else hgt
    mem.x <- .my.memberDend(x)
    yTop <- hgt + edge.root
    if(center) { x1 <- 0.5 ; x2 <- mem.x + 0.5 }
    else       { x1 <- 1   ; x2 <- mem.x }
    xlim <- c(x1 - 1/2, x2 + 1/2)
    ylim <- c(0, yTop)
    if (horiz) {## swap and reverse direction on `x':
    xl <- xlim; xlim <- rev(ylim); ylim <- xl
    tmp <- xaxt; xaxt <- yaxt; yaxt <- tmp
    }
    plot(0, xlim = xlim, ylim = ylim, type = "n", xlab = xlab, ylab = ylab,
     xaxt = xaxt, yaxt = yaxt, frame.plot = frame.plot, ...)
    if(is.null(dLeaf))
        dLeaf <- .75*(if(horiz) strwidth("w") else strheight("x"))

    if (edge.root) {
### FIXME: the first edge + edgetext is drawn here, all others in plotNode()
### -----  maybe use trick with adding a single parent node to the top ?
    x0 <- my.plotNodeLimit(x1, x2, x, center)$x
    if (horiz)
        segments(hgt, x0, yTop, x0)
    else segments(x0, hgt, x0, yTop)
    if (!is.null(et <- attr(x, "edgetext"))) {
        my <- mean(hgt, yTop)
        if (horiz)
        text(my, x0, et)
        else text(x0, my, et)
    }
    }
    my.plotNode(x1, x2, x, type = type, center = center, leaflab = leaflab,
             dLeaf = dLeaf, nodePar = nodePar, edgePar = edgePar, horiz = horiz)
}

### the work horse: plot node (if pch) and lines to all children
my.plotNode <-
    function(x1, x2, subtree, type, center, leaflab, dLeaf,
         nodePar, edgePar, horiz = FALSE)
{
    inner <- !is.leaf(subtree) && x1 != x2
    yTop <- attr(subtree, "height")
    bx <- my.plotNodeLimit(x1, x2, subtree, center)
    xTop <- bx$x
    usrpar <- par("usr");

    ## handle node specific parameters in "nodePar":
    hasP <- !is.null(nPar <- attr(subtree, "nodePar"))
    if(!hasP) nPar <- nodePar

    if(getOption("verbose")) {
    cat(if(inner)"inner node" else "leaf", ":")
    if(!is.null(nPar)) { cat(" with node pars\n"); str(nPar) }
    cat(if(inner)paste(" height", formatC(yTop),"; "),
        "(x1,x2)= (",formatC(x1,wid=4),",",formatC(x2,wid=4),")",
        "--> xTop=", formatC(xTop, wid=8),"\n", sep="")
    }

    Xtract <- function(nam, L, default, indx)
    rep(if(nam %in% names(L)) L[[nam]] else default,
        length.out = indx)[indx]
    asTxt <- function(x) # to allow 'plotmath' labels:
    if(is.character(x) || is.expression(x) || is.null(x)) x else as.character(x)

    i <- if(inner || hasP) 1 else 2 # only 1 node specific par

    if(!is.null(nPar)) { ## draw this node
    pch <- Xtract("pch", nPar, default = 1:2,    i)
    cex <- Xtract("cex", nPar, default = c(1,1),     i)
    col <- Xtract("col", nPar, default = par("col"), i)
    bg <- Xtract("bg", nPar, default = par("bg"), i)
    points(if (horiz) cbind(yTop, xTop) else cbind(xTop, yTop),
           pch = pch, bg = bg, col = col, cex = cex)
    }

    if(leaflab == "textlike")
        p.col <- Xtract("p.col", nPar, default = "white", i)
    lab.col <- Xtract("lab.col", nPar, default = par("col"), i)
    lab.cex <- Xtract("lab.cex", nPar, default = c(1,1), i)
    lab.font <- Xtract("lab.font", nPar, default = par("font"), i)
    if (is.leaf(subtree)) {
    ## label leaf
    if (leaflab == "perpendicular") { # somewhat like plot.hclust
        if(horiz) {
                X <- yTop + dLeaf * lab.cex
                Y <- xTop; srt <- 0; adj <- c(0, 0.5)
        }
        else {
                Y <- yTop - dLeaf * lab.cex
                X <- xTop; srt <- 90; adj <- 1
        }
            nodeText <- asTxt(attr(subtree,"label"))
        text(X, Y, nodeText, xpd = TRUE, srt = srt, adj = adj,
                 cex = lab.cex, col = lab.col, font = lab.font)
    }
    }
    else if (inner) {
    segmentsHV <- function(x0, y0, x1, y1) {
        if (horiz)
        segments(y0, x0, y1, x1, col = col, lty = lty, lwd = lwd)
        else segments(x0, y0, x1, y1, col = col, lty = lty, lwd = lwd)
    }
    for (k in 1:length(subtree)) {
        child <- subtree[[k]]
        ## draw lines to the children and draw them recursively
        yBot <- attr(child, "height")
        if (getOption("verbose")) cat("ch.", k, "@ h=", yBot, "; ")
        if (is.null(yBot))
        yBot <- 0
        xBot <-
        if (center) mean(bx$limit[k:(k + 1)])
        else bx$limit[k] + .my.midDend(child)

        hasE <- !is.null(ePar <- attr(child, "edgePar"))
        if (!hasE)
        ePar <- edgePar
        i <- if (!is.leaf(child) || hasE) 1 else 2
        ## define line attributes for segmentsHV():
        col <- Xtract("col", ePar, default = par("col"), i)
        lty <- Xtract("lty", ePar, default = par("lty"), i)
        lwd <- Xtract("lwd", ePar, default = par("lwd"), i)
        if (type == "triangle") {
        segmentsHV(xTop, yTop, xBot, yBot)
        }
        else { # rectangle
        segmentsHV(xTop,yTop, xBot,yTop)# h
        segmentsHV(xBot,yTop, xBot,yBot)# v
        }
        vln <- NULL
        if (is.leaf(child) && leaflab == "textlike") {
        nodeText <- asTxt(attr(child,"label"))
        if(getOption("verbose"))
            cat('-- with "label"',format(nodeText))
        hln <- 0.6 * strwidth(nodeText, cex = lab.cex)/2
        vln <- 1.5 * strheight(nodeText, cex = lab.cex)/2
        rect(xBot - hln, yBot,
             xBot + hln, yBot + 2 * vln, col = p.col)
        text(xBot, yBot + vln, nodeText, xpd = TRUE,
                     cex = lab.cex, col = lab.col, font = lab.font)
        }
        if (!is.null(attr(child, "edgetext"))) {
        edgeText <- asTxt(attr(child, "edgetext"))
        if(getOption("verbose"))
            cat('-- with "edgetext"',format(edgeText))
        if (!is.null(vln)) {
            mx <-
            if(type == "triangle")
                (xTop+ xBot+ ((xTop - xBot)/(yTop - yBot)) * vln)/2
            else xBot
            my <- (yTop + yBot + 2 * vln)/2
        }
        else {
            mx <- if(type == "triangle") (xTop + xBot)/2 else xBot
            my <- (yTop + yBot)/2
        }
        ## Both for "triangle" and "rectangle" : Diamond + Text

                p.col <- Xtract("p.col", ePar, default = "white", i)
                p.border <- Xtract("p.border", ePar, default = par("fg"), i)
                ## edge label pars: defaults from the segments pars
                p.lwd <- Xtract("p.lwd", ePar, default = lwd, i)
                p.lty <- Xtract("p.lty", ePar, default = lty, i)
                t.col <- Xtract("t.col", ePar, default = col, i)
                t.cex <- Xtract("t.cex", ePar, default =  1,  i)
                t.font<- Xtract("t.font",ePar, default= par("font"), i)
                t.shift <- Xtract("t.shift", ePar, default =  0.01,  i)

        vlm <- strheight(c(edgeText,"h"), cex = t.cex)/2
        hlm <- strwidth (c(edgeText,"m"), cex = t.cex)/2
        hl3 <- c(hlm[1], hlm[1] + hlm[2], hlm[1])
                #polygon(mx+ c(-hl3, hl3), my + sum(vlm)*c(-1:1,1:-1),
                #        col = p.col, border= p.border, lty = p.lty, lwd = p.lwd)
        #text(mx, my, edgeText, cex = t.cex, col = t.col, font = t.font)
                if(horiz) {
                  text(my, mx+t.shift*abs(usrpar[3]-usrpar[4]), edgeText, cex = t.cex, col = t.col, font = t.font)
                } else {
                  text(mx+t.shift*abs(usrpar[2]-usrpar[1]), my, edgeText, cex = t.cex, col = t.col, font = t.font)
                }
        }
        my.plotNode(bx$limit[k], bx$limit[k + 1], subtree = child,
             type, center, leaflab, dLeaf, nodePar, edgePar, horiz)
    }
    }
}

my.plotNodeLimit <- function(x1, x2, subtree, center)
{
    ## get the left borders limit[k] of all children k=1..K, and
    ## the handle point `x' for the edge connecting to the parent.
    inner <- !is.leaf(subtree) && x1 != x2
    if(inner) {
    K <- length(subtree)
    mTop <- .my.memberDend(subtree)
    limit <- integer(K)
    xx1 <- x1
    for(k in 1:K) {
        m <- .my.memberDend(subtree[[k]])
        ##if(is.null(m)) m <- 1
        xx1 <- xx1 + (if(center) (x2-x1) * m/mTop else m)
        limit[k] <- xx1
    }
    limit <- c(x1, limit)
    } else { ## leaf
    limit <- c(x1, x2)
    }
    mid <- attr(subtree, "midpoint")
    center <- center || (inner && !is.numeric(mid))
    x <- if(center) mean(c(x1,x2)) else x1 + (if(inner) mid else 0)
    list(x = x, limit = limit)
}

.my.memberDend <- function(x) {
    r <- attr(x,"x.member")
    if(is.null(r)) {
        r <- attr(x,"members")
        if(is.null(r)) r <- 1:1
    }
    r
}

.my.midDend <- function(x)
    if(is.null(mp <- attr(x, "midpoint"))) 0 else mp


## original Andy Liaw; modified RG, MM :
my.heatmap <- function (x, Rowv=NULL, Colv=if(symm)"Rowv" else NULL,
          distfun = dist, hclustfun = hclust,
          reorderfun = function(d,w) reorder(d,w),
          add.expr, symm = FALSE, revC = identical(Colv, "Rowv"),
          scale = c("row", "column", "none"), na.rm=TRUE,
          margins = c(5, 5), ColSideColors, RowSideColors,
          cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc),
          labRow = NULL, labCol = NULL, main = NULL, xlab = NULL, ylab = NULL,
          keep.dendro = FALSE,
          verbose = getOption("verbose"), imageSize=4, imageVSize=imageSize,imageHSize=imageSize,lasCol=2, lasRow=2, respect=F, ...)
{
    scale <- if(symm && missing(scale)) "none" else match.arg(scale)
    if(length(di <- dim(x)) != 2 || !is.numeric(x))
        stop("'x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if(nr <= 1 || nc <= 1)
        stop("'x' must have at least 2 rows and 2 columns")
    if(!is.numeric(margins) || length(margins) != 2)
        stop("'margins' must be a numeric vector of length 2")

    doRdend <- !identical(Rowv,NA)
    doCdend <- !identical(Colv,NA)
    ## by default order by row/col means
    if(is.null(Rowv)) Rowv <- rowMeans(x, na.rm = na.rm)
    if(is.null(Colv)) Colv <- colMeans(x, na.rm = na.rm)

    ## get the dendrograms and reordering indices

    if(doRdend) {
        if(inherits(Rowv, "dendrogram"))
            ddr <- Rowv
        else {
            hcr <- hclustfun(distfun(x))
            ddr <- as.dendrogram(hcr)
            if(!is.logical(Rowv) || Rowv)
                ddr <- reorderfun(ddr, Rowv)
        }
        if(nr != length(rowInd <- order.dendrogram(ddr)))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else rowInd <- 1:nr

    if(doCdend) {
        if(inherits(Colv, "dendrogram"))
            ddc <- Colv
        else if(identical(Colv, "Rowv")) {
            if(nr != nc)
                stop('Colv = "Rowv" but nrow(x) != ncol(x)')
            ddc <- ddr
        }
        else {
            hcc <- hclustfun(distfun(if(symm)x else t(x)))
            ddc <- as.dendrogram(hcc)
            if(!is.logical(Colv) || Colv)
                ddc <- reorderfun(ddc, Colv)
        }
        if(nc != length(colInd <- order.dendrogram(ddc)))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else colInd <- 1:nc

    ## reorder x
    x <- x[rowInd, colInd]

    labRow <-
        if(is.null(labRow))
            if(is.null(rownames(x))) (1:nr)[rowInd] else rownames(x)
        else labRow[rowInd]
    labCol <-
        if(is.null(labCol))
            if(is.null(colnames(x))) (1:nc)[colInd] else colnames(x)
        else labCol[colInd]

    if(scale == "row") {
        x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
        sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if(scale == "column") {
        x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
        sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }

    ## Calculate the plot layout
    lmat <- rbind(c(NA, 3), 2:1)
    lwid <- c(if(doRdend) 1 else 0.05, imageHSize)
    lhei <- c((if(doCdend) 1 else 0.05) + if(!is.null(main)) 0.2 else 0, imageVSize)
    if(!missing(ColSideColors)) { ## add middle row to layout
        if(!is.character(ColSideColors) || length(ColSideColors) != nc)
            stop("'ColSideColors' must be a character vector of length ncol(x)")
        lmat <- rbind(lmat[1,]+1, c(NA,1), lmat[2,]+1)
        lhei <- c(lhei[1], 0.2, lhei[2])
    }
    if(!missing(RowSideColors)) { ## add middle column to layout
        if(!is.character(RowSideColors) || length(RowSideColors) != nr)
            stop("'RowSideColors' must be a character vector of length nrow(x)")
        lmat <- cbind(lmat[,1]+1, c(rep(NA, nrow(lmat)-1), 1), lmat[,2]+1)
        lwid <- c(lwid[1], 0.2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
    if(verbose) {
        cat("layout: widths = ", lwid, ", heights = ", lhei,"; lmat=\n")
        print(lmat)
    }

    ## Graphics `output' -----------------------

    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = respect)
    ## draw the side bars
    if(!missing(RowSideColors)) {
        par(mar = c(margins[1],0, 0,0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    }
    if(!missing(ColSideColors)) {
        par(mar = c(0.5,0, 0,margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    ## draw the main carpet
    par(mar = c(margins[1], 0, 0, margins[2]))
    if(!symm || scale != "none")
        x <- t(x)
    if(revC) { # x columns reversed
        iy <- nr:1
        ddr <- rev(ddr)
        x <- x[,iy]
    } else iy <- 1:nr

    image(1:nc, 1:nr, x, xlim = 0.5+ c(0, nc), ylim = 0.5+ c(0, nr),
          axes = FALSE, xlab = "", ylab = "", ...)
    axis(1, 1:nc, labels= labCol, las= lasCol, line= -0.5, tick= 0, cex.axis= cexCol)
    if(!is.null(xlab)) mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels= labRow, las= lasRow, line= -0.5, tick= 0, cex.axis= cexRow)
    if(!is.null(ylab)) mtext(ylab, side = 4, line = margins[2] - 1.25,las=lasRow)
    if (!missing(add.expr))
        eval(substitute(add.expr))

    ## the two dendrograms :
    par(mar = c(margins[1], 0, 0, 0))
    if(doRdend)
        my.plot.dendrogram(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    else frame()

    par(mar = c(0, 0, if(!is.null(main)) 1 else 0, margins[2]))
    if(doCdend)
        my.plot.dendrogram(ddc,               axes = FALSE, xaxs = "i", leaflab = "none")
    else if(!is.null(main)) frame()

    ## title
    if(!is.null(main)) title(main, cex.main = 1.5*op[["cex.main"]])

    invisible(list(rowInd = rowInd, colInd = colInd,
                   Rowv = if(keep.dendro && doRdend) ddr,
                   Colv = if(keep.dendro && doCdend) ddc ))
}
