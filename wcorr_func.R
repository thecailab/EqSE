
weight.calc<-function(count, K=2, epsilon=1000, nc=1, celltype=NULL){

  library("zinbwave")
  library("DESeq2")

  if(is.null(celltype)){
    celltype<-rep("cell1",ncol(count))
    covar<-1
  }else{
    covar<-"cluster"
  }

  coldata<-data.frame(cluster=factor(celltype))
  rownames(coldata) <- colnames(count)

  #zinb.wave
  exp.obj <- DESeqDataSetFromMatrix(countData = count, colData=coldata, design=as.formula(paste0("~",covar)))
  zinb.fit <- zinbFit(exp.obj, K=K, epsilon=epsilon,BPPARAM=BiocParallel::MulticoreParam(nc))
  zinb.wave <- zinbwave(exp.obj, fitted_model = zinb.fit, K = 2,
                           epsilon=1000,
                           observationalWeights = TRUE)
  w<- assay(zinb.wave, "weights")

  return(w)
}


wcorr.calc<-function(dat.exp, dat.histone, w.exp, w.histone, alpha=3, method="spearman", celltype=NULL, type="GT"){

  library("wCorr")

  #match matrix
  if(type=="GT"){
    in.histone<-which(sub(".+_","",rownames(dat.histone)) %in% rownames(dat.exp))
    dat.histone<-dat.histone[in.histone,,drop=FALSE]
    w.histone<-w.histone[in.histone,,drop=FALSE]

    match.exp<-match(sub(".+_","",rownames(dat.histone)),rownames(dat.exp))
    dat.exp<-dat.exp[match.exp,,drop=FALSE]
    w.exp<-w.exp[match.exp,,drop=FALSE]
  }

  w2<-(w.histone*w.exp)^alpha

  if(is.null(celltype)){
    celltype<-factor(rep("cell1",ncol(dat.exp)))
  }else{
    celltype<-factor(celltype)
  }
  celltype.level<-levels(celltype)

  out<-list()

  for(j in 1:length(celltype.level)){
    celltype.j<-celltype.level[j]
    pos<-which(celltype==celltype.j)

    out.m<-matrix(NA, nrow(dat.exp),2)

    for(i in 1:nrow(dat.exp)){
      exp.i<-dat.exp[i,]
      histone.i<-dat.histone[i,]
      w2.i<-w2[i,]

      out.m[i,1]<-weightedCorr(exp.i[pos],histone.i[pos],weight=w2.i[pos], method=method)
      out.m[i,2]<-length(which(w2.i[pos]>0.5))
      colnames(out.m)<-paste(c("Rho","n"),celltype.j,sep="_")

      out[[celltype.j]]<-out.m
    }
  }

  out<-do.call("cbind",out)
  out<-cbind(rownames(dat.histone),rownames(dat.exp),out)
  colnames(out)[1:2]<-c("histone.peak","gene")

  return(out)
}


wcorr.plot<-function(dat.exp, dat.histone, w.exp, w.histone, alpha=3, celltype=NULL, cell=NULL, peak=NULL, color="red", pch=20, cex=1){

  #match matrix

  if(is.null(peak)){
    stop("the peak need to be specified")
  }

  if(is.null(cell)){
    pos<-1:ncol(dat.exp)
  }else{
    if(is.null(celltype)){
      stop("the celltype information need to be provided")
    }
    pos<-which(celltype==cell)
  }

  row.histone<-which(rownames(dat.histone)==peak)
  histone.i<-dat.histone[row.histone,pos]
  w.histone.i<-w.histone[row.histone,pos]

  row.exp<-which(rownames(dat.exp)==sub(".+_","",peak))[1]
  exp.i<-dat.exp[row.exp,pos]
  w.exp.i<-w.exp[row.exp,pos]

  w2.i<-(w.histone.i*w.exp.i)^alpha

  cat.i<-as.numeric(cut(w2.i,breaks = 10))

  color.i<-colorRampPalette(c("grey",color))(10)[cat.i]

  plot(histone.i,exp.i, pch=pch, cex=cex, col=color.i,xlab="Histone abundance",ylab="Expression", cex.lab=1.5, cex.axis=1.3)
  abline(lm(exp.i~histone.i,weight=w2.i),col = color.i,lwd=2)
}
