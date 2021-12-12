require(mgcv)
require(Chicago)
require(gamlss)
require(gamlss.add)
require(ggplot2)
require(dplyr)
require(BSgenome.Hsapiens.UCSC.hg19)

ChiCAGOPreProcessing<-function(ChiCAGOData,restrictionFeatures,Dmax=2e6){
  
  #Preprocessed ChiCAGO Subsetted for Storage can be accepted as input
  if(class(ChiCAGOData)[1]=="chicagoData"){
    processed=attr(ChiCAGOData,'x')%>%dplyr::filter(otherEndChr==baitChr & abs(distSign)<=Dmax)
  }else{
    processed=ChiCAGOData%>%dplyr::filter(otherEndChr==baitChr & abs(distSign)<=Dmax)
  }
  
  #Add information for each end
  processed=data.frame(processed,restrictionFeatures[processed$otherEndID,c('chr','start','end')])%>%data.table
  
  processed$gc=restrictionFeatures[processed$otherEndID]$gc
  processed$map=restrictionFeatures[processed$otherEndID]$map
  processed$gene=restrictionFeatures[(processed$baitID)]$gene
  
  #Demarcate Sides
  processed$side='Left'
  processed$side[processed$distSign>0]='Right'
  processed$side=factor(processed$side)
  
  processed$baitID=as.factor(processed$baitID)
  processed$baitID_side=factor(paste0(processed$baitID,processed$side))
  
  #Calculate Distance from the bait
  processed$dist=abs(processed$distSign)
  
  processed
}

ChiCModelRun<-function(processed,mode='baits',model='nb',
                       allBaits=F,baitIDs=NA,
                       extraVars=NA,k=4,cores=1,
                       outputDir='./'){
  dimK<<-k
  
  #Model Terms
  if(model=='nb'){
    mu.form=as.formula(N~s(log(dist),by=side,k=dimK,bs='ts')+gc+map+side+otherEndLen+rep)
  }else if(model=='nb_vardisp'){
    mu.form=as.formula(N~ga(~s(log(dist),by=side,k=dimK,bs='ts')+gc+map+side+otherEndLen+rep,select=T))
    sigma.form=as.formula(~ga(~s(log(dist),k=dimK,bs='ts')+gc+map+side+otherEndLen+rep,select=T))
  }
  #Simplify model for some cutoff based on available data
  nMin=ceiling(quantile((processed%>%group_by(baitID)%>%summarise(nEnds=n()))$nEnds,probs=0.1))
  
  #Generate Replicate Information Strings
  nRep=sum(grepl("N\\.",colnames(processed)))
  repString=paste0('N.',c(1:nRep))
  
  if(is.na(baitIDs)){
    baitIDs=unique(processed$baitID)
  }
  
  #Primary loop through all baits
  results=mclapply(baitIDs,function(x){
    baitData=processed%>%dplyr::filter(baitID==x)
    
    if(nrow(baitData)<nMin){
      sigma.form=as.formula(~1)
    }
    
    #Need to expand out data frames per replicate
    dat=generateChICReplicateFrame(baitData,repString)
    
    #Run model
    if(model=='nb'){
      results=ChiCModel(dat,mu.form)
    }else if(model=='nb_vardisp'){
      results=ChiCModelVarDisp(dat,mu.form,sigma.form)
    }
    results$gene=unique(baitData$gene)
    
    rm(baitData)
    rm(dat)
    
    results
  },mc.cores=cores)
  
  return(rbindlist(results))
}

ChiCModelVarDisp<-function(dat,mu.form,sigma.form){
  #Standard modeling variables
  vars=c('baitChr','start','end','side','otherEndLen','dist','gc','map','nTrans','rep','baitID_side','baitID')
  
  dat<<-dat
  
  #Preliminary fits check if we need variable disperion modes
  init.mod=tryCatch({
    gamlss(mu.form,
           sigma.formula = sigma.form,
           family=NBI,data=dat,method=mixed(50,50))
  },error=function(cond){
    gamlss(mu.form,
           family=NBI,data=dat,method=mixed(50,50))
  })
  
  init.mod2=gamlss(mu.form,
                   family=NBI,data=dat,method=mixed(50,50))
  
  #VC tests on fits 
  if(myVC.test(init.mod,init.mod2)==2){
    sigma.form=as.formula(~1)
    init.mod=init.mod2
  }
  
  #Accept-reject at 0.975
  new.dat<<-remove_outliers_nb_varDisp(dat,init.mod)
  
  #One final fit with Poisson fall back if error is thrown
  mod=tryCatch({
    gamlss(mu.form,
           sigma.formula = sigma.form,
           family=NBI,data=new.dat,method=mixed(50,50))
  },
  error=function(cond) {
    gamlss(mu.form,
           family=PO,data=new.dat,method=mixed(50,50))
  })
  
  #Check for convergence, else increase iterations
  if(!mod$converged){
    mod=gamlss(mu.form,
               family=NBI,data=new.dat,method=mixed(50,50))
  }
  mod<<-mod
  
  #Predict on all terms
  dat$mu=predict(mod,newdata = dat[,..vars],type='response')
  
  if(mod$family[[1]]=='NBI'){
    dat$disp=predict(mod,newdata = dat[,..vars],type='response',what='sigma')
    dat$pval=pnbinom(dat$N-1,mu=dat$mu,size=1/dat$disp,lower.tail=F)
  }else if(mod$family[[1]]=='PO'){
    dat$disp=NA
    dat$pval=ppois(dat$N-1,lambda = dat$mu, lower.tail=F)
  }
  
  dat$qval=p.adjust(dat$pval,method='BH')
  
  #Pool using Fisher's method and compute 
  dat<-fishersMid(dat)
  
  rm(new.dat)
  rm(init.mod)
  rm(mod)
  
  dat
}

ChiCModel<-function(dat,mu.form){
  #Standard modeling variables
  vars=c('baitChr','start','end','side','otherEndLen','dist','gc','map','nTrans','rep','baitID_side','baitID')
  
  dat<<-dat
  
  #Preliminary fits check if we need variable disperion mode
  init.mod=gam(mu.form,family=nb(),data=dat)
  
  #Accept-reject at 0.975
  new.dat<<-remove_outliers_nb(dat,init.mod)
  
  #One final fit with Poisson fall back if error is thrown
  mod=tryCatch({
    gam(mu.form,family=nb(),data=new.dat)
  },
  error=function(cond) {
    gamlss(mu.form,family=poisson(),data=new.dat)
  })
  
  mod<<-mod
  
  #Predict on all terms
  dat$mu=predict(mod,newdata = dat[,..vars],type='response')
  
  if(mod$family$family!='poisson'){
    dat$disp=exp(mod$family$getTheta())
    dat$pval=pnbinom(dat$N-1,mu=dat$mu,size=dat$disp,lower.tail=F)
  }else if(mod$family$family=='poisson'){
    dat$disp=NA
    dat$pval=ppois(dat$N-1,lambda = dat$mu, lower.tail=F)
  }
  
  dat$qval=p.adjust(dat$pval,method='BH')
  
  #Pool using Fisher's method and compute 
  dat<-fishersMid(dat)
  
  rm(new.dat)
  rm(init.mod)
  rm(mod)
  
  dat
}

generateChICReplicateFrame<-function(baitData,repString){
  #Generate primary dataframe
  extractVars=c('baitChr','start','end','side','otherEndLen','dist','gc','map','nTrans','baitID_side','baitID',repString)
  baitData=data.frame(baitData)
  data=baitData[,extractVars]
  
  #Add replicate level information and filter for fragments close to bait
  data=rbindlist(lapply(c(1:length(repString)),function(x){
    data$rep=x
    data$N=data[,repString[x]]
    data
  }))%>%dplyr::filter(dist>5000)
  
  data$rep=factor(data$rep)
  
  data
}

plotResults<-function(results,baitIDs=NA,outputDir="./"){
  #Option to plot only specified baits
  if(is.na(baitIDs)){
    baitIDs=unique(results$baitID)
  }
  #Loop through all baits into pdf file
  for(i in baitIDs){
    baitData=results%>%dplyr::filter(baitID==i)
    if(is.na(unique(baitData$gene))){
      pdf(paste0(outputDir,unique(baitData$baitID),'Results.pdf'))
    }else{
      pdf(paste0(outputDir,unique(baitData$gene),'Results.pdf'))
    }
    
    print(plot(baitData$start,baitData$N))
    print(plot(baitData$start,baitData$mu))
    print(plot(baitData$start,baitData$disp))
    print(plot(baitData$start,-log(baitData$pfish,10)))
    print(plot(baitData$start,-log(baitData$padj_fish,10)))
    dev.off()
  }
}
##########################

#Utilities
fishersMid<-function(x){
  #Compute pooled p-values
  
  x=x%>%group_by(dist,side)%>%dplyr::filter(sum(N)>=0& dplyr::n()==length(unique(x$rep)))%>%arrange((dist))
  x=x%>%data.table
  x$key=paste0(x$baitChr,x$start,x$end)
  x0=x%>%data.table
  setkey(x,key)
  setkey(x0,key)
  
  temp=x%>%group_by(baitChr,start,end,dist,side)%>%
    summarise(pfish=sum(log(pval)),count=dplyr::n())%>%
    dplyr::filter(count>=length(unique(x$rep)))%>%
    mutate(pfishVal=pfish*-2)
  temp$pfish=pchisq(temp$pfishVal,df=temp$count*2,lower.tail=F)
  
  temp$padj_fish=p.adjust(temp$pfish,method='BH')
  temp$key=paste0(temp$baitChr,temp$start,temp$end)
  temp=temp%>%data.table
  setkey(temp,key)
  
  x0$pfish=1
  x0[temp$key,]$pfish=rep(temp[temp$key]$pfish,each=length(unique(x0$rep)))
  
  x0$padj_fish=1
  x0[temp$key,]$padj_fish=rep(temp[temp$key]$padj_fish,each=length(unique(x0$rep)))
  
  x0
}
getGC<-function(ranges,BSgenome=BSgenome.Hsapiens.UCSC.hg19){
  #Calculate GC content for each position in GRanges object
  freqs <- alphabetFrequency(getSeq(BSgenome, ranges))
  gc <- (freqs[,'C'] + freqs[,'G'])/rowSums(freqs)
}
generateFeatureMap<-function(rMapFile,baitMapFile,BSgenome=BSgenome.Hsapiens.UCSC.hg19,mapFile=NA){
  #Bait Map File Read
  baitMap=fread(baitMapFile)
  colnames(baitMap)=c('chr','start','end','baitID','gene')
  baitMap$baitID=as.character(baitMap$baitID)
  setkey(baitMap,baitID)
  
  #Restriction Map Read
  rMap=fread(rMapFile)
  colnames(rMap)=c('chr','start','end','rID')
  rMap$chr=paste0('chr',rMap$chr)
  setkey(rMap,rID)
  
  #Add GC content
  ranges=makeGRangesFromDataFrame(rMap[,1:3])
  
  rMap$gc=getGC(ranges)
  
  #Take sums over mappability scores
  mapScores=import.bw(mapFile,which=ranges)
  overlaps=findOverlaps(ranges,mapScores)
  rMap$map=0
  rMap$map[rMap$rID%in%(unique(queryHits(overlaps)))]=(aggregate(mapScores,overlaps,score=sum(score)))$score
  
  rMap$gene='None'
  rMap[as.numeric(baitMap$baitID),]$gene=baitMap$gene
  
  rMap
}
myVC.test<-function (obj1, obj2, sig.lev = 0.05){
  #Sourced from gamlss::VC.test, modified to provide a decision rule
  if (!is.gamlss(obj1))
    stop(paste("Object 1 is not a gamlss object", "\n", ""))
  if (!is.gamlss(obj2))
    stop(paste("Object 2 is not a gamlss object", "\n", ""))
  l1 <- logLik(obj1)[1]
  l2 <- logLik(obj2)[1]
  if (l1 == l2){
    return(2)
  }
  
  p1 <- obj1$df.fit
  p2 <- obj2$df.fit
  n <- obj1$N
  pobj1 <- predictAll(obj1)
  pobj2 <- predictAll(obj2)
  fname1 <- obj1$family[1]
  fname2 <- obj2$family[1]
  nopar1 <- eval(gamlss.family(fname1))$nopar
  nopar2 <- eval(gamlss.family(fname2))$nopar
  pdfName1 <- paste("d", fname1, sep = "")
  pdfName2 <- paste("d", fname2, sep = "")
  binomial1 <- fname1 %in% .gamlss.bi.list
  binomial2 <- fname2 %in% .gamlss.bi.list
  model1 <- deparse(substitute(obj1))
  model2 <- deparse(substitute(obj2))
  if (binomial1) {
    linc1 <- switch(nopar1, eval(call(pdfName1, pobj1$y,
                                      bd = obj1$bd, mu = pobj1$mu, log = TRUE)), eval(call(pdfName1,
                                                                                           pobj1$y, bd = obj1$bd, mu = pobj1$mu, sigma = pobj1$sigma,
                                                                                           log = TRUE)), eval(call(pdfName1, pobj1$y, bd = obj1$bd,
                                                                                                                   mu = pobj1$mu, sigma = pobj1$sigma, nu = pobj1$nu,
                                                                                                                   log = TRUE)), eval(call(pdfName1, pobj1$y, bd = obj1$bd,
                                                                                                                                           mu = pobj1$mu, sigma = pobj1$sigma, nu = pobj1$nu,
                                                                                                                                           tau = pobj1$tau, log = TRUE)))
  }
  else {
    linc1 <- switch(nopar1, eval(call(pdfName1, pobj1$y,
                                      mu = pobj1$mu, log = TRUE)), eval(call(pdfName1,
                                                                             pobj1$y, mu = pobj1$mu, sigma = pobj1$sigma, log = TRUE)),
                    eval(call(pdfName1, pobj1$y, mu = pobj1$mu, sigma = pobj1$sigma,
                              nu = pobj1$nu, log = TRUE)), eval(call(pdfName1,
                                                                     pobj1$y, mu = pobj1$mu, sigma = pobj1$sigma,
                                                                     nu = pobj1$nu, tau = pobj1$tau, log = TRUE)))
  }
  if (binomial2) {
    linc2 <- switch(nopar1, eval(call(pdfName2, pobj2$y,
                                      bd = obj2$bd, mu = pobj2$mu, log = TRUE)), eval(call(pdfName2,
                                                                                           pobj2$y, bd = obj2$bd, mu = pobj2$mu, sigma = pobj2$sigma,
                                                                                           log = TRUE)), eval(call(pdfName2, pobj2$y, bd = obj2$bd,
                                                                                                                   mu = pobj2$mu, sigma = pobj2$sigma, nu = pobj2$nu,
                                                                                                                   log = TRUE)), eval(call(pdfName2, pobj2$y, bd = obj2$bd,
                                                                                                                                           mu = pobj2$mu, sigma = pobj2$sigma, nu = pobj2$nu,
                                                                                                                                           tau = pobj2$tau, log = TRUE)))
  }
  else {
    linc2 <- switch(nopar2, eval(call(pdfName2, pobj2$y,
                                      mu = pobj2$mu, log = TRUE)), eval(call(pdfName2,
                                                                             pobj2$y, mu = pobj2$mu, sigma = pobj2$sigma, log = TRUE)),
                    eval(call(pdfName2, pobj2$y, mu = pobj2$mu, sigma = pobj2$sigma,
                              nu = pobj2$nu, log = TRUE)), eval(call(pdfName2,
                                                                     pobj2$y, mu = pobj2$mu, sigma = pobj2$sigma,
                                                                     nu = pobj2$nu, tau = pobj2$tau, log = TRUE)))
  }
  li12 <- linc1 - linc2
  w <- sqrt(var(li12) * (n - 1))
  vt <- (l1 - l2 - (p1 - p2)/2 * log(n))/w
  li12b <- li12 - (p1 - p2)/(2 * n) * log(n)
  b <- sum(li12b > 0)
  
  if(vt > qnorm(1 - sig.lev/2) & b >= n/2){
    return(1)
  }else{
    return(2)
  }
}

remove_outliers_nb_varDisp <- function(dat, mod, cut=0.975) {
  #Accept and reject during null estimation
  dat$q <- qNBI(cut, sigma = mod$sigma.fv, mu = mod$mu.fv)
  new.dat <- dat %>% dplyr::filter(.data$N <= .data$q)
  return(new.dat)
}

remove_outliers_nb <- function(dat, mod, cut=0.975) {
  #Accept and reject during null estimation
  dat$q <- qnbinom(cut, mu = mod$fitted.values,size=exp(mod$family$getTheta()))
  new.dat <- dat %>% dplyr::filter(.data$N <= .data$q)
  return(new.dat)
}



