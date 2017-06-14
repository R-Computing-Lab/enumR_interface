  #-------------------------------------
  # Library
  #-------------------------------------
library(psych)
library(GAIPE)
library(GPArotation)
  library(beepr)
  #enumR functions, will load to package later
 # source("C:\\Dropbox\\Lab\\zSoftware\\Github\\enumR\\R\\enumR.R")
  enumR<-function(x, nfactors = 8, rotation = "oblimin", diagonal = FALSE, method = "mle",samplesize=NULL,use="pairwise",cor="cov",truemodel=FALSE,...){
    require(psych)
    require(GAIPE)
    n=nfactors; fm=method
    # cl <- match.call()
    if (rotation == "oblimin") {
      if (!requireNamespace("GPArotation")) {
        stop("You must have GPArotation installed to use oblimin rotation")
      }
    }
    
    old_rotation = rotation
    
    if(truemodel){
      complexrow <- function(x, c) {
        n = length(x)
        temp <- x
        x <- rep(0, n)
        for (j in 1:c) {
          locmax <- which.max(abs(temp))
          x[locmax] <- sign(temp[locmax]) * max(abs(temp))
          temp[locmax] <- 0
        }
        return(x)
      }
      complexmat <- function(x, c) {
        nrows <- dim(x)[1]
        ncols <- dim(x)[2]
        for (i in 1:nrows) {
          x[i, ] <- complexrow(x[i, ], c)
        }
        return(x)
      }
      map <- function(x, n) {
        nvar <- dim(x)[2]
        min.partial <- rep(NA, n)
        e <- eigen(x)
        evect <- e$vectors
        comp <- evect %*% diag(sqrt(e$values))
        if (n >= nvar) {
          n1 <- nvar - 1
        } else {
          n1 <- n
        }
        
        c11.star <- x - comp[, 1:n1] %*% t(comp[, 1:n1])
        d <- diag(1/sqrt(diag(c11.star)))
        rstar <- d %*% c11.star %*% d
        diag(rstar) <- 0
        min.partial <- sum(rstar * rstar)/(nvar * (nvar -1))
        
        return(min.partial)
      }
      if (dim(x)[2] < n)
        n <- dim(x)[2]
      complexfit <- array(0, dim = c(1, 1))
      complexresid <- array(0, dim = c(1, 1))
      vss.df <- data.frame(dof = 0, chisq = NA, prob = NA,sqresid = NA, fit = NA, RMSEA = NA,RMSEA_lower=NA,RMSEA_upper=NA, BIC = NA, SABIC = NA, null.model = NA, null.dof = NA, null.chisq = NA,complex = NA, eChisq = NA, SRMR = NA, eCRMS = NA, eBIC = NA, TLI=NA,fa.fit=NA,fa.fit.off=NA, objective=NA,ENull=NA,eSABIC=NA,null.chisq1=NA,TLI.m=NA, NFI.m=NA, CFI.m=NA,RMSEA_lower.m=NA,RMSEA_upper.m=NA)
      
      if (!is.matrix(x))
        x <- as.matrix(x)
      
      if (is.null(samplesize)) {
        message("samplesize was not specified and was arbitrarily set to 1000.  This only affects the chi square values.")
        samplesize <- 1000}
      
      map.values <- map(x, n)
      if (n > dim(x)[2]) {
        n <- dim(x)[2]}
      
      
      PHI <- diag(n)
      if (n < 2) {
        (rotation = "none")
      }
      else {
        rotation = old_rotation
      }
      
      f <- fa(x, n, rotate = rotation, n.obs = samplesize, warnings = FALSE,
              fm = fm, scores = "none", cor = cor, ...)
      
      original <- x
      sqoriginal <- original * original
      totaloriginal <- sum(sqoriginal) - diagonal *
        sum(diag(sqoriginal))
      
      
      load <- as.matrix(f$loadings)
      model <- load %*% PHI %*% t(load)
      residual <- original - model
      sqresid <- residual * residual
      totalresid <- sum(sqresid) - diagonal * sum(diag(sqresid))
      fit <- 1 - totalresid/totaloriginal
      
      vss.df[ "dof"] <- f$dof
      vss.df[ "chisq"] <- f$STATISTIC
      vss.df[ "prob"] <- f$PVAL
      vss.df["eChisq"] <- f$chi
      vss.df["ENull"] <- f$ENull
      vss.df["SRMR"] <- f$rms
      vss.df["eRMS"] <- f$rms
      vss.df["eCRMS"] <- f$crms
      vss.df["eBIC"] <- f$EBIC
      vss.df["eSABIC"] <- f$ESABIC
      
      vss.df["TLI"]<-f$TLI
      vss.df["TLI.m"]<-(f$null.chisq/f$null.dof - f$STATISTIC/f$dof)/(f$null.chisq/f$null.dof - 1)
      
      vss.df["NFI.m"]<-(f$null.chisq - f$STATISTIC)/(f$null.chisq)
      
      vss.df["CFI.m"]<-(f$null.chisq-f$null.dof - f$STATISTIC-f$dof)/(f$null.chisq-f$null.dof)
      
      vss.df["CFI.m"]<-(f$null.chisq-f$null.dof - f$STATISTIC-f$dof)/(f$null.chisq-f$null.dof)
      
      vss.df["null.model"] = f$null.model
      vss.df["null.dof"] = f$null.dof
      
      vss.df["fa.fit"] = f$fit
      vss.df["fa.fit.off"] = f$fit.off
      vss.df["objective"] = f$objective
      if (!is.null(f$chisq)) {
        vss.df["null.chisq1"] = f$chisq
        vss.df["null.chisq"] = f$null.chisq
      }
      if (!is.null(f$RMSEA[1])) {
        vss.df[ "RMSEA"] <- f$RMSEA[1]
        vss.df[ "RMSEA_upper"] <- f$RMSEA[3]
        vss.df[ "RMSEA_lower"] <- f$RMSEA[2]
        
        RMSEA.m<-CI.RMSEA(f$RMSEA[1],N=samplesize,df=f$dof,clevel=1-f$RMSEA[4])
        
        vss.df[ "RMSEA_upper.m"] <- RMSEA.m$Upper.CI
        vss.df[ "RMSEA_lower.m"] <- RMSEA.m$Lower.CI
        
      }
      else {
        vss.df[ "RMSEA"] <- NA
      }
      
      if (!is.null(f$BIC)) {
        vss.df[ "BIC"] <- f$BIC
      }
      else {
        vss.df[ "BIC"] <- NA
      }
      if (!is.null(f$SABIC)) {
        vss.df[ "SABIC"] <- f$SABIC
      }
      else {
        vss.df[ "SABIC"] <- NA
      }
      if (!is.null(f$complexity)) {
        vss.df[ "complex"] <- mean(f$complexity)
      }
      else {
        vss.df[ "complex"] <- NA
      }
      
      vss.df[ "sqresid"] <- totalresid
      vss.df[ "fit"] <- fit
      for (c in 1:n) {
        simpleload <- complexmat(load, c)
        model <- simpleload %*% PHI %*% t(simpleload)
        residual <- original - model
        sqresid <- residual * residual
        totalsimple <- sum(sqresid) - diagonal * sum(diag(sqresid))
        simplefit <- 1 - totalsimple/totaloriginal
        complexresid[c] <- totalsimple
        complexfit[c] <- simplefit
      }
      
    }else{
      complexrow <- function(x, c) {
        n = length(x)
        temp <- x
        x <- rep(0, n)
        for (j in 1:c) {
          locmax <- which.max(abs(temp))
          x[locmax] <- sign(temp[locmax]) * max(abs(temp))
          temp[locmax] <- 0
        }
        return(x)
      }
      complexmat <- function(x, c) {
        nrows <- dim(x)[1]
        ncols <- dim(x)[2]
        for (i in 1:nrows) {
          x[i, ] <- complexrow(x[i, ], c)
        }
        return(x)
      }
      map <- function(x, n) {
        nvar <- dim(x)[2]
        min.partial <- rep(NA, n)
        e <- eigen(x)
        evect <- e$vectors
        comp <- evect %*% diag(sqrt(e$values))
        if (n >= nvar) {
          n1 <- nvar - 1
        } else {
          n1 <- n
        }
        for (i in 1:n1) {
          c11.star <- x - comp[, 1:i] %*% t(comp[, 1:i])
          d <- diag(1/sqrt(diag(c11.star)))
          rstar <- d %*% c11.star %*% d
          diag(rstar) <- 0
          min.partial[i] <- sum(rstar * rstar)/(nvar * (nvar -1))
        }
        return(min.partial)
      }
      if (dim(x)[2] < n)
        n <- dim(x)[2]
      complexfit <- array(0, dim = c(n, n))
      complexresid <- array(0, dim = c(n, n))
      vss.df <- data.frame(dof = rep(0, n), chisq = NA, prob = NA,sqresid = NA, fit = NA, RMSEA = NA,RMSEA_lower=NA,RMSEA_upper=NA, BIC = NA, SABIC = NA, null.model = NA, null.dof = NA, null.chisq = NA,complex = NA, eChisq = NA, SRMR = NA, eCRMS = NA, eBIC = NA, TLI=NA,fa.fit=NA,fa.fit.off=NA, objective=NA,ENull=NA,eSABIC=NA,null.chisq1=NA,TLI.m=NA, NFI.m=NA, CFI.m=NA,RMSEA_lower.m=NA,RMSEA_upper.m=NA)
      
      if (!is.matrix(x))
        x <- as.matrix(x)
      
      if (is.null(samplesize)) {
        message("samplesize was not specified and was arbitrarily set to 1000.  This only affects the chi square values.")
        samplesize <- 1000}
      
      map.values <- map(x, n)
      if (n > dim(x)[2]) {
        n <- dim(x)[2]}
      
      for (i in 1:n) {
        PHI <- diag(i)
        if (i < 2) {
          (rotation = "none")
        }
        else {
          rotation = old_rotation
        }
        
        f <- fa(x, i, rotate = rotation, n.obs = samplesize, warnings = FALSE,
                fm = fm, scores = "none", cor = cor, ...)
        if (i == 1) {
          original <- x
          sqoriginal <- original * original
          totaloriginal <- sum(sqoriginal) - diagonal *
            sum(diag(sqoriginal))
        }
        
        load <- as.matrix(f$loadings)
        model <- load %*% PHI %*% t(load)
        residual <- original - model
        sqresid <- residual * residual
        totalresid <- sum(sqresid) - diagonal * sum(diag(sqresid))
        fit <- 1 - totalresid/totaloriginal
        
        vss.df[i, "dof"] <- f$dof
        vss.df[i, "chisq"] <- f$STATISTIC
        vss.df[i, "prob"] <- f$PVAL
        vss.df[i, "eChisq"] <- f$chi
        vss.df[i, "ENull"] <- f$ENull
        vss.df[i, "SRMR"] <- f$rms
        vss.df[i, "eRMS"] <- f$rms
        vss.df[i, "eCRMS"] <- f$crms
        vss.df[i, "eBIC"] <- f$EBIC
        vss.df[i, "eSABIC"] <- f$ESABIC
        
        vss.df[i,"TLI"]<-f$TLI
        vss.df[i,"TLI.m"]<-(f$null.chisq/f$null.dof - f$STATISTIC/f$dof)/(f$null.chisq/f$null.dof - 1)
        
        vss.df[i,"NFI.m"]<-(f$null.chisq - f$STATISTIC)/(f$null.chisq)
        
        vss.df[i,"CFI.m"]<-(f$null.chisq-f$null.dof - f$STATISTIC-f$dof)/(f$null.chisq-f$null.dof)
        
        vss.df[i,"CFI.m"]<-(f$null.chisq-f$null.dof - f$STATISTIC-f$dof)/(f$null.chisq-f$null.dof)
        
        vss.df[i,"null.model"] = f$null.model
        vss.df[i,"null.dof"] = f$null.dof
        
        vss.df[i,"fa.fit"] = f$fit
        vss.df[i,"fa.fit.off"] = f$fit.off
        vss.df[i,"objective"] = f$objective
        if (!is.null(f$chisq)) {
          vss.df[i,"null.chisq1"] = f$chisq
          vss.df[i,"null.chisq"] = f$null.chisq
        }
        if (!is.null(f$RMSEA[1])) {
          vss.df[i, "RMSEA"] <- f$RMSEA[1]
          vss.df[i, "RMSEA_upper"] <- f$RMSEA[3]
          vss.df[i, "RMSEA_lower"] <- f$RMSEA[2]
          
         # RMSEA.m<-CI.RMSEA(f$RMSEA[1],N=samplesize,df=f$dof,clevel=1-f$RMSEA[4])
          
       #   vss.df[i, "RMSEA_upper.m"] <- RMSEA.m$Upper.CI
      #    vss.df[i, "RMSEA_lower.m"] <- RMSEA.m$Lower.CI
          
        }
        else {
          vss.df[i, "RMSEA"] <- NA
        }
        
        if (!is.null(f$BIC)) {
          vss.df[i, "BIC"] <- f$BIC
        }
        else {
          vss.df[i, "BIC"] <- NA
        }
        if (!is.null(f$SABIC)) {
          vss.df[i, "SABIC"] <- f$SABIC
        }
        else {
          vss.df[i, "SABIC"] <- NA
        }
        if (!is.null(f$complexity)) {
          vss.df[i, "complex"] <- mean(f$complexity)
        }
        else {
          vss.df[i, "complex"] <- NA
        }
        
        vss.df[i, "sqresid"] <- totalresid
        vss.df[i, "fit"] <- fit
        for (c in 1:i) {
          simpleload <- complexmat(load, c)
          model <- simpleload %*% PHI %*% t(simpleload)
          residual <- original - model
          sqresid <- residual * residual
          totalsimple <- sum(sqresid) - diagonal * sum(diag(sqresid))
          simplefit <- 1 - totalsimple/totaloriginal
          complexresid[i, c] <- totalsimple
          complexfit[i, c] <- simplefit
        }
      }
    }
    vss.stats <- data.frame(vss.df, cfit = complexfit, cresidual = complexresid,map.values)
    
    #class(vss.results) <- c("psych", "vss")
    if(truemodel){
      return(vss.stats[1,])
    }else{return(vss.stats)}
    #return(f)
  }
  


#source("C:\\Dropbox\\Lab\\zSoftware\\Github\\enumR\\R\\GenData.R")
  GenData = function(fmodel,effect,samplesize,inames=NULL) {   # define a general function in terms of a factor model and an effects matrix
    
    numberofvariables = dim(fmodel)[1]        #problem size determined by input to the function
    if(is.null(inames)){inames=paste0("item_",1:numberofvariables)}
    numberoflatent = dim(fmodel)[2]
    tmodel = t(fmodel)      #transpose of model
    communality=diag(fmodel%*%tmodel)       #find how much to weight true scores and errors given the measurement model
    uniqueness=1-communality
    errorweight=diag(sqrt(uniqueness))      #how much to weight the errors
    
    latentscores = matrix(rnorm(samplesize*(numberoflatent)),samplesize) #create true scores for the latent variables
    latentscores = latentscores%*%effect
    truescores = latentscores%*%tmodel
    error = matrix(rnorm(samplesize*(numberofvariables)),samplesize)  #create normal error scores
    error = error%*%errorweight
    observedscore = truescores+error
    observedscore = data.frame(observedscore)
    names(observedscore)<-inames
    return(observedscore) }       #end of function mes
  
  
#source("C:\\Dropbox\\Lab\\zSoftware\\Github\\enumR\\R\\GenStructure.R")
  GenStructure = function(nfactors=5,rfactors=0,r_norm=FALSE,r_norm_sd=.015){
    if(!r_norm){
      Structure =matrix(rep(rfactors,nfactors^2),ncol=nfactors)
    }else{Structure =diag(nfactors)
    Structure[upper.tri(Structure, diag=FALSE)] <- rnorm(mean=rfactors,sd=r_norm_sd,n=(nfactors^2-nfactors)/2)
    Structure <- round(Structure + t(Structure) - diag(diag(Structure)) ,digits=2)}
    return(Structure)
  }
#source("C:\\Dropbox\\Lab\\zSoftware\\Github\\enumR\\R\\GenFactorMatrix.R") 
  GenFactorMatrix = function(nfactors = 5, items= c(5,5,5,5,5), itemsR=c(2,2,2,2,2), loading=.5,loading_norm=FALSE,loading_norm_sd=.025){ # Generate Factor Loadings
    
    Loadings = matrix(rep(0,sum(items)*nfactors),
                      ncol=nfactors,byrow=TRUE)
    if(!loading_norm){
      for(i in 1:nfactors){
        Loadings[(sum(items[1:i-1])+1):sum(items[1:i]),i]<-rep(loading,items[i])
        
        if(itemsR[i]>0){ #reverse score items
          Loadings[(sum(items[1:i-1])+1):(sum(items[1:i-1])+itemsR[i]),i]<- -rep(loading,itemsR[i])}
      }
    } else{
      for(i in 1:nfactors){
        Loadings[(sum(items[1:i-1])+1):sum(items[1:i]),i]<-round(rnorm(items[i], mean = loading,sd=loading_norm_sd),digits=2)
        
        if(itemsR[i]>0){ #reverse score items
          Loadings[(sum(items[1:i-1])+1):(sum(items[1:i-1])+itemsR[i]),i]<- - Loadings[(sum(items[1:i-1])+1):(sum(items[1:i-1])+itemsR[i]),i]
        }
      }
    }
    return(Loadings)
  }
library(beepr)
resultsA=resultsB=NULL
  paramNames <- c("seed",
                  "ndatasets",
                  "patternmatrix",
                  "effectmatrix",
                  "nfactors",
                  "loading",
                  "items",
                  "items_p_f",
                  "itemsR_p_f",
                  "itemsR",
                  "loading_norm",
                  "loading_norm_sd",
                  "rfactors",
                  "r_norm" ,
                  "r_norm_SD",
                  "samplesize",
                  "method",
                  "rotation",
                  "custom_item",
                  "f1_items",
                  "f2_items",
                  "f3_items",
                  "f4_items",
                  "f5_items",
                  "f6_items",
                  "f7_items",
                  "f8_items",
                  "f9_items",
                  "f10_items",
                  "f1_itemsR",
                  "f2_itemsR",
                  "f3_itemsR",
                  "f4_itemsR",
                  "f5_itemsR",
                  "f6_itemsR",
                  "f7_itemsR",
                  "f8_itemsR",
                  "f9_itemsR",
                  "f10_itemsR")

  simulate_enumR <- function(seed=12345,
                             ndatasets=2,
                             patternmatrix=NULL,
                             effectmatrix=NULL,
                             nfactors =5,
                             loading=.5,
                             items=NULL,
                             items_p_f=5,
                             itemsR_p_f=2,
                             itemsR=NULL,
                             loading_norm=FALSE,
                             loading_norm_sd=.05,
                             # effects matrix
                             rfactors=0,
                             r_norm =FALSE,
                             r_norm_SD =.015,
                             samplesize=300,
                             method="mle",
                             rotation="oblimin",
                             custom_item=FALSE,
                             f1_items=NULL,
                             f2_items = NULL,
                             f3_items = NULL,
                             f4_items = NULL,
                             f5_items = NULL,
                             f6_items = NULL,
                             f7_items = NULL,
                             f8_items = NULL,
                             f9_items = NULL,
                             f10_items = NULL,
                             f1_itemsR= NULL,
                             f2_itemsR = NULL,
                             f3_itemsR = NULL,
                             f4_itemsR = NULL,
                             f5_itemsR = NULL,
                             f6_itemsR = NULL,
                             f7_itemsR = NULL,
                             f8_itemsR = NULL,
                             f9_itemsR = NULL,
                             f10_itemsR = NULL
  ){
  
    ###prep
    set.seed(seed)
    require(psych)
    require(GPArotation)
    require(GAIPE);require(audio)
    
    simulation<-ndatasets
    
    # cleanup items
    if(is.null(items)){
      if(custom_item){
        items=c(f1_items,
                f2_items,
                f3_items,
                f4_items,
                f5_items,
                f6_items,
                f7_items,
                f8_items,
                f9_items,
                f10_items)
        
        itemsR=c(f1_itemsR,
                 f2_itemsR,
                 f3_itemsR,
                 f4_itemsR,
                 f5_itemsR,
                 f6_itemsR,
                 f7_itemsR,
                 f8_itemsR,
                 f9_itemsR,
                 f10_itemsR)
        
      }else{
        items=rep(items_p_f,nfactors)
        itemsR=rep(itemsR_p_f,nfactors)
      }
      }

    
    ## Results
meganames<-c( "dof","chisq","prob","sqresid","fit", "RMSEA","RMSEA_lower","RMSEA_upper","BIC","SABIC","null.model","null.dof","null.chisq","complex","eChisq","SRMR","eCRMS","eBIC","TLI","fa.fit","fa.fit.off","objective","ENull","eSABIC","null.chisq1" , "TLI.m","NFI.m","CFI.m","RMSEA_lower.m","RMSEA_upper.m","eRMS","cfit.1","cfit.2","cfit.3","cfit.4","cfit.5","cfit.6","cfit.7","cfit.8","cfit.9","cresidual.1","cresidual.2","cresidual.3","cresidual.4","cresidual.5","cresidual.6","cresidual.7","cresidual.8","cresidual.9","map.values","rfactors","samplesize","loadings","dataset","error","factor")
    mega<-data.frame(matrix(rep(NA,length(meganames)),nrow=1))
    names(mega)=meganames
    mega_null<-mega
    mega_null$error<- 1
    
 #   percentiles<-data.frame(c(1))
 #   names(percentiles)<-"ds"
    
    
    if(is.null(patternmatrix)){
      fmodel=GenFactorMatrix(nfactors=nfactors, items= items[1:nfactors], itemsR=itemsR[1:nfactors], loading=loading,loading_norm=loading_norm,loading_norm_sd=loading_norm_sd)
    }else{fmodel=patternmatrix}
    if(is.null(effectmatrix)){
      smodel=GenStructure(nfactors=nfactors,rfactors=rfactors,r_norm=r_norm,r_norm_sd=r_norm_sd)
    }else{smodel=effectmatrix}
    ## Gen Data
    for(d in 1:simulation) {
      data=GenData(fmodel=fmodel,effect=smodel,samplesize=samplesize)
      cov<-cov(data)
      test<-	try(results<- enumR(cov,n=nfactors,samplesize=samplesize,covar=TRUE,truemodel = TRUE,method = method), silent=FALSE)
      if ('try-error' %in% class(test)) {mega[d,]<- mega_null
      } else {mega[d,]<- results
      mega$factor[d]<-nfactors
      mega$error[d]<- 0
      }
      mega$dataset[d]<-d
      
      print(d)
    }
    mega$rfactors<-rfactors
    mega$samplesize<-samplesize
    mega$loadings<-loading
    
    ds_true<-mega[mega$error==0,]
 
    
    return(ds_true)
  }
  percentile_nav <- function(ds_true) {
  #  wait(length(ds_true)>0,timeout=NA)
    percentiles<-data.frame(c(1))
    names(percentiles)<-"RMSEA_.95"
    percentiles$RMSEA_.95<-quantile(ds_true$RMSEA,.95,na.rm=TRUE)
    percentiles$MAP_.95<-quantile(ds_true$map.values,.95,na.rm=TRUE)
    #percentiles$chisq_.95<-quantile(ds_true$chisq,.95,na.rm=TRUE)
    #percentiles$BIC_.95<-quantile(ds_true$BIC,.95,na.rm=TRUE)
    percentiles$SRMR_.95<-quantile(ds_true$SRMR,.95,na.rm=TRUE)
    #percentiles$SABIC_.95<-quantile(ds_true$SABIC,.95,na.rm=TRUE)
    #percentiles$TLI_.95<-quantile(ds_true$TLI,.95,na.rm=TRUE)
   # percentiles$chisq_sig_.95<-quantile(ds_true$prob,.95,na.rm=TRUE)
    #percentiles$RMSEA_.05<-quantile(ds_true$RMSEA,.05,na.rm=TRUE)
    #percentiles$MAP_.05<-quantile(ds_true$map.values,.05,na.rm=TRUE)
    #percentiles$chisq_.05<-quantile(ds_true$chisq,.05,na.rm=TRUE)
    #percentiles$BIC_.05<-quantile(ds_true$BIC,.05,na.rm=TRUE)
    #percentiles$SRMR_.05<-quantile(ds_true$SRMR,.05,na.rm=TRUE)
    #percentiles$SABIC_.05<-quantile(ds_true$SABIC,.05,na.rm=TRUE)
    percentiles$TLI_.05<-quantile(ds_true$TLI,.05,na.rm=TRUE)
    percentiles$chisq_sig_.05<-quantile(ds_true$prob,.05,na.rm=TRUE)
   # print(percentiles)
    beep()
   return(percentiles)

  }

# Define server logic required to generate and plot a random distribution
#
# Idea and original code by Pierre Chretien
# Small updates by Michael Kapler
#
function(input, output, session) {

	getParams <- function(prefix) {
		input[[paste0(prefix, "_recalc")]]

		params <- lapply(paramNames, function(p) {
			input[[paste0(prefix, "_", p)]]
		})
		names(params) <- paramNames
		params
	}

  # Function that generates scenarios and computes NAV. The expression
  # is wrapped in a call to reactive to indicate that:
  #
  #  1) It is "reactive" and therefore should be automatically
  #     re-executed when inputs change
  #
  resultsA <- eventReactive(input$a_recalc, {do.call(simulate_enumR, getParams("a"))})
  resultsB <- reactive(do.call(simulate_enumR, getParams("b")))

  # Expression that plot NAV paths. The expression
  # is wrapped in a call to renderPlot to indicate that:
  #
  #  1) It is "reactive" and therefore should be automatically
  #     re-executed when inputs change
  #  2) Its output type is a plot
  #
 # output$a_print <- renderPrint({
#  percentile_nav(resultsA())
 # })
  
  output$a_table= renderTable({ head(percentile_nav(resultsA()))},  
              spacing = 'xs',  
              rownames = FALSE,  
              digits = 3)  
 # output$b_table= renderTable({ head(percentile_nav(resultsB()))},  
                           #   spacing = 'xs',  
                          #    rownames = TRUE,  
                          #    digits = 3)  

}
