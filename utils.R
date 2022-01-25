###################
# Utils functions #
###################
library(Formula)
library(LMest)

# Reference:
# - Bartolucci, F., Pandolfi, S., and Pennoni, F. (2017). LMest: An R Package for Latent Markov Models 
#   for Longitudinal Categorical Data. Journal of Statistical Software, 81(4):1-38.


# Function lmestDecoding.Vlogit() is a modified version of function lmestDecoding() contained
# in LMest package that given
#   - the estimated LM model with multinomial logit parametrization (est)
#   - the subject-unit for the decoding (sequence)
# returns:
#   - path prediction Ul obtained by local decoding
#   - path prediction Ug obtained by local decoding
#   - matrix V of posterior probabilities in Equation (4).
lmestDecoding.Vlogit <- function(est, sequence = NULL)
{
  newdata = NULL
  formula = NULL
  index = NULL
  if(is.null(newdata))
  {
    newdata <- est$data
    # id <- attributes(est)$id
    # tv <- attributes(est)$time
    # tv.which <- attributes(est)$whichtv
    # id.which <- attributes(est)$whichid
    id <- newdata$patid
    tv <- newdata$cycno
    tv.which <- 2
    id.which <- 1
    data.new <- newdata[,-c(tv.which,id.which)]
    
    if(is.null(formula))
    {
      formula = attributes(est)$responsesFormula
    }
    #formula = attributes(est)$latentFormula
  }else{
    if(is.null(index))
    {
      id <- attributes(est)$id
      tv <- attributes(est)$time
      tv.which <- attributes(est)$whichtv
      id.which <- attributes(est)$whichid
      data.new <- newdata
    }else{
      id.which <- which(names(newdata) == index[1])
      tv.which <- which(names(newdata) == index[2])
      
      if(is.null(index))
      {
        stop("id and time must be provided")
      }
      if(length(index) !=2)
      {
        stop("id and time must be provided")
      }
      
      if(length(id.which) == 0)
      {
        stop("the id column does not exist")
      }
      
      if(length(tv.which) == 0)
      {
        stop("the time column does not exist")
      }
      id <- newdata[,id.which]
      tv <- newdata[,tv.which]
      data.new <- newdata[,-c(tv.which,id.which)]
    }
  }
  
  if(is.character(id) | is.factor(id))
  {
    warning("conversion of id colum in numeric. The id column must be numeric")
    id <- as.numeric(id)
  }
  if(is.character(tv) | is.factor(tv))
  {
    warning("conversion of time column in numeric. The time column must be numeric")
    tv <- as.numeric(tv)
  }
  
  
  data.new <- newdata[,-c(tv.which,id.which), drop = FALSE]
  if(is.null(formula))
  {
    Y <- data.new
    Xmanifest <- NULL
    Xinitial <- NULL
    Xtrans <- NULL
  }else{
    temp <-  getResponses(data = data.new,formula = formula)
    Y <- temp$Y
    Xmanifest <- temp$X
    Xinitial <- NULL
    Xtrans <- NULL
  }
  
  # if(!is.null(latentFormula) & !is.null(latentFormula[[2]]))
  # {
  #   temp <- getLatent(data = data,latent = latentFormula, responses = responsesFormula)
  #   Xinitial <- temp$Xinitial
  #   Xtrans <- temp$Xtrans
  # }
  tmp <-  long2matrices.internal(Y = Y, id = id, time = tv, yv = rep(1,max(id)),
                                 Xinitial = Xinitial, Xmanifest = Xmanifest, Xtrans = Xtrans)
  
  #model <- tmp$model
  #Xinitial <- tmp$Xinitial
  #Xmanifest <- tmp$Xmanifest
  #Xtrans <- tmp$Xtrans
  Y <- tmp$Y
  yv <- tmp$freq
  if(min(Y,na.rm=T)>0){
    for(i in 1:dim(Y)[3])
    {
      Y[,,i] <- Y[,,i] - min(Y[,,i],na.rm = TRUE)
    }
  }
  if(!is.null(sequence))
  {
    Y <- Y[sequence,,, drop = FALSE]
    yv <- yv[sequence]
  }
  
  
  
  ## Start Computation
  miss = any(is.na(Y))
  if(miss){
    R = 1 * (!is.na(Y))
    Y[is.na(Y)] = 0
  }else{
    R = NULL
  }
  
  if(dim(est$Psi)[3]==1){
    if(is.vector(Y)){
      Y = t(Y)
      if(miss) R = t(R)
    }
    #Y <- Y[,,,drop = TRUE]
    n = nrow(Y); TT = ncol(Y)
  }else{
    if(is.matrix(Y)){
      Y = array(Y,c(1,dim(Y)))
      if(miss) R = array(R,c(1,dim(R)))
    }
    n = dim(Y)[1]; TT = dim(Y)[2]; r = dim(Y)[3]
  }
  piv = est$Piv[sequence,]; Pi = est$PI[,,sequence,]; Psi = est$Psi
  #piv = est$piv; Pi = est$Pi; Psi = est$Psi
  k = length(piv)
  out =  complk(Y,R,rep(1,n),piv,Pi,Psi,k)
  Phi = out$Phi; L = out$L; pv = out$pv
  V = array(0,c(n,k,TT)) # Matrix of posterior probabilities in Equation (4)
  Yvp = matrix(1/pv,n,k)
  M = matrix(1,n,k)
  V[,,TT] = Yvp*L[,,TT] # prod_{l=1}^{TT} phi_{y^(l)|u} * pi^(1) * prod_{m=2}^{TT} tau^{m} / P(Y=y)
  for(t in seq(TT-1,2,-1)){
    M = (Phi[,,t+1]*M)%*%t(Pi[,,t+1])
    V[,,t] = Yvp*L[,,t]*M
  }
  M = (Phi[,,2]*M)%*%t(Pi[,,2])
  V[,,1] = Yvp*L[,,1]*M
  # Local decoding
  Ul = matrix(0,n,TT)
  for(i in 1:n) for(t in 1:TT) Ul[i,t] = which.max(V[i,,t])
  if(n==1) Ul = as.vector(Ul)
  # Global decoding (Viterbi)
  R = L; Ug = matrix(0,n,TT)
  for(i in 1:n) for(t in 2:TT) for(u in 1:k) R[i,u,t] = Phi[i,u,t]*max(R[i,,t-1]*Pi[,u,t])
  if(n==1) Ug[,TT] = which.max(R[,,TT])
  else Ug[,TT] = apply(R[,,TT],1,which.max)
  for(i in 1:n) for(t in seq(TT-1,1,-1)) Ug[i,t] = which.max(R[i,,t]*Pi[,Ug[i,t+1],t+1])
  if(n==1) Ug = as.vector(Ug)
  
  out = list(Ul=Ul,Ug=Ug)
  return(list(V=V,Ul=Ul,Ug=Ug))
}

# Auxiliary functions
getResponses <- function(data, formula)
{
  #data <- data.frame(data)
  if(is.null(formula))
  {
    Y <- data
    X <- NULL
  }else{
    formula <- Formula(formula)
    ll <- length(formula)
    Y <- model.part(formula, data = model.frame(formula, data = data,na.action = NULL), lhs = 1)
    Y <- data.matrix(Y)
    
    X <- NULL
    if(ll[2] != 0)
    {
      X <- model.matrix(formula, model.frame(formula = formula,data,na.action = NULL))
      X <- data.matrix(X)
    }
  }
  out <- list(Y = Y,
              X = X)
  return(out)
}


long2matrices.internal <- function(Y, id, time, yv = NULL,
                                   Xinitial = NULL, Xmanifest = NULL, Xtrans = NULL, cont = FALSE)
{
  # preliminaries
  idu = unique(id)
  n = length(idu)
  TT = max(time)
  
  XXinitial = NULL
  XXmanifest = NULL
  XXtrans = NULL
  init = FALSE
  manifest = FALSE
  trans = FALSE
  if(!is.null(Xinitial)){
    if(dim(Xinitial)[2] == 0)
    {
      Xinitial <- NULL
    }else{
      Xinitial = as.matrix(Xinitial)
      nxInitial = ncol(Xinitial)
      XXinitial <- matrix(NA,n,nxInitial)
    }
    init = TRUE
  }
  
  if(!is.null(Xtrans)){
    
    if(dim(Xtrans)[2] == 0)
    {
      Xtrans <- NULL
    }else{
      Xtrans = as.matrix(Xtrans)
      nxTrans = ncol(Xtrans)
      XXtrans <- array(NA, c(n,TT,nxTrans))
    }
    trans = TRUE
  }
  
  if(!is.null(Xmanifest)){
    Xmanifest = as.matrix(Xmanifest)
    nxMan = ncol(Xmanifest)
    XXmanifest = array(NA,c(n,TT,nxMan))
    manifest = TRUE
  }
  
  if(isTRUE(init) | isTRUE(trans))
  {
    if(isTRUE(manifest))
    {
      model <- "LMlatentManifest"
      stop("covariates on both Latent and Manifest are not allowed",call. = FALSE)
    }else{
      model <- "LMlatent"
      if(cont)
      {
        model = "LMlatentcont"
      }
    }
    
  }else if (isTRUE(manifest))
  {
    model <- "LMmanifest"
    if(ncol(Y) > 1){
      warning("multivariate data are not allowed; only the first response variable is considered", call. = FALSE)
    }
  }else{
    
    model <- "LMbasic"
    if(isTRUE(cont))
    {
      model = "LMbasiccont"
    }
  }
  Y = as.matrix(Y)
  ny = ncol(Y)
  # create matrices
  freq <- NULL
  if(model == "LMbasic" | model == "LMbasiccont")
  {
    if(is.null(yv) && !cont)
    {
      temp <- aggr_data_long(data = Y, id = id, time = time, NAs = 999)
      freq = temp$freq
      id <- temp$Y[,1]
      time <- temp$Y[,2]
      Y = as.matrix(temp$Y[,-c(1,2)])
    }else{
      freq = yv
      id <- id
      time <- time
    }
    
    idu = unique(id)
    n = length(idu)
    
    YY = array(NA,c(n,TT,ny))
    for(i in 1:n){
      ind = which(id==idu[i])
      
      tmp = 0
      for(t in time[ind]){
        tmp=tmp+1
        YY[i,t,] = Y[ind[tmp],]
      }
    }
  }else if(model == "LMlatent" | model == "LMlatentcont")
  {
    YY = array(NA,c(n,TT,ny))
    for(i in 1:n){
      ind = which(id==idu[i])
      timeid <- time[ind]
      if(!is.null(Xinitial))
      {
        timeid1 <- ind[timeid==1]
        if(!length(timeid1)==0)
        {
          XXinitial[i,] = Xinitial[timeid1,]
        }
        
      }
      tmp = 0
      for(t in timeid){
        tmp=tmp+1
        indTemp <- ind[tmp]
        if(!length(indTemp)==0)
        {
          
          if(!is.null(Xtrans))
          {
            XXtrans[i,t,] = Xtrans[indTemp,]
          }
          YY[i,t,] = Y[indTemp,]
        }
      }
    }
    XXtrans <- XXtrans[,-1,, drop = FALSE]
    freq = rep(1,nrow(YY))
  }else if(model == "LMmanifest")
  {
    YY = array(NA,c(n,TT,ny))
    for(i in 1:n){
      ind = which(id==idu[i])
      
      tmp = 0
      for(t in time[ind]){
        tmp=tmp+1
        if(!is.null(Xmanifest))
        {
          XXmanifest[i,t,] = Xmanifest[ind[tmp],]
        }
        YY[i,t,] = Y[ind[tmp],]
      }
    }
    freq = rep(1,nrow(YY))
    
  }
  
  
  # output
  out = list(Y = YY,
             Xinitial = XXinitial,
             Xmanifest = XXmanifest,
             Xtrans = XXtrans,
             model = model,
             freq = freq)
  return(out)
  
}


complk <- function(S,R,yv,piv,Pi,Psi,k){
  # Preliminaries
  sS = dim(S)
  ns = sS[1]
  TT = sS[2]
  if(length(sS)==2) r = 1 else r = sS[3]
  if(r==1){
    if(is.matrix(S)) S = array(S,c(dim(S),1))
    if(is.matrix(R)) R = array(R,c(dim(R),1))
  }
  miss = !is.null(R)
  # Compute log-likelihood
  Phi = array(1,c(ns,k,TT)); L = array(0,c(ns,k,TT))
  if(miss){
    for(j in 1:r) Phi[,,1] = Phi[,,1]*(Psi[S[,1,j]+1,,j]*R[,1,j]+(1-R[,1,j]))
  }else{
    for(j in 1:r) 
      Phi[,,1] = Phi[,,1]*Psi[S[,1,j]+1,,j] 
    # Phi[,,1] = MATRICE n x k
    # prod_{j=1}^{r} P(Y_j^(1)=y_j | U^(1)=u) = phi_{y^(1)|u} u=1,...,k
  }
  # L[,,1] = matrice n x k
  L[,,1] = Phi[,,1]%*%diag(piv) # piv^(1)*phi_{y^(1)|u} u=1,...,k
  for(t in 2:TT){
    if(miss){
      for(j in 1:r) Phi[,,t] = Phi[,,t]*(Psi[S[,t,j]+1,,j]*R[,t,j]+(1-R[,t,j]))
    }else{
      for(j in 1:r) Phi[,,t] = Phi[,,t]*Psi[S[,t,j]+1,,j]
      # Phi[,,t] = MATRICE n x k
      # prod_{j=1}^{r} P(Y_j^(t)=y_j | U^(t)=u) = phi_{y^(t)|u} u=1,...,k
    }
    L[,,t] = Phi[,,t]*(L[,,t-1]%*%Pi[,,t])
    # prod_{l=1}^{t} phi_{y^(l)|u} * pi^(1) * prod_{m=2}^{t} tau^{m}
  }
  if(ns==1) pv = sum(L[1,,TT])
  else pv = rowSums(L[,,TT]) # P(Y=y)
  lk = sum(yv*log(pv))
  out = list(lk=lk,Phi=Phi,L=L,pv=pv)
}


