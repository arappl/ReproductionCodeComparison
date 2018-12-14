# ## A.1 ESTIMATION ACCURACY #################################################
rm(list = ls())
if (!require("parallel")) install.packages("parallel")
if (!require("JM")) install.packages("JM")
if (!require("joineRML")) install.packages("joineRML")
if (!require("assertthat")) install.packages("assertthat")

# path <- "/home/arappl/Dokumente/04_Vergleich_No4/"
# setwd()
path <- getwd()
sub_dir <- "Results"
dir.create(file.path(path, sub_dir))
sub_dir <- "JM_code_snippets"
dir.create(file.path(path, sub_dir))

# ## A.4 FUNCTIONS ###########################################################
simJM <- function(n, n_i, alpha,
                  beta, betas, betals, betatimeind, lambda, betat = 0, noninf = 0, noninfs = 0, noninfls = 0) {
  
  id <- rep(1:n, each = n_i)
  day_in_year <- sample(1:365, n*n_i, replace = TRUE)
  time <- rep(seq(0,(n_i-1)*365, 365), n)
  time <- time + day_in_year
  for(i in 1:n){time[id==i] = time[id==i] - min(time[id==i])}
  time <- time/(n_i*365)
  first = rep(c(1, rep(0,n_i-1)), n)
  rint <- rnorm(n = n, mean =  0, sd = sqrt(2))
  rslope <- rnorm(n = n, mean = 0, sd = sqrt(0.1))
  R_mean <- cbind(rep(rint, each = n_i), rep(rslope, each = n_i))
  eta_ls_mean <- rowSums(cbind(1,time)*R_mean)
  
  if(all(betals == 0) == FALSE){
    if(betatimeind != 0){
      betatime <- betals[betatimeind]
      betals <- betals[-betatimeind]
    }else{betatime <- 0}
    Xls <- matrix(nrow = n*n_i, ncol = length(betals))
    if(ncol(Xls)>0){ ## ~~ 
      for(j in 1:ncol(Xls)){
        sa <- rnorm(n)
        Xls[,j] <- rep(sa, each = n_i)
        Xls[,j] <- (Xls[,j]-mean(Xls[,j]))/sd(Xls[,j])
      }
    }
  }else{Xls =matrix(1,nrow=n*n_i, ncol=1); betatime = 0} 
  eta_ls_mean = eta_ls_mean + Xls%*%betals + betatime*time
  
  X <- matrix(nrow = n*n_i, ncol = length(beta))
  X[,1] <- 1
  if(ncol(X)>1){
    for(j in 2:ncol(X)){
      X[,j] <- rnorm(n*n_i, 0, 1)
      X[,j] <- (X[,j] - mean(X[,j]))/sd(X[,j])
    }}
  eta_l_mean <- X%*%beta + betat*time
  
  
  y <- rnorm(n*n_i, eta_ls_mean + eta_l_mean, sqrt(0.5))
  
  Xs <- matrix(nrow = n, ncol = length(betas))
  for(j in 1:ncol(Xs)){
    Xs[,j] = rnorm(n)
    Xs[,j] <- (Xs[,j] - mean(Xs[,j]))/sd(Xs[,j])
  }
  eta_s_mean = Xs%*%betas
  
  u = runif(n)
  T_surv = (log((-log(1-u)*alpha*(betatime+R_mean[first==1,2]))/(lambda*exp(eta_s_mean))+exp(alpha*(R_mean[first==1,1] + Xls[first==1, , drop=FALSE]%*%betals)))-alpha*(R_mean[first==1,1] + Xls[first==1, , drop=FALSE]%*%betals))/(alpha*(betatime+R_mean[first==1,2]))
  T_surv[is.nan(T_surv)] = 2
  time_mat <- matrix(nrow = n_i, data = time)
  
  delta = rep(1, n)
  for(i in 1:n){
    if(!is.na(T_surv[i])){
      if(T_surv[i]>max(time_mat[,i])){T_surv[i] = max(time_mat[,i]); delta[i] = 0}
      else if(which.max(time_mat[,i]>T_surv[i])<=n_i){time_mat[which.max(time_mat[,i]>T_surv[i]):n_i,i] = 42}
    }
  }
  
  time_zero = which(as.vector(time_mat)==42)
  
  if (length(time_zero == 0) > 0) {
    id <- id[-time_zero]
    y <- y[-time_zero]
    X <- X[-time_zero, ,drop=FALSE]
    Xls <- Xls[-time_zero, , drop=FALSE]
    time <- time[-time_zero]
    R_mean <- R_mean[-time_zero, ]
  }
  
  if(betat!=0){X <- cbind(X, time)}
  
  if (betatimeind != 0) {
    if(is.matrix(Xls)){
      if(betatimeind == ncol(Xls)+1){
        Xls <- cbind(Xls, time)
      }else if(betatimeind == 1){ 
        Xls <- cbind(time,Xls)  
      }else{
        Xls <- cbind(Xls[,c(1:(betatimeind-1))], time, Xls[,c(betatimeind:ncol(Xls))])
      }
    }else{
      if(betatimeind == 1){
        Xls <- cbind(time, Xls)
      }else{
        Xls <- cbind(Xls, time)
      }
    }
  }
  
  if(noninf > 0 | noninfls>0 | noninfs>0){
    for(i in 1:max(c(noninf, noninfls, noninfs))){
      if(i <=noninf){
        X <- cbind(X,rnorm(nrow(X)))}
      if(i <=noninfs){
        Xs <- cbind(Xs,rnorm(nrow(Xs)))}
      if(i<=noninfls){
        Xls <- cbind(Xls, rnorm(length(unique(id)))[id])}
    }
  }
  
  X <- X[,-1, drop = FALSE]
  if(length(betals) != 0 && betals == 0){Xls <- Xls[,-1]}

  return(list(
    #### longitudinal outcome
    "y" = y,
    #### longitudinal predictor fixed effect covariates
    "X" = X,
    #### survival predictor fixed effect covariates
    "Xs" = Xs,
    #### shared predictor fixed effect covariates
    "Xls"=Xls,
    #### shared predictor for survival part
    # "Xls_un" = Xls_un,
    # ### Values for the random effects
    "R_mean"= R_mean,
    #### id indicator
    "id"=id,
    #### measurement times
    "T_long" = time,
    #### Event times
    "T_surv" = T_surv,
    #### Censoring indicator
    "delta" = delta))
}

cvini <- function(x, df, time.effect = time.effect, k = 10, grid1 = NULL, grid2 = NULL, grid3 = NULL){

  if (exists("rel")) {
    if (is.null(grid1)) {grid1 <- c(seq(3, 27, 3), seq(30, 360, 30))}
    if (is.null(grid2)) {grid2 <- c(seq(3, 27, 3), seq(30, 360, 30))}
    if (is.null(grid3)) {grid3 <- c(seq(3, 27, 3), seq(30, 360, 30))}
  } else {
    if (is.null(grid1)) {grid1 <- c(seq(3, 27, 3), seq(30, 360, 30))}
    if (is.null(grid2)) {grid2 <- c(seq(3, 27, 3), seq(30, 360, 30))}
    if (is.null(grid3)) {grid3 <- c(seq(3, 27, 3), seq(30, 360, 30))}
  }
  
  ## 3.1 split data set into k=10 folds --------------------------------------------------------- ##
  # set.seed(1258149)
  set.seed(seedmatrix[x,...])
  m <- length(unique(df$id))/k # size of each subset
  vec <- rep(seq(1:k),each=m)
  testset <- sample(vec, replace = F)
  
  ## 3.2 cross validate each fold via prediction (likelihood) ---------------------------------- ##
  ## inner loop through all k folds of data set (actual cross validation)
  resarr <- indarr <- array(0,c(length(grid1),length(grid2),length(grid3)),
                            dimnames = list(grid1 ,grid2, grid3))
  
  # grid1=ml, grid2=ms, grid3=mls
  
  ## homebrew subset function for lists
  sset <- function(d, filter){
    cand_un <- which(unique(df$id) %in% unique(df$id)[filter])
    cand <- which(df$id %in% unique(df$id)[filter])
    
    if(is.matrix(d)){
      if(dim(d)[1]==length(unique(df$id))){d[cand_un,,drop=FALSE]}else{d[cand,,drop=FALSE]}
    }else{
      if(length(d)==length(unique(df$id))){d[cand_un]}else{d[cand]}
    }
  }
  
  for(kk in seq(1:k)){
    sdf <- lapply(df, sset, testset!=kk)
    kdf <- lapply(df, sset, testset==kk)
    cv.res <- tryCatch({cvres2(data.train = sdf, data.pred = kdf, grid1 = grid1, grid2 = grid2, grid3 = grid3,
                               time.effect = time.effect)},
                       error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
    
    indarr <- indarr + cv.res$indarr
    resarr <- resarr + cv.res$likarr
    
    print(kk)
  }
  
  ## 3.3 combine results and find best set of iterations -------------------------------------- ##
  meanarr <- resarr/(k - indarr) # divide by the number of folds, which worked
  best <- which(meanarr == max(meanarr), arr.ind = T)
  best_iter <- c(grid1[best[1]], grid2[best[2]], grid3[best[3]])
  names(best_iter) <- c("ml","ms","mls")
  return(best_iter)
  return(list("best" = best_iter, "indarr" = indarr))
}

like_cv = function(y, betal, betas, betals, Xl, Xs, Xls,
                   betat, delta, T_surv, alpha, lambda, sigma2, id) {
  first = rep(FALSE, length(id))
  for(i in unique(id)){
    first[which.max(id==i)] = TRUE
  }
  if(is.null(Xl)){
    Xl = as.matrix(rep(1, length(y)))
  }
  etals_m = (as.matrix(Xls)%*%as.matrix(betals))[first==1]
  etals = Xls%*%betals
  etal = cbind(1, Xl)%*%betal
  etas = Xs%*%betas
  
  # if(length(betat) != 0) {
    if (betat != 0) {
      integral = lambda*exp(etas)*((exp(alpha*etals_m + alpha*betat*T_surv) - exp(alpha*etals_m))*(1/(alpha*betat)))
    }else{
      integral = lambda*T_surv*exp(etas + alpha*etals_m)
    }
  # } else {
  #   integral = lambda*T_surv*exp(etas + alpha*etals_m)
  # } 
  
  surv = delta*(log(lambda) + etas + alpha*etals_m + alpha*betat*T_surv) - integral
  long = log(1/sqrt(2*pi*sigma2)) - (y - etals - etal)^2/(2*sigma2)
  like = sum(surv) + sum(long)
  return(like)
}

cvres2 <- function(data.train, data.pred, grid1 = grid1, grid2 = grid2, grid3 = grid3, time.effect = time.effect, alpha = alpha){
  
  likarr <- indarr <- array(0,c(length(grid1), length(grid2), length(grid3)),
                            dimnames = list(grid1 ,grid2, grid3))
  # grid1=ml, grid2=ms, grid3=mls
  
  
  for(ml_akt in grid1){
    for(ms_akt in grid2){
      for(mls_akt in grid3){
        # print(c(ml_akt, ms_akt, mls_akt))
        mod <- tryCatch({JMboost(y = data.train$y, Xl = data.train$X, Xs = data.train$Xs, Xls = data.train$Xls, delta = data.train$delta,
                                 T_long = data.train$T_long, T_surv = data.train$T_surv, id = data.train$id,
                                 mstop_l = ml_akt, mstop_s = ms_akt, mstop_ls = mls_akt, time.effect = time.effect, alpha = alpha, verbose = F)},
                        error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        
        like <-  tryCatch({like_cv(y=data.pred$y,  Xl=data.pred$X, Xs=data.pred$Xs, Xls=data.pred$Xls, delta=data.pred$delta,
                                   T_surv = data.pred$T_surv, id=data.pred$id, time.effect = time.effect, betal=as.matrix(c(mod$int, mod$betal)),
                                   betas=as.matrix(mod$betas), betat = mod$betat, betals=as.matrix(mod$betals), alpha=mod$alpha, lambda=mod$lambda,
                                   sigma2=mod$sigma2)},
                          error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        
        if (is.null(mod)) {
          likarr[as.character(ml_akt), as.character(ms_akt), as.character(mls_akt)] <-  0
          indarr[as.character(ml_akt), as.character(ms_akt), as.character(mls_akt)] <-  1
        } else {
          likarr[as.character(ml_akt), as.character(ms_akt), as.character(mls_akt)] <-  like
        }
        
      }
    }
  }
  return(list("likarr" = likarr, "indarr" = indarr))
}

JMboost = function(y, Xl = NULL, Xs = NULL, Xls = NULL, delta, T_long, T_surv, id, alpha=.001, lambda=.5,
                   mstop_l, mstop_s, mstop_ls, time.effect = TRUE, ny=.1, verbose=FALSE){
  
  sigma2 = var(y)
  mstop = max(mstop_l, mstop_s, mstop_ls)
  n = length(y)
  N = length(unique(id))
  
  ###############################################################
  Xr = matrix(ncol=2*N, nrow = n, data=0)
  Xr1 = matrix(ncol=2*N, nrow = N, data=0)
  unid = order(unique(id))
  id = rep(unid, as.vector(table(id)))
  for(i in 1:N){
    Xr[which(id==as.character(i)),i] = 1
    Xr[which(id==as.character(i)),N+i] = T_long[which(id==as.character(i))]
    Xr1[i,i] = 1; Xr1[i,N+i] = T_surv[i]
  }
  XrA = rbind(Xr[,1:N], Xr1[,1:N])
  lambdaran = mboost:::df2lambda(XrA, 4, weights=1)[2]
  XrAt = t(XrA)
  SA = solve(XrAt%*%XrA + lambdaran*diag(N))%*%XrAt
  XrB = rbind(Xr[,c((N+1):(2*N))], Xr1[,c((N+1):(2*N))])
  lambdaran = mboost:::df2lambda(XrB, 4, weights=1)[2]
  XrBt = t(XrB)
  SB = solve(XrBt%*%XrB + lambdaran*diag(N))%*%XrBt
  gamma0 = rep(0, N)
  gamma1 = rep(0, N)
  ###############################################################
  ### set starting values based on chosen sets of covariates
  int = 0
  betal = 0
  betas = 0
  betat = 0
  betals = 0
  if(is.null(Xl)) {pl = 0; Xl = 0} else {pl = ncol(Xl); betal = rep(0, pl)}
  if(is.null(Xs)) {ps = 0; Xs = 0} else {ps = ncol(Xs); betas = rep(0, ps)}
  if(is.null(Xls)) {
    pls = 0; Xls = 0; Xls_un = 0
  } else {
    pls = ncol(Xls)
    first = rep(FALSE, n)
    for(i in 1:n){
      first[which.max(id==i)] = TRUE
    }
    Xls_un = as.matrix(Xls[first==1,])
    betals = rep(0, pls)
  }
  
  ### define storing matrices/vectors
  GAMMA0 = matrix(0, ncol=mstop, nrow=N)
  GAMMA1 = matrix(0, ncol=mstop, nrow=N)
  BETAL = matrix(0, ncol=mstop, nrow=pl)
  BETAS = matrix(0, ncol=mstop, nrow=ps)
  BETALS = matrix(0, ncol=mstop, nrow=pls)
  BETAT = rep(0, mstop)
  INT = rep(0, mstop)
  ALPHA = rep(0, mstop)
  LAMBDA = rep(0, mstop)
  SIGMA2 = rep(0, mstop)
  
  for(m in 1:mstop){
    ###############################################################
    #### S1 #######################################################
    ###############################################################
    if(m <= mstop_l && pl > 0){
      etal = as.vector(int + Xl%*%betal)
      etas = as.vector(Xs%*%betas)
      etals = as.vector(Xr%*%c(gamma0, gamma1) + as.vector(Xls%*%betals) + T_long*betat)
      etals_un = as.vector(gamma0 + gamma1*T_surv + as.vector(Xls_un%*%betals) + betat*T_surv)
      ###################COMPUTING THE GRADIENT######################
      u = (y - etal - etals)/sigma2
      ##################/COMPUTING THE GRADIENT######################
      fits = matrix(0, 3, pl + as.numeric(time.effect))
      if(pl>0){
        for(i in 1:pl){
          fit = mylm(u, Xl[,i])
          fits[1,i] = fit$int
          fits[2,i] = fit$slp
          fits[3,i] = fit$RSS
        }
      }else if(!time.effect){
        int = int + nyl*mean(u)
      }
      if(time.effect){
        fit = mylm(u, T_long)
        fits[1, pl + 1] = fit$int
        fits[2, pl + 1] = fit$slp
        fits[3, pl + 1] = fit$RSS
      }
      best = which.min(fits[3,])
      ####################COMPUTING STEPLENGTH#######################
      phi = function(v){sum(-dnorm(y, mean=(etal+v*(fits[1,best]+Xl[,best]*fits[2,best])+etals), sd=sqrt(sigma2), log=TRUE))}
      phi = Vectorize(phi)
      steps = 1:10
      y_steps = as.vector(phi(steps))
      ny_l = 0.1*steps[which.min(y_steps)]
      ###################/COMPUTING STEPLENGTH#######################
      if(best==pl+1){
        betat = betat + ny_l*fits[2,best]
        int = int + ny_l*fits[1,best]
      }else{
        betal[best] = betal[best] + ny_l*fits[2,best]
        int = int + ny_l*fits[1,best]
      }
    }
    INT[m] = int
    BETAT[m] = betat
    BETAL[,m] = betal
    
    ###############################################################
    #### S2 #######################################################
    ###############################################################
    if(m<=mstop_s && ps>0){
      etal = as.vector(int + Xl%*%betal)
      etas = as.vector(Xs%*%betas)
      etals = as.vector(Xr%*%c(gamma0, gamma1) + as.vector(Xls%*%betals) + T_long*betat)
      etals_un = as.vector(gamma0 + gamma1*T_surv + as.vector(Xls_un%*%betals) + betat*T_surv)
      ###################COMPUTING THE GRADIENT######################
      if(time.effect){
        if(sum(gamma1)!=0 || betat!=0){
          u = delta - lambda*exp(etas) * (exp(alpha*etals_un) - exp(alpha*(gamma0 + Xls_un%*%betals)))/(alpha*(betat+gamma1))
        }else{
          u = delta - lambda*exp(etas + alpha*etals_un)*T_surv
        }
      }else{
        if(sum(gamma1)!=0){
          u = delta - lambda*exp(etas) * (exp(alpha*etals_un) - exp(alpha*(etals_un - gamma1*T_surv)))/(alpha*gamma1)
        }else{
          u = delta - lambda*exp(etas + alpha*etals_un)*T_surv
        }
      }
      ##################/COMPUTING THE GRADIENT######################
      fits = matrix(0, 3, ps)
      for(i in 1:ps){
        fit = mylm(u, Xs[,i])
        fits[1,i] = fit$int
        fits[2,i] = fit$slp
        fits[3,i] = fit$RSS
      }
      best = which.min(fits[3,])
      bestmod = c(fits[1,best],fits[2,best])
      ####################COMPUTING STEPLENGTH#######################
      if(time.effect){
        if(sum(gamma1)!=0 || betat!=0){
          integral = lambda*(exp(alpha*etals_un) - exp(alpha*(etals_un - (gamma1 + betat)*T_surv)))/(alpha*(betat+gamma1))
        }else{
          integral = lambda*exp(alpha*etals_un)*T_surv
        }
      }else{
        if(sum(gamma1)!=0){
          integral = lambda*(exp(alpha*etals_un) - exp(alpha*(etals_un - gamma1*T_surv)))/(alpha*gamma1)
        }else{
          integral = lambda*exp(alpha*etals_un)*T_surv
        }
      }
      update = Xs[,best]*bestmod[2]
      phi = function(v){-sum(delta*(etas + v*update + alpha*etals_un) - lambda*exp(etas + v*update)*integral)}
      phi = Vectorize(phi)
      steps = 1:10
      y_steps = as.vector(phi(steps))
      ny_s = 0.1*steps[which.min(y_steps)]
      ###################/COMPUTING STEPLENGTH#######################
      betas[best] <- betas[best] + ny_s*fits[2,best]
      lambda = lambda * exp(ny_s*fits[1,best])
    }
    LAMBDA[m] = lambda
    BETAS[,m] = betas
    
    ###############################################################
    #### S3 #######################################################
    ###############################################################
    if(m<=mstop_ls && pls>0){
      etal = as.vector(int + Xl%*%betal)
      etas = as.vector(Xs%*%betas)
      etals = as.vector(Xr%*%c(gamma0, gamma1) + as.vector(Xls%*%betals) + T_long*betat)
      etals_un = as.vector(gamma0 + gamma1*T_surv + as.vector(Xls_un%*%betals) + betat*T_surv)
      ###################COMPUTING THE GRADIENT######################
      u_l = (y - etal - etals)/sigma2
      if(time.effect){
        if(sum(gamma1)!=0 || betat!=0){
          u_s = delta*alpha - lambda*exp(etas)*(exp(alpha*etals_un) - exp(alpha*(gamma0 + Xls_un%*%betals)))/(gamma1 + betat)
        }else{
          u_s = delta*alpha - alpha*lambda*exp(etas)*exp(alpha*etals_un)*T_surv
        }
      }else{
        if(sum(gamma1)!=0){
          #u_s = delta*alpha - lambda*exp(etas)*(exp(alpha*etals_un) - exp(alpha*(etals_un-gamma1*T_surv)))/gamma1
          u_s = delta*alpha - alpha*lambda*exp(etas + alpha*(Xls_un%*%betals + Xr[,1:N]%*%gamma0))*((exp(gamma1*T_surv) - 1)/gamma1)
          #u = delta - lambda*exp(etas) * (exp(alpha*etas) - exp(alpha*(etals_un - gamma1*T_surv)))/(alpha*gamma1)
        }else{
          u_s = delta*alpha - alpha*exp(etas)*lambda*exp(alpha*etals_un)*T_surv
        }
      }
      u = c(u_l,u_s)
      ##################/COMPUTING THE GRADIENT######################
      fits = matrix(0, 3, pls + 1 + as.numeric(time.effect))
      
      # fixed effects
      if (pls > 0){ 
        for(i in 1:pls){
          fit = mylm(u, c(Xls[,i], Xls_un[,i]))
          fits[1, i] = fit$int
          fits[2, i] = fit$slp
          fits[3, i] = fit$RSS
        }
      } 
      
      # random effects 
      rfit = rbind(SA, SB)%*%u
      fits[3, pls + 1] = sum((u-(rbind(Xr, Xr1)%*%rfit))^2)
      
      # fixed time effect
      if(time.effect){
        fit = mylm(u, c(T_long, T_surv))
        fits[1, pls+2] = fit$int
        fits[2, pls+2] = fit$slp
        fits[3, pls+2] = fit$RSS
      }
      
      best = which.min(fits[3,])
      if(m == 1) {best <- pls + 1} #force update random effects in first iteration
      
      ## ################################# ##     
      # Update on random effects
      if(best == pls+1){
        gamma0 = gamma0 + ny*rfit[c(1:length(gamma0))]
        gamma1 = gamma1 + ny*rfit[c((length(gamma0)+1):(2*length(gamma0)))]
        lambda = lambda*exp(ny*alpha*fits[1,best])
      } else if (best == pls+2) { # update on time-effect
        bestmod = c(fits[1, best], fits[2, best])
        # COMPUTING STEPLENGTH
        phi = function(v){
          update = v*bestmod[2]
          updint = v*bestmod[1]
          if(sum(gamma1)!=0 || betat !=0){
            integral = lambda*exp(etas)*(exp(alpha*(etals_un + update*T_surv)) - exp(alpha*(etals_un - (gamma1 + betat)*T_surv)))/(alpha*(betat + update + gamma1))
          }else{
            integral = lambda*exp(etas)*exp(alpha*(etals_un + update*T_surv))*T_surv
          }
          long = -sum((y - etal - updint - etals - v*betat*T_long)/sigma2)
          surv = -sum(delta*(log(lambda) + etas + alpha*(etals_un + betat*T_surv)) - integral)
          return(long + surv)
        }
        phi = Vectorize(phi)
        steps = 1:10
        y_steps = as.vector(phi(steps))
        ny_ls = 0.1*steps[which.min(y_steps)]
        if(sum(is.nan(y_steps))>0){ny_ls = .1}
        
        # UPDATE EFFECT
        betat = betat + ny_ls*fits[2, best]
        int = int + ny_ls*fits[1, best]
        lambda = lambda*exp(ny_ls*alpha*fits[1,best])
      } else { # update on fixed effects
        bestmod = c(fits[1, best], fits[2, best])
        # COMPUTING STEPLENGTH
        phi = function(v){
          update = v*bestmod[2]
          updint = v*bestmod[1]
          if(sum(gamma1)!=0){
            integral = lambda*exp(etas)*exp(alpha*(etals_un + update*Xls_un[,best])) - exp(alpha*(etals_un - gamma1))/(alpha*gamma1)
          }else{
            integral = lambda*exp(etas)*exp(alpha*(etals_un + update*Xls_un[,best]))*T_surv
          }
          long = -sum((y - etal - etals - v*Xls[,best])/sigma2)
          surv = -sum(delta*(log(lambda) + etas + alpha*(etals_un + update*Xls_un[,best])) - integral)
          return(long + surv)
        }
        
        # UPDATE EFFECT
        betals[best] = betals[best] + ny_ls*fits[2,best]
        int = int + ny_ls*fits[1,best]
        lambda = lambda*exp(ny_ls*alpha*fits[1,best])
      } 
    }
    
    INT[m] = int
    BETALS[,m] = betals
    GAMMA0[,m] = gamma0
    GAMMA1[,m] = gamma1
    BETAT[m] = betat
    LAMBDA[m] = lambda
    
    ###############################################################
    #### S4 #######################################################
    ###############################################################
    etal = as.vector(int + Xl%*%betal)
    etas = as.vector(Xs%*%betas)
    etals = as.vector(Xr%*%c(gamma0, gamma1) + as.vector(Xls%*%betals) + T_long*betat)
    etals_un = as.vector(gamma0 + gamma1*T_surv + as.vector(Xls_un%*%betals) + betat*T_surv)
    
    sigma2 = optimize(long.risk, y=y, etal=etal, etals=etals, interval=c(0,100))$minimum
    #sigma2 = sum((y - etal - etals)^2)/(length(y)-pl-pls-1)
    
    # lambda = optimize(surv.risk, alpha=alpha, etals=etals_un, etas=etas, delta=delta, gamma1=gamma1,
    #                   betals=betals, betatimeind=betatimeind, T_surv=T_surv, interval=c(lambda-.1*lambda,lambda+.1*lambda))$minimum
    
    oi.min = min(alpha - 0.1*abs(alpha), alpha - 0.1)
    oi.max = max(alpha + 0.1*abs(alpha), alpha + 0.1)
    optim.int = c(oi.min, oi.max)
    alpha = optimize(surv.risk, lambda=lambda, etals=etals_un, etas=etas, delta=delta, gamma1=gamma1,
                     betals=betals, betat = betat, time.effect = time.effect, T_surv=T_surv, interval=optim.int)$minimum
    
    ALPHA[m] = alpha
    SIGMA2[m] = sigma2
    
    if(verbose){
      print(m)
    }
    
  }
  
  structure(list(GAMMA0 = GAMMA0, GAMMA1 = GAMMA1, BETAL = BETAL, BETAS = BETAS,
                 BETALS = BETALS, BETAT=BETAT, INT=INT, ALPHA = ALPHA, LAMBDA = LAMBDA, SIGMA2 = SIGMA2,
                 gamma0 = gamma0, gamma1 = gamma1, betal = betal, betas = betas, betals = betals,
                 betat = betat, int = int, alpha = alpha, lambda = lambda,  sigma2 = sigma2))
}

long.risk = function(y, etal, etals, sigma2){sum(-dnorm(y, mean = (etal + etals), sd = sqrt(sigma2), log = TRUE))}

surv.risk = function(alpha, lambda, etas, etals_un, T_surv, delta, gamma1, betals, betat, time.effect){
  if(time.effect){
    if(sum(gamma1)!=0 || betat!=0){
      integral = lambda*exp(etas)*(exp(alpha*etals_un) - exp(alpha*(etals_un - (gamma1 + betat)*T_surv)))/(alpha*(betat+gamma1))
    }else{
      integral = lambda*exp(etas)*exp(alpha*etals_un)*T_surv
    }
  }else{
    if(sum(gamma1)!=0){
      integral = lambda*exp(etas)*(exp(alpha*etals_un) - exp(alpha*(etals_un - gamma1*T_surv)))/(alpha*gamma1)
    }else{
      integral = lambda*exp(etas)*exp(alpha*etals_un)*T_surv
    }
  }
  
  risk = -sum(delta*(log(lambda) + etas + alpha*etals_un) - integral)
  return(risk)
}

mylm = function(y,x){
  X = cbind(1, x)
  beta = solve(t(X) %*% X) %*% t(X) %*% y
  RSS = sum((y - X %*% beta)^2)
  return(list("int" = beta[1], "slp" = beta[2], "RSS" = RSS))
}

JMboostc = function(y, Xl = NULL, Xs = NULL, Xls = NULL, delta, T_long, T_surv, id, alpha=.001, lambda=.5,
                    mstop_l, mstop_s, mstop_ls, time.effect = TRUE, ny=.1, eps = 10^(-6), verbose=FALSE){
  
  eps = eps
  
  sigma2 = var(y)
  mstop = max(mstop_l, mstop_s, mstop_ls)
  n = length(y)
  N = length(unique(id))
  
  ###############################################################
  Xr = matrix(ncol=2*N, nrow = n, data=0)
  Xr1 = matrix(ncol=2*N, nrow = N, data=0)
  unid = order(unique(id))
  id = rep(unid, as.vector(table(id)))
  for(i in 1:N){
    Xr[which(id==as.character(i)),i] = 1
    Xr[which(id==as.character(i)),N+i] = T_long[which(id==as.character(i))]
    Xr1[i,i] = 1; Xr1[i,N+i] = T_surv[i]
  }
  XrA = rbind(Xr[,1:N], Xr1[,1:N])
  lambdaran = mboost:::df2lambda(XrA, 4, weights=1)[2]
  XrAt = t(XrA)
  SA = solve(XrAt%*%XrA + lambdaran*diag(N))%*%XrAt
  XrB = rbind(Xr[,c((N+1):(2*N))], Xr1[,c((N+1):(2*N))])
  lambdaran = mboost:::df2lambda(XrB, 4, weights=1)[2]
  XrBt = t(XrB)
  SB = solve(XrBt%*%XrB + lambdaran*diag(N))%*%XrBt
  gamma0 = rep(0, N)
  gamma1 = rep(0, N)
  ###############################################################
  ### set starting values based on chosen sets of covariates
  int = 0
  betal = 0
  betas = 0
  betat = 0
  betals = 0
  if(is.null(Xl)) {pl = 0; Xl = 0} else {pl = ncol(Xl); betal = rep(0, pl)}
  if(is.null(Xs)) {ps = 0; Xs = 0} else {ps = ncol(Xs); betas = rep(0, ps)}
  if(is.null(Xls)) {
    pls = 0; Xls = 0; Xls_un = 0
  } else {
    pls = ncol(Xls)
    first = rep(FALSE, n)
    for(i in 1:n){
      first[which.max(id==i)] = TRUE
    }
    Xls_un = as.matrix(Xls[first==1,])
    betals = rep(0, pls)
  }
  
  ### define storing matrices/vectors
  GAMMA0 = matrix(0, ncol=mstop, nrow=N)
  GAMMA1 = matrix(0, ncol=mstop, nrow=N)
  BETAL = matrix(0, ncol=mstop, nrow=pl)
  BETAS = matrix(0, ncol=mstop, nrow=ps)
  BETALS = matrix(0, ncol=mstop, nrow=pls)
  BETAT = rep(0, mstop)
  INT = rep(0, mstop)
  ALPHA = rep(0, mstop)
  LAMBDA = rep(0, mstop)
  SIGMA2 = rep(0, mstop)
  LIKE <- rep(0, mstop)
  
  m <- 1
  cond <- 1
  
  while(m<=2 || cond > eps){
    
    # for(m in 1:3571){
    ###############################################################
    #### S1 #######################################################
    ###############################################################
    if(pl > 0){
      etal = as.vector(int + Xl%*%betal)
      etas = as.vector(Xs%*%betas)
      etals = as.vector(Xr%*%c(gamma0, gamma1) + as.vector(Xls%*%betals) + T_long*betat)
      etals_un = as.vector(gamma0 + gamma1*T_surv + as.vector(Xls_un%*%betals) + betat*T_surv)
      ###################COMPUTING THE GRADIENT######################
      u = (y - etal - etals)/sigma2
      ##################/COMPUTING THE GRADIENT######################
      fits = matrix(0, 3, pl + as.numeric(time.effect))
      if(pl>0){
        for(i in 1:pl){
          fit = mylm(u, Xl[,i])
          fits[1,i] = fit$int
          fits[2,i] = fit$slp
          fits[3,i] = fit$RSS
        }
      }else if(!time.effect){
        int = int + nyl*mean(u)
      }
      if(time.effect){
        fit = mylm(u, T_long)
        fits[1, pl + 1] = fit$int
        fits[2, pl + 1] = fit$slp
        fits[3, pl + 1] = fit$RSS
      }
      best = which.min(fits[3,])
      ####################COMPUTING STEPLENGTH#######################
      phi = function(v){sum(-dnorm(y, mean=(etal+v*(fits[1,best]+Xl[,best]*fits[2,best])+etals), sd=sqrt(sigma2), log=TRUE))}
      phi = Vectorize(phi)
      steps = 1:10
      y_steps = as.vector(phi(steps))
      ny_l = 0.1*steps[which.min(y_steps)]
      ###################/COMPUTING STEPLENGTH#######################
      if(best==pl+1){
        betat = betat + ny_l*fits[2,best]
        int = int + ny_l*fits[1,best]
      }else{
        betal[best] = betal[best] + ny_l*fits[2,best]
        int = int + ny_l*fits[1,best]
      }
    }
    INT[m] = int
    BETAT[m] = betat
    BETAL[,m] = betal
    
    ###############################################################
    #### S2 #######################################################
    ###############################################################
    if(ps > 0){
      etal = as.vector(int + Xl%*%betal)
      etas = as.vector(Xs%*%betas)
      etals = as.vector(Xr%*%c(gamma0, gamma1) + as.vector(Xls%*%betals) + T_long*betat)
      etals_un = as.vector(gamma0 + gamma1*T_surv + as.vector(Xls_un%*%betals) + betat*T_surv)
      ###################COMPUTING THE GRADIENT######################
      if(time.effect){
        if(sum(gamma1)!=0 || betat!=0){
          u = delta - lambda*exp(etas) * (exp(alpha*etals_un) - exp(alpha*(gamma0 + Xls_un%*%betals)))/(alpha*(betat+gamma1))
        }else{
          u = delta - lambda*exp(etas + alpha*etals_un)*T_surv
        }
      }else{
        if(sum(gamma1)!=0){
          u = delta - lambda*exp(etas) * (exp(alpha*etals_un) - exp(alpha*(etals_un - gamma1*T_surv)))/(alpha*gamma1)
        }else{
          u = delta - lambda*exp(etas + alpha*etals_un)*T_surv
        }
      }
      ##################/COMPUTING THE GRADIENT######################
      fits = matrix(0, 3, ps)
      for(i in 1:ps){
        fit = mylm(u, Xs[,i])
        fits[1,i] = fit$int
        fits[2,i] = fit$slp
        fits[3,i] = fit$RSS
      }
      best = which.min(fits[3,])
      bestmod = c(fits[1,best],fits[2,best])
      ####################COMPUTING STEPLENGTH#######################
      if(time.effect){
        if(sum(gamma1)!=0 || betat!=0){
          integral = lambda*(exp(alpha*etals_un) - exp(alpha*(etals_un - (gamma1 + betat)*T_surv)))/(alpha*(betat+gamma1))
        }else{
          integral = lambda*exp(alpha*etals_un)*T_surv
        }
      }else{
        if(sum(gamma1)!=0){
          integral = lambda*(exp(alpha*etals_un) - exp(alpha*(etals_un - gamma1*T_surv)))/(alpha*gamma1)
        }else{
          integral = lambda*exp(alpha*etals_un)*T_surv
        }
      }
      update = Xs[,best]*bestmod[2]
      phi = function(v){-sum(delta*(etas + v*update + alpha*etals_un) - lambda*exp(etas + v*update)*integral)}
      phi = Vectorize(phi)
      steps = 1:10
      y_steps = as.vector(phi(steps))
      ny_s = 0.1*steps[which.min(y_steps)]
      ###################/COMPUTING STEPLENGTH#######################
      betas[best] <- betas[best] + ny_s*fits[2,best]
      lambda = lambda * exp(ny_s*fits[1,best])
    }
    LAMBDA[m] = lambda
    BETAS[,m] = betas
    
    ###############################################################
    #### S3 #######################################################
    ###############################################################
    if(pls > 0){
      etal = as.vector(int + Xl%*%betal)
      etas = as.vector(Xs%*%betas)
      etals = as.vector(Xr%*%c(gamma0, gamma1) + as.vector(Xls%*%betals) + T_long*betat)
      etals_un = as.vector(gamma0 + gamma1*T_surv + as.vector(Xls_un%*%betals) + betat*T_surv)
      ###################COMPUTING THE GRADIENT######################
      u_l = (y - etal - etals)/sigma2
      if(time.effect){
        if(sum(gamma1)!=0 || betat!=0){
          u_s = delta*alpha - lambda*exp(etas)*(exp(alpha*etals_un) - exp(alpha*(gamma0 + Xls_un%*%betals)))/(gamma1 + betat)
        }else{
          u_s = delta*alpha - alpha*lambda*exp(etas)*exp(alpha*etals_un)*T_surv
        }
      }else{
        if(sum(gamma1)!=0){
          #u_s = delta*alpha - lambda*exp(etas)*(exp(alpha*etals_un) - exp(alpha*(etals_un-gamma1*T_surv)))/gamma1
          u_s = delta*alpha - alpha*lambda*exp(etas + alpha*(Xls_un%*%betals + Xr[,1:N]%*%gamma0))*((exp(gamma1*T_surv) - 1)/gamma1)
          #u = delta - lambda*exp(etas) * (exp(alpha*etas) - exp(alpha*(etals_un - gamma1*T_surv)))/(alpha*gamma1)
        }else{
          u_s = delta*alpha - alpha*exp(etas)*lambda*exp(alpha*etals_un)*T_surv
        }
      }
      u = c(u_l,u_s)
      ##################/COMPUTING THE GRADIENT######################
      fits = matrix(0, 3, pls + 1 + as.numeric(time.effect))
      
      # fixed effects
      if (pls > 0){ 
        for(i in 1:pls){
          fit = mylm(u, c(Xls[,i], Xls_un[,i]))
          fits[1, i] = fit$int
          fits[2, i] = fit$slp
          fits[3, i] = fit$RSS
        }
      } 
      
      # random effects 
      rfit = rbind(SA, SB)%*%u
      fits[3, pls + 1] = sum((u-(rbind(Xr, Xr1)%*%rfit))^2)
      
      # fixed time effect
      if(time.effect){
        fit = mylm(u, c(T_long, T_surv))
        fits[1, pls+2] = fit$int
        fits[2, pls+2] = fit$slp
        fits[3, pls+2] = fit$RSS
      }
      
      best = which.min(fits[3,])
      if(m == 1) {best <- pls + 1} #force update random effects in first iteration
      
      ## ################################# ##     
      # Update on random effects
      if(best == pls+1){
        gamma0 = gamma0 + ny*rfit[c(1:length(gamma0))]
        gamma1 = gamma1 + ny*rfit[c((length(gamma0)+1):(2*length(gamma0)))]
        lambda = lambda*exp(ny*alpha*fits[1,best])
      } else if (best == pls+2) { # update on time-effect
        bestmod = c(fits[1, best], fits[2, best])
        # COMPUTING STEPLENGTH
        phi = function(v){
          update = v*bestmod[2]
          updint = v*bestmod[1]
          if(sum(gamma1)!=0 || betat !=0){
            integral = lambda*exp(etas)*(exp(alpha*(etals_un + update*T_surv)) - exp(alpha*(etals_un - (gamma1 + betat)*T_surv)))/(alpha*(betat + update + gamma1))
          }else{
            integral = lambda*exp(etas)*exp(alpha*(etals_un + update*T_surv))*T_surv
          }
          long = -sum((y - etal - updint - etals - v*betat*T_long)/sigma2)
          surv = -sum(delta*(log(lambda) + etas + alpha*(etals_un + betat*T_surv)) - integral)
          return(long + surv)
        }
        phi = Vectorize(phi)
        steps = 1:10
        y_steps = as.vector(phi(steps))
        ny_ls = 0.1*steps[which.min(y_steps)]
        if(sum(is.nan(y_steps))>0){ny_ls = .1}
        
        # UPDATE EFFECT
        betat = betat + ny_ls*fits[2, best]
        int = int + ny_ls*fits[1, best]
        lambda = lambda*exp(ny_ls*alpha*fits[1,best])
      } else { # update on fixed effects
        bestmod = c(fits[1, best], fits[2, best])
        # COMPUTING STEPLENGTH
        phi = function(v){
          update = v*bestmod[2]
          updint = v*bestmod[1]
          if(sum(gamma1)!=0){
            integral = lambda*exp(etas)*exp(alpha*(etals_un + update*Xls_un[,best])) - exp(alpha*(etals_un - gamma1))/(alpha*gamma1)
          }else{
            integral = lambda*exp(etas)*exp(alpha*(etals_un + update*Xls_un[,best]))*T_surv
          }
          long = -sum((y - etal - etals - v*Xls[,best])/sigma2)
          surv = -sum(delta*(log(lambda) + etas + alpha*(etals_un + update*Xls_un[,best])) - integral)
          return(long + surv)
        }
        
        # UPDATE EFFECT
        betals[best] = betals[best] + ny_ls*fits[2,best]
        int = int + ny_ls*fits[1,best]
        lambda = lambda*exp(ny_ls*alpha*fits[1,best])
      } 
    }
    
    INT[m] = int
    BETALS[,m] = betals
    GAMMA0[,m] = gamma0
    GAMMA1[,m] = gamma1
    BETAT[m] = betat
    LAMBDA[m] = lambda
    
    ###############################################################
    #### S4 #######################################################
    ###############################################################
    etal = as.vector(int + Xl%*%betal)
    etas = as.vector(Xs%*%betas)
    etals = as.vector(Xr%*%c(gamma0, gamma1) + as.vector(Xls%*%betals) + T_long*betat)
    etals_un = as.vector(gamma0 + gamma1*T_surv + as.vector(Xls_un%*%betals) + betat*T_surv)
    
    sigma2 = optimize(long.risk, y=y, etal=etal, etals=etals, interval=c(0,100))$minimum
    #sigma2 = sum((y - etal - etals)^2)/(length(y)-pl-pls-1)
    
    # lambda = optimize(surv.risk, alpha=alpha, etals=etals_un, etas=etas, delta=delta, gamma1=gamma1,
    #                   betals=betals, betatimeind=betatimeind, T_surv=T_surv, interval=c(lambda-.1*lambda,lambda+.1*lambda))$minimum
    
    oi.min = min(alpha - 0.1*abs(alpha), alpha - 0.1)
    oi.max = max(alpha + 0.1*abs(alpha), alpha + 0.1)
    optim.int = c(oi.min, oi.max)
    alpha = optimize(surv.risk, lambda=lambda, etals=etals_un, etas=etas, delta=delta, gamma1=gamma1,
                     betals=betals, betat = betat, time.effect = time.effect, T_surv=T_surv, interval=optim.int)$minimum
    
    ALPHA[m] = alpha
    SIGMA2[m] = sigma2
    
    if(betat != 0){
      integral = lambda*exp(etas)*((exp(alpha*etals_un + alpha*betat*T_surv) - exp(alpha*etals_un))*(1/(alpha*betat)))
    }else{
      integral = lambda*T_surv*exp(etas + alpha*etals_un)
    }
    surv = delta*(log(lambda) + etas + alpha*etals_un + alpha*betat*T_surv) - integral
    long = log(1/sqrt(2*pi*sigma2)) - (y - etals - etal)^2/(2*sigma2)
    like = sum(surv) + sum(long)
    
    LIKE[m] <- like
    
    m <- m + 1
    cond <- abs((LIKE[m - 1] - LIKE[m - 2]) / LIKE[m - 1])
    
    # if (m == mstop || m%%1000 == 0) {
    #   intermediate <- list(GAMMA0 = GAMMA0, GAMMA1 = GAMMA1, BETAL = BETAL, BETAS = BETAS,
    #                        BETALS = BETALS, BETAT=BETAT, INT=INT, ALPHA = ALPHA, LAMBDA = LAMBDA, SIGMA2 = SIGMA2,
    #                        gamma0 = gamma0, gamma1 = gamma1, betal = betal, betas = betas, betals = betals,
    #                        betat=betat, int=int, alpha = alpha, lambda = lambda,  sigma2 = sigma2)
    #   save(intermediate, file = paste0(path, "01_Datenbeispiel/03_Converge_Zwischenspeicher/Dat_conv_m", m, ".RData"))
    #   rm("intermediate")
    # }
    
    if(verbose){
      # if(m%%1000 == 0){print(m)}
      print(m)
    }
    
  }
  conv.iter <- max(which(INT != 0)) + 1
  
  structure(list(GAMMA0 = GAMMA0, GAMMA1 = GAMMA1, BETAL = BETAL, BETAS = BETAS,
                 BETALS = BETALS, BETAT=BETAT, INT=INT, ALPHA = ALPHA, LAMBDA = LAMBDA, SIGMA2 = SIGMA2,
                 gamma0 = gamma0, gamma1 = gamma1, betal = betal, betas = betas, betals = betals,
                 betat=betat, int=int, alpha = alpha, lambda = lambda,  sigma2 = sigma2, 
                 conv.iter = conv.iter))
}

estFun <- function(x, n, n_i, alpha, beta, betas, betals, betatimeind, lambda, betat = 0, noninf = 0, noninfs = 0, noninfls = 0, method = "boost",
                   k = 10, grid1 = NULL, grid2 = NULL, grid3 = NULL) { # method = c("boost", "converge") 
  repeat{
    if ( (noninf + noninfs) > 100 ) {
      set.seed(541689 + x)
    } else {
    # simulate data with specific seed, if fits fail, change the seed (i <-  i + 1)
      i <- i + 1
      seed <- 8469 + x*1000 + i
      set.seed(seed)
    }
    dat <- simJM(n = n, n_i = n_i, alpha = alpha, beta = beta, betas = betas, betals = betals, betatimeind = betatimeind, lambda = lambda, 
                 betat = betat, noninf = noninf, noninfs = noninfs, noninfls = noninfls)
    
    df <- with(dat, data.frame(id, y, X, Xls))
    df.id <- with(dat, data.frame("id" = unique(id), T_surv, delta, Xs))
    
    # +++++ JM fit +++++ # 
    if ( (noninf + noninfs) > 100 ) {
      fitJM <- "too many variables"
    } else {
      fitLME <- tryCatch({ do.call("lme", list(fixed = as.formula(paste0("y~", paste(names(df)[3 : ncol(df)], collapse = "+"))), 
                                    data = df, random = ~ time | id, na.action = na.exclude, control = lmeControl(opt='optim'))) },
                         error = function(e){cat("ERROR :", conditionMessage(e), "/n")})
      fitSURV <- coxph(Surv(time = T_surv, event=delta) ~ . -id, data = df.id, x = TRUE, na.action = na.exclude)
      fitJM <-  tryCatch({ jointModel(fitLME, fitSURV, timeVar = "time", method = "piecewise-PH-aGH") }, 
                         error = function(e){cat("ERROR :",conditionMessage(e), "/n")})
    }
    # +++++ joineRML fit +++++ # 
    # if fitJM returnd no error, fit joineRML
    if (!is.error(fitJM) && !is.null(fitJM)) {
      fitJRML <- tryCatch({mjoint(formLongFixed = list(as.formula(paste0("y~", paste(names(df)[3 : ncol(df)], collapse="+")))),
                                  formLongRandom = list(~ time | id),
                                  formSurv = Surv(time = T_surv, event = delta) ~ . -id,
                                  data = list(df),
                                  survData = df.id, # bei merge(df, df.id) ergeben sich doppelte Variablennamen, muessten sonst extra im code beruecksichtigt werden
                                  timeVar = "time")},
                          error=function(e){cat("ERROR :", conditionMessage(e), "/n")})
      # if fitJRML is also no error, break the repeat loop, else increase i and try again
      if (!is.error(fitJRML) && !is.null(fitJRML)) { break }
    }
  }

  # +++++ JMboost 'fit' part +++++ #
  if (betatimeind == 0) {
    time.effect <- FALSE
  } else {
    time.effect <- TRUE
  }

  if (method == "converge") {
    # +++++ JMboost 'fit' +++++ #
    fitJMb <- tryCatch({JMboostc(y = dat$y, Xl = dat$X, Xs = dat$Xs, Xls = dat$Xls, time.effect = time.effect,
                                 delta = dat$delta, T_long = dat$T_long, T_surv = dat$T_surv, id = dat$id,
                                 mstop_l = 100000, mstop_ls = 100000, mstop_s = 100000, verbose = FALSE)},
                       error=function(e){cat("ERROR :", conditionMessage(e), "/n")})

    return(list(fitJM, fitJRML, fitJMb, seed))

  } else {
    # +++++ JMboost Cross-validation (cv) +++++ #
    # cv_result <- cvini(x = x, df = dat, time.effect = time.effect)
    # best_iter <- cv_result[[1]]
    best_iter <- cvini(x = x, df = dat, time.effect = time.effect, k = k)

    # +++++ JMboost 'fit' +++++ #
    fitJMb <- tryCatch({JMboost(y = dat$y, Xl = dat$X, Xs = dat$Xs, Xls = dat$Xls, time.effect = time.effect,
                                delta = dat$delta, T_long = dat$T_long, T_surv = dat$T_surv, id = dat$id,
                                mstop_l = best_iter[1], mstop_s = best_iter[2], mstop_ls = best_iter[3], verbose = FALSE)},
                       error=function(e){cat("ERROR :", conditionMessage(e), "/n")})

    return(list(fitJM, fitJRML, fitJMb, best_iter, seed))
  }
}

failFun <- function(z, n, n_i, alpha,
                    beta, betas, betals, betatimeind, lambda, betat = 0, noninf = 0, noninfs = 0, noninfls = 0){
  set.seed(541689 + z)
  dat <- simJM(n = n, n_i = n_i, alpha = alpha, beta = beta, betas = betas, betals = betals, betatimeind = betatimeind, lambda = lambda, 
               betat = betat, noninf = noninf, noninfs = noninfs, noninfls = noninfls)
  noninfn <- max(noninf, noninfs, noninfls)
  
  Z <- lapply(1:noninfn, function(nvar, data = dat) {
    df <- with(data, data.frame(id, y, X[, 1:(5 + nvar)], Xls[, 1:(6 + nvar)]))
    df.id <- with(data, data.frame("id" = unique(id), T_surv, delta, Xs[, 1:nvar]))
    
    # +++++ JM fit +++++ # 
    gc(TRUE)
    fitLME <- tryCatch({ do.call("lme", list(fixed = as.formula(paste0("y~", paste(names(df)[3 : ncol(df)], collapse = "+"))), 
                                             data = df, random = ~ time | id, na.action = na.exclude, control = lmeControl(opt='optim'))) },
                       error = function(e){cat("ERROR :", conditionMessage(e), "/n")})
    fitSURV <- coxph(Surv(time = T_surv, event=delta) ~ . -id, data = df.id, x = TRUE, na.action = na.exclude)
    fitJM <-  tryCatch({ jointModel(fitLME, fitSURV, timeVar = "time", method = "piecewise-PH-aGH") }, 
                       error = function(e){cat("ERROR :",conditionMessage(e), "/n")})
    
    fitJRML <- tryCatch({ mjoint(formLongFixed = list(as.formula(paste0("y~", paste(names(df)[3:ncol(df)], collapse = "+")))),
                                 formLongRandom = list(~ time | id),
                                 formSurv = Surv(time = T_surv, event = delta)~ . -id,
                                 data = list(df),
                                 survData = df.id,
                                 timeVar = "time") }, 
                        error=function(e){cat("ERROR :",conditionMessage(e), "/n")})
    
    return(c( ifelse(is.null(fitJM), 0, 1), 
              ifelse(is.null(fitJRML), 0, 1)))
  }
  )
  Z <- matrix(unlist(Z), ncol = length(Z))
  rownames(Z) <- c("JM", "joineRML")
  return(Z)
}


## ACCURACY ############################

n = 100
n_i = 5
alpha = 0.1
lambda = 0.6

noninf <- noninfs <- noninfls <- 0

# w/o fixed time effect
betas = 0.1
beta = c(1.5,-0.5, 0.7, 1.3, 0.3, 0.5)
betat = 0.4
betals = c(0.9, 0.3, -1, 0.2, -0.4)
betatimeind = 0

i <- 0

AM1 <- lapply(1:100, estFun, method = "converge", n = n, n_i = n_i, alpha = alpha, beta = beta[1], betas = betas, betals = 0, betatimeind = betatimeind, lambda = lambda,
                        betat = betat, noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(AM1,file=paste0(path, "Results/", "AM1_results.RData"))

AM3 <- lapply(1:100, estFun, method = "converge", n = n, n_i = n_i, alpha = alpha, beta = beta, betas = betas, betals = 0, betatimeind = betatimeind, lambda = lambda,
              betat = betat, noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(AM3,file=paste0(path, "Results/", "AM3_results.RData"))

AM5 <- lapply(1:100, estFun, method = "converge", n = n, n_i = n_i, alpha = alpha, beta = beta, betas = betas, betals = betals, betatimeind = betatimeind, lambda = lambda,
              betat = betat, noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(AM5,file=paste0(path, "Results/", "AM5_results.RData"))


## w/ fixed time effect
betas = 0.1
beta = c(1.5,-0.5, 0.7, 1.3, 0.3, 0.5)
betals = c(0.4, 0.9, 0.3, -1, 0.2, -0.4)
betatimeind = 1

i <- 0

AM2 <- lapply(1:1, estFun, method = "converge", n = n, n_i = n_i, alpha = alpha, beta = beta[1], betas = betas, betals = betals[1], betatimeind = betatimeind, lambda = lambda,
              noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(AM2,file=paste0(path, "Results/", "AM2_results.RData"))

AM4 <- lapply(1:1, estFun, method = "converge", n = n, n_i = n_i, alpha = alpha, beta = beta, betas = betas, betals = betals[1], betatimeind = betatimeind, lambda = lambda,
              noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(AM4,file=paste0(path, "Results/", "AM4_results.RData"))

AM6 <- lapply(1:1, estFun, method = "converge", n = n, n_i = n_i, alpha = alpha, beta = beta, betas = betas, betals = betals, betatimeind = betatimeind, lambda = lambda,
              noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(AM6,file=paste0(path, "Results/", "AM6_results.RData"))

## VARIABLE SELECTION ##############################
## V1 ##############################################
n = 100
n_i = 5
alpha = 0.1
lambda = 0.6

noninf <- noninfs <- noninfls <- 5
# 
## w/o fixed time effect
betas = 0.1
beta = c(1.5,-0.5, 0.7, 1.3, 0.3, 0.5)
betat = 0.4
betals = c(0.9, 0.3, -1, 0.2, -0.4)
betatimeind = 0

i <- 0

V1M1 <- lapply(1:100, estFun, n = n, n_i = n_i, alpha = alpha, beta = beta[1], betas = betas, betals = 0, betatimeind = betatimeind, lambda = lambda,
                 betat = betat, noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(V1M1,file=paste0(path, "Results/", "V1M1_results.RData"))

V1M3 <- lapply(1:100, estFun, n = n, n_i = n_i, alpha = alpha, beta = beta, betas = betas, betals = 0, betatimeind = betatimeind, lambda = lambda,
              betat = betat, noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(V1M3,file=paste0(path, "Results/", "V1M3_results.RData"))

V1M5 <- lapply(1:100, estFun, n = n, n_i = n_i, alpha = alpha, beta = beta, betas = betas, betals = betals, betatimeind = betatimeind, lambda = lambda,
              betat = betat, noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(V1M5,file=paste0(path, "Results/", "V1M5_results.RData"))


## w/ fixed time effect
betas = 0.1
beta = c(1.5,-0.5, 0.7, 1.3, 0.3, 0.5)
betals = c(0.4, 0.9, 0.3, -1, 0.2, -0.4)
betatimeind = 1

i <- 0

V1M2 <- lapply(1:100, estFun, n = n, n_i = n_i, alpha = alpha, beta = beta[1], betas = betas, betals = betals[1], betatimeind = betatimeind, lambda = lambda,
              noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(V1M2,file=paste0(path, "Results/", "V1M2_results.RData"))

V1M4 <- lapply(1:100, estFun, n = n, n_i = n_i, alpha = alpha, beta = beta, betas = betas, betals = betals[1], betatimeind = betatimeind, lambda = lambda,
              noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(V1M4,file=paste0(path, "Results/", "V1M4_results.RData"))

V1M6 <- lapply(1:100, estFun, n = n, n_i = n_i, alpha = alpha, beta = beta, betas = betas, betals = betals, betatimeind = betatimeind, lambda = lambda,
              noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(V1M6,file=paste0(path, "Results/", "V1M6_results.RData"))

         ## ############################################################# ##
         ## ####################### Evaluation V1 ####################### ##
         ## ############################################################# ##

## V2 ##############################################
n = 100
n_i = 5
alpha = 0.1
lambda = 0.6

noninf <- noninfs <- noninfls <- 50
# 
## w/o fixed time effect
betas = 0.1
beta = c(1.5,-0.5, 0.7, 1.3, 0.3, 0.5)
betat = 0.4
betals = c(0.9, 0.3, -1, 0.2, -0.4)
betatimeind = 0

i <- 0

V2M1 <- lapply(1:100, estFun, n = n, n_i = n_i, alpha = alpha, beta = beta[1], betas = betas, betals = 0, betatimeind = betatimeind, lambda = lambda,
               betat = betat, noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(V2M1,file=paste0(path, "Results/", "V2M1_results.RData"))

V2M3 <- lapply(1:100, estFun, n = n, n_i = n_i, alpha = alpha, beta = beta, betas = betas, betals = 0, betatimeind = betatimeind, lambda = lambda,
              betat = betat, noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(V2M3,file=paste0(path, "Results/", "V2M3_results.RData"))

V2M5 <- lapply(1:100, estFun, n = n, n_i = n_i, alpha = alpha, beta = beta, betas = betas, betals = betals, betatimeind = betatimeind, lambda = lambda,
              betat = betat, noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(V2M5,file=paste0(path, "Results/", "V2M5_results.RData"))


## w/ fixed time effect
betas = 0.1
beta = c(1.5,-0.5, 0.7, 1.3, 0.3, 0.5)
betals = c(0.4, 0.9, 0.3, -1, 0.2, -0.4)
betatimeind = 1

i <- 0

V2M2 <- lapply(1:100, estFun, n = n, n_i = n_i, alpha = alpha, beta = beta[1], betas = betas, betals = betals[1], betatimeind = betatimeind, lambda = lambda,
              noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(V2M2,file=paste0(path, "Results/", "V2M2_results.RData"))

V2M4 <- lapply(1:100, estFun, n = n, n_i = n_i, alpha = alpha, beta = beta, betas = betas, betals = betals[1], betatimeind = betatimeind, lambda = lambda,
              noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(V2M4,file=paste0(path, "Results/", "V2M4_results.RData"))

V2M6 <- lapply(1:100, estFun, n = n, n_i = n_i, alpha = alpha, beta = beta, betas = betas, betals = betals, betatimeind = betatimeind, lambda = lambda,
              noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(V2M6,file=paste0(path, "Results/", "V2M6_results.RData"))

          ## ############################################################# ##
          ## ####################### Evaluation V2 ####################### ##
          ## ############################################################# ##

## V3 ##############################################
n = 100
n_i = 5
alpha = 0.1
lambda = 0.6

noninf <- noninfs <- noninfls <- 250
# 
## w/o fixed time effect
betas = 0.1
beta = c(1.5,-0.5, 0.7, 1.3, 0.3, 0.5)
betat = 0.4
betals = c(0.9, 0.3, -1, 0.2, -0.4)
betatimeind = 0

i <- 0

V3M1 <- lapply(1:100, estFun, n = n, n_i = n_i, alpha = alpha, beta = beta[1], betas = betas, betals = 0, betatimeind = betatimeind, lambda = lambda,
               betat = betat, noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(V3M1,file=paste0(path, "Results/", "V3M1_results.RData"))

V3M3 <- lapply(1:100, estFun, n = n, n_i = n_i, alpha = alpha, beta = beta, betas = betas, betals = 0, betatimeind = betatimeind, lambda = lambda,
               betat = betat, noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(V3M3,file=paste0(path, "Results/", "V3M3_results.RData"))

V3M5 <- lapply(1:100, estFun, n = n, n_i = n_i, alpha = alpha, beta = beta, betas = betas, betals = betals, betatimeind = betatimeind, lambda = lambda,
               betat = betat, noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(V3M5,file=paste0(path, "Results/", "V3M5_results.RData"))


## w/ fixed time effect
betas = 0.1
beta = c(1.5,-0.5, 0.7, 1.3, 0.3, 0.5)
betals = c(0.4, 0.9, 0.3, -1, 0.2, -0.4)
betatimeind = 1

i <- 0 

V3M2 <- lapply(1:100, estFun, n = n, n_i = n_i, alpha = alpha, beta = beta[1], betas = betas, betals = betals[1], betatimeind = betatimeind, lambda = lambda,
               noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(V3M2,file=paste0(path, "Results/", "V3M2_results.RData"))

V3M4 <- lapply(1:100, estFun, n = n, n_i = n_i, alpha = alpha, beta = beta, betas = betas, betals = betals[1], betatimeind = betatimeind, lambda = lambda,
               noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(V3M4,file=paste0(path, "Results/", "V3M4_results.RData"))

V3M6 <- lapply(1:100, estFun, n = n, n_i = n_i, alpha = alpha, beta = beta, betas = betas, betals = betals, betatimeind = betatimeind, lambda = lambda,
               noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(V3M6,file=paste0(path, "Results/", "V3M6_results.RData"))

            ## ############################################################# ##
            ## ####################### Evaluation V3 ####################### ##
            ## ############################################################# ##

## Variable Failure ##############################################

n = 100
n_i = 5
alpha = 0.1
lambda = 0.6

noninf <- noninfs <- noninfls <- 250
# 
## w/o fixed time effect
betas = 0.1
beta = c(1.5,-0.5, 0.7, 1.3, 0.3, 0.5)
betat = 0.4
betals = c(0.9, 0.3, -1, 0.2, -0.4)
betatimeind = 0

VFM1 <- lapply(1:100, failFun, n = n, n_i = n_i, alpha = alpha, beta = beta[1], betas = betas, betals = 0, betatimeind = betatimeind, lambda = lambda,
               betat = betat, noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(VFM1,file=paste0(path, "Results/", "VFM1_results.RData"))

VFM3 <- lapply(1:100, failFun, n = n, n_i = n_i, alpha = alpha, beta = beta, betas = betas, betals = 0, betatimeind = betatimeind, lambda = lambda,
               betat = betat, noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(VFM3,file=paste0(path, "Results/", "VFM3_results.RData"))

VFM5 <- lapply(1:100, failFun, n = n, n_i = n_i, alpha = alpha, beta = beta, betas = betas, betals = betals, betatimeind = betatimeind, lambda = lambda,
               betat = betat, noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(VFM5,file=paste0(path, "Results/", "VFM5_results.RData"))

## w/ fixed time effect
betas = 0.1
beta = c(1.5,-0.5, 0.7, 1.3, 0.3, 0.5)
betals = c(0.4, 0.9, 0.3, -1, 0.2, -0.4)
betatimeind = 1

VFM2 <- lapply(1:100, failFun, n = n, n_i = n_i, alpha = alpha, beta = beta[1], betas = betas, betals = betals[1], betatimeind = betatimeind, lambda = lambda,
               noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(VFM2,file=paste0(path, "Results/", "VFM2_results.RData"))

VFM4 <- lapply(1:100, failFun, n = n, n_i = n_i, alpha = alpha, beta = beta, betas = betas, betals = betals[1], betatimeind = betatimeind, lambda = lambda,
               noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(VFM4,file=paste0(path, "Results/", "VFM4_results.RData"))

VFM6 <- lapply(1:100, failFun, n = n, n_i = n_i, alpha = alpha, beta = beta, betas = betas, betals = betals, betatimeind = betatimeind, lambda = lambda,
               noninf = noninf, noninfs = noninfs, noninfls = noninfls)
save(VFM6,file=paste0(path, "Results/", "VFM6_results.RData"))

## Prediction Precision ##############################################
if (!require("JM")) install.packages("JM")
if (!require("parallel")) install.packages("parallel")

predFun <- function (x, estres) {
  tryCatch({
    set.seed(541689 + x)

    df <- simJM(n = n, n_i = n_i, alpha = alpha, beta = beta, betas = betas, betals = betals, betatimeind = betatimeind, lambda = lambda, 
                 betat = betat, noninf = noninf, noninfs = noninfs, noninfls = noninfls)
 
    # first = rep(FALSE, length(df$y))
    # for(i in 1:length(df$y)){
    #   first[which.max(df$id==i)] = TRUE
    # }
    # Xls_un <-  as.matrix(df$Xls[first==1,])
    # true.etas <-  as.vector(df$Xs[,1] %*% betas)
    # true.etals_un <- as.vector(as.vector(Xls_un%*%betals) + betat*T_surv)
    # true.risk <- lambda * exp(true.etas + alpha * true.etals_um)
    # 
    # etals <- df$Xls[first == 1, 2:(ncol(df$Xls) - noninfls)] %*% betals[2:6] + df$T_surv * betals[1] + rowSums(cbind(1, df$T_surv) * df$R_mean[first == 1,])
    # 
    # 
    # y <- df$y
    # Xls <- cbind(df$X, df$Xls)
    # Xs <- df$Xs
    # delta <- df$delta
    # T_surv <- df$T_surv
    # id <- df$id
    # true.risk <- df$true.risk
    # 
    # res <- estres[[x]][[1]] 
    # betal <- res$coefficients$betas[1]
    # betals <- res$coefficients$betas[-1]
    # betas <- res$coefficients$gammas
    # alpha <- res$coefficients$alpha
    # sigma2 <- res$coefficients$sigma
    # 
    # lambda0_i <- numeric(length(unique(df$id)))
    # for (i in 1:length(lambda0_i)) {
    #   baseHaz <- cbind(c(0, res$control$knots), res$coefficients$xi)
    #   pos <- max(which((df$T_surv[i] - baseHaz[,1]) > 0))
    #   lambda0_i[i] <- baseHaz[pos, 2]
    # }
    # 
    # etal <- as.matrix(rep(1, length(y))) %*% betal
    # etals <- Xls %*% betals
    # etas <- Xs %*% betas
    # 
    # first  <- rep(FALSE, length(id))
    # for(i in unique(id)){
    #   first[which.max(id == i)] = TRUE
    # }
    # 
    # betatimeind <- which(colnames(Xls) == "time")
    # beta_t <- betals[betatimeind]
    # etals_m <- (as.matrix(Xls[,-betatimeind]) %*% as.matrix(betals[-betatimeind]))[first==1] + T_surv * beta_t
    # 
    # y.hat <- etal + etals
    # MSEP.y <- 1/length(y) * sum((y - y.hat)^2) 
    # 
    # risk.hat <- lambda0_i * exp(etas + alpha * etals_m)
    # MSEP.risk <- 1/length(y) * sum((true.risk - risk.hat)^2)
    # 
    # 
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
