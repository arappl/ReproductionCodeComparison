
# ## FUNCTIONS ###########################################################
cvini <- function(x, n, n_i, alphasim, beta, betas, betals, betatimeind, lambda, betat = 0, noninf = 0, noninfs = 0, noninfls = 0, 
                  alphaini = .001, k = 10, grid1 = NULL, grid2 = NULL, grid3 = NULL){

  set.seed(seeds[x])
  df <- simJM(n = n, n_i = n_i, alphasim = alpha, beta = beta, betas = betas, betals = betals, betatimeind = betatimeind, lambda = lambda, 
              betat = betat, noninf = noninf, noninfs = noninfs, noninfls = noninfls)
  
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
  set.seed(1258149)
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
                               betatimeind = betatimeind, alpha = alphaini)},
                       error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    indarr <- indarr + cv.res$indarr
    resarr <- resarr + cv.res$likarr
    
    print(kk)
  }
  
  ## 3.3 combine results and find best set of iterations -------------------------------------- ##
  meanarr <- resarr/(k - indarr) # divide by the number of folds, which worked
  best <- which(meanarr == max(meanarr), arr.ind = T)
  best.iter <- c(grid1[best[1]], grid2[best[2]], grid3[best[3]])
  names(best.iter) <- c("ml","ms","mls")
  return(list("best" = best.iter, "indarr" = indarr))
}

like_cv = function(y, betal, betas, betals, Xl, Xs, Xls,
                   betatimeind, delta, T_surv, alpha, lambda, sigma2, id) {
  first = rep(FALSE, length(id))
  for(i in unique(id)){
    first[which.max(id==i)] = TRUE
  }
  if(is.null(Xl)){
    Xl = as.matrix(rep(1, length(y)))
  }
  beta_t = betals[betatimeind]
  etals_m = (as.matrix(Xls[,-betatimeind])%*%as.matrix(betals[-betatimeind]))[first==1]
  etals = Xls%*%betals
  etal = cbind(1,Xl)%*%betal
  etas = Xs%*%betas
  
  if(length(beta_t) != 0) {
    if (beta_t != 0) {
      integral = lambda*exp(etas)*((exp(alpha*etals_m + alpha*beta_t*T_surv) - exp(alpha*etals_m))*(1/(alpha*beta_t)))
    }else{
      integral = lambda*T_surv*exp(etas + alpha*etals_m)
    }
  } else {
    integral = lambda*T_surv*exp(etas + alpha*etals_m)
  } 
  
  surv = delta*(log(lambda) + etas + alpha*etals_m + alpha*beta_t*T_surv) - integral
  long = log(1/sqrt(2*pi*sigma2)) - (y - etals - etal)^2/(2*sigma2)
  like = sum(surv) + sum(long)
  return(like)
}

cvres2 <- function(data.train, data.pred, grid1 = grid1, grid2 = grid2, grid3 = grid3, betatimeind, alpha = alpha){
  
  likarr <- indarr <- array(0,c(length(grid1), length(grid2), length(grid3)),
                            dimnames = list(grid1 ,grid2, grid3))
  # grid1=ml, grid2=ms, grid3=mls
  
  
  for(ml_akt in grid1){
    for(ms_akt in grid2){
      for(mls_akt in grid3){
        # print(c(ml_akt, ms_akt, mls_akt))
        mod <- tryCatch({JMboost(y = data.train$y, Xl = data.train$X, Xs = data.train$Xs, Xls = data.train$Xls, delta = data.train$delta,
                                 T_long = data.train$T_long, T_surv = data.train$T_surv, id = data.train$id,
                                 mstop_l = ml_akt, mstop_s = ms_akt, mstop_ls = mls_akt, betatimeind = betatimeind, alpha = alpha, verbose = F)},
                        error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        
        like <-  tryCatch({like_cv(y=data.pred$y,  Xl=data.pred$X, Xs=data.pred$Xs, Xls=data.pred$Xls, delta=data.pred$delta,
                                   T_surv = data.pred$T_surv, id=data.pred$id, betatimeind=betatimeind, betal=as.matrix(mod$betal),
                                   betas=as.matrix(mod$betas) ,betals=as.matrix(mod$betals), alpha=mod$alpha, lambda=mod$lambda,
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
  # return(diaarr)
  return(list("likarr" = likarr, "indarr" = indarr))
}

simJM <- function(n, n_i, alphasim,
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
  }else{Xls =matrix(1,nrow=n*n_i, ncol=1); betatime = 0} #  matrix(rep(rnorm(n), each = n_i), ncol=1, nrow=n*n_i); betatime = 0}
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
  T_surv = (log((-log(1-u)*alphasim*(betatime+R_mean[first==1,2]))/(lambda*exp(eta_s_mean))+exp(alphasim*(R_mean[first==1,1] + Xls[first==1, , drop=FALSE]%*%betals)))-alphasim*(R_mean[first==1,1] + Xls[first==1, , drop=FALSE]%*%betals))/(alphasim*(betatime+R_mean[first==1,2]))
  T_surv[is.nan(T_surv)] = 2
  time_mat <- matrix(nrow = n_i, data = time)
  
  delta = rep(1, n)
  for(i in 1:n){
    if(!is.na(T_surv[i])){
      if(T_surv[i]>max(time_mat[,i])){T_surv[i] = max(time_mat[,i]); delta[i] = 0}
      # if(T_surv[i]>1){T_surv[i] = max(time_mat[,i]); delta[i] = 0} # funktioniert nicht mit dieser Parameterkombination, es fallen Individuen schon bei der ersten Beobachtung raus
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
      }else if(betatimeind == 1){ # neue Bedingung hinzugefuegt
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
  
  # if(is.null(Xls)==FALSE){
  #   first = rep(FALSE, length(y))
  #   for(i in 1:length(y)){
  #     first[which.max(id==i)] = TRUE
  #   }
  #   Xls_un = Xls[first==1,]
  #   Xls_un[,betatimeind] = T_surv
  # }
  # 
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
                 betat=betat, int=int, alpha = alpha, lambda = lambda,  sigma2 = sigma2))
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
    
    if (m == mstop || m%%1000 == 0) {
      intermediate <- list(GAMMA0 = GAMMA0, GAMMA1 = GAMMA1, BETAL = BETAL, BETAS = BETAS,
                           BETALS = BETALS, BETAT=BETAT, INT=INT, ALPHA = ALPHA, LAMBDA = LAMBDA, SIGMA2 = SIGMA2,
                           gamma0 = gamma0, gamma1 = gamma1, betal = betal, betas = betas, betals = betals,
                           betat=betat, int=int, alpha = alpha, lambda = lambda,  sigma2 = sigma2)
      save(intermediate, file = paste0(path, "01_Datenbeispiel/03_Converge_Zwischenspeicher/Dat_conv_m", m, ".RData"))
      rm("intermediate")
    }
    
    if(verbose){
      if(m%%1000 == 0){print(m)}
      # print(m)
    }
    
  }
  conv.iter <- max(which(fitJMbcon$INT != 0)) + 1
  
  structure(list(GAMMA0 = GAMMA0, GAMMA1 = GAMMA1, BETAL = BETAL, BETAS = BETAS,
                 BETALS = BETALS, BETAT=BETAT, INT=INT, ALPHA = ALPHA, LAMBDA = LAMBDA, SIGMA2 = SIGMA2,
                 gamma0 = gamma0, gamma1 = gamma1, betal = betal, betas = betas, betals = betals,
                 betat=betat, int=int, alpha = alpha, lambda = lambda,  sigma2 = sigma2, 
                 conv.iter = conv.iter))
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
