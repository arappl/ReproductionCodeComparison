estfun <- function(x, model, n, n_i, alpha, beta, betas, betals, betatimeind, lambda, betat = 0, noninf = 0, noninfs = 0, noninfls = 0, method = "boost",
                   k = 10, grid1 = NULL, grid2 = NULL, grid3 = NULL) # method = c("boost", "converge")
  {
  repeat{
    # simulate data with specific seed, if fits fail, change the seed (i <-  i + 1)
    i <- i + 1
    seed <- 8469 + x*1000 + i
    set.seed(seed)
    dat <- simJM(n = n, n_i = n_i, alpha = alpha, beta = beta, betas = betas, betals = betals, betatimeind = betatimeind, lambda = lambda, 
                 betat = betat, noninf = noninf, noninfs = noninfs, noninfls = noninfls)
    
    df <- with(dat, data.frame(id, y, X, Xls))
    df.id <- with(dat, data.frame("id" = unique(id), T_surv, delta, Xs))
    
    # little workaround to be able to run JM flexibly across various models: jointModel() resorts to the formula call of lme(). If 
    # this formula is passed down through a function to lme() (like in mjoint()), the lme output stores a string of this function. When 
    # jointModel() tries to convert this string into a formula again it will run into an error. Therefore we'll store a code snippet with 
    # the model specific formula in a seperate file and source it back in. Thus jointModel() can access a correct formula and we can only 
    # one function instead of six (or more respectively) versions of it :).
    if (x == 1) { 
      file_ini <- file(paste0(path, "/", sub_dir, "/", "Model_", model, ".txt"))
      
      writeLines(c(
        "lmeIni <- function(df = df){",
        paste0("fitLME <- tryCatch({ lme(", paste0("y~", paste(names(df)[3 : ncol(df)], collapse = "+")),","),
        "                             data = df, random = ~ time | id, na.action = na.exclude, control = lmeControl(opt='optim')) }, ",
        paste0("           error=function(e){cat(","'ERROR :'", ", conditionMessage(e),"," '/n'",")}) "),
        "}"
      )
      , file_ini)
      
      close(file_ini)
    }
    
    # +++++ JM fit +++++ # 
    source(paste0(path, "/", sub_dir, "/", "Model_", model, ".txt"))
    fitLME <- lmeIni(df = df)
    fitSURV <- coxph(Surv(time = T_surv, event=delta) ~ . -id, data = df.id, x = TRUE, na.action = na.exclude)
    fitJM <-  tryCatch({ jointModel(fitLME, fitSURV, timeVar = "time", method="piecewise-PH-aGH") }, 
                       error=function(e){cat("ERROR :",conditionMessage(e), "/n")})
    
    
    # +++++ joineRML fit +++++ # 
    # if fitJM returnd no error, fit joineRML
    if (!is.error(fitJM) && !is.null(fitJM)) {
      fitJRML <- tryCatch({mjoint(formLongFixed = list(as.formula(paste0("y~", paste(names(df)[3:(ncol(df)-3)], collapse="+")))),
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