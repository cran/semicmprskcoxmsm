######## Get the initial Value #######
initial_fit_em_weights <- function(data,X1,X2,event1,event2,w,Trt){

  data$X1 <- data[,X1]
  data$X2 <- data[,X2]
  data$delta1 <- data[,event1]
  data$delta2 <- data[,event2]
  data$A <- data[,Trt]

  event1 <- Surv(time = data$X1, event = data$delta1)
  data$X2.1 <- ifelse((1-data$delta1) * data$delta2 == 1, data$X2, data$X1)
  event2 <- Surv(time = data$X2.1, event = (1-data$delta1) * data$delta2 )

  data.1 <- data[which(data$delta1==1),]
  w.1 <- w[which(data$delta1==1)]
  event3 <- Surv(time = data.1$X1, time2 = data.1$X2, event = (data.1$delta1) * (data.1$delta2))

  ### Fit the weighted cox model:
  fit1 <- coxph(event1 ~ A, data = data, ties = "efron", weights = w)
  fit2 <- coxph(event2 ~ A, data = data, ties = "efron", weights = w)
  fit3 <- coxph(event3 ~ A, data = data.1, ties = "efron", weights = w.1)

  OUT <- list(event1 = event1,
              event2 = event2,
              event3 = event3,
              fit1 = fit1,
              fit2 = fit2,
              fit3 = fit3)
  return(OUT)
}

get_hazard <- function(fit){
  cumbasehaz <- basehaz(fit,centered = FALSE)
  basehaz <- diff(cumbasehaz$hazard)
  basehaz1 <- cumbasehaz[-1,]
  basehaz1$basehaz <- basehaz
  lambda <- basehaz1[which(basehaz1$basehaz!=0),c("time","basehaz")]
  if(nrow(lambda)<fit$nevent){
    a <- c(cumbasehaz[1,2],cumbasehaz[1,1])
    lambda <- rbind(a,lambda)
  }
  Lambda <- stepfun(cumbasehaz$time,c(0,cumbasehaz$hazard))
  OUT <- list(lambda,Lambda)
  return(OUT)
}

initial_lambda_em <- function(OUT){
  fit1 <- OUT$fit1
  fit2 <- OUT$fit2
  fit3 <- OUT$fit3

  lambda01_func <- get_hazard(fit1)
  lambda02_func <- get_hazard(fit2)
  lambda03_func <- get_hazard(fit3)

  lambda1 <- lambda01_func[[1]]
  Lambda1 <- lambda01_func[[2]]

  lambda2 <- lambda02_func[[1]]
  Lambda2 <- lambda02_func[[2]]

  lambda3 <- lambda03_func[[1]]
  Lambda3 <- lambda03_func[[2]]

  OUT <- list(lambda1 = lambda1, Lambda1 = Lambda1,
              lambda2 = lambda2, Lambda2 = Lambda2,
              lambda3 = lambda3, Lambda3 = Lambda3)

}

OUT_em_weights <- function(data,X1,X2,event1,event2,w,Trt){

  ### Get the initial estimate of beta
  fit_0 <- initial_fit_em_weights(data,X1,X2,event1,event2,w,Trt)
  beta1 <- fit_0$fit1$coefficients
  beta2 <- fit_0$fit2$coefficients
  beta3 <- fit_0$fit3$coefficients

  ### Get the initial estimate of lambda0 and Lambda0
  basehaz_0 <- initial_lambda_em(OUT = fit_0)
  lambda1 <- basehaz_0$lambda1
  lambda2 <- basehaz_0$lambda2
  lambda3 <- basehaz_0$lambda3
  Lambda1 <- basehaz_0$Lambda1
  Lambda2 <- basehaz_0$Lambda2
  Lambda3 <- basehaz_0$Lambda3

  ### Get the event
  event1 <- fit_0$event1
  event2 <- fit_0$event2
  event3 <- fit_0$event3

  RES <- list(beta1=beta1,beta2=beta2,beta3=beta3,
              lambda1 = lambda1,lambda2 = lambda2,lambda3 = lambda3,
              Lambda1 = Lambda1,Lambda2 = Lambda2,Lambda3 = Lambda3,
              event1 = event1, event2 = event2, event3 = event3)

  return(RES)

}

get_hazard_offset_weights <- function(fit,data,time1 = NULL,time2,w){

  beta <- fit$coefficients
  base_haz_est <- get_hazard(fit)[[1]]
  base_haz_est$basehaz <- NA

  if(is.null(time1)){
    for (j in 1:nrow(base_haz_est)) {
      base_haz_est$basehaz[j] <- w[which(data[,time2]==base_haz_est$time[j])]/sum( as.numeric(base_haz_est$time[j] <= data[,time2]) * exp(beta*data[,"A"] + fit$offset) * w  )
    }}else{
      for (j in 1:nrow(base_haz_est)) {
        base_haz_est$basehaz[j] <- w[which(data[,time2]==base_haz_est$time[j])]/sum( as.numeric( data[,time1] < base_haz_est$time[j] & base_haz_est$time[j] <= data[,time2] ) *
                                                                                       exp(beta*data[,"A"] + fit$offset) * w  )
      }
    }

  ### Calculate the cumulative baseline hazard:
  cum_base_haz <- as.data.frame(basehaz(fit)$time)
  colnames(cum_base_haz) <- "time"
  cum_base_haz$cum_haz <- NA
  cum_base_haz$base_haz <- NA
  for (j in 1:nrow(cum_base_haz)) {
    if(cum_base_haz$time[j] %in% base_haz_est$time){
      cum_base_haz$base_haz[j] <- base_haz_est$basehaz[which(base_haz_est$time==cum_base_haz$time[j])]
    }else{
      cum_base_haz$base_haz[j] <- 0
    }
  }
  cum_base_haz$cum_haz <- cumsum(cum_base_haz$base_haz)

  Lambda <- stepfun(cum_base_haz$time,c(0,cum_base_haz$cum_haz))
  OUT <- list(lambda = base_haz_est,Lambda = Lambda,cum_base_haz = cum_base_haz)
  return(OUT)

}

em_illness_death_phmm_weight <- function(data,X1,X2,event1,event2,w,Trt,
                                         EM_initial,sigma_2_0){
  rule <- gaussHermiteData(100)
  diff1 <- 1
  diff2 <- 1
  em.n <- 0

  data$X1 <- data[,X1]
  data$X2 <- data[,X2]
  data$delta1 <- data[,event1]
  data$delta2 <- data[,event2]
  data$A <- data[,Trt]
  data$X2.1 <- ifelse((1-data$delta1) * data$delta2 == 1, data$X2, data$X1)

  w.1 <- w[which(data[,event1]==1)]

  beta1 <- EM_initial$beta1
  beta2 <- EM_initial$beta2
  beta3 <- EM_initial$beta3
  sigma_2 <- sigma_2_0

  lambda1 <- EM_initial$lambda1
  lambda2 <- EM_initial$lambda2
  lambda3 <- EM_initial$lambda3

  Lambda1 <- EM_initial$Lambda1
  Lambda2 <- EM_initial$Lambda2
  Lambda3 <- EM_initial$Lambda3

  event1 <- EM_initial$event1
  event2 <- EM_initial$event2
  event3 <- EM_initial$event3

  beta1_est<- NULL
  beta1_est[1] <- beta1
  beta2_est<- NULL
  beta2_est[1] <- beta2
  beta3_est<- NULL
  beta3_est[1] <- beta3
  sigma_2_est<- NULL
  sigma_2_est[1] <- sigma_2

  log_lik_y <- NULL


  while (diff1 > 10^-5 & diff2 > 10^-3) {
    em.n <- em.n + 1
    data$L_y_theta <- NA
    data$log.eeb <- 0
    data$eb2 <- NA
    data$Lam1 <- NA
    data$Lam2 <- NA
    data$Lam3 <- NA
    data$lam01 <- NA
    data$lam02 <- NA
    data$lam03 <- NA

    ###### Get the Q function out: Calculate logE(e^b), E(b^2), E(b)
    for (i in 1:nrow(data)) {

      A_i <- data$A[i]
      delta1_i <- data$delta1[i]
      delta2_i <- data$delta2[i]
      data$Lam1[i] <-  Lambda1(data$X1[i])
      data$Lam2[i]  <- Lambda2(data$X1[i])
      data$Lam3[i]  <- Lambda3(data$X2[i]) - Lambda3(data$X1[i])
      data$lam01[i] <- ifelse(delta1_i==1,
                              lambda1$basehaz[which(lambda1$time==data$X1[i])],0)
      data$lam02[i] <- ifelse((1-delta1_i)*delta2_i==1,
                              lambda2$basehaz[which(lambda2$time==data$X2[i])],0)
      data$lam03[i] <- ifelse(delta1_i*delta2_i==1,
                              lambda3$basehaz[which(lambda3$time==data$X2[i])],0)

      ##### First: calculate f(y_i)
      f_b_y_theta_new <- function(b){
        ( data$lam01[i] * exp(beta1 * A_i + b) )^{delta1_i} *              exp(-data$Lam1[i] * exp(beta1 * A_i + b) ) *
          ( data$lam02[i] * exp(beta2 * A_i + b) )^{(1-delta1_i)*delta2_i} * exp(-data$Lam2[i] * exp(beta2 * A_i + b) ) *
          ( data$lam03[i] * exp(beta3 * A_i + b) )^{delta1_i*delta2_i} *     exp(-data$Lam3[i] * exp(beta3 * A_i + b) ) *
          ( exp( -(b^2) / (2*(sigma_2)) ) / sqrt( 2*pi*(sigma_2) ) )
      }

      f_y_theta_i <- aghQuad(f_b_y_theta_new,muHat = 0, sigmaHat = 1, rule)

      ##### Second: calculate E(e^b | y)
      E_eb_formu <- function(b){
        exp(b) * ( data$lam01[i] * exp(beta1 * A_i + b) )^{delta1_i} *              exp(-data$Lam1[i] * exp(beta1 * A_i + b) ) *
          ( data$lam02[i] * exp(beta2 * A_i + b) )^{(1-delta1_i)*delta2_i} * exp(-data$Lam2[i] * exp(beta2 * A_i + b) ) *
          ( data$lam03[i] * exp(beta3 * A_i + b) )^{delta1_i*delta2_i} *     exp(-data$Lam3[i] * exp(beta3 * A_i + b) ) *
          ( exp( -(b^2) / (2*(sigma_2)) ) / sqrt( 2*pi*(sigma_2) ) ) / f_y_theta_i
      }

      E_eb_i <- aghQuad(E_eb_formu,muHat = 0, sigmaHat = 1, rule)

      ###### Third: Calculate E(b^2 | y)
      E_b2_formu <- function(b){
        (b^2) * ( data$lam01[i] * exp(beta1 * A_i + b) )^{delta1_i} *              exp(-data$Lam1[i] * exp(beta1 * A_i + b) ) *
          ( data$lam02[i] * exp(beta2 * A_i + b) )^{(1-delta1_i)*delta2_i} * exp(-data$Lam2[i] * exp(beta2 * A_i + b) ) *
          ( data$lam03[i] * exp(beta3 * A_i + b) )^{delta1_i*delta2_i} *     exp(-data$Lam3[i] * exp(beta3 * A_i + b) ) *
          ( exp( -(b^2) / (2*(sigma_2)) ) / sqrt( 2*pi*(sigma_2) ) ) / f_y_theta_i
      }

      E_b2_i <- aghQuad(E_b2_formu,muHat = 0, 1, rule)

      data$L_y_theta[i] <- f_y_theta_i
      data$log.eeb[i] <- log(E_eb_i)
      data$eb2[i] <- E_b2_i

      ###### Stop the E step

    }

    if(sum(is.nan(data$log.eeb))>0 | sum(is.na(data$log.eeb))>0){
      diff1 <- 10^{-7}
      diff2 <- 10^{-7}
    }else{

      log_lik_y[em.n] <- sum(log(data$L_y_theta))
      #### M-step:

      #print(paste("# of EM:",em.n))

      data.1 <- data[which(data$delta1==1),]

      fit1 <- coxph(event1 ~ A + offset(log.eeb),data = data,ties = "efron",weights = w)
      fit2 <- coxph(event2 ~ A + offset(log.eeb),data = data,ties = "efron",weights = w)
      fit3 <- coxph(event3 ~ A + offset(log.eeb),data = data.1,ties = "efron",weights = w.1)

      beta1_new <- coef(fit1)
      beta2_new <- coef(fit2)
      beta3_new <- coef(fit3)

      lambda0_est1 <- get_hazard_offset_weights(fit = fit1, data = data, time1 = NULL, time2 = "X1",w = w)
      lambda0_est2 <- get_hazard_offset_weights(fit = fit2, data = data, time1 = NULL, time2 = "X2.1",w = w)
      lambda0_est3 <- get_hazard_offset_weights(fit = fit3, data = data.1, time1 = "X1", time2 = "X2",w = w.1)

      lambda1 <- lambda0_est1$lambda
      Lambda1 <- lambda0_est1$Lambda
      lambda2 <- lambda0_est2$lambda
      Lambda2 <- lambda0_est2$Lambda
      lambda3 <- lambda0_est3$lambda
      Lambda3 <- lambda0_est3$Lambda

      sigma_new <- sum(w*data$eb2) / sum(w)

      if(em.n==1){
        diff1 <- diff2 <- max(abs((beta1-beta1_new)),abs((beta2-beta2_new)),
                              abs((beta3-beta3_new)),abs((sigma_2-sigma_new)))
      }else{
        diff1 <-  abs( ( log_lik_y[em.n]  - log_lik_y[em.n-1] ) / log_lik_y[em.n] )

        diff2 <-  max( abs((beta1-beta1_new)), abs((beta2-beta2_new)),
                       abs((beta3-beta3_new)), abs((sigma_2-sigma_new)) )
      }

      beta1 <- beta1_new
      beta2 <- beta2_new
      beta3 <- beta3_new
      sigma_2 <- sigma_new

      beta3_est[em.n+1] <- beta3
      beta2_est[em.n+1] <- beta2
      beta1_est[em.n+1] <- beta1
      sigma_2_est[em.n+1] <- sigma_2


    }

  }

  res1 <- list(beta1 = beta1_est,beta2 = beta2_est, beta3 = beta3_est,
               Lambda01 = Lambda1,Lambda02 = Lambda2,Lambda03 = Lambda3,
               sigma_2 = sigma_2_est,
               loglik = log_lik_y,
               em.n = em.n,
               data = data)

  return(res1)

}

var_em_illness_death_phmm <- function(data,sigma_2_0,VARS.){

  beta1_boot <- NULL
  beta2_boot <- NULL
  beta3_boot <- NULL
  sigma_2_boot <- NULL

  Lam01_1_boot <- NULL
  Lam02_1_boot <- NULL
  Lam03_1_boot <- NULL

  n <- nrow(data)

  for (b in 1:100) {

    id.boot <- sample(1:n,n,replace = TRUE)
    data_boot <- data[id.boot,]
    data_boot <- data_boot[order(data_boot$id),]
    a <- unique(data_boot$id[duplicated(data_boot)])
    data_boot$X1_new <- data_boot$X1
    data_boot$X2_new <- data_boot$X2
    data_boot$X1_new[data_boot$id %in% a] <- jitter(data_boot$X1[data_boot$id %in% a])
    data_boot$X2_new[which(data_boot$X2==data_boot$X1)] <- data_boot$X1_new[which(data_boot$X2==data_boot$X1)]
    a1 <- which(data_boot$X2!=data_boot$X1)
    data_boot$X2_new[a1] <- data_boot$X2[a1] + data_boot$X1_new[a1] -  data_boot$X1[a1]
    data_boot$X1 <- data_boot$X1_new
    data_boot$X2 <- data_boot$X2_new
    data_boot$X1 <- data_boot$X1*10^4
    data_boot$X2 <- data_boot$X2*10^4
    data_boot$X2.1 <- ifelse((1-data_boot$delta1) * data_boot$delta2 == 1, data_boot$X2, data_boot$X1)

    ## recalculate the propensity score
    ### get the propensity score
    PS_OUT1 <- doPS(data = data_boot,
                    Trt = "A",
                    Trt.name = 1,
                    VARS. = VARS.,
                    logistic = TRUE)
    w <- PS_OUT1$Data$ipw_ate_stab
    data_boot$w <- PS_OUT1$Data$ipw_ate_stab

    EM_initial <- OUT_em_weights(data = data_boot,
                                 X1 = "X1",
                                 X2 = "X2",
                                 event1 = "delta1",
                                 event2 = "delta2",
                                 w = w,
                                 Trt = "A")

    res_boot <- try(em_illness_death_phmm_weight(data = data_boot,
                                                 X1 = "X1",
                                                 X2 = "X2",
                                                 event1 = "delta1",
                                                 event2 = "delta2",
                                                 w = w,
                                                 Trt = "A",
                                                 EM_initial = EM_initial,sigma_2_0 = sigma_2_0),TRUE)

    if(inherits(res_boot,"try-error") | res_boot$em.n==1){
      beta1_boot[b] <- NA
      beta2_boot[b] <- NA
      beta3_boot[b] <- NA
      sigma_2_boot[b] <- NA
      Lam01_1_boot[b] <- NA
      Lam02_1_boot[b] <- NA
      Lam03_1_boot[b] <- NA

    }else{

      beta1_boot[b] <- res_boot$beta1[res_boot$em.n]
      beta2_boot[b] <- res_boot$beta2[res_boot$em.n]
      beta3_boot[b] <- res_boot$beta3[res_boot$em.n]
      sigma_2_boot[b] <- res_boot$sigma_2[res_boot$em.n]

      Lam01_1_boot[b] <- res_boot$Lambda01(10^4)
      Lam02_1_boot[b] <- res_boot$Lambda02(10^4)
      Lam03_1_boot[b] <- res_boot$Lambda03(10^4)

    }
  }

  return(list(beta1_boot = beta1_boot,beta2_boot = beta2_boot, beta3_boot = beta3_boot,
              Lam01_1_boot = Lam01_1_boot,Lam02_1_boot = Lam02_1_boot,Lam03_1_boot = Lam03_1_boot,
              sigma_2_boot = sigma_2_boot))


}

# ##### Get the CIF function:
# cif_est_general <- function(res1,t1_star = t1_star){
#
#   if(res1$em.n==1){
#     a1 <- a2 <- a3 <- b1 <- b2 <- b3 <- Inf
#   }else{
#
#     ### Estimated the basline cumulative hazard of each event at each time
#     ## \Lambda_j(t;a) = \Lambda_0j(t)*exp(\beta * a)
#     bashaz.1 <- list()
#     for (i in 1:3) {
#       coef.1 <- res1[[i]][res1$em.n]
#       a <- res1[[i+3]]$cum_base_haz
#       a <- a[,c("time","cum_haz")]
#       a$hazard.A <- a$cum_haz
#       a$hazard.B <- a$cum_haz*exp(coef.1)
#       bashaz.1[[i]] <- a
#     }
#
#     ### Estimated Survival function:
#     ## S(t;a) = exp^(-(\Lambda_1(t;a) + \Lambda_1(t;a)))
#     bashaz.2 <- bashaz.1[c(1,2)]
#     S.all <- matrix(NA,nrow = nrow(bashaz.2[[1]]),ncol = 3)
#     S.all <- as.data.frame(S.all)
#     colnames(S.all) <- c("time","S.all.A","S.all.B")
#     S.all$time <- bashaz.2[[1]]$time
#     S.all$S.all.A <- exp( - (Reduce("+",bashaz.2)$hazard.A) )
#     S.all$S.all.B <- exp( - (Reduce("+",bashaz.2)$hazard.B) )
#
#
#     ### Calculate the CIF: F_1(t;a) = sum( S(t;a)*\delta(Lambda_1(t;a)) )
#     cif.1 <- matrix(NA,nrow = nrow(S.all),ncol = 3)
#     cif.1 <- as.data.frame(cif.1)
#     colnames(cif.1) <- c("time","cif.A","cif.B")
#     cif.1$time <- S.all$time
#     bashaz.2 <- bashaz.1[[1]]
#     for (i in 1:nrow(cif.1)){
#       if(i==1){
#         cif.1$cif.A[i] <- S.all$S.all.A[i] * bashaz.2$hazard.A[i]
#         cif.1$cif.B[i] <- S.all$S.all.B[i] * bashaz.2$hazard.B[i]
#       }else{
#         cif.1$cif.A[i] <- S.all$S.all.A[i] * (bashaz.2$hazard.A[i] - bashaz.2$hazard.A[i-1])
#         cif.1$cif.B[i] <- S.all$S.all.B[i] * (bashaz.2$hazard.B[i] - bashaz.2$hazard.B[i-1])
#       }
#     }
#     cif.1$cif.A.1 <- cumsum(cif.1$cif.A)
#     cif.1$cif.B.1 <- cumsum(cif.1$cif.B)
#     cif.1 <- cif.1[,c(1,4,5)]
#     colnames(cif.1) <- c("time","cif.A","cif.B")
#
#     ### Calculate the CIF: F_2(t;a) = sum( S(t;a)*\delta(Lambda_2(t;a)) )
#     cif.2 <- matrix(NA,nrow = nrow(S.all),ncol = 3)
#     cif.2 <- as.data.frame(cif.2)
#     colnames(cif.2) <- c("time","cif.A","cif.B")
#     cif.2$time <- S.all$time
#     bashaz.2 <- bashaz.1[[2]]
#     for (i in 1:nrow(cif.2)){
#       if(i==1){
#         cif.2$cif.A[i] <- S.all$S.all.A[i] * bashaz.2$hazard.A[i]
#         cif.2$cif.B[i] <- S.all$S.all.B[i] * bashaz.2$hazard.B[i]
#       }else{
#         cif.2$cif.A[i] <- S.all$S.all.A[i] * (bashaz.2$hazard.A[i] - bashaz.2$hazard.A[i-1])
#         cif.2$cif.B[i] <- S.all$S.all.B[i] * (bashaz.2$hazard.B[i] - bashaz.2$hazard.B[i-1])
#       }
#     }
#     cif.2$cif.A.1 <- cumsum(cif.2$cif.A)
#     cif.2$cif.B.1 <- cumsum(cif.2$cif.B)
#     cif.2 <- cif.2[,c(1,4,5)]
#     colnames(cif.2) <- c("time","cif.A","cif.B")
#
#     ## median of T2 for the third model: 4438
#     ### Calculate the CIF: F_12(t2,4438;a) = 1- exp( sum( \delta( Lambda_12(t2;a) -  Lambda_12(4438;a)  ) )  )
#     bashaz.2 <- bashaz.1[[3]]
#     #t1_star <- median(bashaz.2$time)
#     #t1_star <- 2000
#     cif.3 <- matrix(NA,nrow = nrow(bashaz.2),ncol = 3)
#     cif.3 <- as.data.frame(cif.3)
#     colnames(cif.3) <- c("time","cif.A","cif.B")
#     cif.3$time <- bashaz.2$time
#     Lambda12_a <- stepfun(bashaz.2$time,c(0,bashaz.2$hazard.A))
#     Lambda12_b <- stepfun(bashaz.2$time,c(0,bashaz.2$hazard.B))
#     for (i in 1:nrow(cif.3)){
#       if(cif.3$time[i]<=t1_star){
#         cif.3$cif.A[i] <- 0
#         cif.3$cif.B[i] <- 0
#       }else{
#         cif.3$cif.A[i] <- 1- exp(  -(  bashaz.2$hazard.A[i] - Lambda12_a(t1_star)  )  )
#         cif.3$cif.B[i] <- 1- exp(  -(  bashaz.2$hazard.B[i] - Lambda12_b(t1_star)  )  )
#       }
#     }
#
#     ## use step function to estimate the cif of the event time
#     a1 <- stepfun(cif.1$time,c(0,cif.1$cif.A))
#     a2 <- stepfun(cif.2$time,c(0,cif.2$cif.A))
#     a3 <- stepfun(cif.3$time,c(0,cif.3$cif.A))
#
#     ## use step function to estimate the cif of the event time
#     b1 <- stepfun(cif.1$time,c(0,cif.1$cif.B))
#     b2 <- stepfun(cif.2$time,c(0,cif.2$cif.B))
#     b3 <- stepfun(cif.3$time,c(0,cif.3$cif.B))
#
#   }
#
#   return(list(cif1a=a1,cif2a=a2,cif3a=a3,
#               cif1b=b1,cif2b=b2,cif3b=b3))
#
# }
#
#
# plot_est_cif <- function(cif.data,
#                          color = color,
#                          linetype = linetype){
#   ggplot(cif.data, aes_(x = ~time, y = ~cif.control)) +
#     geom_line(aes(colour=color[1],linetype=linetype[1]),size=1.4)+
#     geom_line(aes_(x=~time, y=~cif.exposure,colour=color[2],linetype=linetype[2]),size=1.4)+
#     ylim(c(0,1))+
#     xlab("Time") +
#     ylab("Cumulative incidence function")+
#     scale_color_manual(name = "Group", values = color, labels = c("Control","Exposure"))+
#     scale_linetype_manual(name = "Group", values = linetype, labels = c("Control","Exposure"))+
#     theme(
#       plot.title =   element_text(size=28, face="bold",hjust = 0.5),
#       axis.title.x = element_text(size=28, face="bold"),
#       axis.title.y = element_text(size=28, face="bold"),
#       axis.text.x = element_text(face="bold",size=22),
#       axis.text.y = element_text(face="bold",size=22),
#       legend.title = element_text(size=22),
#       legend.text = element_text(size=22),
#       legend.position= c(0.8,0.8))
# }
#
