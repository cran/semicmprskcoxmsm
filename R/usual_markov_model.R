usual_illness_death_weight <- function(data,X1,X2,event1,event2,w,Trt){


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

  beta1_est <- fit1$coefficients[1]
  beta2_est <- fit2$coefficients[1]
  beta3_est <- fit3$coefficients[1]

  sd_beta1 <- sqrt(vcov(fit1))
  sd_beta2 <- sqrt(vcov(fit2))
  sd_beta3 <- sqrt(vcov(fit3))

  Lambda1 <- get_hazard(fit1)
  Lambda2 <- get_hazard(fit2)
  Lambda3 <- get_hazard(fit3)

  res1 <- list(beta1 = beta1_est,beta2 = beta2_est, beta3 = beta3_est,
               sd_beta1 = sd_beta1,sd_beta2 = sd_beta2,sd_beta3 = sd_beta3,
               Lambda01 = Lambda1,Lambda02 = Lambda2,Lambda03 = Lambda3)

  return(res1)

}

cif_est_usual <- function(data,X1,X2,event1,event2,w,Trt,
                           t1_star = t1_star){

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

  cumbasehaz1 <- basehaz(fit1,centered = FALSE)
  cumbasehaz2 <- basehaz(fit2,centered = FALSE)
  cumbasehaz3 <- basehaz(fit3,centered = FALSE)

  cumbasehaz1$hazard.A <- cumbasehaz1$hazard
  cumbasehaz2$hazard.A <- cumbasehaz2$hazard
  cumbasehaz3$hazard.A <- cumbasehaz3$hazard

  cumbasehaz1$hazard.B <- cumbasehaz1$hazard * exp(fit1$coefficients[1])
  cumbasehaz2$hazard.B <- cumbasehaz2$hazard * exp(fit2$coefficients[1])
  cumbasehaz3$hazard.B <- cumbasehaz3$hazard * exp(fit3$coefficients[1])

  bashaz.1 <- list(cumbasehaz1,cumbasehaz2,cumbasehaz3)

  ### Estimated Survival function:
  ## S(t;a) = exp^(-(\Lambda_1(t;a) + \Lambda_1(t;a)))
  bashaz.2 <- bashaz.1[c(1,2)]
  S.all <- matrix(NA,nrow = nrow(bashaz.2[[1]]),ncol = 3)
  S.all <- as.data.frame(S.all)
  colnames(S.all) <- c("time","S.all.A","S.all.B")
  S.all$time <- bashaz.2[[1]]$time
  S.all$S.all.A <- exp( - (Reduce("+",bashaz.2)$hazard.A) )
  S.all$S.all.B <- exp( - (Reduce("+",bashaz.2)$hazard.B) )

  ### Calculate the cif_u: F_1(t;a) = sum( S(t;a)*\delta(Lambda_1(t;a)) )
  cif_u.1 <- matrix(NA,nrow = nrow(S.all),ncol = 3)
  cif_u.1 <- as.data.frame(cif_u.1)
  colnames(cif_u.1) <- c("time","cif_u.A","cif_u.B")
  cif_u.1$time <- S.all$time
  bashaz.2 <- bashaz.1[[1]]
  for (i in 1:nrow(cif_u.1)){
    if(i==1){
      cif_u.1$cif_u.A[i] <- S.all$S.all.A[i] * bashaz.2$hazard.A[i]
      cif_u.1$cif_u.B[i] <- S.all$S.all.B[i] * bashaz.2$hazard.B[i]
    }else{
      cif_u.1$cif_u.A[i] <- S.all$S.all.A[i] * (bashaz.2$hazard.A[i] - bashaz.2$hazard.A[i-1])
      cif_u.1$cif_u.B[i] <- S.all$S.all.B[i] * (bashaz.2$hazard.B[i] - bashaz.2$hazard.B[i-1])
    }
  }
  cif_u.1$cif_u.A.1 <- cumsum(cif_u.1$cif_u.A)
  cif_u.1$cif_u.B.1 <- cumsum(cif_u.1$cif_u.B)
  cif_u.1 <- cif_u.1[,c(1,4,5)]
  colnames(cif_u.1) <- c("time","cif_u.A","cif_u.B")

  ### Calculate the cif_u: F_2(t;a) = sum( S(t;a)*\delta(Lambda_2(t;a)) )
  cif_u.2 <- matrix(NA,nrow = nrow(S.all),ncol = 3)
  cif_u.2 <- as.data.frame(cif_u.2)
  colnames(cif_u.2) <- c("time","cif_u.A","cif_u.B")
  cif_u.2$time <- S.all$time
  bashaz.2 <- bashaz.1[[2]]
  for (i in 1:nrow(cif_u.2)){
    if(i==1){
      cif_u.2$cif_u.A[i] <- S.all$S.all.A[i] * bashaz.2$hazard.A[i]
      cif_u.2$cif_u.B[i] <- S.all$S.all.B[i] * bashaz.2$hazard.B[i]
    }else{
      cif_u.2$cif_u.A[i] <- S.all$S.all.A[i] * (bashaz.2$hazard.A[i] - bashaz.2$hazard.A[i-1])
      cif_u.2$cif_u.B[i] <- S.all$S.all.B[i] * (bashaz.2$hazard.B[i] - bashaz.2$hazard.B[i-1])
    }
  }
  cif_u.2$cif_u.A.1 <- cumsum(cif_u.2$cif_u.A)
  cif_u.2$cif_u.B.1 <- cumsum(cif_u.2$cif_u.B)
  cif_u.2 <- cif_u.2[,c(1,4,5)]
  colnames(cif_u.2) <- c("time","cif_u.A","cif_u.B")

  ## median of T2 for the third model: 4438
  ### Calculate the cif_u: F_12(t2,4438;a) = 1- exp( sum( \delta( Lambda_12(t2;a) -  Lambda_12(4438;a)  ) )  )
  bashaz.2 <- bashaz.1[[3]]
  # t1_star <- median(bashaz.2$time)
  # t1_star <- 2000
  cif_u.3 <- matrix(NA,nrow = nrow(bashaz.2),ncol = 3)
  cif_u.3 <- as.data.frame(cif_u.3)
  colnames(cif_u.3) <- c("time","cif_u.A","cif_u.B")
  cif_u.3$time <- bashaz.2$time
  Lambda12_a <- stepfun(bashaz.2$time,c(0,bashaz.2$hazard.A))
  Lambda12_b <- stepfun(bashaz.2$time,c(0,bashaz.2$hazard.B))
  for (i in 1:nrow(cif_u.3)){
    if(cif_u.3$time[i]<=t1_star){
      cif_u.3$cif_u.A[i] <- 0
      cif_u.3$cif_u.B[i] <- 0
    }else{
      cif_u.3$cif_u.A[i] <- 1- exp(  -(  bashaz.2$hazard.A[i] - Lambda12_a(t1_star)  )  )
      cif_u.3$cif_u.B[i] <- 1- exp(  -(  bashaz.2$hazard.B[i] - Lambda12_b(t1_star)  )  )
    }
  }

  ## use step function to estimate the cif of the event time
  a1 <- stepfun(cif_u.1$time,c(0,cif_u.1$cif_u.A))
  a2 <- stepfun(cif_u.2$time,c(0,cif_u.2$cif_u.A))
  a3 <- stepfun(cif_u.3$time,c(0,cif_u.3$cif_u.A))

  ## use step function to estimate the cif of the event time
  b1 <- stepfun(cif_u.1$time,c(0,cif_u.1$cif_u.B))
  b2 <- stepfun(cif_u.2$time,c(0,cif_u.2$cif_u.B))
  b3 <- stepfun(cif_u.3$time,c(0,cif_u.3$cif_u.B))

  return(list(cif1a=a1,cif2a=a2,cif3a=a3,
              cif1b=b1,cif2b=b2,cif3b=b3))

}
