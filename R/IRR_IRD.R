conditional_cif_b <- function(res1,
                              t1_star,
                              b){

  if(res1$em.n==1){
    a1 <- a2 <- a3 <- b1 <- b2 <- b3 <- Inf
  }else{

    ### Estimated the basline cumulative hazard of each event at each time
    ## \Lambda_j(t;a) = \Lambda_0j(t)*exp(\beta * a)
    bashaz.1 <- list()
    for (i in 1:3) {
      coef.1 <- res1[[i]][res1$em.n]
      a <- res1[[i+3]]$cum_base_haz
      a <- a[,c("time","cum_haz")]
      a$hazard.A <- a$cum_haz*exp(b)
      a$hazard.B <- a$cum_haz*exp(coef.1+b)
      bashaz.1[[i]] <- a
    }

    ### Estimated Survival function:
    ## S(t;a) = exp^(-(\Lambda_1(t;a) + \Lambda_2(t;a)))
    bashaz.2 <- bashaz.1[c(1,2)]
    S.all <- matrix(NA,nrow = nrow(bashaz.2[[1]]),ncol = 3)
    S.all <- as.data.frame(S.all)
    colnames(S.all) <- c("time","S.all.A","S.all.B")
    S.all$time <- bashaz.2[[1]]$time
    S.all$S.all.A <- exp( - (Reduce("+",bashaz.2)$hazard.A) )
    S.all$S.all.B <- exp( - (Reduce("+",bashaz.2)$hazard.B) )

    ### Calculate the CIF: F_1(t;a) = sum( S(t;a)*\delta(Lambda_1(t;a)) )
    cif.1 <- matrix(NA,nrow = nrow(S.all),ncol = 3)
    cif.1 <- as.data.frame(cif.1)
    colnames(cif.1) <- c("time","cif.A","cif.B")
    cif.1$time <- S.all$time
    bashaz.2 <- bashaz.1[[1]]
    cif.1$cif.A <- S.all$S.all.A * c(0,diff(bashaz.2$hazard.A))
    cif.1$cif.B <- S.all$S.all.B * c(0,diff(bashaz.2$hazard.B))
    cif.1$cif.A.1 <- cumsum(cif.1$cif.A)
    cif.1$cif.B.1 <- cumsum(cif.1$cif.B)
    cif.1 <- cif.1[,c(1,4,5)]
    colnames(cif.1) <- c("time","cif.A","cif.B")

    ### Calculate the CIF: F_2(t;a) = sum( S(t;a)*\delta(Lambda_2(t;a)) )
    cif.2 <- matrix(NA,nrow = nrow(S.all),ncol = 3)
    cif.2 <- as.data.frame(cif.2)
    colnames(cif.2) <- c("time","cif.A","cif.B")
    cif.2$time <- S.all$time
    bashaz.2 <- bashaz.1[[2]]
    cif.2$cif.A <- S.all$S.all.A * c(0,diff(bashaz.2$hazard.A))
    cif.2$cif.B <- S.all$S.all.B * c(0,diff(bashaz.2$hazard.B))
    cif.2$cif.A.1 <- cumsum(cif.2$cif.A)
    cif.2$cif.B.1 <- cumsum(cif.2$cif.B)
    cif.2 <- cif.2[,c(1,4,5)]
    colnames(cif.2) <- c("time","cif.A","cif.B")

    ## median of T2 for the third model: 4438
    ### Calculate the CIF: F_12(t2,4438;a) = 1- exp( sum( \delta( Lambda_12(t2;a) -  Lambda_12(4438;a)  ) )  )
    bashaz.2 <- bashaz.1[[3]]
    #t1_star <- median(bashaz.2$time)
    #t1_star <- 2000
    cif.3 <- matrix(NA,nrow = nrow(bashaz.2),ncol = 3)
    cif.3 <- as.data.frame(cif.3)
    colnames(cif.3) <- c("time","cif.A","cif.B")
    cif.3$time <- bashaz.2$time
    Lambda12_a_t1_star <- stepfun(bashaz.2$time,c(0,bashaz.2$hazard.A))(t1_star)
    Lambda12_b_t1_star <- stepfun(bashaz.2$time,c(0,bashaz.2$hazard.B))(t1_star)
    cif.3$cif.A <- ifelse(cif.3$time<=t1_star, 0 , 1- exp(  - (  bashaz.2$hazard.A - Lambda12_a_t1_star  )  ) )
    cif.3$cif.B <- ifelse(cif.3$time<=t1_star, 0 , 1- exp(  - (  bashaz.2$hazard.B - Lambda12_b_t1_star  )  ) )

    ## use step function to estimate the cif of the event time
    a1 <- stepfun(cif.1$time,c(0,cif.1$cif.A))
    a2 <- stepfun(cif.2$time,c(0,cif.2$cif.A))
    a3 <- stepfun(cif.3$time,c(0,cif.3$cif.A))

    ## use step function to estimate the cif of the event time
    b1 <- stepfun(cif.1$time,c(0,cif.1$cif.B))
    b2 <- stepfun(cif.2$time,c(0,cif.2$cif.B))
    b3 <- stepfun(cif.3$time,c(0,cif.3$cif.B))

  }

  return(list(a1=a1,a2=a2,a3 = a3,
              b1=b1,b2=b2,b3 = b3,
              cif.1 = cif.1, cif.2 = cif.2, cif.3 = cif.3))

}

individual_RR_RD <- function(dat1,
                             res1,
                             t1_star = t1_star,
                             t){

  dat1$RR1 <- NA
  dat1$RR2 <- NA
  dat1$RR3 <- NA

  dat1$RD1 <- NA
  dat1$RD2 <- NA
  dat1$RD3 <- NA

  n <- nrow(dat1)

  for (i in 1:n) {
    a <- conditional_cif_b( res1,t1_star,b = dat1$b_hat[i] )
    dat1$RR1[i] <- a$b1(t) /  a$a1(t)
    dat1$RR2[i] <- a$b2(t) /  a$a2(t)
    dat1$RR3[i] <- a$b3(t) /  a$a3(t)

    dat1$RD1[i] <- a$b1(t) -  a$a1(t)
    dat1$RD2[i] <- a$b2(t) -  a$a2(t)
    dat1$RD3[i] <- a$b3(t) -  a$a3(t)

  }

  return(dat1)

}

bayesian_boot_irrd <- function(dat2,
                               B,
                               sigma_2_0,
                               EM_initial,
                               varlist,
                               t1_star,t){
  n <- nrow(dat2)
  seed.boot <- sample(1:100000,B)
  RD1_boot <- matrix(NA,nrow = n,ncol = B)
  RD2_boot <- matrix(NA,nrow = n,ncol = B)
  RD3_boot <- matrix(NA,nrow = n,ncol = B)
  RR1_boot <- matrix(NA,nrow = n,ncol = B)
  RR2_boot <- matrix(NA,nrow = n,ncol = B)
  RR3_boot <- matrix(NA,nrow = n,ncol = B)


  for (b in 1:B) {
    print(paste("Bootstrap:",b))
    set.seed(seed.boot[b])
    u <- rexp(n)
    w <- u / mean(u)

    ## re-calculate the weights
    ps1 <- doPS(data = dat2,
                       Trt = "A",
                       Trt.name = 1,
                       VARS. = varlist,
                       logistic = FALSE, w = w)
    w_ate <- ps1$Data$ipw_ate_stab
    w <- w * w_ate

    res_boot <- em_illness_death_phmm_weight(data = dat2,
                                             EM_initial = EM_initial,
                                             sigma_2_0 = sigma_2_0, w = w)

    dat2_boot <- res_boot$data
    dat0 <- individual_RR_RD(dat1 = dat2_boot,
                             res1 = res_boot,
                             t1_star = t1_star,t = t)

    RD1_boot[,b] <- dat0$RD1
    RD2_boot[,b] <- dat0$RD2
    RD3_boot[,b] <- dat0$RD3

    RR1_boot[,b] <- dat0$RR1
    RR2_boot[,b] <- dat0$RR2
    RR3_boot[,b] <- dat0$RR3

  }

  return(RD1_boot = RD1_boot, RD2_boot = RD2_boot, RD3_boot = RD3_boot,
         RR1_boot = RR1_boot, RR2_boot = RR2_boot, RR3_boot = RR3_boot)

}




