####### Generate the data ########
sim_cox_msm_semicmrsk <- function(beta1,beta2,beta3,sigma_2,
                                 alpha0,alpha1,alpha2,alpha3,
                                 n,Cens){


  Lambda01_inv <- function(a){
    if(a >= 0 & a <= (2-2*exp(-3))){
      out <- log(2/(2-a))
    }else{
      out <- 4 + (a-2)/ (2*exp(-3))
    }
    return(out)
  }

  Lambda01 <- function(t){
    if(t>=0 & t<=3){
      out <- 2 - 2*exp(-t)
    }else{
      out <- 2 - 2*exp(-3) + (t-3) * 2 * exp(-3)
    }
    return(out)
  }

  inv_logit <- function(x){
    exp(x)/( 1+exp(x) )
  }

  #### 2 uniform distributed r.v. U1 & U2
  U1 <- runif(n,0,1)
  U2 <- runif(n,0,1)

  #### confounder Z: Z1, Z2, Z3 that depends on U1 & U2
  Z1 <- U1 + U2 + rnorm(n,0,1)
  Z2 <- U1 + U2 + rnorm(n,0,1.5)
  Z3 <- U1 + U2 + rnorm(n,0,1.8)

  #### Treatment assignment A that depends on Z1,Z2,Z3 which follows a logistic regression
  #### A ~ Bernoulli(pA), where pA ~ logit^{-1}(\alpha_0 + \alpha_1*Z1 + \alpha_2*Z2 + \alpha_3*Z3)
  pA <- inv_logit(alpha0+alpha1*Z1+alpha2*Z2+alpha3*Z3)
  A <- NULL
  for (i in 1:n) {
    A[i] <- rbinom(1,1,prob = pA[i])
  }

  ### Frailty model
  b <- rnorm(n,0,sqrt(sigma_2))
  ### Calculate the P(T1 = +\infty)
  p_T1_infty <- exp(beta2*A + b) / ( exp(beta1*A + b) + exp(beta2*A + b) )
  ### Generate the indicator for whether T1=+\infty for this subject
  T1_inf_ind <- NULL
  for (i in 1:n) {
    T1_inf_ind[i] <- rbinom(1,1,p_T1_infty[i])
  }
  ### Generate the survival time
  T1 <- NULL
  T2 <- NULL
  for (i in 1:n) {
    #### T1 = infty situation, only generate T2 based on U1 ###
    if(T1_inf_ind[i]==1){
      T1[i] <- Inf
      a <- - log(U2[i]) / ( exp(beta1*A[i]+b[i]) +  exp(beta2*A[i]+b[i]) )
      T2[i] <- Lambda01_inv(a)
    }
    #### T1 < infty situation, generate T1 first, then generate T2 based on T1 ####
    if(T1_inf_ind[i]==0){
      T1[i] <- Lambda01_inv( - log(U1[i]) / ( exp(beta1*A[i]+b[i]) +  exp(beta2*A[i]+b[i]) ) )
      T2[i] <- Lambda01_inv( - ( log(U2[i]) / ( 2*exp(beta3*A[i]+b[i]) ) ) + Lambda01(T1[i])   )
    }
  }

  # ### get the propensity score
  # fit.1 <- glm(A ~ Z1 + Z2 + Z3,family = binomial(link = "logit"))
  # ps1 <- predict(fit.1,type = "response")
  # w <- A / ps1 + (1-A) / (1-ps1)

  data0 <- as.data.frame(cbind(c(1:n),A,Z1,Z2,Z3,b,pA,T1,T2,T1_inf_ind))
  colnames(data0)[1] <- "id"

  data0$C <- Cens
  data0$X2 <- pmin(data0$T2,data0$C)
  data0$X1 <- pmin(data0$T1,data0$X2)
  data0$delta1 <- ifelse(data0$T1==data0$X1,1,0)
  data0$delta2 <- ifelse(data0$T2==data0$X2,1,0)
  data <- data0

  OUT <- list(data0 = data0,rate=length(which(data0$delta1==0 & data0$delta2==0))/n)
  return(OUT)
}
