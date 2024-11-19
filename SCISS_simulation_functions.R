#all functions used in the SCISS
#regular setting----
Generate_y_x <- function(Graph,p,n,N,c,seed,type,hd=FALSE)
{
  q <- length(Graph[1,])
  #randomness
  set.seed(seed)
  supervise_num <- sort(sample.int(n+N,n))
  
  #generate y
  G <- Graph
  diag(G) <- 0
  mu <- diag(Graph)
  Data <- IsingSampler(n+N,G,thresholds=mu,method="CFTP")
  S_y <- Data
  S_y_supervise <- S_y[supervise_num,]
  
  #generate x
  S_x <- array(1,c(n+N,p+1))
  if(type=="gaussian"){
    for(i in 1:(n+N))
    {
      noise <- rnorm(p,0,1)
      for(j in 1:q)
      {
        S_x[i,j+1] <- S_y[i,j]*c[j,j]
        for(k in 1:q)
        {
          if(k!=j)
          {
            S_x[i,j+1] <- S_x[i,j+1]+c[j,k]*S_y[i,j]*S_y[i,k]
          }
        }
      }
      S_x[i,2:(p+1)] <- S_x[i,2:(p+1)]+noise
    }
  }else if(type=="poisson"){
    for(i in 1:(n+N))
    {
      for(j in 1:q)
      {
        lambda <- S_y[i,j]*c[j,j]
        for(k in 1:q)
        {
          if(k!=j)
          {
            lambda <- lambda+c[j,k]*S_y[i,j]*S_y[i,k]
          }
        }
        S_x[i,j+1] <- rpois(1,lambda)
      }
    }
  }

  if(hd==TRUE){
    S_x_attach <- sample(c(0, 2, 4), size = (n+N)*300*p,
                         replace = TRUE, prob = c(0.7, 0.1, 0.2))
    
    S_x_attach <- array(S_x_attach,c(n+N,300*p))

    S_x <- cbind(S_x,S_x_attach)
  }
  
  S_x_supervise <- S_x[supervise_num,]
  S_x_unsupervise <- S_x[-supervise_num,]
  
  List <- list(S_y_supervise=S_y_supervise,S_x_supervise=S_x_supervise,
               S_x_unsupervise=S_x_unsupervise)
  return(List)
}

A_y_generation <- function(q)
{
  A_y <- array(0,c(2^q,q))
  col_sum <- 1
  for(i in 1:q){
    combination <- combn(q, i)
    for(t in 1:ncol(combination))
    {
      for(k in 1:i){
        y_order <- combination[k,t]
        A_y[col_sum+t,y_order] <- 1
      }
    }
    col_sum <- col_sum+ncol(combination)
  }
  return(A_y)
}


g <- function(a)
{
  out <- exp(a)/(1+exp(a))
  return(out)
}

dg <- function(a)
{
  out <- exp(a)/((1+exp(a))^2)
  return(out)
}

Hessian_est <- function(theta_j,S_y,j)
{
  n <- length(S_y[,1])
  q <- length(S_y[1,])
  
  Hessian_est_single <- array(0,dim <- c(q,q))
  for(i in 1:n) 
  {y_setminus_j <- S_y[i,]
  y_setminus_j[j] <- 1
  Hessian_est_single <- 
    Hessian_est_single +(1/n)*y_setminus_j%*%t(y_setminus_j)*as.numeric(dg(theta_j%*%y_setminus_j))}
  return(Hessian_est_single)
}

S_estimation <- function(inv_Hessian,y1,theta_supervise,j)
{
  q <- length(y1[1,])
  n <- length(y1[,1])
  s <- array(0,c(n,q))
  for(i in 1:n){
    y_s <- y1[i,]
    y_s[j] <- 1
    y_j <- y1[i,j]
    s[i,] <- inv_Hessian%*%y_s*as.numeric(y_j-g(theta_supervise[,j]%*%y_s))
  }
  return(s)
}


SL <- function(y1)
{
  n <- length(y1[,1])
  q <- length(y1[1,])
  theta_supervise <- array(0,c(q,q))
  
  for (j in 1:q){
    y_j1 <- y1
    y_j1[,j] <- 1
    fit <- glm.fit(x = model.matrix(~ 0+y_j1),
                   y = y1[,j],
                   family = binomial(link = "logit"))
    theta_supervise[,j] <- fit$coefficients
  }
  theta_supervise <- (theta_supervise+t(theta_supervise))/2
  
  
  s <- array(0,c(n,q,q))
  for(j in 1:q){
    inv_Hessian_est <- solve(Hessian_est(theta_supervise[,j],y1,j))
    s[,j,] <- S_estimation(inv_Hessian_est,y1,theta_supervise,j)
  }
  
  L<-list(s=s,theta_SL=theta_supervise)
  return(L)
}


m_estimation <- function(inv_Hessian,A_y,P_y_line_z_eta,M,j,theta_supervise)
{
  q <- length(A_y[1,])
  m <- array(0,c(M,q))
  for(i in 1:M){
    S_condition_j_0 <- rep(0,q)
    for(a in 1:2^q){y_s <- A_y[a,]
    y_s[j] <- 1
    y_j <- A_y[a,j]
    S_specified_j_0 <- inv_Hessian%*%y_s*as.numeric(y_j-g(t(theta_supervise[,j])%*%y_s))
    S_condition_j_0 <- S_condition_j_0+S_specified_j_0*as.numeric(P_y_line_z_eta[i,a])
    }
    m[i,] <- S_condition_j_0
  }
  return(m)
}


P_y_line_z_eta_cal <- function(y1,x1,x2,hd=FALSE,theta_supervise,type_conmodel,PoS_model="gaussian")
{
  q <- length(y1[1,])
  p <- length(x1[1,])-1
  n <- length(x1[,1])
  N <- length(x2[,1])
  
  if(type_conmodel=="Aug"){
    gamma_supervise <- array(0,c(q,q*(p+1)))
    psi_j <- array(0,c(n,q*(p+1)))
    psi_sj <- rep(0,(p+1)*q)
    for(j in 1:q){
      for(i in 1:n){
        for(t in 1:q){
          if(t==j){
            psi_sj[((t-1)*(p+1)+1):(t*(p+1))] <- x1[i,]
          }
          else{psi_sj[((t-1)*(p+1)+1):(t*(p+1))] <- x1[i,]*y1[i,t]}
          psi_j[i,]<-psi_sj
        }}
      ifelse(hd==TRUE,fit <- glmnet(psi_j,
                                       as.matrix(y1[,j]),alpha=1,family="binomial",lambda=n^{-3/4},intercept=FALSE),
             fit <- glmnet(psi_j,
                              y1[,j],alpha=0,family="binomial",lambda=n^{-3/4},intercept=FALSE))
      gamma_supervise[j,] <- coef(fit)[2:(q*(p+1)+1)] 
    }
    
    P_all <- array(0,n+N)
    P_a_log <- matrix(0,n+N,2^q)
    P_a <- matrix(0,n+N,2^q)
    a <- 1:2^q
    for(i in 1:n){
      for(j in 1:q){
        P_a_log[i,a] <- P_a_log[i,a]+as.numeric(gamma_supervise[j,((j-1)*(p+1)+1):(j*(p+1))]%*%
                                                  x1[i,])*A_y[a,j]
        if(j!=q){
          for(j1 in (j+1):q)
          {P_a_log[i,a] <- P_a_log[i,a]+as.numeric(gamma_supervise[j,((j1-1)*(p+1)+1):(j1*(p+1))]%*%
                                                     x1[i,])*A_y[a,j]*A_y[a,j1]
          }}
      }
      P_all[i] <- sum(exp(P_a_log[i,a]))
      P_a[i,a] <- exp(P_a_log[i,a])/P_all[i]
    }
    
    for(i in 1:N){
      for(j in 1:q){
        P_a_log[i+n,a] <- P_a_log[i+n,a]+as.numeric(gamma_supervise[j,((j-1)*(p+1)+1):(j*(p+1))]%*%
                                                      x2[i,])*A_y[a,j]
        if(j!=q){
          for(j1 in (j+1):q)
          {P_a_log[i+n,a] <- P_a_log[i+n,a]+as.numeric(gamma_supervise[j,((j1-1)*(p+1)+1):(j1*(p+1))]%*%
                                                         x2[i,])*A_y[a,j]*A_y[a,j1]
          }}
      }
      P_all[i+n] <- sum(exp(P_a_log[i+n,a]))
      P_a[i+n,a]<- exp(P_a_log[i+n,a])/P_all[i+n]
    }
    return(P_a)
  }else if(type_conmodel=="PoS"){
    
    if(hd==TRUE){
      beta <- array(0,c(q,2))
      significance_x_y <- array(0,c(q,length(x1[1,])-1))
      for(j in 1:q)
      {
        Y_j <- y1[,j]
        significance_all <- rep(0,p)
        for(k in 1:(length(x1[1,])-1)){
          if(PoS_model=="gaussian")
          {fit <- glm(x1[,k+1]~Y_j,family="gaussian")
          significance_all[k] <- summary(fit)$coefficients[2,4]}
          if(PoS_model=="binomial")
          {fit <- glm(x1[,k+1]~Y_j,family='binomial')
          significance_all[k] <- summary(fit)$coefficients[2,4]}
          if(PoS_model=="poisson")
          {fit <- glm(x1[,k+1]~Y_j,family='poisson')
          significance_all[k] <- summary(fit)$coefficients[2,4]}
        }
        significance_x_y[j,] <- significance_all
      }
      
      used_significance_index <- c() 
      result_index <- c() 
      
      for (j in 1:q) {
        row_order <- order(significance_x_y[j, ])
        for (k in row_order) {
          if (!(k %in% result_index)) {
            result_index <- c(result_index, k)
            used_significance_index <- c(used_significance_index, k)
            break
          }
        }
      }
      
      x1_significance <- x1[,c(1,result_index+1)]
      x2_significance <- x2[,c(1,result_index+1)]
      
      for(j in 1:q){
        Y_j <- y1[,j]
        if(PoS_model=="gaussian")
        {fit <- glm(x1_significance[,j+1]~Y_j,family="gaussian")}
        if(PoS_model=="binomial")
        {fit <- glm(x1_significance[,j+1]~Y_j,family='binomial')}
        if(PoS_model=="poisson")
        {fit <- glm(x1_significance[,j+1]~Y_j,family='poisson')}
        
        beta[j,]<-coef(fit)
      }
      
      #calculate P_ylinew
      P_ylinew <- array(0,2^q)
      P_ylinew_all <- 0
      a <- 1:2^q
      for(j1 in 1:q)
      {
        P_ylinew[a] <- P_ylinew[a]+theta_supervise[j1,j1]*A_y[a,j1]
        if(j1 < q)
        {
          for(j2 in (j1+1):q)
          {
            P_ylinew[a] <- P_ylinew[a]+theta_supervise[j1,j2]*A_y[a,j1]*A_y[a,j2]
          }
        }
      }
      P_ylinew_all <- sum(exp(P_ylinew))
      
      #calculate P(x,y|w)
      P_xylinew <- matrix(0,n+N,2^q)
      for(i in 1:(n+N)){
        if(i<=n){
          for(a in 1:2^q)
          {P_xylinew[i,a] <- exp(P_ylinew[a])/P_ylinew_all
          for(j1 in 1:q)
          {Y_j <- c(1,1)
          Y_j[2] <- A_y[a,j1]
          if(PoS_model=="gaussian")
          {P_xylinew[i,a] <- P_xylinew[i,a]*dnorm(x1_significance[i,j1+1],Y_j%*%beta[j1,],0.1)}
          if(PoS_model=="binomial")
          {P_xylinew[i,a] <- P_xylinew[i,a]*dbinom(x1_significance[i,j1+1],1,prob=exp(Y_j%*%beta[j1,])/(1+exp(Y_j%*%beta[j1,])) )}
          if(PoS_model=="poisson")
          {P_xylinew[i,a] <- P_xylinew[i,a]*dpois(x1_significance[i,j1+1],exp(Y_j%*%beta[j1,]),log=FALSE)}
          }
          }
        }else{
          for(a in 1:2^q)
          {P_xylinew[i,a] <- exp(P_ylinew[a])/P_ylinew_all
          
          for(j1 in 1:q)
          {Y_j<-c(1,1)
          Y_j[2]<-A_y[a,j1]
          if(PoS_model=="gaussian")
          {P_xylinew[i,a]<-P_xylinew[i,a]*dnorm(x2_significance[i-n,j1+1],Y_j%*%beta[j1,],0.1)}
          if(PoS_model=="binomial")
          {P_xylinew[i,a]<-P_xylinew[i,a]*dbinom(x2_significance[i-n,j1+1],1,prob=exp(Y_j%*%beta[j1,])/(1+exp(Y_j%*%beta[j1,])) )}
          if(PoS_model=="poisson")
          {P_xylinew[i,a]<-P_xylinew[i,a]*dpois(x2_significance[i-n,j1+1],exp(Y_j%*%beta[j1,]),log=FALSE)}
          }
          }
        }
      }    
      
      #calculate P(y|z(x,w))
      P_ylinexw <- matrix(0,n+N,2^q)
      for(i in 1:(n+N)){
        if(sum(P_xylinew[i,]==array(0,2^q))==2^q)
        {P_ylinexw[i,] <- array(1/2^q,2^q)}else{
          for(a in 1:2^q)
          {P_ylinexw[i,a] <- P_xylinew[i,a]/sum(P_xylinew[i,])}
        }
      }
      return(P_ylinexw)
    }else{
      beta<-array(0,c(q,2))
      for(j in 1:q)
      {
        Y_j <- y1[,j]
        if(PoS_model=="gaussian")
        {fit <- glm(x1[,j+1]~Y_j,family="gaussian")}
        if(PoS_model=="binomial")
        {fit <- glm(x1[,j+1]~Y_j,family='binomial')}
        if(PoS_model=="poisson")
        {fit <- glm(x1[,j+1]~Y_j,family='poisson')}
        beta[j,]<-coef(fit)
      }
      
      #calculate P_ylinew
      P_ylinew <- array(0,2^q)
      P_ylinew_all <- 0
      a <- 1:2^q
      for(j1 in 1:q)
      {
        P_ylinew[a] <- P_ylinew[a]+theta_supervise[j1,j1]*A_y[a,j1]
        if(j1<q)
        {
          for(j2 in (j1+1):q)
          {
            P_ylinew[a] <- P_ylinew[a]+theta_supervise[j1,j2]*A_y[a,j1]*A_y[a,j2]
          }
        }
      }
      P_ylinew_all <- sum(exp(P_ylinew))
      
      #calculate P(x,y|w)
      P_xylinew <- matrix(0,n+N,2^q)
      for(i in 1:(n+N)){
        if(i<=n){
          for(a in 1:2^q)
          {P_xylinew[i,a]<-exp(P_ylinew[a])/P_ylinew_all
          for(j1 in 1:q)
          {Y_j<-c(1,1)
          Y_j[2] <- A_y[a,j1]
          if(PoS_model=="gaussian")
          {P_xylinew[i,a] <- P_xylinew[i,a]*dnorm(x1[i,j1+1],Y_j%*%beta[j1,],0.1)}
          if(PoS_model=="binomial")
          {P_xylinew[i,a] <- P_xylinew[i,a]*dbinom(x1[i,j1+1],1,prob=exp(Y_j%*%beta[j1,])/(1+exp(Y_j%*%beta[j1,])) )}
          if(PoS_model=="poisson")
          {P_xylinew[i,a] <- P_xylinew[i,a]*dpois(x1[i,j1+1],exp(Y_j%*%beta[j1,]),log=FALSE)}
          }
          }
        }else{
          for(a in 1:2^q)
          {P_xylinew[i,a] <- exp(P_ylinew[a])/P_ylinew_all
          
          for(j1 in 1:q)
          {Y_j<-c(1,1)
          Y_j[2]<-A_y[a,j1]
          if(PoS_model=="gaussian")
          {P_xylinew[i,a]<-P_xylinew[i,a]*dnorm(x2[i-n,j1+1],Y_j%*%beta[j1,],0.1)}
          if(PoS_model=="binomial")
          {P_xylinew[i,a]<-P_xylinew[i,a]*dbinom(x2[i-n,j1+1],1,prob=exp(Y_j%*%beta[j1,])/(1+exp(Y_j%*%beta[j1,])) )}
          if(PoS_model=="poisson")
          {P_xylinew[i,a]<-P_xylinew[i,a]*dpois(x2[i-n,j1+1],exp(Y_j%*%beta[j1,]),log=FALSE)}
          }
          }
        }
      }    
      
      #calculate P(y|z(x,w))
      P_ylinexw <- matrix(0,n+N,2^q)
      for(i in 1:(n+N)){
        if(sum(P_xylinew[i,]==array(0,2^q))==2^q)
        {P_ylinexw[i,]=array(1/2^q,2^q)}else{
          for(a in 1:2^q)
          {P_ylinexw[i,a]<-P_xylinew[i,a]/sum(P_xylinew[i,])}
        }
      }
    }
    return(P_ylinexw)
    }
    #end of PoS

}


SCISS <- function(y1,x1,x2,hd_Aug=FALSE,hd_PoS=FALSE,PoS_model="gaussian")
{
  n <- length(y1[,1])
  N <- length(x2[,1])
  q <- length(y1[1,])
  p <- length(x1[1,])-1
  
  T_S <- SL(y1)
  theta_supervise <- T_S$theta_SL
  s <- T_S$s
  
  if(hd_Aug==FALSE){
    P_a_Aug <- P_y_line_z_eta_cal(y1,x1,x2,theta_supervise,type_conmodel = "Aug")
    m_Aug <- array(0,c(n+N,q,q))
    theta_ss_Aug <- array(0,c(q,q))
    
    for(j in 1:q){
      Hessian <- Hessian_est(theta_supervise[,j],y1,j)
      inv_Hessian <- solve(Hessian)
      
      m_Aug[,j,] <- m_estimation(inv_Hessian,A_y,P_a_Aug,n+N,j,theta_supervise)
      
      theta_ss_Aug[,j] <- theta_supervise[,j]-n^{-1}*colSums(m_Aug[1:n,j,])
      theta_ss_Aug[,j] <- theta_ss_Aug[,j]+N^{-1}*colSums(m_Aug[((n+1):(n+N)),j,])
    }
    
    theta_ss_Aug <- (t(theta_ss_Aug)+theta_ss_Aug)/2
  }
  
  
  if(hd_Aug==TRUE){
    P_a_Aug <- P_y_line_z_eta_cal(y1,x1,x2,theta_supervise,type_conmodel = "Aug",hd=TRUE)
    m_Aug <- array(0,c(n+N,q,q))
    theta_ss_Aug <- array(0,c(q,q))
    
    for(j in 1:q){
      Hessian <- Hessian_est(theta_supervise[,j],y1,j)
      inv_Hessian <- solve(Hessian)
      
      m_Aug[,j,] <- m_estimation(inv_Hessian,A_y,P_a_Aug,n+N,j,theta_supervise)
      
      theta_ss_Aug[,j] <- theta_supervise[,j]-n^{-1}*colSums(m_Aug[1:n,j,])
      theta_ss_Aug[,j] <- theta_ss_Aug[,j]+N^{-1}*colSums(m_Aug[((n+1):(n+N)),j,])
    }
    
    theta_ss_Aug <- (t(theta_ss_Aug)+theta_ss_Aug)/2
  }
  
  if(hd_PoS==FALSE){
    P_a_PoS <- P_y_line_z_eta_cal(y1,x1,x2,theta_supervise=theta_supervise,hd=FALSE,
                                  type_conmodel="PoS",PoS_model=PoS_model)
    m_PoS <- array(0,c(n+N,q,q))
    theta_ss_PoS <- array(0,c(q,q))
    
    for(j in 1:q){
      Hessian <- Hessian_est(theta_supervise[,j],y1,j)
      inv_Hessian <- solve(Hessian)
      
      m_PoS[,j,] <- m_estimation(inv_Hessian,A_y,P_a_PoS,n+N,j,theta_supervise)
      
      theta_ss_PoS[,j] <- theta_supervise[,j]-n^{-1}*colSums(m_PoS[1:n,j,])
      theta_ss_PoS[,j] <- theta_ss_PoS[,j]+N^{-1}*colSums(m_PoS[((n+1):(n+N)),j,])
    }
    
    theta_ss_PoS <- (t(theta_ss_PoS)+theta_ss_PoS)/2
  }

  if(hd_PoS==TRUE){
    P_a_PoS <- P_y_line_z_eta_cal(y1,x1,x2,theta_supervise=theta_supervise,hd=TRUE,type_conmodel="PoS",PoS_model=PoS_model)
    m_PoS <- array(0,c(n+N,q,q))
    theta_ss_PoS <- array(0,c(q,q))
    
    for(j in 1:q){
      Hessian <- Hessian_est(theta_supervise[,j],y1,j)
      inv_Hessian <- solve(Hessian)
      
      m_PoS[,j,] <- m_estimation(inv_Hessian,A_y,P_a_PoS,n+N,j,theta_supervise)
      
      theta_ss_PoS[,j] <- theta_supervise[,j]-n^{-1}*colSums(m_PoS[1:n,j,])
      theta_ss_PoS[,j] <- theta_ss_PoS[,j]+N^{-1}*colSums(m_PoS[((n+1):(n+N)),j,])
    }
    
    theta_ss_PoS <- (t(theta_ss_PoS)+theta_ss_PoS)/2
  }
    
    L <- list(theta_SL=theta_supervise,theta_SCISS_Aug=theta_ss_Aug,theta_SCISS_PoS=theta_ss_PoS,
              m_Aug=m_Aug,m_PoS=m_PoS,s=s)
    return(L)
}


density_ratio<-function(x1,x2,y)
{
  n <- length(y[,1])
  N <- length(x2[,1])
  q <- length(y[1,])
  p <- length(x1[1,])-1
  
  YY <- array(0,n+N)
  YY[(n+1):(n+N)] <- 1
  XX <- array(0,dim=c(p+1,n+N))
  XX[,1:n] <- x1
  XX[,(n+1):(n+N)] <- x2
  eta <- glm(YY~t(XX[2:(p+1),]),family = binomial)$coef
  w <- array(0,n)
  for(i in 1:n)
  {
    w[i] <- exp(eta%*%x1[i,])
  }
  w <- w/sum(w)
  theta <- array(0,dim=c(q,q))
  for (j in 1:q){
    y_s1 <- y[,-j]
    if(q==2)
    {fit <- glm(y[,j]~y_s1,weight=w,family="binomial")}
    if(q>2)
    {fit <- glm(y[,j]~y_s1,weight=w,family="binomial")}
    theta[j,j] <- coef(fit)[1]
    if(j==1){theta[,j] <- coef(fit)[1:q]}
    else if(j==q){theta[1:(q-1),j] <- coef(fit)[2:q]}
    else{
      theta[1:(j-1),j] <- coef(fit)[2:j]
      theta[(j+1):q,j] <- coef(fit)[(j+1):q]
    }
  }
  theta <- (theta+t(theta))/2
  return(theta)
}


#result table
result_table <- function(to_path,loopnum,q,n,Graph)
{
  theta_sl <- theta_SCISS_Aug <- theta_SCISS_PoS <- theta_dr_result <- array(0,c(q,q,loopnum))
  
  mse_sl <- rse_sl <- bias_sl <- mse_SCISS_Aug <- rse_SCISS_Aug <- bias_SCISS_Aug <- mse_SCISS_PoS <- 
    rse_SCISS_PoS <- bias_SCISS_PoS <- mse_dr <- rse_dr <- bias_dr <- CI_Aug <- CI_PoS <- array(0,c(q,q))

  CI <- array(0,c(q,q,2))
  alpha <- 0.05

  for(loop in 1:loopnum)
  {
      path <- paste(to_path,loop,'.rda',sep='')
      load(path)

      theta_sl[,,loop] <- SCISS_Aug_Re$theta_SL
      theta_SCISS_Aug[,,loop] <- SCISS_Aug_Re$theta_SCISS_Aug
      theta_SCISS_PoS[,,loop] <- SCISS_PoS_Re$theta_SCISS_PoS
      theta_dr_result[,,loop] <- theta_dr
      
      s <- SCISS_Aug_Re$s
      m_Aug <- SCISS_Aug_Re$m
      est_se <- array(0,c(q,q))
      for(j in 1:q)
      {for(k in 1:q)
        {est_se[j,k] <- sum((s[,j,k]+s[,k,j]-m_Aug[1:n,k,j]-m_Aug[1:n,j,k])^2)/(4*n^2)
        }
      }
      CI[,,1] <- SCISS_Aug_Re$theta_SCISS_Aug-sqrt(est_se)*qnorm(1-alpha/2)
      CI[,,2] <- SCISS_Aug_Re$theta_SCISS_Aug+sqrt(est_se)*qnorm(1-alpha/2)
      for(j in 1:q)
      {for(k in 1:q)
        {
          if(Graph[j,k]>=CI[j,k,1]&Graph[j,k]<=CI[j,k,2]){CI_Aug[j,k] <- CI_Aug[j,k]+1}
        }
      }

      m_PoS <- SCISS_PoS_Re$m
      est_se <- array(0,c(q,q))
      for(j in 1:q)
      {for(k in 1:q)
        {
          est_se[j,k]<-sum((s[,j,k]+s[,k,j]-m_PoS[1:n,k,j]-m_PoS[1:n,j,k])^2)/(4*n^2)
        }
      }
      CI[,,1]<-SCISS_PoS_Re$theta_SCISS_PoS-sqrt(est_se)*qnorm(1-alpha/2)
      CI[,,2]<-SCISS_PoS_Re$theta_SCISS_PoS+sqrt(est_se)*qnorm(1-alpha/2)
      for(j in 1:q)
      {for(k in 1:q)
        {
          if(Graph[j,k]>=CI[j,k,1]&Graph[j,k]<=CI[j,k,2]){CI_PoS[j,k] <- CI_PoS[j,k]+1}
        }
      }
  }
  
  mse_sl <- apply((theta_sl-array(Graph,c(q,q,loopnum)))^2,c(1,2),mean)
  rse_sl <- sqrt(mse_sl)
  bias_sl <- apply(theta_sl-array(Graph,c(q,q,loopnum)),c(1,2),mean)
  
  mse_SCISS_Aug <- apply((theta_SCISS_Aug-array(Graph,c(q,q,loopnum)))^2,c(1,2),mean)
  rse_SCISS_Aug <- sqrt(mse_SCISS_Aug)
  bias_SCISS_Aug <- apply(theta_SCISS_Aug-array(Graph,c(q,q,loopnum)),c(1,2),mean)
  
  mse_SCISS_PoS <- apply((theta_SCISS_PoS-array(Graph,c(q,q,loopnum)))^2,c(1,2),mean)
  rse_SCISS_PoS <- sqrt(mse_SCISS_PoS)
  bias_SCISS_PoS <- apply(theta_SCISS_PoS-array(Graph,c(q,q,loopnum)),c(1,2),mean)

  mse_dr <- apply((theta_dr_result-array(Graph,c(q,q,loopnum)))^2,c(1,2),mean)
  rse_dr <- sqrt(mse_dr)
  bias_dr <- apply(theta_dr_result-array(Graph,c(q,q,loopnum)),c(1,2),mean)
  
  RE_sl_SCISS_Aug<-mse_sl/mse_SCISS_Aug
  RE_sl_SCISS_PoS<-mse_sl/mse_SCISS_PoS
  RE_sl_dr<-mse_sl/mse_dr
  
  #print in latex table format
  A <- array(0,c(q*(q-1)/2,13))
  num <- 1
  for(j in 1:(q-1))
  {
    for(k in (j+1):q)
    {
      A[num,1] <- bias_sl[j,k]
      A[num,2] <- rse_sl[j,k]
      
      A[num,3] <- bias_dr[j,k]
      A[num,4] <- rse_dr[j,k]
      A[num,5] <- RE_sl_dr[j,k]
      
      A[num,6] <- bias_SCISS_Aug[j,k]
      A[num,7] <- rse_SCISS_Aug[j,k]
      A[num,8] <- RE_sl_SCISS_Aug[j,k]
      A[num,9] <- CI_Aug[j,k]/loopnum
      
      A[num,10] <- bias_SCISS_PoS[j,k]
      A[num,11] <- rse_SCISS_PoS[j,k]
      A[num,12] <- RE_sl_SCISS_PoS[j,k]
      A[num,13] <- CI_PoS[j,k]/loopnum
      
      num <- num+1
    }
  }

  return(A)
}


#Intrinsic method----

Intrinsic_data_generation <- function(q=3,n,N,seed)
{
  set.seed(seed)
  if(q==3){
    supervise_num <- sort(sample.int(n+N,n))
    Graph <- matrix(c(-2,0,0,0,-1.5,0,1,1,-1),q,q)
    Graph <- pmax(Graph,t(Graph))
    mu <- c(-2, -1.5, -1)
    Graph_without <- Graph
    diag(Graph_without) <- 0
    Data <- IsingSampler(n = n+N, graph = Graph_without, thresholds = diag(Graph), nIter = 200)
    
    S_y <- Data
    S_x <- array(1,dim=c(n+N,q+1))
    for(i in 1:(n+N))
    {
      for(j in 1:q)
      {
        if(S_y[i,j]==1){S_x[i,j+1]=sample(0:1,1,prob = c(0.6, 0.4))}
        else{S_x[i,j+1]=0}
      }
    }
    
    S_y_supervise <- S_y[1:n,]
    S_x_supervise <- S_x[1:n,]
    S_x_unsupervise <- S_x[(n+1):(n+N),]
  }else if(q==4){
    supervise_num <- sort(sample.int(n+N,n))
    Graph <- matrix(c(-2,0,0,1,0,-1.5,0,1,0,0,-2,1,1,1,1,-1),q,q)
    Graph <- pmax(Graph,t(Graph))
    mu <- c(-2, -1.5, -2,-1)
    Graph_without <- Graph
    diag(Graph_without) <- 0
    Data <- IsingSampler(n = n+N, graph = Graph_without, thresholds = diag(Graph), nIter = 200)
    
    S_y <- Data
    S_x <- array(1,dim=c(n+N,q+1))
    for(i in 1:(n+N))
    {
      for(j in 1:q)
      {
        if(S_y[i,j]==1){S_x[i,j+1]=sample(0:1,1,prob = c(0.6, 0.4))}
        else{S_x[i,j+1]=0}
      }
    }}else if(q==5){
      supervise_num <- sort(sample.int(n+N,n))
      Graph <- matrix(c(-2,0,0,0,1,0,-1.5,0,0,1,0,0,-2,0,1,0,0,0,-1.5,1,1,1,1,1,-1),q,q)
      Graph <- pmax(Graph,t(Graph))
      mu <- c(-2, -1.5, -2,-1.5,-1)
      Graph_without <- Graph
      diag(Graph_without) <- 0
      Data <- IsingSampler(n = n+N, graph = Graph_without, thresholds = diag(Graph), nIter = 200)
      
      S_y <- Data
      S_x <- array(1,dim=c(n+N,q+1))
      for(i in 1:(n+N))
      {
        for(j in 1:q)
        {
          if(S_y[i,j]==1){S_x[i,j+1]=sample(0:1,1,prob = c(0.6, 0.4))}
          else{S_x[i,j+1]=0}
        }
      }
    }
  
  S_y_supervise <- S_y[1:n,]
  S_x_supervise <- S_x[1:n,]
  S_x_unsupervise <- S_x[(n+1):(n+N),]
  
  List <- list(S_y_supervise=S_y_supervise,S_x_supervise=S_x_supervise,S_x_unsupervise=S_x_unsupervise)
  return(List)
}


empirical_variance_intrinsic <- function(eta,y1,x1,theta_supervise,j,k)#calculate the empirical variance
{
  n <- nrow(y1)
  q <- ncol(y1)
  p <- ncol(x1)-1
  
  inv_Hessian_theta_j <- solve(Hessian_est(theta_supervise[,j],y1,j))
  inv_Hessian_theta_k <- solve(Hessian_est(theta_supervise[,k],y1,k))
  
  y_setminus_j <- y_setminus_k <- y1
  y_setminus_j[,j] <- 1
  y_setminus_k[,k] <- 1
  
  s_jk <- rep(0,n)
  s_kj <- rep(0,n)
  
  for(i in 1:n){
    y_j <- y1[i,j]
    y_k <- y1[i,k]
    s_jk[i] <- (inv_Hessian_theta_j%*%y_setminus_j[i,]*as.numeric(y_j-g(theta_supervise[,j]%*%y_setminus_j[i,])))[k,1]
    s_kj[i] <- (inv_Hessian_theta_k%*%y_setminus_k[i,]*as.numeric(y_k-g(theta_supervise[,j]%*%y_setminus_k[i,])))[j,1]
  }
  
  P_all <- array(0,n)
  P_a <- matrix(0,n,2^q)
  a <- 1:2^q
  for(i in 1:n){
    for(j1 in 1:q){
      P_a[i,a] <- P_a[i,a]+as.numeric(eta[j,((j-1)*(p+1)+1):(j*(p+1))]%*%
                                        x1[i,])*A_y[a,j1]
      if(j1!=q){
        for(j2 in (j1+1):q)
        {P_a[i,a] <- P_a[i,a]+as.numeric(eta[j,((j1-1)*(p+1)+1):(j1*(p+1))]%*%
                                           x1[i,])*A_y[a,j1]*A_y[a,j2]}}
    }
    P_all[i]<-sum(exp(P_a[i,]))}
  
  m_jk <- rep(0,n)
  m_kj <- rep(0,n)
  for(i in 1:n){
    m_jk_every <- 0
    m_kj_every <-0
    for(a in 1:2^q)
    {
      y_setminus_a_j <- y_setminus_a_k <- A_y[a,]
      y_setminus_a_j[j] <- 1
      y_setminus_a_k[k] <- 1
      y_a_j <- A_y[a,j]
      y_a_k <- A_y[a,k]
      m_jk_every <- (inv_Hessian_theta_j%*%y_setminus_a_j*as.numeric(y_a_j-g(t(theta_supervise[,j])
                                                                             %*%y_setminus_a_j))*as.numeric(exp(P_a[i,a])/P_all[i]))[k]
      m_kj_every <- (inv_Hessian_theta_k%*%y_setminus_a_k*as.numeric(y_a_k-g(t(theta_supervise[,k])
                                                                             %*%y_setminus_a_k))*as.numeric(exp(P_a[i,a])/P_all[i]))[j]
      m_jk[i] <- m_jk[i]+m_jk_every
      m_kj[i] <- m_kj[i]+m_kj_every
    }
  }
  
  empirical_variance <- mean((s_jk+s_kj-m_jk-m_kj)^2)-(mean(s_jk+s_kj-m_jk-m_kj))^2  
  return(empirical_variance)  
}


m_estimation_intri <- function(inv_Hessian,x1,A_y,P_y_line_z_eta,M,theta_supervise,j,k,l)
{
  q <- length(A_y[1,])
  n <- length(x1[,1])
  m <- array(0,c(n,q))
  dm <- array(0,c(n,q))
  ddm <- array(0,c(n,q))
  
  for(i in 1:n){
    S_condition_j_0 <- rep(0,q)
    dS_condition_j_0 <- rep(0,q)
    ddS_condition_j_0 <- rep(0,q)
    
    for(a in 1:2^q){
      y_s <- as.matrix(A_y[a,])
      y_s[j] <- 1
      y_j <- A_y[a,j]
      S_specified_j_0 <- inv_Hessian%*%y_s*as.numeric(y_j-g(t(theta_supervise[,j])%*%y_s))
      S_condition_j_0 <- S_condition_j_0+S_specified_j_0*as.numeric(P_y_line_z_eta[i,a])
      if(j==k){
        dS_condition_j_0 <- dS_condition_j_0+S_specified_j_0*as.numeric(P_y_line_z_eta[i,a])*A_y[a,j]*x1[i,l+1]
        ddS_condition_j_0 <- dS_condition_j_0+S_specified_j_0*as.numeric(P_y_line_z_eta[i,a])*A_y[a,j]*x1[i,l+1]^2
      }else{
        dS_condition_j_0 <- dS_condition_j_0+S_specified_j_0*as.numeric(P_y_line_z_eta[i,a])*A_y[a,j]*A_y[a,k]*x1[i,l+1]
        ddS_condition_j_0 <- dS_condition_j_0+S_specified_j_0*as.numeric(P_y_line_z_eta[i,a])*A_y[a,j]*A_y[a,k]*x1[i,l+1]^2
      }
    }
    m[i,] <- S_condition_j_0
    dm[i,] <- dS_condition_j_0
    ddm[i,] <- ddS_condition_j_0
  }
  List <- list(m=m,dm=dm,ddm=ddm)
  return(List)
}


Intrinsic_iteration <- function(y1,x1,x2,j,k,theta_supervise,upper_times=5,threshold=6)
{
  n <- length(y1[,1])
  q <- length(y1[1,])
  p <- length(x1[1,])-1
  
  gamma_start <- array(0,c(q,q*(p+1)))
  psi_j <- array(0,c(n,q*(p+1)))
  psi_setminus_j <- rep(0,(p+1)*q)
  
  for(j1 in 1:q){
    for(i in 1:n){
      for(t in 1:q){
        if(t==j1){
          psi_setminus_j[((t-1)*(p+1)+1):(t*(p+1))] <- x1[i,]
        }
        else{psi_setminus_j[((t-1)*(p+1)+1):(t*(p+1))] <- x1[i,]*y1[i,t]}
        psi_j[i,] <- psi_setminus_j
      }}
    fit <- glmnet(x <- psi_j, 
                  y<- y1[,j1],family="binomial",lambda=n^{-3/4},alpha=0,intercept=FALSE)
    gamma_start[j1,]<-coef(fit)[2:(q*(p+1)+1)] 
  }
  
  inv_Hessian_theta_j <- solve(Hessian_est(theta_supervise[,j],y1,j))
  inv_Hessian_theta_k <- solve(Hessian_est(theta_supervise[,k],y1,k))
  
  y_setminus_j <- y_setminus_k <- y1
  y_setminus_j[,j] <- 1
  y_setminus_k[,k] <- 1
  
  s_jk <- S_estimation(inv_Hessian_theta_j,y1,theta_supervise,j)[,k]
  s_kj <- S_estimation(inv_Hessian_theta_k,y1,theta_supervise,k)[,j]
  
  gamma_t1 <- array(0,c(q,q*(p+1)))
  for(iter_time in 1:upper_times){
    if(iter_time==1){gamma_t <- gamma_start}
    for(j1 in 1:q){
      for(j2 in 1:q){
        P_all <- array(0,n)
        P_a_log <- matrix(0,n,2^q)
        P_a <- matrix(0,n,2^q)
        a <- 1:2^q
        for(i in 1:n){
          for(k1 in 1:q){
            P_a_log[i,a] <- P_a_log[i,a]+as.numeric(gamma_t[k1,((k1-1)*(p+1)+1):(k1*(p+1))]%*%
                                                      x1[i,])*A_y[a,k1]
            if(k1!=q){
              for(k2 in (k1+1):q)
              {P_a_log[i,a] <- P_a_log[i,a]+as.numeric(gamma_t[k1,((k2-1)*(p+1)+1):(k2*(p+1))]%*%
                                                         x1[i,])*A_y[a,k1]*A_y[a,k2]}}
          }
          P_all[i] <- sum(exp(P_a_log[i,]))
          P_a[i,a] <- exp(P_a_log[i,a])/P_all[i]
        }
        
        for(j3 in 1:p){
          all_m_jk_estimation <- m_estimation_intri(inv_Hessian_theta_j,x1,A_y,P_a,n,theta_supervise,j1,j2,j3)
          all_m_kj_estimation <- m_estimation_intri(inv_Hessian_theta_k,x1,A_y,P_a,n,theta_supervise,j2,j1,j3)
          m_jk <- all_m_jk_estimation$m[,k]
          m_kj <- all_m_kj_estimation$m[,j]
          dm_jk <- all_m_jk_estimation$dm[,k]
          dm_kj <- all_m_kj_estimation$dm[,j]
          ddm_jk <- all_m_jk_estimation$ddm[,k]
          ddm_kj <- all_m_jk_estimation$ddm[,j]
          
          dsig <- 2*mean((s_jk+s_kj-m_jk-m_kj)*(-dm_jk-dm_kj))-
            2*mean(s_jk+s_kj-m_jk-m_kj)*mean(-dm_jk-dm_kj)
          ddsig <- 2*mean((dm_jk+dm_kj)^2+
                            (s_jk+s_kj-m_jk-m_kj)*(-ddm_jk-ddm_kj))-
            2*(mean(dm_jk+dm_kj))^2-2*mean(s_jk+s_kj-m_jk-m_kj)*mean(-ddm_jk-ddm_kj)
          gamma_t1[j1,p*(j2-1)+j3] <- gamma_t[j1,p*(j2-1)+j3]-dsig/ddsig
        }
      }
    }
    if(max(abs(gamma_t1))>threshold){
      break
    }else if(iter_time>1)
    {
      if(empirical_variance_intrinsic(gamma_t1,y1,x1,theta_supervise,j,k)>
         empirical_variance_intrinsic(gamma_t,y1,x1,theta_supervise,j,k)){
        break
      }else{
        gamma_t <- gamma_t1
      }
    }
  }
  
  empirical_variance_begin <- empirical_variance_intrinsic(gamma_start,y1,x1,theta_supervise,j,k)
  empirical_variance_final <- empirical_variance_intrinsic(gamma_t,y1,x1,theta_supervise,j,k)
  
  List <- list(gamma_start=gamma_start,gamma_t=gamma_t,empirical_variance_begin=empirical_variance_begin,
               empirical_variance_final=empirical_variance_final)
  return(List)
}

intrinsic_construct_with_gamma <- function(q,p,x1,x2,y1,gamma,theta_supervise)
{
  n <- nrow(x1)
  N <- nrow(x2)
  
  P_all <- array(0,n+N)
  P_a_log <- matrix(0,n+N,2^q)
  P_a <- matrix(0,n+N,2^q)
  a <- 1:2^q
  for(i in 1:n){
    for(j in 1:q){
      P_a_log[i,a] <- P_a_log[i,a]+as.numeric(gamma[j,((j-1)*(p+1)+1):(j*(p+1))]%*%
                                                x1[i,])*A_y[a,j]
      if(j!=q){
        for(j1 in (j+1):q)
        {P_a_log[i,a] <- P_a_log[i,a]+as.numeric(gamma[j,((j1-1)*(p+1)+1):(j1*(p+1))]%*%
                                                   x1[i,])*A_y[a,j]*A_y[a,j1]
        }}
    }
    P_all[i] <- sum(exp(P_a_log[i,a]))
    P_a[i,a] <- exp(P_a_log[i,a])/P_all[i]
  }
  
  for(i in 1:N){
    for(j in 1:q){
      P_a_log[i+n,a] <- P_a_log[i+n,a]+as.numeric(gamma[j,((j-1)*(p+1)+1):(j*(p+1))]%*%
                                                    x2[i,])*A_y[a,j]
      if(j!=q){
        for(j1 in (j+1):q)
        {P_a_log[i+n,a] <- P_a_log[i+n,a]+as.numeric(gamma[j,((j1-1)*(p+1)+1):(j1*(p+1))]%*%
                                                       x2[i,])*A_y[a,j]*A_y[a,j1]
        }}
    }
    P_all[i+n] <- sum(exp(P_a_log[i+n,a]))
    P_a[i+n,a]<- exp(P_a_log[i+n,a])/P_all[i+n]
  }
  
  m <- array(0,c(n+N,q,q))
  theta_SCISS <- array(0,c(q,q))
  
  for(j in 1:q){
    Hessian <- Hessian_est(theta_supervise[,j],y1,j)
    inv_Hessian <- solve(Hessian)
    
    m[,j,] <- m_estimation(inv_Hessian,A_y,P_a,n+N,j,theta_supervise)
    
    theta_SCISS[,j] <- theta_supervise[,j]-n^{-1}*colSums(m[1:n,j,])
    theta_SCISS[,j] <- theta_SCISS[,j]+N^{-1}*colSums(m[((n+1):(n+N)),j,])
  }
  
  theta_SCISS <- (t(theta_SCISS)+theta_SCISS)/2
  
  List <- list(theta_SCISS=theta_SCISS,m=m) 
  return(List)
}

#construct result
Intrinsic_SCISS <- function(y1,x1,x2)
{
  q <- ncol(y1) 
  n <- nrow(x1)
  N <- nrow(x2)
  p <- ncol(x1)-1
  
  theta_supervise <- SL(y1)$theta_SL 
  Empirical_Variance_begin <- Empirical_Variance_intri <- array(0,c(q,q))
  theta_SCISS_intrinsic <- theta_SCISS_begin <- array(0,c(q,q))
  
  for(j in 1:q){
    for(k in 1:q){
      Intrinsic_jk_result <- Intrinsic_iteration(y1,x1,x2,j,k,theta_supervise,upper_times=5) 
      Empirical_Variance_begin[j,k] <- Intrinsic_jk_result$empirical_variance_begin
      Empirical_Variance_intri[j,k] <- Intrinsic_jk_result$empirical_variance_final
      gamma_begin <- Intrinsic_jk_result$gamma_start
      gamma_final <- Intrinsic_jk_result$gamma_t
      
      theta_SCISS_begin_result <- intrinsic_construct_with_gamma(q,p,x1,x2,y1,gamma_begin,theta_supervise) 
      theta_SCISS_intri_result <- intrinsic_construct_with_gamma(q,p,x1,x2,y1,gamma_final,theta_supervise) 
      
      theta_SCISS_begin[j,k] <- theta_SCISS_begin_result$theta_SCISS[j,k]
      theta_SCISS_intrinsic[j,k] <- theta_SCISS_intri_result$theta_SCISS[j,k]
    }
  }
  
  for(j in 1:q){
    for(k in 1:q){
      if(Empirical_Variance_intri[j,k] < Empirical_Variance_intri[k,j]){
        theta_SCISS_intrinsic[k,j] <- theta_SCISS_intrinsic[j,k]
        Empirical_Variance_intri[k,j] <- Empirical_Variance_intri[j,k]}
      if(Empirical_Variance_begin[j,k] < Empirical_Variance_begin[k,j]){
        theta_SCISS_begin[k,j] <- theta_SCISS_begin[j,k]
        Empirical_Variance_begin[k,j] <- Empirical_Variance_begin[j,k]
      }
    }
  }
  
  List <- list(theta_SCISS_begin=theta_SCISS_begin,theta_SCISS_intri=theta_SCISS_intrinsic,
               Empirical_Variance_begin=Empirical_Variance_begin,
               Empirical_Variance_intri=Empirical_Variance_intri) 
  return(List)
}


Intrinsic_result_analysis <- function(to_path,q,Graph,loopnum)
{
  theta_SCISS_begin <- theta_SCISS_intri <- 
    Begin_SCISS_EAV <- Intrinsic_SCISS_EAV <- array(0,c(loopnum,q,q))
  
  for(loop in 1:loopnum)
  {
    path <- paste(to_path,loop,'.rda',sep='')
    load(path)
    
    theta_SCISS_begin[loop,,] <- Intrinsic_SCISS_Re$theta_SCISS_begin
    theta_SCISS_intri[loop,,] <- Intrinsic_SCISS_Re$theta_SCISS_intri
    Begin_SCISS_EAV[loop,,] <- Intrinsic_SCISS_Re$Empirical_Variance_begin
    Intrinsic_SCISS_EAV[loop,,] <- Intrinsic_SCISS_Re$Empirical_Variance_intri
  }
  
  EAV_begin <- apply(Begin_SCISS_EAV,c(2,3),mean)
  EAV_intri <- apply(Intrinsic_SCISS_EAV,c(2,3),mean)
  mse_begin <- apply((theta_SCISS_begin-array(Graph,c(loopnum,q,q)))^2,c(2,3),mean)
  mse_intri <- apply((theta_SCISS_intri-array(Graph,c(loopnum,q,q)))^2,c(2,3),mean)
  
  List <- list(EAV_begin=EAV_begin,EAV_intri=EAV_intri)
  return(List)
}

#Optimal SL method----
#get the conditional probabillity of y_j|y_\j,w;theta
cdprob <- function(y,theta_est)
{
  n <- length(y[,1])
  q <- length(y[1,])
  prob <- array(0,c(n,q))
  for(j in 1:q){
    y_s <- y
    y_s[,j] <- 1
    for(i in 1:n){
      prob[i,j] <- 1/(1+exp(-theta_est[j,]%*%y_s[i,]))
    }
  }
  return(prob)
}

#get the initialization of estimate
initialize_est <- function(y)
{
  n <- length(y[,1])
  q <- length(y[1,])
  est_ini <- array(0,c(q,q))
  for(j in 1:q){
    for (j in 1:q){
      y_j1 <- y
      y_j1[,j] <- 1
      fit <- glm.fit(x = model.matrix(~ 0+y_j1),
                     y = y[,j],
                     family = binomial(link = "logit"))
      est_ini[,j] <- fit$coefficients
    }
    
    # est_ini[j,j] <- log(sum(y[,j])/(n-sum(y[,j])))
  }
  return(est_ini)
}

loss <- function(y,theta_est)
{
  n <- length(y[,1])
  q <- length(y[1,])
  l1 <- array(0,c(n,q))
  for(j in 1:q){
    y_s <- y
    y_s[,j] <- 1
    for(i in 1:n){
      l1[i,j] <- log(1+exp(theta_est[j,]%*%y_s[i,]))-theta_est[j,]%*%y_s[i,]*y[i,j]#+ or -?
    }
  }
  los <- sum(l1)/n
  return(los)
}

#use coordinate shooting to get the iterative estimate

Joint_SL <- function(y)
{
  n <- length(y[,1])
  q <- length(y[1,])
  est_new <- initialize_est(y)
  
  delta <- 10^-3
  delta1 <- 10^-3
  max_ite <- 10
  
  prob_temp <- cdprob(y,est_new)
  w_temp <- prob_temp*(1-prob_temp)
  
  los <- loss(y,est_new)
  loss_diff <- 1
  ite <- 0
  dist <- 10
  
  while((dist>=delta)&(loss_diff>=delta1)&(ite<=max_ite)){
    dist <- 0
    est_old <- est_new
    
    for(j in 1:q){
      
      for(k in 1:q){
        
        if(k>j){
          t1 <- sum(y[,k]*(y[,j]-prob_temp[,j])+y[,j]*(y[,k]-prob_temp[,k]))/n
          t2 <- sum(y[,k]^2*w_temp[,j]+y[,j]^2*w_temp[,k])/n
          
          est_new[j,k] <- est_old[j,k]-t1/t2
          est_new[k,j] <- est_new[j,k]
          
          prob_temp <- cdprob(y,est_new)
          w_temp <- prob_temp*(1-prob_temp)
          
        }else if(k==j){
          
          t1 <- sum(y[,j]-prob_temp[,j])/n
          t2 <- sum(prob_temp[,j]*(1-prob_temp[,j]))/n
          
          est_new[j,j] <- est_old[j,j]-t1/t2
          
          prob_temp <- cdprob(y,est_new)
          w_temp <- prob_temp*(1-prob_temp)
        }
      }
    }
    
    dist <- sum(abs(est_new-est_old))/length(est_new)
    
    los_old <- los
    los <- loss(y,est_new)
    loss_diff <- los_old-los
    
    ite <- ite+1
  }
  
  s<-array(0,c(q,q,n))
  for(j in 1:q){
    inv_Hessian_est <- inv(Hessian_est(est_new[,j],y,j))
    for(i in 1:n){
      y_s<-y[i,]
      y_s[j]<-1
      y_j<-y[i,j]
      s[j,,i]<-inv_Hessian_est%*%y_s*as.numeric(y_j-g(est_new[,j]%*%y_s))
    }
  }
  
  List <- list(Optimal_SL=est_new, s=s)
  return(List)
}

