#library----
library(MASS)
library(purrr)
library(stats)
library(matlib)
library(glmnet)
library(IsingSampler)

#function----
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

Hessian_est<-function(n,theta_j,S_y,j)
{
  Hessian_est_single <- array(0,dim <- c(q,q))
  for(i in 1:n) 
  {y_setminus_j <- t(S_y)[,i]
  y_setminus_j[j]<- 1
  Hessian_est_single <- 
    Hessian_est_single +(1/n)*y_setminus_j%*%t(y_setminus_j)*as.numeric(dg(theta_j%*%y_setminus_j))}
  return(Hessian_est_single)
}

t_s<-function(y1,q)
{
  theta_supervise=array(0,c(q,q))
  for (j in 1:q){
    y_s<-t(y1)
    y_s1<-y_s[-j,]
    if(q==2)
    {fit<-glm(t(y1)[j,]~y_s1,family="binomial")}
    if(q>2)
    {fit<-glm(t(y1)[j,]~t(y_s1),family="binomial")}
    theta_supervise[j,j]<-coef(fit)[1]
    if(j==1){theta_supervise[,j]=coef(fit)[1:q]}
    else if(j==q){theta_supervise[1:(q-1),j]=coef(fit)[2:q]}
    else{
      theta_supervise[1:(j-1),j]=coef(fit)[2:j]
      theta_supervise[(j+1):q,j]=coef(fit)[(j+1):q]
    }
  }
  theta_supervise=(theta_supervise+t(theta_supervise))/2
  s<-array(0,c(q,q,n))
  for(i in 1:n)
  {
    for(j in 1:q)
    {
      y_s=t(y1)[,i]
      y_s[j]=1
      y_j=t(y1)[j,i]
      s[j,,i]=inv(Hessian_est(n,theta_supervise[,j],y1,j))%*%y_s*as.numeric(y_j-g(theta_supervise[,j]%*%y_s))
    }
  }
  L=list(s=s,theta_SL=theta_supervise)
  return(L)
}

SCISS_Aug<-function(n,N,y1,x1,x2,q,p)
{
  A_y<-array(0,c(q,2^q))
  col_sum<-1
  for(i in 1:q){
    combination<-combn(q, i)
    for(t in 1:ncol(combination))
    {
      for(k in 1:i){
        y_order<-combination[k,t]
        A_y[y_order,col_sum+t]<-1
      }
    }
    col_sum=col_sum+ncol(combination)
  }
  T_S<-t_s(y1,q)
  theta_supervise<-T_S$theta_SL
  s<-T_S$s
  
  gamma_supervise<-array(0,c(q,q*(p+1)))
  psi_j<-array(0,c(q*(p+1),n))
  psi_sj<-rep(0,(p+1)*q)
  for(j in 1:q){
    for(i in 1:n){
      for(t in 1:q){
        if(t==j){
          psi_sj[((t-1)*(p+1)+1):(t*(p+1))]<-x1[,i]
        }
        else{psi_sj[((t-1)*(p+1)+1):(t*(p+1))]<-x1[,i]*(t(y1)[t,i])}
        psi_j[,i]=psi_sj
      }}
    fit<-cv.glmnet(t(psi_j), 
                t(y1)[j,],alpha=0,family="binomial",lambda=c(1:n^{1/2})/n,intercept=FALSE)
    gamma_supervise[j,]<-coef(fit)[2:(q*(p+1)+1)] 
  }
  
  theta_ss<-array(0,dim=c(q,q))
  m<-array(0,dim=c(q,q,n+N))
  for(j in 1:q){
    Hessian <- Hessian_est(n,theta_supervise[,j],y1,j)
    inv_Hessian<-solve(Hessian)
    for(i in 1:n){
      P_all=0
      P_a=array(0,2^q)
      for(a in 1:2^q)
      {
        for(j_1 in 1:q)
        {
          P_a[a]<-P_a[a]+gamma_supervise[j_1,((j_1-1)*(p+1)+1):(j_1*(p+1))]%*%
            x1[,i]*A_y[j_1,a]
        }
        for(j_1 in 1:q)
        {
          if(j_1==q){}
          else{
            for(j_2 in (j_1+1):q)
            {P_a[a]<-P_a[a]+gamma_supervise[j_1,((j_2-1)*(p+1)+1):(j_2*(p+1))]%*%
              x1[,i]*A_y[j_1,a]*A_y[j_2,a]
            }}
        }
        P_all<-P_all+exp(P_a[a])
      }
      S_condition_j_0=array(0,dim=c(q,1))
      for(a in 1:2^q)
      {
        y_s<-A_y[,a]
        y_s[j]<-1
        y_j<-A_y[j,a]
        S_specified_j_0<-inv_Hessian%*%y_s*as.numeric(y_j-g(t(theta_supervise[,j])%*%y_s))
        S_condition_j_0 <-S_condition_j_0+S_specified_j_0*as.numeric(exp(P_a[a])/P_all)
      }
      m[j,,i]<-S_condition_j_0
    }
    
    for(i in 1:N)
    {
      P_all=0
      P_a=array(0,dim=2^q)
      for(a in 1:2^q)
      {
        for(j_1 in 1:q)
        {
          P_a[a]<-P_a[a]+gamma_supervise[j_1,((j_1-1)*(p+1)+1):(j_1*(p+1))]%*%x2[,i]*A_y[j_1,a]
        }
        for(j_1 in 1:q)
        {
          if(j_1==q)
          {}
          else{
            for(j_2 in (j_1+1):q)
            {
              P_a[a]<-P_a[a]+gamma_supervise[j_1,((j_2-1)*(p+1)+1):(j_2*(p+1))]%*%
                x2[,i]*A_y[j_1,a]*A_y[j_2,a]}}
        }
        P_all<-P_all+exp(P_a[a])
      }
      S_condition_j_0<-array(0,dim=c(q,1))
      for(a in 1:2^q)
      {
        y_s<-A_y[,a]
        y_s[j]<-1
        y_j<-A_y[j,a]
        S_specified_j_0<-inv_Hessian%*%y_s*
          as.numeric(y_j-g(t(theta_supervise[,j])%*%y_s))
        S_condition_j_0 <-S_condition_j_0+S_specified_j_0*as.numeric(exp(P_a[a])/P_all)
      }
      m[j,,i+n]<-S_condition_j_0
    }
    theta_ss[,j]<-theta_supervise[,j]-n^{-1}*rowSums(m[j,,1:n])
    theta_ss[,j]<-theta_ss[,j]+N^{-1}*rowSums(m[j,,((n+1):(n+N))])
  }
  theta_ss<-(t(theta_ss)+theta_ss)/2
  
  L=list(theta_SL=theta_supervise,theta_SCISS_Aug=theta_ss,m=m,s=s)
  return(L)
}

SCISS_PoS<-function(n,N,y1,x1,x2,q,p,TYPE)
{
  A_y<-array(0,c(q,2^q))
  col_sum<-1
  for(i in 1:q){
    combination<-combn(q, i)
    for(t in 1:ncol(combination))
    {
      for(k in 1:i){
        y_order<-combination[k,t]
        A_y[y_order,col_sum+t]<-1
      }
    }
    col_sum=col_sum+ncol(combination)
  }
  T1<-t_s(y1,q)
  theta_supervise<-T1$theta_SL
  s<-T1$s
  #beta----
  beta<-array(0,c(q,2))
  for(j in 1:q)
  {
    YY<-t(y1)
    Y_j<-YY[j,]
    if(TYPE=="gaussian")
    {fit<-glm(x1[j+1,]~Y_j,family="gaussian")}
    if(TYPE=="binomial")
    {fit<-glm(x1[j+1,]~Y_j,family='binomial')}
    if(TYPE=="poisson")
    {fit<-glm(x1[j+1,]~Y_j,family='poisson')}
    beta[j,]<-coef(fit)
  }
  
  #P_ylinew----
  P_ylinew<-array(0,2^q)
  P_ylinew_all<-0
  for(a in 1:2^q)
  {
    for(j1 in 1:q)
    {
      P_ylinew[a]<-P_ylinew[a]+theta_supervise[j1,j1]*A_y[j1,a]
      if(j1<q)
      {
        for(j2 in (j1+1):q)
        {
          P_ylinew[a]<-P_ylinew[a]+theta_supervise[j1,j2]*A_y[j1,a]*A_y[j2,a]
        }
      }
    }
    P_ylinew_all<-P_ylinew_all+exp(P_ylinew[a])
  }
  
  
  theta_ss<-array(0,dim=c(q,q))
  m<-array(0,c(q,q,n+N))
  for(j in 1:q){
    Hessian <- Hessian_est(n,theta_supervise[,j],y1,j)
    inv_Hessian<-solve(Hessian)
    #1-n----
    for(i in 1:n){
      #P_xylinew----
      P_xylinew<-array(0,2^q)
      for(a in 1:2^q)
      {
        P_xylinew[a]<-exp(P_ylinew[a])/P_ylinew_all
        for(j1 in 1:q)
        {
          Y_j<-c(1,1)
          Y_j[2]<-A_y[j1,a]
          if(TYPE=="gaussian")
          {P_xylinew[a]=P_xylinew[a]*dnorm(x1[j1+1,i],Y_j%*%beta[j1,],0.1)}
          if(TYPE=="binomial")
          {P_xylinew[a]=P_xylinew[a]*dbinom(x1[j1+1,i],1,prob=exp(Y_j%*%beta[j1,])/(1+exp(Y_j%*%beta[j1,])) )}
          if(TYPE=="poisson")
          {
            P_xylinew[a]=P_xylinew[a]*dpois(x1[j1+1,i],exp(Y_j%*%beta[j1,]),log=FALSE)
        }
        }
      }
      
      #p_ylinexw----
      P_ylinexw<-array(0,2^q)
      if(sum(P_xylinew==array(0,2^q))==2^q)
      {
        P_xylinexw=array(1/2^q,2^q)
      }
      if(sum(P_xylinew==array(0,2^q))!=2^q){
        for(a in 1:2^q)
        {
          P_ylinexw[a]<-P_xylinew[a]/sum(P_xylinew)
        }
      }
      
      S_condition_j_0=array(0,dim=c(q,1))
      for(a in 1:2^q)
      {
        y_s<-A_y[,a]
        y_s[j]<-1
        y_j<-A_y[j,a]
        S_specified_j_0<-inv_Hessian%*%y_s*as.numeric(y_j-g(theta_supervise[,j]%*%y_s))
        S_condition_j_0 <-S_condition_j_0+S_specified_j_0*P_ylinexw[a]
      }
      m[j,,i]<-S_condition_j_0
    }
    
    #(n+1):M----
    for(i in 1:N){
      #P_xylinew----
      P_xylinew<-array(0,2^q)
      for(a in 1:2^q)
      {
        P_xylinew[a]<-exp(P_ylinew[a])/P_ylinew_all
        for(j1 in 1:q)
        {
          Y_j<-c(1,1)
          Y_j[2]<-A_y[j1,a]
          if(TYPE=="gaussian")
          {P_xylinew[a]=P_xylinew[a]*dnorm(x2[j1+1,i],Y_j%*%beta[j,],0.1)}
          if(TYPE=="binomial")
          {P_xylinew[a]=P_xylinew[a]*dbinom(x2[j1+1,i],1,prob=exp(Y_j%*%beta[j,])/(1+exp(Y_j%*%beta[j,])) )}
          if(TYPE=="poisson"){P_xylinew[a]=P_xylinew[a]*dpois(x2[j1+1,i],exp(Y_j%*%beta[j,]),log=FALSE)}
        }
      }
      
      #p_ylinexw----
      P_ylinexw<-array(0,2^q)
      if(sum(P_xylinew==array(0,2^q))==2^q)
      {
        P_xylinexw=array(1/2^q,2^q)
      }
      if(sum(P_xylinew==array(0,2^q))!=2^q){
        for(a in 1:2^q)
        {
          P_ylinexw[a]<-P_xylinew[a]/sum(P_xylinew)
        }
      }
      
      S_condition_j_0=array(0,dim=c(q,1))
      for(a in 1:2^q)
      {
        y_s<-A_y[,a]
        y_s[j]<-1
        y_j<-A_y[j,a]
        S_specified_j_0<-inv_Hessian%*%y_s*as.numeric(y_j-g(t(theta_supervise[,j])%*%y_s))
        S_condition_j_0 <-S_condition_j_0+S_specified_j_0*P_ylinexw[a]
      }
      m[j,,i+n]<-S_condition_j_0
    }
    
    theta_ss[,j]<-theta_supervise[,j]-n^{-1}*rowSums(m[j,,1:n])
    theta_ss[,j]<-theta_ss[,j]+N^{-1}*rowSums(m[j,,((n+1):(n+N))])
  }
  
  theta_ss<-(t(theta_ss)+theta_ss)/2
  
  # w<-array(0,c(q,q))
  # theta_ss_com<-array(0,c(q,q))
  # m_com<-array(0,c(q,q))
  # for(j in 1:q)
  # {
  #   for(k in 1:q)
  #   {
  #     D_sl<-mean((s[j,k,]+s[k,j,])^2)/4/n
  #     D_ssl<-mean((s[j,k,]+s[k,j,]-m[j,k,1:n]-m[k,j,1:n])^2)/4/n
  #     cov_sl_ssl<-cov(theta_supervise[j,],theta_ss[j,])
  #     ww<-optimize(fw,lower=0,upper=1,tol=0.000000001)
  #     w[j,k]<-ww$minimum
  #     m_com[j,k]<-ww$objective
  #     theta_ss_com[j,k]<-w[j,k]*theta_supervise[j,k]+(1-w[j,k])*theta_ss[j,k]
  #   }
  # }
  
  L<-list(theta_SL=theta_supervise,theta_SCISS_PoS=theta_ss,m=m,s=s)
  return(L)
}


density_ratio<-function(x1,x2,y,p,q,n,N)#q=2
{
  YY=array(0,n+N)
  YY[(n+1):(n+N)]=1
  XX=array(0,dim=c(p+1,n+N))
  XX[,1:n]=x1
  XX[,(n+1):(n+N)]=x2
  eta=glm(YY~t(XX[2:(p+1),]),family = binomial)$coef
  w=array(0,n)
  for(i in 1:n)
  {
    w[i]=exp(eta%*%x1[,i])
  }
  w=w/sum(w)
  theta=array(0,dim=c(q,q))
  for (j in 1:q){
    y_s<-y
    y_s1<-y_s[-j,]
    if(q==2)
    {fit<-glm(y[j,]~y_s1,weight=w,family="binomial")}
    if(q>2)
    {fit<-glm(y[j,]~t(y_s1),weight=w,family="binomial")}
    theta[j,j]<-coef(fit)[1]
    if(j==1){theta[,j]=coef(fit)[1:q]}
    else if(j==q){theta[1:(q-1),j]=coef(fit)[2:q]}
    else{
      theta[1:(j-1),j]=coef(fit)[2:j]
      theta[(j+1):q,j]=coef(fit)[(j+1):q]
    }
  }
  theta=(theta+t(theta))/2
  return(theta)
}

