#setting----
q <- 3
p <- q
n <- 200
N <- 10000
M <- n+N

A_y<-array(0,dim=c(q,2^q))
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


G<-matrix(c(0,0.3,-0.6,0.3,0,0.4,-0.6,0.4,0),q,q)
mu=c(0.1,-0.3,0.2)
Graph<-G
diag(Graph)<-mu

cc<-matrix(c(2.5,0.2,0.5,0.2,2.5,0.5,0.5,0.5,2.5),q,q)

#replication----
for(loop in 1:500)
{
  set.seed(loop)
  supervise_num<-sort(sample.int(M,n))
  #generate data
  Data <- IsingSampler(M,G,thresholds=mu,method="CFTP")
  S_y<-t(Data)
  S_y_supervise<-S_y[,supervise_num]
  S_x<-array(1,c(p+1,M))
  
  for(i in 1:M)
  {
    noise<-rnorm(p,0,1)
    for(j in 1:q)
    {
      S_x[j+1,i]<-S_y[j,i]*cc[j,j]
      for(k in 1:q)
      {
        if(k!=j)
        {
          S_x[j+1,i]<-S_x[j+1,i]+cc[j,k]*S_y[j,i]*S_y[k,i]
        }
      }
    }
    S_x[2:(p+1),i]<-S_x[2:(p+1),i]+noise
  }
  
  S_x_supervise<-S_x[,supervise_num]
  S_x_unsupervise<-S_x[,-supervise_num]

  SCISS_Aug_Re<-SCISS_Aug(n,N,t(S_y_supervise),S_x_supervise,S_x_unsupervise,q,p)
  SCISS_PoS_Re<-SCISS_PoS(n,N,t(S_y_supervise),S_x_supervise,S_x_unsupervise,q,p,"gaussian")
  theta_dr<-density_ratio(S_x_supervise,S_x_unsupervise,S_y_supervise,p,q,n,N)
  path<-paste(loop,'.rda',sep='')
  save(G,mu,SCISS_Aug_Re,SCISS_PoS_Re,theta_dr,file =path)
}





