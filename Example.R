#One example for running SCISS under regular setting

#library the package and load the functions----
source("F:/0SSL_Ising/1revision/Githubcode_redo/SCISS_simulation_functions.R")
source("F:/0SSL_Ising/1revision/Githubcode_redo/library.R")

#setting----
q <- 3
p <- q
n <- 200
N <- 10000

Graph2 <- matrix(c(0.1,0.3,-0.6,0.3,-0.3,0.4,-0.6,0.4,0.2),3,3)

cc2 <- matrix(c(2.5,0.2,0.5,0.2,2.5,0.5,0.5,0.5,2.5),3,3)

loopnum <- 500
A_y <- A_y_generation(3)
save_path <- 'F:/0SSL_Ising/1revision/Githubcode_redo/Github/result/q=3/'
for(loop in 1:loopnum)
{
  Generation <- Generate_y_x(Graph2,p,n,N,cc2,loop,"poisson")
  S_y_supervise <- Generation$S_y_supervise
  S_x_supervise <- Generation$S_x_supervise
  S_x_unsupervise <- Generation$S_x_unsupervise
  
  SCISS_Re <- SCISS(S_y_supervise,S_x_supervise,S_x_unsupervise,PoS_model="poisson")
  theta_dr <- density_ratio(S_x_supervise,S_x_unsupervise,S_y_supervise)
  
  path <- paste(save_path,loop,'.rda',sep='')
  
  save(Graph=Graph2,SCISS_Re,theta_dr,file=path)
}

#result analysis----
load_path <- 'F:/0SSL_Ising/1revision/Githubcode_redo/Github/result/q=3/'
result <- result_table(load_path,loopnum,q,n,Graph2)
print(result)
