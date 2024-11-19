This code is for implementation of project *Robust and Efficient Semi-supervised Learning for Ising Model*. It is based on R language.
# Before start  
Ensure the following packages have been installed:  
*MASS*, *stats*, *purrr*, *glmnet*, *IsingSampler*.  

# Introduction to main functions  
The multiple results return in a list. Use $ operator to output each result.  
### SCISS(y1, x1, x2, hd_Aug, hd_PoS, PoS_model)
**Description**    
The function is used for running both SCISS-Aug and SCISS-PoS methods for regular and high-dimensional setting.  
**Argument**  
y1: observed outcomes  
x1: auxiliary features for observed outcomes  
x2: auxiliary features for unobserved outcomes  
hd_Aug: logical; If TRUE, run SCISS-Aug code for high-dimensional setting, if FALSE, run for low-dimensional setting. Default is FALSE.  
hd_PoS: logical; If FALSE, run SCISS-PoS code for high-dimensional setting, if FALSE, run for low-dimensional setting. Default is FALSE.   
PoS_model: one of \{"gaussian","poisson","binomial"\}.  The model chosen for SCISS-PoS.  
**Return Value**  
theta_SL: supervised-learnin(SL) estimator   
theta_SCISS_Aug: SCISS-Aug estimator  
theta_SCISS_PoS: SCISS-PoS estimator  
s: score function for supervised-learning  
m_Aug=m_Aug: imputed score function in SCISS-Aug  
m_PoS=m_PoS: imputed score function in SCISS-PoS  

### Generate_y_x(Graph, p, n, N, c, seed, type, hd=FALSE)  
**Description**    
The function is used for generating regular and high-dimensional data for simulation.    
**Argument**  
Graph: true value of Ising parameter matrix
p: dimensionality of auxiliary feature
n: number of data with observed outcome    
N: number of data without observed outcome  
c: coefficient matrix used to decide the effect of interaction terms of y to x  
seed: seed used to keep reproducible in generation  
type: one of \{"gaussian", "poisson"\}, used to generate the data corresponding to the type  
hd: logical; If TRUE, generate high-dimensional data, if FALSE, generate low-dimensional data  
**Return Value**  
S_y_supervise: observed outcome by generation  
S_x_supervise: auxiliary features for observed outcome by generation   
S_x_unsupervise: auxiliary features for unobserved outcome by generation  

### density_ratio(x1, x2, y)
**Description**    
The function is used for running density-ratio(DR). We use it as a benchmark to SCISS.  
**Argument**  
x1: auxiliary features for observed outcomes    
x2: auxiliary features for unobserved outcomes    
y: observed outcomes     
**Return Value**  
theta: density ratio estimator  

### result_table(to_path,loopnum,q,n,Graph) 
**Description**    
The function is used for analyzing the result of running replications of SCISS function.     
**Argument**  
to_path: path that stores the result  
loopnum:  the number of replications  
q: dimensionality of outcome  
n: number of data with observed outcome    
Graph: true value of Ising parameter matrix   
**Return Value**  
result:  analysis result, each row is the result for one Ising parameter, columns are arranged in this order  
[1] bias of SL estimator [2] root square error(RSE) of SL estimator [3] bias of DR estimator [4] RSE of DR estimator [5] relative efficiency(RE) of DR estimator
[6] bias of SCISS-Aug estimator [7] RSE of SCISS-Aug estimator  [8] RE of SCISS-Aug estimator   
[9] coverage probability(CP) of SCISS-Aug estimator  [10] bias of SCISS-PoS estimator   [11] RSE of SCISS-PoS estimator  
[12] RE of SCISS-PoS estimator [13] CP of SCISS-PoS estimator  

### Intrinsic_data_generation(q,n,N,seed)
**Description**    
The function is used for generating data for intrinsic method. To prove that intrinsic method works well when the model is misspecified, we design the misspecified data in purpose.     
**Argument**  
q: dimensionality of outcome    
n: number of data with observed outcome      
N: number of data without observed outcome   
seed: seed used to keep reproducible in generation    
**Return Value**  
S_y_supervise: observed outcome by generation  
S_x_supervise: auxiliary features for observed outcome by generation   
S_x_unsupervise: auxiliary features for unobserved outcome by generation  

### Intrinsic_SCISS(y1,x1,x2)
**Description**    
The function is used for running intrinsic method.   
**Argument**  
y1: observed outcomes  
x1: auxiliary features for observed outcomes  
x2: auxiliary features for unobserved outcomes    
**Return Value**  
theta_SCISS_begin: SCISS estimator before intrinsic updation   
theta_SCISS_intri: SCISS-INTR estimator  
Empirical_Variance_begin: empirical variance before intrinsic updation  
Empirical_Variance_intri: empirical variance after intrinsic updation    


### Joint_SL(y)
**Description**    
The function is used for running "optimal" SL method.  
**Argument**  
y: observed outcome   
**Return Value**  
Optimal_SL: optimal SL estimator   
s: score function in optimal SL method  

# Example     
See Example.R for an example. This example includes 500 replications of running SCISS method for regular configuration and analysis of the result. 


