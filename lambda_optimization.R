library(tidyverse)


data = read.table("A1_co2.txt",sep=" ",head=TRUE)
N_obs_te = 20
N_obs_tr = nrow(data) - N_obs_te #slice-to index for train-set
train <- data[0:N_obs_tr, ] #collect train dataframe
test <- data[(nrow(data)-19):(nrow(data)), ] #collect test dataframe. stupid non-zero-indexed slicing
p = 12 #period of oscillation. I dont understand if this is a fit-parameter or we are supposed to "guess" it beforehand
omega = 2*pi/p

#time-scale discrete steps j 
t_train = seq(1,length(train$time))
t_test = seq(length(train$time)+1,length(data$time))
t_full = seq(1,length(data$time))

construct_L = function(p){
  #Function handle for constructing L analogous to Madsen eq. (3.86). Need update to include time-term, bad notation in book.
  L = diag(4)*1
  L[1,2] = 0 
  L[2,2] = 1
  L[2,1] = 1
  L[3,3] = cos(2*pi/p)
  L[4,4] = cos(2*pi/p)
  L[4,3] = -sin(2*pi/p)
  L[3,4] = sin(2*pi/p)
  L #implicit returnation
  L #implicit returnation
}

get_fj <- function(j,p){
  Fj = c(1,j,sin(j*2*pi/p),cos(j*2*pi/p))
  Fj
}

#Filter data using local linear trend on the whole test-set 
#Find F718, h718, theta718

get_FNj <- function(lambda,N,p){ #has been checked with maple with N=2
  #Function to generate FN(j)
  seqsum = matrix(0,nrow=4,ncol=4)
  for (j in 0:(N-1)){
    fj = get_fj(j,p)
    seqsum = seqsum + lambda^j*fj%*%t(fj)
  }
  seqsum
}

get_HNj <- function(lambda,Y,p){ #has been checked with maple for Y = co2[717:718]
  seqsum = matrix(0,nrow=4,ncol=1)
  N = length(Y)
  for (j in 0:(N-1)){
    seqsum = seqsum + lambda^j*get_fj(j,p)*Y[N-j]
  }
  seqsum 
}

generate_design_matrix <- function(N){
  #generate design-matrix equal to xN in book
  #inputs: N: number of observations
  design_matrix = matrix(0,nrow=N,ncol=4)
  for (j in N:1){
    j_vals = N-j
    design_matrix[j, ] = get_fj(j_vals,12)
  }
  design_matrix
}



construct_sigma_inv <- function(N,lambda){
  vec = c()
  for (i in 1:N){
    vec = c(vec,lambda^(N-i))
  }
  Sigma = diag(vec)
  Sigma
}

#Use params for one-step prediction errors. Need transciency steps, using 10 
optimize_lambda <- function(lambda) {
  transciency_period = 10
  FNj = matrix(0,4,4)
  HNj = matrix(0,4,1)
  L = construct_L(12)
  L_inv = solve(L)
  #for (j in 0:9){ #Ends at FN10 and HN10 
  #Y = train$co2[j+1] #using up to next point in sequence
  #Old lines: seem to work to well, but should not work well? 
  #xN = generate_design_matrix(j)
  #FNj = FNj + get_FNj(lambda,j,12)
  #HNj = HNj + get_HNj(lambda,Y,12)
  #FNj = FNj + lambda^j*get_fj(j,12)%*%t(get_fj(j,12))
  #HNj = HNj + lambda*L_inv%*%HNj + get_fj(0,12)*Y
  #}
  HNj = get_HNj(lambda,train$co2[1:10],10)
  FNj = get_FNj(lambda,10,12)
  params_est_10 = solve(FNj)%*%HNj #PARAMS10 after transciency 
  #Ends at params estimated at 10 
  start_idx = transciency_period + 1
  #now filter rest of data 
  params_est_N = params_est_10
  one_step_pred_err = c() 
  L = construct_L(12)
  L_inv = solve(L)
  collective_error = c()
  sqrt_factor = c() #for generating PI's
  for (j in 1:10){
    collective_error[j] = 0
    sqrt_factor[j] = 0
  }
  
  
  thetas = c()
  for (j in 11:(length(train$co2)-1)){#j is current N+1. makes prediction for point j and constructs new FNj and HNj
    #print(c("Generating prediction for point",j))
    ynew = get_fj(1,12)%*%params_est_N #predict Y(j) | Y(N-1)
    filtered_preds_onestep = c(filtered_preds_onestep,ynew) #save prediction
    one_step_pred_err = c(one_step_pred_err, train$co2[j]-ynew) #save error
    err_t = train$co2[1:j]-generate_design_matrix(j)%*%params_est_N
    thetas = c(thetas,params_est_N[1])
    T = 0
    for (t in 0:(j-2)){
      T = T+lambda^t
    }
    last_FN = get_FNj(lambda,j-1,12)
    sqrt_factor = c(sqrt_factor,sqrt(1+t(get_fj(1,12))%*%solve(last_FN)%*%get_fj(1,12)))
    collective_error[j] = t(err_t)%*%construct_sigma_inv(j,lambda)%*%err_t/(T-length(params_est_N))  
    #update params (recursive update when YN+1 available)
    
    #FNj = FNj + lambda^(j-1)*get_fj(j,12)%*%t(get_fj(j,12))
    #HNj = lambda*L_inv%*%HNj + get_fj(0,12)*train$co2[j+1]
    FNj = get_FNj(lambda,j,12) #FN(j) for point  #remember j = N+1
    HNj = get_HNj(lambda,train$co2[1:j],12) #HN(j)
    #FNj = FNj + lambda^(j-1) + get_fj((j-1),12)%*%t(get_fj((j-1),12)) #FN(j)
    #HNj = lambda*L_inv%*%HNj + get_fj(0,12)*train$co2[j] #HN(j)
    params_est_N = solve(FNj)%*%HNj #params(j) ie params N + 1 
  }
  #get last point without parameter-update 
  #print(c("Generating prediction for point: ",length(train$co2)))
  ynew = get_fj(1,12)%*%params_est_N
  thetas = c(thetas,params_est_N[1])
  filtered_preds_onestep = c(filtered_preds_onestep,ynew)
  err_t = c(err_t,train$co2[718]-ynew)
  
  #return sum of squared errs
  sum_out = sum(err_t)
  sum_out
  
}

filtered_preds_onestep = c()
lambda = seq(0.9146912192, 0.9872585449,0.001)
error_holder = c()
lambdas = c()
for (i in 1:length(lambda)){
  print(c(i,"/",length(lambda),":",lambda[i]))
  SE = optimize_lambda(lambda[i])
  error_holder = c(error_holder,SE)
  lambdas = c(lambdas,lambda[i])
}

error_df = data.frame(lambdas,error_holder)
ggplot() + geom_point(mapping=aes(x=lambdas,y=error_holder)) + xlab("Lambda-value") + ylab("Sum of errors")
