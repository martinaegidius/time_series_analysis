#data-prep
data = read.table("A1_co2.txt",sep=" ",head=TRUE)
N_obs_te = 20
N_obs_tr = nrow(data) - N_obs_te #slice-to index for train-set
train <- data[0:N_obs_tr, ] #collect train dataframe
test <- data[(nrow(data)-19):(nrow(data)), ] #collect test dataframe. stupid non-zero-indexed slicing

#exercise 1: plot
plot(train$co2 ~train$time,col="steelblue",pch=20,cex=0.6,xlab="Decimal time of measurement",ylab="CO2-level [ppm] ") #makes plot-object and adds points
points(test$co2 ~test$time,col="red",pch=20,cex=0.6) #adds points onto existing plot-object
legend("topleft",c("Train-set","Test-set"),col=c("steelblue","tomato"),lty=2,bty='n',lwd=2)

library(tidyverse)
#redo exc 1 plot. Need to figure out how to add two data-sources
theme_update(legend.position='top') 

ggplot() + 
  geom_point(data=train,mapping=aes(color="train-set",x=time, y=co2),size=0.6,shape=19) + 
  geom_point(data=test,mapping=aes(color="test-set",x=time, y=co2),size=0.6,shape=19) +
  scale_color_manual(name="Mauna Loa observations", 
                     breaks = c("train-set","test-set"),
                     values=c('train-set'="blue","test-set"="red")) + 
  xlab("Decimal time") + 
  ylab("CO2 [PPM]") + 
  theme(legend.text=element_text(size=12)) + 
  theme(axis.text = element_text(face="bold",size=11))

  


three_year_sample = data[11:47-1, ]
lines(three_year_sample$co2~three_year_sample$time)
seasons = c("winter","spring","summer","fall") #make for seasons
seasons = rep(seasons,each=3)
seasons = rep(seasons,times=c(3))
seasons = c(seasons,"winter")
#seasons = cbind(seasons,c(1))
print(seasons)
print(nrow(three_year_sample))
print(length(seasons))
seasonal_coded = cbind(three_year_sample,seasons) #cat season class
#still has todo
ggplot(data=three_year_sample) + 
  geom_line(mapping=aes(x=time, y=co2),linestyle="dashed") + 
  geom_point(mapping=aes(x=time, y=co2,color=seasons),size=3,shape=19) + 
  xlab("Decimal time") + 
  ylab("CO2 [PPM]") + 
  theme(legend.text=element_text(size=12)) + 
  theme(axis.text = element_text(face="bold",size=11))
#we notice a underlying seasonal trend within the range of a year and an underlying systematic increase, possibly linear. Plotting a period of two three years, we notice that CO2-levels tend to rise after the middle of fall-periods and decrease after the end of spring. Thus, there seems to be a systematic oscillation with period 1 year. Thus, we expect to need some sort of seasonal model with an underlying added monthly term of increase

#Question 1.2. We construct a new X-matrix based on the model-formulation
p = 1 #period of oscillation. I dont understand if this is a fit-parameter or we are supposed to "guess" it beforehand
omega = 2*pi/p
#create Vandermond
X = cbind(1,train$time,sin(omega*train$time),cos(omega*train$time)) #should probably be transpose
X_FULL = cbind(1,data$time,sin(omega*data$time),cos(omega*data$time))
X_TEST = cbind(1,test$time,sin(omega*test$time),cos(omega*test$time))
#solve OLS problem
X = cbind(1,train$time,sin(omega*train$time),cos(omega*train$time)) #construct design matrix
(OLS <- solve(t(X)%*%X)%*%t(X)%*%train$co2) #solve for OLS-estimate of theta
pred_vals = OLS[1]+OLS[2]*train$time+OLS[3]*sin(omega*train$time)+OLS[4]*cos(train$time)
#plot fitted values together with the data
ggplot(data=train) + 
  geom_point(mapping=aes(color="Observed",x=time, y=co2),size=0.6,shape=19) + 
  geom_point(mapping=aes(color="OLS-fit",x=time, y=pred_vals),size=0.6,shape=19) + 
  scale_color_manual(name="OLS-model on train-data", 
                     breaks = c("Observed","OLS-fit"),
                     values=c('Observed'="purple","OLS-fit"="cyan")) + 
  xlab("Decimal time") + 
  ylab("CO2 [PPM]") + 
  theme(legend.text=element_text(size=12)) + 
  theme(axis.text = element_text(face="bold",size=11))

#make OLS with lm to also get estimates of uncertainty of parameters
OLS_mlf = lm(train$co2 ~ train$time + sin(omega*train$time)+cos(omega*train$time))
summary(OLS_mlf)
confint(OLS_mlf,level=0.9) #gives CI for parameters. Hereby found 

y_pred = fitted(OLS_mlf)#extract fitted values as vector
ggplot(data=train) + geom_point(mapping=aes(x=time, y=co2),size=0.6,shape=19,color="blue") + geom_point(mapping=aes(x=time, y=y_pred),size=0.6,shape=19,color="red")# want: connect scatter-points and change axis labels
#still same - luckily 

#question 1.4: relaxation algorithm
construct_matrix <- function(N,rho){
  #Constructs square covariance matrix of form given in assignment. Every element is filled with 1,rho,rho^2,rho^3...rho^N
  holder_mat <- matrix(1,nrow=N,ncol=N)
  for (i in 1:N) {
    for (j in 1:N) {
      holder_mat[i,j] <- rho^(abs(j-i))
    }
  }
  holder_mat #implicit return of function
}

relaxation_algorithm <- function(X,y,num_iter,used_lag){
  N = length(y)
  Sigmay = diag(N)
  for (i in 1:num_iter){
    param_est = solve(t(X)%*%solve(Sigmay)%*%X)%*%(t(X)%*%solve(Sigmay)%*%y)
    #error is a constant term. Therefore we need no further advanced estimate of var/covar-structure <- yes we do, we need the sammenhæng between rho and variance/covariance
    err_t <- y - X%*%param_est #get residuals 
    #get lag 1 autocorrelation between residuals
    rho_est = acf(err_t)$acf[used_lag]
    print(c("rho_est: ",rho_est))
    Sigmay = construct_matrix(N,rho_est)
  }
  #uncertainty calculation - get variance of estimator
  #sigma_hat = t((y-X%*%param_est))%*%solve(Sigmay)%*%(y-X%*%param_est)/(N-length(param_est))
  #V_est = sigma_hat*solve(t(X)%*%solve(Sigmay)%*%X)
  #print(V_est)
  output = list(param_est,Sigmay,err_t) #implicit return
  output
}

auto_correlation_function <- function(vals,k){
  #calculates the k-lag autocorrelation function (ACF) manually 
  #probably a bit overkill, use ACF instead 
  holder_mat = matrix(1,nrow=length(vals),ncol=length(vals))
  dims = length(vals)
  t = k+1
  for (i in 1:dims){
    for (j in 1:dims){
      #holder_mat[i,j] = 
        #https://blog.quantinsti.com/autocorrelation-autocovariance/
    }
  }
}

out_list = relaxation_algorithm(X,train$co2,5,2) #3 seems to give best increase of performance in SSE compared to OLS
params_WLS = out_list[1][[1]]
Sigmay = out_list[2][[1]]
Sigmay_inv = solve(Sigmay)
t_err = out_list[3][[1]]

print(params_WLS)
#TODO need find uncertainty 
#using unbiased estimator of variance
sigma_hat_sq = (t(t_err)%*%Sigmay_inv%*%t_err)/(length(train$co2)-length(params_WLS))
v_hat = sigma_hat_sq[1]*solve(t(X)%*%Sigmay_inv%*%X) #probably diag is variance 

#construct model predictions for trainset and test-set 
WLS_preds = X_FULL%*%params_WLS
OLS_preds = X_FULL%*%OLS
ggplot(data=data) +
  geom_point(mapping=aes(x=time,y=co2),size=0.6,shape=19,color="blue") + 
  geom_point(mapping=aes(x=X_FULL[,2],y=WLS_preds),size=0.6,shape=19,color="green") + 
  geom_point(mapping=aes(x=X_FULL[,2], y=OLS_preds),size=0.6,shape=19,color="red") 

res_OLS = sum((data$co2-OLS_preds)^2)
res_WLS = sum((data$co2-WLS_preds)^2)
print(c("SSE OLS: ",res_OLS,"SSE WLS: ",res_WLS,"diff: ",res_OLS-res_WLS))
#generate train and test error: 
train_err_OLS = sum((train$co2-X%*%OLS)^2)
train_err_WLS = sum((train$co2-X%*%params_WLS)^2)
test_err_OLS = sum((test$co2-X_TEST%*%OLS)^2)
test_err_WLS = sum((test$co2-X_TEST%*%params_WLS)^2)
print(c("||| Train errors: ||| OLS: ",train_err_OLS," WLS: ",train_err_WLS))
print(c("||| Test errors: ||| OLS: ",test_err_OLS," WLS: ",test_err_WLS))

#thought: what happens if you do all of this instead on decimal time on numerated time-steps? 
# t_train = seq(1,length(train$time))
# t_test = seq(length(train$time)+1,length(data$time))
# t_full = seq(1,length(data$time))
# 
# #Trend model. Now p=12, as it needs to be relative to stepsize j 
# 
# construct_L = function(p,j){
#   #Function handle for constructing L analogous to Madsen eq. (3.86). Need update to include time-term, bad notation in book.
#   L = diag(3)*cos(2*pi/p)
#   L[1,1] = L[1,1]/cos(2*pi/p)
#   L[3,2] = -sin(2*pi/p)
#   L[2,3] = sin(2*pi/p)
#   L #implicit returnation
# }
# p = 12 
# t_train = seq(1,length(train$time))
# t_test = seq(1,length(test$time))
# t_full = seq(1,length(data$time))
# F0 = c(1,0,1)
# L = construct_L(p)
# 
# generate_f_vector <- function(L,f0,n){
#   holder_mat = matrix(0,nrow=3,ncol=n)
#   holder_mat[,1] = f0
#   for (i in 2:n){
#     res = L%*%f0
#     holder_mat[,i] = res
#     f0 = res
#     }
#   holder_mat #works but calculate F1 and F2 theoretical and check. Remember to use it reversed
# }
# 
# generate_sigma_matrix <- function(lambda,n){
#   #not usefull, as matrix becomes computationally singular
#   holder_mat = diag(n)
#   for (i in 1:n){
#     holder_mat[i,i] = 1/lambda^(n-i)
#   }
#   holder_mat #works but calculate F1 and F2 theoretical and check. Remember to use it reversed
# }
# 
# generate_lambda_vector <- function(lambdaV,n){
#   lambda = rep(1,n)
#   #should be filled using some sort of apply function but i dont have time
#   #thus stupid slow fill procedure
#   for (i in 1:n){
#     lambda[i] = 1/(lambdaV^(n-i))
#   }
#   lambda
# }
# 
# #check fn[,2], fn[,3] with maple: gives correct 
# F_mat = t(generate_f_vector(L,F0,length(t_train))) #f(j)
# #need reverse of F_mat: 
# X_lt = F_mat[nrow(F_mat):1,] #generate local trend design matrix
# lambda = 0.9
# Sigmay_lt = generate_sigma_matrix(lambda,length(t_train)) #is computationally singular. Use computing formula (3.100) instead
# #FN = t(X_lt)%*%solve(Sigmay_lt)%*%X_Lt #numerically unstable. Computational procedure instead
# #calculate sum by hand (could also be done more intelligently)
# seqsum = 0
# lambda_vec = generate_lambda_vector(lambda,length(t_train))
# 
# get_FN <- function(lambdaV,X_lt,N){
#   #takes X_lt as F which is in same format as book
#   #Function to produce FN = (t(XN)%*%XN)3.100)
#   seqsum = matrix(0,nrow=3,ncol=3)
#   for (j in 0:(N-1)){ #error subscript out of bounds at some point
#     prod = (lambdaV^j)*(X_lt[N-j, ]%*%t(X_lt[N-j, ]))
#     seqsum = seqsum + prod
#     if(j<2){
#       print(prod)
#       print(j)
#     }
#   }
#   seqsum
# }
# get_HN <- function(lambdaV,X_lt,Y,N){
#   seqsum = 0
#   for (j in 0:(N-1)){ #error subscript out of bounds at some point
#     prod = (lambdaV^j)*(X_lt[N-j, ])*Y[N-j]
#     seqsum = seqsum + prod
#     if(j<2){
#       print(prod)
#       print(j)
#     }
#   }
#   seqsum
# }
# 
# get_local_trend_params <- function(FN,HN){
#   theta_est = solve(FN)%*%HN
#   theta_est
# }
# 
# FN = get_FN(lambda,X_lt,length(t_train))
# HN = get_HN(lambda,X_lt,train$co2,length(t_train))
# params_lt_N = get_local_trend_params(FN,HN)
# 
# run_local_trend_model <- function(x,y,lambda,F0,transciency_period){
#   #Run iterative local trend model for inputs. 
#     #Args: 
#       #x: time-value 
#       #y: measured y-vals. In i'th step of algorithm, y[i] is first used at the end for updating params, but is included in error calculation. Should probably not be, as we simulate "not knowing" the point
#       #lambda: forgetting-factor 
#       #F0: initial value
#       #transciency-period: the first number of p observations used for generating initial parameter estimates
#       #L: the L-matrix which governs updating the design-matrix to N points
#     #Returns: 
#       #out_list: a list containing
#         #[1][[1]]: vector of next-in-series predictions
#         #[2][[1]]: unbiased estimator of variance of next-in-series-predictions
#   L_inv = solve(L)
#   F_mat = t(generate_f_vector(L,F0,(transciency_period))) #f(j)
#   #need reverse of F_mat: 
#   X_lt = F_mat[nrow(F_mat):1,] #generate local trend design matrix
#   #use the first 10 points as "transient": i guess means use for first estimate
#   FN = get_FN(lambda,X_lt,transciency_period)
#   HN = get_HN(lambda,X_lt,y,transciency_period)
#   print("GENERATED MATRICES SUCCESSFULLY. HN; FN")
#   print(FN)
#   print(HN)
#   params_tra = get_local_trend_params(FN,HN)
#   tval <- transciency_period+1
#   print(c("TVAL IS: ", tval))
#   for (i in tval:length(x)){#should be dynamic but does not work
#     #estimate 
#     print(c("Estimating point: ",i))
#     #updating f(j) -> f(l)
#     F_mat = rbind(F_mat,t(L%*%F_mat[dim(F_mat)[1], ])) #F_mat is f(0)...f(N), so first row is 1,0,1
#     X_lt = F_mat[nrow(F_mat):1,] #generate local trend design matrix
#     
#     if (i<12) { #get all values for data 1:11 and establish vectors to hold data
#       print("Went into first loop")
#       step = 1
#       print(dim(X_lt))
#       preds = F_mat[1:(i-1), ] %*% params_tra
#       pred_p1 = F_mat[i, ] %*% params_tra
#       print("got to here")
#       preds_v = c(pred_p1)
#       err_t = y[1:(i-1)]-preds
#       sigma_sq_est = sum(err_t^2)/(i-length(params_tra))
#       sigma_sq_est_v = c(sigma_sq_est)
#       pred_err_pi = y[i] - pred_p1
#       pred_err_pi_v = c(pred_err_pi)
#     }
#     else{ #update for 1:11+i (ie 1:N+1)
#       print("WENT INTO SECOND LOOP")
#       #step = step + 1 
#       print("Here 1")
#       #print(dim(X_lt))
#       pred_p1 = F_mat[i, ]%*% params_step_i
#       print("Here 2")
#       preds_N = F_mat[1:(i-1), ] %*% params_step_i
#       preds_v = c(preds_v,pred_p1)
#       err_t = y[1:(i-1)]-preds_N
#       sigma_sq_est = sum(err_t^2)/(i-length(params_step_i))
#       sigma_sq_est_v = c(sigma_sq_est_v,sigma_sq_est)
#       pred_err_pi = y[i]-pred_p1
#       pred_err_pi_v = c(pred_err_pi_v,pred_err_pi)
#     }
#     if(i<=length(x)){ #FOR UPDATING MATRICES WHEN Y(N+1) IS AVAILABLE (ie after preds)
#       #New values will be used for params in next iteration. 
#       print("WENT INTO MATRIX-UPDATE-LOOP")
#       print("UPDATING FN")
#       FN = FN+F_mat[(i-1), ]%*%t(F_mat[(i-1), ]) #f(-N) is always going to be equal to non reversed matrix F(i)
#       #Update HN
#       print("UPDATING HN")
#       HN = (L_inv%*%HN)+F_mat[1, ]*y[i]
#       print("UPDATING PARAMS")
#       params_step_i = get_local_trend_params(FN,HN)
#     }  
#     
#     print(c("Predicted value: ",pred_p1," True value: ",y[i]))
#   }
#   out_list = list(preds_v,sigma_sq_est_v)
#   out_list #implicit return
# }  
# #Q1.2 "filter the data with the chosen model. I understand it as: build model, estimate params from full sequence ("filter") and make predictions for sequence
# params_lt = solve(FN)%*%HN 
# 
# #change of understanding: i think we should run algo iteratively and predict all points as this is basically the model 
# out_list = run_local_trend_model(t_train,train$co2,lambda,c(1,0,1),10,L)
# model_preds = out_list[1][[1]]
# model_vars = out_list[2][[1]]
# transient_values = train$co2[1:10]
# model_preds = c(transient_values,model_preds)
# 
# filter_historic_preds = X_lt%*%params_lt
# ggplot(data=train) +
#   geom_point(mapping=aes(x=time,y=co2),size=0.6,shape=19,color="blue") + 
#   geom_point(mapping=aes(x=time,y=model_preds),size=0.6,shape=19,color="green") +  #MUST BE WRONG
#   geom_point(mapping=aes(x=time,y=WLS_preds[1:718]),size=0.4,shape=19,color="red")
#   
# 
# 

