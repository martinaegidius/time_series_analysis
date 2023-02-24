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

FN = get_FNj(0.9,718,12)
HN = get_HNj(0.9,train$co2,12)
params_est_N = solve(FN)%*%HN

#Use params for filtering data at end-point

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

xN = generate_design_matrix(length(train$co2))
filter_preds = xN%*%params_est_N
filtered_df = data.frame(train$time,filter_preds,train$co2)
ggplot(data=filtered_df) + 
  geom_point(mapping=aes(x=train$time,y=filter_preds,color="Filtered values"),size=0.6) + 
  geom_point(mapping=aes(x=train$time,y=train$co2,color="Observed values"),size=0.6) + 
  scale_color_manual(name="Filtered data-set", 
                     breaks = c("Filtered values","Observed values"),
                     values=c("Filtered values"="cyan","Observed values"="orange")) + 
  xlab("Decimal time") + 
  ylab("CO2-concentration [ppm]") + 
  theme(legend.text=element_text(size=12)) + 
  theme(axis.text = element_text(face="bold",size=11)) 
  



#Use params for one-step prediction errors. Need transciency steps, using 10 
filtered_preds_onestep = c()
lambda = 0.9
transciency_period = 10
FNj = matrix(0,4,4)
HNj = matrix(0,4,1)
L = construct_L(12)
L_inv = solve(L)
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

construct_sigma_inv <- function(N,lambda){
  vec = c()
  for (i in 1:N){
    vec = c(vec,lambda^(N-i))
  }
  Sigma = diag(vec)
  Sigma
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
  
  FNj = get_FNj(0.9,j-1,12) #FN(j) for point  #remember j = N+1
  HNj = get_HNj(0.9,train$co2[1:j],12) #HN(j)
  params_est_N = solve(FNj)%*%HNj #params(j) ie params N + 1 
  if(params_est_N[2] < 0){
    print(c("Negative value occured in theta :",j+1,params_est_N[2]))
  }
}
#get last point without parameter-update 
print(c("Generating prediction for point: ",length(train$co2)))
ynew = get_fj(1,12)%*%params_est_N
thetas = c(thetas,params_est_N[1])
filtered_preds_onestep = c(filtered_preds_onestep,ynew)
err_t = c(err_t,train$co2[718]-ynew)
#filtered preds_onestep holds predictions for t=11:718

#also need last var
T = 0
for (t in 0:(length(train$co2)-1)){
  T = T+lambda^t
}
sqrt_factor = c(sqrt_factor,sqrt(1+t(get_fj(1,12))%*%solve(get_FNj(0.9,length(train$co2)-1,12))%*%get_fj(1,12)))
collective_error[length(train$co2)] = t(err_t)%*%construct_sigma_inv(length(train$co2),lambda)%*%err_t/(T-length(params_est_N))  
###end last var calculation
train_transcient_time = train$time[11:length(train$time)]
train_transcient_co2 = train$co2[11:length(train$co2)]
stds = sqrt(collective_error)

#make plot one step prediction error
one_step_pred_errtmp = c(rep(NA,10),one_step_pred_err,train$co2[718]-ynew)
one_step_errors = data.frame(train$time,one_step_pred_errtmp,stds)
g17 = ggplot(data=one_step_errors) +
  geom_point(mapping=aes(x=train$time,y=one_step_pred_errtmp,color="Errors"),size=0.9) + 
  scale_color_manual(name="One-step errors", 
                     breaks = c("Errors"),
                     values=c("Errors"="red")) + 
  xlab("Decimal time") + 
  ylab("Error") + 
  theme(legend.text=element_text(size=12)) + 
  theme(axis.text = element_text(face="bold",size=11)) + 
  ylim(-5,3)

g18 = ggplot(data=one_step_errors) +
  geom_point(mapping=aes(x=train$time,y=stds,color="Std"),size=0.9) + 
  scale_color_manual(name="One-step errors", 
                     breaks = c("Std"),
                     values=c("Std"="red")) + 
  xlab("Decimal time") + 
  ylab("Standard deviation") + 
  theme(legend.text=element_text(size=12)) + 
  theme(axis.text = element_text(face="bold",size=11)) + 
  ylim(1,2)

grid.arrange(g17,g18)

#get PI's 
df_v = (10-length(params_est_N)):(length(train$co2)-1-length(params_est_N))
t_vec = qt(0.975,df_v) #0.975 = 95% CI 
stds = stds[11:length(stds)]
sqrt_factor = sqrt_factor[11:length(sqrt_factor)]
lower_pi = filtered_preds_onestep - t_vec*stds*sqrt_factor
higher_pi = filtered_preds_onestep + t_vec*stds*sqrt_factor

#plot one-step predictions and training-data
filtered_stepwise_df = data.frame(train_transcient_time,filtered_preds_onestep,train_transcient_co2,stds,lower_pi,higher_pi)
g1 = ggplot(data=filtered_stepwise_df) + 
  geom_point(mapping=aes(x=train_transcient_time,y=filtered_preds_onestep,color="One step filtered values"),size=0.6) + 
  geom_point(mapping=aes(x=train_transcient_time,y=train_transcient_co2,color="Observed values"),size=0.6) + 
  geom_line(mapping=aes(x=train_transcient_time,y=higher_pi,color="95% PI")) + 
  geom_line(mapping=aes(x=train_transcient_time,y=lower_pi)) + 
  scale_color_manual(name="Filtered data-set", 
                     breaks = c("One step filtered values","Observed values","95% PI"),
                     values=c("One step filtered values"="cyan","Observed values"="orange","95% PI"="black")) + 
  xlab("Decimal time") + 
  ylab("CO2-concentration [ppm]") + 
  theme(legend.text=element_text(size=12)) + 
  theme(axis.text = element_text(face="bold",size=11)) + 
  ylim(310,417)

show(g1)

#plot local estimate of std
g2 = ggplot(data=filtered_stepwise_df) + 
  geom_point(mapping=aes(x=train_transcient_time,y=stds),size=0.6) + 
  ylim(1,2.1)

require(gridExtra)
grid.arrange(g1,g2,nrow=2)

#generate predictions in test-data
test_preds = c()
last_FN = get_FNj(0.9,718,12)
FN_inv = solve(last_FN)
test_sqrt_factor = c()
#Not assuming more data-points arrive
for (j in 1:20){
  thetas = c(thetas,params_est_N[1])
  test_preds = c(test_preds,get_fj(j,12)%*%params_est_N)
  test_sqrt_factor = c(test_sqrt_factor,sqrt(1+t(get_fj(j,12))%*%FN_inv%*%get_fj(j,12)))
}
t_test = qt(0.975,718-4)
test_lower_pi = test_preds - t_test*stds[718-10]*test_sqrt_factor
test_higher_pi = test_preds + t_test*stds[718-10]*test_sqrt_factor

test_vals_tab = c(test$co2[1],test$co2[2],test$co2[6],test$co2[12],test$co2[20])
pred_vals_tab = c(test_preds[1],test_preds[2],test_preds[6],test_preds[12],test_preds[20])

test_table = cbind(test_vals_tab,pred_vals_tab,test_vals_tab-pred_vals_tab)

test_df = data.frame(test$time,test_preds,test_lower_pi,test_higher_pi,test$co2)
filtered_stepwise_df = subset(filtered_stepwise_df,train_transcient_time >= 2010,select=c(train_transcient_time,filtered_preds_onestep,train_transcient_co2,stds,lower_pi,higher_pi))
g3 = ggplot() + 
  geom_point(data=test_df,mapping=aes(x=test$time,y=test_preds,color="Forecast"),size=2) + 
  geom_line(data=test_df,mapping=aes(x=test$time,y=test_lower_pi,color="95% PI")) + 
  geom_line(data=test_df,mapping=aes(x=test$time,y=test_higher_pi,color="95% PI")) + 
  geom_point(data=test_df,mapping=aes(x=test$time,y=test$co2,color="Observed test"),size=2) + 
  geom_point(data=filtered_stepwise_df, mapping=aes(x=train_transcient_time,y=filtered_preds_onestep,color="Train fit"),size=1.5) + 
  geom_point(data=filtered_stepwise_df,mapping=aes(x=train_transcient_time,y=train_transcient_co2,color="Observed train"),size=1.5) + 
  geom_line(data=filtered_stepwise_df,mapping=aes(x=train_transcient_time,y=higher_pi,color="95% PI")) + 
  geom_line(data=filtered_stepwise_df,mapping=aes(x=train_transcient_time,y=lower_pi)) + 
  scale_color_manual(name="Forecast data-set", 
                     breaks = c("Forecast","95% PI","Observed test", "Train fit", "Observed train"),
                     values=c("Forecast"="cyan","Observed test"="orange","95% PI"="black","Observed train"="red","Train fit"="purple")) + 
  xlab("Decimal time") + 
  ylab("CO2-concentration [ppm]") + 
  theme(legend.text=element_text(size=12)) + 
  theme(axis.text = element_text(face="bold",size=11)) + 
  ylim(380,416)

show(g3)
# 
# #assuming future points may be used (they arrive before next month)
# test_preds = c()
# FN = get_FNj(0.9,718,12)
# test_sqrt_factor = c()
# one_step_pred_err = c()
# collective_error = c()
# sqrt_factor = c()
# for (j in 1:20){
#   ynew = get_fj(1,12)%*%params_est_N
#   test_preds = c(test_preds,ynew)
#   one_step_pred_err = c(one_step_pred_err, test$co2[j]-ynew) #save error
#   err_t = data$co2[1:(718+j)]-generate_design_matrix(718+j)%*%params_est_N
#   T = 0
#   for (t in 0:(718+j-2)){
#     T = T+lambda^t
#   }
#   last_FN = get_FNj(lambda,(718+j-1),12)
#   sqrt_factor = c(sqrt_factor,sqrt(1+t(get_fj(1,12))%*%solve(last_FN)%*%get_fj(1,12)))
#   collective_error[j] = t(err_t)%*%construct_sigma_inv(718+j,lambda)%*%err_t/(T-length(params_est_N))  
#   #update params (recursive update when YN+1 available)
#   
#   #FNj = FNj + lambda^(j-1)*get_fj(j,12)%*%t(get_fj(j,12))
#   #HNj = lambda*L_inv%*%HNj + get_fj(0,12)*train$co2[j+1]
#   FNj = get_FNj(0.9,718+j,12) #FN(j) for point  #remember j = N+1
#   HNj = get_HNj(0.9,data$co2[1:(718+j)],12) #HN(j)
#   #FNj = FNj + lambda^(j-1) + get_fj((j-1),12)%*%t(get_fj((j-1),12)) #FN(j)
#   #HNj = lambda*L_inv%*%HNj + get_fj(0,12)*train$co2[j] #HN(j)
#   params_est_N = solve(FNj)%*%HNj #params(j) ie params N + 1 
# }
# 
# 
# df_v = (718-length(params_est_N)):(length(data$co2)-1-length(params_est_N))
# t_vec = qt(0.975,df_v) #0.975 = 95% CI 
# stds = sqrt(collective_error)
# lower_pi = test_preds - t_vec*stds*sqrt_factor
# higher_pi = test_preds + t_vec*stds*sqrt_factor
# 
# test_df_arrival = data.frame(test$time,test_preds,test_lower_pi,test_higher_pi,test$co2)
# g4 = ggplot() + 
#   geom_point(data=test_df_arrival,mapping=aes(x=test$time,y=test_preds,color="Forecast"),size=2) + 
#   geom_line(data=test_df_arrival,mapping=aes(x=test$time,y=test_lower_pi,color="95% PI")) + 
#   geom_line(data=test_df_arrival,mapping=aes(x=test$time,y=test_higher_pi,color="95% PI")) + 
#   geom_point(data=test_df_arrival,mapping=aes(x=test$time,y=test$co2,color="Observed test"),size=2) + 
#   geom_point(data=filtered_stepwise_df, mapping=aes(x=train_transcient_time,y=filtered_preds_onestep,color="Train fit"),size=0.6) + 
#   geom_point(data=filtered_stepwise_df,mapping=aes(x=train_transcient_time,y=train_transcient_co2,color="Observed train"),size=0.6) + 
#   geom_line(data=filtered_stepwise_df,mapping=aes(x=train_transcient_time,y=higher_pi,color="95% PI")) + 
#   geom_line(data=filtered_stepwise_df,mapping=aes(x=train_transcient_time,y=lower_pi)) + 
#   scale_color_manual(name="Forecast data-set", 
#                      breaks = c("Forecast","95% PI","Observed test", "Train fit", "Observed train"),
#                      values=c("Forecast"="cyan","Observed test"="orange","95% PI"="black","Observed train"="red","Train fit"="purple")) + 
#   xlab("Decimal time") + 
#   ylab("CO2-concentration [ppm]") + 
#   theme(legend.text=element_text(size=12)) + 
#   theme(axis.text = element_text(face="bold",size=11)) + 
#   ylim(380,416)
# show(g4)
# 
# 
# #plot three year subset
# filtered_stepwise_df = subset(filtered_stepwise_df,train_transcient_time >= 1980 | train_transcient_time < 1983,select=c(train_transcient_time,filtered_preds_onestep,train_transcient_co2))
# ggplot(data=filtered_stepwise_df) + 
#   geom_point(mapping=aes(x=train_transcient_time,y=filtered_preds_onestep,color="One step filtered values"),size=0.6) + 
#   geom_point(mapping=aes(x=train_transcient_time,y=train_transcient_co2,color="Observed values"),size=0.6) + 
#   scale_color_manual(name="Filtered data-set", 
#                      breaks = c("One step filtered values","Observed values"),
#                      values=c("One step filtered values"="cyan","Observed values"="orange")) + 
#   xlab("Decimal time") + 
#   ylab("CO2-concentration [ppm]") + 
#   theme(legend.text=element_text(size=12)) + 
#   theme(axis.text = element_text(face="bold",size=11)) + 
#   ylim(300,400)
# 
# 
# 

#plot mean values at all time-steps 
thetas = c(rep(NA,10),thetas)
meandf = data.frame(data$time,data$co2,thetas)
ggplot(data=meandf) +
  geom_point(mapping=aes(x=data$time,y=data$co2,color="Observed"),size=1.5) + 
  geom_line(mapping=aes(x=data$time,y=thetas,color="Estimated mean"),size=1.5) + 
  scale_color_manual(name="Means", 
                     breaks = c("Observed", "Estimated mean"),
                     values=c("Observed"="red","Estimated mean"="cyan")) + 
  xlab("Decimal time") + 
  ylab("CO2-concentration [ppm]") + 
  theme(legend.text=element_text(size=12)) + 
  theme(axis.text = element_text(face="bold",size=11)) + 
  ylim(300,416)
