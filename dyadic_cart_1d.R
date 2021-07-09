#n = 2^l, the number of leaves of the binary tree.
#depth of the tree is l + 1. Root is at depth 1. 
#I will encode a node by its depth and position from left to right. 
#At depth i there are 2^{i - 1} nodes. 
#Node (i,j) represents the interval [n/(2^{i - 1})]*(j - 1) + 1 to n/(2^{i - 1})*j.


###Bottom up pass#########
#initialization
###opt,split,sum,sumsq,size are all lists indexed by i the depth of the tree.
###For each i they are a vector.
dyadic_1d = function(l,y,lambda){
  n = 2^l
  opt = list()
  split = list()
  sum = list()
  sumsq = list()
  size = list()
  var = list()
  partition = list()
  ##partition is a list of lists. Each vertex has a list which denotes the optimal
  ##partition of that vertex. A partition is denoted by a list of two tuples 
  ##representing the end points of the intervals.
  for (i in 1:l){
    opt[[i]] = rep(0,2^{i - 1})
    split[[i]] = rep(0,2^{i - 1})
    sum[[i]] = rep(0,2^{i - 1})
    sumsq[[i]] = rep(0,2^{i - 1})
    var[[i]] = rep(0,2^{i - 1})
    depth = l + 2 - i
    size[[i]] = rep(2^{l-i+1},2^{i - 1})
  }
  
  
  opt[[l + 1]] = rep(lambda,2^{l})
  split[[l + 1]] = rep(0,2^{l})
  sum[[l + 1]] = y
  sumsq[[l + 1]] = y^2
  size[[l + 1]] = rep(1,n)
  var[[l + 1]] = rep(0,n)
  
  ##defining the current list. It is a list of lists. 
  currentlist = list()
  for (i in 1:n){
    currentlist[[i]] = list(c(i,i))
    #currentlist[[i]][1] = c(i,i)
  }
  ####Figure out how to remove for loop here. Or merge this with the bottom
  ##up pass.
  
  ##bottom up pass
  for (i in 2:(l + 1)){
    depth = l + 2 - i
    templist = currentlist
    currentlist = NULL
    currentlist = list()
    for (j in 1:(2^{depth - 1})){
      childdep = depth + 1
      lchild = 2*(j - 1) + 1
      rchild = 2*(j - 1) + 2
      sum[[depth]][j] = sum[[depth + 1]][lchild] + sum[[depth + 1]][rchild]
      sumsq[[depth]][j] = sumsq[[depth + 1]][lchild] + sumsq[[depth + 1]][rchild]
      var[[depth]][j] = sumsq[[depth]][j] - (sum[[depth]][j])^2/size[[depth]][j]
      temp = opt[[depth + 1]][lchild] + opt[[depth + 1]][rchild]
      opt[[depth]][j] = min(var[[depth]][j] + lambda,temp)
      currentlist[[j]] = list()
      currentlist[[j]][[1]] = c((n/(2^{depth - 1}))*(j - 1) + 1,(n/(2^{depth - 1}))*(j)) 
      if (var[[depth]][j] + lambda > temp){
        split[[depth]][j] = 1
        #currentlist[[j]][[1]] = templist[[2*(j - 1) + 1]]
        #currentlist[[2]] = templist[[2*(j - 1) + 2]]
        currentlist[[j]] = c(templist[[2*(j - 1) + 1]],templist[[2*(j - 1) + 2]])
      }
      
    }
  }
  finalpart = currentlist[[1]]
  output = rep(0,n)
  for (i in 1:length(finalpart)){
    output[finalpart[[i]][1]:finalpart[[i]][2]] = mean(y[finalpart[[i]][1]:finalpart[[i]][2]])
  }
  return(list(finalpart,output))
}
#########################Bottom up pass done######################
##Main thing: The optimal partition is being computed bottom up. 
##Every step of the for loop creates the optimal partition at each node.
##The optimal partition is encoded via a list. The list is a list of two
##tuples.
##################MSE experiment###########



ind = function(x,a,b){
  if (x > a && x <= b){
    return(1)
  }
  else {
    return(0)
  }
}

f = function(x){
  return(ind(x,0,0.4))
}


f2 = function(x){
  a1 = 2*ind(x,0,0.4)
  a2 = 4*ind(x,0.4,0.6)
  a3 = ind(x,0.6,0.8)
  a4 = 4*ind(x,0.8,1)
  return(a1 + a2 + a3 + a4)
}

f3 = function(x) {
  return(sin(x))
}

f4 = function(x) {
  return(x^2)
}

f5 = function(x){
  a1 = 2*ind(x,0.2,0.4)
  a2 = 5*ind(x,0.4,0.6)
  a3 = ind(x,0.6,0.8)
  a4 = 4*ind(x,0.8,0.9)
  a5 = 3*ind(x, 0.9, 1)
  return(a1 + a2 + a3 + a4 +a5)
}

f6 = function(x) {
  a1 = .3*ind(x, 0, 0.2)
  a2 = .2*ind(x, 0.2, 0.4)
  a3 = .25*ind(x, 0.4, 0.5)
  a4 = .1*ind(x, 0.5, 0.6)
  a5 = .35*ind(x, 0.6, 0.75)
  a6 = .2*ind(x, 0.75, 0.9)
  a7 = .3*ind(x, 0.9, 1)
  
  return(a1+a2+a3+a4+a5+a6+a7)
}

f7 = function(x) {
  return(1/x)
}

mse = function(iter){
  grid = seq(5,20,by = 1)
  lambda = 
    for (l in grid){
      n = 2^l
      theta = sapply(seq(1:n)/n,f2)
      for (j in 1:iter){
        y = theta + rnorm(2^l,0,sigma) #theta is signal, rnorm is noise
        #ans = dyadic_1d(l,y,lambda[l])[[2]]
        ans = dyadic_1d(l,y,lambda)[[2]]
        mse[j] = mean((theta - ans)^2)
      }
      return(mse)
    }
}

#####Algorithm Background Code for Two-Fold Cross Validation----

#this function takes a vector y, and uses the odd observations, and fills the even elements with the average of neighbors
crossval_odd = function(y) { # y is the list of observations
  
  n = length(y)
  y_odd = vector(mode = "numeric", length = n) # list of odd indexed observations from y
  
  for (i in 1:n) {
    if (i %% 2 != 0) { # odd indexes
      y_odd[i] = y[i]
    } else if (i %% 2 == 0 && i != n) { # even indexes
      y_odd[i] = (y[i-1] + y[i+1])/2 # fill in missing observations with average of neighboring observations
      
    } else if (i == n && i %% 2 == 0) { #make last entry same as previous average
      y_odd[i] = y_odd[i-2] 
    }
  }
  
  return(y_odd)
} 

cv_y_odd = crossval_odd(y)

#Even observations
#does the same thing as crossval_odd, with even and odd reversed
crossval_even = function(y) { # y is the list of observations
  
  n = length(y)
  y_even = vector(mode = "numeric", length = n) # list of even indexed observations from y
  
  for (i in 1:n) {
    if (i == 1) {
      y_even[i] = (y[2] + y[4]) / 2 #fills in first entry, makes same as third
    }
    else {
      if (i %% 2 != 1) { # even indexes
        y_even[i] = y[i]
      } 
      else if (i %% 2 == 1 && i != n) { # even indexes
        y_even[i] = (y[i-1] + y[i+1]) / 2 # fill in odd indexes with average of neighboring observations
        
      } 
      #else if (i == n && i %% 2 == 1) { #if last entry and odd length, make last entry an average of the first and second to last observation
      # y_even[i] = (y[1] + y[i-1])/2
      #} #this won't ever enter because in our case y is always even length
    }
  }
  
  return(y_even)
} 


#HOW MANY LAMBDAS?
# n = sample size = 2^l
# lambda grid = {1, 2^1, 2^2, ... , 2^log(n)}
# so number of lambdas will be log(n) (rounded down to nearest int) + 1

#gets lambdas by powers of 2 until log(n)
get_lambdas = function(y) {
  lambdas = c(0)
  for (i in 0:log(length(y), base = 2)) {
    lambdas[i+1] = 2^i 
  }
  return(lambdas)
}


#function that spits out theta values for each lambda
#will return a vector for each lambda, calls like "fitteddatapoint = theta_vector[[lambda]][datapoint]"
#or can get entire thetahat with thetahat = theta_vector[lambda]
create_theta_vector = function(l, y) {
  # y can be either y_even or y_odd
  n = 2^l
  theta_vector = vector(mode = "numeric", length = length(lambdas))
  
  for (i in 1:length(lambdas)) {
    theta_vector[i] = as.vector(dyadic_1d(l, y, lambdas[i])[2]) #creates vector of theta values
  } 
  
  return(theta_vector)
}

#function that minimizes prediction error, returns final fit after combining even and odd
minimize_pe = function(y, l) { 
  #initiliazations
  pe_odd = c(0)
  pe_even = c(0)
  min_index = 1
  best_lambda_odd = 1
  fit_even = c(0)
  min_index_even = 1
  best_lambda_even = 1
  fit_odd = c(0)
  final_fit_odd = c(0)
  final_fit_even = c(0)
  final_fit = c(0)
  
  
  #pe_odd uses even observations and vice versa
  for(i in 1:length(lambdas)) { #length is same as #of lambdas
    pe_odd[i] = sum((y[seq(2,length(y),2)] - theta_hat_odd[[i]][seq(2,length(y),2)])^2) #sums squared difference of even observations of y and even observations of thetahat
    #y_even = y[c(TRUE, FALSE)] #even observations
    #y_even = y_even[i]
    #pe_odd[i] = sum(y_even^2)
  }
  
  for(i in 1:length(lambdas)) {
    pe_even[i] = sum((y[seq(1,length(y),2)] - theta_hat_even[[i]][seq(1,length(y),2)])^2) #same as a above but reversed
    #y_odd = y[c(FALSE, TRUE)] #odd observations
    #y_odd = y_odd[i]
    #pe_even[i] = sum(y_odd^2)
  }
  
  min_index_even = which.min(pe_even) # returns index of smallest prediction error
  best_lambda_even = lambdas[min_index_even] # returns lambda which has the smallest error
  fit_even = theta_hat_even[[min_index_even]] # final fit for even observations
  
  # repeat process with odd and even switched
  
  
  min_index_odd = which.min(pe_odd) # returns index of smallest error
  best_lambda_odd = lambdas[min_index_odd] # lambda which has the smallest error
  fit_odd = theta_hat_odd[[min_index_odd]] # final fit for odd observations
  
  # now just combine odd and even for final fit
  
  final_fit_odd = fit_odd[seq(1,length(fit_odd),2)] #just the odd obsv
  final_fit_even = fit_even[seq(2,length(fit_even),2)] #just the even obsv
  final_fit = c(rbind(final_fit_odd, final_fit_even)) #combining even and odd observations
  
  return(final_fit)
}

##Run the Algorithm----
#l = 10
#n = 2^l
#sigma = 0.4
#theta = sapply(seq(1:n)/n,f4)
#y = theta + rnorm(2^l,0,sigma); plot(y)

#lambdas = get_lambdas(y); #lambdas
lambdas = c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3, 4, 5, 6, 7, 8, 9)
cv_y_odd = crossval_odd(y);# cv_y_odd
cv_y_even = crossval_even(y);# cv_y_even
theta_hat_even = create_theta_vector(l, cv_y_even);# theta_hat_even
theta_hat_odd = create_theta_vector(l, cv_y_odd);# theta_hat_odd
best_fit = minimize_pe(y,l);# best_fit
###plotting----
plot(x, y, main = "Best Fit mapped onto Y") #original function is black
lines(seq(1,n,1),best_fit, type = "p", col = "red") #fit is red


####Estimating CDFs----

make_t_grid = function(y) {
  t_grid = vector(length = length(y))
  y_sorted = sort(y, decreasing = FALSE)
  
  for (i in 1:length(y)) {
    t_grid[i] = y_sorted[i]
  }
  
  return(t_grid)
}
make_new_data = function(y, t_grid) { #makes new vector, 1 if y <= t, 0 if not
  
  w = list()
  temp_y = c(0)
  
  for (t in 1:length(t_grid)) {
    for (i in 1:length(y)) {
      if (y[i] <= t_grid[t]) {
        temp_y[i] = 1
      } else {
        temp_y[i] = 0
      }
      
      
    }
    w[[t]] = (temp_y)
  }
  return(w)
}

two_normals = function(f1, f2, f3, f4) { # returns y generated from 2 normal dists
  
  coin = c(TRUE, FALSE)
  y = c(0)
  x = runif(n, 0 ,1)
  
  mean1 = sapply(x, f1)
  sigma1 = sapply(x, f2)
  norm1 = rnorm(n, mean1, sigma1)
  
  mean2 = sapply(x, f3)
  sigma2 = sapply(x, f4)
  norm2 = 7+rnorm(n, mean2, sigma2)
  
  actualcdf = (pnorm(norm1, mean1, sigma1) + pnorm(norm2, mean2, sigma2)) / 2
  
  for(i in 1:length(x)) {
    result = sample(coin, size = 1, prob = c(0.5, 0.5))
    if (result == TRUE) {
      y[i] = norm1[i]
    } else {
      y[i] = norm2[i]
    }
  }
  return(y)
}

find_avg_cdf = function(givenx) {
  
  cdf_index = which.min(abs(givenx - x))#finds index in x vector that the given value is closest to
  
  if (cdf_index == 1) {
    averagecdf = (w_matrix[,1]+w_matrix[,2]) / 2
  }
  
  else if (cdf_index == n) {
    averagecdf = (w_matrix[,n]+w_matrix[,n-1]) / 2
  }
  
  else {
    
    averagecdf = (w_matrix[ ,cdf_index + 1] + w_matrix[ ,cdf_index - 1]) / 2
    
  }
  averagecdf = sort(averagecdf)
  return(averagecdf)
  
}

x_big = x
y_big = y
w_matrix_big = w_matrix
find_avg_cdf_big = function(givenx) {
  
  cdf_index = which.min(abs(givenx - x_big))#finds index in x vector that the given value is closest to
  
  if (cdf_index == 1) {
    averagecdf = (w_matrix_big[,1]+w_matrix_big[,2]) / 2
  }
  
  else if (cdf_index == n) {
    averagecdf = (w_matrix_big[,n]+w_matrix_big[,n-1]) / 2
  }
  
  else {
    
    averagecdf = (w_matrix_big[ ,cdf_index + 1] + w_matrix_big[ ,cdf_index - 1]) / 2
    
  }
  averagecdf = sort(averagecdf)
  return(averagecdf)
  
}

#for any given x
random_x = function(anyx) {
  
  averagecdf_small = find_avg_cdf(anyx)
  averagecdf_big = find_avg_cdf_big(anyx)
  plot(t_grid, averagecdf_small, main = paste("cdf at x = ", anyx), xlab = "t values", ylim = c(0,1), type = "l", col = "red")
  lines(t_grid,averagecdf_big, col = "blue")
  #add legend
  mse = mean((averagecdf_big-averagecdf_small)^2)
  message("mse = ", mse)
}

random_x_two_normals = function(anyx) {

  avgcdf = find_avg_cdf(anyx)
  act_cdf = actualcdf # actualcdf defined in line 363
  act_cdf = sort(act_cdf)
 
  plot(t_grid, act_cdf, main = paste("cdf at x = ", anyx), xlab = "t values", ylim = c(0,1), type = "l", col = "blue")
  lines(t_grid, sort(avgcdf), col = "red")
  
  mse = mean((act_cdf-avgcdf)^2)
  message("mse = ", mse)
  
  legend("bottomright", legend=c("original", paste("x = ", anyx)),
         col=c("blue", "red"), lty=1, cex=0.65)
}

threeplots = function(x1, x2, x3) {
  
  avgcdf1 = find_avg_cdf(x1)
  avgcdf1_big = find_avg_cdf_big(x1)
  avgcdf2 = find_avg_cdf(x2)
  avgcdf2_big = find_avg_cdf_big(x2)
  avgcdf3 = find_avg_cdf(x3)
  avgcdf3_big = find_avg_cdf(x3)
  
  plot(t_grid, avgcdf1, main = paste("cdfs at x = ", x1, "x = ", x2, "x = ", x3), xlab = "t", ylab = "P(y<=t)", ylim = c(0,1), type = "l", col = "red", lty = 2)
  lines(t_grid,avgcdf1_big, col = "red")
  lines(t_grid,avgcdf2, col = "blue", lty = 2)
  lines(t_grid,avgcdf2_big, col = "blue")
  lines(t_grid,avgcdf3, col = "green", lty = 2)
  lines(t_grid,avgcdf3_big, col = "green")
  
  
  legend("bottomright", legend=c(paste(x1), paste(x2), paste(x3)),
         col=c("red", "blue", "green"), lty=1, cex=0.65)
  
  message("average mse = ", avg_mse(2))
  
}

plot_int = function(x4, lowerbound, upperbound) {
  
  avgcdf = find_avg_cdf(x4) 
  avgcdf_big = find_avg_cdf_big(x4)
  
  minlower = which.min(abs(avgcdf - lowerbound))
  leftbound = t_grid[minlower]
  
  minupper = which.min(abs(avgcdf - upperbound))
  rightbound = t_grid[minupper]
  
  plot(t_grid, avgcdf, type = "l", col = "blue", main = paste((upperbound-lowerbound)*100, "% Prediction Interval"), lty = 2, xlab = "t", ylab = "P(y<=t)", ylim = c(0,1))
  lines(t_grid, avgcdf_big, col = "blue")
  rect(xleft = leftbound, ybottom = lowerbound, xright = rightbound, ytop = upperbound, border = NA, col = "#FF01AC22")
  message("interval for x = ", x4, ": (", leftbound,",",rightbound,")")
  
}

##to run----
#l = 11
#n = 2^l
#x = runif(n, min = 0, max = 1)


#mean_y = sapply(x, f4)
#sigma_y = sapply(x, f6)
#y = rnorm(2^l,mean_y, sigma_y); plot(y) #make sure if you change the function here
#that you change it everywhere else (random_x, avg_mse, plot_int, threeplots)

#prob_y = sapply(x, f4)
#y = rbinom(n, 10, prob_y)
xldata <- read.csv(file = "actualhd.csv")
heartdisease <- as.matrix(xldata)
l = 11

y = as.double(heartdisease[,10])
x = as.double(heartdisease[,13])
y = y[513:2560]
x = x[513:2560]
n = length(x)
#shape_y = sapply(x, f4)
#scale_y = sapply(x, f2)
#y = rgamma(n, shape_y, scale_y)

y = y[order(x)]
x = x[order(x)]

#fit_cdf = function(y) { #doesnt work as function atm so it's commented out

t_grid = make_t_grid(y) #choose between automated t_grid, or a specific t value or set of t's
#t_grid = c(1.5)

w_matrix = matrix(nrow = length(y), ncol = length(t_grid))
w = make_new_data(y, t_grid)

#lambdas = c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3, 4, 5, 6, 7, 8, 9) #choose between automated lambdas, or a specific set of lambdas
lambdas = get_lambdas(y)

for (t in 1:length(t_grid)) {
  
  cv_w_odd = crossval_odd(w[[t]]); #cv_w_odd
  cv_w_even = crossval_even(w[[t]]); #cv_w_even
  theta_hat_even = create_theta_vector(l, cv_w_even); #theta_hat_even
  theta_hat_odd = create_theta_vector(l, cv_w_odd); #theta_hat_odd
  best_fit = minimize_pe(w[[t]],l); #best_fit
  
  w_matrix[t, ] = best_fit #each row represents yhat based on t_grid[t], #matrix[t, ] = best fit for t'th entry in t_grid, matrix [ ,X] is the cdf of xX, # so matrix [ ,5] is the cdf of x5, based on each t in t_grid
}

for (j in 1:length(x)) {
  
  w_matrix[,j] = sort(w_matrix[,j])
  
}

#return(w_matrix)
#}


plot_t = function(tval) {
  
  binarytvec = c(0)
  
  for (i in 1:length(y)) {
    if (y[i] <= tval) {
      binarytvec[i] = 1
    } else {
      binarytvec[i] = 0
    }
  }
  
  cv_btv_odd = crossval_odd(binarytvec); 
  cv_btv_even = crossval_even(binarytvec); 
  theta_hat_even = create_theta_vector(l, cv_btv_even); 
  theta_hat_odd = create_theta_vector(l, cv_btv_odd); 
  best_fit = minimize_pe(binarytvec,l)
  plot(x, binarytvec, main = paste("t = ", tval), xlab = "x", col = "blue")
  lines(x, best_fit, type = "p", col = "red")
}

avg_mse = function(temp) {
  
  sum = 0
  
  
  for (i in 1:length(x)) {
    
    sum = sum + mean((find_avg_cdf_big(x[i]) - w_matrix[,i])^2)
    
  }
  
  avg_mse = sum/length(x)
  message("average mse = ", avg_mse)
  return(avg_mse)
  
}

# running two normals function
l = 11
n = 2^l
x = runif(n, 0, 1)
f1 = f
f2 = f2
f3 = f3
f4 = f4
y = two_normals(f1,f2,f3,f4)

y = y[order(x)]
x = x[order(x)]

t_grid = make_t_grid(y)

w_matrix = matrix(nrow = length(y), ncol = length(t_grid))
w = make_new_data(y, t_grid)
lambdas = get_lambdas(y)

for (t in 1:length(t_grid)) {
  
  cv_w_odd = crossval_odd(w[[t]]); #cv_w_odd
  cv_w_even = crossval_even(w[[t]]); #cv_w_even
  theta_hat_even = create_theta_vector(l, cv_w_even); #theta_hat_even
  theta_hat_odd = create_theta_vector(l, cv_w_odd); #theta_hat_odd
  best_fit = minimize_pe(w[[t]],l); #best_fit
  
  w_matrix[t, ] = best_fit #each row represents yhat based on t_grid[t], #matrix[t, ] = best fit for t'th entry in t_grid, matrix [ ,X] is the cdf of xX, # so matrix [ ,5] is the cdf of x5, based on each t in t_grid
}

for (j in 1:length(x)) {
  
  w_matrix[,j] = sort(w_matrix[,j])
  
}

random_x_two_normals(0.5)
random_x(32)
plot_t(0) #doesn't work atm, works line by line
plot_int(30, .05, .95)
threeplots(21, 27, 33)
