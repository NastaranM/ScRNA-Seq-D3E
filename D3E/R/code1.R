rm(list = ls())

setwd("~/LiUMaterial/ResearchProject/NSCserver")
require(tidyverse)
require(ggplot2)
require(tibble)
require(magrittr)
require(twosamples)
require(matrixStats)

##########################################################################

params = list("alpha" = 0, "beta" = 0, "gamma" = 0, "c" = 0)
bio_params = list("size" = NULL,"freq"=NULL,"duty"=NULL)


##########################################################################



##########################################################################

set_sample_function = function(func){
  sample_function = func

  draw = function(value,p,max_steps = 1000, save_to_sample = FALSE){
    left_limit = 0
    right_limit = Inf
    x0 = value
    sample = c()
    w = abs(value/2)
    f = sample_function

    logPx = f(x0,p)
    log_slice = log(runif(1)) + logPx
    x_left = x0 - runif(1) * w
    x_right = x_left + w
    if (x_left < left_limit){
      x_left = 0
    }
    if(x_right > right_limit){
      x_right = right_limit
    }
    v = runif(1)
    j = floor(max_steps * v)
    k = max_steps - 1- j

    while (j > 0 & log_slice < f(x_left,p) & x_left - w > left_limit){
      j = j-1
      x_left = x_left - w
    }

    while (k > 0 & log_slice < f(x_right,p) & x_right + w > right_limit){
      k = k -1
      x_right = x_right + w
    }

    n = 10000
    while(1){
      n = n-1
      if(n < 1){
        print("Warning: Cannot find a new value.")
        return(x0)
      }
      x1 = (x_right - x_left) * runif(1) + x_left
      if(log_slice <= f(x1,p)){
        break
      }
      if(x1 < x0){
        x_left = x1
      }
      else{
        x_right  = x1
      }
    }
    value = x1
    if(save_to_sample){
      sample = append(sample,x1)
    }
    return(list("x1" =x1, "sample" = sample, "value" = value))
  }
}

mean_ = function(sample){
  if(length(sample)>0){
  result = mean(sample)}else{
    result = NULL
  }
  return(result)
}

#########################################################################
hyper_para = list("k_alpha" = 1, "theta_alpha" = 100,
                  "k_beta" = 1, "theta_beta" = 100,
                  "k_gamma" = 1)
fc = function(x,p) (params$alpha -1)*log(x)+(params$beta - 1) + pi * log(x) - params$gamma*x
fgamma = function(x,p) (hyper_para$k_gamma-1)*log(x)-
  x/max(p)+log(x)*sum(p) - x*sum(params$c)
falpha = function(x,p) (hyper_para$k_alpha-1)*log(x) -
  x/ hyper_para$theta_alpha+length(p)*log(gamma(x+params$beta)) -
  log(gamma(x))*sum(log(params$c))
fbeta = function(x,p)  (hyper_para$k_beta-1)*log(x) -
  x/ hyper_para$theta_beta+length(p)*log(gamma(x+params$alpha)) -
  log(gamma(x))*sum(log(1-params$c))

##########################################################################


# normalization

geometric_mean = function(data){
  try(if (any(is.na(data)))
    stop("Please remove NA values"))
  try(if(is.null(dim(data)))
  stop("The dimensions of the data is not correct"))
  n = dim(data)[1]
  m = dim(data)[2]
  geomean = exp(colMeans(data))
  return(geomean)
}



weight_normalization = function(data){
  try(if(is.null(dim(data)))
  stop("The dimensions of the data is not correct"))
  try(if (any(is.na(data)))
    stop("Please remove NA values"))
  weight = colMedians(data+1)/geometric_mean(data)
  return(weight)}

normalized_data = function(data, weights){
  try(if(is.null(dim(data)))
    stop("The dimensions of the data is not correct"))
  try(if(length(weights)==0)
    stop("Input length is not correct"))
  try(if(any(is.na(data)))
    stop("Please remove NA values"))
  try(if(any(is.na(weights)))
    stop("Please remove NA values"))
  checking = which(weights == 0)
  if(length(checking > 0)){
    weights = weights[-checking]
    newdata = data[,-checking]}else{
      newdata = data}
  stopifnot(dim(newdata)[2] == length(weights))
  normdata = newdata/weights
  return(normdata)
    }


##########################################################################


kolmogorov_smirnov_test = function(x,y){
  result = try(ks_test(x,y),silent = TRUE)
  if(class(result) == "try-error"){
    result = c(-1,-1)
    }
  return(result[2])
}

anderson_darling_test = function(x,y){
  result = try(ad_test(x,y),silent = TRUE)
  if(class(result) == "try-error"){
    result = c(-1,-1)
  }
  return(result[2])
}
cramer_von_masis = function(x,y){
  result = try(cvm_test(x,y),silent = TRUE)
  if(class(result) == "try-error"){
    result = c(-1,-1)
  }
  return(result[2])
}

##########################################################################


distribution_test = function(x,y,method){
  if(method == 1){
    return(kolmogorov_smirnov_test(x,y))
  }
  else if(method == 2){
    return(anderson_darling_test(x,y))
  }
  else if(method == 3){
    return(0)}
  else{return(cramer_von_masis(x,y))}
}

##########################################################################

rand_Poisson_beta = function(params,n){
  x = rbeta(n,params$alpha,params$beta)
  p = rpois(1,x * params$gamma)
  return(p)
}

##########################################################################

get_para_moments = function(p){

  rm1 = sum(p) / length(p)
  for(x in p){
    rm2 = sum(x^2)/length(p)
    rm3 = sum(x^3)/length(p)}

  fm1 = rm1
  fm2 = rm2 - rm1
  fm3 = rm3 - 3*rm2 + 2*rm1

  r1 = fm1
  r2 = fm2 / fm1
  r3 = fm3 / fm2

  alpha1 = 2*r1 * (r3 - r2) / (r1*r2 - 2*r1*r3 + r2*r3)
  beta1 = 2 * (r2 - r1) * (r1 - r3) * (r3 - r2) / ((r1*r2 - 2*r1*r3 + r2*r3) *(r1 - 2*r2 + r3))
  gamma1 = (-r1*r2 + 2*r1*r3 - r2*r3) / (r1 - 2*r2 +r3)

  params = list("alpha" = alpha1,"beta"= beta1, "gamma" = gamma1, "c" = 0)
  for(elem in params){
  if(is.numeric(elem) & !is.infinite(elem) & !is.nan(elem)){
    params = params
      }
  else{params = list("alpha" = -1,"beta"= -1, "gamma" = -1, "c" = -1)}
  }
return(params)
}

##########################################################################

get_para_bayesion = function(p, nIter = 1000){
  para_fit = get_para_moments(p)

  # hyperparameters

  hyper_para = list("k_alpha" = 1, "theta_alpha" = 100,
                 "k_beta" = 1, "theta_beta" = 100,
                 "k_gamma" = 1, "theta_gamma" = max(p))

  if(para_fit$alpha > 0 & para_fit$beta > 0 & para_fit$gamma > 0){
    cvalues = list()
    resc = list("c" = NULL, "c_sample" = NULL, "c_value" = NULL)
    func1 = set_sample_function(fc)
    for(i in 1:length(p)){
      cvalues[[i]] = c(func1(0.5,p))
    }
    c1 = c()
    c_s = c()
    c_v = c()
    for(i in 1:length(p)){
      c1[i] = cvalues[[i]]$x1
      c_s[i] = cvalues[[i]]$sample
      c_v[i] = cvalues[[i]]$value
    }
    resc$c = c1
    resc$c_sample = c_s
    resc$c_value = c_v
    c_mean = mean_(resc$c_sample)

    names(resc$x1) = NULL
    names(resc$sample) = NULL
    names(resc$value) = NULL


    resgamma = list("gamma" = NULL, "gamma_sample" = NULL, "gamma_value" = NULL)
    funcgamma = set_sample_function(fgamma)
    resgamma = funcgamma(para_fit$gamma,p,max_steps = 1000)
    gamma_mean= mean_(resgamma$sample)


    names(resgamma$x1) = NULL
    names(resgamma$sample) = NULL
    names(resgamma$value) = NULL


    resalpha = list("alpha" = NULL, "alpha_sample" = NULL, "alpha_value" = NULL)
    funcalpha = set_sample_function(falpha)
    resalpha = funcalpha(para_fit$alpha,p,max_steps = 1000)
    alpha_mean= mean_(resalpha$sample)

    names(resalpha$x1)=NULL
    names(resalpha$sample) = NULL
    names(resalpha$value) = NULL


    resbeta = list("beta" = NULL, "beta_sample" = NULL, "beta_value" = NULL)
    funcbeta = set_sample_function(fbeta)
    resbeta = funcbeta(para_fit$beta,p,max_steps = 1000)
    beta_mean = mean_(resbeta$sample)

    names(resbeta$x1)=NULL
    names(resbeta$sample) = NULL
    names(resbeta$value) = NULL

    params = list("alpha" = resalpha$x1, "beta" = resbeta$x1,"gamma" = resgamma$x1, "c" = resc$c )
    mean_params = list("alpha" = alpha_mean, "beta" = beta_mean,"gamma" = gamma_mean, "c" = c_mean)
    sample_params = list("alpha" = resalpha$sample, "beta" = resbeta$sample,"gamma" = resgamma$sample, "c" = resc$c_sample)
    }else{
      cvalues = list()
      resc = list("c" = NULL, "c_sample" = NULL, "c_value" = NULL)
      func1 = set_sample_function(fc)
      for(i in 1:length(p)){
        cvalues[[i]] = c(func1(0.5,p))
      }
      c1 = c()
      c_s = c()
      c_v = c()
      for(i in 1:length(p)){
        c1[i] = cvalues[[i]]$x1
        c_s[i] = cvalues[[i]]$sample
        c_v[i] = cvalues[[i]]$value
      }
      resc$c = c1
      resc$c_sample = c_s
      resc$c_value = c_v
      c_mean = mean_(resc$sample)

      names(resc$x1) = NULL
      names(resc$sample) = NULL
      names(resc$value) = NULL


      resgamma = list("gamma" = NULL, "gamma_sample" = NULL, "gamma_value" = NULL)
      funcgamma = set_sample_function(fgamma)
      xg = mean(p)+1e6
      resgamma = funcgamma(xg,p,max_steps = 1000)
      gamma_mean= mean_(resgamma$sample)


      names(resgamma$x1) = NULL
      names(resgamma$sample) = NULL
      names(resgamma$value) = NULL



      resalpha = list("alpha" = NULL, "alpha_sample" = NULL, "alpha_value" = NULL)
      funcalpha = set_sample_function(falpha)
      resalpha = funcalpha(0.5,p,max_steps = 1000)
      alpha_mean= mean_(resalpha$sample)

      names(resalpha$x1)=NULL
      names(resalpha$sample) = NULL
      names(resalpha$value) = NULL


      resbeta = list("beta" = NULL, "beta_mean"= NULL,"beta_sample" = NULL, "beta_value" = NULL)
      funcbeta = set_sample_function(fbeta)
      resbeta = funcbeta(0.5,p,max_steps = 1000)
      beta_mean = mean_(resbeta$sample)

      names(resbeta$x1)=NULL
      names(resbeta$sample) = NULL
      names(resbeta$value) = NULL

      params = list("alpha" = resalpha$x1, "beta" = resbeta$x1,"gamma" = resgamma$x1, "c" = resc$c )
      mean_params = list("alpha" = alpha_mean, "beta" = beta_mean,"gamma" = gamma_mean, "c" = c_mean)
      sample_params = list("alpha" = resalpha$sample, "beta" = resbeta$sample,"gamma" = resgamma$sample, "c" = resc$c_sample)
    }
  bio_params = list("size" = NULL,"freq"=NULL,"duty"=NULL)
  save = FALSE
  for(i in 1:nIter){
    if(i > floor(nIter/2)){
      save = TRUE
      bio_params = list("size" = NULL,"freq"=NULL,"duty"=NULL)
      bio_params$size = append( bio_params$size,params$gamma/params$beta )
      bio_params$freq = append(bio_params$freq, params$alpha*params$beta/(params$alpha+params$beta))
      bio_params$duty = append(bio_params$duty,params$alpha/(params$alpha+params$beta))
    }
      cvalues = list()
      resc = list("c" = NULL, "c_sample" = NULL, "c_value" = NULL)
      func1 = set_sample_function(fc)
      for(i in 1:length(p)){
        cvalues[[i]] = c(func1(params$c[i],p,save_to_sample = save))
      }
      c1 = c()
      c_s = c()
      c_v = c()
      for(i in 1:length(p)){
        c1[i] = cvalues[[i]]$x1
        c_s[i] = cvalues[[i]]$sample
        c_v[i] = cvalues[[i]]$value
      }
      resc$c = c1
      resc$c_sample = c_s
      resc$c_value = c_v
      c_mean = mean_(resc$c_sample)

      names(resc$x1) = NULL
      names(resc$sample) = NULL
      names(resc$value) = NULL




       resgamma = list("gamma" = NULL,"gamma_sample" = NULL, "gamma_value" = NULL)
       funcgamma = set_sample_function(fgamma)
       resgamma = funcgamma(params$gamma,p,max_steps = 1000,save_to_sample = save)
       gamma_mean = mean_(resgamma$sample)

       names(resgamma$x1)=NULL
       names(resgamma$sample) = NULL
       names(resgamma$value) = NULL


       resalpha = list("alpha" = NULL,"alpha_sample" = NULL, "alpha_value" = NULL)
       funcalpha = set_sample_function(falpha)
       resalpha = funcalpha(params$alpha,p,max_steps = 1000,save_to_sample = save)
       alpha_mean = mean_(resalpha$sample)

       names(resalpha$x1)=NULL
       names(resalpha$sample) = NULL
       names(resalpha$value) = NULL


       resbeta = list("beta" = NULL, "beta_sample" = NULL, "beta_value" = NULL)
       funcbeta = set_sample_function(fbeta)
       resbeta = funcbeta(params$beta,p,max_steps = 1000,save_to_sample = save)
       beta_mean = mean_(resbeta$sample)

       names(resbeta$x1)=NULL
       names(resbeta$sample) = NULL
       names(resbeta$value) = NULL

       params = list("alpha" = resalpha$x1, "beta" = resbeta$x1,"gamma" = resgamma$x1, "c" = resc$c )
       mean_params = list("alpha" = alpha_mean, "beta" = beta_mean,"gamma" = gamma_mean, "c" = c_mean)
       sample_params = list("alpha" = resalpha$sample, "beta" = resbeta$sample,"gamma" = resgamma$sample, "c" = resc$c_sample)

  }
  return(list("params"=params, "mean_params" = mean_params,"sample_params" = sample_params,"bio_params" = bio_params))
  }



##########################################################################

goodness_of_fit = function(p, params, mean_params){
  if(!is.null(mean_params$alpha) & !is.null(mean_params$beta) & !is.null(mean_params$gamma)){
    alpha1 = mean_params$alpha
    beta1 = mean_params$beta
    gamma1 = mean_params$gamma
  }else{
    alpha1 = params$alpha
    beta1 = params$beta
    gamma1 = params$gamma
  }
  params = list("alpha" = alpha1, "beta" = beta1,"gamma" = gamma1, "c" = 0)
  pr = rand_Poisson_beta(params = params,n = length(p))
  return(cramer_von_masis(pr,p))
}

##########################################################################

# Find log of likelihood for sample 'p' with parameters 'params'
# doing poisson-beta sampling 'n' times


log_likelihood = function(p,params,mean_params,n){
  if(!is.null(mean_params$alpha) & !is.null(mean_params$beta) & !is.null(mean_params$gamma)){
    alpha1 = mean_params$alpha
    beta1 = mean_params$beta
    gamma1 = round(mean_params$gamma)+1
  }else{
    if(!is.null(params$alpha) & !is.null(params$beta) & !is.null(params$gamma)){
    alpha1 = params$alpha
    beta1 = params$beta
    gamma1 = round(params$gamma)+1}else{
      return(cat("cannot estimate the parameter"))
    }
  }

  if(n * params$gamma> 1e9){
    n = round(1e9/params$gamma)
  cat("Status1","\n","Reduced Poisson-Beta sample to ", n," ", ".")}
  pval = c()
  for(item in p){
    x = rbeta(n,params$alpha, params$beta)
    pTemp = 0
    for(i in 1:n){
      pTemp = pTemp + rpois(item, params$gamma*x[i])
    }
    pval = append(pval, pTemp/n)
  }
  pval = as.array(pval)
  pval[which(pval %in% 0)] = 1/n
  return(sum(log(pval)))
}

##########################################################################


# Perform likelihood ratio test fitted parameters

likelihood_ratio = function(p, params1, params2, mean_params1, mean_params2,n){
  if(!is.null(mean_params1$alpha) & !is.null(mean_params1$beta) & !is.null(mean_params1$gamma)){
    alpha1 = mean_params1$alpha
    beta1 = mean_params1$beta
    gamma1 = round(mean_params1$gamma)+1
  }else{
    alpha1 = params1$alpha
    beta1 = params1$beta
    gamma1 = round(params1$gamma)+1
  }

  if(!is.null(mean_params2$alpha) & !is.null(mean_params2$beta) & !is.null(mean_params2$gamma)){
    alpha2 = mean_params2$alpha
    beta2 = mean_params2$beta
    gamma2 = round(mean_params2$gamma)+1
  }else{
    alpha2 = params2$alpha
    beta2 = params2$beta
    gamma2 = round(params2$gamma)+1
  }

  for(x in round(p)){
    sum1 = log_likelihood(x,params1,mean_params1,n)
    sum2 = log_likelihood(x,params2,mean_params2,n)
  }
  ratio = sum1 - sum2
  if(ratio >=0){
    cat("status1","\n","Could not perform likelihood ratio test.")
    return(NaN)
  }
  pval = dchisq(-2*ratio,3)
  return(pval)
}


##########################################################################



control_split = function(normalized_dataset, cellname){
  try(if(any(is.na(normalized_dataset)))
    stop("Please remove NA values"))
  try(if(is.null(normalized_dataset))
    stop("The dimensions of the data is not correct"))
  try(if(any(is.na(colnames(normalized_dataset))))
             stop("Please enter column names"))
  try(if(any(is.na(cellname)))
    stop("Please provide a valid cellname"))
  name = colnames(normalized_dataset)
  type_of_cells = name[grep(cellname, name)]
  try(if(length(type_of_cells) == 0)
    stop("Please enter a cellname from the type of cells"))
  selected_cell = colnames(normalized_dataset) == type_of_cells
  selection = normalized_dataset[,selected_cell, drop = FALSE]
  n = dim(selection)[2]
  m = floor(n/2)
  data1 = selection[,1:m]
  colnames(data1) = type_of_cells[1:m]
  rownames(data1) = rownames(normalized_dataset)

  data2 = selection[,(m+1):n]
  colnames(data2) = type_of_cells[(m+1):n]
  rownames(data2) = rownames(normalized_dataset)

  return(list(x = data1,y = data2))
}



comparison_data = function(normalized_dataset, cellname){
  try(if(any(is.na(normalized_dataset)))
    stop("Please remove NA values"))
  try(if(is.null(normalized_dataset))
    stop("The dimensions of the data is not correct"))
  try(if(any(is.na(colnames(normalized_dataset))))
    stop("Please enter column names"))
  try(if(any(is.na(cellname)))
    stop("Please provide a valid cellname"))
  name = colnames(normalized_dataset)
  type_of_cells = name[grep(cellname, name)]
  try(if(length(type_of_cells) == 0)
    stop("Please enter a cellname from the type of cells"))
  selected_cell = colnames(normalized_dataset) == type_of_cells
  selection = normalized_dataset[,selected_cell, drop = FALSE]
  colnames(selection) = type_of_cells
  rownames(selection) = rownames(normalized_dataset)
  return(selection)}

finding_threshold = function(x,pvals, fdr = 0.01){
  names(pvals) = rownames(x)
  pvalues = sort(pvals,decreasing = TRUE)
  n = length(pvalues)
  keeper = c()
  if(n>1){
    for(i in 2:n+1){
      if(pvalues[i-1] <= fdr * i/n)
        cat("Significance threshold with given fdr is ", pvalues[i-1])
      keeper = append(keeper, pvalues[i-1])
    }
    return(list("pvalues" = pvalues, "significant" = keeper))
  }
  else{stop("Please enter p-values")}
}


