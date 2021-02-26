
check_cramer_von_mises = function(pvalue, comment){
  if(pvalue == -1){
    return(comment)
  }
}

parameter_estimation = function(p1,p2,mod, method,n){
  try(if(!(mod %in% c(0,1)))stop("select a mode of analysis from 0 and 1"))
  try(if (!(method %in% c(0,1,2,3)))stop("select a method between 0 to 3 for distributio test"))

  result_para = list("gof1" = NULL, "gof2" = NULL)
  difference = distribution_test(p1,p2,method)
  if(difference == -1){
    cat("Could not estimate p-value. Further analysis aborted.")
  }
  if(mod== 0){
    params1 = get_para_moments(p1)
    params2 = get_para_moments(p2)

    if(params1$alpha <= 0 | params1$beta <= 0 | params1$gamma <= 0 |params2$alpha <= 0 | params2$beta <= 0 | params2$gamma <= 0){
      cat("Could not estimate parameters.")
      params1 = list("alpha" = NULL, "beta" = NULL,"gamma" = NULL)
      params2 = list("alpha" = NULL, "beta" = NULL,"gamma" = NULL)
      mean_params1 = list("alpha" = NULL, "beta" = NULL, "gamma" = NULL)
      mean_params2 = list("alpha" = NULL, "beta" = NULL, "gamma" = NULL)
      bio_params1 = list("size" = NULL,"freq"=NULL,"duty"=NULL)
      bio_params2 = list("size" = NULL,"freq"=NULL,"duty"=NULL)
    }else{

      bio_params1 = list("size" = NULL,"freq"=NULL,"duty"=NULL)
      bio_params2 = list("size" = NULL,"freq"=NULL,"duty"=NULL)
      bio_params1$size = params1$gamma/params1$beta
      bio_params1$freq = 1/params1$alpha
      bio_params1$duty = params1$alpha/(params1$alpha+params1$beta)

      bio_params2$size = params2$gamma/params2$beta
      bio_params2$freq = 1/params2$alpha
      bio_params2$duty = params2$alpha/(params2$alpha+params2$beta)

      if(method == 3){
        difference = likelihood_ratio(p = p2,params1 = params1,
                                      params2 = params2,mean_params1 = mean_params1,
                                      mean_params2 = mean_params2,n)

        if(difference == -1){
          cat("Could not estimate p-value. Further analysis aborted.")
        }
        gof1 = goodness_of_fit(p = p1,params = params1,mean_params = mean_params1)
        check_cramer_von_mises(gof1,comment = 'Goodnes of fit for a first cell type')

        gof2 = goodness_of_fit(p = p2,params = params2,mean_params = mean_params2)
        check_cramer_von_mises(gof1,comment = 'Goodnes of fit for a second cell type')
        result_para = list("gof1" = gof1, "gof2" = gof2)
      }
    }
  }else if(mod == 1){
    params1 = get_para_bayesion(p1)
    params2 = get_para_bayesion(p2)
    params1$alpha = params1$params$alpha
    params1$beta = params1$params$beta
    params1$gamma = params1$params$gamma

    mean_params1 = list("alpha" = NULL, "beta" = NULL, "gamma" = NULL)
    mean_params2 = list("alpha" = NULL, "beta" = NULL, "gamma" = NULL)

    mean_params1$alpha = params1$mean_params$alpha
    mean_params1$beta = params1$mean_params$beta
    mean_params1$gamma = params1$mean_params$gamma

    params2$alpha = params2$params$alpha
    params2$beta = params2$params$beta
    params2$gamma = params2$params$gamma

    mean_params2$alpha = params2$mean_params$alpha
    mean_params2$beta = params2$mean_params$beta
    mean_params2$gamma = params2$mean_params$gamma

    bio_params1 = list("size" = NULL, "freq" = NULL, "duty" = NULL)
    bio_params2 = list("size" = NULL, "freq" = NULL, "duty" = NULL)

    bio_params1$size = params1$bio_params$size
    bio_params1$freq = params1$bio_params$freq
    bio_params1$duty = params1$bio_params$duty

    bio_params2$size = params2$bio_params$size
    bio_params2$freq = params2$bio_params$freq
    bio_params2$duty = params2$bio_params$duty

    bio_params1_mean = list("size" = mean(bio_params1$size),
                            "freq" = mean(bio_params1$freq),
                            "duty" = (bio_params1$duty))
    bio_params2_mean = list("size" = mean(bio_params2$size),
                            "freq" = mean(bio_params2$freq))

    size_p = cramer_von_masis(bio_params1$size, bio_params2$size)
    check_cramer_von_mises(pvalue = size_p,comment = "for size p")

    freq_p = cramer_von_masis(bio_params1$freq, bio_params2$freq)
    check_cramer_von_mises(freq_p, comment = "for freq p")

    duty_p = cramer_von_masis(bio_params1$duty, bio_params2$duty)
    check_cramer_von_mises(duty_p,comment = "for duty p")

    if(method == 3){
      difference = likelihood_ratio(p = p2,params1 = params1,
                                    params2 = params2,mean_params1 = mean_params1,
                                    mean_params2 = mean_params2,n)

      if(difference == -1){
        cat("Could not estimate p-value. Further analysis aborted.")
      }
      pvalue = append(pvalue,difference)
      gof1 = goodness_of_fit(p = p1,params = params1,mean_params = mean_params1)
      check_cramer_von_mises(gof1,comment = 'Goodnes of fit for a first cell type')

      gof2 = goodness_of_fit(p = p2,params = params2,mean_params = mean_params2)
      check_cramer_von_mises(gof1,comment = 'Goodnes of fit for a second cell type')
      result_para = list("gof1" = gof1, "gof2" = gof2)
    }
  }


  if(mod == 0){
    result = list('params1' = params1,
                  'params2' = params2,
                  "bio_params1" = bio_params1,
                  "bio_params2" = bio_params2,
                  "log_size_ratio" = log2(bio_params1$size/bio_params2$size),
                  "log_size_freq" = log2(bio_params1$freq/bio_params2$freq),
                  "log_size_duty" = log2(bio_params1$duty/bio_params2$duty),
                  "para-est" = result_para,
                  "pvalue" = difference,"mean_p1" = mean(p1),
                  "variation_p1" = var(p1),"mean_p2" = mean(p2),
                  "variation_p2" = var(p2))}else if(mod == 1){
                    result = list('params1_mean' = mean_params1,
                                  'params2_mean' =  mean_params2,
                                  "bio_params1_mean" = list("size" = mean(bio_params1$size),
                                                            "freq" = mean(bio_params1$freq),
                                                            "duty" = (bio_params1$duty)),
                                  "bio_params2_mean" = list("size" = mean(bio_params2$size),
                                                            "freq" = mean(bio_params2$freq),
                                                            "duty" = (bio_params2$duty)),
                                  "log_size_ratio" = log2(bio_params1_mean$size/bio_params2_mean$size),
                                  "log_size_freq" = log2(bio_params1_mean$freq/bio_params2_mean$freq),
                                  "log_size_duty" = log2(bio_params1_mean$duty/bio_params2_mean$duty),
                                  "para-est" = result_para,
                                  "pvalue" = difference,"size" = size_p,
                                  "t" = duty_p, "f" = freq_p,"mean_p1" = mean(p1),
                                  "variation_p1" = var(p1),"mean_p2" = mean(p2),
                                  "variation_p2" = var(p2))}
  return(result)
}


analysis = function(FUN = parameter_estimation,x,y,mod,method, n=2){
  benchmark = list(NULL)

  for(i in 1:nrow(x)){
    print(i)
    benchmark[[i]] = parameter_estimation(x[i,], y[i,],mod = mod,method = method, n = n)
  }
  return(benchmark)
}

threshold = function(benchmark,x,fdr = 1){
  pvals = c()
  for(i in 1:nrow(x)){
    pvals[i] = benchmark[[i]]$pvalue

  }
  names(pvals) = rownames(x)
  pvalues = sort(pvals,decreasing = TRUE)
  n = length(pvalues)
  keeper = c()
  if(n>1){
    for(i in 2:n+1){
      print(i)
      if(pvalues[i-1] <= fdr * i/n)
        cat("Significance threshold with given fdr is ", pvalues[i-1])
      keeper = append(keeper, pvalues[i-1])
    }}
  return(list("pvalues" = pvalues, "significant" = keeper))
}



