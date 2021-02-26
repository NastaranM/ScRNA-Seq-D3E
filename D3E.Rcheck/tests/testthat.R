library(testthat)
library(D3E)

data = matrix(c(1,2,3,5,6,8,3,1),nrow = 2)

test_that("geometric_mean",{
  expect_error(geometric_mean(NA))
  expect_error(geometric_mean(NULL))
  expect_condition(geometric_mean("ab"))
  expect_condition(geometric_mean(2))
  expect_condition(geometric_mean(c(2,3)))
  expect_equal(geometric_mean(data) ,c(4.481689,54.598150,1096.633158,7.389056))
})


test_that("weight_normalization",{
  expect_error(weight_normalization(NA))
  expect_error(geometric_mean(NULL))
  expect_condition(geometric_mean("ab"))
  expect_condition(geometric_mean(2))
  expect_condition(geometric_mean(c(2,3)))
  expect_equal(weight_normalization(data), c(0.557825400,0.091578194,0.007295056,0.406005850))
})

test_that("normalized_data",{
  expect_error(normalized_data(NA))
  expect_error(normalized_data(NULL))
  expect_condition(normalized_data("ab"))
  expect_condition(normalized_data(2))
  expect_condition(normalized_data(c(2,3)))
  expect_error(normalized_data(data))
  expect_error(normalized_data(data,data))
  expect_error(normalized_data(data,"a"))
  expect_error(normalized_data(data,c(1,23,4)))
  expect_equal(normalized_data(data,c(1,2,4,5)),
               matrix(c(1,0.75,6,0.75,1,1.00,4,0.20),nrow = 2, byrow = TRUE))
})

test_that("kolmogorov_smirnov_test",{
  expect_equal(kolmogorov_smirnov_test("a"),-1)
  expect_equal(kolmogorov_smirnov_test("a",1),-1)
  expect_equal(kolmogorov_smirnov_test(1),-1)
  expect_equal(kolmogorov_smirnov_test(c(1,2)),-1)

})

test_that("anderson_darling_test",{
  expect_equal(anderson_darling_test("a"),-1)
  expect_equal(anderson_darling_test("a",1),-1)
  expect_equal(anderson_darling_test(1),-1)
  expect_equal(anderson_darling_test(c(1,2)),-1)

})


test_that("cramer_von_masis",{
  expect_equal(cramer_von_masis("a"),-1)
  expect_equal(cramer_von_masis("a",1),-1)
  expect_equal(cramer_von_masis(1),-1)
  expect_equal(cramer_von_masis(c(1,2)),-1)

})

test_that("distribution_test",{
  expect_error(distribution_test("a"))
  expect_error(distribution_test("a",1))
  expect_error(distribution_test(1))
  expect_error(distribution_test(c(1,2)))
  expect_error(distribution_test(NA))
  expect_error(distribution_test(NULL))
  expect_error(distribution_test(Inf))

})

test_that("rand_Poisson_beta",{
  expect_error(rand_Poisson_beta(data))
  expect_error(rand_Poisson_beta(data,1))
  expect_error(rand_Poisson_beta(1))
  expect_error(rand_Poisson_beta(NA))
  expect_error(rand_Poisson_beta(NULL))
  expect_error(rand_Poisson_beta(Inf))
  expect_condition(rand_Poisson_beta())
})

test_that("get_para_moments",{
  expect_error(get_para_moments(NULL))
  expect_condition(get_para_moments())
  expect_error(get_para_moments("a"))
  expect_error(get_para_moments("a",1))
})

test_that("get_para_bayesian",{
  expect_error(get_para_moments(NULL))
  expect_condition(get_para_moments())
  expect_error(get_para_moments("a"))
  expect_error(get_para_moments("a",1))
})


test_that("goodness_of_fit",{
  expect_error(goodness_of_fit(NULL))
  expect_condition(goodness_of_fit())
  expect_error(goodness_of_fit("a"))
  expect_error(goodness_of_fit("a",1))
  expect_error(goodness_of_fit(data))
  expect_error(goodness_of_fit(data,data,data,data))
  expect_error(goodness_of_fit(data,2,"A"))
  expect_error(goodness_of_fit(data,2,4))
  expect_error(goodness_of_fit(data,2,3,10,1))

})

test_that("log_likelihood",{
  expect_error(log_likelihood(NULL))
  expect_condition(log_likelihood())
  expect_error(log_likelihood("a"))
  expect_error(log_likelihood("a",1))
  expect_error(log_likelihood(data))
  expect_error(log_likelihood(data,data,data,data))
  expect_error(log_likelihood(data,2,"A"))
  expect_error(log_likelihood(data,2,4))
  expect_error(log_likelihood(2))
  expect_error(log_likelihood(data,2,3,10,1))
})

test_that("likelihood_ratio",{
  expect_error(likelihood_ratio(NULL))
  expect_condition(likelihood_ratio())
  expect_error(likelihood_ratio("a"))
  expect_error(likelihood_ratio("a",1))
  expect_error(likelihood_ratio(data))
  expect_error(likelihood_ratio(data,data,data,data))
  expect_error(likelihood_ratio(data,2,"A"))
  expect_error(likelihood_ratio(data,2,4))
  expect_error(likelihood_ratio(2))
  expect_error(likelihood_ratio(data,2,3,10,1))
  expect_error(likelihood_ratio(0))
})

test_that("control_split",{
  expect_error(control_split(5))
  expect_condition(control_split())
  expect_error(control_split(c(1,3,5),"q"))
  expect_error(control_split(data,1))
  expect_error(control_split(NA))
  expect_error(control_split(Inf))
  expect_error(control_split(NULL))

})

test_that("comparison_data",{
  expect_error(comparison_data(5))
  expect_condition(comparison_data())
  expect_error(comparison_data(c(1,3,5),"q"))
  expect_error(comparison_data(NA))
  expect_error(comparison_data(Inf))
  expect_error(comparison_data(NULL))
})

test_that("finding_threshold",{
  expect_equal(finding_threshold(data,c(0.5,0.5)),list("pvalues" = c(0.5,0.5), "significant" = 0.5))
  expect_error(finding_threshold(data))
  expect_error(finding_threshold(NA))
  expect_error(finding_threshold(NULL))
  expect_error(finding_threshold(Inf))
  expect_error(finding_threshold())
  expect_error(finding_threshold(2,2))
})

test_that("check_cramer_von_mises",{
  expect_equal(check_cramer_von_mises(-1,"hi"),"hi")
  expect_equal(check_cramer_von_mises("hi"),NULL)
})

test_that("parameter_estimation",{
  expect_error(parameter_estimation(1))
  expect_error(parameter_estimation(1,2))
  expect_error(parameter_estimation(data,1,1))
  expect_error(parameter_estimation(NA))
  expect_error(parameter_estimation(NULL))
  expect_error(parameter_estimation(Inf))
  expect_error(parameter_estimation(data,data,2,1))
  expect_error(parameter_estimation(data,data,2,5))
})
