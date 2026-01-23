local_edition(3)

piecewise_linearise <- function(x, y, m_left = 0, m_right = 0) {
  function(z) {
    n <- length(x)
    
    if (z <= x[1]) return(y[1] - m_left * (x[1]-z))
    if (x[n] <= z) return(y[n] + m_right * (z-x[n]))
    
    x_upp = x[z<=x][1]
    x_low = tail(x[x<z], n=1)
    y_upp = y[z<=x][1]
    y_low = tail(y[x<z], n=1)
    
    y_low + (y_upp-y_low)/(x_upp-x_low)*(z-x_low)
  }
}

compute_rmst <- function(time, status, tau, alpha=0.05) {
  ft <- survival::survfit(Surv(time, status)~1) 
  
  t <- ft$time 
  d <- ft$n.event
  n <- ft$n.risk
  
  # Calculate Kaplan-Meier estimator for survîval function
  surv <- cumprod(1-d/n)
  
  # Assure that the first time entry is 0
  if (t[1] != 0) {
    t <- c(0, t)
    surv <- c(1, surv)
  }
  
  # Integrate survival rate to get RMST function
  t_diff <- diff(t)
  rmst_arr <- c(0, cumsum(t_diff*surv[1:(length(surv)-1)]))
  
  # Transform into piecewise linear function to obtain rmst for tau values that
  # are not contained in the time array t.
  rmst_fun <- piecewise_linearise(t, rmst_arr, m_right = surv[length(surv)])
  rmst <- rmst_fun(tau)
  
  # Use Greenwood's formula to approximate rmst's variance
  # (E.g. see [Lawrence J, Qiu J, Bai S, Hung HMJ. Difference in restricted mean
  # survival time: small sample distribution and asymptotic relative efficiency.
  # Stat Biopharm Res. 2019;11:61-66])
  var_fac <- ifelse(n-d==0, 0, d/(n*(n-d)) )
  rmst_var <- sum((rmst-rmst_arr[t<=tau])^2 * var_fac[t<=tau])
  
  
  list(
    rmst = rmst,
    sd = sqrt(rmst_var),
    CIlow = rmst - qnorm(1-alpha/2) * sqrt(rmst_var),
    CIupp = rmst + qnorm(1-alpha/2) * sqrt(rmst_var)
  )
}

test_that("Test rmst2 without covariates given.", {
  # Cmpute values for rmst and its sd via survRM2::rmst2.
  res <- rmst2(time, status, arm, tau)
  
  
  # arm==0 computations and data extractions via own rmst computations
  time0 <- time[arm==0] 
  status0 <- status[arm==0]
  
  arm0_data <- compute_rmst(time0, status0, tau)
  
  rmst0_cmp <- arm0_data$rmst
  rmst0.sd_cmp <- arm0_data$sd
  rmst0.CIlow_cmp <- arm0_data$CIlow
  rmst0.CIupp_cmp <- arm0_data$CIupp
  
  rmtl0_cmp <- tau - rmst0_cmp
  rmtl0.sd_cmp <- rmst0.sd_cmp
  rmtl0.CIlow_cmp <- rmtl0_cmp - (rmst0_cmp - rmst0.CIlow_cmp)
  rmtl0.CIupp_cmp <- rmtl0_cmp + (rmst0_cmp - rmst0.CIlow_cmp)
  
  # extract data from results for arm==0
  rmst0_res <- res$RMST.arm0$rmst
  rmtl0_res <- res$RMST.arm0$rmtl
  
  # compare if equal 
  expect_equal(rmst0_cmp, rmst0_res[[1]], tolerance = 1e-6) 
  expect_equal(rmst0.sd_cmp, rmst0_res[[2]], tolerance = 1e-6) 
  expect_equal(rmst0.CIlow_cmp, rmst0_res[[3]], tolerance = 1e-6) 
  expect_equal(rmst0.CIupp_cmp, rmst0_res[[4]], tolerance = 1e-6) 
  
  expect_equal(rmtl0_cmp, rmtl0_res[[1]], tolerance = 1e-6) 
  expect_equal(rmtl0.sd_cmp, rmtl0_res[[2]], tolerance = 1e-6) 
  expect_equal(rmtl0.CIlow_cmp, rmtl0_res[[3]], tolerance = 1e-6) 
  expect_equal(rmtl0.CIupp_cmp, rmtl0_res[[4]], tolerance = 1e-6) 
  
  
  # arm==1 computations and data extractions via own rmst computations
  time1 <- time[arm==1] 
  status1 <- status[arm==1]
  
  arm1_data <- compute_rmst(time1, status1, tau)
  
  rmst1_cmp <- arm1_data$rmst
  rmst1.sd_cmp <- arm1_data$sd
  rmst1.CIlow_cmp <- arm1_data$CIlow
  rmst1.CIupp_cmp <- arm1_data$CIupp
  
  rmtl1_cmp <- tau - rmst1_cmp
  rmtl1.sd_cmp <- rmst1.sd_cmp
  rmtl1.CIlow_cmp <- rmtl1_cmp - (rmst1_cmp - rmst1.CIlow_cmp)
  rmtl1.CIupp_cmp <- rmtl1_cmp + (rmst1_cmp - rmst1.CIlow_cmp)
  
  # extract data from results for arm==1
  rmst1_res <- res$RMST.arm1$rmst
  rmtl1_res <- res$RMST.arm1$rmtl
  
  # compare if equal 
  expect_equal(rmst1_cmp, rmst1_res[[1]], tolerance = 1e-6) 
  expect_equal(rmst1.sd_cmp, rmst1_res[[2]], tolerance = 1e-6) 
  expect_equal(rmst1.CIlow_cmp, rmst1_res[[3]], tolerance = 1e-6) 
  expect_equal(rmst1.CIupp_cmp, rmst1_res[[4]], tolerance = 1e-6) 
  
  expect_equal(rmtl1_cmp, rmtl1_res[[1]], tolerance = 1e-6) 
  expect_equal(rmtl1.sd_cmp, rmtl1_res[[2]], tolerance = 1e-6) 
  expect_equal(rmtl1.CIlow_cmp, rmtl1_res[[3]], tolerance = 1e-6) 
  expect_equal(rmtl1.CIupp_cmp, rmtl1_res[[4]], tolerance = 1e-6) 
  
  
  alpha <- 0.05 
  
  # Check RMST difference
  rmst_diff_cmp <- rmst1_cmp - rmst0_cmp
  rmst_diff.sd_cmp <- sqrt(rmst0.sd_cmp^2 + rmst1.sd_cmp^2)
  rmst_diff.CIlow_cmp <- rmst_diff_cmp - qnorm(1-alpha/2) * rmst_diff.sd_cmp
  rmst_diff.CIupp_cmp <- rmst_diff_cmp + qnorm(1-alpha/2) * rmst_diff.sd_cmp
  rmst_diff.pval_cmp <- pnorm(-abs(rmst_diff_cmp)/rmst_diff.sd_cmp)*2
  
  expect_equal(rmst_diff_cmp, res$unadjusted.result[1])
  expect_equal(rmst_diff.CIlow_cmp, res$unadjusted.result[4])
  expect_equal(rmst_diff.CIupp_cmp, res$unadjusted.result[7])
  expect_equal(rmst_diff.pval_cmp, res$unadjusted.result[10])
  
  
  # Check RMST ratio
  rmst_logratio_cmp <- log(rmst1_cmp) - log(rmst0_cmp)
  # Use multivariate delta-method to compute its sd:
  # Var(log(X)-log(Y)) = Var(X)/x^2 + Var(Y)/y^2 for (X,Y) close to (x,y) and X,Y independent.
  rmst_logratio.sd_cmp <- sqrt(rmst0.sd_cmp^2/rmst0_cmp^2 + rmst1.sd_cmp^2/rmst1_cmp^2)
  rmst_logratio.CIlow_cmp <- rmst_logratio_cmp - qnorm(1-alpha/2) * rmst_logratio.sd_cmp
  rmst_logratio.CIupp_cmp <- rmst_logratio_cmp + qnorm(1-alpha/2) * rmst_logratio.sd_cmp
  rmst_logratio.pval_cmp <- pnorm(-abs(rmst_logratio_cmp)/rmst_logratio.sd_cmp)*2
  
  expect_equal(exp(rmst_logratio_cmp), res$unadjusted.result[2])
  expect_equal(exp(rmst_logratio.CIlow_cmp), res$unadjusted.result[5])
  expect_equal(exp(rmst_logratio.CIupp_cmp), res$unadjusted.result[8])
  expect_equal(rmst_logratio.pval_cmp, res$unadjusted.result[11])
  
  
  # Check RMTL ratio
  rmtl_logratio_cmp <- log(rmtl1_cmp) - log(rmtl0_cmp)
  rmtl_logratio.sd_cmp <- sqrt(rmtl0.sd_cmp^2/rmtl0_cmp^2 + rmtl1.sd_cmp^2/rmtl1_cmp^2)
  rmtl_logratio.CIlow_cmp <- rmtl_logratio_cmp - qnorm(1-alpha/2) * rmtl_logratio.sd_cmp
  rmtl_logratio.CIupp_cmp <- rmtl_logratio_cmp + qnorm(1-alpha/2) * rmtl_logratio.sd_cmp
  rmtl_logratio.pval_cmp <- pnorm(-abs(rmtl_logratio_cmp)/rmtl_logratio.sd_cmp)*2
  
  expect_equal(exp(rmtl_logratio_cmp), res$unadjusted.result[3])
  expect_equal(exp(rmtl_logratio.CIlow_cmp), res$unadjusted.result[6])
  expect_equal(exp(rmtl_logratio.CIupp_cmp), res$unadjusted.result[9])
  expect_equal(rmtl_logratio.pval_cmp, res$unadjusted.result[12])
  
})
