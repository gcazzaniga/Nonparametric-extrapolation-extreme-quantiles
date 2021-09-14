### Load required functions ###################################################

source("2_functions.R")

### Main #######################################################################

## Initialize dataframes and list for final memorization of results
quantiles <- data.frame(Sample_length = integer(),
                        Return_period = integer(),
                        N_sim = integer(),
                        HTS_quantile = double(),
                        KDE_quantile = double(),
                        SC0_quantile = double(),
                        SC7_quantile = double(),
                        Theoretical_quantile = double(),
                        HTS_error = double(),
                        KDE_error = double(),
                        SC0_error = double(),
                        SC7_error = double()) 
parameters <- data.frame(Sample_length = integer(),
                         N_sim = integer(),
                         HTS_w = double(),
                         HTS_h = double(),
                         SC0_EVI = double(),
                         SC0_b1 = double(), 
                         SC0_b2 = double(),  
                         SC0_R2 = double(),
                         SC7_EVI = double(),
                         SC7_b1 = double(), 
                         SC7_b2 = double(),  
                         SC7_R2 = double()) 
invalid_data <-  data.frame(Sample_length = integer(),
                            N_sim = integer(),
                            R2_over_0 = logical(),
                            R2_over_0.7 = logical()) 
synthetic_dataset <- list()


for (n in sample_len){
  message(sprintf("Starting analysis for sample length %d...", n))
  ## Initialize dataframe and variables for simulations with sample length n ####
  
  dataset_all <- data.frame(matrix(nrow = n, ncol = n_sim_max))
  colnames(dataset_all) <- as.character(seq(1, n_sim_max, 1))
  
  n_NA <- 0
  n_sim <- 1
  
  ## Compute results for each simulation ####
  
  while (n_sim < (n_sim_max + 1)) {
    
    SC0 <-  data.frame(matrix(nrow = 1, ncol = 4))
    colnames(SC0) <- c("EVI", "b1", "b2", "R2")
    SC7 <-  data.frame(matrix(nrow = 1, ncol = 4))
    colnames(SC7) <- c("EVI", "b1", "b2", "R2")
    
    # Generate synthetic dataset ####
    
    dataset <- sample_gen(distribution, n)
    
    # Compute order statistic ####
    
    Y <- sort(dataset, decreasing = TRUE)
    
    # Probability levels ####
    
    u <- 1 - (1 / return_period)
    u <- u[which(u >= (n / (n+1)))]
    
    # Theoretical quantiles ####
    
    Y_t <- quantile_gen(distribution, u)
    
    ## Quantiles with Hutson (2002) ####
    
    if (distribution$name == "Uniform"){
      
      # data transformation to have upper unbounded data
      lower <- distribution$par1
      upper <- distribution$par2
      dataset_log <- log((dataset - lower) / (upper - dataset)) 
      Y_log <- sort(dataset_log, decreasing = TRUE)
      Y_h_log <- Y_log[1] - (Y_log[1] - Y_log[2]) * log((n + 1) * (1 - u))
      Y_h <- (upper * exp(Y_h_log) + lower) / (1 + exp(Y_h_log))
    } else {
      Y_h <- Y[1] - (Y[1] - Y[2]) * log((n + 1) * (1 - u))
    }
    
    ## Quantiles with kernel density smoothing ####
    
    # Parameters optimization
    h0 <- bw.nrd0(dataset)
    w0 <- 0.5
    log <- capture.output(res <- sceua(lcv_obj, pars = c(h0, w0), 
                                       lower = c(0, 0), upper = c(5, 1),
                                       data = dataset))
    h <- res$par[1]
    w <- res$par[2]
    
    # For Gamma and uniform distribution the PDF is transformed and normalized 
    # in order to have a positive support and a unitary area
    if (distribution$name == "Gamma"){
      
      xl <- 0.1
      xu <- (2 * max(Y_t))
      x <- seq(xl, xu, by = 0.1)
      area <- integrate(KDE,lower = 0, upper = Inf, subdivisions = 200L, data = dataset, h = h, w = w)
      f_x <- KDE(x, data = dataset, h = h, w = w) / area$value 
      
    } else if (distribution$name ==  "Cauchy"){
      
      xl <- min(dataset)
      xu <- max(1000, 2 * max(dataset))
      x <- seq(xl, xu, by = 0.1)
      f_x <- KDE(x, data = dataset, h = h, w = w)
      
    } else if (distribution$name == "Uniform") {
      
      xl <- distribution$par1
      xu <- distribution$par2
      x <- seq(xl, xu, 0.001)
      area <- integrate(KDE, lower = xl, upper = xu, data = dataset, h = h, w = w)
      f_x <- KDE(x, data = dataset, h = h, w = w) / area$value
    }
    
    F_x <- pdf2cdf(f_x, x = x, normalize = TRUE, expand = FALSE)
    F_x_fun <- approxfun(F_x$y, F_x$x)
    Y_k <- F_x_fun(u) 
    
    ## Quantiles with  Scholz (1995) ####
    
    gamma <- 0.50
    prob <- 1 - (1:n - 1 / 3) / (n + 1 / 3)
    data_median <- median(dataset)
    Y_tilde <- Y - data_median
    
    # Define the k values that are tested
    K1 <- floor(max(6, 1.3 * sqrt(n)))
    K2 <- 2 * floor(log10(n) * sqrt(n))
    if ((K2-K1+1) %% 2 == 0) {     
      n_k <- K2 - K1   
    } else {
      n_k <- K2 - K1 + 1
    }
    k_trial <- data.frame(matrix(ncol = 5, nrow = n_k ))
    colnames(k_trial) <- c("k", "R2", "b1", "b2", "c")
    
    # Parameters' estimation for each k 
    for (index in 1:n_k) {
      
      k_trial$k[index] <- K1 + index - 1
      k <- k_trial$k[index]
      
      # EVI estimation
      Y_tilde_k <- Y_tilde[k]
      Y_log_ratio <- log(Y_tilde[1:(k-1)] / Y_tilde_k)
      M1_k <- (1 / (k-1)) * sum(Y_log_ratio)
      M2_k <- (1 / (k-1)) * sum((Y_log_ratio) ^ 2)
      c_hat_k <- M1_k + 1 - 0.5 * (1 - ((M1_k) ^ 2 / M2_k)) ^ (-1)
      k_trial$c[index] <- c_hat_k
      
      # Covariance matrix of order statistic  
      S <- matrix(NA, nrow = k, ncol = k)
      
      for (i in 1:k) {
        for (j in 1:i) {
          S[i, j] <- i ^ (-c_hat_k - 1) * j ^ (-c_hat_k)
          S[j, i] <- S[i, j]
        }
      }
      
      # Least squares solution 
      fc_k <- (((-n * log(prob)) ^ (-c_hat_k)) - 1) / c_hat_k
      X <- cbind(rep(1, k), fc_k[1:k])
      
      if (check(S)){   # Least squares are solved only if S is invertible
        b <- solve(t(X) %*% solve(S) %*% X) %*% (t(X) %*% solve(S) %*% (Y[1:k]))
        k_trial$b1[index] <- b[1]
        k_trial$b2[index] <- b[2]
        est <- (b[1] + fc_k[1:k] * b[2])
        obs <- Y[1:k]
        TSS <- sum((obs - mean(obs)) ^ 2)
        RSS <- sum((obs-est) ^ 2)
        k_trial$R2[index] <- 1-(RSS/TSS)
      } else {
        k_trial$R2[index]<-NA
      }
    }
    
    # Compute optimum EVI and linear regression parameters for SC0
    threshold <- 0
    k_trial<- k_trial[which(k_trial$R2 > threshold), ]
    
    # Memorize if no optimum EVI is found for SC0
    if (sum(which(k_trial$R2 > threshold)) == 0) {
      SC0_check_R2 <- FALSE   
      SC7_check_R2 <- FALSE # If no optimum EVI is found for SC0, no optimum EVI will be found also for SC7
    } else {
      SC0_check_R2 <- TRUE
      c_mean <- mean(k_trial$c, na.rm = T)
      SC0$EVI <- c_mean
      b1_mean <- mean(k_trial$b1, na.rm = T)
      b2_mean <- mean(k_trial$b2, na.rm = T)
      SC0$b1 <- b1_mean
      SC0$b2 <- b2_mean
      fc_k <- (((- n * log(prob)) ^ (- c_mean)) - 1) / c_mean
      est <- (b1_mean + fc_k[1:k] * b2_mean)
      obs <- Y[1:k]
      TSS <- sum((obs - mean(obs)) ^ 2)
      RSS <- sum((obs - est) ^ 2)
      SC0$R2 <- 1 - (RSS / TSS)
    }
    
    if (isTRUE(SC0_check_R2)){ #
      # Compute optimum EVI and linear regression parameters for SC7
      threshold <- 0.7
      k_trial<- k_trial[which(k_trial$R2 > threshold), ]
      
      # Memorize if no optimum EVI is found for SC7
      if (sum(which(k_trial$R2 > threshold)) == 0) {
        SC7_check_R2 <- FALSE
      } else {
        SC7_check_R2 <- TRUE
        c_mean <- mean(k_trial$c, na.rm = T)
        SC7$EVI <- c_mean
        b1_mean <- mean(k_trial$b1, na.rm = T)
        b2_mean <- mean(k_trial$b2, na.rm = T)
        SC7$b1 <- b1_mean
        SC7$b2 <- b2_mean
        fc_k <- (((- n * log(prob)) ^ (- c_mean)) - 1) / c_mean
        est <- (b1_mean + fc_k[1:k] * b2_mean)
        obs <- Y[1:k]
        TSS <- sum((obs - mean(obs)) ^ 2)
        RSS <- sum((obs - est) ^ 2)
        SC7$R2 <- 1 - (RSS / TSS)
        
        # Calculate extrapolated quantiles for SC0 and SC7 only if optimum EVIs are found
        Y_s_0 <- SC0$b1 + SC0$b2 * ((((- n * log(u)) ^ (- SC0$EVI)) - 1) / SC0$EVI)
        Y_s_7 <- SC7$b1 + SC7$b2 * ((((- n * log(u)) ^ (- SC7$EVI)) - 1) / SC7$EVI)
      }
    }
    
    ## Memorize results ####
    
    if (isTRUE(SC0_check_R2) && isTRUE(SC7_check_R2)){
      
      # Save all extrapolated quantiles for the n-th synthetic dataset
      sim_quantiles <- data.frame(Sample_length = rep(n, times = length(return_period)), 
                                  Return_period = return_period, 
                                  N_sim = rep(n_sim, times = length(return_period)), 
                                  HTS_quantile = c(rep(NA, times = length(return_period) - length(u)), Y_h),  
                                  KDE_quantile = c(rep(NA, times = length(return_period) - length(u)), Y_k), 
                                  SC0_quantile = c(rep(NA, times = length(return_period) - length(u)), Y_s_0), 
                                  SC7_quantile = c(rep(NA, times = length(return_period) - length(u)), Y_s_7), 
                                  Theoretical_quantile = c(rep(NA, times = length(return_period) - length(u)), Y_t))
      
      # Calculate relative error (from Scholz & Tjoelker, 1995) on the extrapolated quantiles for the n-th synthetic dataset
      sim_err <- 100 * (sim_quantiles[ , 4:7] - sim_quantiles$Theoretical_quantile) / (sim_quantiles$Theoretical_quantile - median(dataset))
      names(sim_err) <- c("HTS_error", "KDE_error", "SC0_error", "SC7_error")  
      
      # Save optimized parameters (from HTS and KDE) for the n-th synthetic dataset
      sim_parameters <- data.frame(Sample_length = n, 
                                   N_sim = n_sim, 
                                   HTS_w = w,  
                                   HTS_h = h,
                                   SC0_EVI = SC0$EVI, 
                                   SC0_b1 = SC0$b1, 
                                   SC0_b2 = SC0$b2,
                                   SC0_R2 = SC0$R2,
                                   SC7_EVI = SC7$EVI, 
                                   SC7_b1 = SC7$b1, 
                                   SC7_b2 = SC7$b2,
                                   SC7_R2 = SC7$R2)
      
      # Save results in dataframes wich include all simulations
      quantiles <- rbind(quantiles,cbind(sim_quantiles, sim_err))
      parameters <- rbind(parameters, sim_parameters)
      dataset_all[[as.character(n_sim)]] <- dataset
      n_sim <- n_sim + 1
    } else {
      # Save simulations which provide invalid data (when no optimum EVI with SC0 or SC7 is found)
      n_NA <- n_NA + 1
      sim_invalid_data <- data.frame(Sample_length = n, 
                                     N_sim = n_NA, 
                                     R2_over_0 = SC0_check_R2, 
                                     R2_over_0.7 = SC7_check_R2)
      invalid_data <- rbind(invalid_data, sim_invalid_data)
    }
  }
  synthetic_dataset[[paste("Sample_length_",as.character(n))]] <- dataset_all
}

## Calculate median error, 25-th percentile (Q1) and 75-th percentile (Q3) with respect to sample length, return period, and method ####
quantiles_melt <- melt(quantiles, id.vars = c("Sample_length", "Return_period"), measure.vars = c("HTS_error", "KDE_error", "SC0_error", "SC7_error"))
names(quantiles_melt)[3:4] <- c("Method", "Error") 
Median_error <- aggregate(quantiles_melt$Error, by = list(quantiles_melt$Sample_length, quantiles_melt$Return_period, quantiles_melt$Method), function(y) median(y, na.rm = T))
Q1_error <- aggregate(quantiles_melt$Error, by = list(quantiles_melt$Sample_length, quantiles_melt$Return_period, quantiles_melt$Method), function(y) quantile(y, 0.25, na.rm = T))
Q3_error <- aggregate(quantiles_melt$Error, by = list(quantiles_melt$Sample_length, quantiles_melt$Return_period, quantiles_melt$Method), function(y) quantile(y, 0.75, na.rm = T))
quartiles_error <- data.frame(Sample_length = Median_error$Group.1,
                              Return_period = Median_error$Group.2,
                              Method = Median_error$Group.3,
                              Median_error = Median_error$x,
                              Q1_error = Q1_error$x,
                              Q3_error = Q3_error$x)

## Save results in .rds files ####
saveRDS(quantiles, file = paste0("EstimatedQuantiles_and_Errors_", distribution$name, ".rds"))
saveRDS(parameters, file = paste0("OptimumParameters_", distribution$name, ".rds"))
saveRDS(synthetic_dataset, file = paste0("SyntheticDataSet_", distribution$name, ".rds"))
saveRDS(invalid_data, file = paste0("InvalidData_", distribution$name, ".rds"))
saveRDS(quartiles_error, file = paste0("QuartilesErrors_", distribution$name, ".rds"))