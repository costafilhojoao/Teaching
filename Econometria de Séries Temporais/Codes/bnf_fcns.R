####################################################################################################
# Replication Files for Kamber, Morley & Wong:
#
# "Intuitive and Reliable Estimates of the Output Gap from a Beveridge-Nelson Filter"
# The Review of Economics and Statistics
#
# R conversion of some of the MATLAB codes originally written by Ben Wong
#
####################################################################################################
# MATLAB codes converted to R by Luke Hartigan, 2017
# Additional R codes and wrapper class 'bnf' written by Luke Hartigan, 2017
# Updated by James Morley, 2022
# Additional code changes by Cristian, 2022/2023/2024
####################################################################################################

# Helper functions used by the main functions below

# R emulation of MATLAB eye(N)
eye <- function(n)
{
    return (diag(x = 1, nrow = n))
}

# R emulation of MATLAB ones(M, N)
ones <- function(m, n)
{
    return (matrix(data = 1, nrow = m, ncol = n))
}

# R emulation of MATLAB zeros(M, N)
zeros <- function(m, n)
{
    return (matrix(data = 0, nrow = m, ncol = n))
}

# R emulation of MATLAB repmat(A, M, N)
repmat <- function(mat, m, n)
{
    return (kronecker(matrix(data = 1, nrow = m, ncol = n), mat))
}

# Compute the square of the variable 'x'
square <- function(x)
{
    return (x * x)
}

# Returns a suitably transformed, lagged and/or differenced data.
# OPTIONS:
# take_log = take the natural logarithm? Default is no (FALSE)
#
# dcode = option to specify how y is differenced:
#   dcode == "nd" [No difference, i.e., Level] -- Default --
#   dcode == "d1" [1st Difference, i.e., (1 - B)y]
#   dcode == "d4" [4th Difference, i.e., (1 - B^4)y, use with quarterly data]
#   dcode == "d12" [12th Difference, i.e., (1 - B^12)y, use with monthly data]
#
# pcode = option to specify if percentages are computed:
#   pcode == "np" [no change] -- Default --
#   pcode == "p1" [multiply by 100]
#   pcode == "p4" [multiply by 400, annualised quarterly rate]
#   pcode == "p12" [multiply by 1200, annualised monthly rate]
#
# NB: Inspired by the GAUSS procs originally written by Mark Watson:
# http://www.princeton.edu/~mwatson/wp.html
#
transform_series <- function(y, take_log = FALSE, 
    dcode = c("nd", "d1", "d4", "d12"), pcode = c("np", "p1", "p4", "p12"))
{
    # Save 'ts' attributes if 'y' is a 'ts' object
    if (is.ts(y)) {
        y_ts <- TRUE
        endy <- end(y)
        freqy <- frequency(y)
    } else {
        y_ts <- FALSE
    }  

    # Preliminary stuff  
    ty <- as.matrix(y)

    # Log transformation    
    if (take_log) {
        if (any(ty < 0, na.rm = TRUE)) {
            ty[which(ty < 0)] <- NaN
        }
        ty <- log(as.matrix(ty))
    }

    # Difference transform if requested
    dcode <- match.arg(dcode)

    switch (EXPR = dcode,
        "nd" = { ty }, # do nothing
        "d1" = { ty <- diff(x = ty, lag = 1); },
        "d4" = { ty <- diff(x = ty, lag = 4) },
        "d12" = { ty <- diff(x = ty, lag = 12) }
    )

    # Convert to percentage if requested
    pcode <- match.arg(pcode)

    switch (EXPR = pcode,
        "np" = { ty }, # do nothing
        "p1" = { ty <- ty * 100.0; },
        "p4" = { ty <- ty * 400.0 },
        "p12" = { ty <- ty * 1200.0 }
    )

    # Add 'ts' attribute to 'ty' if is a 'ts' object
    if (y_ts) {
        ty <- ts(ty, end = endy, frequency = freqy)
    }

    return (ty)
}

# Demean relative to the break point(s) provided in 'breaks' 
# The values in 'breaks' are assumed to be in ascending order
# NB: Works even if no breaks given, then full sample mean is used instead
piecewise_demean <- function(y, breaks = c())
{
    # Preliminary stuff
    y <- as.matrix(y)
    tobs <- nrow(y)

    # Allocate storage for the demeaned series
    demeaned_y <- matrix(data = NA_real_, nrow = tobs, ncol = 1)
    br <- c(0, breaks, tobs)
    nb <- length(breaks)

    # Number of break points 'nb' should not be greater than T-1
    if (nb > (tobs - 1)) {
        stop(sprintf("Number of break points \'%d\' not feasible with T - 1 = %d\n",
        nb, tobs - 1))
    }

    # Last break point should not be greater than number of observations 'T'
    if (max(breaks) > tobs) {
        stop(sprintf("Largest break point \'%d\' greater than time series length \'%d\'\n",
        breaks[nb], tobs))
    }
    
    for (i in  1:(nb + 1)) {
        j <- br[i] + 1
        k <- br[i + 1]
        demeaned_y[j:k] <- scale(x = y[j:k], center = TRUE, scale = FALSE)
    }
    
    return (demeaned_y)
}


# Demean relative to a backward rolling window 'wind'
rolling_demean <- function(y, y_cycle, wind)
{
    # Preliminary stuff
    y <- as.matrix(y)
    y_cycle <- as.matrix(y_cycle)
    tobs <- nrow(y)
    windp1 <- wind + 1
    windm1 <- wind - 1

    if (wind > tobs) {
        stop(sprintf("Window length \'%d\' greater than time series length \'%d\'\n",
            wind, tobs))
    }

    # Allocate storage for the demeaned series
    demeaned_y <- matrix(data = NA_real_, nrow = tobs, ncol = 1)

    demeaned_y[1:wind] <- t(t(y[1:wind])) - mean(y[1:wind]) + mean(y_cycle[1:wind])

    for (i in windp1:tobs) {
        demeaned_y[i] <- y[i] - mean(y[(i - windm1):i]) + mean(y_cycle[(i - windm1):i])
    }

    return (demeaned_y)
}


# Estimate an OLS reduced form VAR
# INPUTS:
# y      Data(T x N)
# p      lags
# nc     put TRUE if you don't want a constant (default is FALSE)
#
# OUTPUTS:
# A      VAR Coefficients
# SIGMA  Reduced form variance covariance matrix
# U      Time Series of Residuals
# invXX  Inverse of the symmetric square matrix X'X
# X      Design matrix
olsvar <- function(y, p, nc = FALSE)
{
    # Preliminary stuff
    y <- as.matrix(y)
    tobs <- nrow(y) - p                         # length of time series
    Y <- y[(p + 1):nrow(y), , drop = FALSE]     # cut away first p lags
    
    if (!nc) {
        X <- rep(1, tobs)
    } else {
        X <- c();
    }

    for (i in 1:p) {
        Z <- y[((p + 1) - i):(nrow(y) - i), , drop = FALSE]
        X <- cbind(X, Z)
    }

    # Do the QR decomposition
    qr_xmat <- qr(X)

    # Get the OLS coefficients
    A <- qr.coef(qr = qr_xmat, y = Y)

    # Get the OLS residuals
    U <- qr.resid(qr = qr_xmat, y = Y)
    
    # Get the OLS VCov matrix
    # df is total sample-lags-number of regressors
    SIGMA <- sum(square(U)) / (tobs - nrow(A))

    invXX <- qr.solve(crossprod(x = X, y = NULL),tol = 1e-30)

    # Collect results as a list object
    olsvar <- list()
    olsvar$A <- A
    olsvar$SIGMA <- SIGMA
    olsvar$U <- U
    olsvar$invXX <- invXX
    olsvar$X <- X

    return (olsvar)
}

# Estimate an AR or VAR in Dickey-Fuller form without a constant
#
# INPUTS:
# y      Data(T x N)
# p      lags
# nc     put TRUE if you don't want a constant (default is FALSE)
#
# OUTPUTS:
# A      VAR Coefficients
# SIGMA  Reduced form variance covariance matrix
# U      Time Series of Residuals
# invXX  Inverse of the symmetric square matrix X'X
# X      Design matrix
olsvar_df <- function(y, p, nc = FALSE)
{
    # Preliminary stuff
    y <- as.matrix(y)
    tobs <- nrow(y) - p                         # length of time series
    Y <- y[(p + 1):nrow(y), , drop = FALSE]     # cut away first p lags
    
    if (!nc) {
        X <- rep(1, tobs);
    } else {
        X <- c();
    }

    # Put in first lag on LHS
    Z <- y[p:(nrow(y) - 1), , drop = FALSE]
    X <- cbind(X, Z)
    Y <- Y[2:nrow(Y), , drop = FALSE]
    X <- X[2:nrow(X), , drop = FALSE]

    # Now do second differences
    dy <- diff(y)
    for (i in 1:(p - 1)) {
        Z <- dy[((p + 1) - i):(nrow(dy) - i), , drop = FALSE]
        X <- cbind(X, Z)
    }

    # Do the QR decomposition
    qr_xmat <- qr(X)

    # Get the OLS coefficients
    A <- qr.coef(qr = qr_xmat, y = Y)

    # Get the OLS residuals
    U <- qr.resid(qr = qr_xmat, y = Y)
    
    # Get the OLS VCov matrix
    # df is total sample-lags-number of regressors
    SIGMA <- sum(square(U)) / (tobs - nrow(A));

    invXX <- qr.solve(crossprod(x = X, y = NULL),tol = 1e-30)

    # Collect results as a list object
    olsvar_df <- list()
    olsvar_df$A <- A
    olsvar_df$SIGMA <- SIGMA
    olsvar_df$U <- U
    olsvar_df$invXX <- invXX
    olsvar_df$X <- X

    return (olsvar_df)

}

# Estimate an OLS reduced form VAR and adjust estimated variance of errors for outliers
# INPUTS:
# y      Data(T x N)
# p      lags
# temp_dummies	dates of outliers
# nc     put TRUE if you don't want a constant (default is FALSE)
#
# OUTPUTS:
# A      VAR Coefficients
# SIGMA  Reduced form variance covariance matrix
# U      Time Series of Residuals
# invXX  Inverse of the symmetric square matrix X'X
# X      Design matrix
olsvar_outliers <- function(y, p, temp_dummies, nc = FALSE)
{
    # Preliminary stuff
    y <- as.matrix(y)
    tobs <- nrow(y) - p                         # length of time series
    Y <- y[(p + 1):nrow(y), , drop = FALSE]     # cut away first p lags
    
    if (!nc) {
        X <- rep(1, tobs)
    } else {
        X <- c();
    }

    for (i in 1:p) {
        Z <- y[((p + 1) - i):(nrow(y) - i), , drop = FALSE]
        X <- cbind(X, Z)
    }

    # Do the QR decomposition
    qr_xmat <- qr(X)

    # Get the OLS coefficients
    A <- qr.coef(qr = qr_xmat, y = Y)

    # Get the OLS residuals
    U <- qr.resid(qr = qr_xmat, y = Y)
    
    temp_dummies <- temp_dummies - p;
    U_adj <- t(t(U[-temp_dummies]))
    
    # Get the OLS VCov matrix
    # df is total sample-lags-number of regressors
    SIGMA <- sum(square(U_adj)) / (nrow(U_adj) - nrow(A));

    invXX <- qr.solve(crossprod(x = X, y = NULL),tol = 1e-30)

    # Collect results as a list object
    olsvar_outliers <- list()
    olsvar_outliers$A <- A
    olsvar_outliers$SIGMA <- SIGMA
    olsvar_outliers$U <- U
    olsvar_outliers$invXX <- invXX
    olsvar_outliers$X <- X

    return (olsvar_outliers)
}


# Backcast observations based on reversible process
#
# Note: This function depends on functions implemented above
#
# INPUTS:
# y           Time series to implement BN Filter in first differences
# p           AR lags
# delta       Signal to noise ratio
#
# OUTPUTS:
# BN_cycle    The estimated cycle
# BN_cycle_se The estimated cycle standard error
# aux_out     Auxiliary output (AR coefficients & residuals )
backcast <- function(y, p, delta){
    # Do a few preliminary things
    y_flip <- rev(y)
    y_flip <- as.matrix(y_flip)
	y <- rev(y_flip)
	y <- as.matrix(y)
	nr_y = nrow(y) 

    # Compute rho from delta
    rho <- -(1.0 - sqrt(delta)) / sqrt(delta)

    # Re-parametrise model and subtract second differences multiplied by 
    # sum of AR coefficients to fix delta. 
    # Adjust the LHS variable to take note of this

    # Backcast data 
    augmented_y <- rbind(rbind(zeros(p+2,1), y),zeros(p+2,1)) 
    augmented_y_iter <- cbind(ones(nrow(augmented_y),1),augmented_y)
    
    iter <- 2
    end_aug = p+2
    while (iter < 100 && abs(mean(augmented_y_iter[1:end_aug,iter]-augmented_y_iter[1:end_aug,iter-1]))>0.00001*sd(y)){
    	
    	augmented_y_flip <- rev(augmented_y)
    	augmented_y_flip <- as.matrix(augmented_y_flip)

    # Get data matrix from untransformed model. Estimate without a constant
    end_row = p+2+nr_y
    tmp_olsvar_df <- olsvar_df(y = t(t(augmented_y[2:end_row])), p = p, nc = TRUE)
    tmp_olsvar_df_flip <- olsvar_df(y = t(t(augmented_y_flip[2:end_row])), p = p, nc = TRUE)
    
    X <- tmp_olsvar_df$X[ , , drop = FALSE]
    nr_X = nrow(X)
    nc_X = ncol(X)
    X_flip <- tmp_olsvar_df_flip$X[ , , drop = FALSE]
    

    # Get regressors, variance and lags of untransformed AR(p) to perform BN later.
    # Estimate AR by OLS with (for sig_e) and without a constant
    end_row = p+2
    tmp_olsvar <- olsvar(y = rbind(t(t(augmented_y[3:end_row])),y), p = p, nc = TRUE)
    sig2_ols <- tmp_olsvar$SIGMA
    tmp_olsvar_flip <- olsvar(y = rbind(t(t(augmented_y_flip[3:end_row])),y_flip), p = p, nc = TRUE)
    
    X_untransformed <- tmp_olsvar$X[ , , drop = FALSE]
    nr_X_untransformed = nrow(X_untransformed)
    nc_X_untransformed = ncol(X_untransformed)
    X_flip_untransformed <- tmp_olsvar_flip$X[ , , drop = FALSE]
    
    start_row = p+3
    end_row = p+2+nr_y   
    reparameterised_y <- t(t(augmented_y[start_row:end_row, , drop = FALSE])) - 
        rho * X[,1]
    reparameterised_y_flip <- t(t(augmented_y_flip[start_row:end_row, , drop = FALSE])) - 
        rho * X_flip[,1]
 
    #Calculate other preliminaries before estimation

    # Create Data for transformed model to estimate
    Y <- reparameterised_y
    X <- X[1:nr_X,2:nc_X]
    Y_flip <- reparameterised_y_flip
    X_flip <- X_flip[1:nr_X,2:nc_X]

    # Set other parameters for conjugate priors
    A_prior <- zeros(p - 1, 1)     # constant and (p-1) phi_star coefficients
    # Uninformative prior on constant centred on zero
    V_prior <- diag(drop(repmat(0.5, 1, p - 1) / (square(seq(from = 1, to = p - 1, by = 1)))))

    # Calculate Posteriors
    V_post <- qr.solve(qr.solve(V_prior,tol = 1e-30) + (1.0 / sig2_ols) * crossprod(x = X, y = NULL),tol = 1e-30)
    A_post <- V_post %*% (qr.solve(V_post,tol = 1e-30) %*% A_prior + t((1.0 / sig2_ols) * X) %*% Y)
    
    # Add rho back into the posterior of estimated parameters
    A_post <- rbind(rho, A_post)

    # Map the parameters in DF-form back to AR-form
    AR_params <- zeros(p, 1)

    # Last lag
    AR_params[p] <- -A_post[p]
    rseq <- seq(from = p - 1, to = 2, by = -1)

    for (i in rseq) {
        AR_params[i] <- -A_post[i] - sum(AR_params[i:p])
    }

    # Calculate AR(1) term
    AR_params[1] <- rho - sum(AR_params[2:p]) 
 

    # p out of sample forecasts here
    #Add current period in to incoporate current information set to do BN decomposition
       end_col = nc_X_untransformed-1
    X <-rbind(X_untransformed[2:nr_X_untransformed,1:nc_X_untransformed],cbind(y[nr_y,1],t(X_untransformed[nr_X_untransformed,1:end_col])))
    X_flip <- rbind(X_flip_untransformed[2:nr_X_untransformed,1:nc_X_untransformed], cbind(y_flip[nr_y,1], t(X_flip_untransformed[nr_X_untransformed,1:end_col])))


    Companion <- rbind(t(AR_params), cbind(eye(p - 1), zeros(p - 1, 1)))


    j <- 1
    while (j < p+3){
        Companion_j <- eye(nrow(Companion))
        for(i in 1:j-1){Companion_j <- Companion %*% Companion_j}
        X_j = Companion_j %*% t(t(X[nrow(X),1:ncol(X)]))
        overwrite_row = p+2+nrow(y)+j
        augmented_y[overwrite_row] = t(X_j[1,1:ncol(X_j)])
        X_j = Companion_j %*% t(t(X_flip[nrow(X_flip),1:ncol(X_flip)]))
        augmented_y_flip[overwrite_row] = t(X_j[1,1:ncol(X_j)])
    j = j+1
    }


    end_row0 = p+2
    start_row = p+2+nr_y+1
    end_row = p+2+nr_y+p+2

    augmented_y[1:end_row0] = t(t(rev(t(t(augmented_y_flip[start_row:end_row])))))

    augmented_y_iter = cbind(augmented_y_iter, augmented_y)

    iter = iter + 1
    }
     end_row = p+2+nr_y
    augmented_y = t(t(augmented_y[1:end_row]))

    # Return backcast data
    return (augmented_y)
}



# Implement the BN Filter imposing a selected 'delta'
#
# Note: This function depends on functions implemented above
#
# INPUTS:
# y           Time series to implement BN Filter in first differences
# p           AR lags
# delta       Signal to noise ratio
#
# OUTPUTS:
# BN_cycle    The estimated cycle
# BN_cycle_se The estimated cycle standard error
# aux_out     Auxiliary output (AR coefficients & residuals )
BN_Filter <- function(y, p, delta, dynamic_bands, ib, window, dummies, adjust_bands, compute_stderr = TRUE)
{
    # Do a few preliminary things
    y <- as.matrix(y)
    
    if (adjust_bands){
    tmp_olsvar_outliers <- olsvar_outliers(y = y, p = p, temp_dummies = dummies, nc = FALSE)
    sig2_ols_c = tmp_olsvar_outliers$SIGMA
    }
    else{
    tmp_olsvar <- olsvar(y = y, p = p, nc = FALSE)
    sig2_ols_c = tmp_olsvar$SIGMA	
    }	
    

    # Compute rho from delta
    rho <- -(1.0 - sqrt(delta)) / sqrt(delta)

    # Re-parametrise model and subtract second differences multiplied by 
    # sum of AR coefficients to fix delta. 
    # Adjust the LHS variable to take note of this

    # Backcast data 
    if (ib){
    	augmented_y <- backcast(y,p,delta)
    }
    else{
    	augmented_y <- rbind(zeros(p + 2, 1), y)
    }	

    # Get data matrix from untransformed model. Estimate without a constant
    tmp_olsvar_df <- olsvar_df(y = augmented_y[2:nrow(augmented_y)], p = p, nc = TRUE)
    X <- tmp_olsvar_df$X[ , , drop = FALSE]

    # Get regressors, variance and lags of untransformed AR(p) to perform BN later.
    # Estimate AR by OLS with (for sig_e) and without a constant
    tmp_olsvar <- olsvar(y = augmented_y[3:nrow(augmented_y)], p = p, nc = TRUE)
    sig2_ols = tmp_olsvar$SIGMA
    X_untransformed <- tmp_olsvar$X[ , , drop = FALSE]
    
    reparameterised_y <- augmented_y[(p + 3):nrow(augmented_y), , drop = FALSE] - 
        rho * X[ ,1]

    # Calculate the other preliminaries before estimation
    Y <- reparameterised_y
    X <- X[ ,2:ncol(X)]

    # Set other parameters for conjugate priors
    A_prior <- zeros(p - 1, 1)     # constant and (p-1) phi_star coefficients
    # Uninformative prior on constant centred on zero
    V_prior <- diag(drop(repmat(0.5, 1, p - 1) / (square(seq(from = 1, to = p - 1, by = 1)))))

    # Calculate Posteriors
    V_post <- qr.solve(qr.solve(V_prior,tol = 1e-30) + (1.0 / sig2_ols) * crossprod(x = X, y = NULL),tol = 1e-30)
    A_post <- V_post %*% (qr.solve(V_post,tol = 1e-30) %*% A_prior + t((1.0 / sig2_ols) * X) %*% Y)
    
    # Add rho back into the posterior of estimated parameters
    A_post <- rbind(rho, A_post)

    # Map the parameters in DF-form back to AR-form
    AR_params <- zeros(p, 1)

    # Last lag
    AR_params[p] <- -A_post[p]
    rseq <- seq(from = p - 1, to = 2, by = -1)

    for (i in rseq) {
        AR_params[i] <- -A_post[i] - sum(AR_params[i:p])
    }

    # Calculate AR(1) term
    AR_params[1] <- rho - sum(AR_params[2:p])

    # BN Decomposition starts here
    if (ib){
    X <- rbind(X_untransformed[1:nrow(X_untransformed), , drop = FALSE],
            cbind(y[nrow(y), 1, drop = FALSE], 
                X_untransformed[nrow(X_untransformed), 
                1:(ncol(X_untransformed) - 1), drop = FALSE]))
    }
    else {
    X <- rbind(X_untransformed[2:nrow(X_untransformed), , drop = FALSE],
            cbind(y[nrow(y), 1, drop = FALSE], 
                X_untransformed[nrow(X_untransformed), 
                1:(ncol(X_untransformed) - 1), drop = FALSE])) 	
    }            

    # Compute the BN cycle
    Companion <- rbind(t(AR_params), cbind(eye(p - 1), zeros(p - 1, 1)))
    phi <- -Companion %*% qr.solve(eye(p) - Companion)    
    BN_cycle <- phi %*% t(X)
    BN_cycle <- t(BN_cycle[1, , drop = FALSE])
    
    # Compute the change in trend
    Delta_eta <- (1/(1-rho))*(y - X_untransformed %*% AR_params)
    varDeltaBNtrend <- (1/nrow(y))*t(Delta_eta) %*% Delta_eta

    # Collect results into list object
    result <- list()
    result$BN_cycle <- BN_cycle
    result$aux_out <- list()
    result$aux_out$AR_coeff <- AR_params
    result$aux_out$residuals <- y - X_untransformed %*% AR_params
    result$aux_out$Delta_eta<- Delta_eta
    result$aux_out$varDeltaBNtrend <- varDeltaBNtrend


    # Compute the BN cycle standard error if specified
    if (compute_stderr) {
        A <- Companion
        big_A <- qr.solve(eye(square(nrow(A))) - (A %x% A))
        #var_y <- big_A[1, 1] * sig2e
        vecQ <- zeros(p^2, 1)
        
        ind_vec <- cbind(1.0, zeros(1, p - 1))
        Phi <- (A %*% qr.solve(eye(p) - A))
        
        tobs <- nrow(y)
        
        
  if (dynamic_bands){   
  	
    windp1 <- window + 1
    windm1 <- window - 1

    # Allocate storage for the demeaned series
    BN_cycle_se_t <- matrix(data = NA_real_, nrow = tobs, ncol = 1)
    
        y_temp <- y[1:window]
    	y_temp <- as.matrix(y_temp)
        
        if (adjust_bands){
    	tmp_olsvar_outliers <- olsvar_outliers(y = y_temp, p = p, temp_dummies = dummies, nc = FALSE)
    	sig2_ols_c_t = tmp_olsvar_outliers$SIGMA
    	}
    	else{
    	tmp_olsvar <- olsvar(y = y_temp, p = p, nc = FALSE)
    	sig2_ols_c_t = tmp_olsvar$SIGMA	
   		}	       
    
        vecQ[1,1] = sig2_ols_c_t
        vecSigma_X = big_A %*% vecQ
        Sigma_X_alt = matrix(vecSigma_X,p,p)
        
        BN_cycle_se <- as.numeric(sqrt(ind_vec %*% Phi %*% Sigma_X_alt %*% t(Phi) %*% t(ind_vec)))

    BN_cycle_se_t[1:window] <- BN_cycle_se * ones(window,1)

    	for (i in windp1:tobs) {
    	
    		i_start <- i-window+1
    	    y_temp <- y[i_start:i]
    	    y_temp <- as.matrix(y_temp)
    	    window_dummies <- dummies - i + window        	
        	if (adjust_bands){
    		tmp_olsvar_outliers <- olsvar_outliers(y = y_temp, p = p, temp_dummies = window_dummies, nc = FALSE)
    		sig2_ols_c_t = tmp_olsvar_outliers$SIGMA
    		}
    		else{
    		tmp_olsvar <- olsvar(y = y_temp, p = p, nc = FALSE)
    		sig2_ols_c_t = tmp_olsvar$SIGMA	
   			}	
    
        	vecQ[1,1] = sig2_ols_c_t
        	vecSigma_X = big_A %*% vecQ
        	Sigma_X_alt = matrix(vecSigma_X,p,p)

        	BN_cycle_se <- as.numeric(sqrt(ind_vec %*% Phi %*% Sigma_X_alt %*% t(Phi) %*% t(ind_vec)))
 	
        	BN_cycle_se_t[i] <- BN_cycle_se
    	}
         	
   } 
        else{       
        vecQ[1,1] = sig2_ols_c
        vecSigma_X = big_A %*% vecQ
        Sigma_X_alt = matrix(vecSigma_X,p,p)

        BN_cycle_se <- as.numeric(sqrt(ind_vec %*% Phi %*% Sigma_X_alt %*% t(Phi) %*% t(ind_vec)))
        BN_cycle_se_t <- BN_cycle_se * ones(tobs,1)
        }
        
      if (ib){
      	BN_cycle_se_t <- rbind(BN_cycle_se_t[1],BN_cycle_se_t)
      }  
        
      result$BN_cycle_se <- BN_cycle_se_t    
    }     

    # Return results
    return (result)
}



# Pick the signal to noise ratio 'delta' by doing a grid search 
# from 'd0' in increments of 'dt' to find the first local maximum of
# the amplitude to noise ratio (delta_select=1) or local minimum of variance
# of trend shocks (delta_select=2).
#
# Note: This function depends on functions implemented above
#
# INPUTS:
# y          	Time series
# p          	AR lag order
# d0         	Start grid search here, default is 0.005
# dt         	Increment grid by this amount, default is d0
# option		Selection option (1=max amp-to-noise, 2=min var(trend shocks))
#
# OUTPUTS:
# delta_grid
select_delta <- function(y, p, ib, delta_select, d0, dt)
{
  # Initialise the difference to enter the loop
  diff_t <- 1;

  # Initialise the amplitude to noise ratio
  delta_grid <- d0
  tmp <- BN_Filter(y, p, delta_grid, dynamic_bands, ib, window, dummies, F, F)
  old_amp_to_noise <- var(tmp$BN_cycle) / mean(square(tmp$aux_out$residuals))
  old_precisionDeltaBNtrend <- tmp$aux_out$varDeltaBNtrend

  while (diff_t > 0) {
    delta_grid <- delta_grid + dt
    tmp <- BN_Filter(y, p, delta_grid, dynamic_bands, ib, window, dummies, F, F)
    new_amp_to_noise <- var(tmp$BN_cycle) / mean(square(tmp$aux_out$residuals))
    new_precisionDeltaBNtrend <- tmp$aux_out$varDeltaBNtrend

    if (delta_select == 1){
    	diff_t <- new_amp_to_noise - old_amp_to_noise
    }
    
     	if (delta_select == 2){
    	diff_t <- -(new_precisionDeltaBNtrend - old_precisionDeltaBNtrend)
    }

    if (diff_t > 0) {
        old_amp_to_noise <- new_amp_to_noise
        old_precisionDeltaBNtrend <- new_precisionDeltaBNtrend
    }
  }

  delta_grid <- delta_grid - dt

  return (delta_grid)
}



# bnf using an automatically determined 'delta'
bnf <- function(y, 
                delta_select = 2,
                fixed_delta = 0.01,
                p = 12, 
                d0 = 0.005, 
                dt = 0.0005,              
                demean = c("nd", "sm", "pm", "dm"),
                iterative = 100, 
                dynamic_bands = T,
                adjust_bands = F,
                dummies = c(293, 294),
                window = 40,
                ib = T,
                ...)
  # @delta_select: set to 0 if use fixed delta, 1 if max amp-to-noise, 2 if min var(trend shocks)
  # @fixed_delta: set a fixed delta to be used if delta_select=0, or for delta
  # @p is lag order for estimation
  # @d0: min_delta: start point of grid search for delta if delta_select!=0
  # @dt: delta increment
  # @demean = "nd", "sm", "dm", or "pm", where "nd" = no drift, "sm" = sample mean, "dm" = dynamic demeaning, "pm" = structural breaks, set as 'breaks = c(100, 237)'
  # @iterative: set to >1 for max number of iterations for iterative dynamic demeaning
  # @dynamic_bands: set to T for dynamic error bands, F for fixed standard error bands
  # @adjust_bands: set to T to adjusts for outlier observations when calculating error bands, set outliers as 'dummies = c(293, 294)'
  # @window: rolling window length for dynamic demeaning and/or dynamic error bands (e.g., 40 is 10 years for quarterly data)
  # @ib: set to F if no iterative backcasting as in KMW2018 (just unconditional mean), set to T if iterative backcasting
  # @varargs (...): passed into error bands and piecewise_demean function
{

    if(demean != 'dm' && iterative != 0) stop("Set @iterative to 0 if @demean is not 'dm'.")
    
  
    # Save 'ts' attributes if 'y' is a 'ts' object
    if (is.ts(y)) {
        y_ts <- TRUE
        endy <- end(y)
        freqy <- frequency(y)
    } else {
        y_ts <- FALSE
    }

    # Do a few preliminary things
    y <- as.matrix(y)

    # Compute the first-difference and then demean y
    dy <- diff(x = y, lag = 1)
    
    # Demean dy
    demean <- match.arg(demean)

    if (demean == "sm") {
        demeaned_dy <- scale(x = dy, center = TRUE, scale = FALSE)
        demean_method <- "Sample mean"
    } else if (demean == "pm") {
        demeaned_dy <- piecewise_demean(y = dy, ...)
        demean_method <- "Piecewise mean"
    } else if (demean == "dm") {
        demeaned_dy <- rolling_demean(y = dy, y_cycle = zeros(nrow(dy),1), wind = window)
        demean_method <- "Rolling mean"
    } else { # demean == "nd"
        demeaned_dy <- dy
        demean_method <- "No drift"
    }
    
    

    # Automatically compute delta to pass to BN_Filter
    if (delta_select > 0){
    delta <- select_delta(demeaned_dy, p, ib, delta_select, d0, dt)
    }
    else {
    	delta <- fixed_delta
    	}
    tmp <- BN_Filter(demeaned_dy, p, delta, dynamic_bands, ib, window, dummies, adjust_bands)
    cycle <- tmp$BN_cycle
    
    DeltaBNcycle <- diff(x = cycle, lag = 1)
    
    if (iterative > 0 && demean == "dm"){
    	
    	cycle_iter <- cbind(zeros(nrow(cycle),1),cycle)
    	iter <- 1
    	while (iter < iterative && sd(cycle_iter[,ncol(cycle_iter)]-cycle_iter[,ncol(cycle_iter)-1])>0.001*sd(demeaned_dy)){
      	
      	demeaned_dy <- rolling_demean(y = dy, y_cycle = DeltaBNcycle, wind = window)
      	
        if (delta_select > 0){
          delta <- select_delta(demeaned_dy, p, ib, delta_select, d0, dt)
        }
        else{
        	delta <- fixed_delta
        }
      	
        tmp <- BN_Filter(demeaned_dy, p, delta, dynamic_bands, ib, window, dummies, adjust_bands)
        cycle <- tmp$BN_cycle
        DeltaBNcycle <- diff(x = cycle, lag = 1)
        
        cycle_iter <- cbind(cycle_iter,cycle)		
      		
      	iter = iter + 1
    	}    	
    }
    
    colnames(cycle) <- "Cycle"
    cycle_se <- tmp$BN_cycle_se
    
    # Add 'ts' attribute to 'y' and 'cycle' if 'ts' object was originally passed in
    if (y_ts) {
        y <- ts(data = y, end = endy, frequency = freqy)
        cycle <- ts(data = cycle, end = endy, frequency = freqy)
    }

    # Return result invisibly
    result <- list()
    result$call <- match.call()
    result$y <- y
    result$cycle <- cycle
    result$trend <- y - cycle
    result$cycle_se <- cycle_se
    result$delta <- delta
    result$demean_method <- demean_method
    result$iterative <- iterative
    result$cycle_ci <- round(qnorm(p = 0.05 / 2.0, lower.tail = FALSE), 2) * cycle_se
    if (iterative > 0){
    	result$iterations <- ncol(cycle_iter)-1
    }
    class(result) <- "bnfClass"
    invisible (result)
}

# Print method for class "bnf"
print.bnfClass <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
    cat("\nCall: ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
    
    cat("\n\tIntuitive and Reliable Estimates of the Output Gap from a Beveridge-Nelson Filter\n")
    
    cat(sprintf("\nAutomatic Delta: %2.4f\n", x$delta))

    cat(sprintf("\nDemean method: %s\n", x$demean_method))
    
    if (x$iterative > 0){
    	cat(sprintf("\nIterations: %s\n", x$iterations))
    }
    
    cat("\nCycle:\n")
    print(x$cycle)
    
}

# Plot method for class "bnf"
plot.bnfClass <- function(x,  main = "BN Filter Cycle", 
    plot_ci = TRUE, col = "blue", lwd = 2, ...) 
{
    # All par settings which could be changed
    oldpar <- par(no.readonly = TRUE) 
    
    # Reset plotting device to original settings on exiting
    on.exit(par(oldpar))
   
    # Set up new plotting device
    newpar <- par(mar = c(0, 0, 0, 0), oma = c(4, 4, 3, 2) + 0.1, 
        lty = 1, lwd = 1, xaxs = 'r', yaxs = 'r', mgp = c(3, 1, 0), 
        bty = 'o', tcl = -0.5, cex.main = 1.2, font.main = 2, 
        cex.lab = 1.1, font.lab = 1, cex.axis = 1.1, font.axis = 1)

    # Establish the panel layout -- 2 x 1
    layout(matrix(data = c(1, 2), nrow = 1, ncol = 1))

   

    # Preliminary values - cycle
    if (is.ts(x$cycle)) {
        st <- start(x$cycle)[1] + (start(x$cycle)[2] / 12.0)
        et <- end(x$cycle)[1] + (end(x$cycle)[2] / 12.0)
        #x_axis <- seq(from = st + 0.5, to = et + 0.5, length.out = length(x$cycle))
        x_axis <- seq(from = st + 1/12.0, to = et + 1/12.0, length.out = length(x$cycle))
    } else {
        x_axis <- seq(from = 1, to = length(x$cycle), by = 1)
    }

    c_ub <- ceiling(ifelse(plot_ci, ceiling(max(x$cycle) + x$cycle_ci), max(x$cycle)))
    c_lb <- floor(ifelse(plot_ci, floor(min(x$cycle) - x$cycle_ci), min(x$cycle)))

    # Establish the cycle plotting region
    plot(x$cycle, type = 'n', col = NA, 
        ylim = c(c_lb, c_ub), yaxt = 'n', xaxt = 'n',
        main = "", ylab = "", xlab = "")
    abline(h = 0, col = 1, lwd = 1, lty = 1)
    lines(x$cycle, col = col, lwd = lwd)

    # 95% Confidence Intervals (if requested)
    if (plot_ci) { 
        polygon(x = c(x_axis, rev(x_axis)), y = c((x$cycle - x$cycle_ci), rev(x$cycle + x$cycle_ci)),
            col = adjustcolor(col, alpha.f = 0.2), border = NA)
    }

    # Mark the LHS y-axis
    axis(side = 2, las = 1)

    # Mark the x-axis with dates if available
    if (is.ts(x$cycle)) {
        times <- time(x$cycle)
        ticks <- floor(seq(from = times[1], to = times[length(times)], 
            by = ceiling(length(times) / (10 * frequency(x$cycle)))))
        axis(side = 1, at = ticks)
    } else {
        axis(side = 1)
    }    

    # Individual title for the sub-panel
    mtext(text = "Cycle", side = 3, line = -1.5, outer = FALSE, 
        las = 1, cex = 1.1, font = 2)  
    box(which = "plot", col = 1)

    # Main title
    mtext(text = main, side = 3, line = 1, outer = TRUE, 
        las = 1, cex = 1.5, font = 2)

}

# EOF
