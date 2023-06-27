# Import dependencies
library(dplyr)
library(scales)
library(GoFKernel)
library(coda)
library(mvtnorm)

# STORAGE FOR THE FUNCTIONS USED IN THE ANALISYS (FROM SCRATCK)





# ======================================================================================================================
# ONE-DIMENSIONAL MCMC

# This function is ment to return the sequence of the samples for a determined function
random_steps_1D = function (func_wanted, theta_init, n_samples, sigma, print_accept=FALSE) {

    # I create the function to extract the random number
    generation_s = function (x0, std) {

        # Like in the generation of points I use inverse of cumulative distribution
        CDF = function (x) {pnorm(x, mean = x0, sd = std)}
        inv_CDF = inverse(CDF, lower = -Inf, upper = Inf)

        # I create the cumulative function (for the normal distribution I have it for free in R)
        random = runif(as.integer(length(x0)))
        new_value = inv_CDF(random)

        return(new_value)
    }

    # And the one to check the value of the quantiles
    evaluate_Q = function(x0, std, point) {

        # For the case of the normal distribution
        dx = dnorm(point, mean = x0, sd = std)
        sx = dnorm(x0, mean = point, sd = std)

        return(sx/dx)
    }

    # Initilalizing the parameters
    current_theta = theta_init
    current_function = func_wanted(theta_init)
    samples = matrix(data = NA, nrow = n_samples, ncol = length(theta_init) + 1)

    # For statistical purposes
    accepted = 0

    # Evolution loop
    for (n in 1:n_samples) {

        # Take a guessed new theta (s in the slides) and evaluate its probability
        guessed_theta = generation_s(current_theta, sigma)
        guessed_function = func_wanted(guessed_theta)

        # Acceptance conditions
        Q_ratio = evaluate_Q(current_theta, sigma, guessed_theta)
        rho = guessed_function/current_function * Q_ratio
        # cat(guessed_theta, guessed_function, Q_ratio, rho)

        # And then update the general conditions
        if (rho > runif(1)) {
            current_theta = guessed_theta
            current_function = guessed_function
            accepted = accepted + 1
        } # else they remain unchanged and then loaded direcctly in the solutions vector

        # Saving the generated samples because R doesn't accept two elements in returns
        samples[n,] = unlist(c(current_function, current_theta))

    }

    if(print_accept) {
        cat("Acceptance rate = ", round(accepted/n_samples*100, 5), '%\n')
    }

    return(samples)
}








# ======================================================================================================================
# MULTI-DIMENSIONAL MCMC

# This function is ment to return the sequence of the samples for a determined function
random_steps_multi = function (func_wanted, theta_init, n_samples, sigma, print_accept=FALSE) {

    # I create the function to extract the random number
    generation_s_multi = function (x0, std) {

        # I have to loop on the dimensions, in order to use the CDF on all the dimensions
        dimension = as.integer(length(x0))
        new_value = rep(0, dimension)

        for (dim in 1:dimension) {

            # Like in the generation of points I use inverse of cumulative distribution
            CDF = function (x) {pnorm(x, mean = x0[dim], sd = std[dim])}
            inv_CDF = inverse(CDF, lower = -Inf, upper = Inf)

            # I create the cumulative function (for the normal distribution I have it for free in R)
            random = runif(as.integer(length(x0[dim])))
            new_value[dim] = inv_CDF(random)

        }

        return(new_value)
    }

    # And the one to check the value of the quantiles
    evaluate_Q_multi = function(x0, std, point) {

        # For the case of the normal distribution
        dx = dnorm(point, mean = x0, sd = std)
        sx = dnorm(x0, mean = point, sd = std)

        return(prod(sx/dx))
    }

    # Initilalizing the parameters
    current_theta = theta_init
    current_function = func_wanted(theta_init)
    samples = matrix(data = NA, nrow = n_samples, ncol = length(theta_init) + 1)

    # For statistical purposes
    accepted = 0

    # Evolution loop
    for (n in 1:n_samples) {

        # Take a guessed new theta (s in the slides) and evaluate its probability
        guessed_theta = generation_s_multi(current_theta, sigma)
        guessed_function = func_wanted(guessed_theta)

        # Acceptance conditions
        Q_ratio = evaluate_Q_multi(current_theta, sigma, guessed_theta)
        rho = guessed_function/current_function * Q_ratio

        # And then update the general conditions
        if (rho > runif(1)) {
            current_theta = guessed_theta
            current_function = guessed_function
            accepted = accepted + 1
        } # else they remain unchanged and then loaded direcctly in the solutions vector

        # Saving the generated samples because R doesn't accept two elements in returns
        samples[n,] = unlist(c(current_function, current_theta))

    }

    if(print_accept) {
        cat("Acceptance rate = ", round(accepted/n_samples*100, 5), '%\n')
    }

    return(samples)
}









# ======================================================================================================================
# MULTI-DIMENSIONAL MCMC WITH GIBBS SAMPLING

# This function is ment to return the sequence of the samples for a determined function
random_steps_gibbs = function (func_wanted, theta_init, n_samples, sigma, print_accept=FALSE) {

    # I create the function to extract the random number
    generation_s_gibbs = function (x0, std, dim) {

        # Like in the generation of points I use inverse of cumulative distribution
        CDF = function (x) {pnorm(x, mean = x0[dim], sd = std[dim])}
        inv_CDF = inverse(CDF, lower = -Inf, upper = Inf)

        # I create the cumulative function (for the normal distribution I have it for free in R)
        random = runif(1)
        new_value = inv_CDF(random)

        return(new_value)
    }

    # And the one to check the value of the quantiles
    evaluate_Q_gibbs = function(x0, std, point) {

        # For the case of the normal distribution
        dx = dnorm(point, mean = x0, sd = std)
        sx = dnorm(x0, mean = point, sd = std)

        return(prod(sx/dx))
    }

    # Initilalizing the parameters
    current_theta = theta_init
    current_function = func_wanted(theta_init)
    dimensions = length(theta_init)
    samples = matrix(data = NA, nrow = n_samples, ncol = dimensions + 1)

    # For statistical purposes
    accepted = 0

    # Evolution loop
    for (n in 1:n_samples) {

        guessed_theta = current_theta

        # I then can loop on the dimensions of the distribution to allpy the gibbs sampling
        for (dim in 1:dimensions) {

            # Take a guessed new theta (s in the slides) and evaluate its probability
            guessed_theta[dim] = generation_s_gibbs(current_theta, sigma, dim)
            guessed_function = func_wanted(guessed_theta)

            # Acceptance conditions
            Q_ratio = evaluate_Q_gibbs(current_theta, sigma, guessed_theta)
            rho = guessed_function/current_function * Q_ratio
            # cat(guessed_theta, "\t", guessed_function, "\t", Q_ratio, "\t", rho, '\n')

            # And then update the general conditions
            if (rho > runif(1)) {
                current_theta[dim] = guessed_theta[dim]
                current_function = guessed_function
                accepted = accepted + 1
            } # else they remain unchanged and then loaded direcctly in the solutions vector

        }

        # Saving the generated samples because R doesn't accept two elements in returns
        samples[n,] = unlist(c(current_function, current_theta))

    }

    if(print_accept) {
        cat("Acceptance rate = ", round(accepted/(n_samples*dimensions)*100, 5), '%\n')
    }

    return(samples)
}









# ======================================================================================================================
# MULTI-DIMENSIONAL MCMC MVTNORM

# This function is ment to return the sequence of the samples for a determined function
random_steps_mvtnorm = function (func_wanted, theta_init, n_samples, sigma, print_accept=FALSE) {

    # I create the function to extract the random number
    generation_s_mvtnorm = function (x0, cov) {

        # I use the library method to generate the new point
        new_value = rmvnorm(1, mean = x0, sigma = cov, method = c("eigen", "svd", "chol"), pre0.9_9994 = FALSE, checkSymmetry = FALSE)

        return(new_value)
    }

    # And the one to check the value of the quantiles
    evaluate_Q_mvtnorm = function(x0, cov, point) {

        # For the case of the normal distribution
        sx = dmvnorm(point, mean = x0, sigma = cov, log = FALSE)
        dx = dmvnorm(x0, mean = point, sigma = cov, log = FALSE)

        return(sx/dx)
    }

    # Initilalizing the parameters
    current_theta = theta_init
    current_function = func_wanted(theta_init)
    samples = matrix(data = NA, nrow = n_samples, ncol = length(theta_init) + 1)

    # For statistical purposes
    accepted = 0

    # Evolution loop
    for (n in 1:n_samples) {

        # Take a guessed new theta (s in the slides) and evaluate its probability
        guessed_theta = generation_s_mvtnorm(current_theta, sigma)
        guessed_function = func_wanted(guessed_theta)

        # Acceptance conditions
        Q_ratio = evaluate_Q_mvtnorm(current_theta, sigma, guessed_theta)
        rho = guessed_function/current_function * Q_ratio
        # cat(guessed_theta, guessed_function, Q_ratio, rho)

        # And then update the general conditions
        if (rho > runif(1)) {
            current_theta = guessed_theta
            current_function = guessed_function
            accepted = accepted + 1
        } # else they remain unchanged and then loaded direcctly in the solutions vector

        # Saving the generated samples because R doesn't accept two elements in returns
        samples[n,] = unlist(c(current_function, current_theta))

    }

    if(print_accept) {
        cat("Acceptance rate = ", round(accepted/n_samples*100, 5), '%\n')
    }

    return(samples)
}










# ======================================================================================================================
# MULTI-DIMENSIONAL MCMC MVTNORM WITH GIBBS SAMPLING

# This function is ment to return the sequence of the samples for a determined function
random_steps_mvtnorm_gibbs = function (func_wanted, theta_init, n_samples, sigma, print_accept=FALSE) {

    # I create the function to extract the random number
    generation_s_mvtnorm_gibbs = function (x0, cov) {

        # I use the library method to generate the new point
        new_value = rmvnorm(1, mean = x0, sigma = cov, method = c("eigen", "svd", "chol"), pre0.9_9994 = FALSE, checkSymmetry = FALSE)

        return(new_value)
    }

    # And the one to check the value of the quantiles
    evaluate_Q_mvtnorm_gibbs = function(x0, cov, point) {

        # For the case of the normal distribution
        sx = dmvnorm(point, mean = x0, sigma = cov, log = FALSE)
        dx = dmvnorm(x0, mean = point, sigma = cov, log = FALSE)

        return(sx/dx)
    }

    # Initilalizing the parameters
    current_theta = theta_init
    current_function = func_wanted(theta_init)
    dimensions = length(theta_init)
    samples = matrix(data = NA, nrow = n_samples, ncol = length(theta_init) + 1)

    # For statistical purposes
    accepted = 0

    # Evolution loop
    for (n in 1:n_samples) {

        guessed_theta = current_theta

        # I then can loop on the dimensions of the distribution to allpy the gibbs sampling
        for (dim in 1:dimensions) {

            # Take a guessed new theta (s in the slides) and evaluate its probability
            guessed_theta[dim] = generation_s_mvtnorm_gibbs(current_theta, sigma)[dim]
            guessed_function = func_wanted(guessed_theta)

            # Acceptance conditions
            Q_ratio = evaluate_Q_mvtnorm_gibbs(current_theta, sigma, guessed_theta)
            rho = guessed_function/current_function * Q_ratio
            # cat(guessed_theta, guessed_function, Q_ratio, rho)

            # And then update the general conditions
            if (rho > runif(1)) {
                current_theta = guessed_theta
                current_function = guessed_function
                accepted = accepted + 1
            } # else they remain unchanged and then loaded direcctly in the solutions vector

            # Saving the generated samples because R doesn't accept two elements in returns
            samples[n,] = unlist(c(current_function, current_theta))
        
        }

    }

    if(print_accept) {
        cat("Acceptance rate = ", round(accepted/(n_samples*dimensions)*100, 5), '%\n')
    }

    return(samples)
}










# ======================================================================================================================
# GRAPHICAL STUFF

show_results = function (mcmc, init, std, step = 0, show_histo = TRUE, show_chain = TRUE, show_lags = TRUE) {

    # check for the dimension of the distribution
    lungh = as.integer(length(mcmc[1,]))
    num_plots = as.integer(show_histo + show_chain + show_lags)
    if (step == 0) {step = as.integer(length(mcmc[,1])/100)}

    # Plotting every dimension of the plot
    par(mfrow=c(lungh-1, num_plots), oma = c(0, 0, 0, 0))
    options(repr.plot.width=20, repr.plot.height=6*(lungh-1))

    for (dim in 2:lungh) {
        # The chain plot
        plot_g = mcmc[seq(0, length(mcmc[,1]), step), dim]
        plot(1:length(plot_g), plot_g, type = 'o', lwd = 1, col = 'black', xlab = 'Iteration', ylab = 'Parameter evolution', 
            main = paste0('Initial value = ', list(init), '     Standard deviation = ', list(std)))

        # Then the histogram
        hist(mcmc[,dim], breaks = as.integer(sqrt(length(mcmc[,1]))), xlab = 'Lambda', ylab = 'Counts',
            main = paste('Histogram of the posterior distribution in the dimension', dim-1))

        # And finally the coda analysis of the chain
        g_chain = as.mcmc(mcmc[,dim])
        lags = seq(1,length(g_chain),length.out = 100)
        auto_g = autocorr(g_chain, lags=lags)
        if (length(lags) > length(auto_g)) {lags = head(lags, -1)}
        plot(lags, auto_g, ylim=c(0,1), pch=7, col="navy", xlab="lag", ylab="ACF", cex=1.5, main = paste("Autocorrelation behaviour with sigma =", list(std)))
    }

}

