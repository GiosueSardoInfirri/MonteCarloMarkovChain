# Import dependencies
library(dplyr)
library(scales)
library(GoFKernel)
library(coda)
library(mvtnorm)
library(base)

# STORAGE FOR THE FUNCTIONS USED IN THE ANALISYS (FROM SCRATCK)





# ======================================================================================================================
# FUNCTION FOR EVOLUTION IN THE MCMC ALGORITHMS

    # ANDRIEU IMPLEMENTATION

        # Series for evolution

# Exponential decreasing series (reaching 1/2 at "halved_step")
gamma_series_exp = function (halved_step = 1000, step = 0) {
    tau = halved_step/log(2.)
    return (exp(-(step/tau)))
}
# Reciprocal function (1/x decreasing and reaching 1/2 at "halved_step")
gamma_series_rec = function (halved_step = 1000, step = 0) {
    return (halved_step/(step + halved_step))
}

        # Mean and covariance update

# Easy rule for the update of the mean value
update_mean_andrieu = function (old_mean, new_value, step, gamma_function = gamma_series_exp, halved_step = 1000) {
    new_mean = old_mean + gamma_function(halved_step = halved_step,step)*(new_value - old_mean)
    return (new_mean)
}

# Covariance matrix update
update_cov_andrieu = function (old_cov, old_mean, new_value, step, gamma_function = gamma_series_exp, halved_step = 1000) {
    vector = as.vector(new_value - old_mean)
    new_cov = old_cov + gamma_function(halved_step = halved_step,step)*(outer(vector, vector, FUN = "*") - old_cov)
    return (new_cov)
}

    # HAARIO IMPLEMENTATION

# Easy rule to update the mean at every step of the "evolution", from step = 1
update_mean_haario = function (old_mean, new_value, step) {
    new_mean = (old_mean*step + new_value)/(step+1)
    return (new_mean)
}

# Recursive rule for the evolution of covariance matrix
update_covariance_haario = function (old_cov, old_mean, new_value, new_mean, step, dimentions = 1, s_d = 0, epsilon = 0.001) {
    if (s_d == 0) {s_d = (2.38)**2/dimentions}

    old_mean_matrix = outer(as.vector(old_mean), as.vector(old_mean), FUN = "*")
    new_mean_matrix = outer(as.vector(new_mean), as.vector(new_mean), FUN = "*")
    new_value_matrix = outer(as.vector(new_value), as.vector(new_value), FUN = "*")

    new_cov = (step-1)*old_cov/step + s_d*(step*old_mean_matrix - (step - 1)*new_mean_matrix + new_value_matrix + epsilon*diag(as.integer(dimentions)))/step
    return (new_cov)
}

    # FOR ALL THE IMPLEMENTATIONS: PARAMETERS FROM DISTRIBUTION

sampled_mean = function (series) {
    dim = length(series[1,])
    mean = rep(0, dim)
    for (i in 1:dim) {
        mean[i] = mean(series[,i])
    }
    return (mean)
}

sampled_covariance = function (series, mean) {
    k = length(series[,1])
    sam_cov = diag(0, length(mean))
    for (i in 1:k) {
        vec = as.vector(series[i,])
        addition = outer(vec, vec, FUN = "*")
        sam_cov = sam_cov + addition
    }
    sam_cov = (sam_cov - (k + 1)*outer(mean, mean, FUN = "*"))/k
    return (as.matrix(sam_cov))
}









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
# MULTI-DIMENSIONAL SIMPLE ADAPTIVE MCMC

# This function is ment to return the sequence of the samples for a determined function
random_steps_simple = function (func_wanted, theta_init, n_samples, sigma, print_accept=FALSE, t_0 = 1000,
                         gamma_function = gamma_series_exp, halved_step = 1000) {

    # Then I create the function to extract the random number
    generation_s_simple = function (x0, cov) {

        # I use the library method to generate the new point
        new_value = rmvnorm(1, mean = x0, sigma = cov, method = c("eigen", "svd", "chol"), pre0.9_9994 = FALSE, checkSymmetry = FALSE)

        return(new_value)
    }

    # Initilalizing the parameters
    current_theta = theta_init
    current_function = func_wanted(theta_init)
    samples = matrix(data = NA, nrow = n_samples, ncol = length(theta_init) + 1)

    # For statistical purposes
    accepted = 0
    mean = rep(0, length(theta_init))

    # Evolution loop
    for (n in 1:n_samples) {

        # Take a guessed new theta (s in the slides) and evaluate its probability
        guessed_theta = generation_s_simple(current_theta, sigma)
        guessed_function = func_wanted(guessed_theta)

        # And then update the general conditions
        current_theta = guessed_theta
        current_function = guessed_function

        # Saving the generated samples because R doesn't accept two elements in returns
        samples[n,] = unlist(c(current_function, current_theta))

        # Here I evaluate the sampled mean and covariance 
        if (n == t_0) {
            mean = sampled_mean(series = head(samples[,-1], t_0))
            sigma = sampled_covariance(series = head(samples[,-1], t_0), mean = mean)
        } else if (n > t_0) { # and then at every step I update it with the previously declared functions
            old_mean = mean
            mean = update_mean_andrieu(old_mean = mean, new_value = current_theta, step = n,
                                       gamma_function = gamma_function, halved_step = halved_step)
            sigma = update_cov_andrieu(old_cov = sigma, old_mean = old_mean, new_value = current_theta,
                                       step = n, gamma_function = gamma_function, halved_step = halved_step)
        }

    }

    if(print_accept) {
        cat("Acceptance rate = ", round(accepted/n_samples*100, 5), '%\n')
        cat("Final mean = ", mean, "\n")
        cat("Final covariance matrix = \n")
        print(sigma)
    }

    return(samples)
}








# ======================================================================================================================
# MULTI-DIMENSIONAL SIMPLE ADAPTIVE MCMC WITH GIBBS SAMPLING

# This function is ment to return the sequence of the samples for a determined function
random_steps_simple_gibbs = function (func_wanted, theta_init, n_samples, sigma, print_accept=FALSE, t_0 = 1000,
                         gamma_function = gamma_series_exp, halved_step = 1000) {

    # Then I create the function to extract the random number
    generation_s_simple_gibbs = function (x0, cov) {

        # I use the library method to generate the new point
        new_value = rmvnorm(1, mean = x0, sigma = cov, method = c("eigen", "svd", "chol"), pre0.9_9994 = FALSE, checkSymmetry = FALSE)

        return(new_value)
    }

    # Initilalizing the parameters
    current_theta = theta_init
    current_function = func_wanted(theta_init)
    dimensions = length(theta_init)
    samples = matrix(data = NA, nrow = n_samples, ncol = length(theta_init) + 1)

    # For statistical purposes
    accepted = 0
    mean = rep(0, length(theta_init))

    # Evolution loop
    for (n in 1:n_samples) {

        guessed_theta = current_theta

        # I then can loop on the dimensions of the distribution to allpy the gibbs sampling
        for (dim in 1:dimensions) {

            # Take a guessed new theta (s in the slides) and evaluate its probability
            guessed_theta[dim] = generation_s_simple_gibbs(current_theta, sigma)[dim]
            guessed_function = func_wanted(guessed_theta)

            # And then update the general conditions
            current_theta = guessed_theta
            current_function = guessed_function    
        }

        # Saving the generated samples because R doesn't accept two elements in returns
        samples[n,] = unlist(c(current_function, current_theta))

        # Here I evaluate the sampled mean and covariance 
        if (n == t_0) {
            mean = sampled_mean(series = head(samples[,-1], t_0))
            sigma = sampled_covariance(series = head(samples[,-1], t_0), mean = mean)
        } else if (n > t_0) { # and then at every step I update it with the previously declared functions
            old_mean = mean
            mean = update_mean_andrieu(old_mean = mean, new_value = current_theta, step = n,
                                       gamma_function = gamma_function, halved_step = halved_step)
            sigma = update_cov_andrieu(old_cov = sigma, old_mean = old_mean, new_value = current_theta,
                                       step = n, gamma_function = gamma_function, halved_step = halved_step)
        }


    }

    if(print_accept) {
        cat("Acceptance rate = ", round(accepted/(n_samples*dimensions)*100, 5), '%\n')
        cat("Final mean = ", mean, "\n")
        cat("Final covariance matrix = \n")
        print(sigma)
    }

    return(samples)
}









# ======================================================================================================================
# MULTI-DIMENSIONAL HAARIO MCMC

# This function is ment to return the sequence of the samples for a determined function
random_steps_haario = function (func_wanted, theta_init, n_samples, sigma, print_accept=FALSE, t_0 = 1000, eps = 0.001) {

    # Then I create the function to extract the random number
    generation_s_haario = function (x0, cov) {

        # I use the library method to generate the new point
        new_value = rmvnorm(1, mean = x0, sigma = cov, method = c("eigen", "svd", "chol"), pre0.9_9994 = FALSE, checkSymmetry = FALSE)

        return(new_value)
    }

    # And the one to check the value of the quantiles
    evaluate_Q_haario = function(x0, cov, point) {

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
    mean = rep(0, length(theta_init))

    # Evolution loop
    for (n in 1:n_samples) {

        # Take a guessed new theta (s in the slides) and evaluate its probability
        guessed_theta = generation_s_haario(current_theta, sigma)
        guessed_function = func_wanted(guessed_theta)

        # Acceptance conditions
        Q_ratio = 1 #evaluate_Q_haario(current_theta, sigma, guessed_theta)
        rho = guessed_function/current_function * Q_ratio
        # cat(guessed_theta, guessed_function, Q_ratio, rho, "\n")

        # And then update the general conditions
        if (rho > runif(1)) {
            current_theta = guessed_theta
            current_function = guessed_function
            accepted = accepted + 1
        } # else they remain unchanged and then loaded direcctly in the solutions vector

        # Saving the generated samples because R doesn't accept two elements in returns
        samples[n,] = unlist(c(current_function, current_theta))

        # Here I evaluate the sampled mean and covariance 
        if (n == t_0) {
            mean = sampled_mean(series = head(samples[,-1], t_0))
            sigma = sampled_covariance(series = head(samples[,-1], t_0), mean = mean)
        } else if (n > t_0) { # and then at every step I update it with the previously declared functions
            old_mean = mean
            mean = update_mean_haario (old_mean = mean, new_value = current_theta, step = n)
            sigma = update_covariance_haario(old_cov = sigma, old_mean = old_mean, new_value = current_theta, new_mean = mean,
                                             step = n, dimentions = length(theta_init), epsilon = eps)
        }

    }

    if(print_accept) {
        cat("Acceptance rate = ", round(accepted/n_samples*100, 5), '%\n')
        cat("Final mean = ", mean, "\n")
        cat("Final covariance matrix = \n")
        print(sigma)
    }

    return(samples)
}









# ======================================================================================================================
# MULTI-DIMENSIONAL HAARIO MCMC WITH GIBBS SAMPLING

# This function is ment to return the sequence of the samples for a determined function
random_steps_haario_gibbs = function (func_wanted, theta_init, n_samples, sigma, print_accept=FALSE, t_0 = 1000, eps = 0.001) {

    # Then I create the function to extract the random number
    generation_s_haario_gibbs = function (x0, cov) {

        # I use the library method to generate the new point
        new_value = rmvnorm(1, mean = x0, sigma = cov, method = c("eigen", "svd", "chol"), pre0.9_9994 = FALSE, checkSymmetry = FALSE)

        return(new_value)
    }

    # And the one to check the value of the quantiles
    evaluate_Q_haario_gibbs = function(x0, cov, point) {

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
    mean = rep(0, length(theta_init))

    # Evolution loop
    for (n in 1:n_samples) {

        guessed_theta = current_theta

        # I then can loop on the dimensions of the distribution to allpy the gibbs sampling
        for (dim in 1:dimensions) {

            # Take a guessed new theta (s in the slides) and evaluate its probability
            guessed_theta[dim] = generation_s_haario_gibbs(current_theta, sigma)[dim]
            guessed_function = func_wanted(guessed_theta)

            # Acceptance conditions
            Q_ratio = 1#evaluate_Q_haario_gibbs(current_theta, sigma, guessed_theta)
            rho = guessed_function/current_function * Q_ratio
            # cat(guessed_theta, guessed_function, Q_ratio, rho, "\n")

            # And then update the general conditions
            if (rho > runif(1)) {
                current_theta[dim] = guessed_theta[dim]
                current_function = guessed_function
                accepted = accepted + 1
            } # else they remain unchanged and then loaded direcctly in the solutions vector
    
        }

        # Saving the generated samples because R doesn't accept two elements in returns
        samples[n,] = unlist(c(current_function, current_theta))

        # Here I evaluate the sampled mean and covariance 
        if (n == t_0) {
            mean = sampled_mean(series = head(samples[,-1], t_0))
            sigma = sampled_covariance(series = head(samples[,-1], t_0), mean = mean)
        } else if (n > t_0) { # and then at every step I update it with the previously declared functions
            old_mean = mean
            mean = update_mean_haario (old_mean = mean, new_value = current_theta, step = n)
            sigma = update_covariance_haario(old_cov = sigma, old_mean = old_mean, new_value = current_theta, new_mean = mean,
                                             step = n, dimentions = length(theta_init), epsilon = eps)
        }


    }

    if(print_accept) {
        cat("Acceptance rate = ", round(accepted/(n_samples*dimensions)*100, 5), '%\n')
        cat("Final mean = ", mean, "\n")
        cat("Final covariance matrix = \n")
        print(sigma)
    }

    return(samples)
}










# ======================================================================================================================
# MULTI-DIMENSIONAL RAO MCMC

# This function is ment to return the sequence of the samples for a determined function
random_steps_AM_rao = function (func_wanted, theta_init, n_samples, sigma, print_accept=FALSE, t_0 = 1000,
                         gamma_function = gamma_series_exp, halved_step = 1000) {

    # Then I create the function to extract the random number
    generation_s_AM_rao = function (x0, cov) {

        # I use the library method to generate the new point
        new_value = rmvnorm(1, mean = x0, sigma = cov, method = c("eigen", "svd", "chol"), pre0.9_9994 = FALSE, checkSymmetry = FALSE)

        return(new_value)
    }

    # And the one to check the value of the quantiles
    evaluate_Q_AM_rao = function(x0, cov, point) {

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
    mean = rep(0, length(theta_init))

    # Evolution loop
    for (n in 1:n_samples) {

        # Take a guessed new theta (s in the slides) and evaluate its probability
        guessed_theta = generation_s_AM_rao(current_theta, sigma)
        guessed_function = func_wanted(guessed_theta)

        # Acceptance conditions
        Q_ratio = 1 #evaluate_Q_AM_rao(current_theta, sigma, guessed_theta)
        rho = guessed_function/current_function * Q_ratio
        # cat(guessed_theta, guessed_function, Q_ratio, rho, "\n")

        # And then update the general conditions
        if (rho > runif(1)) {
            current_theta = guessed_theta
            current_function = guessed_function
            accepted = accepted + 1
        } # else they remain unchanged and then loaded direcctly in the solutions vector

        # Saving the generated samples because R doesn't accept two elements in returns
        samples[n,] = unlist(c(current_function, current_theta))

        # Here I evaluate the sampled mean and covariance 
        if (n == t_0) {
            mean = sampled_mean(series = head(samples[,-1], t_0))
            sigma = sampled_covariance(series = head(samples[,-1], t_0), mean = mean)
        } else if (n > t_0) { # and then at every step I update it with the previously declared functions
            old_mean = mean
            mean = update_mean_andrieu(old_mean = mean, new_value = current_theta, step = n,
                                       gamma_function = gamma_function, halved_step = halved_step)
            sigma = update_cov_andrieu(old_cov = sigma, old_mean = old_mean, new_value = current_theta,
                                       step = n, gamma_function = gamma_function, halved_step = halved_step)
        }

    }

    if(print_accept) {
        cat("Acceptance rate = ", round(accepted/n_samples*100, 5), '%\n')
        cat("Final mean = ", mean, "\n")
        cat("Final covariance matrix = \n")
        print(sigma)
    }

    return(samples)
}









# ======================================================================================================================
# MULTI-DIMENSIONAL RAO MCMC WITH GIBBS SAMPLING

# This function is ment to return the sequence of the samples for a determined function
random_steps_AM_rao_gibbs = function (func_wanted, theta_init, n_samples, sigma, print_accept=FALSE, t_0 = 1000,
                         gamma_function = gamma_series_exp, halved_step = 1000) {

    # Then I create the function to extract the random number
    generation_s_AM_rao_gibbs = function (x0, cov) {

        # I use the library method to generate the new point
        new_value = rmvnorm(1, mean = x0, sigma = cov, method = c("eigen", "svd", "chol"), pre0.9_9994 = FALSE, checkSymmetry = FALSE)

        return(new_value)
    }

    # And the one to check the value of the quantiles
    evaluate_Q_AM_rao_gibbs = function(x0, cov, point) {

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
    mean = rep(0, length(theta_init))

    # Evolution loop
    for (n in 1:n_samples) {

        guessed_theta = current_theta

        # I then can loop on the dimensions of the distribution to allpy the gibbs sampling
        for (dim in 1:dimensions) {

            # Take a guessed new theta (s in the slides) and evaluate its probability
            guessed_theta[dim] = generation_s_AM_rao_gibbs(current_theta, sigma)[dim]
            guessed_function = func_wanted(guessed_theta)

            # Acceptance conditions
            Q_ratio = 1#evaluate_Q_AM_rao_gibbs(current_theta, sigma, guessed_theta)
            rho = guessed_function/current_function * Q_ratio
            # cat(guessed_theta, guessed_function, Q_ratio, rho, "\n")

            # And then update the general conditions
            if (rho > runif(1)) {
                current_theta[dim] = guessed_theta[dim]
                current_function = guessed_function
                accepted = accepted + 1
            } # else they remain unchanged and then loaded direcctly in the solutions vector
    
        }

        # Saving the generated samples because R doesn't accept two elements in returns
        samples[n,] = unlist(c(current_function, current_theta))

        # Here I evaluate the sampled mean and covariance 
        if (n == t_0) {
            mean = sampled_mean(series = head(samples[,-1], t_0))
            sigma = sampled_covariance(series = head(samples[,-1], t_0), mean = mean)
        } else if (n > t_0) { # and then at every step I update it with the previously declared functions
            old_mean = mean
            mean = update_mean_andrieu(old_mean = mean, new_value = current_theta, step = n,
                                       gamma_function = gamma_function, halved_step = halved_step)
            sigma = update_cov_andrieu(old_cov = sigma, old_mean = old_mean, new_value = current_theta,
                                       step = n, gamma_function = gamma_function, halved_step = halved_step)
        }


    }

    if(print_accept) {
        cat("Acceptance rate = ", round(accepted/(n_samples*dimensions)*100, 5), '%\n')
        cat("Final mean = ", mean, "\n")
        cat("Final covariance matrix = \n")
        print(sigma)
    }

    return(samples)
}










# ======================================================================================================================
# MULTI-DIMENSIONAL GLOBAL ADAPTIVE MCMC

# This function is ment to return the sequence of the samples for a determined function
random_steps_global = function (func_wanted, theta_init, n_samples, sigma, print_accept=FALSE, t_0 = 1000,
                         gamma_function = gamma_series_exp, halved_step = 1000) {

    # Then I create the function to extract the random number
    generation_s_global = function (x0, cov) {

        # I use the library method to generate the new point
        new_value = rmvnorm(1, mean = x0, sigma = cov, method = c("eigen", "svd", "chol"), pre0.9_9994 = FALSE, checkSymmetry = FALSE)

        return(new_value)
    }

    # And the one to check the value of the quantiles
    evaluate_Q_global = function(x0, cov, point) {

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
    rho_mean = 0 # The mean of the alphas in README
    lambda = 1
    mean = rep(0, length(theta_init))

    # Evolution loop
    for (n in 1:n_samples) {

        # Take a guessed new theta (s in the slides) and evaluate its probability
        used_sigma = exp(lambda)*sigma
        if (!isSymmetric(used_sigma)) {
            print("FALSE")
        }
        guessed_theta = generation_s_global(current_theta, used_sigma)
        guessed_function = func_wanted(guessed_theta)

        # Acceptance conditions
        Q_ratio = 1 #evaluate_Q_global(current_theta, sigma, guessed_theta)
        rho = guessed_function/current_function * Q_ratio
        # cat(guessed_theta, guessed_function, Q_ratio, rho, "\n")

        # And then update the general conditions
        if (rho > runif(1)) {
            current_theta = guessed_theta
            current_function = guessed_function
            accepted = accepted + 1
        } # else they remain unchanged and then loaded direcctly in the solutions vector

        # Saving the generated samples because R doesn't accept two elements in returns
        samples[n,] = unlist(c(current_function, current_theta))

        # Here I evaluate the sampled mean and covariance
        rho_mean = (rho_mean*n + rho)/(n + 1)
        if (n == t_0) {
            mean = sampled_mean(series = head(samples[,-1], t_0))
            sigma = sampled_covariance(series = head(samples[,-1], t_0), mean = mean)
        } else if (n > t_0) { # and then at every step I update it with the previously declared functions
            old_mean = mean
            mean = update_mean_andrieu(old_mean = mean, new_value = current_theta, step = n,
                                       gamma_function = gamma_function, halved_step = halved_step)
            sigma = update_cov_andrieu(old_cov = sigma, old_mean = old_mean, new_value = current_theta,
                                       step = n, gamma_function = gamma_function, halved_step = halved_step)
            lambda = lambda + gamma_function(halved_step = halved_step, step = n)*(rho - rho_mean) # it is the logarithm
        }

    }

    if(print_accept) {
        cat("Acceptance rate = ", round(accepted/n_samples*100, 5), '%\n')
        cat("Final mean = ", mean, "\n")
        cat("Final lambda = ", lambda, "\n")
        cat("Final covariance matrix = \n")
        print(sigma)
    }

    return(samples)
}








# ======================================================================================================================
# MULTI-DIMENSIONAL GLOBAL ADAPTIVE MCMC WITH GIBBS SAMPLING

# This function is ment to return the sequence of the samples for a determined function
random_steps_global_gibbs = function (func_wanted, theta_init, n_samples, sigma, print_accept=FALSE, t_0 = 1000,
                         gamma_function = gamma_series_exp, halved_step = 1000) {

    # Then I create the function to extract the random number
    generation_s_global_gibbs = function (x0, cov) {

        # I use the library method to generate the new point
        new_value = rmvnorm(1, mean = x0, sigma = cov, method = c("eigen", "svd", "chol"), pre0.9_9994 = FALSE, checkSymmetry = FALSE)

        return(new_value)
    }

    # And the one to check the value of the quantiles
    evaluate_Q_global_gibbs = function(x0, cov, point) {

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
    rho_mean = 0 # The mean of the alphas in README
    lambda = 1
    mean = rep(0, length(theta_init))

    # Evolution loop
    for (n in 1:n_samples) {

        # Take a guessed new theta (s in the slides) and evaluate its probability
        used_sigma = exp(lambda)*sigma
        if (!isSymmetric(used_sigma)) {
            print("FALSE")
        }

        guessed_theta = current_theta

        # I then can loop on the dimensions of the distribution to allpy the gibbs sampling
        for (dim in 1:dimensions) {

            # Take a guessed new theta (s in the slides) and evaluate its probability
            guessed_theta[dim] = generation_s_global_gibbs(current_theta, used_sigma)[dim]
            guessed_function = func_wanted(guessed_theta)

            # Acceptance conditions
            Q_ratio = 1 #evaluate_Q_global_gibbs(current_theta, sigma, guessed_theta)
            rho = guessed_function/current_function * Q_ratio
            # cat(guessed_theta, guessed_function, Q_ratio, rho, "\n")

            # And then update the general conditions
            if (rho > runif(1)) {
                current_theta[dim] = guessed_theta[dim]
                current_function = guessed_function
                accepted = accepted + 1
            } # else they remain unchanged and then loaded direcctly in the solutions vector
    
        }

        # Saving the generated samples because R doesn't accept two elements in returns
        samples[n,] = unlist(c(current_function, current_theta))

        # Here I evaluate the sampled mean and covariance
        rho_mean = (rho_mean*n + rho)/(n + 1)
        if (n == t_0) {
            mean = sampled_mean(series = head(samples[,-1], t_0))
            sigma = sampled_covariance(series = head(samples[,-1], t_0), mean = mean)
        } else if (n > t_0) { # and then at every step I update it with the previously declared functions
            old_mean = mean
            mean = update_mean_andrieu(old_mean = mean, new_value = current_theta, step = n,
                                       gamma_function = gamma_function, halved_step = halved_step)
            sigma = update_cov_andrieu(old_cov = sigma, old_mean = old_mean, new_value = current_theta,
                                       step = n, gamma_function = gamma_function, halved_step = halved_step)
            lambda = lambda + gamma_function(halved_step = halved_step, step = n)*(rho - rho_mean) # it is the logarithm
        }

    }

    if(print_accept) {
        cat("Acceptance rate = ", round(accepted/(n_samples*dimensions)*100, 5), '%\n')
        cat("Final mean = ", mean, "\n")
        cat("Final lambda = ", lambda, "\n")
        cat("Final covariance matrix = \n")
        print(sigma)
    }

    return(samples)
}










# ======================================================================================================================
# GRAPHICAL STUFF

show_results = function (mcmc, init, std, step = 0, show_histo = TRUE, show_chain = TRUE, show_lags = TRUE) {

    # check for the dimension of the distribution
    lungh = as.integer(length(mcmc[1,]))
    num_plots = as.integer(show_histo + show_chain + show_lags)
    if (step == 0) {step = as.integer(length(mcmc[,1])/100)} # It is prefixed to 100 total shown steps

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

