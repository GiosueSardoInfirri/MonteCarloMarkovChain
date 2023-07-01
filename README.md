# MonteCarloMarkovChain

This is the final project of the Advanced Statistics for Physics Analysis course. It is an analysis and implementation of some Adaptive Monte Carlo Markov Chains (MCMC)

## Normal multi-dimensional MCMC

I implemented the Metropolis-Hastings algorithm, as done in class with the coda/rjags packages, but using plain code, which allowed me to:

- enlarge the dimensionality of the generated data ('Normal Metropolis' in the file)
- verify wheter the Gibbs sampling improved the performances or not ('Gibbs Sampling' in the file)

The code used, from now on, is the modification of a simple code that implements the normal 1-dimension Metrpolis-Hastings algorithm, which is made of the following functions:

1. A function to generate the new guessed point (enlarged in many dimensions for all the next programs):

```python
generation_s = function (x0, std) {

    CDF = function (x) {pnorm(x, mean = x0, sd = std)}
    inv_CDF = inverse(CDF, lower = -Inf, upper = Inf)

    random = runif(as.integer(length(x0)))
    new_value = inv_CDF(random)

    return(new_value)
}
```

2. A function to evaluate the ratio between the two distribution in the case of asymmetrical sampling functions:

```python
evaluate_Q = function(x0, std, point) {

    dx = dnorm(point, mean = x0, sd = std)
    sx = dnorm(x0, mean = point, sd = std)

    return(sx/dx)
}
```

3. A function to return the sequence of the samples for a determined wanted distribution, starting from an initial parameter. It returnes the chain values:

```R
random_steps = function (func_wanted, theta_init, n_samples, sigma, print_accept=FALSE) {

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

```

## MCMC using "mvtnorm" package

After writing these functions, I found many R packages able to create and sample from gaussian distributions from the mean and the covariance of the multi-dimensional gaussian. I decided to use the ***mvtnorm*** package, because of its similarity to the normal declaration functions in R. An example is given by the code for retrieving the distribution value (to evaluate Q-ratio):

```R
dmvnorm(point, mean = x0, sigma = cov, log = FALSE)
```

And similarly here I added the possibility for the Gibbs sampling in the multi-dimentional case

## Are these improvements worth it?

In order to see if the effort of improving the function is worthy, I run it in order to find the 3-D gaussian with fixed different variance, timed the program and looked at the acceptance rate, finding the results in Table:

| Function              | Acceptance Rate (%) | Time needed (s) | Coda Autocorrelation |
| --------------------- | ------------------- | --------------- | -------------------- |
| 1D MCMC               | 60                  | 30              | low                  |
| 3D MCMC Metropolis    | 48                  | 130             | low                  |
| 3D MCMC Gibbs         | 59                  | 135             | low                  |
| 3D MCMC mvtnorm       | 50                  | 40              | low                  |
| 3D MCMC mvtnorm Gibbs | 57                  | 140             | low                  |

These results have been retrieved running the program three times with different initial point looking for a total of 100_000 points using a very wide gaussian for the sampling (unchanged) in order to present the strenghts and weakness of the three methods. The table shows that:

- 1D is very performing, but has the limit of the dimensions
- 3D Metropolis has a low acceptance rate: the gaussian for the selection of the new point is very wide, resulting in a small acceptance rate
- 3D Gibbs returns to the higher acceptance rate found in the 1D case (since it moves in 1D)
- 3D mvtnorm gives the possibility of having covariance, but at the cost of time speed and acceptance rate (it is done without Gibbs sampling, in order to spare time and computational resources, otherwise they would scale linearly with the number of dimentions)
- 3D mvtnorm Gibbs, as done in the previous case, increases the acceptance rate at the cost of time

All of them have a low Autocorrelation due to the very wide sampling ditribution, but unluckily all of them suffer from the same problem: no covariance in the selecting distribution. We want an algorithm that keeps the autocorrelation low, but increases the acceptance rate, shaping the sampling distribution at will in order to archieve these goals

**Attention:** I cannot use these algorithms for **adaptive MCMC**, because they don't allow to shape properly the multi-dimensional gaussian, therefore I need an algorithm that uses directly the covariance matrix

## Adaptive MCMC

There are a lot of different implementations for an adaptive MCMC, which all have strenghts and weatknesses that I want to analize: all of them do have a change in the sampling distribution parameters along the chain.

These algorithms are not hard to be implemented in order to have an adaptive MCMC, as presented by Haario et al. $[1]$ and well studied in many works. The idea behind them is to sample from a very wide multi-dimentional distribution at the beginning, then to update this covariance matrix, with always small corrections as well as marginalization factors in order not to have an exploding matrix in time.

#### AM algorithms

The logical schema for the easiest adaptive MCMC algorithm is presented as shown by Andrieu et al. $[2]$ :

1. Initialize $X_0, \mu_0, \Sigma_0$ as initial step, initial mean and initial covariance matrix for the distribution
2. Then take a sample $X_t$ from the sample distribution $q(X_0, \dots, X_{t-1})$ (usually it is taken $N(0, \Sigma)$ as symmetric simple equation), weighted with the ratio between the functions and the Q-ratio $\alpha(X_{t-1}, Y) = \text{max} \left( 1, \frac{f(Y)}{f(X_{t-1})} \right)$ as:

$$
X_t \backsim P^{SRWM}_q = \int  q(z + \mu, \Sigma) dz
$$

    where I report a simplification of the real formula, which most of the physical distributions

4. Update the parameters according to the new value taken from the distribution and a non-decreasing non-converging series $\left[ \gamma_i \right]$, as follows:

$$
\begin{cases}
    \mu_{i+1} = \mu_i + \gamma_{i+1} (X_{i+1} - \mu_i) \\
    \Sigma_{i+1} = \Sigma_i + \gamma_{i+1} ((X_{i+1} - \mu_i) \cdot (X_{i+1} - \mu_i)^T - \Sigma_i)
\end{cases}
$$

4. Repeat points 3. and 4. until maximum number of samples

**Useful adaptation:** In many applications (and mine as well), the mean is corrected recursively and similarly the covariance matrix update $C_i = s_d \cdot \Sigma_i$ follows:

$$
C_t = 
\begin{cases}
    C_0 \hspace{41.5 mm} \text{for} \ t \leq t_0 \\
    s_d \left[\text{cov}(x_0, \dots, x_{t-1}) + \epsilon I_d\right] \hspace{5 mm} \text{for} \ t > t_0 
\end{cases}
$$

And it is found to be $s_d = \frac{2.38}{d^2}$ with $d$ dimension of the sampled distribution $[1]$ and $\epsilon$ is a small parameter, much smaller with respect to the typical dimensionality of the distribution.

#### Rao-Blackwellised algorithm

This algorithm is a very slightly modification of the previous one, except for the fact that, in order to choose the new value, it uses a normal sampling distribution $Y \backsim P^{Rao} = N(X_i, \Sigma_i)$ and the value is accepted with probability $\alpha(X_{t-1}, X_t)$, otherwise take $X_t = X_{t-1}$, as done in the Simple MCMC.

## Bibliography

1. Haario, H., Saksman, E., Tamminen, J.: "An adaptive Metropolis algorithm"
2. Christophe Andrieu, Johannes Thoms: "A tutorial on adaptive MCMC"
3. Marko Laine: "Adaptive MCMC Methods with applications in environmental and geophysical models"
