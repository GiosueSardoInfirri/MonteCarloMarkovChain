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

#### Adaptive Metropolis algorithms (AM)

The logical schema for the easiest adaptive MCMC algorithm is presented as shown by Andrieu et al. $[2]$ :

1. Initialize $X_0, \mu_0, \Sigma_0$ as initial step, initial mean and initial covariance matrix for the distribution
2. Then take a sample $X_t$ from the sample distribution $q(X_0, \dots, X_{t-1})$ (usually it is taken $N(0, \Sigma)$ as symmetric simple equation), weighted with the ratio between the functions and the Q-ratio $\alpha(X_{t-1}, Y)$ as:

```math
X_t \backsim P^{SRWM}_q = \int_{A-X_{t-1}} \alpha(X_{t-1}, X_{t-1} + z) q(z \ | \ \mu, \Sigma) dz + \mathbb{I}(X_{t-1} \in A)\int_{X-X_{t-1}} (1 - \alpha(X_{t-1}, X_{t-1} + z)) q(z \ | \ \mu, \Sigma) dz
```

where I report a simplification of the real formula, which most of the physical distributions

4. Update the parameters according to the new value taken from the distribution and a non-decreasing non-converging series $\left[ \gamma_i \right]$, as follows:

$$
\begin{cases}
    \mu_{i+1} = \mu_i + \gamma_{i+1} (X_{i+1} - \mu_i) \\
    \Sigma_{i+1} = \Sigma_i + \gamma_{i+1} ((X_{i+1} - \mu_i) \cdot (X_{i+1} - \mu_i)^T - \Sigma_i)
\end{cases}
$$

4. Repeat points 3. and 4. until maximum number of samples

NB: This implementation is therefore useless, because of the fact that it is hard to perform the computations and obtain a shape for $P^{SRWM}_q$ and then one has to implement some approximations in order to retrieve a consistent chain. In my code I therefore show that this implementation (introducing a semplification) is useless.

#### Rao-Blackwellised algorithm (Rao)

This algorithm is a very slightly modification of the previous one: here, in order to solve the problem for the choice of the new value, it uses a normal sampling distribution $Y \backsim P^{Rao} = N(X_i, \Sigma_i)$ and the value is accepted with probability $\alpha(X_{t-1}, X_t)$, otherwise take $X_t = X_{t-1}$, as done in the Simple MCMC.

In this algorithm one has also the choice of the sequence $\gamma_{i+1}$, therefore I tried to implement some non-decreasing functions, such as an exponential decreasing or the reciprocal function rescaled: one can choose a series at will and optimize the evolution of the chain.

**Useful adaptation (Haario):** In many applications (and I tried it as well), the mean is corrected recursively and similarly the covariance matrix update $C_i = s_d \cdot \Sigma_i$ follows:

$$
C_t = 
\begin{cases}
    C_0 \hspace{41.5 mm} \text{for} \ t \leq t_0 \\
    s_d \left[\text{cov}(x_0, \dots, x_{t-1}) + \epsilon I_d\right] \hspace{5 mm} \text{for} \ t > t_0 
\end{cases}
$$

And it is found to be $s_d = \frac{2.38}{d^2}$ with $d$ dimension of the sampled distribution $[1]$ and $\epsilon$ is a small parameter, much smaller with respect to the typical dimensionality of the distribution.

#### Global adaptive scaling algorithm (GAM)

In the previous implementations, the covariance matrix was able to reshape itself, but it was not so able to change the variance in the many dimentions, leading to some problems in the increase of the acceptance rate. In the global-adaptive scaling algorithm there is the introduction of a global parameter that is able to scale the covariance matrix: $\lambda_i$, which is taken into account when evaluating $Y \backsim P^{GAM} = N(X_i, \lambda_i \cdot \Sigma_i)$

$$
\begin{cases}
    \mu_{i+1} = \mu_i + \gamma_{i+1} (X_{i+1} - \mu_i) \\
    \Sigma_{i+1} = \Sigma_i + \gamma_{i+1} ((X_{i+1} - \mu_i) \cdot (X_{i+1} - \mu_i)^T - \Sigma_i) \\
    \text{log}\left(\lambda_{i+1}\right) = \text{log}\left(\lambda_i\right) + \gamma_{i+1} \left(\alpha(X_i, X_{i+1}) - \bar{\alpha}\right)
\end{cases}
$$

where $\bar{\alpha}$ is the mean of the previous $\alpha \text{'s}$ that have been evaluated for the previous points.
This allows to:

- enlarge the variances if $\alpha(X_{t-1}, Y) > \bar{\alpha}$
- decrease them otherwise

## Results

In order to check the improvements of the introduced implementations, I run all the algorithms (available in "class.hmm") and made them run for some iterations. I only used algorithms from the introduction of the mtvnorm package and using the same number of "burn_in" and "total" steps. I divide then the performances based on the searched distribution and on the hyperparameters used, such as the number of "burn_in" and "total" steps.

I underline the fact that, using a gaussian sampling function, the acceptance rate wanted is the one approximating 68,3%, because of the fact that one should expect to switch to the new point according to $1\sigma$ of the gaussian distributioin.

#### 1. 3D Gaussian with good proportionate sample distribution (```cases/case_proportionate.ipynb```)

For a gaussian in three dimentions one sees that the fastest and most efficient way to retrieve a good result is using the "mvtnorm" function, when initializing with a proper correlation matrix. This is due to the fact that there is no need for the algorithms to improve the sampling distribution, because and the well-initialized distribution is rewarded: this applies always when initializing a distribution with a defined central shape and the right dimentions for the covariance matrix.

| Function      | Acceptance Rate (%) | Time needed (s) | Coda Autocorrelation | Distribution |
| ------------- | ------------------- | --------------- | -------------------- | ------------ |
| mvtnorm       | 49                  | 60              | zero                 | good         |
| mvtnorm Gibbs | 60                  | 190             | zero                 | good         |
| Haario        | 13                  | 40              | zero                 | cresty       |
| Haario Gibbs  | 16                  | 160             | zero                 | cresty       |
| Rao           | 50                  | 35              | first 10%/zero *     | peaked/good  |
| Rao Gibbs     | 60                  | 160             | zero                 | good         |
| Global        | 40                  | 65              | first 20%/zero *     | peaked/good  |
| Global Gibbs  | 55                  | 70              | zero                 | good         |

$* = $ the results strongly depends on the run (it can have a little bit of correlation, which vanishes when using )

#### 2. 3D Gaussian with small sample distribution (```cases/case_small.ipynb```)

In order to see the difficulties of these algorithms, I take a smaller covariance matrix: this is a problem for all the algorithms that don't allow a reshape of the matrix, because they al will result in a high autocorrelation.

| Function      | Acceptance Rate (%) | Time needed (s) | Coda Autocorrelation      | Distribution |
| ------------- | ------------------- | --------------- | ------------------------- | ------------ |
| mvtnorm       | 95                  | 50              | circa 10% constant *      | peaked       |
| mvtnorm Gibbs | <95                 | 160             | circa 10% constant *      | peaked       |
| Haario        | 17                  | 40              | zero                      | cresty       |
| Haario Gibbs  | 20                  | 65              | zero                      | cresty       |
| Rao           | 55                  | 35              | circa 10% at begginning * | good         |
| Rao Gibbs     | 65                  | 60              | zero                      | perfect      |
| Global        | 50                  | 60              | first 20%/zero *          | cresty       |
| Global Gibbs  | 65                  | 100             | zero                      | good         |

#### 3. 2D Gaussian times 1D Cauchy (```cases/case_2g_1c.ipynb```)

In order to stress the algorithms giving them different hard distributions, I decided to take a two-cauchy distribution (with two maxima) in order to see how it would have shaped the covariance distributions, in order to look for these two different maxima in the 3D space.

| Function      | Acceptance Rate (%) | Time needed (s) | Coda Autocorrelation               | Distribution           |
| ------------- | ------------------- | --------------- | ---------------------------------- | ---------------------- |
| mvtnorm       | 75                  | 65              | high at beginning in Cauchy *      | "too big" Cauchy tails |
| mvtnorm Gibbs | 85                  | 160             | high at very beginning in Cauchy * | "too big" Cauchy tails |
| Haario        | 17                  | 35              | zero                               | cresty                 |
| Haario Gibbs  | 35                  | 65              | zero                               | good                   |
| Rao           | 45                  | 30              | first 10%/good *                   | peaked                 |
| Rao Gibbs     | 60                  | 60              | zero                               | good                   |
| Global        | 55                  | 70              | first 20%/good *                   | peaked                 |
| Global Gibbs  | 70                  | 100             | zero                               | good                   |

#### 4. 1D Gaussian times 2D Cauchy (```cases/case_1g_2c.ipynb```)

To check if the algorithm was clever enough, I tried to enlarge the Cauchy distribution in another dimention, resulting in similar performances (even though I expected them to be worse)

| Function      | Acceptance Rate (%) | Time needed (s) | Coda Autocorrelation | Distribution              |
| ------------- | ------------------- | --------------- | -------------------- | ------------------------- |
| mvtnorm       | 80                  | 65              | circa 10% constant * | "too big" Cauchy tails    |
| mvtnorm Gibbs | 85                  | 160             | circa 10% constant * | "too big" Cauchy tails    |
| Haario        | 17                  | 35              | zero                 | cresty                    |
| Haario Gibbs  | 40                  | 70              | zero                 | good                      |
| Rao           | 40                  | 40              | zero                 | cresty                    |
| Rao Gibbs     | 50                  | 60              | zero                 | second peak not resoluted |
| Global        | 45                  | 70              | zero                 | cresty                    |
| Global Gibbs  | 75                  | 90              | zero                 | good                      |

#### 5. 5D Gaussian times 1D Cauchy (```cases/case_5g.ipynb```)

I also tried to increase the number of dimentions and the algorithm seemd to be pretty robust for central distributions (and a little less for those with multiple peaks)

| Function      | Acceptance Rate (%) | Time needed (s) | Coda Autocorrelation               | Distribution              |
| ------------- | ------------------- | --------------- | ---------------------------------- | ------------------------- |
| mvtnorm       | 70                  | 80              | circa 10% constant *               | "too big" Cauchy tails    |
| mvtnorm Gibbs | 85                  | 290             | high at very beginning in Cauchy * | "too big" Cauchy tails    |
| Haario        | 15                  | 50              | zero                               | cresty                    |
| Haario Gibbs  | 30                  | 120             | zero                               | cresty in some dimentions |
| Rao           | 35                  | 40              | zero                               | peaked                    |
| Rao Gibbs     | 50                  | 70              | zero                               | good                      |
| Global        | 50                  | 50              | circa 10% at beginning             | cresty                    |
| Global Gibbs  | 70                  | 90              | zero                               | good                      |

NB: I tried also other distributions, but the generated chains really are dependent on the choice of the initial point and covariance matrix, resulting in very unpleasant distributions

## Further improvements

There are other algorithms that I can experiment in the future, with a little adaptation of the code, which allow to increase the performances, at least theoretically speaking, like:

- Promote the global parameter $\lambda$ to vector in order to shape each gaussian per se, from using all the previous algorithms
- Make the algorithm change only in some localized dimentions per iteration, resulting in more robust results (as done in ML algorithms)
- Use principal components in order to get the best parameters for the covariance of the distribution

Further more, my code can be improved and written in a more efficient way, other than summarized into a unique function where to only choose the wanted algorithm.

## Bibliography

1. Haario, H., Saksman, E., Tamminen, J.: "An adaptive Metropolis algorithm"
2. Christophe Andrieu, Johannes Thoms: "A tutorial on adaptive MCMC"
3. Marko Laine: "Adaptive MCMC Methods with applications in environmental and geophysical models"
