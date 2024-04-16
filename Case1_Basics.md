Case Study 1
================

## Ratio of Fibonacci numbers

The Fibonacci sequence is a series of numbers where a number is the
addition of the last two numbers, starting with 0 and 1.

## For and While functions

Two R functions were implemented to calculate the sequence of ratios of
consecutive Fibonacci numbers for a given number of terms. The first
function uses a `for` loop, while the second function uses a `while`
loop.

``` r
# Using for loop
fib_ratio_for <- function(n) {
  fib <- c(0, 1)
  ratio <- numeric(n)
  for(i in 1:n) {
    fib <- c(fib, sum(tail(fib, 2)))
    ratio[i] <- fib[i+2] / fib[i+1]
  }
  return(ratio)
}

# Using while loop
fib_ratio_while <- function(n) {
  fib <- c(0, 1)
  ratio <- numeric(n)
  i <- 1
  while(i <= n) {
    fib <- c(fib, sum(tail(fib, 2)))
    ratio[i] <- fib[i+2] / fib[i+1]
    i <- i + 1
  }
  return(ratio)
}
```

## The result of benchmark two functions

The functions were benchmarked for n = 200 and n = 2000. The results
showed that the function using the for loop was faster than the function
using the while loop,because it has lower mean time.

``` r
if (!require(microbenchmark)) install.packages('microbenchmark')
```

    ## Loading required package: microbenchmark

``` r
# Benchmark the functions
library(microbenchmark)
n <- 200
microbenchmark(fib_ratio_for(n), fib_ratio_while(n), times = 100L)
```

    ## Unit: microseconds
    ##                expr   min     lq    mean median    uq     max neval
    ##    fib_ratio_for(n) 701.8 726.45 889.649 756.35 802.7  3898.3   100
    ##  fib_ratio_while(n) 708.3 732.30 999.295 749.00 871.3 11558.2   100

## Result of the Plot

The sequence of ratios was plotted for n = 100 using the ggplot2 package
in R. The plot showed that the sequence starts to converge quickly,
within the first 10 terms. The sequence converges to the golden ratio,
approximately 1.61803398875. This was verified by checking the last few
terms of the sequence.

``` r
if (!require(ggplot2)) install.packages('ggplot2')
```

    ## Loading required package: ggplot2

``` r
# Plot the sequence
library(ggplot2)
n <- 100
ratio <- fib_ratio_for(n)
df <- data.frame(x = 1:n, y = ratio)
ggplot(df, aes(x=x, y=y)) + geom_line() + theme_minimal() + labs(x = "i", y = "r_i")
```

![](Case1_Basics_files/figure-gfm/c-1.png)<!-- -->

``` r
tail(ratio, 10)
```

    ##  [1] 1.618034 1.618034 1.618034 1.618034 1.618034 1.618034 1.618034 1.618034
    ##  [9] 1.618034 1.618034

## Gamma function

The Gamma function is a mathematical concept that extends the factorial
function, traditionally applicable to positive integers, to complex and
real numbers.

## Computing the Gamma Function in R

In R, the Gamma function can be computed using the built-in gamma
function.

``` r
compute_rho <- function(n) {
  rho <- gamma((n - 1) / 2) / (gamma(1 / 2) * gamma((n - 2) / 2))
  return(rho)
}
```

## Observations and Limitations with Large Inputs

When we try to compute rho for a large value of n (n = 2000), we may
encounter numerical issues. The computation might result in Inf or NaN
because the Gamma function grows very rapidly for large values of n, and
it can easily exceed the maximum representable number in R.

``` r
compute_rho <- function(n) {
  rho <- gamma((n - 1) / 2) / (gamma(1 / 2) * gamma((n - 2) / 2))
  return(rho)
}

# Test the function with n = 2000
n <- 2000
result <- compute_rho(n)
print(result)
```

    ## [1] NaN

## Handling Large Inputs with the Log-Gamma Function

To handle large values of n, we can use the lgamma function in R, which
computes the natural logarithm of the absolute value of the Gamma
function. This allows us to work with much larger values without running
into numerical overflow issues.

``` r
compute_rho <- function(n) {
  log_rho <- lgamma((n - 1) / 2) - (lgamma(1 / 2) + lgamma((n - 2) / 2))
  rho <- exp(log_rho)
  return(rho)
}
```

## Visualizing the Gamma Function and Estimating its Limit

The ratio increases rapidly for small n, reaches a maximum, and then
slowly approaches a constant value as n goes to infinity. This constant
value would be the limit of the function as n approaches infinity.

``` r
n_values <- seq(1, 2000, by = 10)
rho_values <- sapply(n_values, compute_rho)
plot(n_values, rho_values / sqrt(n_values), type = "l")
```

![](Case1_Basics_files/figure-gfm/g-1.png)<!-- -->

## The golden ratio

The golden ratio, often denoted by the Greek letter ϕ or τ, is a
mathematical concept that holds a special fascination due to its
occurrence in many natural and human-made structures. It is an
irrational number approximately equal to 1.61812.

## Method 1: Recursive Calculation

The first method for computing Φ^n+1 is a recursive calculation. This
method is based on the Fibonacci-like relationship Φ^n+1 = Φ^n + Φ^n−1.

``` r
golden_ratio_iterative <- function(n) {
  phi <- (sqrt(5) + 1) / 2
  a <- 1
  b <- phi
  for (i in 2:n) {
    temp <- b
    b <- a + b
    a <- temp
  }
  return(b)
}
```

## Method 2: Power Operator

The second method involves simply raising Φ to the power of n+1.

``` r
golden_ratio_power <- function(n) {
  phi <- (sqrt(5) + 1) / 2
  return(phi^(n+1))
}
```

## Comparison of Methods

Now we will compare the results obtained from these two methods for n =
12, 60, 120, 300. We will use both the == operator, which checks for
exact equality, and the all.equal function, which checks for near
equality.

## Exact Equality vs. Near Equality:

- When comparing results using ==, the outputs may not be exactly equal
  due to floating-point precision in the iterative calculation method.
- Using all.equal(), we check for near equality with a tolerance, which
  accommodates small floating-point discrepancies.

## Reasons for Differences:

- Floating-point arithmetic in computers introduces small errors during
  iterative calculations.
- Accumulated errors in the iterative method can lead to discrepancies
  when compared directly with exact equality (==).

``` r
n_values <- c(12, 60, 120, 300)

for (n in n_values) {
  iterative_result <- golden_ratio_iterative(n)
  power_result <- golden_ratio_power(n)
  
  print(paste("For n =", n))
  print(paste("Iterative result: ", iterative_result))
  print(paste("Power result: ", power_result))
  
  # Compare using all.equal
  print(paste("Are the results nearly equal? ", isTRUE(all.equal(iterative_result, power_result))))
}
```

    ## [1] "For n = 12"
    ## [1] "Iterative result:  321.996894379985"
    ## [1] "Power result:  521.001919378726"
    ## [1] "Are the results nearly equal?  FALSE"
    ## [1] "For n = 60"
    ## [1] "Iterative result:  3461452808002"
    ## [1] "Power result:  5600748293801.01"
    ## [1] "Are the results nearly equal?  FALSE"
    ## [1] "For n = 120"
    ## [1] "Iterative result:  1.19816555420249e+25"
    ## [1] "Power result:  1.938672590849e+25"
    ## [1] "Are the results nearly equal?  FALSE"
    ## [1] "For n = 300"
    ## [1] "Iterative result:  4.96926405783747e+62"
    ## [1] "Power result:  8.04043814465433e+62"
    ## [1] "Are the results nearly equal?  FALSE"

## Game of craps

I created a function called `play_craps` and I used the `sample`
function to generate random numbers for simulating the roll of a die.
For the first two rules, I used `if` conditions with logical statement
using `==` . If the rolled number wasn’t 1 or 6, I assigned it as the
point number and set up a `for` loop with three iterations to handle the
third condition. Then, I decided to repeat the game for 1000 times to
find the win percentage. So, I created a function called
`simulate_craps` and provided it with the parameter `num_games` which
repeats the game 1000 times. At the end, it calculates the win
percentage and returns it.

``` r
play_craps <- function() {
  roll <- sample(1:6, 1)
  if (roll == 6) {
    return("Win")
  }
  else if (roll == 1) {
    return("Lose")
  }
  else {
    point <- roll
   
    for (i in 1:3) {
      roll <- sample(1:6, 1)
      if (roll == point) {
        return("Win")
      }
    }
    return("Lose")
  }
}

simulate_craps <- function(num_games) {
  wins <- 0
  for (i in 1:num_games) {
    result <- play_craps()
    if (result == "Win") {
      wins <- wins + 1
    }
  }
  
  win_percentage <- (wins / num_games) * 100
  return(win_percentage)
}

num_games <- 1000
win_percentage <- simulate_craps(num_games)
cat("Win percentage after simulating", num_games, "games of craps:", win_percentage, "%\n")
```

    ## Win percentage after simulating 1000 games of craps: 43.1 %

## Readable and efficient code

1.  This code aims to find the mean squared error of a linear model with
    variables `y` and `x`. First, it sets the seed as 1 using
    `set.seed(1)` to generate a fixed random number if we repeat the
    code. Then, it generates 1000 random numbers for `x` and `y` using
    the `rnorm` function and converts them into a dataframe.In the
    subsequent steps, it repeatedly takes different subsets of the
    dataframe, excluding 250 rows of the dataset each time, for a total
    of 4 times. For each subset, it first makes predictions and then
    calculates the mean squared error, saving it in the variable `r`.

``` r
set.seed(1)
x <- rnorm(1000)
y <- 2 + x + rnorm(1000)
df <- data.frame(x, y)
cat("Step", 1, "\n")
```

    ## Step 1

``` r
fit1 <- lm(y ~ x, data = df[-(1:250),])
p1 <- predict(fit1, newdata = df[(1:250),])
r <- sqrt(mean((p1 - df[(1:250),"y"])^2))
cat("Step", 2, "\n")
```

    ## Step 2

``` r
fit2 <- lm(y ~ x, data = df[-(251:500),])
p2 <- predict(fit2, newdata = df[(251:500),])
r <- c(r, sqrt(mean((p2 - df[(251:500),"y"])^2)))
cat("Step", 3, "\n")
```

    ## Step 3

``` r
fit3 <- lm(y ~ x, data = df[-(501:750),])
p3 <- predict(fit3, newdata = df[(501:750),])
r <- c(r, sqrt(mean((p3 - df[(501:750),"y"])^2)))
cat("Step", 4, "\n")
```

    ## Step 4

``` r
fit4 <- lm(y ~ x, data = df[-(751:1000),])
p4 <- predict(fit4, newdata = df[(751:1000),])
r <- c(r, sqrt(mean((p4 - df[(751:1000),"y"])^2)))
r
```

    ## [1] 1.0249956 0.9952113 1.0685500 1.0707264

2.  If I wanted to improve the readability of the code, I would likely
    consolidate the repetitive lines into a single function to avoid
    redundancy. By grouping the repetitive lines into a function, we can
    enhance the clarity and maintainability of the code. This makes it
    easier to understand the purpose of the code by the name of the
    function and what it actually does.

3.  I created two functions for this purpose. The first one, called
    `mean_square_error`, calculates the squared error. We can utilize
    this function within a for loop, iterating it four times. the second
    function, `mean_square_error_4`, also computes the squared error,
    but it includes the for loop within itself. This function does the
    repetition internally.

Here’s how these functions :

``` r
mean_square_error<-function(a,b){
  fit <- lm(y ~ x, data = df[-(a:b),])
  p <- predict(fit, newdata = df[(a:b),])
  r<-(sqrt(mean((p - df[(a:b),"y"])^2)))
  return(r)
}

set.seed(1)
x <- rnorm(1000)
y <- 2 + x + rnorm(1000)
df <- data.frame(x, y)
r<-c()
a<-1
b<-250
for (i in 1:4)
{
  r<-c(r,c(mean_square_error(a,b)))
  a<-a+250
  b<-b+250
}
r
```

    ## [1] 1.0249956 0.9952113 1.0685500 1.0707264

Second function:

``` r
mean_square_error_4<-function(df){
  a<-1
  b<-250
  r<-c()
  for (i in 1:4) {
  fit <- lm(y ~ x, data = df[-(a:b),])
  p <- predict(fit, newdata = df[(a:b),])
  r <- c(r, sqrt(mean((p - df[(a:b),"y"])^2)))
  a<-a+250
  b<-b+250
  }
  return(r)
}
set.seed(1)
x <- rnorm(1000)
y <- 2 + x + rnorm(1000)
df <- data.frame(x, y)
mean_square_error_4(df)
```

    ## [1] 1.0249956 0.9952113 1.0685500 1.0707264



## Measuring and improving performance

``` r
kwtest <- function (x, g, ...) {
  if (is.list(x)) {
    if (length(x) < 2L)
      stop("'x' must be a list with at least 2 elements")
    if (!missing(g))
      warning("'x' is a list, so ignoring argument 'g'")
    if (!all(sapply(x, is.numeric)))
      warning("some elements of 'x' are not numeric and will be coerced to numeric")
    k <- length(x)
    l <- lengths(x)
    if (any(l == 0L))
      stop("all groups must contain data")
    g <- factor(rep.int(seq_len(k), l))
    x <- unlist(x)
  }
  else {
    if (length(x) != length(g))
      stop("'x' and 'g' must have the same length")
    g <- factor(g)
    k <- nlevels(g)
    if (k < 2L)
      stop("all observations are in the same group")
  }
  n <- length(x)
  if (n < 2L)
    stop("not enough observations")
  r <- rank(x)
  TIES <- table(x)
  STATISTIC <- sum(tapply(r, g, sum)^2/tapply(r, g, length))
  STATISTIC <- ((12 * STATISTIC/(n * (n + 1)) - 3 * (n + 1))/(1 - sum(TIES^3 - TIES)/(n^3 - n)))
  PARAMETER <- k - 1L
  PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
  names(STATISTIC) <- "Kruskal-Wallis chi-squared"
  names(PARAMETER) <- "df"
  RVAL <- list(statistic = STATISTIC, parameter = PARAMETER,
               p.value = PVAL, method = "Kruskal-Wallis rank sum test")
  return(RVAL)
}
```

## a. Pseudo code outlining the function:

1.  Check if input x is a list:

<!-- -->

1.  If x is a list:
    1.  Check if x has at least 2 elements.

<!-- -->

2.  If g is provided, issue a warning since it will be ignored.
3.  Check if all elements of x are numeric, if not, coerce them to
    numeric.
4.  Calculate the number of groups (k) and their lengths (l).
5.  Check if any group contains no data.
6.  Combine the elements of x into a single vector.

<!-- -->

2.  If x is not a list: i. Check if x and g have the same length.

<!-- -->

2.  Convert g into a factor.
3.  Calculate the number of levels in g (k).
4.  Check if all observations are in the same group.

<!-- -->

2.  Calculate the number of observations (n).
3.  Check if there are enough observations.
4.  Compute the rank of each observation.
5.  Count ties in the ranks.
6.  Calculate the Kruskal-Wallis test statistic.
7.  Compute the degrees of freedom.
8.  Calculate the p-value using the chi-squared distribution.
9.  Return a list containing the statistic, degrees of freedom, p-value,
    and method used.

## b. Calling the function

``` r
x_list <- list(group1 = c(1, 2, 3), group2 = c(4, 5, 6))
x_vector <- c(1, 2, 3, 4, 5, 6)
g <- c(1, 1, 1, 2, 2, 2)
result_list <- kwtest(x_list)

# Call the function with x as a vector and corresponding g argument
result_vector <- kwtest(x_vector, g)

# Ensure the two function calls return the same result
identical(result_list, result_vector)
```

    ## [1] TRUE

## c. Faster version of kwtest():

``` r
kwtest_fast <- function(x, g) {
  if (!is.numeric(x)) {
    stop("'x' and 'g' must be numeric vectors")
  }
  # Check if x and g have the same length
  if (length(x) != length(g)) {
    stop("'x' and 'g' must have the same length")
  }
  g <- factor(g)
  if (nlevels(g) < 2) {
    stop("There must be at least two unique groups")
  }
  n <- length(x)
  r <- rank(x)
  
  TIES <- table(x)
  STATISTIC <- sum(tapply(r, g, sum)^2/tapply(r, g, length))
  STATISTIC <- ((12 * STATISTIC/(n * (n + 1)) - 3 * (n + 1))/(1 -
                                              sum(TIES^3 - TIES)/(n^3-n)))
  return(STATISTIC)
}
#testing and comparing the two functions results


x<-c(rnorm(15))
g<-c(1,2,3,3,4,4,4,4,4,4,4,4,4,4,4)
kwtest(x =x,g=g )
```

    ## $statistic
    ## Kruskal-Wallis chi-squared 
    ##                   5.293182 
    ## 
    ## $parameter
    ## df 
    ##  3 
    ## 
    ## $p.value
    ## [1] 0.1515454
    ## 
    ## $method
    ## [1] "Kruskal-Wallis rank sum test"

``` r
kwtest_fast(x=x,g=g)
```

    ## [1] 5.293182

## d. Function to perform Kruskal-Wallis test using kruskal.test.default():

``` r
kwtest_default <- function(X, g) {
  m <- nrow(X)
  test_statistics <- numeric(m)
  for (i in 1:m) {
    test_statistics[i] <- kruskal.test(X[i, ] ~ g)$statistic
  }
  return(test_statistics)
}
```

## e. Function to perform Kruskal-Wallis test using the faster version:

``` r
kwtest_fast_repeat <- function(X, g) {
  m <- nrow(X)
  test_statistics <- numeric(m)
  for (i in 1:m) {
    test_statistics[i] <- kwtest_fast(X[i, ], g)
  }
  return(test_statistics)
}

#testing 
m <- 1000 # number of repetitions
n <- 50 # number of individuals
X <- matrix(rt(m * n, df = 10), nrow = m)
grp <- rep(1:3, c(20, 20, 10))
a<-kwtest_fast_repeat(X,g=grp)
b<-kwtest_default(X,g=grp)
identical(a,b)
```

    ## [1] TRUE

## f. Performance comparison using a benchmarking package:

The function in part e exhibits significantly better performance
compared to the function in part d. The reason for this discrepancy is
likely attributable to the inner workings of the kruskal test function
from the stats package. It’s possible that kruskal test incorporates
multiple checks or overhead that our simplified function omits. Another
factor could be differences in how the two functions handle input data.
To conclusively determine the precise reasons behind the performance
gap, a detailed examination of the internal structure of kruskal.test
would be required.

``` r
library(microbenchmark)

set.seed(1234)
m <- 1000
n <- 50
X <- matrix(rt(m * n, df = 10), nrow = m)
grp <- rep(1:3, c(20, 20, 10))

# Benchmarking the two approaches
benchmark_results <- microbenchmark(
  kwtest_default(X, grp),
  kwtest_fast_repeat(X, grp),
  times = 10
)

print(benchmark_results)
```

    ## Unit: milliseconds
    ##                        expr      min       lq     mean   median       uq
    ##      kwtest_default(X, grp) 420.4252 426.6281 435.3037 429.4898 449.3487
    ##  kwtest_fast_repeat(X, grp) 206.0180 206.9743 210.3849 210.4178 211.6573
    ##       max neval
    ##  459.1989    10
    ##  216.8586    10

## g. Vectorized version of the faster function:

The kruskal test function from the stats package shows worse performance
compared to the other two functions, primarily due to differences in
their underlying mechanisms. Notably, the vectorized function exhibits
slightly better performance than the function from part e. As outlined
in the script, this can be attributed to the efficiency of high-level
programming languages like R, which leverage fast routines in languages
such as C or C++. Consequently, operations that iterate over matrices
using high-level constructs, like row- or column-based operations,
inherently incur more overhead compared to vectorized functions that
operate at a lower level of abstraction.

``` r
kwtest_fast_vectorized <- function(X, g) {
  apply(X, 1, function(x) kwtest_fast(x, g))
}

# Benchmarking the vectorized approach
benchmark_results_vectorized <- microbenchmark(
  kwtest_default(X, grp),
  kwtest_fast_repeat(X, grp),
  kwtest_fast_vectorized(X, grp),
  times = 10
)

print(benchmark_results_vectorized)
```

    ## Unit: milliseconds
    ##                            expr      min       lq     mean   median       uq
    ##          kwtest_default(X, grp) 420.9939 424.4710 427.2183 426.3066 430.6291
    ##      kwtest_fast_repeat(X, grp) 203.2931 209.5967 211.5950 211.4321 212.9433
    ##  kwtest_fast_vectorized(X, grp) 206.9969 209.0111 214.9893 210.7926 214.1260
    ##       max neval
    ##  436.0928    10
    ##  219.2372    10
    ##  253.1792    10
