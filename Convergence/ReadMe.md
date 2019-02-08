# <p align='center'>Convergence of Numerical Differentiation Methods</p>

# Content
* Common Methods
* Code Structure
* Sources

# Common Methods
There exist a wide range of numerical differentiation methods with different speed of convergence and accuracy. Finite difference formulas evaluate the function at different real points to approximate derivatives. However, as they converge to the true value of the derivative they fall victim to the round-off and subtractive cancellation error. These errors are due to the floating-point arithmetic of common computer systems, which fails to represent very small numbers with sufficient accuracy. An alternative to finite difference formulas is the complex step approach. As the name suggests, it evaluates functions at complex numbers. Fortunately, this prevents the subtractive cancellation error.

This code illustrates the advantage of the complex step approach over finite difference formulas.
By approximating the first-order derivative of the function

<a href="https://www.codecogs.com/eqnedit.php?latex=f(x)&space;=&space;(x^2)sin(1/x)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?f(x)&space;=&space;(x^2)sin(1/x)" title="f(x) = (x^2)sin(1/x)" /></a>

at the point x = 0.1 with the

* two-point formula \
<a href="https://www.codecogs.com/eqnedit.php?latex=f'(x)&space;\approx&space;\frac{f(x&plus;h)-f(x)}{h}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?f'(x)&space;\approx&space;\frac{f(x&plus;h)-f(x)}{h}" title="f'(x) \approx \frac{f(x+h)-f(x)}{h}" /></a>
* three-point formula \
<a href="https://www.codecogs.com/eqnedit.php?latex=f'(x)&space;\approx&space;\frac{f(x&plus;h)-f(x-h)}{2h}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?f'(x)&space;\approx&space;\frac{f(x&plus;h)-f(x-h)}{2h}" title="f'(x) \approx \frac{f(x+h)-f(x-h)}{2h}" /></a>
* five-point formula \
<a href="https://www.codecogs.com/eqnedit.php?latex=f'(x)&space;=&space;\frac{-f(x&plus;2h)&space;&plus;&space;8f(x&plus;h)&space;-8f(x-h)&space;&plus;f(x-2h)&space;}{12h}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?f'(x)&space;=&space;\frac{-f(x&plus;2h)&space;&plus;&space;8f(x&plus;h)&space;-8f(x-h)&space;&plus;f(x-2h)&space;}{12h}" title="f'(x) = \frac{-f(x+2h) + 8f(x+h) -8f(x-h) +f(x-2h) }{12h}" /></a>
* complex step approach \
<a href="https://www.codecogs.com/eqnedit.php?latex=f'(x)&space;\approx&space;\frac{Im(f(x&plus;ih)}{h}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?f'(x)&space;\approx&space;\frac{Im(f(x&plus;ih)}{h}" title="f'(x) \approx \frac{Im(f(x+ih)}{h}" /></a>,

where h denotes the step size.
The different convergence rates and the collapse of the finite different formulas at very small step sizes h are illustrated in a graph.

<img src="Convergence.jpg" width="600">

The finite difference methods collapse to zero with very small absolute values of the step size h due to the subtractive cancellation error. The complex step approach (purple crosses) meanwhile keeps on converging to the true value of the derivative.
### Code Structure

Define the function of which we want to find the derivative.
```r
f <- function(v){
  y = (v^2)*sin(1/v)
  return(y)
}
```

Define the two-point formula function. It takes the target function, the x-value at which the derivative should be approximated, a starting value for the step size, and a maximum number of iterations as inputs.
It will approximate the derivative at y with h^1, h^2, ... h^max and return a vector of the approximation at these step sizes.
```r
# Two-Point Formula
TwoPointRep <- function(f, x, hstart, max){
  # define function
  func <- f
  # empty matrix for returning values
  v <- 1:max
  mat <- matrix(rep(v,2), ncol = 2)
  # set counter to zero
  count = 0
  # start producing the derivatives
  repeat {
    count = count + 1
    h = hstart^count
    # 2point numerical derivative
    deriv <- (func(x+h)-func(x))/h
    mat[count,2] <- deriv
    if (count == max){
      break
    }
  }
  return(mat)
}
```
The three-point and five-point formula follow the same structure.

For the complex step approach, we first have to create a complex number...
```r
# create complex number
i <- complex(real = 0, imaginary = 1)
```
... and retrieve the approximation for the derivative by only evaluating the imaginary part of the complex number.
```r
deriv <-  Im(func(x+i*h))/h
```
The skeleton for producing a sequence of derivative approximation stays the same.
Call the functions to produce the numerical derivatives ...
```r
Call the functions to produce the numerical derivatives ...
mycomplex <- ComplexStepRep(f, x = 0.1, hstart = 0.01, max = 12) %>% as.data.frame()
mytwo <-TwoPointRep(f, x = 0.1, hstart = 0.01, max = 12) %>% as.data.frame()
mythree <- ThreePointRep(f, x = 0.1, hstart = 0.01, max = 12) %>% as.data.frame()
myfive <- FivePointRep(f, x = 0.1, hstart = 0.01, max = 12) %>% as.data.frame()
```

In the end we bind the data as a data.frame and produce a latex table with xtable() as well as the graph seen above with ggplot().

# Sources
For the complex step approach:
 * Joaquim R. R. A. Martins, Peter Sturdza and Juan J. Alonso. The Complex-Step Derivative Approximation. ACM Transactions on Mathematical Software, Vol.29, No. 3, September 2003, Pages 245–262, 2003.
 * William Squire, George Trapp. Using Complex Variables to Estimate Derivatives of Real Functions. SIAM Rev. Vol. 40, No. 1 p. 110-112, 1998. \\

For finite difference formulas and the structure of the floating-point arithmetic:
* C. Woodford, C. Phillips. Numerical Methods with Worked Examples: MatlabEdition. Springer, p. 119-128, 2012.
* James E. Gentle. Numerical Linear Algebra for Applications in Statistics. George Mason University, p. 1 - 12, 1998


