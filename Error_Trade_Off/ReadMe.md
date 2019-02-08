# The Dilemma of Choosing a Step-Size in Finite Difference Methods

Finite difference formulas in numerical differentiation theoretically deliver more accurate results the smaller their step size becomes. The approximation becomes more exact the closer the points used for the approximation are to the function value at which we want to obtain the derivative. This implies that the absolute step size should converge to zero.
However, finite difference methods use the difference between close function values to obtain the derivative. Common computer systems operate on a 64-bit floating-point arithmetic. This arithmetic aims to represent as wide a range of numbers with as much accuracy as possible. In particular, there is a minimum distance between each number that this system can represent.
Especially with very small absolute numbers near zero, computers cannot represent numbers with sufficient accuracy.
This leads to round-off errors and the subtractive cancellation error.
Take for example the two-point formula:

<a href="https://www.codecogs.com/eqnedit.php?latex=f'(x)&space;\approx&space;\frac{f(x&plus;h)-f(x-h)}{2h}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?f'(x)&space;\approx&space;\frac{f(x&plus;h)-f(x-h)}{2h}" title="f'(x) \approx \frac{f(x+h)-f(x-h)}{2h}" /></a>.

As the step size h approaches zero, the computer will fail to represent the distance between the function values in the enumerator. This is because it becomes so close to zero. At first, the computer will round the distance to the next representable number. Because we divide by another small number this can produce rounding errors. Ultimately, the enumerator will become so small that the computer rounds it to zero - this makes the approximation useless!

This code produces a graph that illustrates the trade-off between large and small absolute step sizes in finite difference formulas. As an alternative, the complex step approach is used. The latter is unaffected by the subtractiv cancellation error and converges to the true derivative value long after the finite difference methods have deteriorated. For this reason the complex step approach is the standard in current numerical differentiation software packages (see for example NumDeriv in R).

<img src="ErrorTradeOff.jpg" width="600">

The dilemma was also illustrated in this way by Martins, Sturdza and Alonsa (2003). 
In this example the derivative of the function

<a href="https://www.codecogs.com/eqnedit.php?latex=f(x)&space;=&space;(x^2)sin(1/x)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?f(x)&space;=&space;(x^2)sin(1/x)" title="f(x) = (x^2)sin(1/x)" /></a>

is approximated at the point x = 0.1.

# Code Structure
Define the function of which we want to estimate the derivative.

```r
f <- function(v){
  y = (v^2)*sin(1/v)
  return(y)
}
```

# Sources
For the complex step approach:
 * Joaquim R. R. A. Martins, Peter Sturdza and Juan J. Alonso. The Complex-Step Derivative Approximation. ACM Transactions on Mathematical Software, Vol.29, No. 3, September 2003, Pages 245–262, 2003.
 * William Squire, George Trapp. Using Complex Variables to Estimate Derivatives of Real Functions. SIAM Rev. Vol. 40, No. 1 p. 110-112, 1998. \\

For finite difference formulas and the structure of the floating-point arithmetic:
* C. Woodford, C. Phillips. Numerical Methods with Worked Examples: MatlabEdition. Springer, p. 119-128, 2012.
* James E. Gentle. Numerical Linear Algebra for Applications in Statistics. George Mason University, p. 1 - 12, 1998
