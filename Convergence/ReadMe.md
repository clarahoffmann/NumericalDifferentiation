# Convergence of Numerical Differentiation Methods

# SOEP Explorer
There exist a wide range of numerical differentiation methods with different speed of convergence and accuracy. Finite difference formulas evaluate the function at different real points to approximate derivatives. However, as they converge to the true value of the derivative they fall victim to the round-off and subtractive cancellation error. These errors are due to the floating-point arithmetic of common computer systems, which fails to represent very small numbers with sufficient accuracy. An alternative to finite difference formulas is the complex step approach. As the name suggests, it evaluates functions at complex numbers. Fortunately this prevents the subtractive cancellation error.

This code illustrates the advantage of the complex step approach over finite difference formulas.
By approximating the first-order derivative of the function
f(x) = (x^2)sin(1/x)


<img src="SoepExplorerBasic.PNG" width="1200">

### Packages

```r
# load packages
if (!require("pacman")) 
  install.packages("pacman"); library("pacman") 

p_load("magrittr", 
       "dplyr", 
       "shiny",
       "plotly",
       "DT",
       "haven")

```

<img src="Convergence.jpg" width="400">
