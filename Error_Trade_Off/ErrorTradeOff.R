###########################################################
#
# The Error Dilemma with Finite Difference Methods
# for Numerical Differentiation
#
# Student Project on Numerical Differentiation
# Author: Clara Hoffmann
#
# based on the illustration of the error dilemma
# by Martins, Sturdza and Alonso, 2003
###########################################################

# set working directory
setwd("/Users/claracharlottehoffmann/Desktop/NumericalIntroductory")

# load packages
if (!require("pacman")) 
  install.packages("pacman"); library("pacman") 
p_load("Rcpp",
       "gganimate",
       "dplyr", 
       "reshape2",
       "ggplot2",
       "scales")
if (!require("transformr")) 
  devtools::install_github("thomasp85/transformr"); 
library("transformr") 

# function of which we want to find
# the derivative value
f <- function(v){
  y = (v^2)*sin(1/v)
  return(y)
}

# Two-Point Formula
TwoPointRep <- function(f, x, hstart, max){
  # Returns a dataframe containing the
  # approximation to the 
  # first-order derivative
  # using the two-point formula
  # 
  # Args:
  #   f:   target function for derivative
  #   x:   point at which we want to 
  #        estimate the derivative
  #   hstart: starting value of the step-size
  #   max: defines the range of
  #        the step-size, which starts
  #        with hstart^1 and ends with
  #        hstart^max. In every iteration
  #        1 is added to the power
  #
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

# the ThreePointRep, FivePointRep and
# ComplexStepRep work in the same way
# as the TwoPointRep
ThreePointRep <- function(f, x, hstart, max){
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
    deriv <- (func(x+h)-func(x-h))/(2*h)
    mat[count,2] <- deriv
    if (count == max){
      break
    }
  }
  return(mat)
}

FivePointRep <- function(f, x, hstart, max){
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
    # 5point numerical derivative
    deriv <-(func(x-2*h)-8*func(x-h)+
               8*func(x+h)-func(x+2*h))/(12*h)
    mat[count,2] <- deriv
    if (count == max){
      break
    }
  }
  return(mat)
}

# create complex number
i <- complex(real = 0, imaginary = 1)

ComplexStepRep <- function(f, x, hstart, max){
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
    # 3point numerical derivative
    deriv <-  Im(func(x+i*h))/h
    mat[count,2] <- deriv
    if (count == max){
      break
    }
  }
  return(mat)
}

# compute analytic derivative value
f.first <- function(v){
  y = 2*sin(1/v)*v-cos(1/v)
  return(y)
}
# true value of derivative
deriv.true <- f.first(0.1)

# approximate with 100 steps
mycomplex <- ComplexStepRep(f, x = 0.1, hstart = 0.5, max = 100) %>% as.data.frame()
mytwo <-TwoPointRep(f, x = 0.1, hstart = 0.5, max = 100) %>% as.data.frame()
mythree <- ThreePointRep(f, x = 0.1, hstart = 0.5, max = 100) %>% as.data.frame()
myfive <- FivePointRep(f, x = 0.1, hstart = 0.5, max = 100) %>% as.data.frame()

names(mycomplex) <- c("iteration", "Complex")
names(mytwo) <- c("iteration", "Two-Point")
names(mythree) <- c("iteration", "Three-Point")
names(myfive) <- c("iteration", "Five-Point")

# produce table for latex
output <- cbind(mytwo[,2], mythree[,2], 
                myfive[,2], mycomplex[,2])%>% as.data.frame  
names(output) <- c("two-point", 
                   "three-point", 
                   "five-point", 
                   "complex")
output.melt <- cbind(mytwo, mythree, myfive, mycomplex)  %>% 
  melt(id = "iteration")
names(output.melt) <- c("iteration", "Method", "value")
output.melt$iteration <- 0.1^output.melt$iteration

# plot relative approximation error over step size
error.melt <- output.melt
error.melt$value <- abs((output.melt$value - 
                           deriv.true)/ deriv.true)
error.plot <- ggplot( data = error.melt, 
                      aes(x = iteration, y = value,
                          group = Method, colour = Method,
                          shape = Method,
                          linetype = Method)) +
  geom_line() + 
  coord_cartesian(clip = 'off') +
  shadow_mark() +
  scale_x_log10(name = "h",
                labels=trans_format('log10',
                                    math_format(10^.x))) +
  scale_y_log10(name = "relative approximation error",
                labels=trans_format('log10',
                                    math_format(10^.x))) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +  
  # delete this option if you want a .pdf or .jpg
  #instead of a .gif :
  transition_manual(rev(iteration), cumulative = T) + 
  ease_aes('sine-in-out')
error.plot 
# save as .gif
anim_save("ErrorTradeOff.gif", animation = last_animation())
# or save as .pdf (enable all gganimate options before)
#ggsave("ErrorTradeOff.jpg", plot = error.plot  ) 
