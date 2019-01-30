# set working directory
setwd("/Users/claracharlottehoffmann/Desktop/NumericalIntroductory")

# load packages
if (!require("pacman")) 
  install.packages("pacman"); library("pacman") 
p_load("dplyr", 
       "reshape2",
       "ggplot2",
       "xtable",
       "scales")

# function of which we want to find
# the derivative value
f <- function(v){
  y = (v^2)*sin(1/v)
  return(y)
}

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
    deriv <-(func(x-2*h)-8*func(x-h)+8*func(x+h)-func(x+2*h))/(12*h)
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

mycomplex <- ComplexStepRep(f, x = 0.1, hstart = 0.01, max = 12) %>% as.data.frame()
mytwo <-TwoPointRep(f, x = 0.1, hstart = 0.01, max = 12) %>% as.data.frame()
mythree <- ThreePointRep(f, x = 0.1, hstart = 0.01, max = 12) %>% as.data.frame()
myfive <- FivePointRep(f, x = 0.1, hstart = 0.01, max = 12) %>% as.data.frame()

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
print(xtable(output, 
             type = "latex"), 
      file = "exampleoutput.tex")

# produce graph
output.melt <- cbind(mytwo, mythree, myfive, mycomplex)  %>% 
melt(id = "iteration")
names(output.melt) <- c("iteration", "Method", "value")
output.melt$iteration <- 0.1^output.melt$iteration

conver.plot <- ggplot( data = output.melt, 
                  aes(x = iteration, y = value,
                      group = Method, colour = Method,
                      shape = Method,
                      linetype = Method)) + 
    geom_line() + geom_point() + 
  scale_x_log10(name = "h",
                labels=trans_format('log10',math_format(10^.x))) + 
  ylim(c(-0.1, 1))  +
  ylab("approximated first-order derivative") + theme_bw()
conver.plot
ggsave("Convergence.jpg", plot = conver.plot ) # save