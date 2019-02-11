###########################################################
#
# Approximating Speed and Velocity with Finite Difference
# Formulas
#
# Student Project on Numerical Differentiation
# Author: Clara Hoffmann
#
# This code is the solution to problem 8.37 in
# Gilat and Subramaniam (2014)
###########################################################

# set working directory
setwd("...")

# load packages
if (!require("pacman")) 
  install.packages("pacman"); library("pacman") 
p_load("dplyr",
       "ggplot2",
       "gridExtra")

# create radar data
data <- cbind( t = c(0,3,6,9,12,15),
      r = c(18.803, 18.861, 18.946, 
            19.042, 19.148, 19.260),
      theta = c(0.7854, 0.7792, 0.7701, 
                0.7594, 0.7477, 0.7350)) %>% 
  as.data.frame()

FirstSecDeriv <- function(df){
  # Returns a dataframe containing the time
  # velocity (first derivative) and acceleration
  # (second derivative)
  # 
  # Args:
  #   df:   dataframe containing radar data 
  #         with the following variables:
  #         df$r : distance towards tracked object in km
  #         df$t : time of tracking in s
  #         df$theta : angle of radar in rad
  # Returns:
  #   dataframe of time, velocity in km/h and
  #   acceleration in km/h/s
  #
  l <- length(df$r)
  # first-order derivatives
  # compute derivative with three-point formula 
  # f(x+h)-f(x-h)/h
  deriv.r <- (lead(df$r)-lag(df$r))/(lead(df$t)-lag(df$t))
  # correct first and last value with two-point formula
  # f(x+h)-f(x)/h
  deriv.r[1] <-  (df$r[2] - df$r[1])/(df$t[2] - df$t[1])
  deriv.r[l] <- (df$r[l-1] - df$r[l])/(df$t[l-1] - df$t[l])
  # repeat the same for the derivative of theta
  deriv.theta <- (lead(df$theta)-lag(df$theta))/
    (lead(df$t)-lag(df$t))
  deriv.theta[1] <-  (df$theta[2] - df$theta[1])/
    (df$t[2] - df$t[1])
  deriv.theta[l] <- (df$theta[l-1] - df$theta[l])/
    (df$t[l-1] - df$t[l])
  # compute velocity
  fd <- (sqrt((deriv.r^2) + 
                (df$r*deriv.theta)^2))*3600
  # second-order derivatives
  # use two- and three-point formula as before
  deriv.r.sec <- (lead(deriv.r) - lag(deriv.r))/
    (lead(df$t)-lag(df$t))
  # replace last and first value
  deriv.r.sec[1] <-  (deriv.r[2] - deriv.r[1])/
    (df$t[2] - df$t[1])
  deriv.r.sec[l] <- (deriv.r[l-1] -  deriv.r[l])/
    (df$t[l-1] - df$t[l])
  # repeat for theta
  deriv.theta.sec <- (lead(deriv.theta) - 
                        lag(deriv.theta))/
    (lead(df$t)-lag(df$t))
  deriv.theta.sec[1] <-  (deriv.theta[2] - 
                            deriv.theta[1])/
    (df$t[2] - df$t[1])
  deriv.theta.sec[l] <- (deriv.theta[l-1] - 
                           deriv.theta[l])/
    (df$t[l-1] - df$t[l])
  # compute acceleration
  sd <-(sqrt(
    (deriv.r.sec + df$r * (deriv.theta^2))^2 +
    (df$r*deriv.theta.sec + 
       2*deriv.r*deriv.theta)^2))*3600
  # bind result as dataframe
  result <- cbind(t = df$t, 
                  first.derivative = fd,
                  second.derivative = sd) %>% 
    as.data.frame
  return(result)
}

# compute derivatives for given data
result <-  FirstSecDeriv(data)
result 

# plot the results
# for velocity
velo.plot <- ggplot(data = result , 
                   aes(x = t, 
                       y = first.derivative)) + 
  geom_point() + 
  geom_line() + theme_bw() +
  xlab("time in seconds") +
  ylab("approx. instantaneous velocity in km/h")  +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) 
velo.plot
# for acceleration
acc.plot <- ggplot(data = result , 
                  aes(x = t, y = second.derivative)) + 
  geom_point() + 
  geom_line() + theme_bw() +
  xlab("time in seconds") +
  ylab("approx. instantaneous acceleration in km/h/s") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

acc.plot
plot <- grid.arrange(velo.plot, acc.plot, nrow = 1)
ggsave("RadarSpeed.jpg", plot = plot ) # save
