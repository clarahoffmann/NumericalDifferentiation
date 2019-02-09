# <p align='center'>Measuring Velocity and Acceleration with Numerical Differentiation</p>
Solution to problem 8.37 in Gilat \& Subramaniam (2014)

# Content
* Using Radars to Measure Speed & Velocity
* Code
* Sources

# Using Radars to Measure Speed & Velocity
Many radars that measures speed of a moving object use numerical differentiation. 
In this appplication, a radar station which measures the speed of passing planes uses numerical differentiation to obtain their speed and acceleration. This example is taken from Gilat \& Subramaniam (2014). Assume the angle and distance towards an aircraft is recorded every three seconds. The velocity and acceleration can be determined using the distance to the airplane (r), the angle of the radar to the ground (theta) and the time (t) by 

<a href="https://www.codecogs.com/eqnedit.php?latex=v&space;=&space;\sqrt{\bigg(\frac{dr}{dt}\bigg)^2&space;&plus;&space;\bigg(r\frac{d\theta}{dt}\bigg)^2}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?v&space;=&space;\sqrt{\bigg(\frac{dr}{dt}\bigg)^2&space;&plus;&space;\bigg(r\frac{d\theta}{dt}\bigg)^2}" title="v = \sqrt{\bigg(\frac{dr}{dt}\bigg)^2 + \bigg(r\frac{d\theta}{dt}\bigg)^2}" /></a> \
 and \
<a href="https://www.codecogs.com/eqnedit.php?latex=a&space;=&space;\sqrt{\Bigg[\frac{d^2r}{dt^2}&space;&plus;&space;r\bigg(\frac{d\theta}{dt}\bigg)^2\Bigg]^2&space;&plus;&space;\Bigg[r\frac{d^2\theta}{dt^2}&space;&plus;&space;2\frac{dr}{dt}\frac{d\theta}{dt}}\Bigg]^2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?a&space;=&space;\sqrt{\Bigg[\frac{d^2r}{dt^2}&space;&plus;&space;r\bigg(\frac{d\theta}{dt}\bigg)^2\Bigg]^2&space;&plus;&space;\Bigg[r\frac{d^2\theta}{dt^2}&space;&plus;&space;2\frac{dr}{dt}\frac{d\theta}{dt}}\Bigg]^2" title="a = \sqrt{\Bigg[\frac{d^2r}{dt^2} + r\bigg(\frac{d\theta}{dt}\bigg)^2\Bigg]^2 + \Bigg[r\frac{d^2\theta}{dt^2} + 2\frac{dr}{dt}\frac{d\theta}{dt}}\Bigg]^2" /></a>.

The approximation of the derivatives with a two-point formula (end nodes) and three-point formula (interior nodes) delivers this graph 
\
<img src="RadarSpeed.jpg" width="500">

# Code 

To obtain the first order derivative just apply the three-point formula for interior nodes 

<a href="https://www.codecogs.com/eqnedit.php?latex=f'(x)&space;\approx&space;\frac{f(x&plus;h)-f(x-h)}{2h}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?f'(x)&space;\approx&space;\frac{f(x&plus;h)-f(x-h)}{2h}" title="f'(x) \approx \frac{f(x+h)-f(x-h)}{2h}" /></a>

for example for the instantaneous change in the distance to the airplane...
```r
deriv.r <- (lead(df$r)-lag(df$r))/(lead(df$t)-lag(df$t))
```
For the end points we have to use the two point formula

<a href="https://www.codecogs.com/eqnedit.php?latex=f'(x)&space;\approx&space;\frac{f(x&plus;h)-f(x)}{h}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?f'(x)&space;\approx&space;\frac{f(x&plus;h)-f(x)}{h}" title="f'(x) \approx \frac{f(x+h)-f(x)}{h}" /></a>

like this

```r
  deriv.r[1] <-  (df$r[2] - df$r[1])/(df$t[2] - df$t[1])
  deriv.r[l] <- (df$r[l-1] - df$r[l])/(df$t[l-1] - df$t[l])
```

For computing the second-order derivative, that is the acceleration, we can just compute the instantaneous change in the first-order derivatives.
For example for the distance to the airplane like this:

```r
deriv.r.sec <- (lead(deriv.r) - lag(deriv.r))/
    (lead(df$t)-lag(df$t))
```

Then we repeat this for the derivations for theta (the angle of the radar) and plug in the derivations in the formula for the velocity and acceleration given above.

Then we finish up with plotting the velocity and acceleration with ggplot().

<img src="RadarSpeed.jpg" width="500">

# Sources
Amos Gilat, Vish Subramaniam. *Numerical Methods for Engineers and Scientists - An Introduction with Applications Using Matlab*. Wiley, 2014. \

You can find some similar problems in \
Jaan Kiusalaas. *Numerical Methods in Engineering with MATLAB*. Cambridge University Press, p. 182 - 195, 2005.
