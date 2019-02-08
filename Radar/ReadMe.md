# <p align='center'>Measuring Velocity and Acceleration with Numerical Differentiation</p>
Solution to problem 8.37 in Gilat \& Subramaniam (2014)

# Content
* Using Radars to Measure Speed & Velocity
* Code
* Sources

# Using Radars to Measure Speed & Velocity
Any radar that measures speed of a moving object uses numerical differentiation. For example, a radar station which measures the speed of passing planes uses numerical differentiation to obtain their speed and acceleration. This example is taken from Gilat \& Subramaniam (2014). Assume the angle and distance towards an aircraft is recorded every three seconds. The velocity and acceleration can be determined using the distance to the airplane (r), the angle of the radar to the ground (theta) and the time (t) by ...

<a href="https://www.codecogs.com/eqnedit.php?latex=v&space;=&space;\sqrt{\bigg(\frac{dr}{dt}\bigg)^2&space;&plus;&space;\bigg(r\frac{d\theta}{dt}\bigg)^2}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?v&space;=&space;\sqrt{\bigg(\frac{dr}{dt}\bigg)^2&space;&plus;&space;\bigg(r\frac{d\theta}{dt}\bigg)^2}" title="v = \sqrt{\bigg(\frac{dr}{dt}\bigg)^2 + \bigg(r\frac{d\theta}{dt}\bigg)^2}" /></a> \
 and \
<a href="https://www.codecogs.com/eqnedit.php?latex=a&space;=&space;\sqrt{\Bigg[\frac{d^2r}{dt^2}&space;&plus;&space;r\bigg(\frac{d\theta}{dt}\bigg)^2\Bigg]^2&space;&plus;&space;\Bigg[r\frac{d^2\theta}{dt^2}&space;&plus;&space;2\frac{dr}{dt}\frac{d\theta}{dt}}\Bigg]^2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?a&space;=&space;\sqrt{\Bigg[\frac{d^2r}{dt^2}&space;&plus;&space;r\bigg(\frac{d\theta}{dt}\bigg)^2\Bigg]^2&space;&plus;&space;\Bigg[r\frac{d^2\theta}{dt^2}&space;&plus;&space;2\frac{dr}{dt}\frac{d\theta}{dt}}\Bigg]^2" title="a = \sqrt{\Bigg[\frac{d^2r}{dt^2} + r\bigg(\frac{d\theta}{dt}\bigg)^2\Bigg]^2 + \Bigg[r\frac{d^2\theta}{dt^2} + 2\frac{dr}{dt}\frac{d\theta}{dt}}\Bigg]^2" /></a>.

The approximation of the derivatives with a two-point formula (end nodes) and three-point formula (interior nodes) delivers the following graph \\
<img src="RadarSpeed.jpg" width="500">

# Code 


# Sources
Amos Gilat, Vish Subramaniam. *Numerical Methods for Engineers and Scientists - An Introduction with Applications Using Matlab*. Wiley, 2014. \

You can find some similar problems in \
Jaan Kiusalaas. *Numerical Methods in Engineering with MATLAB*. Cambridge University Press, p. 182 - 195, 2005.
