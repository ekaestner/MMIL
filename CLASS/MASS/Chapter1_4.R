install.packages('MASS')
library(MASS)

## Chapter 1 ############################################################################################
rm(list=ls())

#
ttt = rnorm(1000)
sss = rnorm(1000)

truehist( c(ttt,sss) , nbins=25)
contour( ddd <- kde2d(ttt,sss))
image(ddd)

#
ttt = seq(1,20,0.5)
www = 1+ttt/2
yyy = ttt+www*rnorm(ttt)

dum = data.frame( ttt , www , yyy)
rm( ttt , www , yyy )

fm = lm( ttt ~ yyy , data=dum)
summary(fm)

fm1 = lm( ttt ~ yyy , data=dum , weight = www)
summary(fm1)

lrf = loess( yyy ~ ttt , dum )

plot( dum$ttt , dum$yyy )
lines( spline( dum$ttt , fitted(lrf)) , col=2)
abline(0,1,lty=3,col=3)
abline(fm,col=4)
abline(fm1,col5,lty=4)

plot( fitted(fm) , resid(fm) , xlab='Fitted' , ylab='Residuals' )

qqnorm(resid(fm))
qqline(resid(fm))

#
attach(hills)
hills
pairs(hills)
plot(dist,time)
identify(dist,time,row.names(hills))
abline(lqs(time~dist),lty=3,col=4)
detach()

