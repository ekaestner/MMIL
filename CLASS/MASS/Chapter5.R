library(MASS)

## Probability Distributions ############################################################################################
# Density functions (continuous) and Functions (discrete)
# Desnity (probability functions), cumulative density function, quantile function (inverse of CDF) = dNAME, cNAME, qNAME, rNAME {dnorm,cnorm,qnorm,rnorm}
qt(.975,11) # 5% critical value for t on 11 degrees of freedom

qqplot # compare two distributions
qqnorm # compare sample with normal distribution

x=10
plot( qt(ppoints(x),9) , sort(ppoints(10)))

x = rnorm(250)
qqnorm(x); qqline(x)

#
ttt = rnorm(5000)
fitdistr(ttt,"normal") # fits standard univariate distributions

#
contam = rnorm( 100 , 0 , (1+2*rbinom(100,1,.05)) )
regul = rnorm(100 , 0 , 1)
qqplot( contam , regul ); qqline( contam , regul )

#
sample(100)
sample(contam)
sample(contam,replace=T)

ttt = c( 1, 2, 3, 4, 5)
ppp = c(.01,.01,.6,.37,.01)
sample( ttt , 100, replace=T, prob=ppp)

## Data Examination ############################################################################################

summary(sample( ttt , 100, replace=T, prob=ppp) )# summary function retruns mean, quartiles, number of missing values
attach(geyser)

var(geyser$waiting) #var,  (variance-cosvariance), cor (correlations), cov.wt, quantile, max, min, range, stem
std(geyser$waiting)
cor(geyser) 
quantile(geyser$duration)
range(waiting)

 # plots include hist, boxplot, stem
boxplot(geyser)
par( mfrow = c(1,2) )
boxplot(chem,sub="chem",range=.5)
boxplot( abbey, sub = "abbey")

par(c(1,1))
library(lattice)
bwplot(type~y | meas, data=fgl,scales=list(x="free"))

## Univariate Statistics ############################################################################################
attach(shoes)

t.test(A, mu=8.88)
t.test(A)$conf.int

wilcox.test(A,mu=10)

var.test(A,B)
t.test(A,B,var.equal=T)
t.test(A,B,var.equal=F)
wilcox.test(A,B)
wilcox.test(A,B,paired=T)

# Full list of classical tests
?binom.test
?friedman.test
?prop.test
?chisq.test
?kruskal.test
?t.test
?cor.test
?mantelhaen.test
?var.test
?fisher.test
?mcnemar.test
?wilcox.test
  
## 5.5 Robust Summaries ############################################################################################
# robust/resistant in relation to outliers
# Mean is only optimal for certain distributions of data
# Different variance estimations: variance, inter-quartile range, median absolute deviation

detach()

sort(chem)

mean(chem)
median(chem)
mad(chem)

unlist(huber(chem))
unlist(hubers(chem))

fitdistr(chem,"t",list(m=3,s=0.5),df=5)

# trimmed mean, removes values on each side of distribution
ttt = c( -10, -10, seq( 0 , 10 , .5))
mean(ttt)
mean(ttt,trim=.1)

## 5.6 Density Estimation ############################################################################################
install.packages('sm')
library(sm)

density # bandwidth specified by bw=

truehist(duration , nbins=25, xlim = c(0.5, 6), ymax = 1.2)
lines( density( duration, width = "SJ", n=256) , lty = 3 )
lines( density(duration, width = "SJ-dpi", n=256), lty=3 )

geyser2 = data.frame( as.data.frame(geyser)[-1,] , pduration = geyser$duration[-299] )
attach(geyser2)
par( mfrow = c(2,2) )
plot( pduration, waiting, xlim = c(.5,6), ylim = c(40, 110), xlab = "previous dueration", ylab = "waiting")
f1 = kde2d( pduration, waiting, n=50, lims=c(.4, 5, 40, 110) )
image( f1, zlim = c(0, .075))
f2 = kde2d( pduration, waiting, n=50, lims=c(.4, 5, 40, 110), h = c(width.SJ(duration), width.SJ(waiting) ))
image( f2, zlim = c(0, .075))
persp( f2, phi = 30, theta=20, d=5)

ttt = seq( 0, 10, .5)
length(ttt)
ttt[-21]

# Local polynomial fitting
install.packages('KernSmooth')
library(KernSmooth)

plot( x=c(0, 35000), y=c(0, .003), type = "n", bty = 'l')
rug(galaxies)
lines(bkde(galaxies,bandwidth=dpik(galaxies)))
lines(locpoly(galaxies, bandwidth=dpik(galaxies)), lty=3 )

locpoly(galaxies, bandwidth=dpik(galaxies))

## 5.7 DBootstrap and Permutation Testing ############################################################################################
gal = galaxies/1000

median(gal)

density( gal, n=1, from=20.833, to=20.834, width='SJ')$y
1/( 2 * sqrt(length(gal)) * .13)

set.seed(101)
m=1000
res=numeric(m)
for (i in 1:m) res[i] = median( sample( gal, replace=T))
gal
sample(gal, replace=T)












