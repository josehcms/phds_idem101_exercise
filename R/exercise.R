rm( list = ls() )
graphics.off()

require(tidyverse);require(data.table)
source('R/MortSmooth_bbase.R')
load('data/sinasc_sp_phds.RData')

bw_dt = sinasc[ , .( bweight, time, moage, moage5, preterm, sex ) ]
rm( sinasc )


head(bw_dt)

### Plot 1 - birth weight by time
bw_dt[ , .( mean = mean( bweight ), n = .N ), time  ] %>%
  ggplot() +
  geom_point( aes( x = time, y = mean ) ) +
  geom_line( aes( x = time, y = mean ) ) +
  theme_classic()

bw_dt[ , .( mean = mean( bweight ), n = .N ), .( time, preterm )  ] %>%
  ggplot() +
  geom_point( aes( x = time, y = mean, color = factor( preterm ) ) ) +
  geom_line( aes( x = time, y = mean, color = factor( preterm ) ) ) +
  facet_wrap( ~ preterm, scales = 'free' ) +
  theme_classic() +
  theme(
    legend.position = 'bottom'
  )

mean_dt = 
  bw_dt[ , .( mean = mean( bweight ), n = .N ), .( time, cat = preterm )  ] %>%
  .[ , cat := ifelse( cat == 1, 'preterm', 'term' ) ]

data1 = bw_dt[ , .( count = .N ), .( bweight, time, preterm ) ]


### Linear model
x1 = data1$time
x2 = data1$preterm
w  = data1$count 
y  = data1$bweight
n <- nrow(data1)
X <- cbind(1, x1, x2)

### Linear weighted
tXWX <- t(X)%*%( w * X )
tXWy <- t(X)%*%( w * y )
betas_wls <- solve(tXWX, tXWy)
betas_wls
ywls <- X %*% betas_wls

x1s = seq( min(x1), max(x1), length = 100 )

resultsLM =
  data.table(
    x1s = rep( x1s, 2 ),
    ypred = c( betas_wls[1] + x1s * betas_wls[2],
               betas_wls[1] + x1s * betas_wls[2] + betas_wls[3] ),
    cat   = c( rep( 'term', 100 ), rep( 'preterm', 100 ) ),
    model = 'Linear'
  )

ggplot( ) +
  geom_point( data = resultsLM, aes( x = x1s, y = ypred, color = cat ) ) +
  geom_line( data = resultsLM, aes( x = x1s, y = ypred, color = cat ) ) +
  geom_point( data = mean_dt, aes( x = time, y = mean, color = cat ) ) +
  geom_line( data = mean_dt, aes( x = time, y = mean, color = cat ) ) +
  facet_wrap( ~ cat, scales = 'free' ) +
  theme_classic() +
  theme(
    legend.position = 'bottom'
  )

### Non-linear
library(MortalitySmooth)
B1 <- MortSmooth_bbase(x=x1, xl=min(x1), xr=max(x1), ndx=24, deg=3)
dim(B1)

# xs <- seq(min(x1), max(x1), length=10000)
# Bs <- MortSmooth_bbase(x=xs, xl=min(x1), xr=max(x1), ndx=12, deg=3)
# dim(Bs)
# matplot(xs, Bs, t="l", lty=1, xlab="time", ylab="B")


B = cbind( 1, B1, x2 )


# add penalty
nb = ncol(B)
nb1 <- ncol(B1)

D <- diff(diag(nb1), diff=2)

tDD <- t(D)%*%D

lambda <- 1.25

P1 <- lambda*tDD

P = matrix( 0, nb, nb )
P[1:nb1+1,1:nb1+1] = P1

tBWB <- t(B)%*%( w * B )

Pr <- 10^-6 * diag(nb)
Pr[1,1]  <- 0
Pr[nb,nb]  <- 0

tBBWpP <- tBWB + P + Pr
tBwy <- t(B)%*%( w * y )

betasPS <- solve(tBBWpP, tBwy)
yPS <- B %*% betasPS

x1s <- seq(min(x1),max(x1),length=100)
B1s <- MortSmooth_bbase(x1s,min(x1),max(x1),ndx=24,deg=3)

### extract betas
beta0 <- betasPS[1]
betas1 <- betasPS[ 1:nb1+1 ]
betas2 <- betasPS[nb]

# f1 <- B1s%*%betas1
# plot(x1s,f1)
my.mu.pre <- beta0+B1s%*%betas1 + betas2 * 1
my.mu.nonpre <- beta0+B1s%*%betas1 + betas2 * 0

resultsGAM1 =
  data.table(
    x1s = rep( x1s, 2 ),
    ypred = c( beta0+B1s%*%betas1 + betas2 * 0,
               beta0+B1s%*%betas1 + betas2 * 1 ),
    cat   = c( rep( 'term', 100 ), rep( 'preterm', 100 ) ),
    model = 'GAM1'
  )

ggplot( ) +
  geom_point( data = resultsLM, aes( x = x1s, y = ypred, color = model ) ) +
  geom_line( data = resultsLM, aes( x = x1s, y = ypred, color = model ) ) +
  geom_point( data = resultsGAM1, aes( x = x1s, y = ypred, color = model ) ) +
  geom_line( data = resultsGAM1, aes( x = x1s, y = ypred, color = model ) ) +
  geom_point( data = mean_dt, aes( x = time, y = mean ), shape = 4 ) +
  geom_line( data = mean_dt, aes( x = time, y = mean ), linetype = 'dotted' ) +
  facet_wrap( ~ cat, scales = 'free' ) +
  theme_classic() +
  theme(
    legend.position = 'bottom'
  )


plot(x1s,my.mu.nonpre)
plot(x1s,my.mu.pre)

########


B = cbind( 1, x2, x2*B1, (!x2)*B1 )

# add penalty
nb = ncol(B)
nb1 <- ncol( (x2*B1) )
nb2 <- ncol( ( (!x2)*B1 ) )

D1 <- diff(diag(nb1), diff=2)
D2 <- diff(diag(nb2), diff=2)

tDD1 <- t(D1)%*%D1
tDD2 <- t(D2)%*%D2

lambda1 <- lambda2 <- 1.2

P1 <- lambda1*tDD1
P2 <- lambda2*tDD2

P = matrix( 0, nb, nb )
P[1:nb1+1 + 1,1:nb1+1 + 1] <- P1
P[1:nb2+nb1+1 + 1,1:nb2+nb1+1 + 1] <- P2

tBWB <- t(B)%*%( w * B )

Pr <- 10^-6 * diag(nb)
Pr[1,1]  <- 0
Pr[2,2]  <- 0

tBBWpP <- tBWB + P + Pr
tBwy <- t(B)%*%( w * y )

betasPS <- solve(tBBWpP, tBwy)
yPS <- B %*% betasPS

x1s <- seq(min(x1),max(x1),length=100)
B1s <- MortSmooth_bbase(x1s,min(x1),max(x1),ndx=24,deg=3)

#x2s <- seq(min(x2*x1),max(x2*x1),length=1000)
#B2s <- MortSmooth_bbase(x2s,min(x2),max(x2),ndx=12,deg=3)

### extract betas
betas0 <- betasPS[1]
betas1 <- betasPS[2]
betas2 <- betasPS[ 1:nb1+1 + 1 ]
betas3 <- betasPS[ 1:nb1+nb2+1 + 1 ]

# f1 <- B1s%*%betas1
# f2 <- B2s%*%betas2

my.mu.pre <- betas0 + 1 * betas1 + (B1s*1) %*% betas2 + (B1s*0) %*% betas3
my.mu.nonpre <- betas0 + 0 * betas1 + (B1s*0) %*% betas2 + (B1s*1) %*% betas3

resultsGAM2 =
  data.table(
    x1s = rep( x1s, 2 ),
    ypred = c( my.mu.nonpre, my.mu.pre ),
    cat   = c( rep( 'term', 100 ), rep( 'preterm', 100 ) ),
    model = 'GAM2'
  )

ggplot( ) +
  geom_point( data = resultsLM, aes( x = x1s, y = ypred, color = model ) ) +
  geom_line( data = resultsLM, aes( x = x1s, y = ypred, color = model ) ) +
  geom_point( data = resultsGAM1, aes( x = x1s, y = ypred, color = model ) ) +
  geom_line( data = resultsGAM1, aes( x = x1s, y = ypred, color = model ) ) +
  geom_point( data = resultsGAM2, aes( x = x1s, y = ypred, color = model ) ) +
  geom_line( data = resultsGAM2, aes( x = x1s, y = ypred, color = model ) ) +
  geom_point( data = mean_dt, aes( x = time, y = mean ), shape = 4 ) +
  geom_line( data = mean_dt, aes( x = time, y = mean ), linetype = 'dotted' ) +
  facet_wrap( ~ cat, scales = 'free' ) +
  theme_classic() +
  theme(
    legend.position = 'bottom'
  )



plot(x1s,my.mu.nonpre,ylim=range(my.mu.nonpre,my.mu.pre))
lines(x1s,my.mu.pre)

par(mfrow=c(1,2))
plot(x1s,my.mu.nonpre)
plot(x1s,my.mu.pre)

length(x1s)
length(my.mu.nonpre)
