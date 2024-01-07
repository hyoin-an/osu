##########################################################################################
## Script adapted from
## Data analysis 1 - Rizhoctonia root rot data set #######################################
## Binomial spatial generalized linear mixed models ######################################
## Author: Wagner Hugo Bonat LEG/IMADA ###################################################
## Date: 18/06/2014 ######################################################################
##########################################################################################

## Load necessary packages
require(geoR)
require(Matrix)
require(bbmle)
require(akima)

## Loading the data set
load("rhizoc.rdata")

## Loading extra functions
source("functionssglmm.r")

## Spherical covariance function
spherical.cov.mat <- function(d, sigma = 1, phi = 1, tau = 1){
  return((sigma * (1 - 3/2 * (d/phi) + 1/2 * (d/phi)^3) + tau * ifelse(d == 0, 1, 0)) * ifelse(d < phi, 1, 0))
}

## Logit function
logit <- function(p){
  log(p/(1-p))
}

## Inverse logit function
inv.logit <- function (etas) {
  exp(etas)/(1+exp(etas))
}

###################################################
#### Data Manipulation and EDA
###################################################
## Change the names
names(Rhizoc) <- c("coord.X","coord.Y", "N.trials", "Y", "Yield")

## Calculate p and logit_p for each location
Rhizoc$p <- Rhizoc$Y / Rhizoc$N.trials
Rhizoc$logit_p <- logit(Rhizoc$p)

## EDA plots
hist(Rhizoc$p, main = 'Histogram of Proportion of Root Disease', xlab = 'p')

png('figures/root_prop_map.png', width = 960, height = 640, res = 120)
plot(coord.Y ~ coord.X, cex = p * 7, data = Rhizoc, main = 'Proportion of Root Disease by Location', xlab = 'x', ylab = 'y')
dev.off()


###################################################
#### Model Fitting
###################################################
## Get initial values
initial = start.values.glgm(Y ~ 1, family="binomial", data= Rhizoc, coords = Rhizoc[,1:2],nugget=TRUE, ntrial=Rhizoc$N.trials)
data.geo <- as.geodata(data.frame(Rhizoc[,1:2],res))
plot(variog(data.geo)) ## Look here a big nugget effect, so in this case my initial point is not good, because the nugget is bigger than partial sill.
initial[4] <- 1.1*initial[2]

## Fitting Binomial SGLMM
fit1 <- glgm(Y ~ 1, cov.model = "spherical", kappa = log(0.5), inits = initial, data=Rhizoc, coords = Rhizoc[,1:2],
             nugget=TRUE, family="binomial", ntrial = Rhizoc$N.trials,
             method.optim = "BFGS", method.integrate = "NR")
fit1[[9]]
LAPLACE <- coef(fit1[[1]])[,1] # Get point estimates from laplace model


## Make predictions and draw map
# Create matrices and data frame that will hold predictions for mapping
im <- interp(Rhizoc$coord.X, Rhizoc$coord.Y, rep(0, nrow(Rhizoc)))
pred.locs <- as.data.frame(which(!is.na(im$z), arr.ind = T))
pred.data <- cbind(pred.locs, t(apply(pred.locs, 1, function(k) c(im$x[k[1]], im$y[k[2]]))))
colnames(pred.data)[3:4] <- c('coord.X','coord.Y')

all.locs <- rbind(Rhizoc[,c('coord.X','coord.Y')], pred.data[,c('coord.X','coord.Y')])

# Setup components for prediction
n <- nrow(Rhizoc)
k <- nrow(pred.data)

D <- dist(all.locs)
cov.mat <- as.matrix(spherical.cov.mat(D, sigma = LAPLACE[2], phi = LAPLACE[3], tau = LAPLACE[4])) + (LAPLACE[2] + LAPLACE[4]) * diag(nrow(all.locs))
Sigma_obs <- cov.mat[1:n, 1:n]
Sigma_pred <- cov.mat[(n+1):(n+k), (n+1):(n+k)]
Sigma_cross <- cov.mat[1:n, (n+1):(n+k)]

B <- solve(Sigma_obs, Sigma_cross)

mu_obs <- rep(1,n) * LAPLACE[1]
mu_pred <- rep(1,k) * LAPLACE[1]

# Estimated BLUP using Universal Kriging Predictor
pred.data$pred_logit <- mu_pred + crossprod(B, (Rhizoc$logit_p - mu_obs))
pred.data$sd <- sqrt(diag(Sigma_pred - crossprod(B, Sigma_cross)))

# Convert back to original scale by taking invesrse logit
pred.data$pred_p <- inv.logit(pred.data$pred_logit)
pred.data$p_upper <- inv.logit(pred.data$pred_logit + qnorm(.975) * pred.data$sd)
pred.data$p_lower <- inv.logit(pred.data$pred_logit - qnorm(.975) * pred.data$sd)


## Map the results
im$pred_logit <- im$pred_p <- im$p_upper <- im$p_lower <- im$sd <- im$z
for(k in 1:nrow(pred.data)){
  cur_loc <- pred.data[k,]
  im$pred_logit[cur_loc$row, cur_loc$col] <- cur_loc$pred_logit
  im$pred_p[cur_loc$row, cur_loc$col] <- cur_loc$pred_p
  im$p_upper[cur_loc$row, cur_loc$col] <- cur_loc$p_upper
  im$p_lower[cur_loc$row, cur_loc$col] <- cur_loc$p_lower
  im$sd[cur_loc$row, cur_loc$col] <- cur_loc$sd
}

pal <- colorRampPalette(c("black", "white"))
Rhizoc$color <- pal(10)[as.numeric(cut(Rhizoc$logit_p, breaks = 10))]

# Predicted logit map
par(mfrow=c(1,1), cex=0.75, mar=c(4,4,1,1), bty="L")
plot(Rhizoc[,c('coord.X','coord.Y')], type = 'n')
with(im, image(x,y,pred_logit, main = "Predicted logit(p)"))
points(Rhizoc[,c('coord.X','coord.Y')], col = Rhizoc$color, pch = 19)

# Predicted p map
# zlim = c(min(im$p_lower, na.rm = T), max(im$p_upper, na.rm = T))


png('figures/p_predict_map.png', width = 960, height = 640, res = 120)
par(cex=0.75, mar=c(4,4,1,1), bty="L")
plot(Rhizoc[,c('coord.X','coord.Y')], type = 'n')
with(im, image(x,y,pred_p, main = "Predicted probability", xlab = '', ylab = ''))
points(Rhizoc[,c('coord.X','coord.Y')], cex = 7 * Rhizoc$p, pch = 19, col = 'black')
dev.off()

# Upper 95% prediction interval
png('figures/p_upper95_map.png', width = 960, height = 640, res = 120)
par(cex=0.75, mar=c(4,4,1,1), bty="L")
plot(Rhizoc[,c('coord.X','coord.Y')], type = 'n')
with(im, image(x,y,p_upper, main = "Upper 95% for p", xlab = '', ylab = ''))
dev.off()

# Lower 95% prediction interval
png('figures/p_lower95_map.png', width = 960, height = 640, res = 120)
par(cex=0.75, mar=c(4,4,1,1), bty="L")
plot(Rhizoc[,c('coord.X','coord.Y')], type = 'n')
with(im, image(x,y,p_lower, main = "Lower 95% for p", xlab = '', ylab = ''))
dev.off()

# Standard Deviation map
png('figures/p_sd_map.png', width = 960, height = 640, res = 120)
par(cex=0.75, mar=c(4,4,1,1), bty="L")
plot(Rhizoc[,c('coord.X','coord.Y')], type = 'n')
with(im, image(x,y,sd, main = "Standard Deviation", xlab = '', ylab = ''))
points(Rhizoc[,c('coord.X','coord.Y')], pch = 19)
dev.off()






