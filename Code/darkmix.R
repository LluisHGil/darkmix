# March 29, 2020
# 
# ================================================================================
# Title: Mixture Models for detection and characterization of Dark Matter halos
# Authors: Lluís Hurtado-Gil, Michael A. Kuhn, Pablo Arnalte-Mur, Eric D. Feigelson, Vicent J. Martínez
# ================================================================================
# 
# Code names: darkmix.R
# Language: R
# Code tested under the following compilers/operating systems: 
#      R version 3.5.2 (2018-12-20) -- "Eggshell Igloo" on Macintosh
# 
# Functions of darkmix.R library. Each function includes a brief explanation of it purpose and the 
# definition of its input and output arguments. A usage example is provided.
#
# System requirements: R version 3.5.2 or later. Spatstat version 1.58-2 or later.
# 
# Calls to external routines: 
#    This code makes use of spatstat and spatstat.utils functions owin, pp3, unique, as.im, 
#         quadscheme, box3, variablesinformula and volume
#    This code makes use of plotrix function draw.circle.
#	   This code makes use of LaplacesDemon functions Blocks and LaplacesDemon
#	   This code makes use of expint functions gamma and gammainc.
#    This code makes use of misc3d function kde3d.
#    This code makes use of numerous Base R functions (in particular the function optim).
#
# List of Functions in this library:
#    quad.3d
#    einasto.model
#    dn.fit
#    const.model
#    centers
#    param2profile
#    profile2param
#    action.initial
#    action.multi
#    action.control
#    action.prior
#    action.draw
#    action.pop
#    mixture.model
#    model.lik
#    int.model
#    mask.freeze
#    mask.model
#    mask.adjust_positions
#    mask.adjust_rotations
#    mask.adjust_radii
#    mask.adjust_sersic
#    mask.adjust_mix
#    mask.adjust_shape
#    AIC.BIC
#    goodness
#    residuals.lambda
#    residuals.mm
#    ready.mcmc
#    run.ld
#    model.post
#    outside.win
#    image.3D
#    ds9.colors
#    american.colors
#    plot.dens
#    plot.mm3d
#    plot.mm2d
#    gen.pattern
#    profile.emp
#    volume.shells
#    mass.ein
#    profile.ind
#    profile.mm
#    inside.shell
#    inside
#    plot.profile
#    profile.sphere
#    einasto.1d
#    around
# 
# Additional comments: The R library defined in darkmix.R will be accessed by the
# finite mixture model fitting example in the article.
# 

#
# quad.3d (Creates a quadrature scheme for a 3D point process)
# 
# The input is a spatial point pattern data as a pp3 object (clust) and the 
# number of tiles per dimension in a vector (ntile, default is c(32, 32, 32)).
# The return is a grid of points as a pp3 object in the window occupied by "clust"
# Usage> quad.3d(clust=clust, ntile)
quad.3d <- function(clust=clust, ntile = c(32, 32, 32)) {
	yrange <- clust$domain$yrange
	zrange <- clust$domain$zrange
	clustyz <- ppp(x=clust$data$z[1], y=clust$data$y[1], window=owin(x=c(zrange[1],zrange[2]), y=c(yrange[1],yrange[2])))
	Qyz <- quadscheme(clustyz, ntile=c(ntile[3],ntile[2]))
	grid.yz <- cbind(x=Qyz$dummy$y, y=Qyz$dummy$x)
	grid.yz <- head(grid.yz,-4)
	grid.size <- dim(grid.yz)[1]
	
	xrange <- clust$domain$xrange
	xlong <- xrange[2] - xrange[1]
	epsx <- xlong/ntile[1]
	
	Q <- rbind()
	for(i in 1:ntile[1]) {
		x <- (i-1)*epsx+xrange[1] + epsx/2
		grid.yzi <- cbind(rep(x,grid.size),grid.yz)
		Q <- rbind(Q,grid.yzi)
	}
	Q.pp3 <- pp3(x=Q[,1], y=Q[,2], z=Q[,3], window=box3(x=xrange, y=yrange, z=zrange))	
	return(Q.pp3)
}

#
# einasto.model (The Einasto profile surface density distribution)
# 
# The inputs are x, y and z coordinates and a vector "param" containing the 
# (x0,y0,z0) coordinates of the halo center, the radius of the Einasto profile (rs) 
# and the Sérsic índex (n). Function returns value -Inf when n>30.
# The output is a value of the model (not normalized) at the location (x,y,z).
# Usage> einasto.model(0,1,2,param=c(0.0, 0.0, 0.0, 4, 3))
#
einasto.model <- function(x1,y1,z1,param=param,param2=NULL) {
  x0=param[1]
  y0=param[2]
  z0=param[3]
  re=param[4]
  n=param[5]
  if(n < 30) {
    dn <- uniroot(dn.fit,c(0,200),n)$root
    dist <- sqrt((x1-x0)^2 + (y1-y0)^2 + (z1-z0)^2)
    res <- exp(-dn*((dist/re)^(1/n)-1))
  }
  else {res <- -Inf}
  res
}

#
# dn.fit (Return the value of the dn term in the Einasto profile evaluation)
# 
# The inputs are a value x where the equation is evaluated and the Sérsic index n. 
# The call of einasto.model estimates a root for this function and takes the returned 
# value as the term dn. This function is not meant to be used by the user
# Usage> dn.fit(8.67,3)
#
dn.fit <- function(x,n) {
  g <- gamma(3*n) - 2*gammainc(3*n,x)
  return(g)
}

# 
# const.model (The constant surface density distribution)
# 
# The inputs are x, y and z coordinates.
# The output is the value 1 regardless of the input coordinates, representing a constant surface density.
# This function is a 3D version of function const.model from Kuhn et al. 2014.
# Usage> const.model(0,1,2,param2=NULL)
#
const.model <- function(x1,y1,z1,param=NULL,param2=param2) { 1.0 }

#
# centers (Finds k candidates for halos centers)
#
# Inputs are the data set (clust), the bandwidth for the kernel density estimation (h, default is 1), 
# the number of tiles per dimension (ntile) and the number of halos that we want to model (k).
# Output is a list of two elements: 1) a matrix of k rows and 4 columns and 2) the resulting 
# density map. Each row contains the x, y, z coordinates of the center and the density at this coordinate.
# Usage> centers(clust, h=0.7, ntile=c(25,25,7), k=6)
#
centers <- function(clust, h=1, ntile = c(32,32,32), nk) {
  
  window <- clust$domain
  xlong <- window$xrange[2] - window$xrange[1]
  ylong <- window$yrange[2] - window$yrange[1]
  zlong <- window$zrange[2] - window$zrange[1]
  epsx2 <- xlong/(2*ntile[1])
  epsy2 <- ylong/(2*ntile[2])
  epsz2 <- zlong/(2*ntile[3])
  lims <- c(window$xrange[1]+epsx2, window$xrange[2]-epsx2, window$yrange[1]+epsy2, window$yrange[2]-epsy2, window$zrange[1]+epsz2, window$zrange[2]-epsz2)
  dens <-  kde3d(clust$data$x, clust$data$y, clust$data$z, h=h, n=ntile, lims=lims)

  T <- prod(ntile)
  dens.mat <- matrix(NA,T,4)
  mode <- matrix(0,1,4)
  l <- 1
  for(i in 1:ntile[1]) {
    for(j in 1:ntile[2]) {
      for(k in 1:ntile[3]) {
        dens.mat[l,] <- c(dens[[1]][i], dens[[2]][j], dens[[3]][k], dens[[4]][i,j,k])
        if((i>1) & (i<ntile[1]) & (j>1) & (j<ntile[2]) & (k>1) & (k<ntile[3])) {
          max <- around(dens,i,j,k)
          if(max == TRUE) {
            mode <- rbind(mode, c(dens[[1]][i], dens[[2]][j], dens[[3]][k], dens[[4]][i,j,k]))
          }
        }
        l <- l + 1
      }
    }
  }
  mode <- mode[order(mode[,4], decreasing=TRUE),]
  return(list(mode[1:nk,], dens.mat))
}

#
# param2profile (Converts an array of parameters into an R data frame)
# 
# The input is an array of parameters in the variable "param" as described
# for the function mixture.model. The assumed model form is
# for a model with "param2" of c(<n+1>,model.einasto,5, ... , model.einasto,5,model.const,0)
# The output is a data.frame object with columns, x, y, z, re, n, and mix. 
# This function is a 3D version of function param2ellipse from Kuhn et al. 2014.
# Usage> profiles <- param2profile(param)
# 
param2profile <-function(param){
  v <- 6
  k <- length(param)/v # i.e. six parameters for Einasto profile, minus one for the first mix, plus one for constant component
  profile <- data.frame(x=rep(NA,k),y=rep(NA,k),z=rep(NA,k),re=rep(NA,k),n=rep(NA,k),mix=rep(NA,k))
  adjust <- 0 # compensate for missing mix[1]
  for (i in 1:k) {
    for(j in 1:v) {
      profile[i,j] <- param[v*(i-1)+j-adjust]  		
    }
    if (i > 1) {profile[i,v] <- param[v*(i-1)+v-adjust]}	
    if (i == 1) {profile[i,v] <- 0.0}	
    adjust <- 1
  }
  print(profile)
  cat("Background mix: "); cat(tail(param,1));
  cat("\n")
  return(list(profile=profile, bmix=tail(param,1)))
}

#
# profile2param (Converts an R data frame into an array of parameters)
# 
# The input is a data.frame object with columns, x, y, z, re, n, 
# and mix, describing a collection of Einasto profiles. The variable
# "bmix" is the log mixing parameter for the background component.
# The output is an array of parameters in the format of the "param" argument
# to the function mixture.model. The assciated "param2" variable for this 
# "param" array has the form c(<n+1>,model.einasto,5, ... , model.einasto,5,model.const,0)
# This function is a 3D version of function ellipse2param from Kuhn et al. 2014.
# Usage> param <- profile2param(profile)
# 
profile2param <-function(profile) {
  param <- c(t(as.matrix(profile[[1]])))
  if (param[6] != 0.0) {print("Error! The first mixture coefficent is not zero")}
  param <- param[-6]
  return(c(param,profile[[2]]))
}

#
# action.initial (creates initial vecter of parameters)
#
# Inputs are an empty vector for the model parameters (param), the variables from function 
# mixture.model containing information regarding the model component and its parameters
# (i, par1, par2, n_param, n_models, p1, model) and the extra argument with information to 
# initialize the model parameters (arg). This information includes the output of function centers
# and candidate values for the radius (rs)s, the Sérsic index (n), the mixture coefficient (mix) 
# and the background mixture coefficient (back). The output is a vector with the initial values of 
# the model parameters (param).
# This function is used as an argument by function mixture.model, not to be used by the user.
# cent <- centers(clust, h, ntile, k)
# arg = list(centers=cent, re=re, n=n, mix=mix, back=back) 
#
action.initial <- function(param=c(),i,par1,par2,n_param,n_models,p1,model=NULL, arg) {
  if(i < n_models) {
    param[p1:(p1+n_param-1)] <- c(arg$centers[i,1:3], arg$re, arg$n)
    if (i > 1) { param[p1+n_param] <- arg$mix }
  }
  else {
    param[p1+n_param] <- arg$mix
  }
  if(i == n_models) {param[p1+n_param] <- arg$back}
  return(param)
}

#
# action.multi (Evaluates a component of the mixture model in data set)
# 
# Inputs are the vector with the model evaluation (output), the variables from function 
# mixture.model containing information regarding the model component and its parameters
# (i, par1, par2, n_param, n_models, p1, model) and the extra argument with the locations
# where the model is to be evaluated (arg). 
# After the second, third, fourth (and so on) model parameters, there is a mixture coefficient 
# indicating the weight of this model relative to the first model component. I.e. For a mixture
# component for the second component log mix = -1, the peak surface density of component 2 
# would be one tenth the peak surface density of component 1.
# The output is a value of the model (not normalized) at the location (x,y,z).
# This function is used as an argument by function mixture.model, not to be used by the user.
#
action.multi <- function(output,i,par1,par2,n_param,n_models,p1,model, arg) {
  # arg <- list(x1 = data[,1], y1 = data[,2], z1 = data[,3])
  param_new=NULL
  if (n_param > 0) { param_new=par1[1:n_param] }
  scale <- 1.0
  if ( i == 1) {output <- scale*model(arg$x1,arg$y1,arg$z1,param=param_new)}
  if (i > 1) { 
    scale <- 10.0^(par1[[n_param+1]]) 
    output <- output + scale*model(arg$x1,arg$y1,arg$z1,param=param_new)
  }
  
  output
}

#
# action.member ()
#
#
#
action.member <- function(output,i,par1,par2,n_param,n_models,p1,model, arg) {
  param_new=NULL
  if (n_param > 0) { param_new=par1[1:n_param] }
  scale <- 1.0
  if ( i == 1) {output[,i] <- scale*model(arg$x1,arg$y1,arg$z1,param=param_new)}
  if (i > 1) { 
    scale <- 10.0^(par1[[n_param+1]]) 
    output[,i] <- scale*model(arg$x1,arg$y1,arg$z1,param=param_new)
  }
  
  output
}

#
# action.control (Retuns TRUE if the proposed set of parameters is not valid)
# 
# Inputs are a logical variable (alert), the variables from function 
# mixture.model containing information regarding the model component and its parameters
# (i, par1, par2, n_param, n_models, p1, model) and the extra argument with information to 
# initialize the model parameters (arg). This information only includes the data window.
# The output is a boolean variable stating wether the parameters satisfies the window conditions
# and some other conditions.
# This function is used as an argument by function mixture.model, not to be used by the user.
# arg = win
#
action.control <- function(alert,i,par1,par2,n_param,n_models,p1,model,arg) {
  if(i == 1) {alert=FALSE}
  param_new=NULL
  if(i < n_models) {
    if(n_param > 0) { param_new=par1[1:n_param] }
    if(outside.win(param_new, arg$window)) {alert=TRUE; print(paste("Center of halo", i, "is outside the window."))}
    if(!is.na(param_new[4])){
      if(param_new[4] <= 0) {alert=TRUE; print(paste("Halo", i, "has negative radius:", param_new[4]))}
    }
    if(!is.na(param_new[5])){
      if(param_new[5] <= 0) {alert=TRUE; print(paste("Halo", i, "has negative Sérsic index:", param_new[5]))}
      if(param_new[5] > 30) {alert=TRUE; print(paste("Halo", i, "has Sérsic index > 30:", param_new[5]))}
    }
  }
  return(alert)
}

#
# action.prior (Returns de logarithm of the prior probabilities)
#
# Inputs are an empty vector for the log priors (logp), the variables from function 
# mixture.model containing information regarding the model component and its parameters
# (i, par1, par2, n_param, n_models, p1, model) and the extra NULL argument (arg). 
# Output is a vector containing the log prior of each parameters. 
# Gaussian priors with large variances are used to simulate nearly flat priors.
# This function is used as an argument by function mixture.model, not to be used by the user.
#
action.prior <- function(logp,i,par1,par2,n_param,n_models,p1,model,arg) {
  param_new=NULL
  if(i < n_models) {
    if (n_param > 0) { param_new=par1[1:n_param] }
    lpx = dnorm(param_new[1], 0, 1000, log=TRUE)
    lpy = dnorm(param_new[2], 0, 1000, log=TRUE)
    lpz = dnorm(param_new[3], 0, 1000, log=TRUE)
    lpr = dnorm(param_new[4], 1, 1000, log=TRUE)  # dlnorm(param_new[4], 1.5, 100, log=TRUE)
    lpn = dnorm(param_new[5], 3, 1000, log=TRUE)  # dlnorm(param_new[5], 2.5, 100, log=TRUE)
    logp[p1:(p1+n_param-1)] <- c(lpx,lpy,lpz,lpr,lpn)
    
    if (i > 1) {
      lpmix = dnorm(par1[n_param+1], 0, 1000, log=TRUE)
      logp[p1+n_param] <- lpmix
    } 
  }
  else {
    logp <- c(logp, dnorm(par1,0.5,1000, log=TRUE))
  }
  return(logp)
}

#
# action.draw (Overplots model radii on an image of the data set)
# 
# Inputs are a NULL output variable (output), the variables from function 
# mixture.model containing information regarding the model component and its parameters
# (i, par1, par2, n_param, n_models, p1, model) and an extra argument with plotting 
# parameters and the data set to be plotted (arg). 
# There is no output. This function plots the spatial point pattern (variables X and Y) 
# with circles of radius r_e.
# This function is used as an argument by function mixture.model, not to be used by the user.
# This function is a 3D version of function draw.model with circles from Kuhn et al. 2014.
# arg = list(mar=c(0,0,0,0), clust=clust, col=2, pch=16, border=2, cex=0.5, lwd=2)
#
action.draw <- function(output=NULL,i,par1,par2,n_param,n_models,p1,model,arg) {
  if(i == 1) {
    par(mar=arg$mar)
    plot(arg$clust$data$x,arg$clust$data$y,cex=arg$cex,pch=arg$pch,axes=TRUE,cex.lab=2,cex.axis=2,xlab=expression(paste("h"^-1, " Mpc", sep="")), ylab=expression(paste("h"^-1, " Mpc", sep="")))	}
  param_new=NULL
  if (n_param > 0) { param_new=par1[1:n_param] }
  if (!is.null(param_new)) {
    text(param_new[1],param_new[2],labels=toString(i),col=arg$col,cex=arg$cexl)
  }
  draw.circle(param_new[1],param_new[2],param_new[4],border=arg$border,lwd=arg$lwd)
}

#
# action.pop (Estimates the number of particles per component)
#
# Inputs are an empty vector where the population number will be estimated (pop), 
# the variables from function mixture.model containing information regarding the model 
# component and its parameters (i, par1, par2, n_param, n_models, p1, model) and an extra 
# argument (arg) with the data set (clust) and the grid (quad). 
# This function is used as an argument by function mixture.model, not to be used by the user.
# Output is the estimated number of points in a model component.
# arg = list(clust=clust, quad=quad)
#
action.pop <- function(pop,i,par1,par2,n_param,n_models,p1,model,arg) {
  N <- length(arg$clust$data$x)
  window <- arg$clust$domain
  low <- c(window$xrange[1],window$yrange[1],window$zrange[1])
  upp <- c(window$xrange[2],window$yrange[2],window$zrange[2])
  eps <- arg$quad$data$z[2] - arg$quad$data$z[1]
  eps3 <- eps*eps*eps
  
  if(i == 1) {pop <- c()}
  param_new=NULL
  if (n_param > 0) { param_new=par1[1:n_param] }
  scale <- 1.0
  if (i > 1) {scale <- 10.0^(par1[[n_param+1]])}
  pop[i] <- sum(scale*model(arg$quad$data$x,arg$quad$data$y,arg$quad$data$z,param=param_new)*eps3)
  if(i == n_models) {
    pop[i] <- scale*volume(window)
    norm <- sum(pop)
    pop <- N*pop/norm
  }
  return(pop)
}

# 
# mixture.model (The finite mixture model surface density distribution)
# 
# This function admits different functions (action) that can be evaluated on every model component
# with different purposes. The inputs are a vector of model parameters (param), a vector describing the 
# component models to be used (param2), the output of the selected action funciton (output), the 
# chosen action function (action) and any necessary additional arguments to function action in a list (arg)
# The vector describing the component models "param2" starts with the total number of components, 
# followed by the function names and the number of parameters to pass to each function.
# This function is a 3D version of function multi.model from Kuhn et al. 2014.
# Usage> mixture.model(param=c(0.0, 0.0, 0.0, 4.0, 3.0, 2.0, 2.0, 2.0, 5.0, 2.0, 1.0),
#        param2=c(2,einasto.model,5,einasto.model,5), output, action.multi, 
#        arg <- list(x1 = data[,1], y1 = data[,2], z1 = data[,3]))
#
mixture.model <- function(param=param,param2=param2,output=NULL,action,arg) {
    n_models <- param2[[1]]
    p1 <- 1
    p2 <- 2
    n1 <- length(param)
    n2 <- length(param2)
  
    i <- 1
    while (i < n_models+1) {
      
      par1 <- param[p1:n1]
      par2 <- param2[p2:n2]
      model <- par2[[1]]
      n_param <- par2[[2]]
      
      output <- action(output,i,par1,par2,n_param,n_models,p1,model,arg)
      
      p2 <- p2 + 2
      if (i > 1) { p1 <- p1 + 1}
      p1 <- p1 + n_param
      i <- i+1
      
    }
    output
}

#
# model.lik (Returns the negative of the log likelihood)
# 
# The inputs are the "param" and "param2" vectors (described in mixture.model), the 
# spatial point pattern data as a pp3 object "clust" and the spatial 
# point pattern quadrature "quad" (described in quad.3d).
# The output is the negative of the log likelihood for the model. (A negative value 
# is given so that optim can find the maximum likelihood by minimization.)
# This function is a 3D version of function model.lik from Kuhn et al. 2014.
# Usage> -model.lik(param,param2=param2,clust=clust,quad)
#
model.lik <- function(param,param2=param2,model,clust=clust,quad,all) {
    n <- length(clust$data$x)
    arg = list(x1=clust$data$x,y1=clust$data$y,z1=clust$data$z,window=clust$domain)
    output=NULL
    if(all==TRUE) {
      alert <- mixture.model(param,param2,output,action.control,arg)
    }
    else {alert <- FALSE}
    if(alert == TRUE) {
      # a <- param2profile(param)
      result <- Inf
    }
    else {
      norm <- int.model(param, param2, quad, model)  
      lprior <- model(param, param2, output=c(), action.prior, arg)
      result <- sum(log(model(param=param,param2=param2,output,action.multi,arg)/norm*n))-n
      result <- result # + sum(lprior)
    }
    -result
}

#
# int.model (Integrates the finite mixture model in a given window)
# 
# Inputs are the "param" and "param2" vectors (described in mixture.model) and 
# the quadrature of point process, a grid of points covering the window W. 
# The output is the numerical integration of the model in the observation window 
# before normalization
# Usage> int.model(model=mixture.model,action=action.multi,param=c(0.0, 0.0, 0.0, 4.0, 3.0, 2.0, 2.0, 2.0, 5.0, 2.0, 1.0),
#        param2=c(2,einasto.model,5,einasto.model,5),quad=quad)
int.model <- function(param, param2, quad, model) {
	eps <- quad$data$z[2] - quad$data$z[1]
	arg = list(x1=quad$data$x,y1=quad$data$y,z1=quad$data$z,window=quad$domain)
	output <- c()
	eval <- model(param, param2, output, action.multi, arg)
	norm <- sum(eval*eps*eps*eps)
	norm	
}

#
# mask.freeze (this optimization procedure will only allow certain parameters to be fit by optim, while the others stay fixed)
# 
# The input are the arrays "param" and "param2", a vector of 1s and 0s "mask" which indicates 
# which parameters in "param" should be fit (1) or frozen (0), and the ppp object clust and 
# quad. The output is an optimized parameter array.
# This function is obtained from Kuhn et al. 2014.
# Usage> param <- mask.freeze(param,param2,mask.adjust_positions(k,v),clust=clust,quad=quad)
# 
mask.freeze <- function(param_input,param2,mask,clust=clust,quad=quad) {
  ind_masked <- which(mask == 0)
  ind_unmasked <- which(mask == 1)
  param_out <- param_input
  ocf <- optim(param_input[ind_unmasked], model.lik, param2=c(param_input[ind_masked],ind_masked), model=mask.model, clust=clust, quad=quad, all=FALSE)
  param_out[ind_unmasked] <- ocf$par
  param_out
}

#
# mask.model (model that is passed to optim by the mask.freeze function)
# 
# The input is the same as other models, x, y and z coordinates, the model "param" parameters, 
# and the model "param2" description.
# The output is the value of the model at position (x,y,z).
# This function is a 3D version of function mask.model from Kuhn et al. 2014.
# Usage> ocf <- optim(param_input[ind_unmasked],model.lik,model=mask.model,clust=clust,param2=c(param_input[ind_masked],ind_masked))
# 
#mask.model <- function(x1,y1,z1,param=NULL,param2=param2) {
mask.model <- function(param,param2=param2, output=NULL, action, arg) {

  n_masked <- length(param2)/2
  n_unmasked <- length(param)
  n_par <- n_masked + n_unmasked
  mask <- rep(1,n_par)

  ind_masked=round(param2[(n_masked+1):(n_masked*2)]) 
  mask[ind_masked] <- 0 
  ind_unmasked=which(mask == 1) 

  param_input <- rep(0.0,n_par)  
  param_input[ind_masked] <- param2[1:n_masked]
  param_input[ind_unmasked] <- param # parametres originals (input)

  n <- length(param_input)/6 # v+1
  param_models <- c(n+1) # k+1
  for (i in 1:n) { param_models=c(param_models,einasto.model,5) } # c(k+1,profile,v)
  param_models <- c(param_models,const.model,0) # c(k+1,profile,v,const profile, 0)
  alert <- mixture.model(param_input,param_models,output,action.control,arg)
  if(alert == FALSE) {
    out <- mixture.model(param=param_input,param2=param_models,output,action.multi,arg)
  }
  else {
    out <- output
  }
  out
}

# 
# mask.adjust_positions (creates an array of 1s and 0s for freezeing and thawing parameters in mask.freeze)
# 
# The input is an integer "k" corresponding to the number of components in the finite mixture model described by the 
# "param2" variable c(<k+1>,model.einasto,v, ... , model.einasto,v,model.const,0). And the input "v" is the number
# of parameters per Einasto profile including the mixture coefficient (v=6) 
# The output is an array of 1s and 0s with length v*(k-1). Position and mixing parameters are thawed.
# This function is a 3D version of function mask.adjust_positions from Kuhn et al. 2014.
# Usage> param <- mask.adjust_positions(8,6)
# 
mask.adjust_positions <- function(k,v) {

  mask <- rep(0,(v+1)*k)
  mask[1:3] <- 1
  if (k > 1) {
     indx <- (1:(k-1))*(v+1)
     indy <- 1+(1:(k-1))*(v+1)
     indz <- 2+(1:(k-1))*(v+1)
     indmix <- v+(1:(k-1))*(v+1)
     mask[indx] <- 1
     mask[indy] <- 1
     mask[indz] <- 1
     mask[indmix] <- 1
  }
  mask[(v+1)*k] <- 1
  mask
}

#
# mask.adjust_radii (creates an array of 1s and 0s for freezeing and thawing parameters in mask.freeze)
# 
# The input is an integer "k" corresponding to the number of components in the finite mixture model described by the 
# "param2" variable c(<k+1>,model.einasto,v, ... , model.einasto,v,model.const,0). And the input "v" is the number
# of parameters per Einasto profile including the mixture coefficient (v=6) 
# The output is an array of 1s and 0s with length v*(k-1). Radii parameters are thawed.
# This function is a 3D version of function mask.adjust_radii from Kuhn et al. 2014.
# Usage> param <- mask.adjust_radii(8,6)
# 
mask.adjust_radii <- function(k,v) {

  mask <- rep(0,(v+1)*k)
  mask[4] <- 1
  if (k > 1) {
     indr <- 3+(1:(k-1))*(v+1)
     mask[indr] <- 1
  }
  mask
}

#
# mask.adjust_sersic (creates an array of 1s and 0s for freezeing and thawing parameters in mask.freeze)
# 
# The input is an integer "k" corresponding to the number of components in the finite mixture model described by the 
# "param2" variable c(<k+1>,model.einasto,v, ... , model.einasto,v,model.const,0). And the input "v" is the number
# of parameters per Einasto profile including the mixture coefficient (v=6) 
# The output is an array of 1s and 0s with length v*(k-1). Sérisc index parameters are thawed.
# This function is a 3D version of function mask.adjust_sersic from Kuhn et al. 2014.
# Usage> param <- mask.adjust_sersic(8,6)
# 
mask.adjust_sersic <- function(k,v) {
  
  mask <- rep(0,(v+1)*k)
  mask[5] <- 1
  if (k > 1) {
     inds <- 4+(1:(k-1))*(v+1)
     mask[inds] <- 1
  }
  mask
}

#
# mask.adjust_mix (creates an array of 1s and 0s for freezeing and thawing parameters in mask.freeze)
# 
# The input is an integer "k" corresponding to the number of components in the finite mixture model described by the 
# "param2" variable c(<k+1>,model.einasto,v, ... , model.einasto,v,model.const,0). And the input "v" is the number
# of parameters per Einasto profile including the mixture coefficient (v=6) 
# The output is an array of 1s and 0s with length v*(k-1). Mixture coefficients are thawed.
# This function is a 3D version of function mask.adjust_mix from Kuhn et al. 2014.
# Usage> param <- mask.adjust_mix(8,6)
# 
mask.adjust_mix <- function(k,v) {

  mask <- rep(0,(v+1)*k)
  # mask[6] <- 1
  if (k > 1) {
     inds <- 5+(1:(k-1))*(v+1)
     mask[inds] <- 1
  }
  mask[(v+1)*k] <- 1
  mask
}

# mask.adjust_shape (creats an array of 1s and 0s for freezeing and thawing parameters in mask.freeze)
#
# The input is an integer "n" corresponding to the number of ellipsoids in a finite mixture model described by the 
# "param2" variable c(<n+1>,model.ell,5, ... , model.ell,5,model.const,0).
# The output is an array of 1s and 0s with length 6*n. Core radius and axis ratio parameters are thawed.
# Usage> param <- mask.freeze(param,mask.adjust_positions(n),clust=clust)
# 
mask.adjust_shape <- function(k,v) {
  mask <- rep(0,(v+1)*k)
  mask[4] <- 1
  mask[5] <- 1
  if (k > 1) {
    indsre <- 3+(1:(k-1))*(v+1)
    indsn <- 4+(1:(k-1))*(v+1)
    indsmix <- 5+(1:(k-1))*(v+1)
    mask[indsre] <- 1
    mask[indsn] <- 1
    mask[indsmix] <- 1
  }
  mask
}

#
# AIC.BIC (Calculates the AIC and BIC for a given set of parameters)
# 
# Inputs are the model parameters (param), the model functions (param2), the data set (clust),
# and the grid (quad).
# Output is the Maximum Log-Likelihod, the Akaike Information Criteron and the Bayesian Information Criterion
# Usage> AIC.BIC(param, param2, clust, quad)
#
AIC.BIC <- function(param,param2,clust,quad)
{
  N <- length(clust$data$x)
  L <- -model.lik(param, param2, mixture.model, clust, quad, all=TRUE);
  AIC <- -2*L + length(param)*2.0; 
  BIC <- -2*L + log(N)*length(param)
  print(paste("Maximum log-likelihood:", L))
  print(paste("AIC:", AIC))
  print(paste("BIC:", BIC))
  return(c(L,AIC,BIC));
}

#
# goodness (Evaluates the quality of a model)
# 
# The imputs are the data set (clust), 
# the grid where the residuals are evaluated (quad.residuals), the number of tiles per axis in the 
# grid (ntile),  the type of residuals (raw, inverse or Pearson; default is raw), the bandwidth for 
# the smoothing (bandwidth, default is 0.7), a distance for truncation of high density regions 
# (epsilon, default is 0.01), the functions that must be plotted (maps; s: absolute errors, 
# Output is a float between 0 and 1. 1 is perfect model fitting and 0 is no correlation between 
# model and data.
# Usage> goodness(fields, fast=TRUE)
#
goodness <- function(param, param2, clust, quad.residuals, ntile = c(32,32,32), residual.type="raw", bandwidth=1) 
{
  N <- length(clust$data$x);
  M <- length(quad.residuals$data$x);
  index <- 1:(N+M)
  window <- quad.residuals$domain
  minx <- window$xrange[1]; maxx <- window$xrange[2]; #epsx <- (maxx-minx)/ntile[1];
  miny <- window$yrange[1]; maxy <- window$yrange[2]; #epsy <- (maxy-miny)/ntile[2];
  minz <- window$zrange[1]; maxz <- window$zrange[2]; #epsz <- (maxz-minz)/ntile[3];
  
  unionquad <- rbind(cbind(clust$data$x, clust$data$y, clust$data$z), cbind(quad.residuals$data$x, quad.residuals$data$y, quad.residuals$data$z))
  unionquad <- cbind(index, unionquad)
  options(scipen = 20)
  write.table(unionquad, file="unionquad.txt", col.names = FALSE, row.names = FALSE)
  command <- paste("voro++ -o", minx, maxx, miny, maxy, minz, maxz, "unionquad.txt", sep=" ")
  system(command)
  
  # x y z residual voronoi lambda
  residuals <- residuals.lambda(N, M, param, param2, quad.residuals, residual.type, truncate)
  
  # Data density
  if(!file.exists("data_smoothed.txt")) {
    smooth.data <- centers(clust, bandwidth, ntile, 1)[[2]]
    write.table(smooth.data, file="data_smoothed.txt", col.names = FALSE, row.names = FALSE)
  }
  smooth.data <- read.table("data_smoothed.txt") # Sigma^*(r)
  
  # Model density
  if(!file.exists("dummy_smoothed.txt")) {
    command <- paste("./Code/kernel_lambda -i residuals_lambda.txt -n", toString(N), "-m", toString(M), "-g", toString(ntile[1]), toString(ntile[2]), toString(ntile[3]), "-w", toString(minx), toString(maxx), toString(miny), toString(maxy), toString(minz), toString(maxz), "-b", toString(bandwidth), "-r dummy_smoothed.txt")
    system(command)
  }
  smooth.dummy <- read.table("dummy_smoothed.txt") # Sigma^dag(r)
  
  fit <- summary(lm(smooth.dummy[,4] ~ smooth.data[,4] - 1))
  R2 <- fit$r.squared
  return(R2)
}

#
# residuals.lambda (Calculates the residuals at every grid location and data set point)
#
# The inputs are the number of data set points (N), the number of grid points (M), the model parameters
# (param), the model functions (param2), the grid (quad.residuals), the type of residuals (residual.type)
# and a distance for truncation of high density regions (epsilon).
# Output is a matrix containing in its first r columns, the locations of every data or grid point,
# the residuals as evaluated on each point, the weight corresponding to each point in a voronoi tesselation, 
# and the model intensity. This matrix is stored in file "residuals_lambda.txt".
# This function is a 3D version of function mpl.prepare from spatstat package.
# This function is called from residuals.mm, not by the user.
#
residuals.lambda <- function(N,M,param,param2,quad.residuals,residual.type="raw", truncate=0)
{
  # Prepare data variables
  voro <- read.table("unionquad.txt.vol") # Quadrature and weights
  M <- dim(voro)[1] - N
  weights <- voro[,5]
  Z <- c(rep(1,N),rep(0,M)) 
  Y <- rep(0,N+M) 
  Y[1:N] <- 1/weights[1:N]
  subset <- rep(TRUE,N+M) 
  marks <- rep(0,N+M) 
  
  #Evaluate model
  arg <- list(x1=voro[,2], y1=voro[,3], z1=voro[,4])
  output <- c()
  A <- mixture.model(param=param,param2=param2,output,action.multi, arg)
  
  subset[A==0] <- FALSE
  B=data.frame(A)
  colnames(B) <- c("B")
  trend <- ~offset(log(B))   #.mpl.Y ~ slope  ??? 
  spv <- package_version(versionstring.spatstat())
  the.version <- list(major = spv$major, minor = spv$minor, release = spv$patchlevel, date = "$Date: 2014/01/24 05:43:38 $")
  
  interaction = NULL; covariates = NULL; eps=NULL; computed=NULL;
  covfunargs = list(); correction = "border"; rbord = 0; use.gam = FALSE; 
  gcontrol = list(); famille = NULL; forcefit = FALSE; nd = NULL; 
  allcovar = FALSE; callstring = ""; precomputed = NULL; 
  savecomputed = FALSE; preponly = FALSE; rename.intercept = TRUE; 
  justQ = FALSE; weightfactor = NULL; vnamebase = c("Interaction", "Interact."); 
  vnameprefix = NULL; warn.illegal = TRUE; warn.unidentifiable = TRUE; 
  weightfactor = NULL; skip.border = FALSE
  
  Q <- voro[,2:4]
  want.trend <- !is.null(trend)  # && !identical.formulae(trend,~1)
  want.inter <- !is.null(interaction) && !is.null(interaction$family)
  trend.formula <- trend
  computed <- list()
  problems <- list()
  names.precomputed <- names(precomputed)
  likelihood.is.zero <- FALSE
  is.identifiable <- TRUE
  
  .mpl <- list(W = weights, Z = Z, Y = Y, MARKS = marks, SUBSET = subset)
  
  glmdata <- data.frame(.mpl.W = .mpl$W, .mpl.Y = .mpl$Y)
  izdat <- .mpl$Z[.mpl$SUBSET]
  ndata <- sum(izdat)
  ndummy <- sum(!izdat)
  
  skip.border <- skip.border && (correction == "border")
  internal.names <- c(".mpl.W", ".mpl.Y", ".mpl.Z", ".mpl.SUBSET", "SUBSET", ".mpl")
  reserved.names <- c("x", "y", "z", "marks", internal.names)
  
  check.clashes <- function(forbidden, offered, where) {
    name.match <- outer(forbidden, offered, "==")
    if (any(name.match)) {
      is.matched <- apply(name.match, 2, any)
      matched.names <- (offered)[is.matched]
      if (sum(is.matched) == 1) {
        return(paste("The variable", sQuote(matched.names), 
                     "in", where, "is a reserved name"))
      }
      else {
        return(paste("The variables", paste(sQuote(matched.names), 
                                            collapse = ", "), "in", where, "are reserved names"))
      }
    }
    return("")
  }
  
  trendvariables <- variablesinformula(trend)
  cc <- check.clashes(internal.names, trendvariables, "the model formula")
  cc <- check.clashes(reserved.names, names(covariates), sQuote("covariates"))
  covariates.df <- B
  needed <- names(covariates.df) %in% trendvariables #Revisar
  covariates.needed <- covariates.df[, needed, drop = FALSE]
  glmdata <- data.frame(glmdata, covariates.needed)
  nbg <- is.na(covariates.needed)
  
  Vnames <- NULL
  IsOffset <- NULL
  
  glmdata <- cbind(glmdata, data.frame(.mpl.SUBSET = .mpl$SUBSET, stringsAsFactors = FALSE))
  trendpart <- paste(as.character(trend), collapse = " ")
  rhs <- trendpart
  fmla <- paste(".mpl.Y ", rhs)
  fmla <- as.formula(fmla)
  trendfmla <- paste(".mpl.Y ", trendpart)
  prep <- list(fmla = fmla, trendfmla = trendfmla, covariates = B, glmdata = glmdata, 
               Vnames = Vnames, IsOffset = IsOffset, problems = problems, 
               likelihood.is.zero = likelihood.is.zero, is.identifiable = is.identifiable, 
               computed = computed)
  
  fmla <- prep$fmla
  glmdata <- prep$glmdata
  problems <- prep$problems
  likelihood.is.zero <- prep$likelihood.is.zero
  is.identifiable <- prep$is.identifiable
  computed <- append(computed, prep$computed)
  IsOffset <- prep$IsOffset
  .mpl.W <- glmdata$.mpl.W
  .mpl.SUBSET <- glmdata$.mpl.SUBSET
  
  gc <- "glm.control"
  gcontrol <- list()
  gcontrol <- do.call(gc, gcontrol)
  FIT <- glm(fmla, family = quasi(link = log, variance = mu), 
             weights = .mpl.W, data = glmdata, subset = .mpl.SUBSET, 
             control = gcontrol, model = FALSE)
  
  environment(FIT$terms) <- sys.frame(sys.nframe())
  co <- FIT$coef
  W <- glmdata$.mpl.W
  SUBSET <- glmdata$.mpl.SUBSET
  Vnames <- prep$Vnames
  
  fitin <- fii()
  rslt <- list(method = "mpl", fitter = if (use.gam) "gam" else "glm", 
               projected = FALSE, coef = co, trend = if (want.trend) trend else NULL, 
               interaction = if (want.inter) interaction else NULL, 
               fitin = fitin, Q = Q, maxlogpl = 0, internal = list(glmfit = FIT, 
                                                                   glmdata = glmdata, Vnames = Vnames, IsOffset = IsOffset, 
                                                                   fmla = fmla, computed = computed), covariates = covariates, 
               covfunargs = covfunargs, correction = correction, rbord = rbord, 
               terms = terms(trend.formula), version = the.version, 
               problems = problems)
  rslt2 <- rslt
  class(rslt2) <- "ppm"
  
  #Lambda
  lambda <- GLMpredict(rslt$internal$glmfit, rslt$internal$glmdata, rslt$coef, changecoef = FALSE)
  
  Z <- c(rep(1,N),rep(0,M))
  if(residual.type == "raw") {residuals <- Z - weights*lambda}  #Raw residuals
  if(residual.type == "inverse") {residuals <- (Z/lambda) - weights*rep(-1,N+M)}  #Inverse residuals
  if(residual.type == "pearson") {pearson <- Z/sqrt(lambda) - weights*sqrt(lambda)}  #Pearson residuals
  
  norm <- int.model(param, param2, quad.residuals, mixture.model)   
  sigma <- N*A/norm
  R <- cbind(voro[,2:4],residuals,weights,sigma) 
  write.table(R,file="residuals_lambda.txt",row.names=FALSE,col.names=FALSE)
  return(R)
}

#
# residuals.mm (Produces the residuals of an estimated model)
# 
# The imputs are the model parameters (param), the model functions (param2), the data set (clust),
# the grid where the residuals are evaluated (quad.residuals), the number of tiles per axis in the 
# grid (ntile), the type of residuals (raw, inverse or Pearson), a distance for truncation of high 
# density regions (epsilon) and a logical variable indicating if the Fast Fourier Transformation 
# is to be used (fast).
# No output is produced. This function is called by functions plot.mm3d and goodness.
# Usage> residuals.mm(param, param2, clust, quad.residuals, ntile, residual.type, epsilon, fast
#
residuals.mm <- function(param, param2, clust, quad.residuals, ntile, residual.type="raw", truncate=0,fast=TRUE)
{
  # Voronoi tesselation
  N <- length(clust$data$x);
  M <- length(quad.residuals$data$x);
  index <- 1:(N+M)
  window <- quad.residuals$domain
  minx <- window$xrange[1]; maxx <- window$xrange[2]; epsx <- (maxx-minx)/ntile[1];
  miny <- window$yrange[1]; maxy <- window$yrange[2]; epsy <- (maxy-miny)/ntile[2];
  minz <- window$zrange[1]; maxz <- window$zrange[2]; epsz <- (maxz-minz)/ntile[3];
  
  unionquad <- rbind(cbind(clust$data$x, clust$data$y, clust$data$z), cbind(quad.residuals$data$x, quad.residuals$data$y, quad.residuals$data$z))
  unionquad <- cbind(index, unionquad)
  options(scipen = 20)
  write.table(unionquad, file="unionquad.txt", col.names = FALSE, row.names = FALSE)
  command <- paste("voro++ -o", minx, maxx, miny, maxy, minz, maxz, "unionquad.txt", sep=" ")
  system(command)
  
  #if(!file.exists("unionquad.txt.vol")) {
  #  print("voro++ failed. Starting alternative.")
  #  write.table(c(N,ntile,minx, maxx, miny, maxy, minz, maxz), file="dimensions.txt",col.names = FALSE, row.names = FALSE)
  #  command <- "./Code/quadscheme"
  #  system(command)
  #  command <- "cat datacat.txt quadrature.txt > unionquad.txt"
  #  system(command)
  #  command <- paste("voro++ -o", minx, maxx, miny, maxy, minz, maxz, "unionquad.txt", sep=" ")
  #  system(command)
  #}
  
  # x y z residual voronoi lambda
  residuals <- residuals.lambda(N, M, param, param2, quad.residuals, residual.type, truncate,fast)
}

#
# ready.mcmc (Prepares the necessary inputs for the MCMC routine)
# 
# Inputs are the spatial point pattern as a pp3 object (clust), the quadrature (quad),
# the finite mixture model expressed as list (param2) and a distance (rbord) for edge corrections.
# The output is a list containing the necessary information to start the MCMC routine.
# Usage> MyData <- ready.mcmc(clust,quad,param2,rbord=0.5)
#
ready.mcmc <- function(clust=clust,quad=quad,param2,rbord=rbord)
{
	#Laplaces parameters
	mon.names="LP"
	k <- param2[[1]]-1
	parm.list <- list()
	nam <- paste("beta", 1, sep = "")
	parm.list[[1]] <-  rep(0,3); names(parm.list)[1] <- "beta1";
	parm.list[[2]] <-  0; names(parm.list)[2] <- "re1";
	parm.list[[3]] <-  0; names(parm.list)[3] <- "n1";
	for(i in 1:(k-1))
	{
		ini <- i*4;
		parm.list[[ini]] <-  rep(0,3); names(parm.list)[ini] <- paste("beta",i+1,sep="");
		parm.list[[ini+1]] <-  0; names(parm.list)[ini+1] <- paste("re",i+1,sep="");
		parm.list[[ini+2]] <-  0; names(parm.list)[ini+2] <- paste("n",i+1,sep="");
		parm.list[[ini+3]] <-  0; names(parm.list)[ini+3] <- paste("mix",i+1,sep="");
	}
	parm.list[[4*(k-1)+4]] <- 0; names(parm.list)[4*(k-1)+4] <- "back"
	parm.names <- as.parm.names(parm.list)
	
	win=box3(xrange=quad$domain$xrange,yrange=quad$domain$yrange,zrange=quad$domain$zrange)
	N <- length(clust$data$x)
	MyData <- list(clust=clust, mon.names=mon.names, parm.names=parm.names, param2=param2, win=win, rbord=rbord, quad=quad, N=N)
	return(MyData)
}

#
# run.ld (Executes the MCMC LaplacesDemon routine)
# 
# Inputs are the initial parameters of the Mixture Model, the necessary information as generated
# by the ready.mcmc() function, the MCMC method and the number of iterations. This function is 
# prepared to use either method "AMWG" and "twalk". For more information please visit:
# https://www.rdocumentation.org/packages/LaplacesDemon/versions/16.1.1/topics/LaplacesDemon
# Usage> Fit <- run.ld(Initial.Values, MyData, "AMWG", 1000)
run.ld <- function(Initial.Values,MyData,method,iter)
{
	set.seed(666)
	k <- MyData$param2[[1]]
	v <- MyData$param2[[3]]
	if(method=="AMWG")
	{
		bloq <- Blocks(Initial.Values,N=k)
		bloq[[1]] <- c(1:5)
		for(i in 2:k)
		{
			ini <- (i-1)*(v+1); las <- ini+v;
			bloq[[i]] <- c(ini:las)
		}
		bloq[[k]] <- c(las+1)
		Fit <- LaplacesDemon(model.post, Data=MyData, Initial.Values, Covar=NULL, Iterations=iter, Status=100, Thinning=2, Algorithm="AMWG", Specs=list(B=NULL,n=0, Periodicity=50)) # B=bloq
	}
	if(method=="twalk")
	{
		l <- length(Initial.Values)
		altri <- Initial.Values + rnorm(l,0,0.1)
		Fit <- LaplacesDemon(model.post, Data=MyData, Initial.Values, Covar=NULL, Iterations=iter, Status=100, Thinning=2, Algorithm="twalk", Specs=list(SIV=altri, n1=4, at=6, aw=1.5))		
	}
	return(Fit)
}

#
# model.post (Returns the log Posterior)
# 
# The inputs are the "param" vector (described in mixture.model) and list MyData with all 
# necessary information.
# The output is the log posterior probability of the model.
# Usage> logpost <- model.post(param, MyData)
#
model.post <- function(parm, MyData)
{
  output=NULL
  arg = MyData$win
  alert <- mixture.model(parm,MyData$param2,output,action.control,arg)
	if(alert==TRUE){LP=-Inf; LL=-Inf;}
	else
	{
		#Log (Prior Densities)
		lprior <- mixture.model(parm, MyData$param2, output=c(), action.prior, arg) #.prior(parm,MyData$param2)
		#Normalization cnt (for the likelihood)
		norm <- int.model(parm,MyData$param2,MyData$quad, mixture.model)
		#Log-likelihood
		LL <- sum(log(MyData$N*mixture.model(param=parm,param2=MyData$param2,output=c(), action.multi, MyData$clust$data)/norm))-MyData$N
		#Log-posterior
		LP <- LL + sum(lprior)
	}
	Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=parm, parm=parm)
  return(Modelout)
}

#
# outside.win (Checks if a component center is inside the window)
# 
# Inputs are the parameters of a halo component ("param_new") and the 
# window "win" of the spatial point pattern
# The output is a boolean avariable being TRUE when the center of a components is outside the window
# Usage> outside.win(c(0,0,0,3,5,0.7),win)
#
outside.win <- function(param_new,win) {
	outside <- FALSE
	if(!is.na(param_new[1])){
	  if(param_new[1] < win$xrange[1]) {outside=TRUE}
	  if(param_new[1] > win$xrange[2]) {outside=TRUE}
	}
	if(!is.na(param_new[2])){
    if(param_new[2] < win$yrange[1]) {outside=TRUE}
	  if(param_new[2] > win$yrange[2]) {outside=TRUE}
	}
	if(!is.na(param_new[3])){
    if(param_new[3] < win$zrange[1]) {outside=TRUE}
    if(param_new[3] > win$zrange[2]) {outside=TRUE}
	}
	return(outside)
}

#
# image.3D (Sums all values with two common coordinates and creates an image)
# 
# The inputs are a matrix with four columns (R), the direction in wich the values are summed (proj.var)
# and the data set window (window). Matrix R contains in its first three columns a grid as generated
# by function quad.3d. The fourth column are the evaluation of a funciton on every grid location
# (data density, model density or residuals). This function sums the fourth column for every location 
# with two coincidental dimensions. If proj.var="Z", then all values with a common X and Y coordinates
# are summed and a two dimensional image is produced.
# Output is an image shaped like one of the faces of the data set window and the matriz used to produce it.
# This function is called by function plot.dens().
#
image.3D <- function(R, proj.var="Z", window=clust$domain)
{
  x <- sort(unique(R[,1])); #x <- x[-c(1,length(x))];
  y <- sort(unique(R[,2])); #y <- y[-c(1,length(y))];
  if(dim(R)[2] > 2) {z <- sort(unique(R[,3]));} #z <- z[-c(1,length(z))];
	
	# X Y
	if(proj.var=="Z") {
	  R.ima <- matrix(NA,length(y),length(x)); 
	  win.ima <- owin(window$xrange,window$yrange);  
	  k <- 1 
	  for(i in x) 
	  {
		  pla <- R[R[,1]==i,] 
		  l <- 1 
		  for(j in y) 
		  {
			  recta <- pla[pla[,2]==j,]
			  R.ima[l,k] <- sum(recta[,4])
			  l <- l + 1
		  }
		  k <- k + 1
	  }
	  ima <- as.im(R.ima,W=win.ima)		
	}	
	
	# X Z
	if(proj.var=="Y") {
	  R.ima <- matrix(NA,length(z),length(x)); 
	  win.ima <- owin(window$xrange,window$zrange); 
	  k <- 1
	  for(i in x)
	  {
		  pla <- R[R[,1]==i,]
		  l <- 1
		  for(j in z)
		  {
			  recta <- pla[pla[,3]==j,]
			  R.ima[l,k] <- sum(recta[,4])
			  l <- l + 1
	  	}
		  k <- k + 1
	  }
	  ima <- as.im(R.ima,W=win.ima)
	}
	
	# Y Z
	if(proj.var=="X") {
	  R.ima <- matrix(NA,length(z),length(y)); 
	  win.ima <- owin(window$yrange,window$zrange); 
	  k <- 1
  	for(i in y)
	  {
		  pla <- R[R[,2]==i,]
		  l <- 1
		  for(j in z)
		  {
			  recta <- pla[pla[,3]==j,]
			  R.ima[l,k] <- sum(recta[,4])
			  l <- l + 1
		  }
		  k <- k + 1
	  }
	  ima <- as.im(R.ima,W=win.ima)
	}
	
	return(list(ima,R.ima))
}

# 
# ds9.colors (creates an intensity colormap mimicking the ds9 hsv scheme)
# 
# The input "n" is the number of steps in the map.
# The output is an array of hexadecimal RGB colors.
# This function is obtained from Kuhn et al. 2014.
# Usage> image(t(volcano)[ncol(volcano):1,],col=ds9.colors(512))
# 
ds9.colors <- function(n) {
 x <- (0:n)*1.0/n
 h <- (270.0/360.0-x+1.0) %% 1.0
 s <- 4.0*x*(1.0-x)
 v <- x^0.33
 color <- hsv(h=h,s=s,v=v)
 color
}

# 
# american.colors (creates an intensity colormap with white as the middle color)
# 
# The input "n" is the number of steps in the map.
# The output is an array of hexadecimal RGB colors.
# This function is obtained from Kuhn et al. 2014.
# Usage> image(t(volcano)[ncol(volcano):1,],col=american.colors(512))
# 
american.colors <- function(n) {
 x <- (0:n)*1.0/n
 h <- round(1-x)*0.667
 v <- pmin(1.0,(1.1-abs(1-2*x))*3)
 s <- abs(2*x-1)
 hsv(h=h,s=s,v=v)
}

#
# plot.dens (plots images of a function evaluated on a grid)
# 
# The inputs are the matrix with four columns (R) as explained in function image.3D, 
# the data set (clust), the color of the image scale (col.plot, gray or color), the color for 
# the data set points (col.point), the scale of the image scale (scale: linear (lin) or logarithmic (log)),
# the direction in which the matrix R will be summed (proj.var: X, Y, Z), the data set window (window),
# a logical variable stating if the image color scale will be symmetric around zero (symmetry), a logical
# variable stating if a figure has to be stored (print), the dimensions of the stored figure (w, h) and
# its name (name).
# There is no output. This function plots an image or stores it.
# Usage> plot.dens(R,clust, col.plot="color",scale="lin",proj.var="X", window, symmetry=TRUE, w=880, h=880, name="figure")
#
plot.dens <- function(R,clust=clust,col.plot="gray",col.point = 1, scale="lin",proj.var="Z", window=clust$domain, symmetry=TRUE, print=FALSE, print.data=FALSE, title='Data', w = 880, h = 880, name="density.png")
{
	#X <- clust$data
	image <- image.3D(R,proj.var, window)
	ima <- image[[1]]
#	ima <- ima / max(ima)
	R.ima <- image[[2]]
#	R.ima <- R.ima / max(R.ima)

	if(proj.var=="Z") {X <- as.matrix(cbind(clust$data$x, clust$data$y))}
	if(proj.var=="Y") {X <- as.matrix(cbind(clust$data$x, clust$data$z))}
	if(proj.var=="X") {X <- as.matrix(cbind(clust$data$y, clust$data$z))}
	
	if(col.plot=="gray") {
		delta <- (log10(1) - log10(0.001))/100
		escala_log <- 10^(log10(0.001) + (c(100:1)-0.5)*delta) 
		escala <- c(100:0)*0.01 
		if(scale=="lin") {col.scale <- gray(escala)}
		if(scale=="log") {col.scale <- gray(escala_log)}
	}
	if(col.plot=="color") {
	  col.scale <- american.colors(512)
	}
	if(symmetry==TRUE) {
	  ddmax <- max(c(max(R.ima),-min(R.ima))); 
	  ddmin <- -ddmax;
	} else {
	  ddmax <- max(R.ima);
	  ddmin <- min(R.ima);
	}
	if(print==TRUE) {png(name, width = w, height = h)}
	plot(ima,main='',box=TRUE,col=col.scale,zlim=c(ddmin,ddmax),ribsep=0.02,ribside='right',ribargs=,cex.axis=2,axes=TRUE)
	title(title,cex.main=3)
	# mtext(text = "Particle density", side = 4, line = 1, cex=2)
	#,xlab=expression(paste("h"^-1, " Mpc", sep="")),ylab=expression(paste("h"^-1, " Mpc", sep="")))
	mtext(expression(paste("h"^-1, " Mpc", sep="")), side=1, line=1.5, at=31, cex=2)
	mtext(expression(paste("h"^-1, " Mpc", sep="")), side=2, line=0.5, at=184, cex=2)
	if(print.data==TRUE) {points(X[,1], X[,2], pch=20,cex=0.2, col=col.point)}
	if(print==TRUE) {dev.off();}
}

#
# plot.mm3d (Plot different figures of our model)
#
# The imputs are the model parameters (param), the model functions (param2), the data set (clust), 
# the grid where the residuals are evaluated (quad.residuals), the number of tiles per axis in the 
# grid (ntile),  the type of residuals (raw, inverse or Pearson; default is raw), the bandwidth for 
# the smoothing (bandwidth, default is 1), a distance for truncation of high density regions 
# (epsilon, default is 0.01), the functions that must be plotted (maps; s: absolute errors, 
# e: relative errors, d: data density, m: model density), a logical variable indicating if the Fast 
# Fourier Transformation is to be used (fast), the points color (col.point), the chromatic scale 
# (linear: "lin", default; or logarithmic: "log"), the projection direction (proj.var, default is Z)
# a logical variable stating if the plots is to be stored as a figure (print) and its dimensions (w, h).
# No output is produced, different plots are made and stored when indicated.
# Usage> plot.mm3d(param, param2, clust, quad.residuals, "raw", bandwidth=1, epsilon=0.01, maps=c("s", "e", "d", "m"), fast, col.point=1, scale="lin", proj.var="Z", print=FALSE, w = 880, h = 880)
#
plot.mm3d <- function(param, param2, clust, quad.residuals, residual.type="raw", bandwidth=1, maps = c("s", "e", "d", "m"), fast=TRUE, col.point=1, scale="lin", proj.var="Z", print=TRUE, print.data=FALSE, w = 880, h = 880, ntile) {# ntile
  
  # Voronoi tesselation
  N <- length(clust$data$x);
  M <- length(quad.residuals$data$x);
  index <- 1:(N+M)
  window <- quad.residuals$domain
  minx <- window$xrange[1]; maxx <- window$xrange[2]; 
  miny <- window$yrange[1]; maxy <- window$yrange[2]; 
  minz <- window$zrange[1]; maxz <- window$zrange[2]; 
  
  unionquad <- rbind(cbind(clust$data$x, clust$data$y, clust$data$z), cbind(quad.residuals$data$x, quad.residuals$data$y, quad.residuals$data$z))
  unionquad <- cbind(index, unionquad)
  options(scipen = 20)
  write.table(unionquad, file="unionquad.txt", col.names = FALSE, row.names = FALSE)
  command <- paste("voro++ -o", minx, maxx, miny, maxy, minz, maxz, "unionquad.txt", sep=" ")
  system(command)

  # x y z residual voronoi lambda
  residuals <- residuals.lambda(N, M, param, param2, quad.residuals, residual.type, truncate)
  
  # Absolute errors
  if("s" %in% maps) {
    if(!file.exists("residuals_smoothed.txt")) {
      # Compute
      if(fast==TRUE) {
        command <- paste("./Code/kernel_fourier_residuals -i residuals_lambda.txt -n", toString(N), "-m", toString(M), "-g", toString(ntile[1]), toString(ntile[2]), toString(ntile[3]), "-w", toString(minx), toString(maxx), toString(miny), toString(maxy), toString(minz), toString(maxz), "-b", toString(bandwidth), "-r residuals_smoothed.txt")
      } else {
        command <- paste("./Code/kernel_absolute_residuals -i residuals_lambda.txt -n", toString(N), "-m", toString(M), "-g", toString(ntile[1]), toString(ntile[2]), toString(ntile[3]), "-w", toString(minx), toString(maxx), toString(miny), toString(maxy), toString(minz), toString(maxz), "-b", toString(bandwidth), "-r residuals_smoothed.txt")
      }
      system(command)
    }
    # Plot
    smooth.absolute <- read.table("residuals_smoothed.txt") # s(r)
    name <- "absolute_errors3d.png"
    plot.dens(smooth.absolute,clust, col.plot="color", col.point,scale,proj.var, window, symmetry=TRUE, print, print.data, title='Smoothed raw residuals', w, h, name) # Abs errors
  }

  # Relative errors
  if("e" %in% maps) {
    if(!file.exists("residuals_smoothed.txt")) {
      # Compute
      if(fast==TRUE) {
        command <- paste("./Code/kernel_fourier_residuals -i residuals_lambda.txt -n", toString(N), "-m", toString(M), "-g", toString(ntile[1]), toString(ntile[2]), toString(ntile[3]), "-w", toString(minx), toString(maxx), toString(miny), toString(maxy), toString(minz), toString(maxz), "-b", toString(bandwidth), "-r residuals_smoothed.txt")
      } else {
        command <- paste("./Code/kernel_absolute_residuals -i residuals_lambda.txt -n", toString(N), "-m", toString(M), "-g", toString(ntile[1]), toString(ntile[2]), toString(ntile[3]), "-w", toString(minx), toString(maxx), toString(miny), toString(maxy), toString(minz), toString(maxz), "-b", toString(bandwidth), "-r residuals_smoothed.txt")
      }
      system(command)
    }
    if(!file.exists("dummy_smoothed.txt")) {
      command <- paste("./Code/kernel_lambda -i residuals_lambda.txt -n", toString(N), "-m", toString(M), "-g", toString(ntile[1]), toString(ntile[2]), toString(ntile[3]), "-w", toString(minx), toString(maxx), toString(miny), toString(maxy), toString(minz), toString(maxz), "-b", toString(bandwidth), "-r dummy_smoothed.txt")
      system(command)
    }
    # Plot
    smooth.absolute <- read.table("residuals_smoothed.txt") 
    smooth.dummy <- read.table("dummy_smoothed.txt")
    smooth.relative <- smooth.absolute[,1:3]
    smooth.relative[,4] <- smooth.absolute[,4]/smooth.dummy[,4]
    write.table(smooth.relative, file="relative_smoothed.txt", col.names=FALSE, row.names=FALSE)
    name <- "relative_errors3d.png"
    plot.dens(smooth.relative,clust, col.plot="color", col.point, scale,proj.var, window, symmetry=TRUE, print, print.data, title='Relative residuals', w, h, name)
  }

  # Data density
  if("d" %in% maps) {
    if(!file.exists("data_smoothed.txt")) {    
      smooth.data <- centers(clust, bandwidth, ntile, 1)[[2]]
      write.table(smooth.data, file="data_smoothed.txt", col.names = FALSE, row.names = FALSE)
    }
    smooth.data <- read.table("data_smoothed.txt")
    name <- "data_density3d.png"
    plot.dens(smooth.data,clust, col.plot="gray", col.point, scale, proj.var, window, symmetry=FALSE, print, print.data, title = 'Data density distribution', w, h, name)
  }

  # Model density
  if("m" %in% maps) {
    if(!file.exists("dummy_smoothed.txt")) {
      command <- paste("./Code/kernel_lambda -i residuals_lambda.txt -n", toString(N), "-m", toString(M), "-g", toString(ntile[1]), toString(ntile[2]), toString(ntile[3]), "-w", toString(minx), toString(maxx), toString(miny), toString(maxy), toString(minz), toString(maxz), "-b", toString(bandwidth), "-r dummy_smoothed.txt")
      system(command)
    }
    smooth.dummy <- read.table("dummy_smoothed.txt")
    name <- "model_density3d.png"
    plot.dens(smooth.dummy, clust, col.plot="gray", col.point, scale, proj.var, window, symmetry=FALSE, print, print.data, title='Model density distribution', w, h, name)
  }
}

#
# plot.mm2d (Plot different figures of our model)
#
# The imputs are the model parameters (param), the model functions (param2), the data set (clust), 
# the grid where the residuals are evaluated (quad.residuals), the number of tiles per axis in the 
# grid (ntile), the bandwidth for the smoothing (bandwidth, default is 1), the functions that must be 
# plotted (maps; s: absolute errors, e: relative errors, d: data density, m: model density), the points 
# color (col.point), the chromatic scale (linear: "lin", default; or logarithmic: "log"), the projection 
# direction (proj.var, default is Z) a logical variable stating if the plots is to be stored as a figure 
# (print) and its dimensions (w, h). The density fields are collapsed along dimension proj.var and 2d 
# smoothings are done using spatstat functions.
# No output is produced, different plots are made and stored when indicated.
# Usage> plot.mm2d(param, param2, clust, quad.residuals, bandwidth=1, epsilon=0.01, maps=c("s", "e", "d", "m"), fast, col.point=1, scale="lin", proj.var="Z", print=FALSE, w = 880, h = 880)
#
plot.mm2d <- function(param, param2, clust, quad.residuals, bandwidth=1, maps = c("s", "e", "d", "m"), col.point=2, scale="lin", proj.var="Z", print=TRUE, print.data=FALSE, w = 880, h = 880) {
  
  x = clust$data$x; y = clust$data$y; z = clust$data$z;
  window <- clust$domain
  minx <- window$xrange[1]; maxx <- window$xrange[2]; 
  miny <- window$yrange[1]; maxy <- window$yrange[2]; 
  minz <- window$zrange[1]; maxz <- window$zrange[2]; 
  if(proj.var=="X") {data.ppp <- ppp(x=y, y=z, window=owin(c(miny,maxy),c(minz,maxz)))}
  if(proj.var=="Y") {data.ppp <- ppp(x=x, y=z, window=owin(c(minx,maxx),c(minz,maxz)))}
  if(proj.var=="Z") {data.ppp <- ppp(x=x, y=y, window=owin(c(minx,maxx),c(miny,maxy)))}
  X <- cbind(data.ppp$x, data.ppp$y)
  
  # Color scale
  delta <- (log10(1) - log10(0.001))/100
  escala_log <- 10^(log10(0.001) + (c(100:1)-0.5)*delta) 
  escala <- c(100:0)*0.01 
  if(scale=="lin") {col.scale <- gray(escala)}
  if(scale=="log") {col.scale <- gray(escala_log)}
  
  # Data density
  if("d" %in% maps) {
    data.dens <- density(data.ppp, sigma=bandwidth); 
    name <- "data_density2d.png"
    if(print==TRUE) {png(name, width = w, height = h)}
    par(mar = c(0,0,0,1))
    plot(data.dens,box=FALSE,main=NULL,col=col.scale,ribsep=0.02,ribside='right',ribargs=,cex.axis=2)
    if(print.data==TRUE) {points(x, y, pch=20,cex=0.2, col=col.point) }
    if(print==TRUE) {dev.off();}
  }
  
  # Model density
  if("m" %in% maps) {
    arg = list(x1=quad.residuals$data$x, y1=quad.residuals$data$y, z1=quad.residuals$data$z)
    A <- mixture.model(param,param2,output=c(),action.multi,arg)
    R <- cbind(arg$x1, arg$y1, arg$z1, A)
    image <- image.3D(R, proj.var, window)
    ima <- image[[1]]
    R.ima <- image[[2]]
    # Point pattern fit
    mlfit <- ppm(data.ppp, ~offset(log(B)), covariates = list(B=ima))
    model.dens <- diagnose.ppm(mlfit,which="marks",sigma=bandwidth,main='',plot.it=FALSE)  
    # Plot
    name <- "model_density2d.png"
    if(print==TRUE) {png(name, width = w, height = h)}
    par(mar = c(0,0,0,1))
    plot(-model.dens$Ydens,box=FALSE,main=NULL,col=col.scale,ribsep=0.02,ribside='right',ribargs=,cex.axis=2) 
    if(print.data==TRUE) {points(x, y, pch=20,cex=0.2, col=col.point) }
    if(print==TRUE) {dev.off();}
  }
  
  # Absolute errors
  if("s" %in% maps) {
    abs.dens <- diagnose.ppm(mlfit,which="smooth",sigma=bandwidth,plot.smooth='image',main='',plot.it=FALSE)
    ddmax <- max(c(max(abs.dens$smooth$Z),-min(abs.dens$smooth$Z)))
    name <- "absolute_errors2d.png"
    if(print==TRUE) {png(name, width = w, height = h)}
    par(mar = c(0,0,0,1))
    plot(abs.dens$smooth$Z,box=FALSE,main='',zlim=c(-ddmax,ddmax),col=american.colors(512),ribsep=0.02,ribside='right',ribargs=,cex.axis=2)
    if(print.data==TRUE) {points(x, y, pch=20,cex=0.2, col=col.point) }
    if(print==TRUE){dev.off();}
  }  
    
  # Relative errors
  if("e" %in% maps) {
    rel.dens <- as.im(-abs.dens$smooth$Z/model.dens$Ydens)
    ddmax <- max(c(max(rel.dens),-min(rel.dens)))
    name <- "relative_errors2d.png"
    if(print==TRUE) {png(name, width = w, height = h)}
    par(mar = c(0,0,0,1))
    plot(rel.dens,box=FALSE,main='',zlim=c(-ddmax,ddmax),col=american.colors(512),ribsep=0.02,ribside='right',ribargs=,cex.axis=2)
    if(print.data==TRUE) {points(x, y, pch=20,cex=0.2, col=col.point) }
    if(print==TRUE){dev.off();}
  }
}

#
# gen.pattern (Generates a realization of the process)
# 
# Inputs are the model parameters (param), the model functions (param2), the number of points per
# component (pop), the data set window (window) and an integer (scl). Argument scl is to be increased
# if the function do not generate as many points as sum(pop).
# Output is a pp3 object with the generated data set.
# Usage> gen.pattern(param, param2, pop, window, scl=2)
#
gen.pattern <- function(param, param2, pop, window, scl=2)
{
  # Fill empty values if existing
  df <- data.frame(matrix(c(1:(k+1), rep(0,k+1)), ncol = 2))
  colnames(df) <- c('Var1', 'Zero')
  colnames(pop) <- c('Var1', 'Freq')
  pop <- merge(df, pop, by='Var1', all=TRUE)
  pop[is.na(pop)] <- 0
  pop <- pop$Freq

  # Maximum distance inside window
  max.dist <- sqrt((window$xrange[2]-window$xrange[1])^2 + (window$yrange[2]-window$yrange[1])^2 + (window$zrange[2]-window$zrange[1])^2)
  k <- param2[[1]]-1
  v <- param2[[3]]
  param.pass <- param[1:v]
  param.pass <- c(param.pass, 0)
  
  # Parameters for generated values
  h <- 1e-3; xlim <- c(0,max.dist); d <- seq(0,max.dist,h)
  
  halo <- c(); 
  for(i in 1:k) {
    M <- round(pop[i])*scl
    # Generate distances
    pro <- profile.sphere(d,param.pass[4],param.pass[5]) # Einasto profile in 3D
    pron <- pro/sum(pro*h) # Normlized
    cum <- pron*0; cum[1] <- pron[1]*h # Create cummulative function
    for(j in 2:length(d)) { cum[j] <- cum[j-1] + pron[j]*h }
    dist <- rep(0,M); epsi <- runif(M,0,1)
    for(j in 1:M) {   # Find distance for epsi value
      ide <- max(which((epsi[j] > cum) == TRUE))
      dist[j] <- (d[ide+1] + d[ide])/2
    }
    # Generate halo
    theta <- acos(1-2*runif(M,0,1))
    phi <- 2*pi*runif(M,0,1)
    x <- dist*sin(theta)*cos(phi) + param.pass[1]
    y <- dist*sin(theta)*sin(phi) + param.pass[2]
    z <- dist*cos(theta) + param.pass[3]
    iner <- inside.shell(x,y,z,window)[[1]]
    min <- min(length(iner[,1]), round(pop[i]))
    if(min != 0) {
      iner <- iner[1:min,]
      halo <- rbind(halo, iner) # Save halo
    }
    # New halo
    p1 <- (v+1)*i
    p2 <- (v+1)*i+v
    param.pass <- param[p1:p2]
  } 
  # Background
  x <- runif(round(pop[k+1]), window$xrange[1], window$xrange[2])
  y <- runif(round(pop[k+1]), window$yrange[1], window$yrange[2])
  z <- runif(round(pop[k+1]), window$zrange[1], window$zrange[2])
  halo <- rbind(halo, cbind(x,y,z))  
  
  if(dim(halo)[1] < sum(pop)) {
    cat(dim(halo)[1])
    cat(" points generated. ")
    cat(sum(pop))
    cat(" expected. Increase 'scl' and repeat.")
  }
  
  new_clust <- pp3(x=halo[,1],y=halo[,2],z=halo[,3],domain=window)
  return(new_clust)
}

#
# profile.emp (Calculates the empyrical profile of the point process centred in a halo)
# 
# Imputs are the estimated parameters of a single model component (halo), the point process (clust) and 
# the radii of the shells where the profile is calculated (spheres).
# Output is a matrix with the middle points of the used shells and the density of each one
# Usage> profile.emp(halo, clust, spheres)
#
profile.emp <- function(halo, clust, spheres) {
	
	# Calculate relative distances of each particle with respect to the halo center
	dist <- sqrt((clust$data$x-halo[1])^2 + (clust$data$y-halo[2])^2 + (clust$data$z-halo[3])^2)
	dist <- sort(dist)
	xlim <- max(spheres)
	dist <- dist[dist < xlim]

	# Density of particles
	shells <- hist(dist,breaks=c(0,spheres),plot=FALSE)
	window <- clust$domain
	vol.shells <- volume.shells(halo, shells$breaks, window, M=1e6)[-1]
	den.shells <- shells$counts/vol.shells
	profile <- cbind(shells$mids, den.shells)
	profile <- profile[!is.na(profile[,2]),]
	profile <- profile[!is.infinite(profile[,2]),]

	return(profile)
}

#
# volume.shells (Calculates the volume of the concentric shells where the dataset empyric profile is calcualted)
# 
# Inputs are the parameters of the chosen halo (halo), the radii of each shell
# (spheres), the window data set (window) and the number of points used for the Monte Carlo integreation (M).
# Output is the volume of each shell where the profile is calculated.
# Usage> volume.shell(halo, spheres, window, M=1e6)
#
volume.shells <- function(halo, spheres, window, M = 1e6) {
	
  num.spheres <- length(spheres)
  balls <- c()
  vol.shells <- c()
  i <- 1
  shell.in <- inside(halo, spheres[i], window) 
  while(shell.in == TRUE)
  {
    balls[i] <- (4/3)*pi*spheres[i]^3
    if(i==1) {vol.shells[i] <- balls[i]}
    else {
      vol.shells[i] <- balls[i] - balls[i-1]
    }
    i <- i + 1
    shell.in <- inside(halo, spheres[i], window)
  }
  
	Xpop <- runif(M,window$xrange[1], window$xrange[2])
	Ypop <- runif(M,window$yrange[1], window$yrange[2])
	Zpop <- runif(M,window$zrange[1], window$zrange[2])
	Vol <- volume(window)
	rho <- Vol/M
	dist <- (Xpop-halo[1])^2 + (Ypop-halo[2])^2 + (Zpop-halo[3])^2
	
	spheres2 <- spheres^2
	num.spheres <- length(spheres)
	for(j in i:num.spheres){balls[j] <- length(dist[dist < spheres2[j]])*rho}
	for(j in i:num.spheres){vol.shells[j] <- (balls[j]-balls[j-1])};
	return(vol.shells)
}

#
# mass.ein (Analytical calculation of the mass of an Einasto halo at a radius r)
# 
# Inputs are the parameters of the halo (halo) and radius of the sphere
# centred at the halo where we want to integreate the mass.
# Output is the mass of a halo contained in a sphere of radius r
# Usage> mass.ein(halo=c(0,0,0,1,3), r=2)
#
mass.ein <- function(halo, r) {
  re <- halo[4]
  n <- halo[5]
  dn <- uniroot(dn.fit, c(0, 200), n)$root
  h <- re/(dn^n)
  s <- (dn^n)*r/re  
  G <- gamma(3*n)
  M <- 4*pi*(h^3)*n*G*exp(dn)
  Mr <- M*(1-(gammainc(3*n,s^(1/n))/G))
  Mr
}

#
# profile.ind (Calculates the model profile of a halo as centred in the same halo)
# 
# Inputs are the parameters of the chosen halo (halo), the point process (clust),
# the shells radii (shells) and the model normalizing constant (norm)
# Output is the value of the profile at each shell centered in the halo center
# Usage> profile.ind(halo=c(0,0,0,1,3), clust=clust, shells=shells, norm=norm)
#
profile.ind <- function(halo, clust, shells, norm)
{
  N <- length(clust$data$x)
  window <- clust$domain
  num.spheres <- length(shells$breaks)
  # Volume in the shells
  vol.shells <- volume.shells(halo, shells$breaks, window)[-1]
  # Mass in the shells
  mass.spheres <- mass.ein(halo,shells$breaks)
  mass.shells <- mass.spheres[-1] - mass.spheres[-num.spheres]

  norm.mass <- (N*(10^halo[6])*mass.shells/norm)/vol.shells
  profile <- cbind(shells$mids,norm.mass)
  profile <- profile[!is.na(profile[,2]),]
  profile <- profile[!is.infinite(profile[,2]),]
  
  return(profile)
}

#
# profile.mm (Calculates the model profile of the entire Mixture Model as centred in a halo)
# 
# Inputs are the component where the calculation is centred (comp), the model parameters (param),
# the model functions (param2), the data set (clust), the shells radii (shells), the model normalizing
# constant (norm) and the number of points used for the Monte Carlo integration (M).
# Output is the value of the total profile at each shell centered in the halo center.
# Usage> profile.mm(comp=4, param, param2, clust, shells, norm, M=1e4)
#
profile.mm <- function(comp, param, param2, clust, shells, norm, M = 1e4)
{
  qwe <- c()
  k <- param2[[1]]-1
  v <- param2[[comp*2+1]]
  p1 <- (v+1)*(comp-1)
  p2 <- (v+1)*(comp-1)+v
  halo <- param[p1:p2]
  if(comp==1) {halo <- c(halo, 0)}

  N <- length(clust$data$x)
  window <- clust$domain
  num.spheres <- length(shells$mids)

  profile <- matrix(0, num.spheres, 2)
  profile[,1] <- shells$mids

  for(i in 1:num.spheres)
  {
    rmax <- shells$breaks[i+1]
    rmin <- shells$breaks[i]
    theta = runif(M,0,1)*2.0*pi;
    phi = acos(2.0*runif(M,0,1) - 1.0);
    r = runif(M,rmin^3,rmax^3)^(1/3);
    x = r * sin(phi) * cos(theta) + halo[1];
    y = r * sin(phi) * sin(theta) + halo[2];
    z = r * cos(phi) + halo[3];
    pop <- inside.shell(x,y,z,window); 
    x <- pop[[1]][,1]; y <- pop[[1]][,2]; z <- pop[[1]][,3]; Ms <- pop[[2]];
    for(j in 1:k)
    {
        p1 <- (v+1)*(j-1)
        p2 <- (v+1)*(j-1)+v
        param.pass <- param[p1:p2]
        if(j==1) {param.pass <- c(param.pass, 0)}
        qwe[i] <- sum(einasto.model(x,y,z,param.pass))
        profile[i,2] <- profile[i,2] + sum((10^param.pass[6])*einasto.model(x,y,z,param.pass))/Ms
    }
    # Background
    param.pass <- param[k*(v+1)]
    profile[i,2] <- profile[i,2] + (10^param.pass)
  }
  profile[,2] <- N*profile[,2]/norm
  profile <- profile[!is.na(profile[,2]),]
  profile <- profile[!is.infinite(profile[,2]),]
  
  return(profile)
}

#
# inside.shell (Counts the number of points inside the window)
# 
# Inputs are the coordinates of points (x, y, z) inside the data set window (window).
# Output is a list with the points that lie inside the window and its number.
# Usage> inside.shell(x,y,z,window)
#
inside.shell <- function(x,y,z,window)
{
  pop <- cbind(x,y,z)
  minx <- window$xrange[1]; maxx <- window$xrange[2];
  miny <- window$yrange[1]; maxy <- window$yrange[2];
  minz <- window$zrange[1]; maxz <- window$zrange[2];
  pop <- subset(pop, x>minx & x<maxx & y>miny & y<maxy & z>minz & z<maxz)
  Ms <- length(pop[,1])
  return(list(pop,Ms))
}

#
# inside (Determines if a shell is entirely inside the window)
# 
# Inputs are the halo parameters (halo), the radius of the shell (radius) and the window (window)
# Output is a boolean: TRUE if the shell is entirely inside the window and FALSE
# otherwise
# Usage> inside(halo, radius, window)
#
inside <- function(halo, radius, window) {
  param_new <- halo[1:3]+radius
  plus <- outside.win(param_new,window)
  param_new <- halo[1:3]-radius
  minus <- outside.win(param_new,window)
  return(!plus & !minus)
}

#
# plot.profile (Plots the profiles centred in a halo)
# 
# Inputs are the model parameters (param), the model functions (param2), the data set (clust),
# the component in which the profiles will be centred (comp), the limits where the profile is 
# calculated (xlim), the number of bins (nbin), the title of the produced plot (main) and flag 
# stating if a legend is to be included (flag.legend)
# Output is a matrix containing the produced profiles and a plot.
# Usage> plot.profile(param, param2, clust, 5, xlim=c(0.1, 40), nbin=60, flag.legend=FALSE)
#
plot.profile <- function(param, param2, clust, comp, xlim=c(0.1, 40), nbin=60, main=NULL, flag.legend=TRUE) {
  
  # Extract halo component parameters
  k <- param2[[1]]-1
  v <- param2[[comp*2+1]]
  p1 <- (v+1)*(comp-1)
  p2 <- (v+1)*(comp-1)+v
  halo <- param[p1:p2]
  if(comp==1) {halo <- c(halo, 0)}
  
  # Shell geometry
  spheres <- exp(seq(log(xlim[1]), log(xlim[2]), length.out = nbin))
  shells <- hist(0,breaks=c(0,spheres),plot=FALSE)
  
  # Empyric profile
  emp.pro <- profile.emp(halo, clust, spheres)

  # Single component
  norm <- int.model(param,param2,quad,mixture.model)
  comp.pro <- profile.ind(halo, clust, shells, norm)
  
  # Whole MM
  mm.pro <- profile.mm(comp, param, param2, clust, shells, norm, M = 1e4)

  # Plot
  main <- paste('Component', toString(comp))
  miny <- min(mm.pro[,2]); maxy <- max(c(emp.pro[,2], comp.pro[,2], mm.pro[,2]))
  png("Profile_xy_3D.png", width=880, height=880)
  par(mar=c(6,7,3,0))
  plot(emp.pro[,1], emp.pro[,2], cex=3, pch=16, log="xy", ylim=c(miny, maxy), 
       type="b",xlab="", ylab="", main=main, cex.axis=3, cex.main=3)
  title(xlab=expression(paste("h"^-1, " Mpc", sep="")), line=4, cex.lab=3)
  title(ylab="Particle number density", line=4, cex.lab=3)
  lines(comp.pro[,1], comp.pro[,2], lwd=4, col=2)
  lines(mm.pro[,1], mm.pro[,2], lwd=4, col=3)
  if(flag.legend==TRUE) {
  legend("topright", legend=c("Observed profile", "One component profile", 
                              "Mixture model profile"),
         col=c("black", "red", "green"), lty=1:2, lwd=2, cex=3)}
  dev.off()
  profiles <- as.data.frame(cbind(emp.pro, comp.pro[,2], mm.pro[,2])) 
  colnames(profiles) <- c("r", "Empyrical", "Individual", "Mixture")
  return(profiles)
}
  
#
# profile.sphere (Calculates the integral of the Einasto profile in a shell)
#
# Inputs are the distances where the profile is to be calculated (x), the Einasto radius (rs)
# and the Sérsic index (n).
# Output is the integral of the one dimensional Einasto profile along a shell 
# Usage> profile.sphere(1,1,3)
#
profile.sphere <- function(x,rs,n)
{
  #Perfil einasto
  ein <- einasto.1d(x,rs,n)
  #Product
  ret <- ein*4*pi*x*x
  ret
}

#
# einasto.1d (Calculates the one dimensional Einasto profile)
# 
# Inputs are the distance where the profile is to be calculated (x), the Einasto radius (rs)
# and the Sérsic index (n).
# Output is the value of the Einasto profile
# Usage> einasto.1d(1,1,3)
#
einasto.1d <- function(x,rs,n)
{
  dn <- uniroot(dn.fit,c(0,200),n)$root
  rre <- exp(-dn*((x/rs)^(1/n)-1))
  rre
}
  
#
# around (Checks if a given coordinate is a local maximum)
#
# Inputs are a density estimation as produced by function kde3d from package misc3d, and the 
# location of a coordinate in the file (i,j,k).
# Output is a logical indicating if this point is a local maximum (TRUE) or not (FALSE).
# This function called by function centers().
# Usage> around(dens,10,10,5)
#
around <- function(dens,i,j,k) {
  max <- TRUE
  for(a in -1:1) {
    for(b in -1:1) {
      for(c in -1:1) {
        if(dens[[4]][i,j,k] < dens[[4]][i+a,j+b,k+c]) {max <- FALSE}
      }
    }
  }
  max
}

#
# membership (Calculates the membership probability of each particle)
# 
# Inputs are the parameters (param), the model functions (param2), the mixture.model (model),
# the data set (clust) and a logical variable stating if the plot is to be stored as a figure (print) 
# and its dimensions (w, h).
# Output is a dataframe with as many columns as components + 1 and as many rows as particles.
# Each column contains the probability that particle i belongs to components j. The last column
# is the component assigned to particle i according to a multinomial distribution.
# Usage> membership(param, param2, mixture.model, clust, FALSE, 0, TRUE, 880, 880)
#
membership <- function(param, param2, model, clust, threshold=0, ci=FALSE, sigma=0, num_sam=100, print=TRUE, w=880, h=880) {
  n <- length(clust$data$x)
  k <- param2[[1]]
  arg = list(x1=clust$data$x,y1=clust$data$y,z1=clust$data$z)
  output=matrix(NA,n,param2[[1]])
  prob <- mixture.model(param=param,param2=param2,output,action.member,arg)
  prob <- t(apply(prob, 1, function(x)(x/sum(x))))
  names <- c("p1")
  for(i in 2:k) {names <- c(names, paste("p",i, sep=""))}
  class <- c()
  for(i in 1:n) {
    class[i] <- which.max(rmultinom(1,1,prob[i,]))  		
    if ((threshold > 0) & (all(prob[i,] < threshold))) {
        class[i] <- k
    } else if(max(prob[i,]) == prob[i,k]) {
        class[i] <- k
    }
  }
  ret <- cbind(prob,class)
  colnames(ret) <- c(names,"class")
  ret <- as.data.frame(ret)
  
  # Error bars
  sigma_sam <- c()
  if (ci) {
  	l <- length(param)
  	member_sam <- matrix(NA,n,num_sam)
  	for(j in 1:num_sam) {
     		p <- param + rnorm(l,0,sigma)
        prob <- mixture.model(param=p,param2=param2,output,action.member,arg)
        prob <- t(apply(prob, 1, function(x)(x/sum(x))))
        for(i in 1:n) {
  	        if (all(prob[i,] < threshold)) {
  		        member_sam[i,j] <- k
  	        }
  	        else {
	            member_sam[i,j] <- which.max(rmultinom(1,1,prob[i,]))  		
  	        }
        }
  	}
  	number_sam <- matrix(NA,k,num_sam)
  	for(i in 1:num_sam) {
  	    number_sam[,i] <- as.matrix(table(member_sam[,i]))
  	}
  	for(j in 1:k) {
  	    sigma_sam[j] <- round(sd(number_sam[j,]))
  	}
  }
  
  # Make plot
  if(print==TRUE) {
    name <- "classification.png"
    png(name, width = w, height = h) 
    par(mar=c(0,0,0,0))
  }
  clust2d <- ppp(x=clust$data$x, y=clust$data$y, 
                 window=owin(xrange=clust$domain$xrange, yrange=clust$domain$yrange),
                 marks=ret$class) 
  colmap <- colourmap(rainbow((5+k)), inputs=1:(5+k))
  sy <- symbolmap(pch=21, bg=colmap, inputs=1:(5+k), cex=2, nsymbols=(5+k))
  plot(clust2d, symap=sy, legend=FALSE, main="", leg.args=list(cex=2, cexl=2))
  axis(1)
  axis(2)
  if(print==TRUE){dev.off()}
  
  return(list(ret,sigma_sam))
}  

