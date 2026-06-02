#' Bandwidth for the PI estimator in the case of the logarithmic transformation.
#'
#' Bandwidth for the PI estimator for estimating \eqn{g}, the true probability density function (pdf), correponding to the logarithmic transformation of data based on replacing \eqn{g} by the reference \eqn{N(\bar y,s_Y^2)} distribution. See Section 7 of the article by Savchuk (2026).
#'
#' @param x vector of the original X-data.
#' @return The bandwidth for the PI estimator for estimating \eqn{g}, the density after the logarithmic transformation, obtained by replacing \eqn{g} by the reference \eqn{N(\bar y,s_Y^2)} density.
#' @references Savchuk, O. (2026, under review). Density estimation for log-transformed data.
#' @seealso \code{\link{KDE_ln_PI}}.
#' @examples
#' # Example. Finding bandwidth for the PI estimator for the log-transformed sample of size n=100 
#' #originated from a gamma distribution 
#' #with the parameters: shape=2, scale=2.
#' xdat=rgamma(100,2,scale=2)
#' paste("PI bandwidth for the PI estimator=", h_ln_PI(xdat))
#' @export
#' @importFrom stats sd
h_ln_PI=function(x){
  #PI bandwidth for KDE_ln_PI based on the normal reference rule for g (see Section 7 of the article)
  #x is a vector of the original X-data
  n=length(x)
  y=log(x)
  mu_Y=mean(y)
  sigma_Y=sd(y)
  sigma_Y*(4*exp(5*mu_Y-7*sigma_Y^2/2)/(2*sigma_Y^2+3))^(1/5)*n^(-1/5)
}


#' The PI kernel density estimate in the case of the logaritjmic transformation.
#'
#' Computing the PI kernel density estimate for the log-transformed data based on the Gaussian kernel.
#'
#' @param h the bandwidth, scalar,
#' @param x vector of the original X-data,
#' @param u vector of the estimation points.
#' @return Vector of density estimates at the corresponding values of the argument (u). 
#' @references Savchuk, O. (2026, under review). Density estimation for log-transformed data.
#' @seealso \code{\link{h_ln_PI}}.
#' @examples
#' #Example(Gamma(shape=2,scale=2) density). Estimating the pdf of the density corresponding to the
#' #logarithmic transformation of sample of size n=100 originated from the gamma density.
#' xdat=rgamma(100,2,scale=2)  #original data
#' n=length(xdat)  #sample size
#' u=seq(-3,5,len=1000) #estimation points
#' dens_lntransform=exp(2*u-exp(u)/2)/4  #density after log-transformation
#' h_PI=h_ln_PI(xdat)
#' dev.new()
#' dens_PI=KDE_ln_PI(h_PI,xdat,u)
#' plot(u,dens_PI,'l',ylim=c(0,1.05*max(dens_PI)),lwd=2,cex.lab=1.7,cex.axis=1.7,cex.main=1.5,
#' main="Gamma dens. after ln-transform. and its PI estimate",xlab="u",ylab="density")
#' lines(u,dens_lntransform,lwd=2,lty="dashed",col="red",bty="n")
#' legend(-3,0.5,legend=c("true density","PI estimate"),lwd=c(2,2),col=c("red","black"),
#' lty=c("dashed","solid"),cex=1.5,bty="n")
#' text(3.25,0.525,paste("n=",n),cex=2)
#' text(3.85,0.425,bquote(h[PI] == .(sprintf("%.4f", h_PI))),cex=2)
#' @export
#' @importFrom "grDevices" "dev.new"
#' @importFrom "graphics" "plot" "title"
#' @importFrom stats dnorm
KDE_ln_PI=function(h,x,u){    #Plug-in estimator in the case of the logarithmic transformation
  #x is a vertor of the original X-data
  # h is the bandwidth(scalar)
  # u is the vector of values where we compute our density estimate
  n=length(x)
  arg=(matrix(exp(u),length(u),1)%*%matrix(1,1,n)-matrix(1,length(u),1)%*%matrix(x,1,n))/h
  (as.numeric(dnorm(arg)%*%matrix(1,n,1)/(n*h)))*exp(u)
}

#' Adaptive normal PI estimate in the case of the logarithmic transformation.
#'
#' Computing the adaptive normal PI kernel density estimate for the log-transformed data based on the Gaussian kernel.
#'
#' @param x vector of the original X-data,
#' @param u vector of the estimation points.
#' @return Vector of density estimates at the corresponding values of the argument (u). 
#' @references Savchuk, O. (2026, under review). Density estimation for log-transformed data.
#' @seealso \code{\link{KDE_ln_PI}}.
#' @examples
#' #Example(skewed bimodal density). Estimating the density corresponding to the logarithmic 
#' #transformation of a sample of size n=1000 originated from the skewed bimodal density 
#' #(density #8 of Marron and Wand (1992)).
#' #Density after the logarithmi transformation
#' g_skew=function(u)
#'   3/4*dnorm(u)+1/4*dnorm(u,mean=3/2,sd=1/3) 
#' skew_realiz=function(n){
#'   # n is the sample size
#'   mix_comp=runif(n)<0.25
#'   dat=numeric(n)
#'   dat[mix_comp]=rnorm(sum(mix_comp),mean=3/2,sd=1/3)
#'   dat[!mix_comp]=rnorm(sum(!mix_comp))
#'  dat
#' }
#' #Plotting an estimate
#' n=1000
#' u=seq(-4,4,len=1000)
#' y_skew=skew_realiz(n)
#' x_skew=exp(y_skew)
#' h_PI=h_ln_PI(x_skew)
#' KDE_skew_PI=KDE_ln_PI(h_PI,x_skew,u)
#' KDE_skew_adapt_norm=KDE_PI_adapt_norm(x_skew,u)
#' dev.new()
#' plot(u,g_skew(u),'l',lwd=2,ylim=c(0,1.05*max(KDE_skew_PI)),col="red",xlab="y",ylab="density",
#' cex.lab=1.7,cex.axis=1.7,cex.main=1.7,main="Skewed bimodal distribution")
#' lines(u,KDE_skew_PI,lty="longdash",lwd=2)
#' lines(u,KDE_skew_adapt_norm,lty="solid",lwd=2,col="blue")
#' legend(-4.5,0.45,legend=c("PI estimate","Adaptive normal","density g"),lty=c("longdash",
#' "solid","solid"),col=c("black","blue","red"),lwd=c(3,3,3),bty="n",cex=1.7)
#' legend(1.25,0.45,legend=paste("n=",n),cex=2,bty="n")
#' @export
#' @importFrom "grDevices" "dev.new"
#' @importFrom "graphics" "plot" "title"
#' @importFrom stats dnorm quantile sd
KDE_PI_adapt_norm=function(x,u){    #PI estimator in the case of the ln transformation based on the PI reference bandwidth based on the N(mu_Y,sigma_Y^2) distribution used for g
  #x is a vertor of the original data
  # u is the vector of values where we compute our density estimate
  n=length(x)
  l=length(u)
  y=log(x)
  Y999=quantile(y,probs=0.999)[[1]] #99.9-th percentile
  Y001=quantile(y,probs=0.001)[[1]] #0.1-th percentile
  u_new=pmax(pmin(u,Y999),Y001)
  mu_Y=mean(y)
  sigma_Y=sd(y)
  u1=mu_Y-(3*sigma_Y^2+sigma_Y*sqrt(sigma_Y^2+4))/2 #first point of discontinuity for Bfnorm
  u2=mu_Y-(3*sigma_Y^2-sigma_Y*sqrt(sigma_Y^2+4))/2 #second point of discontinuity for Bfnorm
  u11=(u1-0.5*sigma_Y)
  u12=(u1+0.5*sigma_Y)
  u21=(u2-0.5*sigma_Y)
  u22=(u2+0.5*sigma_Y)
  f=function(y)
    exp((y-mu_Y)^2/(10*sigma_Y^2))*exp(y)/((abs((y-mu_Y)^2+3*(sigma_Y^2)*(y-mu_Y)-sigma_Y^2+2*sigma_Y^4))^(2/5))
  Const=(sigma_Y^(9/5))/(2^(1/10))  #constant
  func=ifelse(u_new<u11,f(u_new),
         ifelse(((u_new>=u11)&(u_new<=u12)),(f(u11)+(f(u12)-f(u11))*(u_new-u11)/sigma_Y),
		   ifelse(((u_new>u12)&(u_new<u21)),f(u_new),
		     ifelse(((u_new>=u21)&(u_new<=u22)),(f(u21)+(f(u22)-f(u21))*(u_new-u21)/sigma_Y),f(u_new)))))
  band=Const*func*n^(-1/5)
  C=matrix(band,l,1)%*%matrix(1,1,n)  #matrix of bandwidths
  arg=(matrix(exp(u),l,1)%*%matrix(1,1,n)-matrix(1,l,1)%*%matrix(x,1,n))/C
  (as.numeric(dnorm(arg)%*%matrix(1,n,1)))*exp(u)/(n*band) #array of size l
}

#' The PI kernel density estimate in the case of the power transformation.
#'
#' Computing the PI kernel density estimate for data after a power transformation based on the Gaussian kernel.
#'
#' @param alpha the exponent of the power transformation (alpha>0),
#' @param h the bandwidth, scalar,
#' @param x vector of the original X-data,
#' @param u vector of the estimation points.
#' @return Vector of density estimates at the corresponding values of the argument (u). 
#' @references Savchuk, O., Schick A. (2013). Density estimation for power transformations. \emph{Journal of Nonparametric Statistics}, 25(3), 545-559.
#' @seealso \code{\link{KDE_ln_PI}}.
#' @examples
#' #Example(original data is sandard normal, n=1000 and (1)alpha=0.5 (Case 1), (2) alpha=2 (Case 2))
#' n=1000
#' xdat=rnorm(n)
#' h_P=bw.SJ(xdat) #Sheather-Jones PI for original x-data
#'
#' #Case1: alpha=0.5
#' alpha=0.5
#' u=seq(-2.5,2.5,len=1000)
#' g_0.5=1/alpha/sqrt(2*pi)*abs(u)^(1/alpha-1)*exp(-abs(u)^(2/alpha)/2) #true density
#' PI_power_est_0.5=KDE_power_PI(alpha,h_P,xdat,u)
#' dev.new()
#' plot(u,g_0.5,lwd=2,'l',col="red",ylim=c(0,1.125*max(PI_power_est_0.5)),xlab="y",ylab="density",
#' cex.axis=1.7,cex.main=1.7,cex.lab=1.7,main=expression(paste("Norm density, ",alpha,"=0.5")))
#' lines(u,PI_power_est_0.5,lwd=2,lty="dashed")
#' legend(-2.75,1.2*max(max(PI_power_est_0.5),0.52),legend=c("true density","PI estimate"),
#' lty=c("solid","dashed"),col=c("red","black"),lwd=c(2,2),bty="n",cex=1.5)
#' legend(0.75,0.5,legend=paste("n=",n),cex=2,bty="n")
#'
#' #Case 2: alpha=2
#' alpha=2
#' y=sign(xdat)*abs(xdat)^alpha
#' u=seq(-1.5,1.5,len=1000)
#' g_2=1/alpha/sqrt(2*pi)*abs(u)^(1/alpha-1)*exp(-abs(u)^(2/alpha)/2) #true density
#' PI_power_est_2=KDE_power_PI(alpha,h_P,xdat,u)
#' dev.new()
#' plot(u,g_2,lwd=2,'l',col="red",ylim=c(0,1.1*max(PI_power_est_2)),xlab="y",ylab="density",
#' cex.axis=1.7,cex.main=1.7,cex.lab=1.7,main=expression(paste("Norm density, ",alpha,"=2")))
#' lines(u,PI_power_est_2,lwd=2,lty="dashed")
#' legend(0,2,legend=c("true density","PI estimate"),lty=c("solid","dashed"),col=c("red","black"),
#' lwd=c(2,2),cex=1.5,bty="n")
#' legend(-1.5,5,legend=paste("n=",n),cex=2,bty="n")
#' @export
#' @importFrom "grDevices" "dev.new"
#' @importFrom "graphics" "plot" "title"
#' @importFrom stats dnorm
KDE_power_PI=function(alpha,h,x,u){
    # alpha is the exponent of the power transformation
    # h is the bandwidth(scalar)
    # x is a vertor of the original X-data
    # u is the vector of values where we compute the density estimate
    n=length(x)
    arg=(matrix((sign(u)*abs(u)^(1/alpha)),length(u),1)%*%matrix(1,1,n)-matrix(1,length(u),1)%*%matrix(x,1,n))/h
    vec=dnorm(arg)%*%matrix(1,n,1)
    as.numeric(diag(abs(u)^(1/alpha-1))%*%vec/(n*h*alpha))
}
