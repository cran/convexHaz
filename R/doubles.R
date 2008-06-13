#################################################
# hpar.F, suppF: put support and weights together

hpar.F <- function(alpha, nu, mu) {
	return(list(alpha=alpha, nu=nu, mu=mu))
	}
 
supp.F <- function(constant, tau, eta) {
	return(list(constant = constant, tau = tau, eta = eta))
	}

# basis functions and their integrals

dp <- function(x, d, p = 1) (x > d) * (x - d) ^ p

g0  <- 1
g1  <- function(x, t)         ifelse(x > t, x - t, 0)  # NB. x is the tau
g2  <- function(x, t)         ifelse(t > x, t - x, 0)  # NB. x is the eta
G0  <- function(t)            t
G1  <- function(x, t)         x ^ 2 / 2 - dp(x, t, 2) / 2 
GG1 <- function(x, t)         max(x, t) * min(x, t) ^ 2 / 2 - min(x, t) ^ 3 / 6
G2  <- function(x, t)         dp(t, x, 2) / 2
GG2 <- function(x, t)         dp(t, x, 3) / 6



# hazard function
h <- function(tt, supp, hpar) {
      	k <- length(supp$tau)
	m <- length(supp$eta)
	T1<- if(supp$constant>0){supp$constant*hpar$alpha} else 0
	T2<- if(k>0) sum(hpar$nu*sapply(supp$tau, g1, t=tt)) else 0
	T3<- if(m>0) sum(hpar$mu*sapply(supp$eta, g2, t=tt)) else 0
      	return(T1+T2+T3)  
	} 


# cumulative hazard function
H <- function(tt, supp, hpar) {
	k <- length(supp$tau)
	m <- length(supp$eta)
	T1<- if(supp$constant>0){supp$constant*hpar$alpha*tt} else 0
	T2<- if(k>0) sum(hpar$nu*sapply(supp$tau, G1, t=tt)) else 0
	T3<- if(m>0) sum(hpar$mu*sapply(supp$eta, G2, t=tt)) else 0
      	return(T1+T2+T3)  
      }  


# integral of cumulative hazard function
HH <- function(tt, supp, hpar) {
      	k <- length(supp$tau)
	m <- length(supp$eta)
	T1<- if(supp$constant>0){supp$constant*hpar$alpha*0.5*(tt)^2} else 0
	T2<- if(k>0) sum(hpar$nu*sapply(supp$tau, GG1, t=tt)) else 0
	T3<- if(m>0) sum(hpar$mu*sapply(supp$eta, GG2, t=tt)) else 0
	return(T1+T2+T3)
	}  



# empirical distribution function
FFn <- function(t, xx){
	xx <- sort(xx)
 	return(sum(t >= xx)/length(xx))
	}

# empirical cumulative hazard function
HHn <- function(t, xx){
	xx <- sort(xx) 
	return(sum ((t >= xx)/(length(xx):1)))
	}

# integral of empirical cumulative hazard function
YYn <- function(t, xx){
	xx <- sort(xx)
	return(sum((t >= xx)*(t-xx)/(length(xx):1)))
	}
