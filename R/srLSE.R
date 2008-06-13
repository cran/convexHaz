# support reduction algorithm for LSE

srLSE <- function(x, a, TT, M = 100, GRIDLESS=0, tol=1e-8, max.loop=250, print=0, solve.tol=.Machine$double.eps) {

	# check that tolerance isn't too small 
	if (tol < 1e-13){
		tol <- 1e-13
		if (print==1) cat("Tolerance too small.  Set to 1e-13.", "\n")
		}

	xx <- sort(x)        
	
	a=min(a,TT)
	GRID <- seq(0, TT, length = M + 1)
	GRID1 <- GRID[GRID <= a]
	GRID2 <- GRID[GRID >= a]

	PROBLEM <- FALSE
    
	supp <- supp.F(0, numeric(), numeric())       		# initialize support
	hpar <- hpar.F(numeric(), numeric(), numeric())		# initialize weights
          
	for (i in 1:max.loop) {
 
		if (PROBLEM == TRUE)  break

		# set the initial estimate 
		if (i == 1) {
			supp$constant <- 1
			hpar$alpha <- HHn(TT, xx)/TT  # this is the same as min.phi.F, but don't have to worry about warnings
                 	}
 
		if(print == 1){
			cat("iteration #", i, "\n")
			cat(" support: tau = ", supp$tau, "c = ", supp$constant, "eta = ", supp$eta, "\n")
			cat(" weights: nu = ", hpar$nu, "alpha = ", hpar$alpha, "mu = ", hpar$mu, "\n")
			}

		# calculate directional derivatives for each point in GRID
		DirDer <- dir.der.phi.F(supp, hpar, GRID1, GRID2, TT, xx)


		if (DirDer$min > - tol) {
			if(print==1) cat("min of derivatives:", DirDer$min ,"\n")
			break		# NPLSE found, exit loop
			}

		if(print == 1){
			cat("min of derivatives:", DirDer$min, "with case:", DirDer$case, "and new point:", if(DirDer$case==1)GRID1[DirDer$index] else if(DirDer$case==2)GRID2[DirDer$index] else 0, "\n")
			}  
 
 
		# check that proposed support point isn't already in support (else abort) 
		if (DirDer$case == 1 && GRID1[DirDer$index] %in% supp$tau  ||  DirDer$case == 2 &&
                       GRID2[DirDer$index] %in% supp$eta ||
                       DirDer$case == 0 && supp$constant > 0) {
                        PROBLEM <- TRUE
                        if(print==1) cat("new point is already in support\n")
                        next
                 		}

 
		# add new support point which corresponds the minDD
		if(DirDer$case == 1) {
			k <- length(supp$tau)
			supp$tau[k + 1] <- GRID1[DirDer$index]
			}
		if(DirDer$case == 2) {
			m <- length(supp$eta)
			supp$eta[m + 1] <- GRID2[DirDer$index]
			}
		if(DirDer$case == 0) {
			supp$constant <- 1
			}

		# if GRIDLESS=1, augment GRID (NB. second case must be first)
		if ((DirDer$case==2)&(GRIDLESS==1)){
			loc		<-	DirDer$index+length(GRID1)- length(intersect(GRID1, GRID2))
			newpoints 	<- 	0.5*c(GRID[loc]+GRID[min(loc+1, length(GRID))], GRID[loc]+GRID[max(1,loc-1)])
			GRID		<- 	sort(union(GRID, newpoints))
			GRID1 	<- 	GRID[GRID <= a]
			GRID2 	<- 	GRID[GRID >= a]
			}
		if ((DirDer$case==1)&(GRIDLESS==1)){
			loc		<-	DirDer$index
			newpoints 	<- 	0.5*c(GRID[loc]+GRID[min(loc+1, length(GRID))], GRID[loc]+GRID[max(1,loc-1)])
			GRID		<- 	sort(union(GRID, newpoints))
			GRID1 	<- 	GRID[GRID <= a]
			GRID2 	<- 	GRID[GRID >= a]
			}

		# calculate new weights corresponding to the new support
		wnew <- min.phi.F(supp, xx, TT, solve.tol)

		prop.PROBLEM <- FALSE			
		if(DirDer$case == 1) {
			k	<-	length(supp$tau)
			if (wnew$warn==FALSE)	prop.PROBLEM <- (wnew$newpar$nu[k]<0)
                 	}
            if(DirDer$case == 2) {
			m	<-	length(supp$eta)
			if (wnew$warn==FALSE)	prop.PROBLEM <- (wnew$newpar$mu[m]<0)
				}
            if(DirDer$case == 0) {
                  if (wnew$warn==FALSE)	prop.PROBLEM <- (wnew$newpar$alpha<0)
                 	}


		if ((wnew$warn==TRUE)||(prop.PROBLEM==TRUE)) {

			if (print==1) cat("PROBLEM with linear solve occurred.  Auto fix in progress.", "\n")
			if(DirDer$case == 1) {
                  	k 	<- 	length(supp$tau)-1
                        k0 	<- 	order(abs(supp$tau[1:k]-supp$tau[k+1]))[1]
				supp$tau	<-	supp$tau[-k0]
				hpar$nu  	<-	hpar$nu[-k0]
                 		}
                 	if(DirDer$case == 2) {
                  	m 	<- 	length(supp$eta)-1
                        m0 	<- 	order(abs(supp$eta[1:m]-supp$eta[m+1]))[1]
				supp$eta 	<-	supp$eta[-m0]
				hpar$mu 	<-	hpar$mu[-m0]
                 		}
                 	if(DirDer$case == 0) {
                        PROBLEM <- TRUE
				if (print==1) cat("PROBLEM with linear solve occurred.  Abort.", "\n")
				break
                 		}

			wnew <- min.phi.F(supp, xx, TT, solve.tol)

			if(wnew$warn==TRUE){
				PROBLEM <- TRUE
				if (print==1) cat("PROBLEM with linear solve occurred.  Abort.", "\n")
				break
				}

			}

		wnew <- wnew$newpar 

		# check if all the weights are positive, otherwise, do support reduction step
		SR <- FALSE
		if (min( wnew$nu, wnew$alpha, wnew$mu) < 0) SR <- TRUE

		if (SR == FALSE) {    	# all weights positive, go to the next iterate
			hpar <- wnew 	# change currect estimate		
			next              # continue to the next loop
			}
			
		# "else" do support reduction step
		# first, modify comparison weights
		if(DirDer$case == 0) hpar$alpha<-0       
		if(DirDer$case == 1) hpar$nu<-c(hpar$nu,0)    
		if(DirDer$case == 2) hpar$mu<-c(hpar$mu,0)    
                
		newsupport <- supp.red.F(supp, hpar, wnew, xx, TT, solve.tol) 
		if(newsupport$prob==TRUE){
			PROBLEM <- TRUE
			if (print==1) cat("PROBLEM with linear solve occurred.  Abort.", "\n")
			break
			}
		hpar <- newsupport$weight
		supp <- newsupport$support	  
	
		} # end of main loop in full version

	if (i == max.loop) PROBLEM <- TRUE  # algorithm didn't converge
	if ((PROBLEM == TRUE)&(print == 1))  cat("Problem occurred.")
	if ((PROBLEM == FALSE)&(print == 1)) {  # display support points with corresponding weights
		cat("\n")
		cat("FINAL RESULT: ",i, " iterations", "\n")
		cat("support: tau= ", supp$tau, " c= ", supp$constant, " eta= ", supp$eta, "\n")
		cat("weights: nu= ", hpar$nu, " alpha= ", hpar$alpha, " mu= ", hpar$mu, "\n")
		}

	phi <- phi.F(supp, hpar, xx, TT) 

	return(list(lse=list(supp=supp, mix=hpar), ls=phi, iter=i, conv=(PROBLEM==0)))

	} # end of srLSE()

####################################################################
# the support reduction step as a separate function

supp.red.F <- function(supp, hpar, hparnew, xx, TT, solve.tol) {

	sr.PROBLEM <- FALSE
	m <- length(supp$eta)
	k <-  length(supp$tau)

	test0 <- if(supp$constant != 0) hparnew$alpha / hpar$alpha  else 1  
	test1 <- if(k>0) hparnew$nu / hpar$nu else 1  
	test2 <- if(m>0) hparnew$mu / hpar$mu else 1  

	# vector of minimum [new weight] / [old weight] of three types of support points
	minTest <- c(test0, min(test1), min(test2))

	# find which of eta, constant, tau needs to be removed
	case <- which(minTest == min(minTest, na.rm = T) )[1] - 1
 
	# figure out lambdastar        
	if(case == 0)	lambdastar <- 1 / (1 - test0)               
	if(case == 1) {
		j0 <- which(test1 == min(test1))[1]
		lambdastar <- 1 / (1 - test1[j0])
		}
	if(case == 2) {
		j0 <- which(test2 == min(test2))[1]  
		lambdastar <- 1 / (1 - test2[j0])
		}

	# calculate weighted average new weights alphatilde, nutilde, mutilde
	alphatilde <- (1 - lambdastar) * hpar$alpha + lambdastar * hparnew$alpha
 	nutilde <- rep(0, k)
	nutilde <- (1 - lambdastar) * hpar$nu + lambdastar * hparnew$nu
	mutilde <- rep(0, m)
	mutilde <- (1 - lambdastar) * hpar$mu + lambdastar * hparnew$mu
            
        
	# remove the boundary point where mutilde equals zero
	if(case == 0) {
		supp$constant <- 0  
		alphatilde <- 0                
		}        
	if(case == 1) {
		supp$tau <- supp$tau[-j0]
		nutilde <- nutilde[-j0]        
		}
	if(case == 2) {
		supp$eta <- supp$eta[-j0]
		mutilde <- mutilde[-j0]        
		}
  
	hpartilde<-hpar.F(alphatilde, nutilde, mutilde)	# set up "old" weights

	# find new weights                
	hpar0 <- min.phi.F(supp, xx, TT, solve.tol)
	if( hpar0$warn==TRUE ) { 
		sr.PROBLEM <- TRUE
		hpar0 <- hpar0$newpar
		} else {  hpar0<-hpar0$newpar  }

          
	# call recursively  
	if ( (min(hpar0$alpha, hpar0$nu, hpar0$mu) < 0)&&(sr.PROBLEM==FALSE) ) supp.red.F(supp, hpartilde, hpar0, xx, TT, solve.tol)  else (list(weight = hpar0, support = supp, prob=sr.PROBLEM))

}# end of supportReductionF()



######################################################
# Return the minimum directional derivative, case (0, 1 or 2) and index.

dir.der.phi.F <- function(supp, hpar, GRID1, GRID2, TT, xx) {

	DirDerFn1 <- function(t) 	HH(t, supp, hpar) - YYn(t, xx) 
	DirDerFn2 <- function(t)	(TT - t)*(H(TT, supp, hpar) - HHn(TT, xx))-(HH(TT, supp, hpar) - YYn(TT, xx)) +(HH(t, supp, hpar) - YYn(t, xx))

	DirDer1 <- sapply(GRID1, DirDerFn1)
	minDirDer1 <- min(DirDer1)		
	
	DirDer0 <- H(TT, supp, hpar) - HHn(TT, xx)
	
	DirDer2 <- sapply(GRID2, DirDerFn2)
	minDirDer2 <- min(DirDer2)

	temp <- c(minDirDer1, minDirDer2, DirDer0)
	minDirDer <- min(temp)


	case <- (which( temp == minDirDer ))[1]

	if ( case == 1 ) index <- (which( DirDer1 == minDirDer ))[1] else if
	   ( case == 2 ) index <- (which( DirDer2 == minDirDer ))[1] else
	   { case <- 0 ; index <- 1 }

      return(list(min=minDirDer, case=case, index=index))
 
} # end of dir.der.phi.F

##################################################
# Return the A matrix and b function that goes into min.phi.F

ab.phi.F	<- function(supp, xx, TT){

	hint02 <- function(x)      	(TT-x)^2 /2
	hint01 <- function(x)         x^2/ 2
	hint11 <- function(x, y)     	min(x,y)^2 /6 * (3*max(x,y)-min(x,y))
	hint22 <- function(x, y)    	TT^3/3 - (x+y)*TT^2 /2 + x*y*TT - (max(x,y)^2 /6) * (3*min(x,y)-max(x,y))
	hint12 <- function(x, y)      ifelse(y < x, x^3 /6 - x^2 * y/2 + x*y^2 /2 - y^3 / 6, 0)

	k <- length(supp$tau)
	m <- length(supp$eta)
	c <- supp$constant

	A11 <- if(k>0) 	matrix(sapply(1:(k*k),function(x,i)hint11(x[i,1],x[i,2]),x=matrix(c(rep(supp$tau,each=k), rep(supp$tau,k)),k*k,2)),k,k) else numeric()
	A10 <- if(k*c>0) 	matrix(hint01(supp$tau[1:k]),k,c)  else numeric()
	A12 <- if(k*m>0) 	matrix(0,k,m) else numeric()
	A01 <- if(k*c>0) 	t(matrix(hint01(supp$tau[1:k]),k,c))  else numeric()
	A00 <- if(c>0) 	matrix(TT,c,c) else numeric()
	A02 <- if(c*m>0) 	matrix(hint02(supp$eta[1:m]),c,m) else numeric()
	A21 <- if(k*m>0) 	matrix(0,m,k) else numeric()
	A20 <- if(c*m>0) 	t(matrix(hint02(supp$eta[1:m]),c,m)) else numeric()
	A22 <- if(m>0)  	matrix(sapply(1:(m*m),function(x,i)hint22(x[i,1],x[i,2]),x=matrix(c(rep(supp$eta,each=m), rep(supp$eta,m)),m*m,2)),m,m) else numeric()

	#NB. although A21=t(A12), the two must be programmed separately, so that the empty matrix is read correctly later on

	A1 <- if(k>0) cbind(A11, A10, A12) else numeric()
	A0 <- if(c>0) cbind(A01, A00, A02) else numeric()
	A2 <- if(m>0) cbind(A21, A20, A22) else numeric()
                
	A<-rbind(A1, A0, A2)

	b1<-if(k>0) sapply(supp$tau, YYn, xx=xx) else numeric()
	b0<-if(c>0) HHn(TT, xx) else numeric()
	b2<-if(m>0) sapply(supp$eta, function(t,xx,TT)(TT - t) * HHn(TT, xx) - YYn(TT, xx) + YYn(t, xx),xx=xx, TT=TT) else numeric()

	b<-c(b1,b0,b2)

	return(list(A=A,b=b))

	} # end of ab.phi.F


##################################################
# Return the minimum of criterion for fixed support

min.phi.F <- function(supp, xx, TT, solve.tol){

	k <- length(supp$tau)
	m <- length(supp$eta)
	c <- supp$constant

	Ab	<-	ab.phi.F(supp, xx, TT)

	weights <- 	numeric()
	try(weights<-solve(Ab$A, Ab$b, tol=solve.tol), silent=TRUE) # allow tolerance as input
	if (length(weights)==0){
		# exit with warning
		newpar <- numeric()
		warn	 <- TRUE
		} else { # do as before
			nu <- if(k>0) weights[1:k] else numeric()
			alpha <- if(c>0) weights[k+1] else numeric()
			mu <- if(m>0) weights[(k+c+1):(k+c+m)] else numeric()

			newpar <- hpar.F(alpha, nu, mu)
			warn	 <- FALSE
			}

	return(list(newpar=newpar, warn=warn))
	}# end of min.phi.F
 
####################################################
# Return the criterion for fixed support and weights

phi.F <- function(supp, hpar, xx, TT){

	Ab <- ab.phi.F(supp, xx, TT)
	A <- Ab$A
	b <- Ab$b

	tce <- c(hpar$nu, hpar$alpha, hpar$mu)

	return(Phi <- 0.5 * tce %*% A %*% tce - tce %*% b)
	}
