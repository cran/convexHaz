# support reduction algorithm for the MLE

srMLE <- function(x, a, M = 100, GRIDLESS=0, ini.LSE=0, tol=1e-8, max.inner.loop=250, max.outer.loop=50, print=0,  solve.tol=.Machine$double.eps) {
 
	if (tol < 1e-13){
		tol = 1e-13
		if (print==1) cat("Tolerance too small.  Set to 1e-13.", "\n")
		}

	xx <- sort(x)	# sort original data in increasing order

	a=min(a, max(xx)) # not really necessary, but...
	GRID <- seq(0, max(x), length = M + 1)
      	GRID1 <- GRID[GRID <= a]
      	GRID2 <- GRID[GRID >= a]

	PROBLEM <- FALSE
    
	supp <- supp.F(0, numeric(), numeric())       		# initialize support
	hpar <- hpar.F(numeric(), numeric(), numeric())		# initialize weights

	for(i in 1:max.outer.loop){		

		if (PROBLEM == TRUE) break

		# set the initial estimate for bathtub hazard function
      		if (i == 1) {
			if (ini.LSE==1){
				if (print==1) cat("find LSE", "\n")
				TT	<- xx[floor(0.95*length(xx))]
				h.LSE <- srLSE(xx, a, TT=TT, M, tol, max.loop=max.inner.loop, print=print, GRIDLESS=0)
				if (h.LSE$conv == 0){
					PROBLEM <- TRUE
					if (print==1) cat("Initial LSE failed to converge.  Abort.", "\n")
					break
					}
				supp  <- h.LSE$lse$supp
				hpar	<- h.LSE$lse$mix
				if (print==1) cat("LSE found. Set as initial estimate", "\n") 
				if (print==1) cat("\n")
				} else {
					supp$constant <- 1
					hpar$alpha    <- 1/mean(xx)
					}  
			}
		
      	if(print == 1){
			cat("OUTER iteration #", i, "\n")
            	cat(" support: tau = ", supp$tau, "c = ", supp$constant, "eta = ", supp$eta, "\n")
            	cat(" weights: nu = ", hpar$nu, "alpha = ", hpar$alpha, "mu = ", hpar$mu, "\n")
			}

      	DirDer <- dir.der.psi.F(supp, hpar, GRID1, GRID2, xx)  # calculate directional derivatives 
		if(print == 1){
			cat("min of derivatives: ", DirDer$min, "\n")
			}  

      
      	if (DirDer$min > -tol) break	# NPMLE found, exit outer loop

      	# INNER LOOP:

		cur.supp <-	supp
		cur.par  <-	hpar
		cur.vec  <- sapply(xx, h, supp=cur.supp, hpar=cur.par)
		
		for (j in 1 : max.inner.loop) {  
 
            	if (PROBLEM == TRUE)  break

                 	# set the initial estimate for bathtub hazard function
                 	if (j == 1) {
				inner.supp 	<-	cur.supp
                        inner.par  	<- 	min.psi.approx.F(inner.supp, cur.vec, xx, solve.tol)
				temp.par	<-	cur.par

				if (inner.par$warn==TRUE) {
					
					if (print==1) cat("PROBLEM with linear solve occurred.  Auto fix in progress.", "\n")
					
					tau.diff		<- 	abs(diff(sort(inner.supp$tau)))
					min.diff.tau	<-	if(length(tau.diff)>0) min(tau.diff) else Inf
					eta.diff		<- 	abs(diff(sort(inner.supp$eta)))
					min.diff.eta	<-	if(length(eta.diff)>0) min(eta.diff) else Inf

					if  (min.diff.tau <= min.diff.eta){ # remove one tau
						remove <- order(tau.diff)[1]
						out	 <- which(inner.supp$tau==sort(inner.supp$tau)[remove])
						inner.supp$tau <- inner.supp$tau[-out]
						temp.par$nu	   <- temp.par$nu[-out]							
						} else { # remove one eta
							remove <- order(eta.diff)[1]
							out	 <- which(inner.supp$eta==sort(inner.supp$eta)[remove])
							inner.supp$eta <- inner.supp$eta[-out]
							temp.par$mu	   <- temp.par$mu[-out]
							}


					inner.par <- min.psi.approx.F(inner.supp, cur.vec, xx, solve.tol)

					if(inner.par$warn==TRUE){
						PROBLEM <- TRUE
						if (print==1) cat("PROBLEM with linear solve occurred.  Abort.", "\n")
						break
						}

					}
			
				inner.par<-inner.par$newpar

				ini.SR <- FALSE
				if (min(inner.par$alpha, inner.par$nu, inner.par$mu) < 0) ini.SR <- TRUE

				if (ini.SR == TRUE){
					inner.new <- supp.red.MLE.F(inner.supp, temp.par, inner.par, cur.vec, xx, solve.tol)
					if(inner.new$prob==TRUE){
						PROBLEM <- TRUE
						if (print==1) cat("PROBLEM with linear solve occurred.  Abort.", "\n")
						break
						}
                 			inner.par <- inner.new$weight
                 			inner.supp <- inner.new$support
					}
				}
 

                 	if(print == 1){
				if(j == 1)  cat("INNER LOOP:" , "\n")
		     		cat("iteration #", j, "\n")
                 		cat(" support: tau = ", inner.supp$tau, "c = ", inner.supp$constant, "eta = ", inner.supp$eta, "\n")
                 		cat(" weights: nu = ", inner.par$nu, "alpha = ", inner.par$alpha, "mu = ", inner.par$mu, "\n")
		     	}

 
                 	# calculate directional derivatives for each point in GRID
                	DirDerInner <- dir.der.psi.approx.F(inner.supp, inner.par, cur.vec, GRID1, GRID2, xx)


                 	if (DirDerInner$min > -tol) {
				if(print==1) cat("min of derivatives:",DirDerInner$min ,"\n")
				break	# exit inner loop
				}

                 	if(print == 1){
		     	cat("min of derivatives:", DirDerInner$min, "with case:", DirDerInner$case, "and new point:", if(DirDerInner$case==1)GRID1[DirDerInner$index] else if(DirDerInner$case==2)GRID2[DirDerInner$index] else 0, "\n")
			}  

 
                 	# if new point is already in support, then we have a problem  
                 	if (DirDerInner$case == 1 && GRID1[DirDerInner$index] %in%
                  	inner.supp$tau  ||  DirDerInner$case == 2 &&
                       	GRID2[DirDerInner$index] %in% inner.supp$eta ||
                       	DirDerInner$case == 0 && inner.supp$constant > 0) {
                        	PROBLEM <- TRUE
                        	cat("new point is already in support\n")
                        	next
                 			}
 

                 	# add new support point which corresponds the minDD
                 	if(DirDerInner$case == 1) {
                  	k <- length(inner.supp$tau)
                        inner.supp$tau[k + 1] <- GRID1[DirDerInner$index]
                 		}
                 	if(DirDerInner$case == 2) {
                        m <- length(inner.supp$eta)
                        inner.supp$eta[m + 1] <- GRID2[DirDerInner$index]
                 		}
                 	if(DirDerInner$case == 0) {
                        inner.supp$constant <- 1
                 		}

			# if GRIDLESS=1, augment GRID (NB. second case must be first)
			if ((DirDerInner$case==2)&(GRIDLESS==1)){
				loc		<-	DirDerInner$index+length(GRID1)- length(intersect(GRID1, GRID2))
				newpoints 	<- 	0.5*c(GRID[loc]+GRID[min(loc+1, length(GRID))], GRID[loc]+GRID[max(1,loc-1)])
				GRID		<- 	sort(union(GRID, newpoints))
				GRID1 	<- 	GRID[GRID <= a]
				GRID2 	<- 	GRID[GRID >= a]
				}
			if ((DirDerInner$case==1)&(GRIDLESS==1)){
				loc		<-	DirDerInner$index
				newpoints 	<- 	0.5*c(GRID[loc]+GRID[min(loc+1, length(GRID))], GRID[loc]+GRID[max(1,loc-1)])
				GRID		<- 	sort(union(GRID, newpoints))
				GRID1 	<- 	GRID[GRID <= a]
				GRID2 	<- 	GRID[GRID >= a]
				}

                 	# calculate new weights corresponding to the new support
                 	inner.par.new <- min.psi.approx.F(inner.supp, cur.vec, xx, solve.tol)

			prop.PROBLEM <- FALSE			
			if(DirDerInner$case == 1) {
				k	<-	length(inner.supp$tau)
				if (inner.par.new$warn==FALSE)	prop.PROBLEM <- (inner.par.new$newpar$nu[k]<0)
                 		}
                 	if(DirDerInner$case == 2) {
				m	<-	length(inner.supp$eta)
				if (inner.par.new$warn==FALSE)	prop.PROBLEM <- (inner.par.new$newpar$mu[m]<0)
					}
                 	if(DirDerInner$case == 0) {
                        if (inner.par.new$warn==FALSE)	prop.PROBLEM <- (inner.par.new$newpar$alpha<0)
                 		}


			if ((inner.par.new$warn==TRUE)||(prop.PROBLEM==TRUE)) {

				#cat("first try lm idea:", "\n")
				#ab	<-	ab.psi.F(inner.supp, cur.vec, xx)
				#cat(lm.fit(ab$A, ab$b)$coefficients, "\n")
				#cat(qr.coef(qr(ab$A), ab$b), "\n")
				#cat(solve(ab$A, ab$b, tol=1e-30), "\n")				
					
				if (print==1) cat("PROBLEM with linear solve occurred.  Auto fix in progress.", "\n")
				if(DirDerInner$case == 1) {
                  		k 	<- 	length(inner.supp$tau)-1
                        	k0 	<- 	order(abs(inner.supp$tau[1:k]-inner.supp$tau[k+1]))[1]
					inner.supp$tau	<-	inner.supp$tau[-k0]
					inner.par$nu  	<-	inner.par$nu[-k0]
                 			}
                 		if(DirDerInner$case == 2) {
                  		m 	<- 	length(inner.supp$eta)-1
                        	m0 	<- 	order(abs(inner.supp$eta[1:m]-inner.supp$eta[m+1]))[1]
					inner.supp$eta 	<-	inner.supp$eta[-m0]
					inner.par$mu 	<-	inner.par$mu[-m0]
                 			}
                 		if(DirDerInner$case == 0) {
                        	PROBLEM <- TRUE
					if (print==1) cat("PROBLEM with linear solve occurred.  Abort.", "\n")
					break
                 			}

				inner.par.new <- min.psi.approx.F(inner.supp, cur.vec, xx, solve.tol)

				if(inner.par.new$warn==TRUE){
					PROBLEM <- TRUE
					if (print==1) cat("PROBLEM with linear solve occurred.  Abort.", "\n")
					break
					}

				}

			inner.par.new <- inner.par.new$newpar
 
                 	# check if all the weights are positive, otherwise, do Support Reduction Step
                 	inner.SR <- FALSE
                 	if (min( inner.par.new$nu, inner.par.new$alpha, inner.par.new$mu) < 0) inner.SR <- TRUE  
                 	# cat("support reduction", inner.SR, "\n")  

                 	if (inner.SR == FALSE) {
				inner.par <- inner.par.new
				next
				}                 	
 
                	if(DirDerInner$case == 0) inner.par$alpha<-0       
		     	if(DirDerInner$case == 1) inner.par$nu<-c(inner.par$nu,0)    
                 	if(DirDerInner$case == 2) inner.par$mu<-c(inner.par$mu,0)    
                
                 	newh <- supp.red.MLE.F(inner.supp, inner.par, inner.par.new, cur.vec, xx, solve.tol)
			if(newh$prob==TRUE){
				PROBLEM <- TRUE
				if (print==1) cat("PROBLEM with linear solve occurred.  Abort.", "\n")
				break
				}
                 	inner.par <- newh$weight
                 	inner.supp <- newh$support	  
 
		} # end of inner loop

		# cat("inner iteration", j, "\n")
		if (j==max.inner.loop) PROBLEM <- TRUE

       	if ((PROBLEM == FALSE)&(print == 1)) {	
                	cat("INNER LOOP RESULT: ",j, " iterations", "\n")
	 	    	cat("support: tau= ", inner.supp$tau, " c= ", inner.supp$constant, " eta= ", inner.supp$eta, "\n")
	 	    	cat("weights: nu= ", inner.par$nu, " alpha= ", inner.par$alpha, " mu= ", inner.par$mu, "\n")
		    	}
	
		
		if (PROBLEM==TRUE) break
		# ELSE do ARMIJO STEP

		temp <- function(l){
			l.new  <- join.F(l, cur.supp, inner.supp, cur.par, inner.par)
			l.supp <- l.new$hsupp
			l.par  <- l.new$hpar
			return(psi.F(l.supp, l.par, xx))
			}
		
		# combination Armijo step (fastest)
		if (i==1){ # minimize temp(l)
		 	ll<-0:10/10
			ll.res<-sapply(ll, temp)
			ll.index<-which(ll.res==min(ll.res))[1]
			l<-ll[ll.index]
			}

		if(i>1){ # simple check
			l<-1			
			while (temp(l)-temp(0) > 0 ) { l=0.9*l }
			}

			hnew  <- join.F(l, cur.supp, inner.supp, cur.par, inner.par)
			supp <- hnew$hsupp
			hpar  <- hnew$hpar

       	if ((PROBLEM == FALSE)&(print == 1)) {	
                	cat("ARMIJO step lambda:", l,"\n")
			cat("\n")
		    	}

	} # end of outer loop

	if (i==max.outer.loop) PROBLEM <- TRUE

	psi <- psi.F(supp, hpar, xx)

      if ((PROBLEM == FALSE)&(print == 1)) {	
		# display support points with corresponding weights
		cat("\n")
            cat("FINAL RESULT: ",i, " iterations", "\n")
	 	cat("support: tau= ", supp$tau, " c= ", supp$constant, " eta= ", supp$eta, "\n")
	 	cat("weights: nu= ", hpar$nu, " alpha= ", hpar$alpha, " mu= ", hpar$mu, "\n")
		cat("log likelihood: ", -psi, "\n")
		}

      # if (PROBLEM == TRUE) cat("Problem occurred", "\n")

	return(list(mle=list(supp=supp, mix=hpar), llh=-psi*length(xx), conv=(PROBLEM==0)))

} # end of SRAlgF()


####################################################################
# support reduction step of the algorithm

supp.red.MLE.F <- function(supp, hpar, hparnew, cur.vec, xx, solve.tol) {

	m <- length(supp$eta)
	k <-  length(supp$tau)
	sr.PROBLEM <- FALSE

	test0 <- if(supp$constant != 0) hparnew$alpha / hpar$alpha  else 1  
	test1 <- if(k>0) hparnew$nu / hpar$nu else 1  
	test2 <- if(m>0) hparnew$mu / hpar$mu else 1  

	# vector of minimum [new weight] / [old weight] of three types of support points
	minTest <- c(test0, min(test1), min(test2))

 	# find which of eta, constant, tau to remove
	case <- which(minTest == min(minTest, na.rm = T) )[1] - 1
 
	# figure out lambdastar        
	if(case == 0) 	lambdastar <- 1 / (1 - test0)               
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

	hpartilde<-hpar.F(alphatilde, nutilde, mutilde)

	hpar0 <- min.psi.approx.F(supp, cur.vec, xx, solve.tol)
	if( hpar0$warn==TRUE ) { 
		sr.PROBLEM <- TRUE
		hpar0 <- hpar0$newpar
		} else {  hpar0<-hpar0$newpar  }
          
	# call recursively  
	if ( (min (hpar0$alpha, hpar0$nu, hpar0$mu) < 0)&&(sr.PROBLEM==FALSE) ) {
			supp.red.MLE.F(supp, hpartilde, hpar0, cur.vec, xx, solve.tol)} else ( list(weight = hpar0, support = supp, prob=sr.PROBLEM) )

}# end of supp.red.MLE.F

#################################################
# add two hazard functions

join.F <- function(l, supp1, supp2, hpar1, hpar2){

	# joint supports
	l.constant  	<-	if ((supp1$constant==1)|(supp2$constant==1)) 1
	l.tau 		<-	union(supp1$tau, supp2$tau)
	l.eta 		<-	union(supp1$eta, supp2$eta)
			
	# add mixtures
	h1		<- 	if(supp1$constant==1) hpar1$alpha else 0
	h2		<-	if(supp2$constant==1) hpar2$alpha else 0
	l.alpha	<-	(1-l)*h1+l*h2

	l.nu 		<-	rep(0, length(l.tau))
	if(length(l.nu)>0){   	
		for(i in 1:length(l.nu)){
			h1 <- if(l.tau[i] %in% supp1$tau) hpar1$nu[which(supp1$tau==l.tau[i])] else 0
			h2 <- if(l.tau[i] %in% supp2$tau) hpar2$nu[which(supp2$tau==l.tau[i])] else 0
			l.nu[i] <- (1-l)*h1 + l*h2
			}
		}  

	l.mu 		<-	rep(0, length(l.eta))  
	if(length(l.mu)>0){	
		for(i in 1:length(l.mu)){
			h1 <- if(l.eta[i] %in% supp1$eta) hpar1$mu[which(supp1$eta==l.eta[i])] else 0
			h2 <- if(l.eta[i] %in% supp2$eta) hpar2$mu[which(supp2$eta==l.eta[i])] else 0
			l.mu[i] <- (1-l)*h1 + l*h2
			}
		}	  

	# delete mixture and support if mixture is zero
	if (l.alpha == 0) l.constant = 0
	j1 <- which(l.nu == 0)
	if (length(j1) > 0){
		l.tau <- setdiff(l.tau, l.tau[j1])
		l.nu  <- setdiff(l.nu, l.nu[j1])
		}
	j2 <- which(l.mu == 0)
	if (length(j2) > 0){
		l.eta <- setdiff(l.eta, l.eta[j2])
		l.mu  <- setdiff(l.mu, l.mu[j2])
		}
				
	l.supp 	<-	supp.F(l.constant ,l.tau , l.eta)
	l.par  	<-	hpar.F(l.alpha, l.nu, l.mu)	

	return(list(hpar=l.par, hsupp=l.supp))
	} # end of join.F



#####################################################################
# directional derivative functions 

mat1.F <-function(GRID, xx, k){
	xx <- sort(xx)
	MM <- length(GRID)
	N  <- length(xx)
	mat1.f <- matrix(0, nrow=N, ncol=MM)
	for(i in 1:N){
		mat1.f[i,]	<- ifelse(GRID >= xx[i], (GRID-xx[i])^k, 0)
		}
	return(mat1.f)
	}

mat2.F <-function(GRID, xx, k){
	xx <- sort(xx)
	MM <- length(GRID)
	N  <- length(xx)
	mat2.f <- matrix(0, nrow=N, ncol=MM)
	for(i in 1:N){
		mat2.f[i,]	<- ifelse(GRID <= xx[i], (xx[i]-GRID)^k, 0)
		}
	return(mat2.f)
	}



######################################################
# Return the minimum directional derivatives (true and approximate), case (0, 1 or 2) and index.

dir.der.psi.F <- function(supp, hpar, GRID1, GRID2, xx) {

	MM1	<- 	length(GRID1)
	MM2	<-	length(GRID2)
	N	<- 	length(xx)

	h.vec		<-	sapply(xx, h, supp=supp, hpar=hpar)
	h.mat1 		<- 	matrix(h.vec, N, MM1)
	h.mat2 		<- 	matrix(h.vec, N, MM2) 

	mat.g1  	<-	mat1.F(GRID1, xx, 1)
	mat.G1  	<- 	matrix(GRID1^2, N, MM1, byrow=TRUE)/2-mat1.F(GRID1, xx, 2)/2
	mat.g2  	<-	mat2.F(GRID2, xx, 1) 
	mat.G2  	<- 	mat2.F(GRID2, xx, 2)/2

	
	DirDer1 <- as.vector((rep(1,N)/N) %*% (mat.G1-mat.g1/h.mat1))
	minDirDer1 <- min(DirDer1)		
	
	DirDer0 <- mean(xx-c(1/h.vec[1:(N-1)],0))  # MOD
	
	DirDer2 <- as.vector((rep(1,N)/N) %*% (mat.G2-diag(c(rep(1,N-1),0)) %*% mat.g2/h.mat2)) # MOD
	minDirDer2 <- min(DirDer2)

	temp <- c(minDirDer1, minDirDer2, DirDer0)
	minDirDer <- min(temp)

      	return(list(min=minDirDer))
 
} # End of function dir.der.psi.F


dir.der.psi.approx.F <- function(supp, hpar, cur.vec, GRID1, GRID2, xx) {

	MM1	<- 	length(GRID1)
	MM2	<-	length(GRID2)
	N	<- 	length(xx)
	
	h.vec	  	<-	sapply(xx, h, supp=supp, hpar=hpar)
	h.mat1		<-	matrix(h.vec, N, MM1)
	cur.mat1	<-	matrix(cur.vec, N, MM1)
	h.mat2		<-	matrix(h.vec, N, MM2)
	cur.mat2	<- 	matrix(cur.vec, N, MM2) 


	mat.G1 	<- 	matrix(GRID1^2, N, MM1, byrow=TRUE)/2-mat1.F(GRID1, xx, 2)/2
	mat.g1	<-	mat1.F(GRID1, xx, 1)
	mat.G2 	<- 	mat2.F(GRID2, xx, 2)/2  
	mat.g2	<-	mat2.F(GRID2, xx, 1)    

	DD1.t1	<-	as.vector((rep(1,N)/N) %*% (mat.G1-2*mat.g1/cur.mat1 + mat.g1*h.mat1/(cur.mat1)^2))
	DD1.t2	<-	as.vector((rep(1,N)/N) %*% (mat.g1/cur.mat1)^2)
	DD1		<- 	DD1.t1/sqrt(DD1.t2)
	DD1[which(DD1=="NaN")] 	<-Inf
	minDirDer1 	<- 	min(DD1, na.rm=TRUE)		
	
	DD0.t1	<-	mean( xx - c(2/cur.vec[1:(N-1)]-h.vec[1:(N-1)]/cur.vec[1:(N-1)]^2,0) ) # MOD 
	DD0.t2	<-	mean(c(1/cur.vec[1:(N-1)]^2,0))   # MOD 
	DirDer0 	<-	if(DD0.t2>0)  DD0.t1/sqrt(DD0.t2) else Inf

	DD2.t1	<-	as.vector((rep(1,N)/N) %*% (mat.G2-diag(c(rep(1,N-1),0)) %*% (2*mat.g2/cur.mat2 - mat.g2*h.mat2/(cur.mat2)^2)))  #MOD
	DD2.t2	<-	as.vector((c(rep(1,N-1),0)/N) %*% (mat.g2/cur.mat2)^2)  #MOD
	DD2		<- 	DD2.t1/sqrt(DD2.t2)
	DD2[which(DD2=="NaN")] 	<-Inf
	minDirDer2 	<- 	min(DD2, na.rm=TRUE)


	temp 		<- 	c(minDirDer1, minDirDer2, DirDer0)
	minDirDer 	<- 	min(temp)

	case 		<- 	(which( temp == minDirDer ))[1]

	if ( case == 1 ) index <- (which( DD1 == minDirDer ))[1] else if
	   ( case == 2 ) index <- (which( DD2 == minDirDer ))[1] else
	   { case <- 0 ; index <- 1 }

      	return(list(min=minDirDer, case=case, index=index))
 
} # end of dir.der.psi.approx.F


##################################################
#   return the A matrix and b vector for the function min.psi

ab.psi.F <-function(supp, cur.vec, xx){

	k <- length(supp$tau)
	m <- length(supp$eta)
	c <- supp$constant
	N <- length(xx)

	mat.g1  <-  mat1.F(supp$tau, xx, 1)
	mat.G1  <-  matrix(supp$tau^2, N, k, byrow=TRUE)/2-mat1.F(supp$tau, xx, 2)/2
	mat.g2  <-  mat2.F(supp$eta, xx, 1)  
	mat.G2  <-  mat2.F(supp$eta, xx, 2)/2

	D.Y1	<- if (k>0) diag(1/cur.vec) %*% mat.g1 			else numeric()
	D.Y0	<- if (c>0) matrix(c(1/cur.vec[1:(N-1)],0), N, 1) 	else numeric()   #MOD
	D.Y2	<- if (m>0) diag(c(1/cur.vec[1:(N-1)],0)) %*% mat.g2 	else numeric()   #MOD 
	

	D.Y	<- cbind(D.Y1, D.Y0, D.Y2)
	A	<- (t(D.Y) %*% D.Y)/N

	b1	<- if (k>0) - as.vector((rep(1,N)/N) %*% (mat.G1-2*(diag(1/cur.vec)%*% mat.g1))) 				else numeric()
	b0	<- if (c>0)   ((N-1)/N)*mean(2/cur.vec[1:(N-1)])-mean(xx) 							else numeric()  #MOD
	b2	<- if (m>0) - as.vector((rep(1,N)/N) %*% (mat.G2-2*(diag(c(1/cur.vec[1:(N-1)],0))%*% mat.g2))) 	else numeric()  #MOD

	b<-c(b1,b0,b2)

return(list(A=A,b=b))
} # end of ab.psi.F


##################################################
# Return the minimum of criterion for fixed support


min.psi.approx.F <- function(supp, cur.vec, xx, solve.tol){
	
	k <- length(supp$tau)
	m <- length(supp$eta)
	c <- supp$constant

	Ab<-ab.psi.F(supp, cur.vec, xx)
      

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
	}# end of min.psi.approx.F
 

####################################################
# Return the criterion for fixed support and weights

psi.F <- function(supp, hpar, xx){ 

	n 	<- length(xx) # MOD 
	T1	<- mean(sapply(xx, H, supp=supp, hpar=hpar))  # MOD
	T2	<- ((n-1)/n)*mean(sapply(xx[1:(n-1)], function(x, supp, hpar) log(h(x, supp, hpar)), supp=supp, hpar=hpar))  # MOD
	return(T1-T2)  # MOD

	} # end of psi.F
