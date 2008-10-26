#############################################################
# find the LSE of a convex hazard using bisection + support reduction


convexLSE <- function(x, TT, M = 100, GRIDLESS=0, tol =1e-5, tol.SR=1e-8, type=1, max.loop=100, max.loop.SR=250, range=c(0, TT), solve.tol=.Machine$double.eps) 
{

	if(type==1) delta.F<- function(x, xold) sum((x-min(x))^2)
   	if(type==2) delta.F<- function(x, xold) sqrt(sum((x - xold)^2))
   	if(type==3) delta.F<- function(x, xold) sort(x)[2]+sort(x)[3]-2*sort(x)[1]
	
	zzz <- file("convexLSE.out", "w")
	cat("phi","antimode","conv","\n", file=zzz)
	
	m 		<- (range[2]-range[1])*(0:4)/4+range[1]
	phi.vec 	<- rep(0,5)
	temp.vec 	<- list()
	

	for(i in 1:5){
		temp 		 <- 	srLSE(x, m[i], TT, M, GRIDLESS, tol, max.loop.SR, print=0, solve.tol)
		temp.vec[i]  <-   list(temp)
		phi.vec[i] <- 	temp$ls
		cat(phi.vec[i], m[i], temp$conv,"\n", file=zzz)
		}

	phi.delta <- delta.F(phi.vec, rep(10,5))


	for(i in 1:max.loop)
	{
	
	if (phi.delta < tol) break

	j0 <- which(phi.vec == min(phi.vec))[1]
	if(length(j0)==5)  j0 <- j0[3]
	if(length(j0)==4)	 j0 <- j0[2]	
	if(length(j0)==3)  j0 <- j0[2]
	if(length(j0)==2)  j0 <- j0[1]
	# if length(j0)==1 do nothing

	if(j0==1) j0 <- 2
	if(j0==5) j0 <- 4	

	index <- (j0-1):(j0+1)

	mtilde <- rep(0,5)
	mtilde[c(1,3,5)] <- m[index]
	mtilde[2] <- (mtilde[1] + mtilde[3])/2
	mtilde[4] <- (mtilde[3] + mtilde[5])/2
	m <- mtilde

	phi.vec.old 		<- phi.vec 
	phi.vec[c(1,3,5)] 	<- phi.vec.old[index]
	temp.vec[c(1,3,5)]	<- temp.vec[index] 

	for(i in 1:2){
		temp		 	<- 	srLSE(x, m[2*i], TT, M, GRIDLESS, tol, max.loop.SR, print=0, solve.tol)
		temp.vec[2*i]     <-	list(temp)
		phi.vec[2*i] 	<- 	temp$ls
		cat(phi.vec[2*i], m[2*i], temp$conv,"\n", file=zzz)
		}

	phi.delta <- delta.F(phi.vec, phi.vec.old)

	} # end loop	

	close(zzz)

	j0 <- which(phi.vec == min(phi.vec))
	if(length(j0)==5)  j0 <- j0[3]
	if(length(j0)==4)	 j0 <- j0[2]	
	if(length(j0)==3)  j0 <- j0[2]
	if(length(j0)==2)  j0 <- j0[1]
	# if length(j0)==1 do nothing

	res <- temp.vec[[j0]]

	return(list(lse=list(supp=res$lse$supp, mix=res$lse$mix), ls=res$ls, antimode=m[j0], conv=(i<max.loop), iter=i))

} # end of convexLSE function

