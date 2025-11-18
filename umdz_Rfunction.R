###############################################################################
#            Unimodal Density Estimation Using Bernstein Polynomials          #
#																	      	  #
#          Copyright(2013): Bradley C. Turnbull, NC State University and	  #	
#                           Dr. Sujit K. Ghosh, NC State University 		  # 
#																			  #	
#                         Version 2.0 - August 25, 2013					      #	
###############################################################################

#Packages needed
require(quadprog)
# HALLO
#Function to determine the optimal number of weights
#using the CN criterion
mOpt.CN <- function(tData, L){
	
#set starting point for m
	n <- length(tData)
	m <- floor(  n^(1/3) ) - 1
	
#set start value for while loop
	logratio <- 1
	
	while( logratio < sqrt(n) ){
		m <- m+1
		
	#construct D matrix using m value 
		B <- NULL
		for(k in 1:m){
			B <- cbind(B, pbeta(tData, shape1=k, shape2=m-k+1))
		}
		Dmat <- t(B) %*% L %*% B		
			
	#take spectral decomposition to find eigenvalues
		spec <- eigen(Dmat, symmetric=TRUE)
		d <- spec$values
		min.eigenValue <- max( min(d), 0 )
		max.eigenValue <- max(d)
		
		logratio <- log10(max.eigenValue) - log10(min.eigenValue)		
		}	
	m-1 #return number of weights
}

#Function to determine the index of the maximum weight
maxWeight <- function(m, Fn, lower, upper){
	#find max of weights
	maxplace <- which.max( Fn( ((1:m)/m)) - 
			        Fn( ((0:(m-1))/m) )  )
	}
	
#Function to generate the constraint matrix 
constraintMat <- function( m, maxplace){
	A <- suppressWarnings(   
	   rbind( rep(1,m), diag(rep(1,m)), 
	   matrix( rep( c(-1,1, rep(0,m - 1)) , maxplace-1), maxplace-1, m, byrow=TRUE), 
	   matrix( rep( c( rep (0,maxplace-1),1,-1,rep(0, m-maxplace)),m-maxplace),
	   m-maxplace,m,byrow=TRUE))
	   )
	Amat <- t(A)
	}

#Function to solve for the weight vector	
solveWeights <- function(m, Fn, lower, upper, Dmat, dvec){	
	#find the location of the maximum weight
		max.place <- maxWeight(m, Fn, lower, upper)
	#make the constraint matrix
		Amat <- constraintMat(m, max.place)	
	#make bvec vector of constraints
		bvec=c(1,rep(0,2*m-1))
	#find weights using solve.QP function
		w.hat = solve.QP(Dmat,dvec,Amat,bvec,meq=1)$solution
	#function to find max of an element and 0
		max0 <- function(x){max(x,0)}
	#make sure no weights are < 0
		w.hat <- sapply( w.hat, max0)
	#make sure the weights sum to 1
		wsum <- sum(w.hat)
		w.hat <- w.hat / wsum
}	
##################################################################################

#The main function which the user will call
umd <- function(data, bound.type="sd", fix.lower=NA, fix.upper=NA, crit="CN", 
															m=NA, warning=TRUE){
#Tranform the data to the [0,1] support

	#delta definition
		n <- length(data)
		delta <- sd(data)/ sqrt(n)
		sort.data <- sort(data)	
	
	#find lower and upper values for transformation
	#use given lower and upper values if given
		if(bound.type=="sd"){
			  lower.calc <- min(data) - delta
		}else if(bound.type=="Carv"){
			  lower.calc <- sort.data[1] - 
			  			(sort.data[2] - sort.data[1])/( (1-0.05)^-2 - 1)
		}	
			
		if( is.na(fix.lower) == TRUE){
			lower <- lower.calc
		}else{ 
			#if fix.lower is given take the larger of 
			#fix.lower and min data point with correction
			lower <- max( lower.calc, fix.lower ) 
			  #warn if fix.lower is less than the data
			  if( min(data) < fix.lower && warning == TRUE){
		         cat("WARNING: data contains values less than fix.lower \n")
			  }	
			}
	
	#upper bound		
		if(bound.type=="sd"){
				upper.calc <- max(data) + delta
		}else if(bound.type=="Carv"){
				upper.calc <- sort.data[n] + 
						(sort.data[n] - sort.data[n-1])/( (1-0.05)^-2 - 1)
		}	
		
		if( is.na(fix.upper) == TRUE){
			upper <- upper.calc	
		}else{ 
			#if fix.upper is given, take the smaller of 
			#fix.upper and max data point with correction
			upper <- min( upper.calc, fix.upper ) 
			  #warn if gUpper is more than data
			  if( max(data) > fix.upper && warning==TRUE){
			   cat("WARNING: data contains values greater than fix.upper \n")
		       }	
			}

	#Make the transformed data set
		tdata <- (data-lower)/(upper-lower)	
	
#Construct the L matrix and make vector of ecdf values
	Fn <- ecdf(tdata)
	ep <- 3/( 8*length(data) )	
	ecdf.vec <- Fn(tdata)	
	
	#construct L
	L <- diag( n / ( ( ecdf.vec + ep)*(1 + ep - ecdf.vec ) ) )
	
	
#Find the optimal number of weights, depends on the given crit
#then return the desired functions
if( crit == "CN"){
	mOpt <- mOpt.CN(tdata, L)
	
	if( is.na(m) == FALSE){
		
		if( mOpt < m && warning==TRUE){
		   cat("WARNING: given number of weights is larger than optimal number,\n"
		       , "\t \t optimal number of weights =", mOpt, "\n")
			}
		if( mOpt > m && warning==TRUE){
		   cat("WARNING: given number of weights is less than optimal number,\n"
		   	   ,"\t \t optimal number of weights =", mOpt, "\n")
			}
				
		#make B, Dmat, and dvec
			B <- NULL
			for(k in 1:m){
				B <- cbind(B, pbeta(tdata, shape1=k, shape2=m-k+1))
			}
			Dmat <- t(B) %*% L %*% B
			dvec <- t(B)%*% L %*% ecdf.vec
			
		#Make sure there are no zero eigen values	
			spec <- eigen(Dmat, symmetric=TRUE)
			d <- spec$values
			Q <- spec$vectors
		
		#find which values are < 10e-6, and set to smallest
		#eigen value >= 10e-6
			if( min(d) < 10e-6){
			tooSmall <- which( d < 10e-6)	
			d[tooSmall] <- d[ min(tooSmall) - 1]
		#Recreate pos. def. Dmat matrix 
			Dmat <- Q %*% diag(d) %*% t(Q)				
			}	
						
		} #if no m provided then use the optimal m
	else{ m <- mOpt 				
		#make B, Dmat, and dvec
			B <- NULL
			for(k in 1:m){
				B <- cbind(B, pbeta(tdata, shape1=k, shape2=m-k+1))
			}
			Dmat <- t(B) %*% L %*% B
			dvec <- t(B)%*% L %*% ecdf.vec

		#Make sure there are no zero eigen values	
			spec <- eigen(Dmat, symmetric=TRUE)
			d <- spec$values
			Q <- spec$vectors
		
		#find which values are < 10e-6, and set to smallest
		#eigen value >= 10e-6
			if( min(d) < 10e-6){
			tooSmall <- which( d < 10e-6)	
			d[tooSmall] <- d[ min(tooSmall) - 1]
		#Recreate pos. def. Dmat matrix 
			Dmat <- Q %*% diag(d) %*% t(Q)
			}
	}
		
#Solve for the weights
	weights <- solveWeights(m,Fn,lower,upper,Dmat,dvec)
	
#use the weights to create the 4 distribution functions
	#pdf
	dumd = function(x){ 
		mix.pdf = function(x){ifelse( (x-lower)*(upper-x) >= 0, 
			sum(weights*dbeta((x-lower)/(upper-lower),1:m,m:1))/(upper-lower),0)}
		sapply(x, mix.pdf)
		}
	#cdf
	pumd = function(x){ 
		mix.cdf= function(j){sum(weights*pbeta( (j-lower)/(upper-lower), 1:m, m:1)) }
		sapply(x, mix.cdf)
		}
	#rand generator
	rumd = function(n=1){
		rsample = function(){
			k <- sample(1:m,size=1,prob=weights)
			return(  lower+( rbeta(1, k, m-k+1)*(upper-lower) )   )  
		}
		replicate(n, rsample() )
	}
	#quantile function
	qumd = function(q){
		mix.quantile = function(q){ g = function(x){ pumd(x)-q }
			return( uniroot(g, interval=c(lower,upper) )$root ) 
			}
		sapply(q, mix.quantile)
	}

#Return list (ends function) 	
return(list(weights=weights,m.hat=m,dumd=dumd,pumd=pumd,rumd=rumd,qumd=qumd,
													lower=lower, upper=upper))
}
else if(crit %in% c("AIC","BIC")) {

#Find Optimal number of weights using AIC and BIC criterion	
	m.test <- seq( 2, ceiling(n/log(n)) )
	
	#lists and vectors to hold the values in the loop
	AICs <- rep(0, length(m.test) )
	BICs <- rep(0, length(m.test) )
	weights.list <- list()
	
	#loop through m values and calculate the AIC
	for( i in 1:length(m.test)){
		
		#make B, Dmat, and dvec for this m
		B <- NULL
		for(k in 1:m.test[i]){
			B <- cbind(B, pbeta(tdata, shape1=k, shape2=m.test[i]-k+1))
		}
		Dmat <- t(B) %*% L %*% B
		dvec <- t(B)%*% L %*% ecdf.vec
			
	#Make sure there are no zero eigen values	
		spec <- eigen(Dmat, symmetric=TRUE)
		d <- spec$values
		Q <- spec$vectors
		if( min(d) < 10e-6){
		tooSmall <- which( d < 10e-6)	
		d[tooSmall] <- d[ min(tooSmall) - 1]
		#Recreate pos. def. Dmat matrix 
		Dmat <- Q %*% diag(d) %*% t(Q)		
		}
		
		#Calculate weights
		weights.list[[i]] <- solveWeights(m.test[i],Fn,lower,upper,Dmat,dvec)
		
		#Make the pdf
		dumd.now = function(x){ 
		  mix.pdf = function(x){ifelse( (x-lower)*(upper-x) >= 0, 
		  sum(weights.list[[i]]*dbeta((x-lower)/(upper-lower),
		  1:m.test[i],m.test[i]:1))/(upper-lower),0)}
		  sapply(x, mix.pdf)
		}
		
		#Calculate AIC and BIC
		AICs[i]<- -2*sum( log(dumd.now(data))) + 2*(m.test[i]-2)
		BICs[i]<- -2*sum( log(dumd.now(data))) + log(n)*(m.test[i]-2)	
	}
	
	#Which AIC or BIC is the smallest
	if(crit == "AIC"){
		best.m <- which.min(AICs)
	}else if(crit == "BIC"){
		best.m <- which.min(BICs)
	}
	
	#Set weights and mOpt 
	weights <- weights.list[[best.m]]
	mOpt <- m.test[best.m]
 	
	#see if user supplied their own m value
	if( is.na(m) == FALSE){
		
		if( mOpt < m && warning==TRUE){
		   cat("WARNING: given number of weights is larger than optimal number,\n"
		       , "\t \t optimal number of weights =", mOpt, "\n")
			}
		if( mOpt > m && warning==TRUE){
		   cat("WARNING: given number of weights is less than optimal number,\n"
		   	   ,"\t \t optimal number of weights =", mOpt, "\n")
			}	
						
		#make B, Dmat, and dvec
			B <- NULL
			for(k in 1:m){
				B <- cbind(B, pbeta(tdata, shape1=k, shape2=m-k+1))
			}
			Dmat <- t(B) %*% L %*% B
			dvec <- t(B)%*% L %*% ecdf.vec
			
		#Make sure there are no zero eigen values	
			spec <- eigen(Dmat, symmetric=TRUE)
			d <- spec$values
			Q <- spec$vectors
		
		#find which values are < 10e-6, and set to smallest
		#eigen value >= 10e-6
			if( min(d) < 10e-6){
			tooSmall <- which( d < 10e-6)	
			d[tooSmall] <- d[ min(tooSmall) - 1]
			#Recreate pos. def. Dmat matrix 
			Dmat <- Q %*% diag(d) %*% t(Q)				
			}
			
		#Calculate weights	
		weights <- solveWeights(m,Fn,lower,upper,Dmat,dvec)						
	}
	else{
	#if no m is provided then just set m to mOpt
	#weights are already set
		m <- mOpt
	}	 
			
	#Calculate items to return to user:
	#pdf
	dumd = function(x){ 
	  mix.pdf = function(x){ifelse( (x-lower)*(upper-x) >= 0, 
	    sum(weights*dbeta((x-lower)/(upper-lower),1:m,m:1))/(upper-lower),0)}
	  sapply(x, mix.pdf)
	}
	#cdf
	pumd = function(x){ 
	  mix.cdf= function(j){sum(weights*pbeta( (j-lower)/(upper-lower), 1:m, m:1)) }
	  sapply(x, mix.cdf)
	}
	#rand generator
	rumd = function(n=1){
	  	rsample = function(){
		 k <- sample(1:m,size=1,prob=weights)
		 return(  lower+( rbeta(1, k, m-k+1)*(upper-lower) )   )  
		}
		replicate(n, rsample() )
	}
	#quantile function
	qumd = function(q){
	  mix.quantile = function(q){ g = function(x){ pumd(x)-q }
  	    return( uniroot(g, interval=c(lower,upper) )$root ) 
	  }
  	  sapply(q, mix.quantile)
	}	
	
#Return list (ends function) 	
return(list(weights=weights,m.hat=m,dumd=dumd,pumd=pumd,rumd=rumd,qumd=qumd,
			lower=lower, upper=upper))
}	
#end the function
}

#Function to calculate standard errors 
umd.se <- function(data, vals, bound.type="sd", fix.lower=NA, 
			fix.upper=NA, crit="CN", m=NA, warning=TRUE, boot.num=100){
#Tranform the data to the [0,1] support

	#delta definition
		n <- length(data)
		delta <- sd(data)/ sqrt(n)
		sort.data <- sort(data)	
	
	#find lower and upper values for transformation
	#use given lower and upper values if given
		if(bound.type=="sd"){
			  lower.calc <- min(data) - delta
		}else if(bound.type=="Carv"){
			  lower.calc <- sort.data[1] - 
			  			(sort.data[2] - sort.data[1])/( (1-0.05)^-2 - 1)
		}	
			
		if( is.na(fix.lower) == TRUE){
			lower <- lower.calc
		}else{ 
			#if fix.lower is given take the larger of 
			#fix.lower and min data point with correction
			lower <- max( lower.calc, fix.lower ) 
			  #warn if fix.lower is less than the data
			  if( min(data) < fix.lower && warning == TRUE){
		         cat("WARNING: data contains values less than fix.lower \n")
			  }	
			}
	
	#upper bound		
		if(bound.type=="sd"){
				upper.calc <- max(data) + delta
		}else if(bound.type=="Carv"){
				upper.calc <- sort.data[n] + 
						(sort.data[n] - sort.data[n-1])/( (1-0.05)^-2 - 1)
		}	
		
		if( is.na(fix.upper) == TRUE){
			upper <- upper.calc	
		}else{ 
			#if fix.upper is given, take the smaller of 
			#fix.upper and max data point with correction
			upper <- min( upper.calc, fix.upper ) 
			  #warn if gUpper is more than data
			  if( max(data) > fix.upper && warning==TRUE){
			   cat("WARNING: data contains values greater than fix.upper \n")
		       }	
			}

	#matrices to hold values
		dumd.est <- matrix(0, boot.num, length(vals) )
		pumd.est <- matrix(0, boot.num, length(vals) )
	
	for( i in 1:boot.num){
		cur.data <- sample(data, n, replace=T )
		
		cur.obj <- umd(cur.data, bound.type=bound.type, fix.lower=lower,
						fix.upper=upper, crit=crit, m=m, warning=F) 
	
		dumd.est[i,] <- cur.obj$dumd(vals)
		pumd.est[i,] <- cur.obj$pumd(vals)
	}

	dumd.se <- apply(dumd.est, 2, sd)
	pumd.se <- apply(pumd.est, 2, sd)

	#return(list(weights=weights,m.hat=m,dumd=dumd,pumd=pumd,rumd=rumd,qumd=qumd))
	return( list(dumd.se=dumd.se, pumd.se=pumd.se, vals=vals, lower=lower, upper=upper ))
}

