expand.breakpoints <- function(dataset, breakpoints, index = "patnum", status
	 = "dead", tevent = "stop", Zvar = T)
{
# Expand <dataset> so that each individual has a row for each unique
# failure time that is <= his/her/its own failure time.
#
# Original version due to Kathy Russell, under the supervision
# of Kenneth Hess, Department of Biomathematics, 
# University of Texas M. D. Anderson Cancer Centre.
#
# Modified to handle problem with cbind when single columns are
# extracted from data frames, and to handle events at time = 0
# JHM, Nov 19, 1998.
#
# ERROR checking
# Require dataset to be of type data.frame
	if((missing(dataset) || !is.data.frame(dataset)
		)) stop("\n   Parameter, dataset, must be a data frame"
			)	
	# Require breakpoints to be a vector of length >= 1
	if((missing(breakpoints)) || (!(is.numeric(
		breakpoints))) || (!(is.vector(
		breakpoints))) || (length(breakpoints) < 
		1))
		stop("\n  Parameter, breakpoints, must be a numeric vector with at least
1 element"
			)
	varset <- names(dataset)
	covset <- varset
	lvset <- 1:length(varset)	
	# Require dataset to have unique names
	if(any(duplicated(varset))) {
		stop("\n Parameter, dataset, must have uniquely defined column  names"
			)
	}
# Require index to match exactly 1 dataset name
	if(length((indexloc <- lvset[varset == index])) !=
		1)
		stop("\n   Parameter, index, must be a character string matching exactly
1 name in the the first paramter, dataset"
			)
	covset <- covset[covset != index]	
	# Require status to match exactly 1 dataset name
	if(length((statusloc <- lvset[varset == status]
		)) != 1)
		stop("\n   Parameter, status, must be a character string matching exactly
1 name in the the first paramter, dataset"
			)
	covset <- covset[covset != status]	
	# Require tevent to match exactly 1 dataset name
	if(length((teventloc <- lvset[varset == tevent]
		)) != 1)
		stop("\n   Parameter, tevent, must be a character string matching exactly
1 name in the the first paramter, dataset"
			)
	covset <- covset[covset != tevent]	
	#*****************
#Begin
#*****************
	n <- nrow(dataset)
	temp.stop <- breakpoints
	if(breakpoints[1] > 0)
		temp.start <- c(0, breakpoints)
	else temp.start <- c(-1, breakpoints)	
	# JHM: Added 19.11.98
	temp.epoch <- 1:(length(breakpoints) + 1)
	n.break <- length(breakpoints) + 1
	temp.status <- rep(0, n.break)
	if(Zvar) {
		Zmat <- diag(rep(1, n.break))
		Z.name <- paste("Z", seq(1, n.break), 
			sep = "")
		Z.add <- vector()
	}
	epoch <- vector()
	ind <- vector()
	Tstart <- vector()
	Tstop <- vector()
	status2 <- vector()
	for(i in 1:n) {
		cur.status <- unlist(dataset[i, status]
			)
		cur.time <- unlist(dataset[i, tevent])
		curlength <- sum(breakpoints < cur.time
			)
		curlp1 <- curlength + 1
		ind <- c(ind, rep(i, curlp1))
		if(Zvar)
			Z.add <- rbind(Z.add, Zmat[1:
				curlp1,  ])
		Tstart <- c(Tstart, temp.start[1:curlp1
			])
		if(curlength > 0) {
			Tstop <- c(Tstop, temp.stop[1:
				curlength], cur.time)
			status2 <- c(status2, 
				temp.status[1:curlength
				], cur.status)
		}
		else {
			Tstop <- c(Tstop, cur.time)
			status2 <- c(status2, 
				cur.status)
		}
		epoch <- c(epoch, temp.epoch[1:curlp1])
	}
	new.dataset <- cbind(dataset[ind, index, drop
		 = F], Tstart, Tstop, status2, epoch, 
		dataset[ind, covset, drop = F])	
	# JHM 19.11.98: Added drop=F (x2)
	if(Zvar)
		new.dataset <- cbind(new.dataset, Z.add
			)
	nam <- c(index, "Tstart", "Tstop", status, 
		"epoch", covset)
	if(Zvar)
		nam <- c(nam, Z.name)
	dimnames(new.dataset) <- list(1:length(Tstop), 
		nam)
	return(new.dataset)
}

#breakpoints <- sort(unique(time))
#dat.exp<-expand.breakpoints(dat,breakpoints,index="id",status="death",tevent="time",Zvar=F)
