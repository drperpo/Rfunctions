
allin <- function(X,R,death,Ftime,iter.max=50,mon=FALSE){
 ### X matrix of td covariates
 ### R matrix of fixed covariates
call <- match.call()
    n <- nrow(X); p <- ncol(X)
    q <- ncol(Ftime);v <- ncol(R)
 
    R <- apply(R, 2, function(x) {
        x - mean(x) })
   r <- 1
    eventtimes <- (1:n)[death == 1]
    k <- length(eventtimes)
    Ftime <- Ftime[eventtimes, ]
    dl <- rep(0,v)
    b <- matrix(rep(0, p), ncol = 1)
    gama <- matrix(runif(q)/10, ncol = 1)             
    gama[1, 1] <- 1
    delta <- 5
    iter <- 0
      k <- length(eventtimes)
	dm <- p*q+ncol(R)   
 while ( delta>1e-4 & iter < iter.max) {
        theta <- b %*% t(gama)
        iter <- iter + 1
        tmp <- sapply(eventtimes,sum.eventsf,X,R,Ftime,theta,dl,n,p,v,eventtimes)
        hi <- unlist(tmp[1, ])
          scores <- matrix(unlist(tmp[2,]),ncol=k) %*% rep.int(1, k) + 
               matrix(t(X[eventtimes,]) %*% Ftime,ncol=1)
      scores <- rbind(scores, matrix(unlist(tmp[3,]),ncol=k)%*%rep.int(1,k)+t(R[eventtimes,]) %*%rep.int(1,k))
      #information matrix
      hess <- -(matrix(matrix(unlist(tmp[4,]),ncol=k) %*% rep(1,k), ncol=dm, byrow=FALSE))
        vecu <- as.vector(scores)
        bc<- c(as.vector(b),dl)
        db <- kronecker(gama, diag(p))
	  db <- rbind(db,t(t(rep(0,v))) %*%rep(1,ncol(db)))
	  db<- cbind(db,t(t(rep(0,dm)))%*%rep(1,v))
        if(v>1){ db[(p*q+1):dm,(p+1):(p+v)] <- diag(rep(gama[1,1],v))}
        if(v==1){db[(p*q+1):dm,(p+1):(p+v)] <- gama[1,1]}
        d1 <- t(db) %*% vecu
        d2 <- -t(db) %*% hess %*% db
        sd2 <- solve(d2)
        bc <- bc + sd2 %*% d1
        b <- matrix(bc[1:p], ncol = r) 
        dl <- bc[(1+p):(p+v)]
        theta <- b %*% t(gama)
         bst <- sqrt(abs(diag(sd2)))
       tmp <- sapply(eventtimes,sum.eventsf,X,R,Ftime,theta,dl,n,p,v,eventtimes)
        hi <- unlist(tmp[1, ])
         scores <- matrix(unlist(tmp[2,]),ncol=k) %*% rep.int(1, k) + 
                matrix(t(X[eventtimes,]) %*% Ftime,ncol=1)
        scores <- rbind(scores, matrix(unlist(tmp[3,]),ncol=k)%*%rep.int(1,k)+t(R[eventtimes,]) %*%rep.int(1,k))
        #information matrix
        hess <- -(matrix(matrix(unlist(tmp[4,]),ncol=k) %*% rep(1,k), ncol=dm, byrow=FALSE))
        vecu <- scores
        parg <- as.vector(t(gama))
        dg <- kronecker(diag(q), b)
        dg <- rbind(dg,t(t(dl)) %*%rep(1,ncol(dg)))
        dg[(p*q+1):(p*q+v),2:q] <- 0
        d1 <- t(dg) %*% vecu
        d2 <- -t(dg) %*% hess %*% dg
	  sd2 <- solve(d2)
        parg <- parg +  sd2%*% d1
        gst <- sqrt(abs(diag(sd2)))
        gama <- matrix(parg, ncol = r, byrow = TRUE)
        gama <- gama/gama[1,1]
        theta2 <- b %*% t(gama)
         delta <-max(abs(theta-theta2))
if(mon==TRUE){cat(iter,"  ", round(dl[1],5),"dif",round(delta,6),"\n")}
}

        bc<- c(as.vector(b),dl)
        db <- kronecker(gama, diag(p))
	  db <- rbind(db,t(t(rep(0,v))) %*%rep(1,ncol(db)))
	  db<- cbind(db,t(t(rep(0,dm)))%*%rep(1,v))
        if(v>1){ db[(p*q+1):dm,(p+1):(p+v)] <- diag(rep(gama[1,1],v))}
        if(v==1){db[(p*q+1):dm,(p+1):(p+v)] <- gama[1,1]}
    	  parg <- as.vector(t(gama))
        dg <- kronecker(diag(q), b)
        dg <- rbind(dg,t(t(dl)) %*%rep(1,ncol(dg)))
        dg[(p*q+1):(p*q+v),2:q] <- 0
            dbg <- cbind(db, dg)
    d2 <- t(dbg) %*% hess %*% dbg
    covtheta <- dbg %*% ginv(d2, tol = 10^-15) %*% t(dbg)
    sterr <- sqrt(abs(diag(covtheta)))
    sterr <- round(sterr,5)



 pvals <- 1-pchisq((b/bst[1:p])^2 ,1)
  lik <- sum(log(hi)) + sum(X[eventtimes, ] * t(theta %*%t(Ftime))) +sum(R[eventtimes, ] %*%t(t(dl)) )
df <- (ncol(X)+ncol(Ftime)-1)+ncol(R)
aic <- -2*lik+2*df

out <- list(f.coef =dl,sterr=sterr,g.st=gst,p.vals=pvals, b.se=bst,theta=theta,betas = b,
gamma=gama,iter=iter.max,hess=hess,logL=lik,d2=d2,h0=hi,AIC=aic,call=call,df =df)
out}



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sum.eventsf <- function(i,X,R,Ftime,theta,dl,n,p,v,eventtimes){ 
#function to compute risk sets, weighted means,part of scores and information
    xj <- matrix(X[i:n,],nrow=(n-i+1))
    rj <- matrix(R[i:n,],nrow=(n-i+1))
    Fi <- Ftime[eventtimes==i,]
    rr <-exp(rj%*%dl+xj %*% theta %*% Fi) #relative risk
    hi <- 1/sum(rr)  #baseline hazard
    xmean <- hi*t(xj)%*%rr        
    rmean <- hi*t(rj)%*%rr  # first dev wrt d        
    xdif <- xj -matrix(rep(xmean,n-i+1),ncol=p, byrow=TRUE)
    rdif <- rj -matrix(rep(rmean,n-i+1),ncol=v, byrow=TRUE)
    XX <- t(xdif)%*% (xdif*matrix(rep(rr,p),ncol=p,byrow=F))
    FF <- Fi %*% t(Fi)
    RR <- t(rdif) %*%(rdif*matrix(rep(rr,v),ncol=v,byrow=FALSE)) #second der wrt dxd
    RB <- t(xdif) %*% (rdif * matrix(rep(rr, v), ncol = v, byrow = FALSE)) #second der wrt dx
    m1 <- matrix(c(kronecker(FF[,1],RB)),ncol=v, byrow = FALSE)
    j1 <- cbind(kronecker(FF,XX),m1)
      j2 <- rbind(j1,cbind(t(m1),RR))
        #  matrix(c(kronecker(FF[,1],RB),RR),nrow=v, byrow = TRUE))
list (hi,                      # hazard hi
     -xmean%*% t(Fi),            # part of term i of score 
         -rmean, 
         hi*j2)}  



tst <- function(fit,b1=NULL){
p <- length(fit$f.coef)
q <- length(fit$b.se)
b2 <- 1-pchisq((fit$f.coef[p]/fit$b.se[q])^2,1)
if(!is.null(b1))
b1 <- fit$p.vals[b1]
out <- list(p1 = b2,p2 = b1)
out
}


#### try more advanced models
sgn <- function(fit){
  p <- nrow(fit$theta)
  nd <- length(fit$b.se)
  1-pchisq((fit$f.coef/fit$b.se[(p+1):nd])^2,1)  
}

