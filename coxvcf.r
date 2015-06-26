## fits reduced and full rank model with 
## time dependent and fixed effects
##  two ways of computing  st.errors, not sure yet which one
## is correct (probably second)
##
library(survival)
library(MASS)
coxvcf <- function (formula = formula(data), X,  Ft, rank, iter.max = 30, data=sys.parent(), theta0) 
{
# formula = survival type Surv(time,death)~x1+x2
# X       = matrix of covariates that have td effects
# Ft      = matrix of time functions, first column should be ones
# rank    = rank of the reduced rank model
# theta0 =  starting values of theta if known

    nar <- nargs()
    call <- match.call()
    m <- match.call(expand = FALSE)
    temp <- c("", "formula", "data")
    m <- m[match(temp, names(m), nomatch = 0)]
    Terms <- terms(formula, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())
    if (NROW(m) == 0) 
        stop("No (non-missing) observations")
    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")) 
        stop("Response must be a survival object")
    R <- model.matrix(Terms, m)[, -1, drop = FALSE]
    time <- Y[, 1]
    ord <- order(time)
    Ftime <- Ft[ord, ]
    n <- nrow(X)
    d <- Y[, 2]
    p <- ncol(X)
    q <- ncol(Ftime)
    v <- ncol(R)
    r <- rank
  
  if (sum(Ft[, 1]) != n) 
        stop("The matrix of time functions should include a constant in the first column")
    nm <- attr(Terms, "term.labels")
    tmp <- paste("f", 1:(q - 1), "(t)", sep = "")

    R <- apply(R, 2, function(x) {
        x - mean(x)
    })
    

     if (r > p | r > q) 
        stop("Rank r cannot exceed the number of covariates or time functions", 
            "\n", "Full rank model for this data is: R=", min(p, 
                q))
 
    eventtimes <- (1:n)[d == 1]
    k <- length(eventtimes)
    Ftime <- Ftime[eventtimes, ]
    dl <- rep(0,v)
    
   if (r == min(p, q)) {
        theta <-matrix(rep(0, p * q), ncol = q)
if(nar==6){theta <- theta0}
         fvc <- full.fitf(eventtimes, X, R, Ftime, theta, dl, iter.max,    n, p, q,v)
        nm <- c(nm, paste(rep(nm, q - 1), rep(tmp, each = p),             sep = ":"))

   aic1 <- fvc$loglik-r*(p+q-r)-v
  aic2 <- 2* (r*(p+q-r)-v)-2*fvc$loglik
 aicc <- aic2+ 2*(r*(p+q-r)-v)*(r*(p+q-r)-v+1)/(n-(r*(p+q-r)-v)-1)
        fit <- list(call = call, theta = fvc$th, ster = fvc$se,  f.coef=fvc$f.coef
, iter = fvc$iter, coef.names = nm, hazard = fvc$hi, time=Y[,1],death=Y[, 2],Ftime=Ft, loglik=fvc$loglik, aic1=aic1,aic2=aic2,aicc=aicc)
        class(fit) <- "coxvc"
        fit
}    

    else {
       b <- matrix(rep(0, r * p), ncol = r)
        gama <- matrix(runif(q * r)/10, ncol = r)
        if (r == 1) {
            gama[1, 1] <- 1
        }
        else {
            gama[, 1] <- 1
        }
      frr <- rr.fitf(eventtimes, X, R,Ftime, b, gama,dl, iter.max,n, p, q,v, r)

   nm <- c(nm, paste(rep(nm, q - 1), rep(tmp, each = p), 
            sep = ":"))
  aic1 <- frr$loglik-r*(p+q-r)-v
  aic2 <- 2* (r*(p+q-r)-v)-2*frr$loglik
 aicc <- aic2+ 2*(r*(p+q-r)-v)*(r*(p+q-r)-v+1)/(n-(r*(p+q-r)-v)-1)
        fit <- list(call = call, b = frr$b, gama = frr$gama, f.coef= frr$f.coef, ster=frr$se, s2=frr$s2,
            theta = frr$th,  rank = r, coef.names = nm, hazard=frr$hi,time=Y[,1],death=Y[, 2],Ftime=Ft,loglik=frr$loglik,aic1=aic1,aic2=aic2,aicc=aicc, hess=frr$hess)
        class(fit) <- "coxrr"
        fit
    }            
}




 rr.fitf <- 
function (eventtimes, X, R, Ftime, b, gama, dl, iter.max, n, p, q, v,r) 
{


    delta <- 5
    iter <- 0
    k <- length(eventtimes)
dm <- p*q+ncol(R)   
 while (delta > 1e-05 && iter < iter.max) {
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
        parb <- as.vector(b)
        db <- kronecker(gama, diag(p))
        db <- rbind(db,t(t(dl)) %*%rep(1,ncol(db)))
        d1 <- t(db) %*% vecu
        d2 <- -t(db) %*% hess %*% db
        parb <- parb + solve(d2) %*% d1
        b <- matrix(parb, ncol = r)
        theta <- b %*% t(gama)

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
        d1 <- t(dg) %*% vecu
        d2 <- -t(dg) %*% hess %*% dg
        parg <- parg + solve(d2) %*% d1
        gama <- matrix(parg, ncol = r, byrow = TRUE)
     theta2 <- b %*% t(gama)

delta <-max(abs(theta-theta2))
cat(iter,"  ", delta,"\n")


}

cri <- 4
iter=0
while( cri > 0.00001 & iter<30){
iter <- iter+1
 tmp <- sapply(eventtimes,sum.eventsf,X,R,Ftime,theta,dl,n,p,v,eventtimes)
 hi <- unlist(tmp[1, ])
 scores <- matrix(unlist(tmp[3,]),ncol=k)%*%rep.int(1,k)+t(R[eventtimes,]) %*%rep.int(1,k)
      #information matrix
 hess <- -(matrix(matrix(unlist(tmp[4,]),ncol=k) %*% rep(1,k), ncol=dm, byrow=FALSE))
hs <- hess[(p*q+1):(p*q+v),(p*q+1):(p*q+v)]
dl2<- dl-solve(hs,scores)
cri <- max(abs(dl-dl2))
dl <- dl2

cat(iter, "   ", cri, "\n")}

    parb <- as.vector(b)
    parg <- as.vector(t(gama))
    parbg <- c(parb, parg)
    db <- kronecker(gama, diag(p))
    dg <- kronecker(diag(q), b)
    dbg <- cbind(db, dg)
    d2 <- t(dbg) %*% hess[1:(p*q),1:(p*q)] %*% dbg
    covtheta <- dbg %*% ginv(d2, tol = 10^-15) %*% t(dbg)
    sterr1 <- sqrt(abs(diag(covtheta)))
    sterr1 <- round(matrix(sterr1, nrow = p), 5)
    ster2 <- sqrt(abs(diag(solve(hs))))
     ster2 <- round(ster2,3)
    sterr <- c(sterr1,ster2)

dl <- as.vector(dl)

 lik <- sum(log(hi)) + sum(X[eventtimes, ] * t(theta %*%t(Ftime)))+sum(R[eventtimes,] *dl)
   s2 <-   round(matrix(sqrt(abs(diag(solve(hess)))), nrow = dm),  3)
  outlist <- list(th = round(theta, 5),  f.coef =dl, b=b, gama=gama, se = sterr, s2= s2, hi = hi, loglik=lik, hess=d2)
    outlist
}


full.fitf <- function (eventtimes, X, R,  Ftime, theta,dl,  iter.max, n, p, q,v) 
{
delta <- 5; iter <-0; k <- length(eventtimes)
dm <- p*q+ncol(R)
  while(delta>1e-5 & iter < iter.max){
      th <- theta 
      iter <- iter+1
      #call function sum.events to compute risk sets etc
      tmp <- sapply(eventtimes,sum.eventsf,X,R,Ftime,theta,dl,n,p,v,eventtimes)
      hi <- unlist(tmp[1,])

      scores <- matrix(unlist(tmp[2,]),ncol=k) %*% rep.int(1, k) + 
                matrix(t(X[eventtimes,]) %*% Ftime,ncol=1)
      scores <- rbind(scores, matrix(unlist(tmp[3,]),ncol=k)%*%rep.int(1,k)+t(R[eventtimes,]) %*%rep.int(1,k))
      #information matrix
      hess <- -(matrix(matrix(unlist(tmp[4,]),ncol=k) %*% rep(1,k), ncol=dm, byrow=F))
             
      coef <- c(matrix(theta,ncol=1),dl)-solve(hess,matrix(t(scores),ncol=1))
      theta <- matrix(coef[1:(p*q)],ncol=q)
      dl <- coef[(p*q+1):dm]
      delta <- max(abs(theta-th))
      cat(iter,"   ", delta,"\n")
              if (iter>iter.max) 
           stop("Likelihood did not converge after ",iter," iterations")}

  sterr <- round(matrix(sqrt(abs(diag(solve(hess)))), nrow = dm), 
        3)

 lik <- sum(log(hi)) + sum(X[eventtimes, ] * t(theta %*%t(Ftime)))+sum(R[eventtimes,] *dl)
    outlist <- list(th = round(theta, 5), se = sterr, f.coef = dl,
        iter = iter, hi = hi,loglik=lik)
    outlist
}

sum.eventsf2 <- function(i,X,R,Ftime,theta,dl,n,p,v,eventtimes){ 
#function to compute risk sets, weighted means,part of scores and information
    xj <- matrix(X[i:n,],nrow=(n-i+1))
    rj <- matrix(R[i:n,],nrow=(n-i+1))
    Fi <- Ftime[eventtimes==i,]
    rr <-exp(rj%*%dl+xj %*% theta %*% Fi) #relative risk
    hi <- 1/sum(rr)  #baseline hazard
    xmean <- hi*t(xj)%*%rr        
    rmean <- hi*t(rj)%*%rr  # first dev wrt d        
    xdif <- xj -matrix(rep(xmean,n-i+1),ncol=p, byrow=T)
    rdif <- rj -matrix(rep(rmean,n-i+1),ncol=v, byrow=T)
    XX <- t(xdif)%*% (xdif*matrix(rep(rr,p),ncol=p,byrow=F))
    FF <- Fi %*% t(Fi)
    RR <- t(rdif) %*%(rdif*matrix(rep(rr,v),ncol=v,byrow=F)) #second der wrt dxd
    RB <- t(xdif) %*% (rdif * matrix(rep(rr, v), ncol = v, byrow = FALSE)) #second der wrt dx
    
    j1 <- cbind(kronecker(FF,XX),kronecker(FF[,1],RB))
    j2 <- rbind(j1,matrix(c(kronecker(FF[,1],RB),RR),nrow=v, byrow = FALSE))
    list (hi,                      # hazard hi
     -xmean%*% t(Fi),            # part of term i of score 
         -rmean, 
         hi*j2)} 

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



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

