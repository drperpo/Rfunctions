
calc.h0 <- function(obj) {
	if (is(obj,"coxph")) {
      ha <- coxph.detail(obj)$haz
	ha*exp(-sum(obj$means*obj$coef))}
      else {
      ha <- obj$haz
      ha * exp(-sum(obj$means * obj$theta))}}
#################################################################################
coxvc <- function (formula = formula(data), Ft, rank, iter.max = 30, data = sys.parent()) 
{
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
    X <- model.matrix(Terms, m)[, -1, drop = FALSE]
    time <- Y[, 1]
    ord <- order(time)
    Ftime <- Ft[ord, ]
    n <- nrow(X)
    d <- Y[, 2]
    p <- ncol(X)
    q <- ncol(Ftime)
    r <- rank
    if (sum(Ft[, 1]) != n) 
        stop("The matrix of time functions should include a constant in the first column")
    nm <- attr(Terms, "term.labels")
    tmp <- paste("f", 1:(q - 1), "(t)", sep = "")

    X <- apply(X, 2, function(x) {
        x - mean(x)
    })
    
    means<-  matrix(apply(matrix(apply(X,2,function(x)x*Ft),ncol=p*q),2,mean),ncol=q,byrow=TRUE)
    
    if (r > p | r > q) 
        stop("Rank r cannot exceed the number of covariates or time functions", 
            "\n", "Full rank model for this data is: R=", min(p, 
                q))
    eventtimes <- (1:n)[d == 1]
    k <- length(eventtimes)
    Ftime <- Ftime[eventtimes, ]
    if (r == min(p, q)) {
        theta <- matrix(rep(0, p * q), ncol = q)
        fvc <- full.fit(eventtimes, X, Ftime, theta, iter.max, 
            n, p, q)
        nm <- c(nm, paste(rep(nm, q - 1), rep(tmp, each = p), 
            sep = ":"))
        fit <- list(call = call, theta = fvc$th, ster = fvc$se, 
            logL = fvc$lik, iter = fvc$iter, coef.names = nm, hazard = fvc$hi, means=means,time=Y[,1],death=Y[, 2],Ftime=Ft)
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
        frr <- rr.fit(eventtimes, X, Ftime, b, gama, iter.max, 
            n, p, q, r)
        nm <- c(nm, paste(rep(nm, q - 1), rep(tmp, each = p), 
            sep = ":"))
        fit <- list(call = call, b = frr$b, gama = frr$gama, 
            theta = frr$th, ster = frr$se, logL = frr$lik, iter = frr$iter, 
            rank = r, coef.names = nm, hazard=frr$hi,means=means,time=Y[,1],death=Y[, 2],Ftime=Ft,hessian=frr$hessian)
        class(fit) <- "coxrr"
        fit
    }
}

#############################################################################################
expand.haz <- function(haz,status,fun=c("baseline","cumulative")){
c <- status
c[status==0]  <- 1
c[status==1] <- 0
h2 <- NULL
it <- 1
k <- 1
while(it < (length(status)+1)){
    if (c[it]==0){h2[it] <- haz[it-k+1]}
    if (c[it]==1 && fun=="cumulative"){
        h2[it] <-  h2[it-1]
        k <- k+1}
    if (c[it]==1 && fun=="baseline"){
        h2[it] <-  0
        k <- k+1}
it <- it+1}
h2}      

full.fit <- 
function (eventtimes, X, Ftime, theta, iter.max, n, p, q) 
{
    delta <- 5
    iter <- 0
    k <- length(eventtimes)
    while (delta > 1e-06 && iter < iter.max) {
        th <- theta
        iter <- iter + 1
        tmp <- sapply(eventtimes, sumevents, X, Ftime, theta, 
            n, p, eventtimes)
        hi <- unlist(tmp[1, ])
        lik <- sum(log(hi)) + sum(X[eventtimes, ] * t(theta %*% 
            t(Ftime)))
        scores <- matrix(unlist(tmp[2, ]), ncol = k) %*% rep.int(1, 
            k) + matrix(t(X[eventtimes, ]) %*% Ftime, ncol = 1)
        hess <- -(matrix(matrix(unlist(tmp[3, ]), ncol = k) %*% 
            rep(1, k), ncol = p * q, byrow = FALSE))
        theta <- matrix(matrix(theta, ncol = 1) - solve(hess, 
            matrix(t(scores), ncol = 1)), ncol = q)
        delta <- max(abs(theta - th))
        if (iter > iter.max) 
            stop("Likelihood did not converge after ", iter, 
                " iterations")
    }
    sterr <- round(matrix(sqrt(abs(diag(solve(hess)))), nrow = p), 
        5)
    outlist <- list(th = round(theta, 5), se = sterr, lik = lik, 
        iter = iter, hi = hi)
    outlist
}
##############################################################################################

plotcoxvc <- function(object, fun=c("survival","effects"),ylab="",xlab="",lwd=1,col=1,lty=1,...){
    if(!inherits(object, "coxvc")  && !inherits(object,"coxrr"))
        stop("First arg must be the result of coxvc")

    if (fun=="survival"){
    haz <- object$haz * exp(-sum(object$means * object$theta))
    surv <- exp(-cumsum(haz*(exp( sum(object$means * object$theta)))))    
    plot(object$time[death==1],surv,ylim=c(0,1),type='s',ylab=ylab,xlab=xlab,lwd=lwd,col=col,lty=lty)}
     if (fun=="effects"){
     th <- object$theta
     p <- nrow(th)
     nms <- object$coef.names[1:nrow(th)]
     Ft <- object$Ftime[death==1,]
     ef <- th%*%t(Ft)
   plot(time[death==1],ef[1,],ylim=c(floor(min(ef)),ceiling(max(ef))),type='l',lwd=lwd,col=col,ylab=ylab,xlab=xlab,lty=lty)     
     i <- 1
                while (i < p) {
                i <- i+1
                lines(time[death==1],ef[i,],lty=i,col=i,lwd=lwd)}
                legend(min(time[death==1]),ceiling(max(ef)),nms,lty=1:p,col=1:p)}

    }
##############################################################################################
print.coxrr <- function(x,...) {
    cat("call:"); cat("\n")
    print(x$call);cat("\n")
    coef <- as.vector(x$theta)
    se <- as.vector(x$ster)
    tmp <- cbind(coef,exp(coef),se)
    dimnames(tmp) <- list(x$coef.names,c("coef","exp(coef)","se(coef)"))
    prmatrix(tmp)
    cat("\n")
    cat("log-likelihood= ",x$logL,", ","Rank=",x$rank,"\n")
    cat("algorithm converged in", x$iter ,"iterations","\n")
   cat("\n")
    cat("Beta :"); cat("\n"); prmatrix(round(x$b,4))   ;cat("\n")
    cat("Gamma:"); cat("\n"); prmatrix(round(x$gama,4));cat("\n") }
############################################################################################
print.coxvc <- function (x, digits =max(options()$digits-4,3),...) 
{
    cat("call:")
    cat("\n")
    print(x$call)
    cat("\n")
    coef <- as.vector(x$theta)
    se <- as.vector(x$ster)
    tmp <- cbind( round(coef,4), round(exp(coef),3), round(se,3), round(coef/se,4),   signif(1 - pchisq((coef/se)^2, 1), digits-1 ))
    dimnames(tmp) <- list(x$coef.names, c("coef", "exp(coef)", 
        "se(coef)", "z", "p"))
    prmatrix(tmp)
    cat("\n")
    cat("log-likelihood= ", x$logL)
    cat("\n")
    cat("\n")
    cat("algorithm converged in", x$iter, "iterations")
}       
############################################################################################
 rr.fit <- 
function (eventtimes, X, Ftime, b, gama, iter.max, n, p, q, r) 
{
    delta <- 5
    iter <- 0
    k <- length(eventtimes)
    while (delta > 1e-06 && iter < iter.max) {
        theta <- b %*% t(gama)
        iter <- iter + 1
        tmp <- sapply(eventtimes, sumevents, X, Ftime, theta, 
            n, p, eventtimes)
        hi <- unlist(tmp[1, ])
        lik <- sum(log(hi)) + sum(X[eventtimes, ] * t(theta %*% 
            t(Ftime)))
        scores <- matrix(unlist(tmp[2, ]), ncol = k) %*% rep.int(1, 
            k) + matrix(t(X[eventtimes, ]) %*% Ftime, ncol = 1)
        hess <- -(matrix(matrix(unlist(tmp[3, ]), ncol = k) %*% 
            rep(1, k), ncol = p * q, byrow = FALSE))
        vecu <- as.vector(scores)
        parb <- as.vector(b)
        db <- kronecker(gama, diag(p))
        d1 <- t(db) %*% vecu
        d2 <- -t(db) %*% hess %*% db
        parb <- parb + solve(d2) %*% d1
        b <- matrix(parb, ncol = r)
        theta <- b %*% t(gama)
        tmp <- sapply(eventtimes, sumevents, X, Ftime, theta, 
            n, p, eventtimes)
        hi <- unlist(tmp[1, ])
        lik2 <- sum(log(hi)) + sum(X[eventtimes, ] * t(theta %*% 
            t(Ftime)))
        scores <- matrix(unlist(tmp[2, ]), ncol = k) %*% rep.int(1, 
            k) + matrix(t(X[eventtimes, ]) %*% Ftime, ncol = 1)
        hess <- -(matrix(matrix(unlist(tmp[3, ]), ncol = k) %*% 
            rep(1, k), ncol = p * q, byrow = FALSE))
        vecu <- scores
        parg <- as.vector(t(gama))
        dg <- kronecker(diag(q), b)
        d1 <- t(dg) %*% vecu
        d2 <- -t(dg) %*% hess %*% dg
        parg <- parg + solve(d2) %*% d1
        gama <- matrix(parg, ncol = r, byrow = TRUE)
	gama <- gama/gama[1,1]
        delta <- lik2 - lik
    }
    parb <- as.vector(b)
    parg <- as.vector(t(gama))
    parbg <- c(parb, parg)
    db <- kronecker(gama, diag(p))
    dg <- kronecker(diag(q), b)
    dbg <- cbind(db, dg)
    d2 <- t(dbg) %*% hess %*% dbg
    covtheta <- dbg %*% ginv(d2, tol = 10^-15) %*% t(dbg)
    sterr <- sqrt(abs(diag(covtheta)))
    sterr <- round(matrix(sterr, nrow = p), 5)
    outlist <- list(b = b, gama = gama, th = round(theta, 6), 
        se = sterr, lik = lik2, iter = iter, hi =hi,hessian=hess)
}

###########################################################################################
       
sumevents <- function(i,X,Ftime,theta,n,p,eventtimes){ 
#function to compute risk sets, weighted means,part of scores and information
    xj <- matrix(X[i:n,],nrow=(n-i+1))
    Fi <- Ftime[eventtimes==i,]
    rr <-exp(xj %*% theta %*% Fi) #relative risk
    hi <- 1/sum(rr)  #baseline hazard
    xmean <- hi*t(xj)%*%rr        
    xdif <- xj -matrix(rep(xmean,n-i+1),ncol=p, byrow=TRUE)
    XX <- t(xdif)%*% (xdif*matrix(rep(rr,p),ncol=p,byrow=FALSE))
    FF <- Fi %*% t(Fi)
    list (hi,                      # hazard hi
         -xmean%*% t(Fi),            # part of term i of score 
         hi*kronecker(FF,XX))} 

summary.coxrr <- function(object,...) {
    cat("call:"); cat("\n")
    print(object$call);cat("\n")
     cat("Beta :"); cat("\n"); prmatrix(round(object$b,4))   ;cat("\n")
    cat("Gamma:"); cat("\n"); prmatrix(round(object$gama,4));cat("\n") }

