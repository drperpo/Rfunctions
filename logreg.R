logreg <- function(y,x,X0,kapa,...){
	fit <- list()
	m <- nrow(X0)
	n <- ncol(X0)
	X <- cbind(rep(1,m),X0)
	n <- n+1
	je <- nargs()
	 	if (je<4){
			kapa <- 0}
	if (je < 5){
    p0 <- (sum(y) + 1) / (sum(x) + 2)
    a <-  rep(0,n)
    a[1] <- log(p0 / (1 - p0))}
	
	R <- kapa * (outer(1:n,1:n,"==")+0);
	R[1, 1] <- 0
	
	rp <-  1
	it <- 0
		while (rp>1e-04){
    		it <-it + 1
    		z <- X %*% a
          a1 <- a
    		p <- 1 / (1 + exp(-z))
    		mu <- x * p
    		w <- mu * (1 - p)
    		W <- rep(w, n)
    		G0 <- t(X) %*% (W * X)
    		G <- G0 + R
    		a <- solve( G, (t(X) %*% (y - mu + w * z)))
    		da <- max(abs(a-a1))
    		rp <- da / max(abs(a))  
			cat(da,"\n")}
fit$p <-p
fit$a <- a
fit$G <-G
fit$G0 <- G0
fit$z <- z
fit}


 bindevnew <-function(y, t, p){
d <- -2 * sum(y * log(p) + (t - y)* log(1 - p))}


kapa <- 5
ff <- logreg(y,t,X,kapa)

#if nargout > 2
dev <- bindevnew(y, t, ff$p)
tr <- sum(diag(solve(ff$G) %*% ff$G0))
aic <- dev + 2 * tr
#end

