scatterplot <- function(formula,data,horizon,plot=TRUE,xlab)
{
    cox1 <- coxph(formula, data)
    y <- cox1[["y"]]
    p <- ncol(y)
    time <- y[,p-1]
    status <- y[,p]
    x <- cox1$linear.predictors
    if (missing(horizon)) horizon <- max(time[status==1])*1.05
    ylim <- c(0,1.05*horizon) # extra space to print R-square
    imputed <- time
    cens <- which(status==0)
    events <- which(status==1)
    ncens <- length(cens)
    if (missing(xlab)) xlab <- "x"
    cx <- coxph(Surv(time,status) ~ x, method="breslow")
    for (i in 1:ncens) {
        xi <- x[cens[i]]
        nd <- data.frame(x=xi)
        sf <- survfit(cx,newdata=nd)
        sf <- data.frame(time=sf$time,surv=sf$surv)
        sf <- sf[!duplicated(sf$surv),]
        if (imputed[cens[i]]>max(sf$time))
            imputed[cens[i]] <- horizon
        if (imputed[cens[i]]<max(sf$time)) {
            sf1 <- sf[sf$time <= imputed[cens[i]],]
            sf2 <- sf[sf$time >= imputed[cens[i]],]
            sf2$surv <- 1-sf2$surv/min(sf1$surv)
            if (sf2$surv[1]==0) sf2 <- sf2[-1,]
            rand <- runif(1)
            survtimes <- c(0,sf2$surv,1)
            idx <- as.numeric(cut(rand,survtimes))
            if (idx>nrow(sf2)) imputed[cens[i]] <- horizon
            if (idx<=nrow(sf2)) imputed[cens[i]] <- sf2$time[idx]
        }
    }
    res <- data.frame(x=x, imputed=imputed)
    attr(res, "horizon") <- horizon
    if (plot) {
        plot(x, imputed, type="n", ylim=ylim, xlab=xlab, ylab="(Imputed) survival time")
        points(x[cens], imputed[cens])
        points(x[events], imputed[events], pch=16)
        abline(lsfit(x,imputed),lty=2)
        lmova <- lm(imputed ~ x)
        r2 <- summary(lmova)$r.squared
        text(max(x),1.05*horizon,paste("R-squared =",round(r2,3)),adj=1)
    }
    return(res)
}




pew <-
  function (time, status, tsurv, survmat, tcens, censmat, width,
            FUN = c("KL", "Brier"), tout)
  {
    ### Data needs to be ordered according to time (asc) and status (desc)
    ord <- order(time, -status)
    time <- time[ord]
    status <- status[ord]
    survmat <- survmat[,ord]
    censmat <- censmat[,ord]
    ### Both tsurv and tcens need to be sorted and need to start with 0
    if (any(!(tsurv==sort(tsurv)))) stop("argument 'tsurv' needs to be sorted")
    if (any(!(tcens==sort(tcens)))) stop("argument 'tcens' needs to be sorted")
    if (min(tsurv)<0) stop("no negative values allowed for 'tsurv'")
    if (min(tcens)<0) stop("no negative values allowed for 'tcens'")
    if (min(tsurv)>0) {
      tsurv <- c(0,tsurv)
      survmat <- rbind(rep(1,ncol(survmat)),survmat)
    }
    if (min(tcens)>0) {
      tcens <- c(0,tcens)
      censmat <- rbind(rep(1,ncol(censmat)),censmat)
    }
    ### Prepare for call to prederrw
    FUN <- match.arg(FUN)
    FUNn <- 2
    if (FUN == "Brier")
      FUNn <- 1
    n <- length(time)
    nsurv <- length(tsurv)
    ncens <- length(tcens)
    nout <- length(tout)
    res <- .C("prederrw", sn = as.integer(n), time = as.double(time),
              status = as.integer(status), snsurv = as.integer(nsurv),
              sncens = as.integer(ncens), snout = as.integer(nout),
              tsurv = as.double(tsurv), survmat = as.double(survmat),
              tcens = as.double(tcens), censmat = as.double(censmat),
              w = as.double(width), tout = as.double(tout), sFUN = as.integer(FUNn),
              err = double(nout), work = double(n))
    res <- data.frame(time = tout, Err = res$err)
    attr(res, "Score") <- FUN
    attr(res, "width") <- width
    return(res)
  }

pewcox <- function(formula, censformula, width, data, censdata,
                   FUN = c("KL", "Brier"), tout, CV = FALSE, progress = FALSE)
{
  cox1 <- coxph(formula, data)
  y <- cox1[["y"]]
  p <- ncol(y)
  time <- y[, p - 1]
  status <- y[, p]
  n <- length(time)
  if (nrow(data) != n)
    stop("missing data in time or status not allowed")
  ord <- order(time, -status)
  time <- time[ord]
  status <- status[ord]
  tt <- sort(unique(time[status == 1]))
  nt <- length(tt)
  if (!CV) {
    progress <- FALSE
    x <- cox1$linear.predictors[ord]
    cox1 <- coxph(Surv(time, status) ~ x)
    # if no covariates, then use Kalbfleisch-Prentice type in survfit
    if (sum(x^2)==0)
      sf <- survfit(cox1, newdata = data.frame(x = x), type="kalbfl")
    else
      sf <- survfit(cox1, newdata = data.frame(x = x))
    tt <- sf$time
    survmat <- sf$surv
  }
  if (CV) {
    x <- cox1$linear.predictors[ord]
    if (sum(x^2)==0) stop("No cross-validation for null model implemented")
    X <- model.matrix(formula, data = data)
    X <- X[, -1, drop = FALSE]
    X <- X[ord, , drop = FALSE]
    if (progress) {
      m <- floor(log10(n)) + 1
      pre <- rep("\b", 2 * m + 1)
      cat("Calculating cross-validated survival curves:\n")
    }
    survmat <- matrix(NA, nt, n)
    for (i in 1:n) {
      if (progress) {
        cat(pre, i, "/", n, sep = "")
        flush.console()
      }
      cmini <- coxph(Surv(time[-i], status[-i]) ~ X[-i,
                                                    , drop = FALSE], method = "breslow")
      xi <- as.vector(X[-i, , drop = FALSE] %*% cmini$coef)
      ximini <- as.numeric(X[i, ] %*% cmini$coef)
      cmini <- coxph(Surv(time[-i], status[-i]) ~ xi, method = "breslow")
      survi <- survfit(cmini, newdata = data.frame(xi = ximini))
      survmat[, i] <- evalstep(survi$time, survi$surv,
                               tt, subst = 1)
    }
    if (progress) {
      cat("\nCalculating prediction error ...")
      flush.console()
    }
  }
  if (tt[1] > 0) {
    tsurv <- c(0, tt)
    survmat <- rbind(rep(1, n), survmat)
  }
  else tsurv <- tt
  nsurv <- length(tsurv)
  ## censoring
  if (missing(censdata)) censdata <- data
  coxcens <- coxph(censformula, censdata)
  ycens <- coxcens[["y"]]
  p <- ncol(ycens)
  tcens <- ycens[, p - 1]
  dcens <- ycens[, p]
  xcens <- coxcens$linear.predictors[ord]
  coxcens <- coxph(Surv(tcens, dcens) ~ xcens)
  # if no covariates, then use Kalbfleisch-Prentice type in survfit
  if (sum(xcens^2)==0)
    sfcens <- survfit(coxcens, newdata = data.frame(xcens = xcens),
                      type="kalbfl")
  else
    sfcens <- survfit(coxcens, newdata = data.frame(xcens = xcens))
  tcens <- sfcens$time
  censmat <- sfcens$surv
  if (tcens[1] > 0) {
    tcens <- c(0, tcens)
    censmat <- rbind(rep(1, n), censmat)
  }
  ncens <- length(tcens)
  ## vector at which prediction error is to be calculated
  if (missing(tout)) {
    tout <- unique(c(tsurv, tcens))
    tout <- c(tout, tout - width)
    tout <- tout[tout >= 0]
    if (min(tout) > 0)
      tout <- c(0, tout)
    tout <- sort(unique(tout))
    nout <- length(tout)
    if (nout > 0) {
      tout <- tout[-nout]
      nout <- nout - 1
    }
  }
  else {
    tout <- sort(tout)
    nout <- length(tout)
  }
  res <- pew(time, status, tsurv, survmat, tcens, censmat,
             width, FUN, tout)
  if (progress)
    cat("\n")
  return(res)
}



pe <- function(time,status,tsurv,survmat,tcens,censmat,FUN=c("KL","Brier"),tout)
{
  FUN <- match.arg(FUN)
  FUNn <- 2
  if (FUN=="Brier") FUNn <- 1
  n <- length(time) # no checks that this coincides with ncols of survmat and censmat
  nsurv <- length(tsurv) # no check that this coincides with nrows of survmat
  ncens <- length(tcens) # no check that this coincides with nrows of censmat
  nout <- length(tout)
  res <- .C('prederr',
            sn = as.integer(n),
            time = as.double(time),
            status = as.integer(status),
            snsurv = as.integer(nsurv),
            sncens = as.integer(ncens),
            snout = as.integer(nout),
            tsurv = as.double(tsurv),
            survmat = as.double(survmat),
            tcens = as.double(tcens),
            censmat = as.double(censmat),
            tout = as.double(tout),
            sFUN = as.integer(FUNn),
            err = double(nout),
            work = double(n)
  )
  res <- data.frame(time=tout,Err=res$err)
  attr(res,"score") <- FUN
  return(res)
}

pecox <- function(formula, censformula, data, censdata,
                  FUN = c("KL", "Brier"), tout, CV = FALSE, progress = FALSE)
{
  cox1 <- coxph(formula, data)
  y <- cox1[["y"]]
  p <- ncol(y)
  time <- y[, p - 1]
  status <- y[, p]
  n <- length(time)
  if (nrow(data) != n)
    stop("missing data in time or status not allowed")
  ord <- order(time, -status)
  time <- time[ord]
  status <- status[ord]
  tt <- sort(unique(time[status == 1]))
  nt <- length(tt)
  if (!CV) {
    progress <- FALSE
    x <- cox1$linear.predictors[ord]
    cox1 <- coxph(Surv(time, status) ~ x)
    # if no covariates, then use Kalbfleisch-Prentice type in survfit
    if (sum(x^2)==0)
      sf <- survfit(cox1, newdata = data.frame(x = x), type="kalbfl")
    else
      sf <- survfit(cox1, newdata = data.frame(x = x))
    tt <- sf$time
    survmat <- sf$surv
  }
  if (CV) {
    x <- cox1$linear.predictors[ord]
    if (sum(x^2)==0) stop("No cross-validation for null model implemented")
    X <- model.matrix(formula, data = data)
    X <- X[, -1, drop = FALSE]
    X <- X[ord, , drop = FALSE]
    if (progress) {
      m <- floor(log10(n)) + 1
      pre <- rep("\b", 2 * m + 1)
      cat("Calculating cross-validated survival curves:\n")
    }
    survmat <- matrix(NA, nt, n)
    for (i in 1:n) {
      if (progress) {
        cat(pre, i, "/", n, sep = "")
        flush.console()
      }
      cmini <- coxph(Surv(time[-i], status[-i]) ~ X[-i,
                                                    , drop = FALSE], method = "breslow")
      xi <- as.vector(X[-i, , drop = FALSE] %*% cmini$coef)
      ximini <- as.numeric(X[i, ] %*% cmini$coef)
      cmini <- coxph(Surv(time[-i], status[-i]) ~ xi, method = "breslow")
      survi <- survfit(cmini, newdata = data.frame(xi = ximini))
      survmat[, i] <- evalstep(survi$time, survi$surv,
                               tt, subst = 1)
    }
    if (progress) {
      cat("\nCalculating prediction error ...")
      flush.console()
    }
  }
  if (tt[1] > 0) {
    tsurv <- c(0, tt)
    survmat <- rbind(rep(1, n), survmat)
  }
  else tsurv <- tt
  nsurv <- length(tsurv)
  ## censoring
  if (missing(censdata)) censdata <- data
  coxcens <- coxph(censformula, censdata)
  ycens <- coxcens[["y"]]
  p <- ncol(ycens)
  tcens <- ycens[, p - 1]
  dcens <- ycens[, p]
  xcens <- coxcens$linear.predictors
  coxcens <- coxph(Surv(tcens, dcens) ~ xcens)
  # if no covariates, then use Kalbfleisch-Prentice type in survfit
  if (sum(xcens^2)==0)
    sfcens <- survfit(coxcens, newdata = data.frame(xcens = xcens),
                      type="kalbfl")
  else
    sfcens <- survfit(coxcens, newdata = data.frame(xcens = xcens))
  tcens <- sfcens$time
  censmat <- sfcens$surv
  if (tcens[1] > 0) {
    tcens <- c(0, tcens)
    censmat <- rbind(rep(1, n), censmat)
  }
  ncens <- length(tcens)
  if (missing(tout))
    tout <- unique(c(tsurv, tcens))
  tout <- sort(tout)
  nout <- length(tout)
  res <- pe(time, status, tsurv, survmat, tcens, censmat, FUN,
            tout)
  if (progress)
    cat("\n")
  return(res)
}




Fwindow <- function(object, width, variance=TRUE, conf.level=0.95)
{
  if (variance)
    sf <- data.frame(time=object$time,surv=object$surv,varHaz=object$std.err^2)
  else sf <- data.frame(time=object$time,surv=object$surv)
  sf$Haz <- -log(sf$surv)
  tt <- c(0,sf$time) # assuming no event at t=0 or t<0
  ttt <- c(0,tt,tt-width)
  ttt <- ttt[ttt >= 0]
  ttt <- sort(unique(ttt))
  ttt <- unique(ttt)
  H <- outer(c(0,sf$Haz),c(0,sf$Haz),"-")
  dimnames(H) <- list(tt,tt)
  tt <- c(tt,Inf)
  idx1 <- as.numeric(cut(ttt,tt,right=FALSE))
  idx2 <- as.numeric(cut(ttt+width,tt,right=FALSE))
  Fw <- diag(H[idx2,idx1])
  nt <- length(Fw)
  Fw[nt] <- Fw[nt-1]
  if (variance) {
    varH <- outer(c(0,sf$varHaz),c(0,sf$varHaz),"-")
    varFw <- diag(varH[idx2,idx1])
    varFw[nt] <- varFw[nt-1]
    ciwidth <- qnorm(1-(1-conf.level)/2)*sqrt(varFw)
    low <- Fw - ciwidth
    up <- Fw + ciwidth
    low[low<0] <- 0 # negative lower values of cum hazard set to zero
  }
  # Return on probability scale
  Fw <- 1 - exp(-Fw)
  if (variance) {
    low <- 1 - exp(-low)
    up <- 1 - exp(-up)
    res <- data.frame(time=ttt,Fw=Fw,low=low,up=up)
  }
  else res <- data.frame(time=ttt,Fw=Fw)
  attr(res,"width") <- width
  return(res)
}


evalstep <- function(time, stepf, newtime, subst=-Inf, to.data.frame=FALSE)
{
  n <- length(time)
  if (is.vector(stepf))
    if (length(stepf) != n)
      stop("arguments 'time' and 'stepf' should have the same length")
  if (is.matrix(stepf) | is.data.frame(stepf))
    if (nrow(stepf) != n)
      stop("argument 'stepf' should have the same number of rows as length of argument 'time'")
  # time should be ordered, not contain duplicates, and not contain +/- Inf
  if (any(!(order(time) == 1:n))) stop("argument 'time' should be ordered")
  if (any(duplicated(time))) stop("argument 'time' should not contain duplicates")
  if (any(is.infinite(time))) stop("(-) infinity not allowed in 'time'")
  idx <- cut(newtime,c(-Inf,time,Inf),right=FALSE)
  idx <- as.numeric(idx)
  if (is.vector(stepf)) res <- c(subst,stepf)[idx]
  if (is.matrix(stepf) | is.data.frame(stepf)) {
    stepf <- rbind(rep(subst,ncol(stepf)),stepf)
    res <- stepf[idx,]
  }
  if (to.data.frame) return(data.frame(newtime=newtime,res=res))
  else return(res)
}


deb <- function(x, method=c("print","cat")) {
  call <- match.call()
  method <- match.arg(method)
  if (method=="print") {
    cat(deparse(call$x), ":\n")
    print(x)
  }
  if (method=="cat") cat(deparse(call$x), ":", x, "\n")
  flush.console()
}



CVPL <- function(formula, data, progress=TRUE, overall=FALSE, shrink=1)
{
  # Extract data (time and status)
  cox1 <- coxph(formula, data)
  y <- cox1[["y"]]
  p <- ncol(y)
  time <- y[, p - 1]
  status <- y[, p]
  n <- length(time)
  if (nrow(data) != n)
    stop("missing data in time or status not allowed")
  ord <- order(time, -status)
  time <- time[ord]
  status <- status[ord]
  # Model matrix
  X <- model.matrix(formula, data = data)
  if (ncol(X)>1) {
    X <- X[, -1, drop = FALSE]
    X <- X[ord, , drop = FALSE]
    X <- t(t(X) - apply(X,2,mean)) # Centering
  }
  if (overall) {
    cfull <- coxph(Surv(time,status) ~ X, data = data, method="breslow")
    pl1 <- cfull$loglik[2]
  }
  if (progress) {
    m <- floor(log10(n)) + 1
    pre <- rep("\b", 2 * m + 1)
  }
  rat <- rep(NA,n)
  res <- 0
  for (i in 1:n) {
    if (progress) {cat(pre,i,"/",n,sep=""); flush.console()}
    # leave out i
    if (!overall) { # estimate coefficients without i
      cmin1 <- coxph(Surv(time[-i], status[-i]) ~ X[-i,], method="breslow")
      rat[-i] <- exp(shrink * cmin1$linear.predictors)
      rat[i] <- exp(sum(X[i,] * shrink * cmin1$coef))
      rcsrat <- rev(cumsum(rev(rat)))
      rat1 <- rat[status==1]
      rcsrat1 <- rcsrat[status==1]
      pl1 <- sum(log(rat1) - log(rcsrat1))
      if (shrink==1) pl2 <- cmin1$loglik[2]
    }
    else # use overall estimate to define rat
      rat <- exp(cfull$linear.predictors)
    if (overall | shrink != 1) {
      # Can't use loglik of cmin1 now
      ratmin1 <- rat[-i]
      rcsratmin1 <- rev(cumsum(rev(ratmin1)))
      rat1 <- ratmin1[status[-i]==1]
      rcsrat1 <- rcsratmin1[status[-i]==1]
      pl2 <- sum(log(rat1) - log(rcsrat1))
    }
    res <- res + (pl1-pl2)
  }
  if (progress) cat(paste(rep("\b", 4*m+3), collapse=""))
  return(res)
}




CVcindex <- function(formula, data, type="single", matrix=FALSE)
{
  cox1 <- coxph(formula, data)
  y <- cox1[["y"]]
  p <- ncol(y)
  time <- y[,p-1]
  status <- y[,p]
  # Center covariates
  # First get design matrix
  X <- model.matrix(formula, data = data)
  X <- X[,-1] # remove intercept
  X <- t(t(X) - apply(X,2,mean))
  cfull <- coxph(Surv(time,status) ~ X, data = data, method="breslow")
  n <- nrow(data)
  if (type=="single" | type=="pair") {
    m <- floor(log10(n))+1 # for formatting progress count
    pre <- rep("\b",2*m+1)
    xmat <- matrix(NA,n,n) # in columns (-i)
    # xmat[j,i] will contain PI_{j,(-i)}
    for (i in 1:n) {
      cat(pre,i,"/",n,sep=""); flush.console()
      # leave out i
      cmin1 <- coxph(Surv(time[-i], status[-i]) ~ X[-i,], method="breslow")
      # evaluate at all j except i
      xmat[-i,i] <- cmin1$linear.predictors
      # evaluate at i
      xmat[i,i] <- sum(X[i,] * cmin1$coef)
    }
    cat("\n")
  }
  if (type=="single") {
    formula <- as.formula("Surv(time,status) ~ x")
    ndata <- data.frame(time=time,status=status,x=diag(xmat))
    res <- cindex(formula=formula, data=ndata)
    if (matrix) res <- list(concordant=res$concordant,total=res$total,cindex=res$cindex,matrix=xmat)
  }
  if (type=="pair") {
    n <- length(time) # check if = length(status) and length(x)
    ord <- order(time,-status)
    time <- time[ord]
    status <- status[ord]
    xmat <- xmat[ord,ord]
    # pairs (i,j) for which the smallest observed time is an event time
    wh <- which(status==1)
    total <- concordant <- 0
    for (i in wh) {
      if (i < n) {
        for (j in ((i+1):n)) {
          if (time[j] > time[i]) { # ties not counted
            total <- total + 2
            if (xmat[j,i] < xmat[i,i]) concordant <- concordant + 1
            if (xmat[j,j] < xmat[i,j]) concordant <- concordant + 1
          }
        }
      }
    }
    if (matrix) res <- list(concordant=concordant,total=total,cindex=concordant/total,matrix=xmat) else res <- list(concordant=concordant,total=total,cindex=concordant/total)
  }
  if (type=="fullpairs") {
    m <- floor(log10(n*(n-1)/2))+1 # for formatting progress count
    pre <- rep("\b",2*m+1)
    # xmat[i,j] will contain PI_{i,(-i,-j)}; xmat[j,i] will contain PI_{j,(-i,-j)}
    xmat <- matrix(NA,n,n)
    cnt <- 0
    for (i in 1:n) {
      if (i < n) {
        for (j in ((i+1):n)) {
          cnt <- cnt+1
          cat(pre,cnt,"/",n*(n-1)/2,sep=""); flush.console()
          # leave out i and j
          cmin2 <- coxph(Surv(time[-c(i,j)], status[-c(i,j)]) ~ X[-c(i,j),], method="breslow")
          # evaluate at i
          xmat[i,j] <- sum(X[i,] * cmin2$coef)
          # evaluate at j
          xmat[j,i] <- sum(X[j,] * cmin2$coef)
        }
      }
    }
    ord <- order(time,-status)
    time <- time[ord]
    status <- status[ord]
    xmat <- xmat[ord,ord]
    # pairs (i,j) for which the smallest observed time is an event time
    wh <- which(status==1)
    total <- concordant <- 0
    for (i in wh) {
      if (i < n) {
        for (j in ((i+1):n)) {
          if (time[j] > time[i]) {# ties not counted
            total <- total + 1
            if (xmat[j,i] < xmat[i,j]) concordant <- concordant + 1
          }
        }
      }
    }
    if (matrix) res <- list(concordant=concordant,total=total,cindex=concordant/total,matrix=xmat) else res <- list(concordant=concordant,total=total,cindex=concordant/total)
  }
  cat("\n")
  attr(res, "type") <- type
  return(res)
}



cutLM <- function(data, outcome, LM, horizon, covs,
                  format = c("wide","long"), id, rtime, right=TRUE)
{
  format <- match.arg(format)
  if (format=="wide") {
    LMdata <- data
    if (!is.null(covs$varying))
      LMdata[[covs$varying]] <- 1 - as.numeric(LMdata[[covs$varying]] > LM)
  } else {
    if (missing(id))
      stop("argument 'id' should be specified for long format data")
    if (missing(rtime))
      stop("argument 'rtime' should be specified for long format data")
    ord <- order(data[[id]],data[[rtime]])
    data <- data[ord,]
    ids <- unique(data[[id]])
    n <- length(ids)
    # initialize LMdata; copy first row of each subject
    LMdata <- data[which(!duplicated(data[[id]])),]
    for (i in 1:n) {
      wh <- which(data[[id]]==ids[i])
      di <- data[wh,]
      idx <- cut(LM,c(data[[rtime]][wh],Inf),right=right,labels=FALSE)
      if (!is.na(idx)) LMdata[i,] <- di[idx,]
      else {
        LMdata[i,] <- di[1,]
        LMdata[[covs$varying]][i] <- NA
        LMdata[[rtime]][i] <- NA
      } 
    }
  }
  LMdata <- LMdata[LMdata[[outcome$time]] > LM,]
  if (format=="long") LMdata <- LMdata[!is.na(LMdata[[id]]),]
  # apply administrative censoring at horizon
  LMdata[outcome$status] <- LMdata[[outcome$status]] *
    as.numeric(LMdata[[outcome$time]] <= horizon)
  LMdata[outcome$time] <- pmin(as.vector(LMdata[[outcome$time]]),horizon)
  LMdata$LM <- LM
  if (format=="long")
    cols <- match(c(id,outcome$time,outcome$status,covs$fixed,covs$varying,rtime,"LM"),
                  names(LMdata))
  else
    cols <- match(c(outcome$time,outcome$status,covs$fixed,covs$varying,"LM"),
                  names(LMdata))
  return(LMdata[,cols])
}


cindex <- function(formula, data)
{
  cox1 <- coxph(formula, data)
  y <- cox1[["y"]]
  p <- ncol(y)
  time <- y[,p-1]
  status <- y[,p]
  x <- cox1$linear.predictors
  n <- length(time)
  ord <- order(time,-status)
  time <- time[ord]
  status <- status[ord]
  x <- x[ord]
  # pairs (i,j) for which the smallest observed time is an event time
  wh <- which(status==1)
  total <- concordant <- 0
  for (i in wh) {
    for (j in ((i+1):n)) {
      if (time[j] > time[i]) {# ties not counted
        total <- total + 1
        if (x[j] < x[i]) concordant <- concordant + 1
        if (x[j] == x[i]) concordant <- concordant + 0.5
      }
    }
  }
  return(list(concordant=concordant,total=total,cindex=concordant/total))
}



AUCw <- function(formula,data,width)
{
  cox1 <- coxph(formula, data)
  y <- cox1[["y"]]
  p <- ncol(y)
  time <- y[,p-1]
  status <- y[,p]
  x <- cox1$linear.predictors
  tt <- sort(unique(time[status==1]))
  ttw <- c(tt,tt-width)
  ttw <- ttw[ttw>0]
  ttw <- sort(unique(ttw))
  ntw <- length(ttw)
  AUCw <- rep(NA,ntw)
  for (j in 1:ntw) {
    twj <- ttw[j]
    ttj <- tt[((tt >= twj) & (tt <= twj+width))]
    ntj <- length(ttj)
    AUCt <- rep(NA,ntj)
    numsum <- denomsum <- 0
    for (i in 1:ntj) {
      ti <- ttj[i]
      # risk set
      Y <- sum(time>=ti)
      R <- which(time>ti) # !!! R(ti) is which(time>=ti), but for the "controls", j=i should be excluded, only valid without ties
      xi <- x[time==ti]
      num <- sum(x[R]<xi) + 0.5*sum(x[R]==xi)
      numsum <- numsum + num
      denomsum <- denomsum + Y-1 # Also only valid in absence of ties
    }
    AUCw[j] <- numsum/denomsum
  }
  res <- data.frame(time=ttw,AUCw=AUCw)
  attr(res, "width") <- width
  return(res)
}


AUC <- function(formula,data,plot=TRUE)
{
  cox1 <- coxph(formula, data)
  y <- cox1[["y"]]
  p <- ncol(y)
  time <- y[,p-1]
  status <- y[,p]
  x <- cox1$linear.predictors
  tt <- sort(unique(time[status==1]))
  nt <- length(tt)
  AUCt <- rep(NA,nt)
  numsum <- denomsum <- 0
  for (i in 1:nt) {
    ti <- tt[i]
    # risk set
    Y <- sum(time>=ti)
    R <- which(time>ti) # !!! R(ti) is which(time>=ti), but for the "controls", j=i should be excluded, only valid without ties
    xi <- x[time==ti]
    num <- sum(x[R]<xi) + 0.5*sum(x[R]==xi)
    AUCt[i] <- num/(Y-1) # Also only valid in absence of ties
    numsum <- numsum + num
    denomsum <- denomsum + Y-1 # Also only valid in absence of ties
  }
  AUC <- numsum/denomsum
  if (plot) {
    plot(tt,AUCt,xlab="Time t",ylab="AUC(t)")
    lines(lowess(data.frame(tt,AUCt)))
    abline(h=0.5,lty=3)
  }
  return(list(AUCt=data.frame(time=tt,AUC=AUCt),AUC=AUC))
}




toleranceplot <- function(formula,data,coverage=0.8,horizon,plot=TRUE,xlab)
{
  cox1 <- coxph(formula, data)
  y <- cox1[["y"]]
  p <- ncol(y)
  time <- y[,p-1]
  status <- y[,p]
  x <- cox1$linear.predictors
  xs <- sort(unique(x))
  nx <- length(xs)
  if (missing(horizon)) horizon <- max(time[status==1])*1.05
  ylim <- c(0,horizon)
  if (missing(xlab)) xlab <- "x"
  if (plot) plot(range(xs),c(0,horizon),type="n",xlab=xlab,ylab="Tolerance interval")
  cx <- coxph(Surv(time,status) ~ x, method="breslow")
  res <- matrix(NA,nx,3)
  for (i in 1:nx) {
    xi <- xs[i]
    nd <- data.frame(x=xi)
    sf <- survfit(cx,newdata=nd)
    sf <- data.frame(time=sf$time,surv=sf$surv)
    low <- max(sf$time[sf$surv>1-(1-coverage)/2])
    up <- min(sf$time[sf$surv<(1-coverage)/2])
    if (is.infinite(up)) up <- horizon
    lines(rep(xi,2),c(low,up),type="l")
    res[i,] <- c(xi,low,up)
  }
  res <- as.data.frame(res)
  names(res) <- c("x","lower","upper")
  attr(res, "coverage") <- coverage
  attr(res, "horizon") <- horizon
  return(res)
}
