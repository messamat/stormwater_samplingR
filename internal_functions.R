#Modify rcompanion transformTukey to return lambda
#This Tukey ladder computation is designed to take in a raster attribute table and repeatedly sample it to determine whether a transformation 
#is required, and if so, performs the transformation
getlambda <- function(x, start = -2, end = 2, int = 0.025, verbose = FALSE, statistic = 1) {
  n = (end - start)/int
  lambda = as.numeric(rep(0, n))
  W = as.numeric(rep(0, n))
  Shapiro.p.value = as.numeric(rep(0, n))
  if (statistic == 2) {
    A = as.numeric(rep(1000, n))
    Anderson.p.value = as.numeric(rep(0, n))
  }
  for (i in (1:n)) {
    lambda[i] = signif(start + (i - 1) * int, digits = 4)
    if (lambda[i] > 0) {
      TRANS = x^lambda[i]
    }
    if (lambda[i] == 0) {
      TRANS = log(x)
    }
    if (lambda[i] < 0) {
      TRANS = -1 * x^lambda[i]
    }
    W[i] = NA
    if (statistic == 2) {
      A[i] = NA
    }
    if (any(is.infinite(TRANS)) == FALSE & any(is.nan(TRANS)) == 
        FALSE) {
      W[i] = signif(shapiro.test(TRANS)$statistic, digits = 4)
      Shapiro.p.value[i] = signif(shapiro.test(TRANS)$p.value, 
                                  digits = 4)
      if (statistic == 2) {
        A[i] = signif(ad.test(TRANS)$statistic, digits = 4)
        Anderson.p.value[i] = signif(ad.test(TRANS)$p.value, 
                                     digits = 4)
      }
    }
  }
  if (statistic == 2) {
    df = data.frame(lambda, W, Shapiro.p.value, A, Anderson.p.value)
  }
  if (statistic == 1) {
    df = data.frame(lambda, W, Shapiro.p.value)
  }
  if (verbose == TRUE) {
    print(df)
  }
  if (statistic == 1) {
    df2 = df[with(df, order(-W)), ]
  }
  if (statistic == 2) {
    df2 = df[with(df, order(A)), ]
  }
  cat("\n")
  print(df2[1, ])
  cat("\n")
  cat("if (lambda >  0){TRANS = x ^ lambda}", "\n")
  cat("if (lambda == 0){TRANS = log(x)}", "\n")
  cat("if (lambda <  0){TRANS = -1 * x ^ lambda}", "\n")
  cat("\n")
  return(df2[1, "lambda"])
}

transformTukey_lambda <- function (x,start = -2, end = 2, int = 0.025, rastertab=T, rep=10, verbose = FALSE, statistic = 1) {
  if (rastertab) { #If analyze a raster attribute table
    lambda_iter <- NULL
    for (iter in 1:rep) {
      x$prob <- x$Count/sum(x$Count)
      xsamp <- sample(x$Value, size=5000, replace=T, prob=x$prob)
      lambda_iter[iter] <- getlambda(x=xsamp)
      lambda <- mean(lambda_iter)
    }
    if (lambda > 0) {
      TRANS = x$Value^lambda
    }
    if (lambda == 0) {
      TRANS = log(x$Value)
    }
    if (lambda < 0) {
      TRANS = -1 * x$Value^lambda
    }
  } else{ #If processing a vector
    if (length(x)<=5000) { 
      lambda <- getlambda(x=x)
    } else {
      for (iter in 1:rep) {
        xsamp <- sample(x, size=5000)
        lambda_iter[iter] <- getlambda(x=xsamp)
        lambda <- mean(lambda_iter)
      }
    }
    if (lambda > 0) {
      TRANS = x^lambda
    }
    if (lambda == 0) {
      TRANS = log(xValue)
    }
    if (lambda < 0) {
      TRANS = -1 * x^lambda
    }
  }
  return(list(TRANS,lambda))
}

#Back transformTukey
backtransformTukey <- function (x, lambda) 
{
  if (lambda > 0) {
    BACKTRANS = x^(1/lambda)
  }
  if (lambda == 0) {
    BACKTRANS = exp(x)
  }
  if (lambda < 0) {
    BACKTRANS = (-x)^(1/lambda)
  }
  return(BACKTRANS)
}

#Function to bin raster attribute table
bin_rastertab<- function(tab, nbins, tukey=T, rep=10, rastertab) {
  df <- tab
  if (tukey) {
    dftrans <- transformTukey_lambda(df, rep=rep, rastertab=rastertab)
    df$valtuk <- dftrans[[1]]
    lambda <- dftrans[[2]]
  } else {
    df$valtuk <- df$Value
    lambda <- 1
  }
  binsize <- (max(df$valtuk)-min(df$valtuk))/nbins
  #minbin = min(df$valtuk)-min(df$valtuk)%%binsize+binsize
  #maxbin = max(df$valtuk)-max(df$valtuk)%%binsize+binsize
  bin <- with(df,seq(min(df$valtuk),max(df$valtuk), binsize))
  binback <- data.frame(bin,binmax=backtransformTukey(bin, lambda))
  setDT(binback)[,binmin:=shift(binmax, 1, type='lag'),]
  binback$id <- as.numeric(row.names(binback))-1
  df$bin <- .bincode(x=df$valtuk, breaks=bin, TRUE, include.lowest = T)
  if (rastertab) {
    df_hist <- setDT(df)[,sum(Count), .(bin)]
  } else {
    df_hist <- setDT(df)[,length(Value), .(bin)]
  }
  colnames(df_hist) <- c('id','count')
  df_hist <- merge(df_hist, binback, by='id', all.x=T, all.y=T)
  return(df_hist[-1,])
}
