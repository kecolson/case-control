ccwc.weights <-function(entry=0, # time of entry to follow up
                        exit, # time of exit from follow up
                        fail, # status on exit (1 = fail, 0 = censored)
                        origin=0, # origin of analysis time scale
                        controls=1, # number of controls to select for each case
                        weights = NULL, # sampling weights
                        match=list(), # list of categorical variables on which to match cases & controls
                        include=list(), # list of other variables to be carried across into case control study
                        data=NULL, # data frame in which to look for input variables
                        silent=FALSE # If FALSE, echos a . to the screen for each case-control set created; otherwise produces no output.
)
{
  # Check arguments
  entry <- eval(substitute(entry), data)
  exit <- eval(substitute(exit), data)
  fail <- eval(substitute(fail), data)
  origin <- eval(substitute(origin), data)
  weights <- eval(substitute(weights), data) #Allows weights to be evaluated from the data frame specified
  
  n <- length(fail)
  if (length(exit)!=n)
    stop("All vectors must have same length")
  if (length(entry)!=1 && length(entry)!=n)
    stop("All vectors must have same length")
  if (length(weights)!=n & !is.null(weights)) #Do we want to add this check for weights?
    stop("All vectors must have same length")
  if (length(origin)==1) {
    origin <- rep(origin, n)
  }
  else {
    if (length(origin)!=n)
      stop("All vectors must have same length")
  }
  
  # Transform times to correct scale
  t.entry <- as.numeric(entry - origin)
  t.exit <- as.numeric(exit - origin)
  # match= argument
  marg <- substitute(match)
  if (mode(marg)=="name") {
    match <- list(eval(marg, data))
    names(match) <- as.character(marg)
  }
  else if (mode(marg)=="call" && marg[[1]]=="list") {
    mnames <- names(marg)
    nm <- length(marg)
    if (!is.null(mnames)) {
      if (nm>1) {
        for (i in 2:nm) {
          if (mode(marg[[i]])=="name")
            mnames[i] <- as.character(marg[[i]])
          else
            stop("illegal argument (match)")
        }
      }
      else {
        for (i in 2:nm) {
          if (mode(marg[[i]])=="name")
            mnames[i] <- as.character(marg[[i]])
          else
            stop("illegal argument (match)")
        }
        mnames[1] <= ""
      }
    }
    names(marg) <- mnames
    match <- eval(marg, data)
  }
  else {
    stop("illegal argument (match)")
  }
  
  m <- length(match)
  mnames <- names(match)
  if (m>0) {
    for (i in 1:m) {
      if (length(match[[i]])!=n) {
        stop("incorrect length for matching variable")
      }
    }
  }
  # include= argument
  iarg <- substitute(include)
  if (mode(iarg)=="name") {
    include <- list(eval(iarg, data))
    names(include) <- as.character(iarg)
  }
  else if (mode(iarg)=="call" && iarg[[1]]=="list") {
    ni <- length(iarg)
    inames <- names(iarg)
    if (ni>1) {
      if (!is.null(inames)) {
        for (i in 2:ni) {
          if (mode(iarg[[i]])=="name")
            inames[i] <- as.character(iarg[[i]])
          else
            stop("illegal argument (include)")
        }
      }
      else {
        for (i in 2:ni) {
          if (mode(iarg[[i]])=="name")
            inames[i] <- as.character(iarg[[i]])
          else
            stop("illegal argument (include)")
        }
        inames[1] <= ""
      }
    }
    names(iarg) <- inames
    include <- eval(iarg, data)
  }
  else {
    stop("illegal argument (include)")
  }
  
  ni <- length(include)
  inames <- names(include)
  if (ni>0) {
    for (i in 1:ni) {
      if (length(include[[i]])!=n) {
        stop("incorrect length for included variable")
      }
    }
  }
  # create group codes using matching variables
  grp <- rep(1,n)
  pd <- 1
  if (m>0) {
    for (im in 1:m) {
      v <- match[[im]]
      if (length(v)!=n)
        stop("All vectors must have same length")
      if (!is.factor(v))
        v <- factor(v)
      grp <- grp + pd*(as.numeric(v) - 1)
      pd <- pd*length(levels(v))
    }
  }
  # Create vectors long enough to hold results
  nn <- (1+controls)*sum(fail!=0) # Case + # of controls * number of cases
  pr <- numeric(nn)
  sr <- numeric(nn)
  tr <- vector("numeric", nn)
  fr <- numeric(nn)
  nn <- 0
  # Sample each group
  if (!silent) {
    cat("\nSampling risk sets: ")
  }
  set <- 0
  nomatch <- 0
  incomplete <- 0
  ties <- FALSE
  fg <- unique(grp[fail!=0])
  failures_map <- (fail!=0) # changed from original function for minor gains in efficiency
  
  for (g in fg) {
    group_map <- (grp==g) # changed from original function for minor gains in efficiency
    # Failure times
    ft <- unique( t.exit[group_map & failures_map] )
    # Create case-control sets
    for (tf in ft) {
      if (!silent) {
        cat(".")
      }
      tf_eq <- (t.exit == tf) # changed from original function for minor gains in efficiency
      tf_gt <- (t.exit > tf) # changed from original function for minor gains in efficiency
      te_lte <- (t.entry <= tf) # changed from original function for minor gains in efficiency
      
      set <- set+1
      case <- group_map & tf_eq & failures_map
      
      ncase <- sum(case)
      if (ncase>0)
        ties <- TRUE
      noncase <- group_map & (te_lte) &
        (tf_eq | tf_gt) & !case
      noncaseweights <- weights[group_map & (te_lte) & # For a given case, pull the noncase weights in the given group
                                  (tf_eq | tf_gt) & !case]
      
      ncont <- controls*ncase
      if (ncont>sum(noncase)) {
        ncont <- sum(noncase)
        if (ncont>0) incomplete <- incomplete + 1
      }
      if (ncont>0) {
        newnn <- nn+ncase+ncont
        sr[(nn+1):newnn] <- set
        tr[(nn+1):newnn] <- tf
        fr[(nn+1):(nn+ncase)] <- 1
        fr[(nn+ncase+1):newnn] <- 0
        pr[(nn+1):(nn+ncase)] <- (1:n)[case]
        pr[(nn+ncase+1):(newnn)] <-
          sample((1:n)[noncase], size=ncont, replace = TRUE, #Should replace be F or T??
                 prob = noncaseweights) #added sampling weights here
        nn <- newnn
      }
      else {
        nomatch <- nomatch + ncase
      }
    }
  }
  if (!silent) {
    cat("\n")
  }
  res <- vector("list", 4+m+ni)
  if (nn>0) {
    res[[1]] <- sr[1:nn]
    res[[2]] <- map <- pr[1:nn]
    res[[3]] <- tr[1:nn] + origin[map]
    res[[4]] <- fr[1:nn]
  }
  if (m>0) {
    for (i in 1:m) {
      res[[4+i]] <- match[[i]][map]
    }
  }
  if (ni>0) {
    for (i in 1:ni) {
      res[[4+m+i]] <- include[[i]][map]
    }
  }
  names(res) <- c("Set", "Map", "Time", "Fail", mnames, inames)
  if (incomplete>0)
    warning(paste(incomplete, "case-control sets are incomplete"))
  if (nomatch>0)
    warning(paste(nomatch, "cases could not be matched"))
  if (ties)
    warning("there were tied failure times")
  data.frame(res)
}