library(parallel)
#library(fastmatch)
#library(ggplot2)
library(Rfast)
library(Rfast2)
library(ape)

myCores <- if (detectCores() > 2L) detectCores() - 2


################################
###### FUNCTION concatenate ####
################################
# This is a modified version of cibnd.DNAbin 
catAAbin <- function(x,
           mss = '-',
           check.names = TRUE,
           fill.with.gaps = TRUE,
           quiet = FALSE)
  {
    obj <- x
    n <- length(obj)
    if (n == 1)
      return(obj[[1]])
    for (i in 1:n)
      if (!is.matrix(obj[[1]]))
        stop("the 'cbind' method for \"DNAbin\" accepts only matrices")
    NR <- unlist(lapply(obj, nrow))
    for (i in 1:n)
      class(obj[[i]]) <- NULL
    if (check.names) {
      NMS <- lapply(obj, rownames)
      for (i in 1:n)
        if (anyDuplicated(NMS[[i]]))
          stop("Duplicated rownames in matrix ",
               i,
               ": see ?cbind.DNAbin")
      nms <- unlist(NMS)
      if (fill.with.gaps) {
        NC <- unlist(lapply(obj, ncol))
        nms <- sort(unique(nms))
        qm <- charToRaw(mss)
        ans <- matrix(qm, length(nms), sum(NC))
        rownames(ans) <- nms
        from <- 1
        for (i in 1:n) {
          to <- from + NC[i] - 1
          k <- match(NMS[[i]], nms)
          ans[k, from:to] <- obj[[i]]
          from <- to + 1
        }
      }
      else {
        tab <- table(nms)
        ubi <- tab == n
        nms <- names(tab)[which(ubi)]
        ans <- obj[[1]][nms, , drop = FALSE]
        for (i in 2:n)
          ans <- cbind(ans, obj[[i]][nms, ,
                                     drop = FALSE])
        if (!quiet && !all(ubi))
          warning("some rows were dropped.")
      }
    }
    else {
      if (length(unique(NR)) > 1)
        stop("matrices do not have the same number of rows.")
      ans <- matrix(unlist(obj), NR)
      rownames(ans) <- rownames(obj[[1]])
    }
    class(ans) <- "AAbin"
    ans
  }


##############################
##### find factor ############
##############################
# Find the closest number to x that is divisible by y
# e.g. diiv(50, 10) #ans 50
# e.g. diiv(54, 10) #ans 50
# e.g. diiv(55, 10) #ans 60
# e.g. diiv(46, 10) #ans 50

diiv <-
  function(x, y) {
    if (x %% y >= y - (x %% y))
      x <- x + (y - (x %% y))
    else
      x <- x - (x %% y)
    x
  }

##############################
######### unwhich ############
##############################
unwhich <- 
  function (whch, dim = max(whch)) {
  y <- array(logical(length(whch)), dim = dim)
  y[whch] <- TRUE
  y
}

##########################################
######### subsample Function #############

subsamG <-
  function(x, y = 0L, z = 0L, v = 0L, seed = NULL) {
    # x: object of class multiAAbin list, that contains the single gene alignments
    # y: The desired subsample size (number of genes per subsample)
    # z: The desired approximate Number of subsamples to generate
    # v: The desired number times each gene to be sampled
    # --Note-- set either z or v, otherwise, v is considered.
    
    # adding seed for reproducibility
    # use provided seed, then revert to the original random seed on exit
    if (!is.null(seed)) {
      old_seed <- .Random.seed
      set.seed(seed)
      on.exit({.Random.see <<- old_seed})
    }
    
    orgx <- x
    u <- sapply(x, ncol)
    w <- length(x)
    
    nt <- cumsum(u)
    nf <- c(1, nt[-(length(nt))] + 1)
    x <- catAAbin(x)
    aa <- charToRaw(c('GPAVLIMCFYWHKRQNEDSTgpavlimcfywhkrqnedst'))
    xgap <- lapply(seq_along(nf), function(z) x[,(nf[z]):(nt[z])])
    # Sequences with less than 20 non-gap positions are considered as missing
    xgap <- sapply(xgap, function(z1) apply(z1, 1, function(z2) sum(z2 %in% aa) >= 20)) 
    
    if (y == 0) y <- w
    
    if (z == 0 & v == 0) {
      z <- 1L
      } else if (v == 0) {
      z <- diiv(z, w)
      } else z <- diiv(w*v/y, w)  #same as# z <- (diiv(v,y) / y) * w
        
    stopifnot(y <= w)

    cnt <- 1
    
    while (cnt <= 10) {
      mb <- array(FALSE, c(z, w))
      r <- sample(w, w)
      
      for (i in seq_len(z)) {
        if (length(r) >= y) {
          mb[i, r[1:y]] <- TRUE
          r <- r[-(1:y)]
        } else {
          r1 <- setdiff(seq_len(w), r)
          r1 <- sample(r1, length(r1))
          
          n <- y - length(r)
          mb[i, c(r, r1[1:n])] <- TRUE
          r <- c(r, r1[-(1:n)])
          r <- sample(r, length(r))
         
        }
      }
      
      inspect <- sum(xgap %*% t(mb) <= 0)
      if (inspect <= 0) {
        message("Sucess after ", cnt,  " attempt(s)")
          return(setNames(list(
          orgx,
          xgap,
          mb), c("alignments", "gapmap", "P1subsets")))
        break
      } else {
        message("Attempt ", cnt, " failed! trying again...");cnt <- cnt + 1}
    }
    
    message("\nYou may want to consider removing gappy taxa, or increase the subset size then try again")
   
  }

########### END ###############################


#########################################################
####### To print the subsets fasta files for P1 #########
writep1subsets <- function(x, sufx, path = NULL) {
  # x: Successful output from the subsamG function
  # sufx: suffix for the subsamples fasta files
  # path: destination path to save the output fasta files
  stopifnot(!is.null(path))
  
  u <- sapply(x[[1]], ncol)
  k <- unlist(mapply(rep_len, seq_along(u), u))
  xd <- catAAbin(x[[1]])
  xr <- split(seq_len(ncol(xd)), k)
      if (dir.exists(path))  {
        message("Directory <", path, "> already exist, choose another name!\n")
        } else {
          dir.create(path = path, recursive = TRUE)
          for (i in seq_len(nrow(x[[3]]))) {
            write.FASTA(xd[, unlist(xr[x[[3]][i,]])], file = sprintf("%s/%s_%05d.fast", path, sufx, i))
          }
        }
}



#############
### NOTE! ###
#############

# Order of the trees should be similar to the fasta files
# using the numbering sequence suffix in the fasta file names.

############################################################    
######### Evaluate the trees from subsamples ###############    

txerr <- function(tr, subg, rm = 0.10, rgn = 0.5, rtx = 0.5, outrgx = "^ub") {
  #tr: Object of class multiphylo, that have the inferred trees from the subsamples.
  #subg: Successful output from the subsamG, from which the trees assigned to tr were inferred
  #rm: percentage to mask
  #rgn: threshold of missing data after masking gene-wise
  #rtx: threshold of missing data after masking taxon-wise
  #outrgx: Regular expression unique to the outgroup
  
  bn <- subg[[1]]
  tb <- subg[[3]]
  tbb <- sweep(tb,2,colSums(tb),`/`) # Normalize tb
  
  tr <- sapply(tr, function(z) {
    cophenetic(compute.brlen(z,1))
    }, simplify = "array" )
  
  igs <- vapply(1:dim(tr)[2], function(z) {
    rowTrimMean(tr[,z,], a = 0.40)
    }, array(0.0, dim(tr)[2]))
  
  ig <- abs(sweep(tr,1:2,igs))^2
  
  og <- ig
  outg <- grep(outrgx,colnames(tr))
  ig[outg,,] <- 0L
  ig[,outg,] <- 0L
  
  # Rescaling twice
  ig <- sweep(ig^2,3,(colMeans(ig,dims = 2)),`/`)
  ig <- sweep(ig^2,3,colMeans(ig,dims = 2),`/`)
  
  ig <- colSums(ig)
  igu <- ig %*% tbb
  su <- which.min(apply(igu,1, function(n) {
    if (sum(n) > 0) sd(n) * mean(n) else Inf }))
  
  kp <-  0.1
  stu <- c(outg,su)
  gfrq <- colSums(tbb > 0)
  frm <- nrow(tbb) - gfrq + 1
  tom <- frm + ceiling(gfrq * kp) + 2
  
  igm <- t(vapply(seq_len(nrow(ig)), function(p) {
    mmt <- ig[p,] * tbb;
    mms <- colSort(mmt);
    vapply(seq_len(ncol(mms)), function(z) {
      mean(mms[frm[z]:tom[z],z])
      }, 0.0)
    }, array(0.0, ncol(tbb))))
    
  igm[is.nan(igm)] <- 0
  igm <- igm * igu
  
  og[-(stu),,] <- 0L
  og[,-(stu),] <- 0L
  
  # Rescaling twice as with the ingroup
  og <- sweep(og^2,3,colMeans(og,dims = 2),`/`)
  og <- sweep(og^2,3,colMeans(og,dims = 2),`/`)
  
  og <- colSums(og)
  ogu <- og %*% tbb

  ogm <- t(vapply(seq_len(nrow(og)), function(p) {
    mmt <- og[p,] * tbb;
    mms <- colSort(mmt);
    vapply(seq_len(ncol(mms)), function(z) {
      mean(mms[frm[z]:tom[z],z])
      }, 0.0)
    }, array(0.0, ncol(tbb))))
  
    ogm[is.nan(ogm)] <- 0
    ogm <- ogm * ogu
  ogm[su,] <- 0
  
  gaps <- charToRaw(c('-?xX'))
  nt <- cumsum(sapply(bn,ncol))
  nf <- c(1, nt[-(length(nt))] + 1)
  bn2 <- catAAbin(bn)
  bn2 <- lapply(seq_along(nf), function(z) bn2[,(nf[z]):(nt[z])])
  qm <- charToRaw('?')
  
  # gaps in or
  notgap0 <- notgap <- sapply(bn2, function(z1) {
   apply(z1, 1, function(z2) !all(z2 %in% gaps))
   })
  
  irc <- round(sum(igm > 0) * rm)
  irm <- arrayInd(order(igm,decreasing = TRUE)[1:irc],dim(igm))
  for (i in seq_len(irc)) 
    bn2[[(irm[i,][2])]][(irm[i,][1]),] <- qm
  
  orc <- round(sum(ogm > 0) * rm)
  orm <- arrayInd(order(ogm,decreasing = TRUE)[1:orc],dim(ogm))
  for (i in seq_len(orc)) 
    bn2[[(orm[i,][2])]][(orm[i,][1]),] <- qm
  

  notgap1 <- notgap <- sapply(bn2, function(z1) {
    apply(z1, 1, function(z2) !all(z2 %in% gaps))
    })
  
  if (!is.null(rgn)) {rpgn <- (rowMeans(notgap) < rgn); notgap[rpgn,] <- FALSE}
  if (!is.null(rtx)) {rptx <- (colMeans(notgap) < rtx); notgap[,rptx] <- FALSE}
  
  notgap2 <- notgap
  notgaps <- lapply(seq_len(ncol(notgap)), function(i) notgap[,i])
  
    bn2 <- mapply(function(a,b) a[b, ], bn2, notgaps)
    bn2 <- bn2[!rptx]
    
    nt <- cumsum(sapply(bn2,ncol))
    nf <- c(1, nt[-(length(nt))] + 1)
    bn2 <- catAAbin(bn2)
    bn2 <- lapply(seq_along(nf), function(z) bn2[,(nf[z]):(nt[z])])
    notgap3 <- sapply(bn2, function(z1) apply(z1, 1, function(z2) !all(z2 %in% gaps)))
    
    message(sum(!rpgn),
            " Taxa, in ",
            sum(!rptx),
            " Partitions ... Original Data:\n",
            dim(igu)[1],
            " Taxa, in ",
            dim(igu)[2],
            " Partition") 


  nigm <- (((100 - 1) / (max(igm[-outg,]) - min(igm[-outg,]))) * (igm[-outg,] - min(igm[-outg,]))) + 1
  nogm <- (((100 - 1) / (max(ogm[ outg,]) - min(ogm[ outg,]))) * (ogm[ outg,] - min(ogm[ outg,]))) + 1
  
  ngm <- rbind(nigm, nogm)
  return(list(catAAbin(bn2), bn2, notgap0, notgap1, notgap2, notgap3, ngm))

}



###########################################
############## prprt #####################

prprt <- function(x, multi = TRUE) {
# This function take trees and create a partition table specifically designed
# for the chops analysis (makeing use of the ape::prop.part function)
  xtip <- sort(unique(unlist(lapply(x, function(z) z[['tip.label']]))))
  
  xprt <- mclapply(x, function(z) {
    if (!is.null(z$node.label)) {
      zspt <- as.numeric(z$node.label[-1])
      if (any(is.na(zspt))) zspt[is.na(zspt)] <- max(zspt, na.rm = TRUE)
      } else zspt <- rep(1,(z$Nnode) - 1)
      
    brlen <- z[["edge.length"]][z[["edge"]][,2] > length(z[["tip.label"]])]
    z <- prop.part(z, check.labels = FALSE)
    whl <- length(attr(z,'labels'))
    hlf <- whl / 2
    
    if (hlf %in% lengths(z)) {
      rmv <- match(hlf,lengths(z))
      z1 <- append(z, list(setdiff(1:whl, z[[rmv]] )))
      attr(z1,'number') <- append(attr(z,'number'), attr(z,'number')[rmv])
      attr(z1, 'class') <- attr(z,'class')
      attr(z1, 'labels') <- attr(z,'labels')
      rm(z)
      z <- z1
      zspt <- append(zspt, zspt[rmv])
      brlen <- append(brlen, brlen[rmv])
    }
    
    y <- lapply(z, function(g) if (length(g) > hlf) setdiff(1:whl,g) else g)
    y <- sapply(seq_along(y), function(g) {
      unwhich(match(attr(z,'labels')[y[[g]]], xtip), length(xtip))
      })
    
    y <- y[,-1]
    attr(y, 'xtip') <- xtip
    attr(y, 'number') <- attr(z,'number')[-1]
    attr(y, 'labels') <- unwhich(match(attr(z,'labels'), xtip), length(xtip))
    attr(y, 'nodelab') <- zspt
    attr(y, 'brlen') <- brlen
    y
  }, mc.cores = myCores)
  xprt
}


# Here I am including the matrices needed for the figures in the output and removing old code


