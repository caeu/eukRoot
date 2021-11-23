# sourcing the protocols required functions
# Packeges required:
# ape, parallel, Rfast, Rfast2
source('scripts/mitoP.R')


dir.create(savepath <-
             'intermediate_data/P2/P2_loci', recursive = TRUE)

files <-
  dir(path = 'intermediate_data/P2/P2_loci',
      pattern = '\\.fast$',
      full.names = TRUE)

sTaxa <-
  scan("intermediate_data/P2/P2_loci/taxaUsed_P2.txt", what = "character")
sTaxa <- paste(sTaxa, collapse = "|")

myfragfiles <-
  lapply(files, function(x) {
    x <-
      as.matrix(read.FASTA(x, type = 'AA'))
    x[sort(grep(sTaxa, labels(x), value = TRUE)), ]
  })

namefiles <-
  gsub('^intermediate_data/P2/P2_loci/|\\.fast$', '', files)
names(myfragfiles) <- namefiles


# To create the datasets before P2 masking (Core1 and Core2)
# unifying the taxaIds
core1set <- lapply(myfragfiles, function(z) {
  updateLabel(
    x = z,
    old = rownames(z),
    new = sub("^((?:[^_]+_){2}[^_]+).*", "\\1", rownames(z))
  )
})

#1 Core1 is just direct concatenation of the initial alignments

write.FASTA(catAAbin(core1set), file = sprintf("%s/%s_%s_%s.fas",
                                     'main_datasets',
                                     'core1',
                                     length(core1set),
                                     ncol(catAAbin(core1set))))

# Write Single alignments
dir.create(core1 <- 'core1')
gaps <- charToRaw(c('-?xX'))
invisible(lapply(seq_along(core1set), function(y){
  z = core1set[[y]]
  z <- updateLabel(
    z,
    old = rownames(z),
    new = sub("^((?:[^_]+_){2}[^_]+).*", "\\1", rownames(z))
        )
  z1 <- apply(z, 1, function(z2) (sum(!(z2 %in% gaps))) > 20 )
  z <- z[z1, ]
  write.FASTA(z, file = sprintf('%s/core1_%02d.fast', core1, y))
}))


#2 Core2 is for alignments that has at least one taxa for each of the discoba
#  16 alignments do not satisfy the condition: 
#  18, 26, 33, 44, 46, 51, 57, 58, 61, 62, 63, 64, 65, 66, 68, 72
core2set <- lapply(core1set, function(z) {
  if (all(any(grepl('^ex-EG', rownames(z))),
          any(grepl('^ex-HL', rownames(z))),
          any(grepl('^ex-JA', rownames(z))))) {
    z } else NULL
})
which(sapply(core2set, is.null)) # print to console which loci to be removed for core2
core2set <- Filter(Negate(is.null), core2set)

write.FASTA(catAAbin(core2set), file = sprintf("%s/%s_%s_%s.fas",
                                     'main_datasets',
                                     'core2',
                                     length(core2set),
                                     ncol(catAAbin(core2set))))


#### Protocol2 starts here

ncolmin <- min(sapply(myfragfiles, ncol))
ncolmax <- max(sapply(myfragfiles, ncol))
minchop <- 2
maxchop <- 9
exv <- c("ex-JA", "ex-HL", "ex-EG")

dir.create(savepath <- "intermediate_data/P2/P2_chops", recursive = TRUE)
chopseq <- list()

writefs <-
  FALSE  # writefs <- TRUE #Change to TRUE to write chop files to disk
for (i in seq_along(myfragfiles)) {
  incol <- ncol(myfragfiles[[i]])
  #ifrg <- round(incol / round(incol / round(sqrt(incol)/3)))
  ifrg <-
    round((incol ^ (1 / 3) - ncolmin ^ (1 / 3)) / (ncolmax ^ (1 / 3) - ncolmin ^
                                                     (1 / 3)) * (maxchop - minchop)) + minchop + 1
  iseq <- cumsum(c(0, rep_len(incol %/% ifrg, ifrg)))
  excess <- incol - iseq[length(iseq)]
  if (excess > 0) {
    iseq <- iseq + c(rep(0, (length(iseq) - excess)), seq_len(excess))
  }
  chopseq[[i]] <- iseq
  qm <- charToRaw('-')
  gaps <- charToRaw(c('-?xX'))
  names(myfragfiles)[i] # file name
  if (writefs)
    write.FASTA(myfragfiles[[i]], file = sprintf("%s/%s_All_00.fas", savepath, names(myfragfiles)[i]))
  excv <- grep("^ex-", labels(myfragfiles[[i]])) # excavata
  for (y in exv) {
    zv <- setdiff(exv, y)
    exca <-
      grepl(paste(zv, collapse = "|"), labels(myfragfiles[[i]], value = TRUE))
    if (writefs)
      if (any(grepl("ex-", labels(myfragfiles[[i]][!(exca),]))))
        write.FASTA(myfragfiles[[i]][!(exca), ],
                    file = sprintf("%s/%s_%s_00.fas", savepath, names(myfragfiles)[i], y))
  }
  for (j in seq_len(ifrg - 1)) {
    crfi <- myfragfiles[[i]]
    crfi[excv, -((iseq[j] + 1):iseq[(j + 2)])] <- qm
    gappy <- apply(crfi, 1, function(z) {
      sum(!(z %in% gaps)) < 20
    })
    crfi <- crfi[!gappy,]
    if (writefs)
      write.FASTA(crfi, file = sprintf("%s/%s_All_%02d.fas", savepath, names(myfragfiles)[i], j))
    for (y in exv) {
      zv <- setdiff(exv, y)
      exca <-
        grepl(paste(zv, collapse = "|"), labels(crfi, value = TRUE))
      if (writefs)
        if (any(grepl("ex-", labels(crfi[!(exca),]))))
          write.FASTA(crfi[!(exca), ],
                      file = sprintf("%s/%s_%s_%02d.fas", savepath, names(myfragfiles)[i], y, j))
    }
  }
}

names(chopseq) <- namefiles
chopseqFiles <- list(chopseq = chopseq, myfragfiles = myfragfiles)

# After inferring 4 trees for each of the files above, save in the directory:
dir.create(treepath <- "intermediate_data/P2/chopTrees")
cntRuns <- 4 # number of trees per chop

# Reading the trees
# (the inferred trees are saved in: "mitop_p2_choptrees.RData")
# to use the already inferred trees jump to step in line 159
 

filesIQ <-
  dir(path = treepath,
      pattern = '\\.treefile$|\\.bestTree$',
      full.names = TRUE)
namefilesIQ <-
  gsub('^intermediate_data/P2/chopTrees/|\\.treefile$|\\.raxml.bestTree$',
       '',
       filesIQ)

mitop_p2_choptrees <-
  setNames(mclapply(filesIQ, function(z) {
    z <-
      unroot(read.tree(z))
    updateLabel(z,
                z[['tip.label']],
                sub("^([^-]{2}-[^\\W_]+_[^\\W_]+_[^\\W_]+).*$", "\\1", z[['tip.label']]))
  }, mc.cores = myCores) , namefilesIQ)


# save(mitop_p2_choptrees, file = 'intermediate_data/P2/chopTrees/mitop_p2_choptrees.RData')

## To reload the already inferred chopTrees:
load('intermediate_data/P2/chopTrees/mitop_p2_choptrees.RData')
namefilesIQ <- names(mitop_p2_choptrees)

unqprt3IQ <- prprt(mitop_p2_choptrees)


{
  xtips <-
    sort(unique(unlist(lapply(mitop_p2_choptrees, function(z)
      z[['tip.label']]))))
  exca <- grepl('ex-', xtips)
  diph <- grepl('di-', xtips)
  ubac <- grepl('ub-', xtips)
  amor <- grepl('am-', xtips)
  notEx <- rbind(amor = amor,
                 diph = diph,
                 ubac = ubac)
}

allinOne <- setNames(mclapply(1:nrow(notEx), function(v) {
  lapply(unqprt3IQ, function(z) {
    y <- attr(z, 'nodelab')
    ax <- attr(z, 'labels')
    xnod <-
      (colSums(exca & z) > 0L) &  # if it has excavata in the clade
      (colSums(notEx[v, ] &
                 z) > 0L) &
      # if it has taxa from the specified group where notEx[v,] = {amor, ubac or diph}
      (colSums((apply(notEx[-v, ], 2, any)) &
                 z) == 0L) # if does not have any taxa from other than the specified group
    if (any(xnod)) {
      m <-
        z[, xnod] # m [affected bipartitions index, the taxa in that bipartition]
      n <-
        y[xnod] # n support values for the affected bipartitions (for future)
      sx <- 1
      if (NCOL(m) == 1) {
        #if it is one then use it as is
        mx <- m
        nx <- n
        
      } else {
        # if it is more than one, reduce to one
        b <- NULL
        
        for (j in seq_len(NCOL(m) - 1)) {
          for (k in (j + 1):NCOL(m)) {
            if (all((m[, j] & m[, k]) == m[, k]) & any(m[, k] & exca)) {
              if (!all((m[, j] &
                        !exca) == (m[, k] &
                                   !exca)))
                sx <- sx + 1 # increase nesting level by one level
              if (all((m[, j] &
                       exca) == (m[, k] & exca)))
                b <- append(b, j)
              else
                b <- append(b, k)
            } else if (all((m[, j] &
                            m[, k]) == m[, j]) & any(m[, j] & exca)) {
              if (!all((m[, j] &
                        !exca) == (m[, k] &
                                   !exca)))
                sx <- sx + 1 # increase nesting level by one level
              if (all((m[, j] &
                       exca) == (m[, k] & exca)))
                b <- append(b, k)
              else
                b <- append(b, j)
            }
          }
        }
        
        
        if (is.null(b)) {
          mx <- m
          nx <- n
        } else {
          b <- unique(b)
          mx <- m[, -(b)]
          nx <- n[-(b)]
        }
        if (NCOL(mx) > 1)
          mx <- apply(mx, 1, any)
      }
    } else {
      mx <- rep(FALSE, length(ax))
      nx <- NA
      sx <- 0
    } #if (is.null(nx)) NULL else list(mx,nx)
    list(cbind(mx, ax), c(length(nx), sx), nx)
  })
}, mc.cores = myCores),
rownames(notEx))


# Create a list, each element represent a gene. each gene is a 7 dimentional array
orgGens <-
  sort(unique(sub("(^P2mito[^_]+)_.*", "\\1", namefilesIQ)))
namesNOTuse <-
  sub("(^P2mito[^_]+_[^_]+_[^_]+)_.*", "\\1", namefilesIQ)

chopsy <- setNames(mclapply(orgGens, function(z) {
  cntChop <- sum(grepl(paste0(z, "_All"), namesNOTuse)) / cntRuns
  perGene <-
    lapply(c("All", "ex-EG", "ex-HL", "ex-JA"), function(y) {
      setNames(lapply(0:(cntChop - 1), function(v) {
        which(namesNOTuse %in% sprintf("%s_%s_%02d", z, y, v))
      }), sprintf("chop%02d", 0:(cntChop - 1)))
    })
  names(perGene) <- c("All", "exEG", "exHL", "exJA")
  perGene
}, mc.cores = myCores), orgGens)


Giant <- setNames(lapply(chopsy, function(z) {
  D1 <- length(z[[1]]) # Fragments count
  D2 <- length(c("amor", "diph", "ubac"))
  D3 <- length(c("All", "exEG", "exHL", "exJA"))
  D4 <- cntRuns
  D5 <- 2
  array(NA, c(D1, D2, D3, D4, D5))
}), orgGens)


for (i in seq_along(Giant)) {
  dimChop <- paste("chop", c(0:(dim(Giant[[i]])[1] - 1)), sep = "")
  dimTree <- c('iqtr', 'iqtrG', 'rxml', 'rxmlG')
  dimnames(Giant[[i]]) <-
    list(dimChop,
         names(allinOne),
         names(chopsy[[1]]),
         dimTree,
         c("NoBP", "Scor"))
}

for (i in seq_along(Giant)) {
  d <- dim(Giant[[i]])
  for (d1 in seq_len(d[1])) {
    for (d2 in seq_len(d[2])) {
      for (d3 in seq_len(d[3])) {
        if (length(chopsy[[i]][[d3]][[d1]]) > 0) {
          for (d4 in seq_len(d[4])) {
            Giant[[i]][d1, d2, d3, d4, ] <-
              allinOne[[d2]] [[chopsy[[i]] [[d3]] [[d1]] [d4]]] [[2]]
          }
        }
      }
    }
  }
}

chopvsgroup <- mclapply(Giant, function(z) {
  dimChop <- paste("chop", c(0:(dim(z)[1] - 1)), sep = "")
  setNames(lapply(seq_len(dim(z)[1]), function(x) {
    sepa <- apply(z[x, , , , 2], c(1, 2) , function(a) {
      if (any(is.na(a)))
        NA
      else if (all(a >= 2))
        median(a)
      else
        0
    })
    sepa
  }), dimChop)
}, mc.cores = myCores)


qm <- charToRaw('x')
gaps <- charToRaw(c('-?xX'))
exv <- c("ex-EG", "ex-HL", "ex-JA")


maskedchops <-
  setNames(mclapply(seq_along(chopvsgroup), function(a1) {
    zbin <- chopseqFiles[[2]][[a1]]
    zlg <- array(0,
                 dim = c(nrow(zbin), ncol(zbin)),
                 dimnames = list(labels(zbin)))
    
    zq <- chopseqFiles[[1]][[a1]]
    zqcnt <-
      length(zq) - 2 # 5 chops with overlaps results in 7 points
    
    for (j in seq_len(zqcnt)) {
      cnf <-
        #thischop[[a1]][[j + 1]] > 0  # because No. 1 is chop0 which is full gene no choppings
        chopvsgroup[[a1]][[j + 1]] > 0  # because No. 1 is chop0 which is full gene no choppings
      
      cnf1 <- cnf[, -1]
      
      cnfC1 <- colSums(cnf1, na.rm = TRUE) > 0
      
      if (any(cnfC1, na.rm = TRUE)) {
        exvdel <- grep(paste(exv[cnfC1], collapse = "|"), labels(zbin))
        zlg[exvdel, (zq[j] + 1):zq[j + 2]] <-
          zlg[exvdel, (zq[j] + 1):zq[j + 2]] - 1
      }
      
      if (any(!cnfC1, na.rm = TRUE)) {
        exvkep <- grep(paste(exv[!cnfC1], collapse = "|"), labels(zbin))
        zlg[exvkep, (zq[j] + 1):zq[j + 2]] <-
          zlg[exvkep, (zq[j] + 1):zq[j + 2]] + 1
      }
    }
    
    zbin[zlg < 0] <- qm
    
    rgappy <-
      apply(zbin, 1, function(z) {
        sum(!(z %in% gaps)) < 2
      })
    
    cgappy <-
      apply(zbin, 2, function(z) {
        sum(!(z %in% gaps)) < 2
      })
    
    if (length(zbin[!rgappy, !cgappy]) > 0)
      mskd <- zbin[!rgappy, !cgappy]
    else
      mskd <- NULL
    
    mskd
    
  }, mc.cores = myCores), orgGens)



notcatbin <-
  setNames(lapply(seq_along(maskedchops), function(a1) {
    zbin <- maskedchops[[a1]]
    
    if (is.null(zbin)) {
      fzbin <- list(NULL, NULL)
      
    } else {
      leng <- nrow(zbin)
      exwant <-
        vapply(exv, function(z) {
          grepl(z, labels(zbin))
        }, logical(leng))
      
      if (any(apply(exwant, 2, any))) {
        w1 <- apply(exwant, 1, any)
        w2 <- match(zbin[w1,], gaps, nomatch = FALSE)
        w3 <- matrix(!w2, ncol = ncol(zbin))
        w4 <- colAny(w3)
        zbin1 <- zbin[, w4]
        
        if (all(apply(exwant, 2, any))) {
          # This is to keep the alingment if it has at least one discoba for each of the three discoba
          #if (any(apply(exwant, 2, any))) { # This is keep the alignment if it has at least one discoba from any of the three discoba
          w1 <- apply(exwant, 1, any)
          w2 <- match(zbin[w1,], gaps, nomatch = FALSE)
          w3 <- matrix(!w2, ncol = ncol(zbin))
          w4 <- colAny(w3)
          zbin2 <- zbin[, w4]
        } else {
          zbin2 <- NULL
        }
        fzbin <- list(zbin1, zbin2)
      } else {
        fzbin <- list(NULL, NULL)
      }
      
      fzbin
    }
  }), orgGens)

notcatCore1 <- lapply(notcatbin, function(z)
  z[[1]])
notcatCore2 <- lapply(notcatbin, function(z)
  z[[2]])

tocatCore1 <- Filter(Negate(is.null), notcatCore1)
tocatCore2 <- Filter(Negate(is.null), notcatCore2)


fastaP2 <-
  function(x,
           id = c('core1_P2', 'core2_P2'),
           pathsign = "main_datasets") {
    id <- match.arg(id)
    if (length(x) > 0) {
      x <- lapply(x, function(z) {
        updateLabel(
          x = z,
          old = rownames(z),
          new = sub("^((?:[^_]+_){2}[^_]+).*", "\\1", rownames(z))
        )
      })
      bincnt <- length(x)
      catbin <- catAAbin(x)
      catlen <- ncol(catbin)
      write.FASTA(catbin,
                  file = sprintf("%s/%s_%s_%s.fas",
                                 pathsign,
                                 id,
                                 bincnt,
                                 catlen))
    } else
      message ('Nothing to print!')
  }

# To print the concatenated fasta files to the
fastaP2(tocatCore1)
fastaP2(tocatCore2, id = "core2")


# To print the single alignments after P2 (core1_p2, core2_p2)
# 22 genes do not meet the requirement for core2
cat(which(!sapply(seq_along(notcatCore1), function(z) identical(notcatCore1[[z]], notcatCore2[[z]]))), sep = ', ')
# 12, 18, 26, 29, 31, 33, 41, 44, 46, 51, 54, 56, 57, 58, 61, 62, 63, 64, 65, 66, 68, 72

dir.create(core1p2 <- 'core1_p2')
gaps <- charToRaw(c('-?xX'))
invisible(lapply(seq_along(tocatCore1), function(y){
  z = tocatCore1[[y]]
  z <- updateLabel(
    z,
    old = rownames(z),
    new = sub("^((?:[^_]+_){2}[^_]+).*", "\\1", rownames(z))
        )
  z1 <- apply(z, 1, function(z2) (sum(!(z2 %in% gaps))) > 20 )
  z <- z[z1, ]
  write.FASTA(z, file = sprintf('%s/core1p2_%02d.fast', core1p2, y))
}))




# Concordance iqtree gcf scf
files <-  dir(path = '../eukroot_concordance', pattern = "cf_scf\\.cf\\.tree$", recursive = TRUE, full.names = TRUE)
filenames <- basename(files)

cftrees <- setNames(lapply(files, function(z) {
  unroot(read.tree(z))
}), filenames)


setNames(lapply(cftrees, function(z) {
  prtnm <- prop.part(z)
  
  myrootex <- grepl('di-|am-', z[["tip.label"]])
  myrootam <- grepl('di-|ex-', z[["tip.label"]])
  myrootdi <- grepl('ex-|am-', z[["tip.label"]])
  
  if (is.monophyletic(z, which(myrootex))) {
    myrt <- "ex"
    mypos <- Position(function(y) identical(y, which(myrootex)), prtnm, nomatch = 0)
    if (mypos > 0) {
      myret <- z[["node.label"]][mypos]
    } else {
      mypos <- Position(function(y) identical(y, which(!myrootex)), prtnm, nomatch = 0)
      myret <- z[["node.label"]][mypos]
    }} else if (is.monophyletic(z, which(myrootam))) {
    myrt = 'am'
    mypos <- Position(function(y) identical(y, which(myrootam)), prtnm, nomatch = 0)
    if (mypos > 0) {
      myret <- z[["node.label"]][mypos]
    } else {
      mypos <- Position(function(y) identical(y, which(!myrootam)), prtnm, nomatch = 0)
      myret <- z[["node.label"]][mypos]
    }} else if (is.monophyletic(z, which(myrootdi))) {
    myrt <- 'di'
    mypos <- Position(function(y) identical(y, which(myrootdi)), prtnm, nomatch = 0)
    if (mypos > 0) {
      myret <- z[["node.label"]][mypos]
    } else {
      mypos <- Position(function(y) identical(y, which(!myrootdi)), prtnm, nomatch = 0)
      myret <- z[["node.label"]][mypos]
    }}
  c(myrt,myret)
}), filenames)
