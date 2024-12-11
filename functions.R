### FUNCTIONS
# pos 8001 in chrM corresponds to pos 1 in shifted chrM
# this function transforms coordinate from original chrM into the shifted
# reference
s8k <-function(x) {
    if (x > 16569) return(NA)
    if (x < 1) return(NA)
    if (x > 8000) return(x - 8000)
    if (x <= 8000) return(x + 8569)
    return(NA)
}
# this one from shifted ref to original
us8k <-function(x) {
    if (x > 16569) return(NA)
    if (x < 1) return(NA)
    if (x < 8570) return(x + 8000)
    if (x >= 8570) return(x - 8569)
    return(NA)
}

# From chatGPT, which is often smarter than me
# UNTESTED
shifted_to_rcrs <- function(positions, shift = 8000, genome_length = 16569) {
  zero_based_positions <- positions - 1
  rcrs_positions <- (zero_based_positions + shift) %% genome_length
  return(rcrs_positions + 1)
}

# Calculate genotype from base A, C, G, T counts
gt <- function(counts) {
    if (length(counts) != 4) {
        stop("The counts argument must be of length 4")
    }
    props <- counts/sum(counts)
    return(props)
}
err_het <- function(freqs) {
    sum(abs(sort(freqs)[3:4] - 0.5))
}

scale2 <- function(x, na.rm=FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)

stars.pvalue <- function(pval) {
    stars <- rep('', length(pval))
    stars[pval<=0.1] <- "."
    stars[pval<=0.05] <- "*"
    stars[pval<=0.01] <- "**"
    stars[pval<=0.001] <- "***"
    return(stars)
}

dummify <- function(m, reduce.binary=FALSE) {
    factor.cols <- sapply(m, function(x) is.factor(x) || is.character(x) || is.logical(x))
    if (sum(factor.cols) == 0) return(m)
    d <- lapply(colnames(m)[factor.cols], function(col.name) {
               mm <- model.matrix.lm(~.+0, as.data.frame(m[,col.name, drop=FALSE]), na.action=na.pass)
               if (ncol(mm)==2 && reduce.binary==TRUE) {
                   mm <- mm[,1,drop=FALSE]
               }
               return(mm)
    }) %>% do.call('cbind', .)
    return(cbind(m[!factor.cols], d))
}

assocVars <- function(M, id.col=1, method=c("pearson", "kendall", "spearman"), use="everything", reduce.binary=FALSE,
                      pvalues=FALSE, adjust="holm", alpha=.05) {
    # Test variance
    M <- as.data.frame(M) %>% mutate_if(is.character, as.factor)
    # Get IDs
    ids <- M[[id.col]]
    M <- M[-id.col]
    # Dummify
    m <- dummify(M, reduce.binary=reduce.binary)
    m <- as.data.frame(m)
    zeroVar <- (apply(mutate_if(M, is.factor, as.integer), 2, var, na.rm=TRUE) == 0)
    apply(mutate_if(M, is.factor, as.integer), 2, sd, na.rm=TRUE)
    if (any(zeroVar)) stop("Variables ", paste(colnames(M)[zeroVar]), " have zero variance")
    rownames(m) <- ids
    if (pvalues == TRUE) {
        cor.mat <- psych::corr.test(m, method=method, use=use, adjust=adjust, alpha=alpha)$p
    } else {
        cor.mat <- cor(m, method=method, use=use)
    }
    diag(cor.mat) <- NA
    return(cor.mat)
}

transf.betareg <- function(y) {
    n.obs <- sum(!is.na(y))
    (y * (n.obs - 1) + 0.5) / n.obs
}

