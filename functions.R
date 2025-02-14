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

circ_arcs <- function(sid, min_del_len = 0, alpha = 0.2) {
    circos.clear()
    col_text <- "grey40"
    circos.par("track.height" = 0.1, gap.degree = 1, cell.padding = c(0, 0, 0, 0), start.degree = 90)
    circos.initialize(as.factor(rep("mtDNA", 16569)), x = 1:16569)
    #circos.initialize("mtDNA", xlim = c(1,16569))
    
    circos.track(ylim = c(0,1), panel.fun = function(x, y) {
        chr = CELL_META$sector.index
        xlim = CELL_META$xlim
        ylim = CELL_META$ylim
        circos.text(mean(xlim), mean(ylim), chr, cex = 1, col = col_text,
                    facing = "bending.inside", niceFacing = TRUE)
    }, bg.col = "grey90", bg.border = FALSE, track.height = 0.06)
    
    # axis
    brk <- c(1, seq(501, 16001, by = 500), 16569)
    circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
        circos.axis(h = "top", major.at = brk, labels.cex = 1, col = col_text,
                    labels.col = col_text, lwd = 0.7, labels.facing = "clockwise")
    }, bg.border = FALSE)
    
    # gc content
    circos.genomicTrack(data = mt_gc_content, ylim = c(0, 1), panel.fun = function(region, value, ...) {
        circos.genomicLines(region, value, type = "l", col = "lightblue", lwd = 0.6, area = TRUE)
        circos.segments(x0 = 0, x1 = 16569, y0 = 0.25, y1 = 0.25, lwd = 0.6, lty = "11", col = "grey60")
        circos.segments(x0 = 0, x1 = 16569, y0 = 0.50, y1 = 0.50, lwd = 0.6, lty = "11", col = "grey60")
        circos.segments(x0 = 0, x1 = 16569, y0 = 0.77, y1 = 0.75, lwd = 0.6, lty = "11", col = "grey60")
    }, track.height = 0.1, bg.border = TRUE)
    circos.yaxis(at = c(0, 1), labels.cex = 0.6, lwd = 0, tick.length = 0, labels.col = col_text, col = col_text)
    
    circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
        Xleft <- mt_genes_df$start
        Xright <- mt_genes_df$end
        Ybottom <- ifelse(mt_genes_df$drop, 0, 0.5)
        Ytop <- ifelse(mt_genes_df$drop, 0.5, 1)
        Col <- case_when(mt_genes_df$gene_biotype == 'protein_coding' ~ "grey50",
                         mt_genes_df$gene_biotype == 'Mt_rRNA' ~ "grey80",
                         mt_genes_df$gene_biotype == 'Mt_tRNA' ~ "white")
        circos.rect(xleft = Xleft, ybottom = Ybottom,
                    xright = Xright, ytop = Ytop,
                    col = Col, border = col_text)
    }, bg.border = NA, track.height = 0.05)
    
    #circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    #    for (i in 1:nrow(mt_genes_df)) {
    #        circos.arrow(x1 = mt_genes_df[i,]$start, x2 = mt_genes_df[i,]$end,
    #                     arrow.position = ifelse(mt_genes_df[i,]$strand == '+', "start", "end"),
    #                     col = ifelse(mt_genes_df[i,]$strand == '+', "red", "green"),
    #                     border = "black")
    #    }
    #}, bg.border = NA, track.height = 0.2)
    
    circos.genomicLabels(mt_genes_df, labels.column = 5, col = col_text, line_lwd = 0.5,
        line_col = "grey80", side = "inside", connection_height = 0.05, labels_height = 0.04)
    
    del_a <- Del %>% mutate(del_size = ifelse(del_start <= del_end, 
                             del_end - del_start, 
                             del_end + (16569 - del_start))) %>%
        filter(del_size >= min_del_len, sample_id == sid) %>%
        mutate(seqnames = factor("mtDNA")) %>%
        select(seqnames, start = del_start, end = del_start, read_length:del_size)
    del_b <- Del %>% mutate(del_size = ifelse(del_start <= del_end, 
                             del_end - del_start, 
                             del_end + (16569 - del_start))) %>%
        filter(del_size >= min_del_len, sample_id == sid) %>%
        mutate(seqnames = factor("mtDNA")) %>%
        select(seqnames, start = del_end, end = del_end, read_length:del_size)
    
    circos.genomicLink(del_a, del_b, col = add_transparency("darkred", 1-alpha), border = NA)
    
    title(paste0("Deletions in ", Metadata[Metadata$sample_id == sid,]$tissue[1]))
    
    circos.clear()
}

circ_del_cov <- function(sid, min_del_len = 0, alpha = 0.5) {
    circos.clear()
    col_text <- "grey40"
    circos.par(circle.margin = 0.2, "track.height" = 0.1, gap.degree = 1, cell.padding = c(0, 0, 0, 0), start.degree = 90)
    circos.initialize(as.factor(rep("mtDNA", 16569)), x = 1:16569)
    
    circos.track(ylim = c(0,1), panel.fun = function(x, y) {
        chr = CELL_META$sector.index
        xlim = CELL_META$xlim
        ylim = CELL_META$ylim
        circos.text(mean(xlim), mean(ylim), chr, cex = 1, col = col_text,
                    facing = "bending.inside", niceFacing = TRUE)
    }, bg.col = "grey90", bg.border = FALSE, track.height = 0.06)
    
    # axis
    brk <- c(1, seq(501, 16001, by = 500), 16569)
    circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
        circos.axis(h = "top", major.at = brk, labels.cex = .8, col = col_text,
                    labels.col = col_text, lwd = 0.7, labels.facing = "clockwise")
    }, bg.border = FALSE)
    
    circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
        Xleft <- mt_genes_df$start
        Xright <- mt_genes_df$end
        Xmid <- mt_genes_df %>% filter(big) %>% mutate(midp = (start+end)/2) %>% pull(midp)
        Ymid <- ifelse(mt_genes_df %>% filter(big) %>% pull(drop), 0.25, 0.75)
        Labels <- mt_genes_df %>% filter(big) %>% pull(gene_name)
        ColText <- ifelse((mt_genes_df %>% filter(big) %>% pull(gene_biotype)) == "protein_coding", "white", "black")
        Ybottom <- ifelse(mt_genes_df$drop, 0, 0.5)
        Ytop <- ifelse(mt_genes_df$drop, 0.5, 1)
        Col <- case_when(mt_genes_df$gene_biotype == 'protein_coding' ~ "grey50",
                         mt_genes_df$gene_biotype == 'Mt_rRNA' ~ "grey85",
                         mt_genes_df$gene_biotype == 'Mt_tRNA' ~ "white")
        circos.rect(xleft = Xleft, ybottom = Ybottom,
                    xright = Xright, ytop = Ytop,
                    col = Col, border = col_text)
        circos.text(x = Xmid, y = Ymid, labels = Labels, cex = .6, 
                    col = ColText,
                    facing = "bending.inside", niceFacing = TRUE)
    }, bg.border = NA, track.height = 0.06)
    
    seqinfo <- Seqinfo(
        seqnames = c("mtDNA"),
        seqlengths = c(16569),
        isCircular = c(TRUE),
        genome = "hg38"
    )

    non_hanging <- Del %>%
        mutate(del_size = ifelse(del_start <= del_end, 
                             del_end - del_start, 
                             del_end + (16569 - del_start))) %>%
        filter(del_size >= min_del_len, sample_id == sid) %>%
        mutate(seqnames = factor("mtDNA")) %>%
        select(seqnames, start = del_start, end = del_end) %>%
        filter(start < end) %>%
        makeGRangesFromDataFrame(seqinfo = seqinfo)

    hanging <- Del %>%
        mutate(del_size = ifelse(del_start <= del_end, 
                             del_end - del_start, 
                             del_end + (16569 - del_start))) %>%
        filter(del_size >= min_del_len, sample_id == sid) %>%
        mutate(seqnames = factor("mtDNA")) %>%
        select(seqnames, start = del_start, end = del_end) %>%
        filter(start > end)
    hang_left <- hanging %>%
        mutate(end = 16569) %>%
        makeGRangesFromDataFrame(seqinfo = seqinfo)
    hang_right <- hanging %>%
        mutate(start = 1) %>%
        makeGRangesFromDataFrame(seqinfo = seqinfo)
    del_cov <- as.numeric(coverage(c(non_hanging, hang_left, hang_right))$mtDNA)
        
    circos.track(factor(rep("mtDNA", 16569)), x = 1:16569, y = del_cov,
                 panel.fun = function(x, y) {
                     circos.lines(x, y, col = add_transparency("brown2", 1-alpha), border = "brown2", lwd = 0.75, type = "l", area = TRUE)
                 }
    )
    
    del_a <- Del %>% mutate(del_size = ifelse(del_start <= del_end, 
                             del_end - del_start, 
                             del_end + (16569 - del_start))) %>%
        filter(del_size >= min_del_len, sample_id == sid) %>%
        mutate(seqnames = factor("mtDNA")) %>%
        select(seqnames, start = del_start, end = del_start, read_length:del_size)
    del_b <- Del %>% mutate(del_size = ifelse(del_start <= del_end, 
                             del_end - del_start, 
                             del_end + (16569 - del_start))) %>%
        filter(del_size >= min_del_len, sample_id == sid) %>%
        mutate(seqnames = factor("mtDNA")) %>%
        select(seqnames, start = del_end, end = del_end, read_length:del_size)
    
    circos.genomicLink(del_a, del_b, col = add_transparency("darkred", 1-alpha), border = NA)
    title(paste0("Deletions in ", Metadata[Metadata$sample_id == sid,]$tissue[1]))
    
    circos.clear()
}
