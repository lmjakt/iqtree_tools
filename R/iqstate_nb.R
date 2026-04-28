
nib.enc <- c("-"=0, A=1, C=2, M=3, G=4, R=5, S=6, V=7, T=8, W=9, Y=10, H=11, K=12, D=13, B=14, N=15)
trans.enc <- sapply(names(nib.enc), function(x){ paste(x, names(nib.enc), sep="..") })
rownames(trans.enc) <- nib.enc


nb.state.read <- function(fn){
    con <- file(fn, open="rb")
    f.size <- file.size(fn)
    magic <- vector(mode='integer', length=2)
    magic <- readBin(con, magic, n=2, size=4)
    seq.l <- readBin(con, "integer", n=1, size=4)
    nodes.n <- readBin(con, "integer", n=1, size=4)
    ## we can read in the full set of data after this.
    ## lets read the sequences as a matrix of ints
    ## with nibble encoding the number of integers used per sequence is:
    seq.bytes <- seq.l %/% 8 + ifelse(seq.l %% 8 > 0, 1, 0)
    nb.seq <- matrix( readBin(con, "integer", n=seq.bytes * nodes.n, size=4), ncol=nodes.n )
    nodes <- readBin(con, "character", n=nodes.n)
    colnames(nb.seq) <- nodes
    close(con)
    nb.seq
}

## convert an integer nibble encoded sequence to a character vector
## nb is an integer vector
nb.seqToChar <- function(nb){
    ## Note that in this encoding 0 indicates a gap, 15 an N.
    ## use X for other uncertainty symbols to begin with.
    nb.ch <- rep("X", 16) 
    nb.ch[ 1+c(0, 1, 2, 4, 8, 15) ] <- c("-", "A", "C", "G", "T", "N")
    nb.unc <- t(sapply(0:7, function(i){
        bitwAnd( bitwShiftR( nb, 28 - i*4 ), 15 )
    }))
    ch <- nb.ch[ nb.unc+1 ]
    paste(ch, collapse="")
}

read.fa <- function(fn){
    lines <- readLines(fn)
    id.i <- grep("^>", lines)
    beg <- id.i + 1
    end <- c(id.i - 1, length(lines))[-1]
    seq <- sapply(seq_along(id.i), function(i){ paste(lines[ beg[i]:end[i] ], collapse="") })
    names(seq) <- sub("^>([^ ]+).*", "\\1", lines[id.i])
    seq
}

state.to.seq <- function(fn){
    lines <- readLines(fn)
    lines <- lines[!grepl("^#", lines)]
    lines <- lines[-1] ## remove the header line
    cols <- strsplit(lines, "\t")
    nucs <- sapply(cols, function(x){ x[c(1,3)] })
    ## Note: nesting tapply inside c() converts the "array"
    ## object to a character object. This makes it more similar
    ## to the other objects obtained by the functions provided here.
    c(tapply( nucs[2,], nucs[1,], paste, collapse=""))
}

## Obtain an 8 bit representation of transition from nib1 to nib2
nibble.trans <- function(nib1, nib2){
    nb.tr <- t(sapply(0:7, function(i){
        b1 <- bitwAnd( bitwShiftR( nib1, 28 - i*4 ), 15 )
        b2 <- bitwAnd( bitwShiftR( nib2, 28 - i*4 ), 15 )
        bitwOr(bitwShiftL(b1, 4), b2)
    }))
    as.integer(nb.tr)
}

## if full, return all 256 possible transitions. Otherwise return
## the counts for (A,C,T,G,-,N)
count.nibble.trans <- function(nib1, nib2, full=FALSE){
    ntr <- nibble.trans(nib1, nib2)
    tr.n <- as.integer(table( c(0:255, ntr) ) - 1)
    names(tr.n) <- trans.enc
    if(full)
        return(tr.n)
    sel.n <- c("-", "A", "C", "G", "T", "N")
    i <- sapply(nib.enc[sel.n], function(i){ i * 16 + nib.enc[sel.n] })
    ## note nib.enc and trans.enc are defined at the top of this file
    tr.n[ i+1 ]
}
