require("ape")
source("iqstate_nb.R")

nb.fn <- list.files("../parse_state_file/test", pattern="nb_states$", full.names=TRUE)
fa.fn <- list.files("../parse_state_file/test", pattern="fas$", full.names=TRUE)
state.fn <- list.files("../parse_state_file/test", pattern="state$", full.names=TRUE)
tree.fn <- list.files("../parse_state_file/test", pattern="treefile$", full.names=TRUE)
tr.n.fn <- list.files("../extract_transitions/test", pattern="tsv$", full.names=TRUE)

fa.seq <- read.fa(fa.fn)
states.seq <- state.to.seq(state.fn)
txt.seq <- c(states.seq, fa.seq)

tree <- read.tree( tree.fn )
plot(tree, show.node.label=TRUE)

nibbles <- nb.state.read(nb.fn)
bin.seq <- apply(nibbles, 2, nb.seqToChar)

### Compare the sequences. Note that bin.seq is likely to be slighly longer.
all(txt.seq == substr( bin.seq[names(txt.seq)], 1, nchar(txt.seq) )) ## TRUE

### The encoding does seem to be correct.

### Count transitions:
## this just gives us a vector of the transitions:
n1.n2 <- nibble.trans(nibbles[,'Node1'], nibbles[,'Node2'])
n1.n2.n <- count.nibble.trans(nibbles[,'Node1'], nibbles[,'Node2'], full=TRUE)

## we can compare it to what extract_transitions gave us.
trans.n <- read.table(tr.n.fn, sep="\t")
colnames(trans.n) <- c("i", "j", "d", "p", "c", trans.enc)

par(mfrow=c(2,1))
barplot( as.numeric(trans.n[4,-(1:5)]), las=2, cex.names=0.8, names.arg=trans.enc)
barplot(n1.n2.n, las=2, cex.names=0.8)

all(trans.n[4,-(1:5)] == n1.n2.n) ## FALSE
## but:

which(trans.n[4,-(1:5)] != n1.n2.n) ## 1
## and the first entry is gap to gap. There will be additional gap to gap
## counts using the function here due to gaps inserted at the end. This
## is thus reasonable and suggests that that the counting is done correctly.

