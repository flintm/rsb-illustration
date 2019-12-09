# Exploration of Cartesian product to identify set of alternative systems, A
# and performance paramter values, Theta.
# Written by Madeleine Flint, 2019-11-20 to 2019-12-05

# Define subsets --------
n <- 4 # number of subsets (i.e., S, F, L, E)
S <- factor(paste0('s',1:2))
nS = length(S)
F <- factor(paste0('f',1:4))
nF = length(F)
L <- factor(paste0('l',1:11), levels = paste0('l',1:11))
nL = length(L)
E <- factor(paste0('e',1:32), levels = paste0('e',1:32))
nE = length(E)
nA = nS*nF*nL*nE
combine <- function(..., prefix = "", sep = "_") {
  paste0(prefix, levels(interaction(..., sep = sep)))
}

# Define A = S x F x L x E ----------
A <- combine(S,F,L,E) # a warning will appear but can be ignored
length(A) == nA
m.A <- matrix(unlist(strsplit(A,"_")), nrow = nA, ncol = n, byrow = TRUE)
head(m.A)

# Define performance parameter values for subsystems, Theta
theta <- factor(seq(from=0,to=1,length.out=5))
nk = length(theta)
nK = nk^4
Theta <- combine(theta,theta,theta,theta)
m.Theta <- matrix(unlist(strsplit(Theta,"_")), nrow = nK, ncol = n, byrow = TRUE)
head(m.Theta)

# Find numbers associated with hypothetical systems in Fig. 1--------
ex.Match <- function(row, example){
	all(row==example)
}
n.ex <- 3 # a, b, c
ex.A <- matrix(c("s1", "f3", "l6", "e7", 
               "s1", "f3", "l11", "e21", 
               "s1", "f4", "l6", "e14"), 
             nrow = n.ex, ncol = n, byrow = TRUE)
ex.Theta <- matrix(c("0", "0", "0", "0", 
                     "1", "0", "0.25", "0.25", 
                     "0", "1", "1", "1"), 
                   nrow = n.ex, ncol = n, byrow = TRUE)
ex.k     <- sapply(1:n.ex, function(i) which(apply(m.A, MARGIN = 1, ex.Match, ex.A[i,])))
ex.kappa <- sapply(1:n.ex, function(i) which(apply(m.Theta, MARGIN = 1, ex.Match, ex.Theta[i,])))

print(matrix(c(ex.k, ex.kappa), nrow = n.ex, ncol = 2))