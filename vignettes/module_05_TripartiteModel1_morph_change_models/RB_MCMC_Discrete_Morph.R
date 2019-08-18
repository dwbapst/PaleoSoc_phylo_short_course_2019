## ----global_options, eval = TRUE, include=TRUE---------------------------


## ----eval = TRUE---------------------------------------------------------
    example <- 1.0


## ---- include=TRUE, eval = TRUE------------------------------------------
    morpho <- readDiscreteCharacterData("data/Cinctans.nex")


## ---- include=TRUE, eval = TRUE------------------------------------------
    taxa <- morpho.names()
    num_taxa <- morpho.size() 
    num_branches <- 2 * num_taxa - 2


## ---- include=TRUE, eval = TRUE------------------------------------------
    mvi = 1
    mni = 1


## ---- include=TRUE, eval = TRUE------------------------------------------
    br_len_lambda ~ dnExp(0.2)
    moves[mvi++] = mvScale(br_len_lambda, weight=2)


## ---- include=TRUE, eval=FALSE-------------------------------------------
## library(ggplot2)
## draws <- rexp(10000, .2)
## hist(draws)


## ---- include=TRUE, eval = TRUE------------------------------------------
    
    phylogeny ~ dnUniformTopologyBranchLength(taxa, branchLengthDistribution=dnExponential(br_len_lambda))
    moves[mvi++] = mvNNI(phylogeny, weight=num_branches/2.0)
    moves[mvi++] = mvSPR(phylogeny, weight=num_branches/10.0)
    moves[mvi++] = mvBranchLengthScale(phylogeny, weight=num_branches)
    tree_length := phylogeny.treeLength()


## ---- include=TRUE, eval = TRUE------------------------------------------
    alpha_morpho ~ dnUniform( 0, 1E6 )
    rates_morpho := fnDiscretizeGamma( alpha_morpho, alpha_morpho, 4 )
    #Moves on the parameters to the Gamma distribution.
    moves[mvi++] = mvScale(alpha_morpho, lambda=1, weight=2.0)


## ---- include=TRUE, eval=FALSE-------------------------------------------
## library(ggplot2)
## alpha_morpho <- runif(1, 0, 1E6 )
## 
## draws <- rgamma(1000, shape = alpha_morpho, rate = alpha_morpho)
## hist(draws)


## ---- include=TRUE, eval = TRUE------------------------------------------
n_max_states <- 7
idx = 1
morpho_bystate[1] <- morpho
for (i in 2:n_max_states) {
    # make local tmp copy of data
    # only keep character blocks with state space equal to size i
    morpho_bystate[i] <- morpho
    morpho_bystate[i].setNumStatesPartition(i)
	# get number of characters per character size wth i-sized states
    nc = morpho_bystate[i].nchar()
    # for non-empty character blocks
    if (nc > 0) {
        # make i-by-i rate matrix
        q[idx] <- fnJC(i)
# create model of evolution for the character block
        m_morph[idx] ~ dnPhyloCTMC( tree=phylogeny,
                                    Q=q[idx],
                                    nSites=nc,
                                    siteRates=rates_morpho,
                                    type="Standard")

        # attach the data
	    m_morph[idx].clamp(morpho_bystate[i])

        # increment counter
        idx = idx + 1
idx
}
}


## ---- include=TRUE, eval = TRUE------------------------------------------
    mymodel = model(phylogeny)


## ---- include=TRUE, eval = TRUE------------------------------------------
    monitors[mni++] = mnModel(filename="output/mk_gamma.log", printgen=10)


## ---- include=TRUE, eval = TRUE------------------------------------------
    monitors[mni++] = mnFile(filename="output/mk_gamma.trees", printgen=10, phylogeny)


## ---- include=TRUE, eval = TRUE------------------------------------------
    monitors[mni++] = mnScreen(printgen=100)


## ---- include=TRUE, eval = TRUE------------------------------------------
    mymcmc = mcmc(mymodel, monitors, moves, nruns=2, combine="mixed")


## ---- include=TRUE, eval = TRUE------------------------------------------
    mymcmc.run(generations=10000, tuningInterval=200)


## ---- include=TRUE, eval = TRUE------------------------------------------
    q()


## ------------------------------------------------------------------------
knitr::purl("module_05_TripartiteModel1_morph_change_models/RB_Discrete_Morphology/RB_MCMC_Discrete_Morph.Rmd")


##     source("scripts/mk_gamma.Rev")

