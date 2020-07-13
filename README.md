<!-- README.md is generated from README.Rmd. Please edit that file -->

SpiecEasi
=========

[![Build
Status](https://travis-ci.org/zdk123/SpiecEasi.svg?branch=pulsar)](https://travis-ci.org/zdk123/SpiecEasi)

Sparse InversE Covariance estimation for Ecological Association and
Statistical Inference

This package will be useful to anybody who wants to infer graphical
models for all sorts of compositional data, though primarily intended
for microbiome relative abundance data (generated from 16S amplicon
sequence data). It also includes a generator for \[overdispersed, zero
inflated\] multivariate, correlated count data. Please see the paper
published in [PLoS Comp
Bio](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004226).

One small point on notation: we refer to the method as “SPIEC-EASI” and
reserve “SpiecEasi” for this package.

TOC
---

1.  [Installation](#installation)
2.  [News](#news)
3.  [Basic Usage](#basic-usage)
4.  [American Gut Data](#analysis-of-american-gut-data)
5.  [Using phyloseq](#working-with-phyloseq)
6.  [Learning latent variable graphical
    models](#learning-latent-variable-graphical-models)
7.  [Cross Domain SPIEC-EASI](#cross-domain-interactions)
8.  [pulsar & batch
    options](#pulsar-parallel-utilities-for-model-selection)
9.  [Troubleshooting](#troubleshooting)

Installation
------------

I assume that all auxiliary packages are already installed - for example
pulsar, huge, MASS, etc. If you get an unexpected error, you may need to
download and install a missing dependency.

From an interactive R session:

    library(devtools)
    install_github("zdk123/SpiecEasi")
    library(SpiecEasi)

News
----

The latest SpiecEasi (version 1.0.0 and higher) now uses the [pulsar
package](https://cran.r-project.org/package=pulsar) for stability-based
model selection. The methods are similar to previous releases, but
contains some additional methods for [speeding up
computations](#pulsar-parallel-utilities-for-model-selection)

The input arguments have changed slightly (but are backwards compatible)
but the data structure that is returned from `spiec.easi` has changed.

The output to spiec.easi-fit models structure can be easily processed
using new getter functions. See `?getOptInd` for usage.

You can revert to the previous release
([0.1.4](https://github.com/zdk123/SpiecEasi/releases/tag/v0.1.4)) to
avoid code-breaking changes.

Basic Usage
-----------

Lets simulate some multivariate data under zero-inflated negative
binomial model, based on (high depth/count) round 1 of the American gut
project, with a sparse network. The basic steps are

1.  load the data and normalize counts to to common scale (min depth)
2.  fit count margins to the model
3.  generate a synthetic network
4.  generate some synthetic data
5.  clr transformation
6.  inverse covariance estimation along a lambda (sparsity) path
7.  stability selection using the StARS criterion
8.  evaluate performance

Obviously, for real data, skip 1-4.

    data(amgut1.filt)
    depths <- rowSums(amgut1.filt)
    amgut1.filt.n  <- t(apply(amgut1.filt, 1, norm_to_total))
    amgut1.filt.cs <- round(amgut1.filt.n * min(depths))

    d <- ncol(amgut1.filt.cs)
    n <- nrow(amgut1.filt.cs)
    e <- d

Synthesize the data

    set.seed(10010)
    graph <- make_graph('cluster', d, e)
    Prec  <- graph2prec(graph)
    Cor   <- cov2cor(prec2cov(Prec))

    X <- synth_comm_from_counts(amgut1.filt.cs, mar=2, distr='zinegbin', Sigma=Cor, n=n)

the main SPIEC-EASI pipeline: Data transformation, sparse inverse
covariance estimation and model selection

    se <- spiec.easi(X, method='mb', lambda.min.ratio=1e-2, nlambda=15)
    # Applying data transformations...
    # Selecting model with pulsar using stars...
    # Fitting final estimate with mb...
    # done

examine ROC over lambda path and PR over the stars index for the
selected graph

    huge::huge.roc(se$est$path, graph, verbose=FALSE)
    stars.pr(getOptMerge(se), graph, verbose=FALSE)
    # stars selected final network under: getRefit(se)

The above example does not cover all possible options and parameters.
For example, other generative network models are available, the
lambda.min.ratio (the scaling factor that determines the minimum
sparsity/lambda parameter) shown here might not be right for your
dataset, and its possible that you’ll want more repetitions (number of
subsamples) for StARS.

Analysis of American Gut data
-----------------------------

Now let’s apply SpiecEasi directly to the American Gut data. Don’t
forget that the normalization is performed internally in the
`spiec.easi` function. Also, we should use a larger number of stars
repetitions for real data. We can pass in arguments to the inner stars
selection function as a list via the parameter `pulsar.params`. If you
have more than one processor available, you can also supply a number to
`ncores`. Also, let’s compare results from the MB and glasso methods as
well as SparCC (correlation).

    se.mb.amgut <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=1e-2,
                              nlambda=20, pulsar.params=list(rep.num=50))
    se.gl.amgut <- spiec.easi(amgut1.filt, method='glasso', lambda.min.ratio=1e-2,
                              nlambda=20, pulsar.params=list(rep.num=50))
    sparcc.amgut <- sparcc(amgut1.filt)
    ## Define arbitrary threshold for SparCC correlation matrix for the graph
    sparcc.graph <- abs(sparcc.amgut$Cor) >= 0.3
    diag(sparcc.graph) <- 0
    library(Matrix)
    sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)
    ## Create igraph objects
    ig.mb     <- adj2igraph(getRefit(se.mb.amgut))
    ig.gl     <- adj2igraph(getRefit(se.gl.amgut))
    ig.sparcc <- adj2igraph(sparcc.graph)

Visualize using igraph plotting:

    library(igraph)
    ## set size of vertex proportional to clr-mean
    vsize    <- rowMeans(clr(amgut1.filt, 1))+6
    am.coord <- layout.fruchterman.reingold(ig.mb)

    par(mfrow=c(1,3))
    plot(ig.mb, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="MB")
    plot(ig.gl, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso")
    plot(ig.sparcc, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="sparcc")

![](https://i.imgur.com/husP6HN.png)

We can evaluate the weights on edges networks using the terms from the
underlying model. SparCC correlations can be used directly, while
SpiecEasi networks need to be massaged a bit. Note that since SPIEC-EASI
is based on penalized estimators, the edge weights are not directly
comparable to SparCC (or Pearson/Spearman correlation coefficients)

    library(Matrix)
    secor  <- cov2cor(getOptCov(se.gl.amgut))
    sebeta <- symBeta(getOptBeta(se.mb.amgut), mode='maxabs')
    elist.gl     <- summary(triu(secor*getRefit(se.gl.amgut), k=1))
    elist.mb     <- summary(sebeta)
    elist.sparcc <- summary(sparcc.graph*sparcc.amgut$Cor)

    hist(elist.sparcc[,3], main='', xlab='edge weights')
    hist(elist.mb[,3], add=TRUE, col='forestgreen')
    hist(elist.gl[,3], add=TRUE, col='red')

![](https://i.imgur.com/xQb3LyF.png)

Lets look at the degree statistics from the networks inferred by each
method.

    dd.gl     <- degree.distribution(ig.gl)
    dd.mb     <- degree.distribution(ig.mb)
    dd.sparcc <- degree.distribution(ig.sparcc)

    plot(0:(length(dd.sparcc)-1), dd.sparcc, ylim=c(0,.35), type='b',
          ylab="Frequency", xlab="Degree", main="Degree Distributions")
    points(0:(length(dd.gl)-1), dd.gl, col="red" , type='b')
    points(0:(length(dd.mb)-1), dd.mb, col="forestgreen", type='b')
    legend("topright", c("MB", "glasso", "sparcc"),
            col=c("forestgreen", "red", "black"), pch=1, lty=1)

![](https://i.imgur.com/mheniLP.png)

Working with phyloseq
---------------------

SpiecEasi includes some convience wrappers to work directly with
`phyloseq` objects.

    library(phyloseq)
    ## Load round 2 of American gut project
    data('amgut2.filt.phy')
    se.mb.amgut2 <- spiec.easi(amgut2.filt.phy, method='mb', lambda.min.ratio=1e-2,
                               nlambda=20, pulsar.params=list(rep.num=50))
    ig2.mb <- adj2igraph(getRefit(se.mb.amgut2),  vertex.attr=list(name=taxa_names(amgut2.filt.phy)))
    plot_network(ig2.mb, amgut2.filt.phy, type='taxa', color="Rank3")

![](https://i.imgur.com/EuIv2Q6.png)

Learning latent variable graphical models
-----------------------------------------

It can be shown that unobserved, latent variables introduce artifacts
into empirical estimates of OTU-OTU associations. These effects can be
removed from the network by treating the inverse covariance selection
problem as a sparse + low-rank decomposition (SPIEC-EASI slr), where the
sparse component are the associations encoded as a conditional
independence graph, and the low-rank components are the isolated latent
effects.

Please see the
[preprint](https://www.biorxiv.org/content/10.1101/2019.12.21.885889v1.full)
and the manuscript [Synapse
project](https://www.synapse.org/#!Synapse:syn20843558) or [github
repository](https://github.com/zdk123/SpiecEasiSLR_manuscript) for more
details.

To demonstrate this in action we can show that removing latent effects
improves a consistency measure between round 1 and round 2 of the
American Gut project data.

First we fit the networks, assuming that there are 10 latent components
in the dataset:

    se1.mb.amgut <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=1e-2,
                              nlambda=20, pulsar.params=list(rep.num=20, ncores=4))
    se2.mb.amgut <- spiec.easi(amgut2.filt.phy, method='mb', lambda.min.ratio=1e-2,
                              nlambda=20, pulsar.params=list(rep.num=20, ncores=4))


    se1.slr.amgut <- spiec.easi(amgut1.filt, method='slr', r=10, lambda.min.ratio=1e-2,
                              nlambda=20, pulsar.params=list(rep.num=20, ncores=4))
    se2.slr.amgut <- spiec.easi(amgut2.filt.phy, method='slr', r=10, lambda.min.ratio=1e-2,
                              nlambda=20, pulsar.params=list(rep.num=20, ncores=4))

Then we compare the consistency between the edge sets within each method
using the Jaccard index.

    otu1 <- colnames(amgut1.filt)
    otu2 <- taxa_names(amgut2.filt.phy)
    edge.diss(getRefit(se1.mb.amgut), getRefit(se2.mb.amgut), 'jaccard', otu1, otu2)
    edge.diss(getRefit(se1.slr.amgut), getRefit(se2.slr.amgut), 'jaccard', otu1, otu2)

Consistency should be a bit better for the slr networks.

Construct the robust PCA from amgut2 data

    X <- se2.slr.amgut$est$data
    L <- se2.slr.amgut$est$resid[[getOptInd(se2.slr.amgut)]]
    pca <- robustPCA(X, L)

We can also check the correlation between AGP meta-data and the latent
factors (scores of the robust PCA).

    age <- as.numeric(as.character(sample_data(amgut2.filt.phy)$AGE))
    bmi <- as.numeric(as.character(sample_data(amgut2.filt.phy)$BMI))
    depth <- colSums(otu_table(amgut2.filt.phy))

    cor(age, pca$scores, use='pairwise')
    cor(bmi, pca$scores, use='pairwise')
    cor(depth, pca$scores, use='pairwise')

Cross domain interactions
-------------------------

SpiecEasi now includes a convenience wrapper for dealing with multiple
taxa sequenced on the same samples, such as 16S and ITS, as seen in
[Tipton, Müller, et.
al. (2018)](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0393-0).
It assumes that each taxa is in it’s own data matrix and that all
samples are in all data matrices in the same order.

Here’s an example run from the [HMP2
project](https://ibdmdb.org/tunnel/public/summary.html) with 16S and
Proteomics data.

    library(phyloseq)
    data(hmp2)
    se.hmp2 <- spiec.easi(list(hmp216S, hmp2prot), method='mb', nlambda=40,
                  lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05))

    dtype <- c(rep(1,ntaxa(hmp216S)), rep(2,ntaxa(hmp2prot)))
    plot(adj2igraph(getRefit(se.hmp2)), vertex.color=dtype+1, vertex.size=9)

![](https://i.imgur.com/86OwfM2.png)

pulsar: parallel utilities for model selection
----------------------------------------------

SpiecEasi is now using the [pulsar
package](https://cran.r-project.org/package=pulsar) as the backend for
performing model selection. In the default parameter setting, this uses
the same [StARS](https://arxiv.org/abs/1006.3316) procedure as previous
versions. As in the previous version of SpiecEasi, we can supply the
`ncores` argument to the pulsar.params list to break up the subsampled
computations into parallel tasks. In this example, we set the random
seed to make consistent comparison across experiments.

    ## Default settings ##
    pargs1 <- list(rep.num=50, seed=10010)
    t1 <- system.time(
    se1 <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=1e-3, nlambda=30,
                  sel.criterion='stars', pulsar.select=TRUE, pulsar.params=pargs1)
    )
    ## Parallel multicore ##
    pargs2 <- list(rep.num=50, seed=10010, ncores=4)
    t2 <- system.time(
    se2 <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=1e-3, nlambda=30,
                  sel.criterion='stars', pulsar.select=TRUE, pulsar.params=pargs2)
    )

We can further speed up StARS using the
[bounded-StARS](https://arxiv.org/abs/1605.07072) (‘bstars’) method. The
B-StARS approach computes network stability across the whole lambda
path, but only for the first 2 subsamples. This is used to build an
initial estimate of the summary statistic, which in turn gives us a
lower/upper bound on the optimal lambda. The remaining subsamples are
used to compute the stability over the restricted path. Since denser
networks are more computational expensive to compute, this can
significantly reduce computational time for datasets with many
variables.

    t3 <- system.time(
    se3 <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=1e-3, nlambda=30,
                   sel.criterion='bstars', pulsar.select=TRUE, pulsar.params=pargs1)
    )
    t4 <- system.time(
    se4 <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=1e-3, nlambda=30,
                   sel.criterion='bstars', pulsar.select=TRUE, pulsar.params=pargs2)
    )

We can see that in addition to the computational savings, the refit
networks are identical.

    ## serial vs parallel
    identical(getRefit(se1), getRefit(se2))
    t1[3] > t2[3]
    ## stars vs bstars
    identical(getRefit(se1), getRefit(se3))
    t1[3] > t3[3]
    identical(getRefit(se2), getRefit(se4))
    t2[3] > t4[3]

### Batch Mode

Pulsar gives us the option of running stability selection in batch mode,
using the [batchtools](https://mllg.github.io/batchtools/) package. This
will be useful to anyone with access to an hpc/distributing computing
system. Each subsample will be independently executed using a
system-specific cluster function.

This requires an external config file which will instruct the batchtools
registry how to construct the cluster function which will execute the
individual jobs. `batch.pulsar` has some built in config files that are
useful for testing purposes (serial mode, “parallel”, “snow”, etc), but
it is advisable to create your own config file and pass in the absolute
path. See the [batchtools
docs](https://mllg.github.io/batchtools/articles/batchtools.html#configuration-file)
for instructions on how to construct config file and template files
(i.e. to interact with a queueing system such as TORQUE or SGE).


    ## bargs <- list(rep.num=50, seed=10010, conffile="path/to/conf.txt")
    bargs <- list(rep.num=50, seed=10010, conffile="parallel")
    ## See the config file stores:
    pulsar::findConfFile('parallel')

    ## uncomment line below to turn off batchtools reporting
    # options(batchtools.verbose=FALSE)
    se5 <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=1e-3, nlambda=30,
                sel.criterion='stars', pulsar.select='batch', pulsar.params=bargs)

Troubleshooting
---------------

A common issue that comes up with when running `spiec.easi` is coming up
with an empty network after running StARS.

For example:

    pargs <- list(seed=10010)
    se <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=5e-1, nlambda=10, pulsar.params=pargs)
    # Warning in pulsar(data = X, fun = match.fun(estFun), fargs = args, seed =
    # 10010, : Optimal lambda may be smaller than the supplied values
    getOptInd(se)
    # [1] 1
    sum(getRefit(se))/2
    # [1] 0

As the warning indicates, the network stability could not be determined
from the lambda path. Looking at the stability along the lambda path,
`se$select$stars$summary`, we can see that the maximum value of the
StARS summary statistic never crosses the default threshold (0.05).

This problem we can fix by lowering `lambda.min.ratio` to explore denser
networks.

    se <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=1e-1, nlambda=10, pulsar.params=pargs)

We have now fit a network, but since we have only a rough, discrete
sampling of networks along the lambda path, we should check how far we
are from the target stability threshold (0.05).

    getStability(se)
    # [1] 0.034032
    sum(getRefit(se))/2
    # [1] 158

To get closer to the mark, we should bump up `nlambda` to more finely
sample of the lambda path, which gives a denser network.

    se <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=1e-1, nlambda=100, pulsar.params=pargs)
    getStability(se)
    # [1] 0.04946882
    sum(getRefit(se))/2
    # [1] 210
