<!-- README.md is generated from README.Rmd. Please edit that file -->


SpiecEasi
=========
[![Build Status](https://travis-ci.org/zdk123/SpiecEasi.svg?branch=pulsar)](https://travis-ci.org/zdk123/SpiecEasi)

Sparse InversE Covariance estimation for Ecological Association and Statistical Inference

This package will be useful to anybody who wants to infer graphical models for all sorts of compositional data, though primarily intended for microbiome relative abundance data (generated from 16S amplicon sequence data). It also includes a generator for [overdispersed, zero inflated] multivariate, correlated count data. Please see the paper published in [PLoS Comp Bio](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004226).

One small point on notation: we refer to the method as "SPIEC-EASI" and reserve "SpiecEasi" for this package.

## TOC ##
[Installation](#installation)
[Basic Usage](#basic-usage)
[American Gut Data](#analysis-of-american-gut-data)
[phyloseq](#working-with-phyloseq)
[Cross Domain SPIEC-EASI](#cross-domain-interactions)

## Installation ##

I assume that all auxiliary packages are already installed - for example huge, MASS, etc. If you get an unexpected error, you may need to download and install a missing dependency.

From an interactive R session:


```r
library(devtools)
install_github("zdk123/SpiecEasi", ref='pulsar')
library(SpiecEasi)
```
## Basic Usage ##

Lets simulate some multivariate data under zero-inflated negative binomial model, based on (high depth/count) round 1 of the American gut project, with a sparse network. The basic steps are

1. load the data and normalize counts to to common scale (min depth)
2. fit count margins to the model
3. generate a synthetic network
4. generate some synthetic data
5. clr transformation
6. inverse covariance estimation along a lambda (sparsity) path
7. stability selection using the StARS criterion
8. evaluate performance

Obviously, for real data, skip 1-4.


```r
data(amgut1.filt)
depths <- rowSums(amgut1.filt)
amgut1.filt.n  <- t(apply(amgut1.filt, 1, norm_to_total))
amgut1.filt.cs <- round(amgut1.filt.n * min(depths))

d <- ncol(amgut1.filt.cs)
n <- nrow(amgut1.filt.cs)
e <- d
```
Synthesize the data

```r
set.seed(10010)
graph <- make_graph('cluster', d, e)
Prec  <- graph2prec(graph)
Cor   <- cov2cor(prec2cov(Prec))

X <- synth_comm_from_counts(amgut1.filt.cs, mar=2, distr='zinegbin', Sigma=Cor, n=n)
```

the main SPIEC-EASI pipeline: Data transformation, sparse inverse covariance estimation and model selection

```r
se <- spiec.easi(X, method='mb', lambda.min.ratio=1e-2, nlambda=15)
# Applying data transformations...
# Selecting model with pulsar using stars...
# Fitting final estimate with mb...
# done
```

examine ROC over lambda path and PR over the stars index for the selected graph

```r
huge::huge.roc(se$est$path, graph, verbose=FALSE)
stars.pr(getOptMerge(se), graph, verbose=FALSE)
# stars selected final network under: se.est$refit$stars
```

The above example does not cover all possible options and parameters. For example, other generative network models are available, the lambda.min.ratio (the scaling factor that determines the minimum sparsity/lambda parameter) shown here might not be right for your dataset, and its possible that you'll want more repetitions (number of subsamples) for StARS.


## Analysis of American Gut data ##


Now let's apply SpiecEasi directly to the American Gut data. Don't forget that the normalization is performed internally in the `spiec.easi` function. Also, we should use a larger number of stars repetitions for real data. We can pass in arguments to the inner stars selection function as a list via the parameter `pulsar.params`. If you have more than one processor available, you can also supply a number to `ncores`. Also, let's compare results from the MB and glasso methods as well as SparCC (correlation).


```r
se.mb.amgut <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=1e-2,
                          nlambda=20, pulsar.params=list(rep.num=50))
se.gl.amgut <- spiec.easi(amgut1.filt, method='glasso', lambda.min.ratio=1e-2,
                          nlambda=20, pulsar.params=list(rep.num=50))
sparcc.amgut <- sparcc(amgut1.filt)
## Define arbitrary threshold for SparCC correlation matrix for the graph
sparcc.graph <- abs(sparcc.amgut$Cor) >= 0.3
diag(sparcc.graph) <- 0
sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)
## Create igraph objects
ig.mb     <- adj2igraph(getRefit(se.mb.amgut))
ig.gl     <- adj2igraph(getRefit(se.gl.amgut))
ig.sparcc <- adj2igraph(sparcc.graph)
```

Visualize using igraph plotting:

```r
library(igraph)
## set size of vertex proportional to clr-mean
vsize    <- rowMeans(clr(amgut1.filt, 1))+6
am.coord <- layout.fruchterman.reingold(ig.mb)

par(mfrow=c(1,3))
plot(ig.mb, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="MB")
plot(ig.gl, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso")
plot(ig.sparcc, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="sparcc")
```

![plot of chunk unnamed-chunk-7](https://i.imgur.com/KOjxkzI.png)

We can evaluate the weights on edges networks using the terms from the underlying model. SparCC correlations can be used directly, while SpiecEasi networks need to be massaged a bit. Note that since SPIEC-EASI is based on penalized estimators, the edge weights are not directly comparable to SparCC (or Pearson/Spearman correlation coefficients)


```r
library(Matrix)
secor  <- cov2cor(getOptCov(se.gl.amgut))
sebeta <- symBeta(getOptBeta(se.mb.amgut), mode='maxabs')
elist.gl     <- summary(triu(secor*getRefit(se.gl.amgut), k=1))
# Error in triu(secor * getRefit(se.gl.amgut), k = 1): unused argument (k = 1)
elist.mb     <- summary(sebeta)
elist.sparcc <- summary(sparcc.graph*sparcc.amgut$Cor)

hist(elist.sparcc[,3], main='', xlab='edge weights')
hist(elist.mb[,3], add=TRUE, col='forestgreen')
hist(elist.gl[,3], add=TRUE, col='red')
```

![plot of chunk unnamed-chunk-8](https://i.imgur.com/f3ompVE.png)

Lets look at the degree statistics from the networks inferred by each method.


```r
dd.gl     <- degree.distribution(ig.gl)
dd.mb     <- degree.distribution(ig.mb)
dd.sparcc <- degree.distribution(ig.sparcc)

plot(0:(length(dd.sparcc)-1), dd.sparcc, ylim=c(0,.35), type='b',
      ylab="Frequency", xlab="Degree", main="Degree Distributions")
points(0:(length(dd.gl)-1), dd.gl, col="red" , type='b')
points(0:(length(dd.mb)-1), dd.mb, col="forestgreen", type='b')
legend("topright", c("MB", "glasso", "sparcc"),
        col=c("forestgreen", "red", "black"), pch=1, lty=1)
```

![plot of chunk unnamed-chunk-9](https://i.imgur.com/39lTxIu.png)


## Working with phyloseq ##

SpiecEasi includes some convience wrappers to work directly with `phyloseq` objects.

```r
library(phyloseq)
## Load round 2 of American gut project
data('amgut2.filt.phy')
se.mb.amgut2 <- spiec.easi(amgut2.filt.phy, method='mb', lambda.min.ratio=1e-2,
                           nlambda=20, pulsar.params=list(rep.num=50))
ig2.mb <- adj2igraph(getOptRefit(se.mb.amgut2),  vertex.attr=list(name=taxa_names(amgut2.filt.phy)))
# Error in getOptRefit(se.mb.amgut2): could not find function "getOptRefit"
plot_network(ig2.mb, amgut2.filt.phy, type='taxa', color="Rank3")
# Error in "igraph" %in% class(graph): object 'ig2.mb' not found
```

## Cross domain interactions ##

SpiecEasi now includes a convenience wrapper for dealing with multiple taxa sequenced on the same samples, such as 16S and ITS, as seen in [Tipton, Müller, et. al. (2018)](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0393-0). It assumes that each taxa is in it's own data matrix and that all samples are in all data matrices in the same order.

Here's an example run from the [HMP2 project](https://ibdmdb.org/tunnel/public/summary.html) with 16S and Proteomics data.


```r
library(phyloseq)
data(hmp2)
se.hmp2 <- spiec.easi(list(hmp216S, hmp2prot), method='mb', nlambda=40,
              lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05))

dtype <- c(rep(1,ntaxa(hmp216S)), rep(2,ntaxa(hmp2prot)))
plot(adj2igraph(getRefit(se.hmp2)), vertex.color=dtype+1, vertex.size=9)
```

![plot of chunk unnamed-chunk-11](https://i.imgur.com/KTajLks.png)


## pulsar: parallel utilities for model selection ##
SpiecEasi is now using the [pulsar package](https://cran.r-project.org/package=pulsar) as the backend for performing model selection. In the default parameter setting, this uses the same [StARS](https://arxiv.org/abs/1006.3316) procedure as previous versions.
As in the previous version of SpiecEasi, we can supply the `ncores` argument to the pulsar.params list to break up the subsampled computations into parallel tasks.
In this example, we set the random seed to make consistent comparison across experiments.

```r
## Default settings ##
pargs1 <- list(rep.num=50, seed=10010)
t1 <- system.time(
se1 <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=1e-3, nlambda=30,
              sel.criterion='stars', pulsar.select=TRUE, pulsar.params=pargs1)
)
# Applying data transformations...
# Selecting model with pulsar using stars...
# Fitting final estimate with mb...
# done
## Parallel multicore ##
pargs2 <- list(rep.num=50, seed=10010, ncores=4)
t2 <- system.time(
se2 <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=1e-3, nlambda=30,
              sel.criterion='stars', pulsar.select=TRUE, pulsar.params=pargs2)
)
# Applying data transformations...
# Selecting model with pulsar using stars...
# Fitting final estimate with mb...
# done
```

We can further speed up StARS using the [bounded-StARS](https://arxiv.org/abs/1605.07072) ('bstars') method. The B-StARS approach computes network stability across the whole lambda path, but only for the first 2 subsamples. This is used to build an initial estimate of the summary statistic, which in turn gives us a lower/upper bound on the optimal lambda. The remaining subsamples are used to compute the stability over the restricted path. Since denser networks are more computational expensive to compute, this can significantly reduce computational time for datasets with many variables.

```r
t3 <- system.time(
se3 <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=1e-3, nlambda=30,
               sel.criterion='bstars', pulsar.select=TRUE, pulsar.params=pargs1)
)
# Applying data transformations...
# Selecting model with pulsar using bstars...
# Fitting final estimate with mb...
# done
t4 <- system.time(
se4 <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=1e-3, nlambda=30,
               sel.criterion='bstars', pulsar.select=TRUE, pulsar.params=pargs2)
)
# Applying data transformations...
# Selecting model with pulsar using bstars...
# Fitting final estimate with mb...
# done
```

We can see that in addition to the computational savings, the refit networks are identical.

```r
## serial vs parallel
identical(getRefit(se1), getRefit(se2))
# [1] TRUE
t1[3] > t2[3]
# elapsed 
#    TRUE
## stars vs bstars
identical(getRefit(se1), getRefit(se3))
# [1] TRUE
t1[3] > t3[3]
# elapsed 
#    TRUE
identical(getRefit(se2), getRefit(se4))
# [1] TRUE
t2[3] > t4[3]
# elapsed 
#    TRUE
```

### Batch Mode ###

Pulsar gives us the option of running stability selection in batch mode, using the [batchtools](https://mllg.github.io/batchtools/) package. This will be useful to anyone with access to an hpc/distributing computing system. Each subsample will be independently executed using a system-specific cluster function.

This requires an external config file which will instruct the batchtools registry how to construct the cluster function which will execute the individual jobs. `batch.pulsar` has some built in config files that are useful for testing purposes (serial mode, "parallel", "snow", etc), but it is advisable to create your own config file and pass in the absolute path. See the [batchtools docs](https://mllg.github.io/batchtools/articles/batchtools.html#configuration-file) for instructions on how to construct config file and template files (i.e. to interact with a queueing system such as TORQUE or SGE).

```r

## bargs <- list(rep.num=50, seed=10010, conffile="path/to/conf.txt")
bargs <- list(rep.num=50, seed=10010, conffile="parallel")
## See the config file stores:
pulsar::findConfFile('parallel')
# [1] "/usr/local/lib/R/3.4/site-library/pulsar/config/batchtools.conf.parallel.R"

## uncomment line below to turn off batchtools reporting
# options(batchtools.verbose=FALSE)
se5 <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=1e-3, nlambda=30,
            sel.criterion='stars', pulsar.select='batch', pulsar.params=bargs1)
# Applying data transformations...
# Selecting model with batch.pulsar using stars...
# Sourcing configuration file '/usr/local/lib/R/3.4/site-library/pulsar/config/batchtools.conf.snow.R' ...
# Created registry in '/private/var/folders/8h/m7c712_948df5tnnk9cxhhjr0000gp/T/RtmpFu1BZj/registry1c57f0ad6a5' using cluster functions 'Socket'
# Adding 50 jobs ...
# Submitting 50 jobs in 50 chunks using cluster functions 'Socket' ...
# Fitting final estimate with mb...
# done
```





<!-- ## Extracting Neighborhood ##

Sometimes we are interested in the neighborhood that surrounds specific taxa/nodes, as seen in figure 4 in Tipton, Müller, et. al. The node(s) of interest and their direct neighbors can be extracted for future plotting purposes. For example, if we want to look at the first 5 nodes of the Am Gut Data.

```r

hubind <- which(rowSums(se.hmp2$refit) > 7)
# Error in base::rowSums(x, na.rm = na.rm, dims = dims, ...): 'x' must be an array of at least two dimensions
hmp2hubs <- extract_hood(se.hmp2, hubind)
# Error in extract_hood(se.hmp2, hubind): could not find function "extract_hood"
```
And to plot, using igraph plotting:

```r
plot(adj2igraph(hmp2hubs), color=dtype[hubind])
``` -->
# Error: attempt to use zero-length variable name
```
# test
