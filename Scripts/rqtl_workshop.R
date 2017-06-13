## R/qtl workshop for CTC in Memphis, TN, 2017-06-13

## 1: load data + summaries

# load R/qtl
library(qtl)

# load example data set from web
sug <- read.cross("csv", "http://rqtl.org", "sug.csv",
                  genotypes=c("CC", "CB", "BB"), alleles=c("C", "B"))

# summary information
summary(sug)

# other summary info
nind(sug)
nmar(sug)
totmar(sug)
nphe(sug)
nchr(sug)

# summary plot
plot(sug)

# individual parts of that plot
plotMap(sug)
plotMissing(sug)
plotPheno(sug, 1)
plotPheno(sug, 2)

## 2: calculate genotype probabilities; QTL genome scan

# calculate QTL genotype probabilities
sug <- calc.genoprob(sug, step=1)

# genome scan of first phenotype by interval mapping
out.em <- scanone(sug)

# plot LOD curves
plot(out.em)

# interactive plot
library(qtlcharts)
iplotScanone(out.em, sug, chr=c(7, 11, 15))

# top LOD score on each chromosome
summary(out.em)

# LOD scores above 4
summary(out.em, threshold=4)

## 3. permutation test
# we'll just do 200 permutations, because it takes a while
operm.em <- scanone(sug, n.perm=200)

# significance thresholds
summary(operm.em)

# alpha = 0.05, 0.2
summary(operm.em, alpha=c(0.05, 0.2))

# histogram of results
plot(operm.em)

# Add threshold to plot
plot(out.em)
add.threshold(out.em, perms=operm.em, lty=2, col="orchid")

## 4. LOD support intervals
# default is to drop 1.5
lodint(out.em, chr=7)

# 2-LOD support interval
lodint(out.em, chr=7, drop=2)

# expand to flanking markers
lodint(out.em, chr=7, drop=2, expandtomarkers=TRUE)

# chr 15
lodint(out.em, chr=15, drop=2, expandtomarkers=TRUE)

# approximate Bayes intervals
bayesint(out.em, chr=7, expandtomarkers=TRUE)
bayesint(out.em, chr=15, expandtomarkers=TRUE)

## 5. Haley-Knott regression 
out.hk <- scanone(sug, method="hk")

# plot the two together
plot(out.em, out.hk, col=c("slateblue", "orchid"), lty=c(1,2))

# another way to make that plot
plot(out.em, col="slateblue")
plot(out.hk, col="orchid", lty=2, add=TRUE)

# plot the differences
plot(out.em - out.hk, ylim=c(-0.5, 0.5), ylab="LOD(EM)-LOD(HK)")

# real advantage is with permutations
operm.hk <- scanone(sug, method="hk", n.perm=1000)

# can also use multiple CPU
operm.hk <- scanone(sug, method="hk", n.perm=1000, n.cluster=8)

## 6. multiple imputation
# first perform the imputations
sug <- sim.geno(sug, step=1, n.draws=256)

# genome scan by imputation
out.imp <- scanone(sug, method="imp")

# plot the three together
plot(out.em, out.hk, out.imp, lty=c(1,2,3))
plot(out.em - out.hk, out.em-out.imp, col=c("slateblue", "orchid"),
     ylim=c(-0.5, 0.5))

## 7. effect plots
# raw phenotype x genotype
(maxpos <- max(out.em))
marker <- find.marker(sug, 7, 47.7)
plot.pxg(sug, marker=marker)

# plot mean phenotype for each genotype
effectplot(sug, mname1=marker)
effectplot(sug, mname1="7@47.7")

# chr 15
max(out.em, chr=15)
effectplot(sug, mname1="15@12")

# and recall that we can use iplotScanone
iplotScanone(out.em, sug, chr=c(7,15))

## 8. genome scan with a different trait + nonparametric
# use pheno.col, as phenotype name
out.hw <- scanone(sug, pheno.col="heart_wt")
# or numeric index
out.hw <- scanone(sug, pheno.col=4)

# plot of lod curves
plot(out.hw)

# heart weight looks perfectly fine to me...
plotPheno(sug, 4)

# but if we wanted to convert to normal quantiles
hw_transf <- nqrank(sug$pheno[,4])

# can give this directly to scanone
out.hw.transf <- scanone(sug, pheno.col=hw_transf)
plot(out.hw.transf, out.hw, col=c("slateblue", "orchid"), lty=1:2)

# non-parametric scan
out.np <- scanone(sug, pheno.col=4, model="np")
plot(out.np, col="green", add=TRUE)

# binary trait: use model="binary"
hw_bin <- as.numeric(sug$pheno[,4] > 120)
out.bin <- scanone(sug, pheno.col=hw_bin, model="binary")
plot(out.bin)

## 9. covariates
cov <- sug$pheno$bw

out.nocov <- scanone(sug, pheno.col="heart_wt")
out.acov <- scanone(sug, pheno.col="heart_wt", addcovar=cov)
plot(out.nocov, out.cov, col=c("slateblue", "orchid"), lty=1:2)

out.icov <- scanone(sug, pheno.col="heart_wt", intcovar=cov)
plot(out.icov - out.acov)

# to do permutations, you need to make sure their properly paired
runif(1, 0, 10^8)
set.seed(22907491)
operm.acov <- scanone(sug, pheno.col="heart_wt", addcovar=cov, 
                      method="hk", n.perm=200)
set.seed(22907491)
operm.icov <- scanone(sug, pheno.col="heart_wt", intcovar=cov,
                      method="hk", n.perm=200)

# plot the permutation results against each other
plot(unclass(operm.acov), unclass(operm.icov))
abline(0,1)

# threshold for the interaction LOD score
summary(operm.icov - operm.acov)

## 10. Split on sex
# this particular cross has all individuals male
# but say we faked it...
sugMF <- sug
sugMF$pheno$sex <- sample(0:1, nind(sug), replace=TRUE)

# now we can split into males and females
sex <- sugMF$pheno$sex
fem <- subset(sugMF, ind=(sex==0))
mal <- subset(sugMF, ind=(sex==1))

# scan for females and males separately
out.fem <- scanone(fem, method="hk")
out.mal <- scanone(mal, method="hk")
plot(out.fem, out.mal, col=c("orange", "purple"))

# test for interaction
out.acov <- scanone(sugMF, addcovar=sex, method="hk")
out.icov <- scanone(sugMF, intcovar=sex, method="hk")
plot(out.icov - out.acov)

# permutations as above, using set.seed() to pair them
set.seed(22907491)
operm.acov <- scanone(sugMF, addcovar=sex,
                      method="hk", n.perm=200)
set.seed(22907491)
operm.icov <- scanone(sugMF, intcovar=sex,
                      method="hk", n.perm=200)

summary(operm.icov - operm.acov)

## 11. X-specific permutation threshold
operm <- scanone(sug, method="hk", perm.Xsp=TRUE, n.perm=1000)
# but these data don't have an X chromosome...
summary(operm)

## 12. Two-dimensional scan
# use a more-coarse step size
sug <- calc.genoprob(sug, step=2.5)

# the 2d scan
out2 <- scantwo(sug, method="hk")

# permutations are really important, but take a very long time
load(url("http://rqtl.org/various.RData"))
summary(operm2)

# summary of scantwo results, with alpha=0.2
summary(out2, perms=operm2, alpha=0.2)

# just chr 7, 12, 15
plot(out2, lower="fv1", upper="int")

# interactive plot
iplotScantwo(out2, sug)

## 13. include marker as covariate
plot(out.hk)
max(out.hk)

# go back to more-fine step size
sug <- calc.genoprob(sug, step=1)

# form covariate matrix for the chr 7 QTL
X <- formMarkerCovar(sug, "7@47.7")

# for an intercross, this is a matrix with 2 columns
head(X)

# scan conditional on the chr 7 QTL
out.c7 <- scanone(sug, method="hk", addcovar=X)
plot(out.hk, out.c7, col=c("slateblue", "orchid"))

## 14. multiple-QTL models
# make a QTL object just with the chr 7 QTL
qtl <- makeqtl(sug, 7, 47.7, what="prob")

# fit the single-QTL model
out.fq <- fitqtl(sug, qtl=qtl, method="hk")
summary(out.fq)

# summary without p-values
summary(out.fq, pvalues=FALSE)

# scan for an additional QTL
out.aq <- addqtl(sug, qtl=qtl, method="hk")
plot(out.aq)

# add to the QTL object
max(out.aq)
qtl <- addtoqtl(sug, qtl=qtl, chr=15, pos=13)

# fit qtl model
out.fq2 <- fitqtl(sug, qtl=qtl, method="hk")
summary(out.fq2, pvalues=FALSE)

# refine the QTL model
rqtl <- refineqtl(sug, qtl=qtl, method="hk")

# plot QTL on the genetic map
plot(rqtl)

# plot LOD profiles
plotLodProfile(rqtl)

# add the single-QTL scan
plot(out.hk, chr=c(7,15), add=TRUE, col="orchid", lty=2)

## 15. stepwise QTL analysis
summary(operm)
out.sq <- stepwiseqtl(sug, method="hk", max.qtl=5, 
                      additive.only=TRUE, 
                      penalties=c(3.45, Inf, Inf))
out.sq

(pen <- calc.penalties(operm2))
out.sqi <- stepwiseqtl(sug, method="hk", max.qtl=5,
                       penalties=pen)
out.sqi
