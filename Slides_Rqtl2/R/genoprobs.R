source("colors.R")

pdf(file="../Figs/genoprobs.pdf", width=9.75, height=6.5, pointsize=24, onefile=TRUE)
par(mar=rep(0.1,4),las=1,fg="white",col="white",col.axis="white",col.lab="white",
    bg=bgcolor,bty="n")


plot(0,0,xlim=c(0,100),ylim=c(0,100),xaxt="n",yaxt="n",xlab="",ylab="",type="n")

x <- seq(0.5, 50, length=8)
text(x, 95, LETTERS[1:8])


file <- "_cache/attie.RData"
if(file.exists(file)) {
    load(file)
} else {
    attieDO <- readRDS("~/Projects/AttieDO/CalcGenoProb/attieDO.rds")

    set.seed(14511670)
    chr <- sample(1:19, 1)

    K <- qtl2geno::n_mar(attieDO)[chr]

    k <- 20
    start <- sample(1:(K-k), 1)

    fg <- attieDO$founder_geno[[chr]][,start - 1 + 1:k]
    g <- attieDO$geno[[chr]][,start-1+1:k]

    save(k, fg, g, file=file)
}

y <- seq(2.5, 87.5, length=k)
for(i in 1:8)
    points(rep(x[i], k), y, pch=21, bg=c("white", bgcolor)[(fg[i,]==3)+1])

gg <- rbind(g[1,], g[1,])
gg[1,gg[1,]==2] <- 1
gg[2,gg[2,]==2] <- 3

for(i in 1:k) {
    if(runif(1) < 0.5)
        gg[1:2,i] <- gg[2:1,i]
}

xx <- x[1:2]+70
for(i in 1:2)
    points(rep(xx[i], k), y, pch=21, bg=c("white", bgcolor)[(gg[i,]==3)+1])

text(mean(xx), 95, "DO-101")

dev.off()
