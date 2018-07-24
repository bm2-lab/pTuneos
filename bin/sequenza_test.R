####read data and normalization
Args <- commandArgs()
snp_file=Args[6]
cnv_file=Args[7]
out_dir=Args[8]
sample_id=Args[9]

library('sequenza')
snp <- read.table(snp_file, header = TRUE, sep = "\t")
cnv <- read.table(cnv_file, header = TRUE, sep = "\t")
seqz.data <- VarScan2seqz(varscan.somatic = snp, varscan.copynumber = cnv)
#seqz.data <- VarScan2seqz(varscan.somatic = snp)

write.table(seqz.data, paste(out_dir,sample_id,"_seqz.txt",sep=""), col.names = TRUE, row.names = FALSE, sep = "\t", quote=FALSE)
#data.file <- system.file("data", "example.seqz.txt.gz", package = "sequenza")
#seqz.data <- read.seqz(data.file)
#gc.stats <- gc.sample.stats(data.file)
gc.stats <- gc.norm(x = seqz.data$depth.ratio, gc = seqz.data$GC.percent)
gc.vect <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
seqz.data$adjusted.ratio <- seqz.data$depth.ratio / gc.vect[as.character(seqz.data$GC.percent)]
pdf(file=paste(out_dir,sample_id,"_GC_content.pdf",sep=""))
par(mfrow = c(1,2), cex = 1, las = 1, bty = 'l')
matplot(gc.stats$gc.values, gc.stats$raw, type = 'b', col = 1, pch = c(1, 19, 1), lty = c(2, 1, 2), xlab = 'GC content (%)', ylab = 'Uncorrected depth ratio')
legend('topright', legend = colnames(gc.stats$raw), pch = c(1, 19, 1))
hist2(seqz.data$depth.ratio, seqz.data$adjusted.ratio, breaks = prettyLog, key = vkey, panel.first = abline(0, 1, lty = 2), xlab = 'Uncorrected depth ratio', ylab = 'GC-adjusted depth ratio')
dev.off()

## Extract the information from the seqz file
test <- sequenza.extract(paste(out_dir,sample_id,"_seqz.txt",sep=""), gz = FALSE)

names(test)
###Plot chromosome view with mutations, BAF, depth ratio and segments
pdf(file=paste(out_dir,sample_id,"BAF.pdf",sep=""))
chromosome.view(mut.tab = test$mutations[[1]], baf.windows = test$BAF[[1]],ratio.windows = test$ratio[[1]], min.N.ratio = 1,segments = test$segments[[1]], main = test$chromosomes[1])
dev.off()

###Inference of cellularity and ploidy
CP.example <- sequenza.fit(test)

sequenza.results(sequenza.extract = test, cp.table = CP.example,sample.id = sample_id, out.dir=paste(out_dir,sample_id,sep=""))

cint <- get.ci(CP.example)

#cp.plot(CP.example)
#cp.plot.contours(CP.example, add = TRUE, likThresh = c(0.95))


pdf(file=paste(out_dir,sample_id,"_cellularity.pdf",sep=""))
par(mfrow = c(2,2))
cp.plot(CP.example)
cp.plot.contours(CP.example, add = TRUE)
plot(cint$values.cellularity, ylab = "Cellularity",xlab = "posterior probability", type = "n")
select <- cint$confint.cellularity[1] <= cint$values.cellularity[,2] & cint$values.cellularity[,2] <= cint$confint.cellularity[2]
polygon(y = c(cint$confint.cellularity[1], cint$values.cellularity[select, 2], cint$confint.cellularity[2]),x = c(0, cint$values.cellularity[select, 1], 0), col='red', border=NA)
lines(cint$values.cellularity)
abline(h = cint$max.cellularity, lty = 2, lwd = 0.5)
plot(cint$values.ploidy, xlab = "Ploidy",ylab = "posterior probability", type = "n")
select <- cint$confint.ploidy[1] <= cint$values.ploidy[,1] & cint$values.ploidy[,1] <= cint$confint.ploidy[2]
polygon(x = c(cint$confint.ploidy[1], cint$values.ploidy[select, 1], cint$confint.ploidy[2]),y = c(0, cint$values.ploidy[select, 2], 0), col='red', border=NA)
lines(cint$values.ploidy)
abline(v = cint$max.ploidy, lty = 2, lwd = 0.5)
dev.off()


cellularity <- cint$max.cellularity
write.table(cellularity, paste(out_dir,sample_id,"_cellularity.txt",sep=""), col.names = FALSE, row.names = FALSE, sep = "\t", quote=FALSE)
ploidy <- cint$max.ploidy
avg.depth.ratio <- mean(test$gc$adj[, 2])

mut.tab <- na.exclude(do.call(rbind, test$mutations))

mut.alleles <- mufreq.bayes(mufreq = mut.tab$F, depth.ratio = mut.tab$adjusted.ratio, cellularity = cellularity, ploidy = ploidy, avg.depth.ratio = avg.depth.ratio)

allele_tab <- (cbind(mut.tab[,c("chromosome","position","F","adjusted.ratio", "mutation")], mut.alleles))
write.table(allele_tab, paste(out_dir,sample_id,"_mutation_allele.txt",sep=""), col.names = TRUE, row.names = FALSE, sep = "\t", quote=FALSE)


seg.tab <- na.exclude(do.call(rbind, test$segments))
cn.alleles <- baf.bayes(Bf = seg.tab$Bf, depth.ratio = seg.tab$depth.ratio, cellularity = cellularity, ploidy = ploidy, avg.depth.ratio = avg.depth.ratio)

seg.tab <- cbind(seg.tab, cn.alleles)
write.table(seg.tab, paste(out_dir,sample_id,"_seg_copynumber.txt",sep=""), col.names = TRUE, row.names = FALSE, sep = "\t", quote=FALSE)

#pyclone_input <- sequenza2PyClone(mut.tab,seg.tab,sample_id,norm.cn = 2)
#write.table(pyclone_input, paste(out_dir,sample_id,"_pyclone_input_sequenza.txt",sep=""), col.names = TRUE, row.names = FALSE, sep = "\t", quote=FALSE)