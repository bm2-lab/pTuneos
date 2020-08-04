####read data and normalization
Args <- commandArgs()
seq_file=Args[6]
out_dir=Args[7]
sample_id=Args[8]

library('sequenza')

#seqz.data <- VarScan2seqz(varscan.somatic = snp)
#seqz.data <- read.seqz("small.seqz.gz")

#test <- sequenza.extract(seq_file, chromosome.list=c((1:22),"X","Y"), verbose = FALSE)

test <- sequenza.extract(seq_file, chromosome.list=paste("chr",c((1:22),"X","Y"),sep=""), verbose = FALSE)

CP <- sequenza.fit(test)

sequenza.results(sequenza.extract = test, cp.table = CP, sample.id = sample_id, out.dir=out_dir)

cint <- get.ci(CP)

cellularity <- cint$max.cellularity
write.table(cellularity, paste(out_dir,paste(sample_id,"_cellularity.txt",sep=""),sep='/'), col.names = FALSE, row.names = FALSE, sep = "\t", quote=FALSE)
ploidy <- cint$max.ploidy
write.table(ploidy, paste(out_dir,paste(sample_id,"_ploidy.txt",sep=""),sep='/'), col.names = FALSE, row.names = FALSE, sep = "\t", quote=FALSE)


