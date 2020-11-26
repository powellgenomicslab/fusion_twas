args = commandArgs(trailingOnly=T)

cluster <- args[1]
chrom <- args[2]

print(sprintf("Cluster %s - Chromosome %s", cluster, chrom))

input_glob <- paste0("/directflow/SCCGGroupShare/projects/annsen/analysis/AMD_scRNA/twas_analysis", "/", "outputs", "/", cluster, "/", "Chr", chrom, "/", "weights", "/", "*.wgt.RDat")
output_file <- paste0("/directflow/SCCGGroupShare/projects/annsen/analysis/AMD_scRNA/twas_analysis", "/", "outputs", "/", cluster, "/", cluster, "_", "Chr", chrom, "_ModelSummary.txt")

lst = Sys.glob(input_glob)
names = gsub(".wgt.RDat","",basename(lst))
N = length(lst)

models = c("top1","blup","enet","bslmm","lasso")
colnames = c( "id" , "nsnps" , "hsq" , "hsq.se" , "hsq.pv" , paste(models,"r2",sep='.') , paste(models,"pv",sep='.') )

mat.snps = matrix(nrow=N,ncol=1)
mat.hsq = matrix(nrow=N,ncol=2)
mat.hsqpv = matrix(nrow=N,ncol=1)
mat.mod.rsq = matrix(nrow=N,ncol=length(models))
mat.mod.pv = matrix(nrow=N,ncol=length(models))

for ( i in 1:N ) {
load(lst[i])
m = match( models , colnames(cv.performance) )
mat.snps[i] = nrow(snps)
mat.hsq[i,] = hsq
mat.hsqpv[i] = hsq.pv
mat.mod.rsq[i,] = cv.performance[1,m]
mat.mod.pv[i,] = cv.performance[2,m]
}

# Build this table properly so we can actually review this
summary_df <- data.frame(id = names, nsnps = mat.snps[,1])
hsq_mat <- cbind(mat.hsq, mat.hsqpv)
hsq_df <- as.data.frame(hsq_mat)
colnames(hsq_df) <- c("hsq", "hsq.se", "hsq.pv")
mod_df <- cbind(mat.mod.rsq, mat.mod.pv)
colnames(mod_df) <- c(paste(models,"r2",sep='.') , paste(models,"pv",sep='.'))
summary_df <- as.data.frame(cbind(summary_df, hsq_df))
summary_df <- as.data.frame(cbind(summary_df, mod_df))
summary_df <- format(summary_df, digits = 3)
write.csv(summary_df, output_file)

#write.table( cbind( names , mat.snps , format(mat.hsq,digits=2) , format(mat.hsqpv,digits=2) , format(round(mat.mod.rsq,3),digits=3) , format(mat.mod.pv,digits=2) ),
#quote=F, col.names=colnames, row.names=F,sep='\t')
best = apply(mat.mod.rsq,1,max,na.rm=T)

options(digits=3)
cat( "Average hsq:" , mean(mat.hsq[,1]) , '(' , sd(mat.hsq[,1])/sqrt(N) , ') \n' , file=stderr() )

cat( "\nAverage crossvalidation R2:\n" , file=stderr() )
write.table( cbind( format(apply(mat.mod.rsq,2,mean,na.rm=T),digits=3) , format(apply(mat.mod.rsq,2,sd,na.rm=T) / sqrt(N),digits=3) ) , quote=F , col.names=c("R2","SE") , row.names=models , sep='\t' , file=stderr() )
cat( "BEST" , format(mean(best),digits=3) , '\n' , sep='\t' , file=stderr() )

cat( "\n% Model is best:\n" , file=stderr() )
for ( i  in 1:length(models) ) {
cat( models[i] , ": " , 100*mean( mat.mod.rsq[,i] == best , na.rm=T ) , "%" , '\n' , sep='' , file=stderr() )
}
