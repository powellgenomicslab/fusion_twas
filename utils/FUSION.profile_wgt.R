args = commandArgs(trailingOnly=T)
library(data.table)
setDTthreads(7)

cluster <- args[1]

# Define output filenames
inputs_dir <- paste0("/directflow/SCCGGroupShare/projects/annsen/analysis/AMD_scRNA/twas_analysis/processed/", cluster)
output_dir <- paste0("/directflow/SCCGGroupShare/projects/annsen/analysis/AMD_scRNA/twas_analysis/outputs/", cluster)

weights_dir <- paste0(output_dir, "/", "WEIGHTS/")
pos_output <- paste0(output_dir, "/", cluster, ".pos")
profile_output <- paste0(output_dir, "/", cluster, ".profile")
profile_err <- paste0(output_dir, "/", cluster, ".err")
weights_list_output <- paste0(output_dir, "/", cluster, ".list")

input_glob <- Sys.glob(paste0(weights_dir, "*.wgt.RDat"))
weights_list_df <- data.table(file = paste0("WEIGHTS/", basename(input_glob)))
fwrite(weights_list_df, weights_list_output, row.names = FALSE, col.names = FALSE)

names = gsub(".wgt.RDat","",basename(input_glob))
N = length(names)

models = c("top1","blup","enet","bslmm","lasso")
colnames = c( "id" , "nsnps" , "hsq" , "hsq.se" , "hsq.pv" , paste(models,"r2",sep='.') , paste(models,"pv",sep='.') )

mat.snps = matrix(nrow=N,ncol=1)
mat.hsq = matrix(nrow=N,ncol=2)
mat.hsqpv = matrix(nrow=N,ncol=1)
mat.mod.rsq = matrix(nrow=N,ncol=length(models))
mat.mod.pv = matrix(nrow=N,ncol=length(models))

for ( i in 1:N ) {
load(input_glob[i])
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
fwrite(summary_df, profile_output, sep = "\t")

#write.table( cbind( names , mat.snps , format(mat.hsq,digits=2) , format(mat.hsqpv,digits=2) , format(round(mat.mod.rsq,3),digits=3) , format(mat.mod.pv,digits=2) ),
#quote=F, col.names=colnames, row.names=F,sep='\t')
best = apply(mat.mod.rsq,1,max,na.rm=T)

# Write out summary stats
output_stream <- file(profile_err)
options(digits=3)
cat( "Average hsq:" , mean(mat.hsq[,1]) , '(' , sd(mat.hsq[,1])/sqrt(N) , ') \n' , file = output_stream) 
close(output_stream)

output_stream <- file(profile_err, open = "a")
cat( "\nAverage crossvalidation R2:\n" , file=output_stream )
write.table( cbind( format(apply(mat.mod.rsq,2,mean,na.rm=T),digits=3) , format(apply(mat.mod.rsq,2,sd,na.rm=T) / sqrt(N),digits=3) ) , quote=F , col.names=c("R2","SE") , row.names=models , sep='\t' , file=output_stream )
cat( "BEST" , format(mean(best),digits=3) , '\n' , sep='\t' , file=output_stream )

cat( "\n% Model is best:\n" , file=output_stream )
for ( i  in 1:length(models) ) {
cat( models[i] , ": " , 100*mean( mat.mod.rsq[,i] == best , na.rm=T ) , "%" , '\n' , sep='' , file=output_stream )
}

close(output_stream)

# Merge gene pos with ensembl gene information
gene_loc_list <- paste0("/directflow/SCCGGroupShare/projects/annsen/analysis/AMD_scRNA/twas_analysis/inputs/", cluster, "/gene_loc/AMD_GA_Disease_", cluster, "_Chr", as.character(1:22), "_GeneLoc.tsv")
gene_df_list <- lapply(gene_loc_list, fread)
gene_df <- rbindlist(gene_df_list)
gene_ref <- fread("/directflow/SCCGGroupShare/projects/annsen/analysis/AMD_scRNA/important_files/GeneReference.tsv")
gene_ref <- gene_ref[, c("ensembl_gene_id", "gene_symbol")]
colnames(gene_df) <- c("ensembl_gene_id", "CHR", "P0", "P1")
gene_output_df <- merge(gene_df, gene_ref, by = "ensembl_gene_id")
gene_output_df <- gene_output_df[, c("ensembl_gene_id", "gene_symbol", "CHR", "P0", "P1")]
colnames(gene_output_df) <- c("WGT", "ID", "CHR", "P0", "P1")
gene_output_df <- gene_output_df[match(names, gene_output_df$WGT), ]
gene_output_df$WGT <- paste0("WEIGHT/",gene_output_df$WGT, ".wgt.RDat")
fwrite(gene_output_df, pos_output, sep = "\t")

print(sprintf("%s FUSION TWAS files complete!", cluster))
