#########################
#   Testing PLINK files using SSD file 
#########################
lbrat.SSD.OneSet_SetIndex = function(SSD.INFO, SetIndex, lbrat_est.obj){

  id1<-which(SSD.INFO$SetInfo$SetIndex == SetIndex)
  if(length(id1) == 0){
    MSG<-sprintf("Error: cannot find set index [%d] from SSD!", SetIndex)
    stop(MSG)
  }
  SetID<-SSD.INFO$SetInfo$SetID[id1]

  try1<-try(Get_Genotypes_SSD(SSD.INFO, SetIndex, is_ID=T),silent = TRUE)
  if(class(try1) != "try-error"){
    G<-try1
    Is.Error<-FALSE
  } else {
    err.msg<-geterrmessage()
    msg<-sprintf("Error to get genotypes of %s: %s",SetID, err.msg)
    stop(msg)
  }
  re<-cat_test(lbrat_est.obj, G)

  return(re)
}

#'LBRAT test or RGMMAT test using SSD format files
#'
#' 
#'
#' @param SSD.INFO SSD format information file, output of function "Open_SSD". The genome wide scan are run set by set.
#' @param lbrat_est.obj ouput from lbrat_est.R
#' @param ... Other options of the LBRAT or RGMMAT test. Deﬁned same as in function "lbrat_test()".
#' 
#' 
#' 
#' @return reults of the LBRAT or RGMMAT test. First column contains batchID, second column contains SNP ID, third column concains prospective P-value and forth column contains retrospective P-value
#' 
#' @export 

lbrat.SSD.All = function(SSD.INFO, lbrat_est.obj, ...){
  N.Set<-SSD.INFO$nSets
  OUT.pvalue<-c()

  for(i in 1:N.Set){
    if(i%%100==0){print(paste0(i," sets finished"))}
    Is.Error<-TRUE
    try1 = try(lbrat.SSD.OneSet_SetIndex(SSD.INFO=SSD.INFO, SetIndex=i, lbrat_est.obj=lbrat_est.obj))

    if(class(try1) != "try-error"){
      re<-try1;
      Is.Error<-FALSE
    } else {

      err.msg<-geterrmessage()
      msg<-sprintf("Error to run GA for %s: %s",SSD.INFO$SetInfo$SetID[i], err.msg)
      warning(msg,call.=FALSE)

    }
    print(paste(i,"-th set",SSD.INFO$SetInfo[i,2],",",SSD.INFO$SetInfo[i,3],"SNPs"))

    if(!Is.Error){
      temp.single<-data.frame(rep(SSD.INFO$SetInfo[i,2],nrow(re)),rownames(re),re$maf, re$pval.pro, re$pval.retro)
      rownames(temp.single)<-NULL
      OUT.pvalue<-rbind(OUT.pvalue,temp.single)

      #OUT.Marker.Test[i]<-re$param$n.marker.test
    }
  }
  colnames(OUT.pvalue) = c('Batch.ID', "SNP.name", "MAF", "Pval.pro", "Pval.retro")
  OUT.pvalue$Batch.ID <- NULL

  return(OUT.pvalue)
}

#' Read in Phenotype and Genotype Files
#'
#' This function loads a phenotype file, a covariates file and a plink fam file
#' it's able to match the phenotype ID with genotype ID and generate input dataframe
#' for L-BRAT/RGMMAT test.
#'
#' @param y.file Path to the input file
#' @param cov.file Path to the covariates file
#' @param plink.file Path to the plink file frefix
#' @return A formatted dataframe used as L-BRAT/RGMMAT test input
#' @export

lbrat_read_phe  <- function(y.file, cov.file, plink.file){

 cat("Checking phenotype file......\n");
 cov <-  try(read.csv(cov.file, header = T, stringsAsFactors=F))
 y  <-  try(read.csv(y.file, header = T, stringsAsFactors=F));

if (class(cov)=="try-error"||class(y)=="try-error")
{
	cat("! Can not open file phenotype files \n");
	return(list(bSuccess=F));
}else{
{	
  cat('From Y file:\n')
	cat("* Individuals =", NROW(y), "\n");
	cat("* Times =", NCOL(y), "\n");
	cat("* Mean =",  colMeans(y, na.rm=T), "\n");
	cat("* SD =",    colSds(y, na.rm=T),"\n");

	
	cat('From Cov file:\n')
	cat("* Individuals =", NROW(cov), "\n");
	
}}

cat("* PLINK.FAM =", fam.file , "\n");

 fam.file = paste(plink.file, '.fam', sep = "");
 fam <- try( read.table(fam.file, header=F, stringsAsFactors=F) );
if (class(fam)=="try-error")
{
	cat("! Can not open file(", fam.file, ")\n");
	return(list(bSuccess=F));
}else{
  cat("* Individuals = ", nrow(fam), "\n")
	cat("  First 5 Items of PLINK .fam:\n");
	print(head(fam, n=5));
}

 #match with PLINK id.
 # g_id <- fam$V2
cat("Start Matching ......\n")
 time  <-  matrix(, ncol = 2)
 matched_id <- Reduce(intersect,list(y[,1],fam$V2,cov[,1]))

 y_match <- y[y[,1]%in%matched_id,];
 cov_match <- cov[cov[,1]%in%matched_id,];
y_order <- y_match[na.omit(match( fam$V2, y_match[,1])),]
cov_order <- cov_match[na.omit(match(fam$V2, cov_match[,1])),]

 if(length(matched_id)<10){
 	cat("! ID MATCH ERROR between phenotype and genotype. \n");
 }else{
	cat("  Matched Individuals=", length(matched_id), "\n");
  cat("Table of # of Obs. ")
  print(table(apply(y_order[,-1],1, function(x) sum(!is.na(x)))))
 }
 

#remove too much missingness
  time_count <- apply(y_order[,-1], 1, function(x)sum(!is.na(x)))
  id.remove <- which(time_count < 1)

  if(length(id.remove)>0){
    cat("! ",length(id.remove), " subjects removed due to missingness.\n")
    y_order<-y_order[-c(id.remove),];
    cov_order<-cov_order[-c(id.remove),];
  cat("  First 5 of Matched Outcome:\n");
  print(head(y_order, n=5));
  cat("  First 5 of Covariates:\n")
  print(head(cov_order, n = 5))
  }

cat("* Number of subjects for analysis: ", nrow(y_order), "\n")



#long format

plink_list <- (fam[fam$V2 %in% y_order[,1],1:2])
write.table(plink_list, file = "ID_to_keep.txt", quote=F, row.names = F, col.names = F)

phe.wide = y_order
cov.wide = cov_order
cat("Start converting ......\n")
 cov_order <- cov_order[,-1]
 y_order <- y_order[,-1]
 y.long  <- c()
 cov.long  <- matrix(, ncol = ncol(cov_order))
colnames(cov.long) <- colnames(cov_order)

 for (i in 1:nrow(y_order)){
  for (j in 1:ncol(y_order)){
    if (!is.na(y_order[i,j]))
    {
      y.long  <-  rbind(y.long, y_order[i,j]);
      time  <-  rbind(time, c(i, j));
      cov.long  <-  rbind(cov.long, cov_order[i,]);
    }
  }
}
 time  <-  time[-1,]
 cov.long  <-  cov.long[-1,]
 cov.long  <-  data.frame(cov.long)
 row.names(cov.long) <- c()
 p0 = list(phe.long=y.long, phe.cov=cov.long, phe.time = time, phe.wide = phe.wide, cov.wide = cov.wide)


cat('Output files -------------------------------------\n')
cat("* Plink keep list save to File:   ID_to_keep.txt.\n")
cat('--------------------------------------------------\n')
cat('Done !\n')
 return(p0)
}

colSds<-function(mat, na.rm=T)
{
	r<-c();
	for(i in 1:dim(mat)[2])
		r <- c(r, sd(mat[,i], na.rm=na.rm));
	return(r);
}


SNPs_maps<-function(bimfile, rsID){
	bim<-fread(bimfile)
	bim.matched<-bim[which(bim$V2 %in% rsID),c(2, 1, 4)]
	colnames(bim.matched)= c('SNP.name', 'chr', 'bp')	
	bim.matched$SNP.name<-as.character(bim.matched$SNP.name)
	bim.matched <- data.frame(bim.matched)
	return(bim.matched)
}



cat_addpos<-function(bimfile, result){

  # n_occur <- data.frame(table(result$SNP.name));
	result$SNP.name<-as.character(result$SNP.name);
	# result$GENE.name<-as.character(result$GENE.name);
        # duplicated_snp<-as.character(n_occur$Var1[n_occur$Freq>1])
        # rna_gene<-c('AC072062.1', 'RP11-325B23.2')
        # result = result[-which(result$GENE.name %in% rna_gene & result$SNP.name %in% duplicated_snp),] 
	rsID<-unique(result$SNP.name);
	pos <-SNPs_maps(bimfile, rsID);
	# result<-result[which(result$Pval.retro<0.05),]
	forplot<-join(pos, result, type = "inner")
        forplot<-forplot[order(forplot$Pval.retro ),]
	return(forplot)
}

gen_batch_SetID  <- function(bimfile, file.setid, batch.size = 1000){
  snp.list  <- fread(bimfile, select = c(2))$V2
  n.snp <- length(snp.list)
  batch.num <-c(rep(1:(n.snp%/%batch.size),each = batch.size), rep(n.snp%/%batch.size + 1, n.snp%%batch.size))
  batch.list  <- paste("batch", batch.num, sep = "")
  batch.setid <- cbind(batch.list, snp.list)
  write.table(batch.setid, file = file.setid, row.names = F, col.names = F, quote = F)
  cat('* Batch SetID file is saved!\n')
}

read.grm  = function(filename){
  mat <- scan(filename, what = numeric())
  ncol <- (sqrt(8 * length(mat) + 1) - 1) / 2
  diag_idx <- cumsum(seq.int(ncol))
  split_idx <- cummax(sequence(seq.int(ncol)))
  split_idx[diag_idx] <- split_idx[diag_idx] - 1
  splitted_rows <- split(mat, f = split_idx)
  mat_full <- suppressWarnings(do.call(rbind, splitted_rows))
  mat_full[upper.tri(mat_full)] <- t(mat_full)[upper.tri(mat_full)]
  return(mat_full)
}


