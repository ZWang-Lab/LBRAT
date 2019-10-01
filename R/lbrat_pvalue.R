#' Calculate prospective and retrospective P-values for GEE or GLMM model
#'
#' This function tests a SNPs for a given SNP set for a given lbrat estimated null model.
#'
#' @param lbrat.est The output of function "lbrat_est()"
#' @param G The genotype matrix, an m*q matrix where m is the number of subjects and q is the total number genetic variants. 
#' @param impute.method choose the iputation method when there is missing genotype. Optional options are: 'random', 'fixed' or 'bestguess'.
#' @param GRM takes m-by-m genetic correlation matrix or kinship matrix.
#' 
#' @return This function returns a dataframe. The row name is the SNP ID, the first column is the prospective score statistics, the second colum is the retrospective score statistics, the third column is the prospective pvalue and the forth column is the restrospective pvalue
#' 
#' @export
lbrat_test <-function(lbrat.est, G, impute.method='fixed', GRM = NULL)
{
    res<-lbrat.est$Y.res; V<-lbrat.est$V; V.inv<-lbrat.est$V.inv;X<-lbrat.est$X;N<-nrow(X)
    m<-lbrat.est$m;time<-lbrat.est$time;mu<-lbrat.est$mu;tau<-lbrat.est$tau;cluster.id<-lbrat.est$cluster.id;
    snp.names<-colnames(G); family = lbrat.est$family;
    if(!is.null(GRM)){
      GRM = cov2cor(GRM)
    }

    if(is.vector(G))
    {
      G[G==9]<-NA;G<-Impute(G,'fixed')
      center.G<-as.matrix(G-mean(G))
      var_g<-t(center.G)%*%(center.G)/(m-1);
      maf<-min(mean(G)/2, 1-mean(G)/2);
      G<-as.matrix(center.G[match(time[,1],cluster.id),1])
    }else
    {
      G[G==9]<-NA
      N_MISS<-sum(is.na(G))
      if(N_MISS>0)
      {
          msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS/nrow(G)/ncol(G))
          warning(msg,call.=F)
          G<-Impute(G,'fixed')
      }
      center.G<-t(t(G)-colMeans(G));
      var_g<-apply(G, 2, var)
      maf<-apply(G, 2, function(x) min(mean(x)/2, 1- mean(x)/2) );
      G<-as.matrix(center.G[match(time[,1],cluster.id),])
    }
    P1<-lbrat.est$P1; P2<-lbrat.est$P2;
    Z<-G-P2%*%(P1%*%G);
    #calculate score;
    n.total<-1;n.rep<-as.numeric(table(time[,1]));tran_res<-rep(0, N);V.pro<-0;V.retro<-0; tran_res_B <-rep(0, m)
    for(i in 1:m)
    {
        ni<-n.rep[i];index<-n.total:(n.total+ni-1);n.total<-n.total + ni
        V.inv.i<-V.inv[[i]];mu.i<-mu[index];G.i<-Z[index,];
        if(family$family == "binomial"){
          if(ni>1){
          delta<-diag(mu.i*(1-mu.i))}else{
          delta<-matrix(mu.i*(1-mu.i))
          }
        }else{
          if(ni>1){
          delta  <- diag(length(mu.i));}else{
          delta <- matrix(1);
          }
          # print('here')
        }
        if(lbrat.est$est.type =="GLMM")
        {
            tran_res[index]<-res[index];
            # tran_res[index]<-V.inv.i%*%res[index];
           if(ni>1){
            V.pro<-V.pro+t(G.i)%*%V.inv.i%*%G.i;}else{
            V.pro<- V.pro + G.i%*%V.inv.i%*%(G.i)
           }
        }else{
            tran_res[index]<-delta%*%V.inv.i%*%res[index];
            CovY <- res[index]%*%t(res[index])
          
            if(ni>1){
            V.pro<-V.pro+t(G.i)%*%delta%*%V.inv.i%*%CovY%*%V.inv.i%*%delta%*%G.i;
            }else{
              V.pro<-V.pro+(G.i)%*%delta%*%V.inv.i%*%CovY%*%V.inv.i%*%delta%*%(G.i);
            }
        }
        
        if(is.null(GRM)){    
          V.retro<-V.retro+sum(tran_res[index])^2*var_g;
        }else{
          tran_res_B[i] = sum(tran_res[index])
        }

    } 
        if(!is.null(GRM)){
          # print(GRM)
          V.retro = c(t(tran_res_B)%*%GRM%*%tran_res_B);
          print(V.retro)
        }

    score = c(t(G)%*%tran_res);
    std.pro<-sqrt(diag(V.pro));std.retro<-sqrt(V.retro)
    score.pro<-score/std.pro;score.retro<-score/std.retro
    pval.pro<-pchisq(score.pro^2,df = 1, lower.tail = F);pval.retro<-pchisq(score.retro^2,df = 1, lower.tail=F)

    result<-cbind(score.pro,score.retro, pval.pro,  pval.retro, maf)
    rownames(result)=snp.names;
    result <- as.data.frame(result)
    return(result)

}

Impute<-function(Z, impute.method){
  p<-dim(Z)[2]
  if(impute.method =="random"){
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-rbinom(length(IDX),2,maf1)
      }
    }
  } else if(impute.method =="fixed"){
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-2 * maf1
      }
    }
  } else if(impute.method =="bestguess") {
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-round(2 * maf1)
      }
    }
  } else {
    stop("Error: Imputation method shoud be \"fixed\", \"random\" or \"bestguess\" ")
  }
  return(as.matrix(Z))
}


