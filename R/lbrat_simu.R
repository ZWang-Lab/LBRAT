#' Simulation for LBRAT test
#'
#' This function use pre-defined parameters to make the simulation data for the LBRAT test (including type I and power test)
#'
#' @param n.sample Numeric, sample size, number of individuals
#' @param n.time Numeric, number of measurements for each individual
#' @param par List, the parameters for the phenotype traits, including covaraites and individual specific time dependent random effects
#' @param time_cov Logical variable, indicating whether time effect is included in phenotypic traits
#' @param snp.count Numeric, number of SNPs
#' @param intercept Logical variable, indicating whether intercept is used in phenotypic traits
#' @param disease.para List, the parameters for disease allele and its effect size for power simulation
#' @param onlypower Logical variable, indicating whether include disease SNPs in the generated SNPs
#' @param phe.model String, the phenotype model, two optional values: 'logistic', 'liability'
#' @param oversampling String, the ascertainment scheme, three optional value: 'random', 'baseline', 'sum'
#' 
#' 
#' 
#' @return A list object is returned to be used as object for LBRAT test
#' @export

lbrat_simu<-function(n.sample=1000, n.time=5, par=list(),
    time_cov = TRUE, snp.count = 1000, intercept=TRUE, disease.para = list(), onlypower = FALSE, phe.model = 'logistic', oversampling = "random")
{

    if(missing(par) || length(par)==0 )
    {
        par <- list(b0 = -2.0, b1 = 0.5, b2= 0.5, btime = 0.2, 
            sig.a = 0.8, sig.b = 0.8, sig.e = 0.8, 
            rho=0.7);
        if(phe.model == "logistic"){
            par$b0 = -2.5;
        }
    }

    if(missing(disease.para)||length(disease.para)==0)
    {
        disease.para<- list(model = "dominance", gamma = 0.34,  p1=0.5, p2=0.1, D=0)
    }
    if(!onlypower){
        snp.mat <- simu_snp(n.sample, snp.count);
    }else{
        snp.mat<- matrix(,nrow = n.sample);
    }

    #### Ascertainment sampling #####

    if(oversampling =="random")
    {
      cat('* No ascertainment, random sampling:\n')
        phe <- simu.binary.phe(n.sample*30, n.time, par, intercept, time_cov, disease.para, phe.model);#changed 12/07
        index <- sample(1:(n.sample*30), n.sample)
        ID.select <- paste("id", index, sep = "")
        phe$y <- phe$y[index,];phe$cov <- phe$cov[which(rownames(phe$cov)%in%ID.select),]; phe$causal.snp<-phe$causal.snp[index,];


      cat('Disease Prevalence: ', apply(phe$y, 2, mean), "\n")


    }else if(oversampling =="baseline")
    {
      cat('* Ascertainment based on baseline:\n')
      phe <- simu.binary.phe(n.sample*30, n.time, par, intercept, time_cov, disease.para, phe.model);
      cat('Disease Prevalence: ', apply(phe$y, 2, mean), "\n")
      index.g1 <-which(phe$y[,1] ==1);
      index.g0 <- which(phe$y[,1]==0);
      index.select.0 <- index.g0[sample(n.sample*0.5)];
      index.select.1 <- index.g1[sample(n.sample*0.5)];
      index.select<-c(index.select.0, index.select.1)
      ID.select <- paste("id", index.select, sep = "")

      phe$y <- phe$y[index.select,];phe$cov <- phe$cov[which(rownames(phe$cov)%in%ID.select),]; phe$causal.snp<-phe$causal.snp[index.select,];
     # }
     cat("After Ascertainment: ", apply(phe$y, 2, mean), "\n")

    }else if(oversampling == "sum"){
      cat('* Ascertainment based on sum:\n')
      phe<- simu.binary.phe(n.sample*50, n.time, par, intercept, time_cov, disease.para, phe.model);
      cat('Disease Prevalence: ', apply(phe$y, 2, mean), "\n")
      row.sum<-apply(phe$y, 1, sum)
      index.g2 <- which(row.sum==n.time);
      index.g0 <- which(row.sum==0);
      index.g1 <- which(row.sum<n.time& row.sum>0);

      index.select.0 <- index.g0[sample(n.sample*0.05)];
      index.select.1 <- index.g1[sample(n.sample*0.9)];
      index.select.2 <- index.g2[sample(n.sample*0.05)];
      index.select <- c(index.select.0, index.select.1, index.select.2);
      ID.select <- paste("id", index.select, sep = "")
     phe$y <- phe$y[index.select,];phe$cov <- phe$cov[which(rownames(phe$cov)%in%ID.select),]; phe$causal.snp<-phe$causal.snp[index.select,];
     # }
     cat("After Ascertainment: ", apply(phe$y, 2, mean), "\n")
    }else{
        print('Error! Choose oversampling methods from random, baseline and sum')
    }
    
 

    colnames(phe$y) <- paste("Y", 1:(NCOL(phe$y)), sep="") ;
    rownames(phe$y) <- paste("ID", 1:NROW(phe$y), sep="");

    colnames(phe$cov) <- paste("X", 1:(NCOL(phe$cov)), sep="") ;
    rownames(phe$cov) <- rep(paste("ID", 1:n.sample, sep=""), each = n.time);

    y.long <- c(t(phe$y));
    y.time <- cbind(rep(1:nrow(phe$y), each =n.time), rep(1:n.time, nrow(phe$y)))
    # cov.long <- phe$cov;
    causal.snp <- phe$causal.snp
    snp.mat <- cbind(snp.mat, causal.snp);
    colnames(snp.mat)[c(ncol(snp.mat)-1, ncol(snp.mat))] <- paste("CAUSAL",1:2, sep="")
    # rownames(cov.long) = NULL

    return(list(phe.wide = phe$y, phe.long = y.long, phe.time = y.time, phe.cov.long=phe$cov, snp.mat = snp.mat))
   
}

 # generate random effect
f.simu<-function( sample, n.time, par)
{
    ncol <- n.time;

    AR1 <- array(0,dim=c(ncol,ncol));
    for(i in 1:ncol)
    for(j in 1:ncol)
        AR1[i,j] <- par$rho^abs(i-j);

    sigma.b <- par$sig.b^2*AR1;
     
    r <- rnorm( sample,  0, par$sig.a ) %*% array(1, dim=c(1,ncol))+
         rmvnorm( sample,  rep(0, ncol), sigma.b ) ;
  
    return(r);
}

simu_snp<-function(n.sample, snp.count = 1000)
{
    file.snp.hap1 <- system.file("extdata", "skat-test-1.hap.gz", package="LBRAT");
    snp.hap <- read.table(file.snp.hap1, header=F);
    p.sel<-sample(3:ncol(snp.hap))

    snp.mat1<-snp.hap[ sample(n.sample), p.sel];
    snp.mat2 <- snp.hap[ sample(n.sample), p.sel];
    snp.mat <- snp.mat1 + snp.mat2 -2;
    #check var(g) !=0       
    maf <- colMeans(snp.mat)/2;m.same <- which( maf==1 |  maf==0 );
    if (length(m.same)>0) snp.mat <- snp.mat[, -m.same, drop=F ];
    #check MAF >0.05
    maf <- colMeans(snp.mat)/2; m.small<-which(maf<0.05|maf>0.95);
    if(length(m.small)>0) snp.mat <- snp.mat[, -m.small, drop=F ];
    snp.mat <-snp.mat[,1:snp.count]
    rownames(snp.mat)<-paste("ID", 1:nrow(snp.mat), sep = "");
    colnames(snp.mat)<-paste("SNP", 1:ncol(snp.mat), sep = "")
    return(snp.mat)
}

#generate null phenotype
simu.binary.phe<-function( n.sample, n.time, par, intercept, time_cov, disease.para, phe.model){


    cov.mat <- cbind( rnorm(n.sample*n.time, 0, 1 ), rep(ifelse(runif(n.sample)>0.5, 0, 1), each = n.time));

    if(intercept){
        mu <- f.simu(n.sample, n.time, par) + matrix(cbind(1, cov.mat )%*%c( par$b0, par$b1, par$b2 ), n.sample, n.time, byrow = TRUE)
    }else{
        mu <- f.simu(n.sample, n.time, par) + cov.mat %*% c(par$b1, par$b2 );
    }

    if(time_cov == TRUE){
        time.effect <- rep(1, n.sample)%*%(t(seq(0, n.time-1)*par$btime));
        mu = time.effect+mu;
    }

    #two causal SNPs per simulation:
    model <- disease.para$model;gamma <-disease.para$gamma;
    p1<-disease.para$p1;p2<-disease.para$p2;D<-disease.para$D;


    hap1<-simu_hap(p1, p2, D, n.sample );
    hap2<-simu_hap(p1, p2, D, n.sample );
    causal.snp<-hap1+hap2;

    if(model =="dominance"){
        f_g<-apply(causal.snp, 1, function(x) x[1]>0 && x[2]>0)
    }else if(model =="additive"){
        f_g<-apply(causal.snp, 1, sum)
    }
    disease.effect<-f_g*gamma;

    mu<-mu + disease.effect;

    #prepare output
    mu = c(t(mu));
    if(phe.model=='logistic')
    {
        print('logistic phenotypes')
        y <- matrix(rbinom(n.sample*n.time, 1, inv.logit(mu)),  n.sample,n.time, byrow = T)
    }else if (phe.model=='liability'){#liability model
        print('liability phenotypes')
        #add sigma.e to the mu;
        # print(head(mu))
        mu <- mu + rnorm( n.sample*n.time,  0, par$sig.e);
        # print(head(mu))
        y<- matrix(ifelse(mu>0, 1, 0), n.sample, n.time, byrow =T)
    }else{
        print('Wrong parameter, please choose from logistic and liability\n')
    }
    causal.snp<-causal.snp;
    # cov.mat <- matrix(cov.mat, n.sample, n.time, byrow =T)

    rownames(cov.mat) <- rep(paste("id", 1:n.sample, sep=""), each = n.time);
    return(list(y=y, cov=cov.mat, causal.snp=causal.snp, phe.model = phe.model));
}


inv.logit=function(x){
    return(exp(x)/(1+exp(x)))
}

logit=function(x){
    return(log(x/(1-x)))
}

simu_hap<-function(p1, p2, D, n.sample){
    if(p1*p2<D){print("invalid setting")}
    aa<-p1*p2+D;bb<-(1-p1)*(1-p2)+D;ab<-p1*(1-p2)-D;ba<-(1-p1)*p2-D
    type<-sample(1:4, n.sample, replace = T, prob = c(aa,ab,ba, bb));
    out = matrix(0, n.sample, 2)
    for(i in 1:n.sample){
        if(type[i] ==1){out[i,]<-c(1,1)}
        if(type[i] ==2){out[i,]<-c(1,0)}
        if(type[i] ==3){out[i,]<-c(0,1)}
        if(type[i] ==4){out[i,]<-c(0,0)}
    }
    return(out)
}

