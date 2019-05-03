#' GEE NULL model estimation
#'
#' This function estimate the parameters and residuals for the NULL model in LBRAT
#'
#' @param y.long Long-formatted phenotype vector 
#' @param time Time covarites matched with phenotype vector
#' @param y.cov Covariate matrix denoting the covariate variables measured at each time
#' @param timecov Logical variable, indicating whether the time fixed effect is estimated 
#' @param corstr String, correlation structure for GEE model, optional values are: 'ar1', 'ind', 'mixture'
#' @param tol Numeric, tolerance for the iterative estimation in when using mixture correlation structure
#' @param max.iter Numeric, the maximum count for the iterative estimation in when using mixture correlation structure
#' 
#' @return This function returns a list object with model parameters and residuals of the NULL GEE model 
#' @export
lbrat_est.gee<-function(y.long, time, y.cov, timecov = TRUE, corstr = "ar1",tol = 10^-6,max.iter = 50)
{
    # just for test
    # y.long = p0$phe.long;
    # time = p0$phe.time;
    # y.cov = p0$phe.cov
    # tol = 10^-6; max.iter = 50
    # corstr = "ar1"
    # timecov= TRUE

    if(length(table(y.long))>2){
      family = gaussian(); scale.fix = F
      cat('Phenotype is continuous, fitting identity link......\n')
    }else{
      family = binomial(); scale.fix = T
      cat('Phenotype is dichotomous, fitting logistic link ......\n')

    }

    Y<-as.matrix(y.long);
    if(timecov == TRUE){
      X=as.matrix(cbind(y.cov, time[,2]-1));
      colnames(X)[ncol(X)]="Time";
    }else 
    {
      X <- as.matrix(y.cov);
    }

    N<-length(y.long)
    X_1 = cbind(rep(1, nrow(X)),X)
    cluster.id<-unique(time[,1]);m<-length(cluster.id)

    if(corstr == "mixture")
    {
      if(family$family == "gaussian"){cat('Error: Mixture Correlation is only available for Binary Outcome\n'); next}
      nullgee0<-geeglm(Y~ X, family=family, id=time[,1], corstr ="unstructured")
          mu.est<-nullgee0$fitted.values;beta.new<-nullgee0$coefficients
      diff <-10;iter <-0
      while(diff > tol && iter<max.iter)
          {
              var.est<-function(x)
              {
                out <- c(0,0);out.mean <- c(0,0);
                n.total<-1;n.rep<-as.numeric(table(time[,1]))
                for(i in 1:m)
                {
                    ni<-n.rep[i];index<-n.total:(n.total+ni-1);n.total<-n.total + ni;
                    y.i = Y[index];mu.i<-mu.est[index];
                    v1<-0.7^abs(outer(time[index,2], time[index,2],"-"));
                    v2<-matrix(1,length(index), length(index));
                    deltasq.inv<-diag((mu.i*(1-mu.i))^(-0.5))
                    sigma.inv<-solve(x[1]*v1+x[2]*v2+(1-x[1]-x[2])*diag(length(index)))
                    # sigma.inv<-solve(x[1]*v2+(1-x[1]-x[2])*v1+x[2]*diag(length(index)))
                    out[1]<-out[1]+t(y.i-mu.i)%*%sigma.inv%*%deltasq.inv%*%sigma.inv%*%(v1-diag(length(index)))%*%deltasq.inv%*%(y.i-mu.i)
                    out[2]<-out[2]+t(y.i-mu.i)%*%sigma.inv%*%deltasq.inv%*%sigma.inv%*%(v2-diag(length(index)))%*%deltasq.inv%*%(y.i-mu.i)
                    out.mean[1]<-out.mean[1]+sum(diag(sigma.inv%*%v1-diag(length(index))));
                     out.mean[2]<-out.mean[2]+sum(diag(sigma.inv%*%v2-diag(length(index))));
                  }
                  out<-out-out.mean;
                  return(out)
              }  

              tau<-nleqslv(c(0.4, 0.4), var.est, jacobian=TRUE)$x
              tau[which(tau<0)]=0;tau[which(tau>1)] =1-tau[which(tau<1)]; if(sum(tau)==0) tau = rep(0.3, 2)

              beta_est<-function(x)
              {   
                  mu<-inv.logit(X_1%*%x);out<-0;
                  n.total<-1;n.rep<-as.numeric(table(time[,1]));phi<-0;
                  for(i in 1:m)
                  {
                      ni<-n.rep[i];index<-n.total:(n.total+ni-1);n.total<-n.total + ni;
                      y.i = Y[index];mu.i<-mu[index];x.i<-X_1[index,];
                      v1<-0.7^abs(outer(time[index,2], time[index,2],"-"));
                      v2<-matrix(1,length(index), length(index));
                      sigma.inv<-solve(tau[1]*v1+tau[2]*v2+(1-tau[1]-tau[2])*diag(length(index)))
                      # sigma.inv<-solve(tau[1]*v2+(1-tau[1]-tau[2])*v1+tau[2]*diag(length(index)))
                      deltasq<-diag((mu.i*(1-mu.i))^(0.5))
                      out<-out+t(x.i)%*%deltasq%*%sigma.inv%*%solve(deltasq)%*%(y.i-mu.i);
                      phi<-phi+t(y.i-mu.i)%*%solve(deltasq)%*%sigma.inv%*%solve(deltasq)%*%(y.i-mu.i);
                  }
                  # print(phi/N)
                  return(out)
              }
              beta.old<-beta.new
              beta.new<-nleqslv(beta.old, beta_est, jacobian=TRUE)$x
              diff = sum((beta.old-beta.new)^2);iter = iter +1;
              # print(diff)
              head(X_1)
              mu.est<-inv.logit(X_1%*%beta.new)
          }
     
    }else if(corstr == "ar1")
    {
      nullgee0<-geeglm(Y~X, family=family, id=time[,1], corstr ="ar1", scale.fix = scale.fix)
          mu.est<-nullgee0$fitted.values;beta.new<-nullgee0$coefficients;
          rho <- as.numeric(summary(nullgee0)$corr[1]);disper = as.numeric(summary(nullgee0)$dispersion[1])
          tau = list(rho = rho, disper = disper);

    }else if(corstr == "ind")
    {
      nullgee0<-geeglm(Y~ X, family=family, id=time[,1], corstr ="independence", scale.fix = scale.fix)
          mu.est<-nullgee0$fitted.values;beta.new<-nullgee0$coefficients;
          disper = as.numeric(summary(nullgee0)$dispersion[1])
          tau = list(disper = disper)
    }
     Y.res<-Y-mu.est;
    n.total<-1;n.rep<-as.numeric(table(time[,1]))
    V<-list();V.inv<-list();P1<-matrix(0, ncol(X_1), N); 
      for (i in 1:m)
        {
          ni<-n.rep[i];index<-n.total:(n.total+ni-1);n.total<-n.total + ni
          mu.i<-mu.est[index];y.res.i <- Y.res[index];covy.i<-y.res.i%*%t(y.res.i)
          if(corstr=="mixture"){
            v1<-0.7^abs(outer(time[index,2], time[index,2],"-"));
            v2<-matrix(1,length(index), length(index));
            sigma<-tau[1]*v1+tau[2]*v2+(1-tau[1]-tau[2])*diag(length(index));
          }else if(corstr =="ar1"){
            sigma<-rho^abs(outer(time[index,2], time[index,2],"-"));
          }else if (corstr =="ind"){
            sigma<-diag(length(index));
          }else{
            cat("Error: incorrect correlation structure!")
          }
          if(family$family=="binomial")
            { 
              if(length(index)>1){delta <- diag((mu.i*(1-mu.i))) }else{delta <- mu.i*(1-mu.i)}
            }else{
              if(length(index)>1){delta <- diag(length(mu.i))}else{delta = 1}
            }
           
          Vi <- sqrt(delta)%*%sigma%*%sqrt(delta);
  
          Vi.inv<-solve(Vi); 
          V[[i]]<-Vi; V.inv[[i]]<-Vi.inv;
          if(length(index)>1){

            # covy.i = Vi # added 12/07/2018
            P1[,index]<-t(X_1[index,])%*%delta%*%Vi.inv%*%covy.i%*%Vi.inv%*%delta;              }else{
              P1[,index]<-t(X_1[index,])*c(covy.i);
          }
        }
        P2<-X_1%*%solve(P1%*%X_1)

    return(list(Y=Y,time=time,X=X_1,cluster.id=cluster.id,m=m,coef = beta.new, tau = tau, est.type = "GEE", P1=P1, P2=P2, V = V, V.inv = V.inv, Y.res=Y.res, mu = mu.est, family = family))
}

#' GLMM NULL model estimation
#'
#' This function estimate the parameters and residuals for the NULL model in RGMMAT test
#'
#' @param y.long Long-formatted phenotype vector 
#' @param time Time covarites matched with phenotype vector
#' @param y.cov Covariate matrix denoting the covariate variables measured at each time
#' @param timecov Logical variable, indicating whether the time fixed effect is estimated
#' 
#' @return This function returns a list object with model parameters and residuals of the NULL GLMM model 
#' @export

lbrat_est.glmm<-function(y.long, time, y.cov, timecov = TRUE)
{
    # just for test
    # y.long = p0$phe.long;
    # y.wide <- p0$phe.wide;
    # time = p0$phe.time;
    # y.cov = p0$phe.cov.long;
    #  tol = 10^-6; max.iter = 50
    #timecov = T
    
    
    Y<-as.matrix(y.long);
    if(timecov == TRUE){
      X=as.matrix(cbind(y.cov, time[,2]-1));
      colnames(X)[ncol(X)]="Time";
    }else 
    {
      X <- as.matrix(y.cov);
    }

    N<- length(y.long);X_1 = cbind(rep(1, nrow(X)),X)
    cluster.id<-unique(time[,1]);m<-length(cluster.id);
    subj = time[,1]; rep = time[,2];
    nullglmm0<-glmer(Y ~ X+(1|subj)+(1|rep), family = binomial);
    beta<-nullglmm0@beta;var.df = as.data.frame(VarCorr(nullglmm0));mu<-nullglmm0@resp$mu;
    tau<-c(var.df$sdcor[2]^2,var.df$sdcor[1]^2)
    # Y.til.res<- (mu*(1-mu))^(-1)*(Y - mu)+logit(mu)-X_1%*%beta;
    Y.res <- (Y-mu);
    #V, V.inv, P1, P2
    n.total<-1;n.rep<-as.numeric(table(time[,1]))
    V<-list();V.inv<-list();P1<-matrix(0, ncol(X_1), N); 
      for (i in 1:m)
        {
          ni<-n.rep[i];index<-n.total:(n.total+ni-1);n.total<-n.total + ni
          mu.i<-mu[index];
          v1<-0.7^abs(outer(time[index,2], time[index,2],"-"));
          v2<-matrix(1,length(index), length(index));
          v3<-diag((mu.i*(1-mu.i))^-1)
          if(ni >1){
            sigma<-tau[1]*v1+tau[2]*v2+v3;
          }else{
            sigma<-(mu.i*(1-mu.i))^-1
            }
          
          Vi<-sigma; Vi.inv<-solve(Vi);  
          V[[i]]<-Vi;V.inv[[i]]<-Vi.inv;
          if(ni>1){
          P1[,index]<-t(X_1[index,])%*%Vi.inv;
          }else{
            P1[,index] <- t(X_1[index,])*c(Vi.inv)
          }
        }
        P2<-X_1%*%solve(P1%*%X_1)

    return(list(Y=Y,time=time,X=X_1,cluster.id=cluster.id,m=m,coef = beta, tau = tau, est.type = "GLMM", P1=P1, P2=P2, V = V, V.inv = V.inv, Y.res=Y.res, mu = mu,family = binomial()))

    #calculate V.inv, P

}



   
