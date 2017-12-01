#put relevant functions in the model
library(ape)
library(optimx)
library(numDeriv)
library(corpcor)
library(MASS)
library(mvtnorm)
library(matrixcalc)
library(Matrix)

#BM model
bm<-function(x.tre,response=response, predictor=predictor){
  NegLogLikeBM<-function(n=n,G=G,response=response,predictor=predictor){
    solveG<-solve(G)
    one<-array(1,c(length(response),1))
    dsm<-cbind(one, matrix(predictor,ncol=1))
    b_est<-pseudoinverse(t(dsm)%*%solveG%*%dsm)%*%t(dsm)%*%solveG%*%response
    sigma_sq_hat<-t(response-dsm%*%b_est)%*%solveG%*%(response-dsm%*%b_est)/n 
    negloglike<-n/2*log(pi)+n/2*log(sigma_sq_hat)+(1/2)*log(det(G))
    negloglike<-negloglike+(1/(2*sigma_sq_hat))*t(response-dsm%*%b_est)%*%solveG%*%(response-dsm%*%b_est)
    return(negloglike)
  }
  
  #tree in newick format with branch length and topology
  G<-vcv(x.tre)
  G<-G/max(G)
  distG<-2*(1-G)
  n<-dim(G)[1]
  one<-array(1,c(n,1))
  solveG<-solve(G)
  
  
  BMmean_predictor<-t(one)%*%solveG%*%predictor
  BMmean_predictor<-BMmean_predictor/(t(one)%*%solveG%*%one)
  BMvcv_predictor<-t(predictor-BMmean_predictor)%*%solveG%*%(predictor-BMmean_predictor)/n
  BMvcv_predictor<-c(BMvcv_predictor)*G
  
  BMmean_response<-t(one)%*%solveG%*%response
  BMmean_response<-BMmean_response/(t(one)%*%solveG%*%one)
  BMvcv_response<-t(response-BMmean_response)%*%solveG%*%(response-BMmean_response)/n
  BMvcv_response<-c(BMvcv_response)*G
  
  dsm<-cbind(one,predictor)#design matrix for initial estimate of b
  b_est<- solve(t(dsm)%*%solveG%*%dsm)%*%t(dsm)%*%solveG%*%response
  sigma_sq_hat<- t(response-dsm%*%b_est)%*%solveG%*%(response-dsm%*%b_est)/n
  est_V<-c(sigma_sq_hat)*G 
  BMvalue<-NegLogLikeBM(n=n,G=G,response=response,predictor=predictor)
  r_sq_value<-r_sq(b_est=b_est,est_V=est_V,dsm=dsm,predictor=predictor,response=response)
  AICc_value<-AICc(n=n,k=3,BMvalue)
  result<-list(BMvalue=BMvalue,sigma_sq=sigma_sq_hat,b_est=b_est,est_V=est_V,r_sq_value=r_sq_value,AICc_value=AICc_value,BMmean_predictor=BMmean_predictor, BMvcv_predictor=BMvcv_predictor,BMmean_response=BMmean_response,BMvcv_response=BMvcv_response)
  return(result)
}#bm

#-----------------------------------------------------------
#OU model
OUvcv<-function(x,G=G,distG=distG){        #vcv for OU, not include sigma_sq
  alp<-x[1]
  A<-(1-exp(-2*alp*G))/(2*alp)
  A<-A*exp(-alp*distG)
  return(A)
}

ou<-function(x.tre,response=response, predictor=predictor){
  
  NegLogLike<-function(x,response=response,predictor=predictor,n=n,one=one,G=G,distG=distG,dsm=dsm){
    alp<-x[1]
    V<-OUvcv(alp,G=G,distG=distG)
    InvOUvcv<-pseudoinverse(V)
    b_est<-pseudoinverse(t(dsm)%*%InvOUvcv%*%dsm)%*%t(dsm)%*%InvOUvcv%*%response
    sigma_sq_hat<-t(response-dsm%*%b_est)%*%InvOUvcv%*%(response-dsm%*%b_est)/n 
    negloglike<-n/2*log(pi)+n/2*log(sigma_sq_hat)+(1/2)*log(det(V))
    negloglike<-negloglike+(1/(2*sigma_sq_hat))*t(response-dsm%*%b_est)%*%InvOUvcv%*%(response-dsm%*%b_est)
    return(negloglike)
    }
  
  
  UniMLEs<-function(x,trait=trait,n=n,one=one,G=G,distG=distG){
    alp<-x[1]
    V<-OUvcv(alp,G=G,distG=distG)
    InvOUvcv<-pseudoinverse(V)
    mu_hat<-t(one)%*%InvOUvcv%*%trait
    mu_hat<-c(mu_hat/(t(one)%*%InvOUvcv%*%one))
    sigma_sq_hat<-t(trait-mu_hat*one)%*%InvOUvcv%*%(trait-mu_hat*one)/n 
    result<-list(V=V,InvOUvcv=InvOUvcv,mu_hat=mu_hat,sigma_sq_hat=sigma_sq_hat)
    }
  
  UniNegLogLike<-function(x,trait=trait,n=n,one=one,G=G,distG=distG){
    alp<-x[1]
    V<-OUvcv(alp,G=G,distG=distG)
    InvOUvcv<-pseudoinverse(V)
    mu_hat<-t(one)%*%InvOUvcv%*%trait
    mu_hat<-c(mu_hat/(t(one)%*%InvOUvcv%*%one))
    sigma_sq_hat<-t(trait-mu_hat*one)%*%InvOUvcv%*%(trait-mu_hat*one)/n 
    negloglike<-n/2*log(pi)+n/2*log(sigma_sq_hat)+(1/2)*log(det(V))
    negloglike<-negloglike+(1/(2*sigma_sq_hat))*t(trait-mu_hat*one)%*%InvOUvcv%*%(trait-mu_hat*one)
    return(negloglike)
  }
  #tree in newick format with branch length and topology
  G<-vcv(x.tre)
  G<-G/max(G)
  distG<-2*(1-G)
  n<-dim(G)[1]
  one<-array(1,c(n,1))
  
  MLEs_predictor<-optimize(UniNegLogLike,c(10^(-10),10),trait=predictor,n=n,one=one,G=G,distG=distG)
  output_predictor<-UniMLEs(MLEs_predictor$minimum,trait=predictor,n=n,one=one,G=G,distG=distG)
  OUmean_predictor<-output_predictor$mu_hat
  OUvcv_predictor<-c(output_predictor$sigma_sq_hat)*output_predictor$V
  
  MLEs_response<-optimize(UniNegLogLike,c(10^(-10),10),trait=response,n=n,one=one,G=G,distG=distG)
  output_response<-UniMLEs(MLEs_response$minimum,trait=response,n=n,one=one,G=G,distG=distG)
  OUmean_response<-output_response$mu_hat
  OUvcv_response<-c(output_response$sigma_sq_hat)*output_response$V
  
  dsm<-cbind(one,matrix(predictor,ncol=1))
  
  old_value<-NegLogLike(1,response=response,predictor=predictor,n=n,one=one,G=G,distG=distG,dsm=dsm)
  new_value<-old_value-1
  attempt<-0
  while(old_value>new_value){
    attempt<-attempt+1
    #print(paste("this is attempt: ", attempt))
    old_value<-new_value
    MLEs<-optimize(NegLogLike,c(10^(-10),100),response=response,predictor=predictor,n=n,one=one,G=G,distG=distG,dsm=dsm) 
    new_value<-MLEs$objective
    #print("old value, new value")
    #print(c(old_value, new_value))
  }
  
  OUvcv_hat<-OUvcv(MLEs$minimum,G=G,distG=distG)
  InvOUvcv_hat<-pseudoinverse(OUvcv_hat)
  b_est<- solve(t(dsm)%*%InvOUvcv_hat%*%dsm)%*%t(dsm)%*%InvOUvcv_hat%*%response
  sigma_sq_hat<- t(response-dsm%*%b_est)%*%InvOUvcv_hat%*%(response-dsm%*%b_est)/n
  est_V<-c(sigma_sq_hat)*OUvcv(MLEs$minimum,G=G,distG=distG)   
  r_sq_value<-r_sq(b_est=b_est,est_V=est_V,dsm=dsm,predictor=predictor,response=response)
  AICc_value<-AICc(n=n,k=4,MLEs$objective)
  result<-list(MLEs=MLEs,sigma_sq=sigma_sq_hat,b_est=b_est,est_V=est_V,r_sq_value=r_sq_value,AICc_value=AICc_value,OUmean_predictor=OUmean_predictor, OUvcv_predictor=OUvcv_predictor,OUmean_response=OUmean_response,OUvcv_response=OUvcv_response)
  return(result)
}#ou

#get sig for optimal 
compute.optimal.sigma<-function(pred.params,predictor=predictor, V=V){
  r<-dim(predictor)[2]
  n<-dim(V)[1]
  one<-array(1,c(n,1))
  inv.V<-pseudoinverse(V)
  sig1<-0
  for(bIndex in 1:r){
    temp.m.pred<-c(t(one)%*%inv.V%*%predictor[,bIndex])
    temp.m.pred<-c(temp.m.pred/(t(one)%*%inv.V%*%one))
    temp.v.pred <- t(predictor[,bIndex] - temp.m.pred*one)%*%inv.V%*%(predictor[,bIndex] - temp.m.pred*one)/n
    sig1<-sig1+pred.params[bIndex]^2*temp.v.pred
  }
  sig1<-sqrt(sig1)
  return(sig1)
}


#-----------------------------------------------------------
#OUBM model
#Input Tree, Trait:response, precitor
#output MLEs, Likelihood, regression coefficient


  
phyloCorrectFactor.OUBM<-function(model.params, pred.params=pred.params, predictor=predictor,x.tre=x.tre){
  alp<-model.params[1]
  sig<-model.params[2]
  #b1<-x[3]  
  
  #alp
  #sig
  #alp<-0.5;sig<-1;b1<-2;b2<-3;b3<-1
  #free.parameters<-c(alp,sig,b1,b2,b3)
  #free.parameters
  #r<- dim(predictor)[2] # b1 to br
  #r
  #pred.params<-x[(2+1):(2+r)] # need to adjust 2 for other models
  #pred.params
  
  rho_y_th<-0
  y_0<-0
  th_0<-0
  
  tree.vcv<-vcv(x.tre)
  sig1<-compute.optimal.sigma(pred.params,predictor=predictor, V=tree.vcv)
  trh<-max(tree.vcv) #tree height
  cov_y_th<--(th_0*(exp(alp*trh) - 1) + y_0)*th_0*exp(-alp*trh) - (sig1*rho_y_th*sig - (alp*exp(alp*trh) - alp)*th_0^2 - alp*th_0*y_0 - sig1^2 - (alp*sig1^2*trh + sig1*rho_y_th*sig - sig1^2)*exp(alp*trh))*exp(-alp*trh)/alp
  var_th<-sig1^2*trh 
  return(c(cov_y_th/var_th))
}




oubm<-function(x.tre,response=response, predictor=predictor){
  CovRes<-function(x,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor){
    alp<-x[1]
    sig<-x[2]
    #b1<-x[3]  
    
    rho_y_th<-0
    y_0<-0
    th_0<-0
    
    
    
    sig1<-compute.optimal.sigma(pred.params,predictor=predictor, V=G)
    
    cov_yi_yj<-array(0,c(n,n))
    cov_y_th<-array(0,c(n,n))
    cov_thi_thj<-array(0,c(n,n))
    b1_hat<-array(0,c(n,n))
    cov_ri_rj<-array(0,c(n,n))  
    
    for(i in 1:n){
      for(j in 1:n){
        
        cov_yi_yj[i,j]<-
          sig1^2*G[i,j]*(exp(1/2*alp*distG[i,j]) - 1)^2*exp(-alp*distG[i,j]) - 2*((th_0*(exp(alp*G[i,j]) - 1) + y_0)*th_0*exp(-alp*G[i,j]) + (sig1*rho_y_th*sig - (alp*exp(alp*G[i,j]) - alp)*th_0^2 - alp*th_0*y_0 - sig1^2 - (alp*sig1^2*G[i,j] + sig1*rho_y_th*sig - sig1^2)*exp(alp*G[i,j]))*exp(-alp*G[i,j])/alp)*(exp(1/2*alp*distG[i,j]) - 1)*exp(-alp*distG[i,j]) - 1/2*(2*(th_0*(exp(alp*G[i,j]) - 1) + y_0)^2*exp(-2*alp*G[i,j]) - (2*sig1*rho_y_th*sig + 2*(alp*exp(2*alp*G[i,j]) - 2*alp*exp(alp*G[i,j]) + alp)*th_0^2 + 4*(alp*exp(alp*G[i,j]) - alp)*th_0*y_0 + 2*alp*y_0^2 - sig1^2 - sig^2 + (2*alp*sig1^2*G[i,j] + 2*sig1*rho_y_th*sig - 3*sig1^2 + sig^2)*exp(2*alp*G[i,j]) - 4*(sig1*rho_y_th*sig - sig1^2)*exp(alp*G[i,j]))*exp(-2*alp*G[i,j])/alp)*exp(-alp*distG[i,j])
        
        cov_y_th[i,j]<-
          -(th_0*(exp(alp*G[i,j]) - 1) + y_0)*th_0*exp(-alp*G[i,j]) - (sig1*rho_y_th*sig - (alp*exp(alp*G[i,j]) - alp)*th_0^2 - alp*th_0*y_0 - sig1^2 - (alp*sig1^2*G[i,j] + sig1*rho_y_th*sig - sig1^2)*exp(alp*G[i,j]))*exp(-alp*G[i,j])/alp
        
        cov_thi_thj[i,j]<-
          sig1^2*G[i,j]
        
        if(distG[i,j]<10^(-10)){b1_hat[i,j]<-0}else{
          b1_hat[i,j]<-   #1-   ((1-exp(-alp*distG[i,j]/2))/ (alp*distG[i,j]/2))
            -((th_0*(exp(alp*distG[i,j]/2) - 1) + y_0)*th_0*exp(-alp*distG[i,j]/2) 
              - ((sig1*exp(alp*distG[i,j]/2) - sig1)*rho_y_th*sig + 
                   (alp*exp(alp*distG[i,j]/2) - alp)*th_0^2 + alp*th_0*y_0 + sig1^2 +
                   (alp*sig1^2*distG[i,j]/2 - sig1^2)*exp(alp*distG[i,j]/2))*exp(-alp*distG[i,j]/2)/alp)/(sig1^2*distG[i,j]/2) # we got problem here
        }
        cov_ri_rj[i,j]<-
          cov_yi_yj[i,j]-2*b1_hat[i,j]*cov_y_th[i,j]+(b1_hat[i,j])^2*cov_thi_thj[i,j]
        
      }#end of j loop
    }#end of i loop
    
    return(cov_ri_rj)
  }#end of CovRes fcn
  
  
  NegLogLike<-function(x,regb=regb,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor,dsm=dsm){
    badval<-(0.5)*.Machine$double.xmax
    alp<-x[1]
    sig<-x[2]  
    b1<-x[3]
    
    p_est<-c(alp,sig,b1)
    #print(p_est)
    V<-CovRes(p_est,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor)
    #V<-as.matrix(nearPD(V)$mat)
    #CovRes(p_est,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor)
    
    regb[2]<-b1
    #put reg b here modified 
    negloglike<-n/2*log(2*pi)+1/2*log(abs(det(V)))
    negloglike<-negloglike+1/2*t(response-dsm%*%regb)%*%pseudoinverse(V)%*%(response-dsm%*%regb)
    #print(negloglike)
    if(alp<0 || alp>100|| sig<0 ||  !is.finite(negloglike) || negloglike <= -1000 ){
      return(badval)
    }
    
    matrix.condition <- kappa(V, exact=TRUE)
    #print(paste("oubm matrix condition: ", log(matrix.condition) ))
    if(log(matrix.condition)>2){ #2 is the precision
      #print("entering the smooth spline method")
      proportions <- seq(from=1, to=0, length.out=101) 
      lnl.vector <- rep(NA, length(proportions))
      max.diff <- 0
      
      
      proportions <- seq(from=1, to=0, length.out=101) 
      lnl.vector <- rep(NA, length(proportions))
      max.diff <- 0
      for(i in sequence(length(proportions))) {
        V.modified.by.proportions<-(1-proportions[i]) * V + proportions[i] * diag(dim(V)[1]) * diag(V)
        try(local.lnl <- (n/2)*log(2*pi)+(1/2)*t(response-dsm%*%regb)%*%pseudoinverse(V.modified.by.proportions)%*%(response-dsm%*%regb) + (1/2)*log(abs(det(V.modified.by.proportions)))) 
        if(i>6) {
          very.local.lnl <- lnl.vector[(i-6):(i-1)]
          max.diff <- NA
          try(max.diff <- max(abs(very.local.lnl[-1] - very.local.lnl[-length(very.local.lnl)]))) #looking locally for jumps in the likelihood)
          current.diff <- NA
          try(current.diff <- abs(local.lnl - lnl.vector[i-1]))
          if(!is.na(max.diff) && !is.na(current.diff)) {
            if(current.diff > 2 * max.diff) {
              #print(paste("breaking after ", i))
              break() #the modified matrix is still poorly conditioned, so stop here  
            }
          }	
        }
        lnl.vector[i] <- local.lnl
      }
      
      proportions<-proportions[which(!is.na(lnl.vector))]
      lnl.vector<-lnl.vector[which(!is.na(lnl.vector))]
      proportions<-proportions[which(is.finite(lnl.vector))]
      lnl.vector<-lnl.vector[which(is.finite(lnl.vector))]
      
      NegLogML <- predict(smooth.spline(proportions, lnl.vector), data.frame(proportions =0.000))$y
      #plot(c(0, proportions), c(NegLogML, lnl.vector), type="n")
      #points(proportions, lnl.vector, pch=20)
      #points(0, NegLogML, col="red")
    }#end of precision
    
    return(negloglike[1])
  }
  #tree in newick format with branch length and topology
  G<-vcv(x.tre)
  #G<-G#/max(G)
  distG<-2*(max(G)-G)
  n<-dim(G)[1]
  one<-array(1,c(n,1))
  solveG<-solve(G)
  
  # BMmean<-t(one)%*%solveG%*%predictor
  # BMmean
  # 
  # BMmean<-BMmean/ c(t(one)%*%solveG%*%one)
  # BMmean
  # 
  # BMmean<-c(BMmean)*one
  # BMsigma<-t(predictor-BMmean)%*%solveG%*%(predictor-BMmean)/n
  # BMvcv<-c(BMsigma)*G
  
  #PUT BM est here  function with return
  
  dsm<-cbind(one,predictor)#design matrix for initial estimate of b
  b_ini<- solve(t(dsm)%*%dsm)%*%t(dsm)%*%response
  #abline(a=b_ini[1],b=b_ini[2])
  #print(b_ini)  
  #p0<-c(rexp(1),sd(predictor),b_ini[2])
  p0<-c(runif(0,20), runif(max(response)- min(response))) #alpha, sigma 
  #print(p0)
  GlobalAttempt<-GlobalOptim(NegLogLike,CovRes,precision=1,p0=p0,regb=b_ini,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor,dsm=dsm,model="OUBM") 
  #print(GlobalAttempt)
  MLEs<-GlobalAttempt$MLEs
  #print(MLEs)
  est_V<-GlobalAttempt$est_V
  dsm<-GlobalAttempt$dsm
  b_est<-GlobalAttempt$b_est
  #print(b_est)
  r_sq_value<-r_sq(b_est=b_est,est_V=est_V,dsm=dsm,predictor=predictor,response=response)
  b_est[2]<-b_est[2]*phyloCorrectFactor.OUBM(MLEs$par,predictor=predictor,x.tre=x.tre)
  #abline(a=b_est[1],b=b_est[2])
  
  AICc_value<-AICc(n=n,k=3,MLEs$value)
  result<-list(MLEs=MLEs,b_est=b_est,est_V=est_V,dsm=dsm,r_sq_value=r_sq_value,AICc_value=AICc_value)
  return(result)
  }#oubm
#-----------------------------------------------------------
#OUOU model
#Input Tree, Trait:response, precitor
#output MLEs, Likelihood, regression coefficient
phyloCorrectFactor.OUOU<-function(x,predictor=predictor,x.tre=x.tre){
  alp<-x[1]
  bet<-x[2]
  sig<-x[3]
  b1<-x[4]
  
  rho_y_th<-0
  y_0<-0
  th_0<-0
  
  G<-vcv(x.tre)
  G<-G/max(G)
  distG<-2*(1-G)
  n<-dim(G)[1]
  one<-array(1,c(n,1))
  solveG<-solve(G)
  
  inv_A_bet<-pseudoinverse(OUvcv(bet,G=G,distG=distG))
  mean_predictor<-c(t(one)%*%inv_A_bet%*%predictor)
  mean_predictor<-c(mean_predictor/t(one)%*%inv_A_bet%*%one)
  var_predictor<-t(predictor-mean_predictor*one)%*%inv_A_bet%*%(predictor-mean_predictor*one)/n
  sig1<-abs(b1*sqrt(abs(var_predictor)))
  
  cov_y_th<-
    -((alp - bet)*y_0*exp(bet*1) + (alp*exp(alp*1) - alp*exp(bet*1))*th_0)*th_0*exp(-alp*1 - 2*bet*1)/(alp - bet) - 1/2*((alp^2 + alp*bet)*sig1^2*exp(alp*1) - 2*(alp^2*bet - bet^3)*th_0*y_0*exp(bet*1) - 2*((alp^2*bet + alp*bet^2)*exp(alp*1) - (alp^2*bet + alp*bet^2)*exp(bet*1))*th_0^2 - (2*(alp*bet - bet^2)*sig1*rho_y_th*sig + (alp^2 - alp*bet)*sig1^2)*exp(alp*1 + 2*bet*1) - 2*(alp*bet*sig1^2 - (alp*bet - bet^2)*sig1*rho_y_th*sig)*exp(bet*1))*exp(-alp*1 - 2*bet*1)/(alp^2*bet - bet^3)
  var_th<-
    -1/2*bet*sig1^2*(exp(-2*bet*1) - 1)*1
  return(c(cov_y_th/var_th))
}


OUmean<-function(x,predictor=predictor,one=one,G=G,distG=distG){
  alp<-x[1]
  V<-OUvcv(alp,G=G,distG=distG)
  InvOUvcv<-pseudoinverse(V)
  mean_predictor<-t(one)%*%V%*%predictor
  mean_predictor<-c(mean_predictor/(t(one)%*%V%*%one))
  return(mean_predictor*one) 
}

ouou<-function(x.tre,response=response, predictor=predictor){
  CovRes<-function(x,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor){
    alp<-x[1]
    bet<-x[2]
    sig<-x[3]
    b1<-x[4]
    
    rho_y_th<-0
    y_0<-0
    th_0<-0
    
    #b1<-0.92636400
    ###LINEAR###
    #linear relationship between optium and predictor
    inv_A_bet<-pseudoinverse(OUvcv(bet,G=G,distG=distG))
    mean_predictor<-c(t(one)%*%inv_A_bet%*%predictor)
    mean_predictor<-c(mean_predictor/t(one)%*%inv_A_bet%*%one)
    var_predictor<-t(predictor-mean_predictor*one)%*%inv_A_bet%*%(predictor-mean_predictor*one)/n
    sig1<-abs(b1*sqrt(abs(var_predictor)))
    #print(sig1)
    
    cov_yi_yj<-array(0,c(n,n))
    cov_y_th<-array(0,c(n,n))
    cov_thi_thj<-array(0,c(n,n))
    b1_hat<-array(0,c(n,n))
    cov_ri_rj<-array(0,c(n,n))  
    
    for(i in 1:n){
      for(j in 1:n){
        
        cov_yi_yj[i,j]<-
          -(2*((alp - bet)*y_0*exp(bet*G[i,j]) + (alp*exp(alp*G[i,j]) - alp*exp(bet*G[i,j]))*th_0)*th_0*exp(-alp*G[i,j] - 2*bet*G[i,j])/(alp - bet) + ((alp^2 + alp*bet)*sig1^2*exp(alp*G[i,j]) - 2*(alp^2*bet - bet^3)*th_0*y_0*exp(bet*G[i,j]) - 2*((alp^2*bet + alp*bet^2)*exp(alp*G[i,j]) - (alp^2*bet + alp*bet^2)*exp(bet*G[i,j]))*th_0^2 - (2*(alp*bet - bet^2)*sig1*rho_y_th*sig + (alp^2 - alp*bet)*sig1^2)*exp(alp*G[i,j] + 2*bet*G[i,j]) - 2*(alp*bet*sig1^2 - (alp*bet - bet^2)*sig1*rho_y_th*sig)*exp(bet*G[i,j]))*exp(-alp*G[i,j] - 2*bet*G[i,j])/(alp^2*bet - bet^3))*(alp*exp(1/2*alp*distG[i,j]) - alp*exp(1/2*bet*distG[i,j]))*exp(-alp*distG[i,j] - 1/2*bet*distG[i,j])/(alp - bet) - 1/2*(2*th_0^2*exp(-2*bet*G[i,j]) - (2*bet*th_0^2 + sig1^2*exp(2*bet*G[i,j]) - sig1^2)*exp(-2*bet*G[i,j])/bet)*(alp*exp(1/2*alp*distG[i,j]) - alp*exp(1/2*bet*distG[i,j]))^2*exp(-alp*distG[i,j] - bet*distG[i,j])/(alp - bet)^2 - 1/2*(((alp^4 + alp^3*bet)*sig1^2*exp(2*alp*G[i,j]) - 2*(alp^4*bet - alp^3*bet^2 - alp^2*bet^3 + alp*bet^4)*y_0^2*exp(2*bet*G[i,j]) - 2*((alp^4*bet + alp^3*bet^2)*exp(2*alp*G[i,j]) - 2*(alp^4*bet + alp^3*bet^2)*exp(alp*G[i,j] + bet*G[i,j]) + (alp^4*bet + alp^3*bet^2)*exp(2*bet*G[i,j]))*th_0^2 - 4*((alp^4*bet - alp^2*bet^3)*exp(alp*G[i,j] + bet*G[i,j]) - (alp^4*bet - alp^2*bet^3)*exp(2*bet*G[i,j]))*th_0*y_0 - 4*(alp^3*bet*sig1^2 - (alp^3*bet - alp^2*bet^2)*sig1*rho_y_th*sig)*exp(alp*G[i,j] + bet*G[i,j]) - (2*(alp^3*bet - alp*bet^3)*sig1*rho_y_th*sig - (alp^3*bet + alp^2*bet^2)*sig1^2 - (alp^3*bet - alp^2*bet^2 - alp*bet^3 + bet^4)*sig^2 + (2*(alp^3*bet - 2*alp^2*bet^2 + alp*bet^3)*sig1*rho_y_th*sig + (alp^4 - 2*alp^3*bet + alp^2*bet^2)*sig1^2 + (alp^3*bet - alp^2*bet^2 - alp*bet^3 + bet^4)*sig^2)*exp(2*alp*G[i,j]))*exp(2*bet*G[i,j]))*exp(-2*alp*G[i,j] - 2*bet*G[i,j])/(alp^4*bet - alp^3*bet^2 - alp^2*bet^3 + alp*bet^4) + 2*((alp - bet)*y_0*exp(bet*G[i,j]) + (alp*exp(alp*G[i,j]) - alp*exp(bet*G[i,j]))*th_0)^2*exp(-2*alp*G[i,j] - 2*bet*G[i,j])/(alp - bet)^2)*exp(-alp*distG[i,j])
        
        #if(alp==bet){
        #  cov_y_th[i,j]<-(sig1^2/(2*bet^3))*(2*bet^3*G[i,j]*exp(bet*G[i,j])+bet^2*exp(bet*G[i,j])-bet^2*exp(3*bet*G[i,j]))
        #}else{
        cov_y_th[i,j]<-
          -((alp - bet)*y_0*exp(bet*G[i,j]) + (alp*exp(alp*G[i,j]) - alp*exp(bet*G[i,j]))*th_0)*th_0*exp(-alp*G[i,j] - 2*bet*G[i,j])/(alp - bet) - 1/2*((alp^2 + alp*bet)*sig1^2*exp(alp*G[i,j]) - 2*(alp^2*bet - bet^3)*th_0*y_0*exp(bet*G[i,j]) - 2*((alp^2*bet + alp*bet^2)*exp(alp*G[i,j]) - (alp^2*bet + alp*bet^2)*exp(bet*G[i,j]))*th_0^2 - (2*(alp*bet - bet^2)*sig1*rho_y_th*sig + (alp^2 - alp*bet)*sig1^2)*exp(alp*G[i,j] + 2*bet*G[i,j]) - 2*(alp*bet*sig1^2 - (alp*bet - bet^2)*sig1*rho_y_th*sig)*exp(bet*G[i,j]))*exp(-alp*G[i,j] - 2*bet*G[i,j])/(alp^2*bet - bet^3)
        #  }
        cov_thi_thj[i,j]<-
          -1/2*bet*sig1^2*(exp(-2*bet*G[i,j]) - 1)*exp(-bet*distG[i,j])
        
        if(distG[i,j]<10^(-10)){
          b1_hat[i,j]<- (2*alp*(-alp/2 +3/2*bet) +alp^2+alp*bet )  / (alp^2-bet^2)
          
        }else{
          b1_hat[i,j]<-
            (2*((alp - bet)*y_0*exp(bet*distG[i,j]/2) + (alp*exp(alp*distG[i,j]/2) - alp*exp(bet*distG[i,j]/2))*th_0)*th_0*exp(-alp*distG[i,j]/2 - 2*bet*distG[i,j]/2)/(alp - bet) - (2*alp*bet*sig1^2*exp(bet*distG[i,j]/2) + (alp^2 - alp*bet)*sig1^2*exp(alp*distG[i,j]/2 + 2*bet*distG[i,j]/2) - (alp^2 + alp*bet)*sig1^2*exp(alp*distG[i,j]/2) + 2*(alp^2*bet - bet^3)*th_0*y_0*exp(bet*distG[i,j]/2) + 2*((alp*bet - bet^2)*sig1*exp(alp*distG[i,j]/2 + 2*bet*distG[i,j]/2) - (alp*bet - bet^2)*sig1*exp(bet*distG[i,j]/2))*rho_y_th*sig + 2*((alp^2*bet + alp*bet^2)*exp(alp*distG[i,j]/2) - (alp^2*bet + alp*bet^2)*exp(bet*distG[i,j]/2))*th_0^2)*exp(-alp*distG[i,j]/2 - 2*bet*distG[i,j]/2)/(alp^2*bet - bet^3))/(2*th_0^2*exp(-2*bet*distG[i,j]/2) - (2*bet*th_0^2 + sig1^2*exp(2*bet*distG[i,j]/2) - sig1^2)*exp(-2*bet*distG[i,j]/2)/bet)
        }    
        cov_ri_rj[i,j]<-
          cov_yi_yj[i,j]-2*b1_hat[i,j]*cov_y_th[i,j]+(b1_hat[i,j])^2*cov_thi_thj[i,j]
        
      }#end of j loop
    }#end of i loop
    
    return(cov_ri_rj)
  }#end of CovRes fcn
  
  
  #NegLogLike(p0,regb=b_ini,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor,dsm=dsm)
  
  
  
  NegLogLike<-function(x,regb=regb,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor,dsm=dsm){
    badval<-(0.5)*.Machine$double.xmax
    alp<-x[1]
    bet<-x[2]
    sig<-x[3]
    b1<-x[4]
    
    p_est<-c(alp,bet,sig,b1)
    V<-CovRes(p_est,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor)
    #print(V)
    #print(c(alp,bet,sig,b1))
    #V<-as.matrix(nearPD(V)$mat)
    
    regb[2]<-b1
    #print("alp==bet?")
    #print(alp==bet)
    negloglike<-n/2*log(2*pi)+1/2*log(abs(det(V)))
    negloglike<-negloglike+1/2*t(response-dsm%*%regb)%*%pseudoinverse(V)%*%(response-dsm%*%regb)
    if(alp<0 ||alp>100||bet<0 ||bet>100|| alp==bet ||sig<0 || !is.finite(negloglike) || negloglike <= -1000){
      return(badval)
    }
    
    matrix.condition <- kappa(V, exact=TRUE)
    #print(paste("ouou matrix condition: ",log(matrix.condition) ))
    if(log(matrix.condition)>2){ #2 is the precision
      proportions <- seq(from=1, to=0, length.out=101) 
      lnl.vector <- rep(NA, length(proportions))
      max.diff <- 0
      
      
      proportions <- seq(from=1, to=0, length.out=101) 
      lnl.vector <- rep(NA, length(proportions))
      max.diff <- 0
      for(i in sequence(length(proportions))) {
        V.modified.by.proportions<-(1-proportions[i]) * V + proportions[i] * diag(dim(V)[1]) * diag(V)
        try(local.lnl <- (n/2)*log(2*pi)+(1/2)*t(response-dsm%*%regb)%*%pseudoinverse(V.modified.by.proportions)%*%(response-dsm%*%regb) + (1/2)*log(abs(det(V.modified.by.proportions)))) 
        if(i>6) {
          very.local.lnl <- lnl.vector[(i-6):(i-1)]
          max.diff <- NA
          try(max.diff <- max(abs(very.local.lnl[-1] - very.local.lnl[-length(very.local.lnl)]))) #looking locally for jumps in the likelihood)
          current.diff <- NA
          try(current.diff <- abs(local.lnl - lnl.vector[i-1]))
          if(!is.na(max.diff) && !is.na(current.diff)) {
            if(current.diff > 2 * max.diff) {
              #print(paste("breaking after ", i))
              break() #the modified matrix is still poorly conditioned, so stop here  
            }
          }	
        }
        lnl.vector[i] <- local.lnl
      }
      
      proportions<-proportions[which(!is.na(lnl.vector))]
      lnl.vector<-lnl.vector[which(!is.na(lnl.vector))]
      
      proportions<-proportions[which(is.finite(lnl.vector))]
      lnl.vector<-lnl.vector[which(is.finite(lnl.vector))]
      
      
      NegLogML <- predict(smooth.spline(proportions, lnl.vector), data.frame(proportions =0.000))$y
    }#end of precision
    #print(negloglike)
    return(negloglike[1])
  }
  #tree in newick format with branch length and topology
  G<-vcv(x.tre)
  G<-G/max(G)
  distG<-2*(1-G)
  n<-dim(G)[1]
  one<-array(1,c(n,1))
  
  dsm<-cbind(one,predictor)#design matrix for initial estimate of b
  b_ini<- solve(t(dsm)%*%dsm)%*%t(dsm)%*%response
  #print(b_ini)
  #p0<-c(rexp(1),rexp(1),sd(response),b_ini[2])
  p0<-c(0.05,0.12,0.1,0.72)
  GlobalAttempt<-GlobalOptim(NegLogLike,CovRes,precision=1,p0=p0,regb=b_ini,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor,dsm=dsm,model="OUOU") 
  
  MLEs<-GlobalAttempt$MLEs
  est_V<-GlobalAttempt$est_V
  dsm<-GlobalAttempt$dsm
  b_est<-GlobalAttempt$b_est
  #print(b_est)
  AICc_value<-AICc(n=n,k=4,MLEs$value)
  r_sq_value<-r_sq(b_est=b_est,est_V=est_V,dsm=dsm,predictor=predictor,response=response)
  b_est[2]<-phyloCorrectFactor.OUOU(MLEs$par,predictor=predictor,x.tre=x.tre)*b_est[2]
  #abline(a=b_est[1],b=b_est[2],col="red")
  
  oumean<-OUmean(MLEs$par[2],predictor=predictor,one=one,G=G,distG=distG)
  ousigma<-OUvcv(MLEs$par[2],G=G,distG=distG)
  
  result<-list(MLEs=MLEs,b_est=b_est,est_V=est_V,r_sq_value=r_sq_value,AICc_value=AICc_value,dsm=dsm,oumean=oumean,ousigma=ousigma)
  return(result)
}#ouou

#-----------------------------------------------------------
#OUBMBM model
#Input Tree, Trait:response, precitor
#output MLEs, Likelihood, regression coefficient
phyloCorrectFactor.OUBMBM<-function(x,predictor=predictor,x.tre=x.tre){
  alp<-x[1]
  xi<-x[2]
  b1<-x[3]
  
  rho_y_th<-0
  rho_y_sig<-0
  rho_th_sig<-0
  y_0<-0
  th_0<-0
  sig_0<-0
  
  G<-vcv(x.tre)
  G<-G/max(G)
  distG<-2*(1-G)
  n<-dim(G)[1]
  one<-array(1,c(n,1))
  solveG<-solve(G)
  
  inv_BMG<-pseudoinverse(G)
  mean_predictor<-c(t(one)%*%inv_BMG%*%predictor)
  mean_predictor<-c(mean_predictor/t(one)%*%inv_BMG%*%one)
  var_predictor<-t(predictor-mean_predictor*one)%*%inv_BMG%*%(predictor-mean_predictor*one)/n
  sig_th<-abs(b1*sqrt(abs(var_predictor)))
  
  cov_y_th<-
    -(th_0*(exp(alp*1) - 1) + y_0)*th_0*exp(-alp*1) - (rho_y_th*sig_0*sig_th - (alp*exp(alp*1) - alp)*th_0^2 - alp*th_0*y_0 - sig_th^2 - (alp*sig_th^2*1 + rho_y_th*sig_0*sig_th - sig_th^2)*exp(alp*1))*exp(-alp*1)/alp  
  var_th<-sig_th^2*1 
  return(c(cov_y_th/var_th))
}
oubmbm<-function(x.tre,response=response, predictor=predictor){
  CovRes<-function(x,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor){
    alp<-x[1]
    #sig_th<-x[2]
    xi <-x[2]
    #b1<-x[3]
    
    rho_y_th<-0
    rho_y_sig<-0
    rho_th_sig<-0
    y_0<-0
    th_0<-0
    sig_0<-0
    
    inv_BMG<-pseudoinverse(G)
    mean_predictor<-c(t(one)%*%inv_BMG%*%predictor)
    mean_predictor<-c(mean_predictor/t(one)%*%inv_BMG%*%one)
    var_predictor<-t(predictor-mean_predictor*one)%*%inv_BMG%*%(predictor-mean_predictor*one)/n
    sig_th<-abs(b1*sqrt(abs(var_predictor)))
    
    #fixed initial values, it will be cancelled eventaully for CovRes
    #or use the OU model estimate  see rmk1 in Jhwueng and Maroulas 2014
    # sig_0<-t(Y-y_0*one)%*%inv_A_alp%*%(Y-y_0*one)/n    #OUOU has not this item
    cov_yi_yj<-array(0,c(n,n))
    cov_y_th<-array(0,c(n,n))
    cov_thi_thj<-array(0,c(n,n))
    b1_hat<-array(0,c(n,n))
    cov_ri_rj<-array(0,c(n,n))  
    
    for(i in 1:n){
      for(j in 1:n){
        
        cov_yi_yj[i,j]<-
          sig_th^2*G[i,j]*(exp(1/2*alp*distG[i,j]) - 1)^2*exp(-alp*distG[i,j]) - 2*((th_0*(exp(alp*G[i,j]) - 1) + y_0)*th_0*exp(-alp*G[i,j]) + (rho_y_th*sig_0*sig_th - (alp*exp(alp*G[i,j]) - alp)*th_0^2 - alp*th_0*y_0 - sig_th^2 - (alp*sig_th^2*G[i,j] + rho_y_th*sig_0*sig_th - sig_th^2)*exp(alp*G[i,j]))*exp(-alp*G[i,j])/alp)*(exp(1/2*alp*distG[i,j]) - 1)*exp(-alp*distG[i,j]) - 1/4*(4*(th_0*(exp(alp*G[i,j]) - 1) + y_0)^2*exp(-2*alp*G[i,j]) - (4*alp*rho_y_th*sig_0*sig_th + 4*alp^2*y_0^2 - 2*alp*sig_0^2 - 2*alp*sig_th^2 + 4*(alp^2*exp(2*alp*G[i,j]) - 2*alp^2*exp(alp*G[i,j]) + alp^2)*th_0^2 + ((2*alp*G[i,j] - 1)*exp(2*alp*G[i,j]) + 1)*xi^2 + 8*(alp^2*exp(alp*G[i,j]) - alp^2)*th_0*y_0 + 2*(2*alp^2*sig_th^2*G[i,j] + 2*alp*rho_y_th*sig_0*sig_th + alp*sig_0^2 - 3*alp*sig_th^2)*exp(2*alp*G[i,j]) - 8*(alp*rho_y_th*sig_0*sig_th - alp*sig_th^2)*exp(alp*G[i,j]))*exp(-2*alp*G[i,j])/alp^2)*exp(-alp*distG[i,j])
        
        cov_y_th[i,j]<-
          -(th_0*(exp(alp*G[i,j]) - 1) + y_0)*th_0*exp(-alp*G[i,j]) - (rho_y_th*sig_0*sig_th - (alp*exp(alp*G[i,j]) - alp)*th_0^2 - alp*th_0*y_0 - sig_th^2 - (alp*sig_th^2*G[i,j] + rho_y_th*sig_0*sig_th - sig_th^2)*exp(alp*G[i,j]))*exp(-alp*G[i,j])/alp
        
        cov_thi_thj[i,j]<-
          sig_th^2*G[i,j]
        if(distG[i,j]<10^(-10)){b1_hat[i,j]<-0}else{
          b1_hat[i,j]<-
            -((th_0*(exp(alp*distG[i,j]/2) - 1) + y_0)*th_0*exp(-alp*distG[i,j]/2) - (rho_y_th*sig_0*sig_th*(exp(alp*distG[i,j]/2) - 1) + ((alp*distG[i,j]/2 - 1)*exp(alp*distG[i,j]/2) + 1)*sig_th^2 + (alp*exp(alp*distG[i,j]/2) - alp)*th_0^2 + alp*th_0*y_0)*exp(-alp*distG[i,j]/2)/alp)/(distG[i,j]/2*sig_th^2)
        }
        cov_ri_rj[i,j]<-
          cov_yi_yj[i,j]-2*b1_hat[i,j]*cov_y_th[i,j]+(b1_hat[i,j])^2*cov_thi_thj[i,j]
        
      }#end of j loop
    }#end of i loop
    
    #print(cov_ri_rj)
    return(cov_ri_rj)
  }#end of CovRes fcn
  
  
  
  NegLogLike<-function(x,regb=regb,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor,dsm=dsm){
    badval<-(0.5)*.Machine$double.xmax
    alp<-x[1]
    #sig_th<-x[2]
    xi <-x[2]
    #b1<-x[3]
    
    
    #regb[2]<-b1
    
    p_est<-c(alp,xi)
    V<-CovRes(p_est,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor)
    
    negloglike<-n/2*log(2*pi)+1/2*log(abs(det(V)))
    negloglike<-negloglike+1/2*t(response-dsm%*%regb)%*%pseudoinverse(V)%*%(response-dsm%*%regb)
    #print(negloglike)
    if(alp<0||alp>100 ||xi<0 || !is.finite(negloglike)|| negloglike <= -1000) {
      negloglike<-badval 
    }
    
    matrix.condition <- kappa(V, exact=TRUE)
   # print(paste("oubmbm matrix condition: ", log(matrix.condition) ))
    if(log(matrix.condition)>2){ #2 is the precision
      proportions <- seq(from=1, to=0, length.out=101) 
      lnl.vector <- rep(NA, length(proportions))
      max.diff <- 0
      
      
      proportions <- seq(from=1, to=0, length.out=101) 
      lnl.vector <- rep(NA, length(proportions))
      max.diff <- 0
      for(i in sequence(length(proportions))) {
        V.modified.by.proportions<-(1-proportions[i]) * V + proportions[i] * diag(dim(V)[1]) * diag(V)
        try(local.lnl <- (n/2)*log(2*pi)+(1/2)*t(response-dsm%*%regb)%*%pseudoinverse(V.modified.by.proportions)%*%(response-dsm%*%regb) + (1/2)*log(abs(det(V.modified.by.proportions)))) 
        if(i>6) {
          very.local.lnl <- lnl.vector[(i-6):(i-1)]
          max.diff <- NA
          try(max.diff <- max(abs(very.local.lnl[-1] - very.local.lnl[-length(very.local.lnl)]))) #looking locally for jumps in the likelihood)
          current.diff <- NA
          try(current.diff <- abs(local.lnl - lnl.vector[i-1]))
          if(!is.na(max.diff) && !is.na(current.diff)) {
            if(current.diff > 2 * max.diff) {
              #print(paste("breaking after ", i))
              break() #the modified matrix is still poorly conditioned, so stop here	
            }
          }	
        }
        lnl.vector[i] <- local.lnl
      }
      
      proportions<-proportions[which(!is.na(lnl.vector))]
      lnl.vector<-lnl.vector[which(!is.na(lnl.vector))]
      
      proportions<-proportions[which(is.finite(lnl.vector))]
      lnl.vector<-lnl.vector[which(is.finite(lnl.vector))]
      
      
      NegLogML <- predict(smooth.spline(proportions, lnl.vector), data.frame(proportions =0.000))$y
      }#end of precision
    
    #V<-as.matrix(nearPD(V)$mat)
    

    #print(negloglike)
    return(negloglike[1])
  }#end of likelihood
  
  #tree in newick format with branch length and topology
  G<-vcv(x.tre)
  G<-G/max(G)
  distG<-2*(1-G)
  n<-dim(G)[1]
  one<-array(1,c(n,1))
  solveG<-solve(G)
  
  ###############################
  OUBMBMmean<-t(one)%*%solveG%*%predictor
  OUBMBMmean<-OUBMBMmean/t(one)%*%solveG%*%one
  OUBMBMmean<-c(OUBMBMmean)*one
  OUBMBMsigma<-t(predictor-OUBMBMmean)%*%solveG%*%(predictor-OUBMBMmean)/n
  OUBMBMvcv<-c(OUBMBMsigma)*G
  ###############################
  
  dsm<-cbind(one,predictor)#design matrix for initial estimate of b
  b_ini<- solve(t(dsm)%*%dsm)%*%t(dsm)%*%response
  #print(b_ini)
  #p0<-c(rexp(1),sd(predictor),b_ini[2])
  p0<-c(0.05,0.3,0.72)
  GlobalAttempt<-GlobalOptim(NegLogLike,CovRes,precision=1,p0=p0,regb=b_ini,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor,dsm=dsm,model="OUBMBM") 
  
  MLEs<-GlobalAttempt$MLEs
  #print(MLEs)
  est_V<-GlobalAttempt$est_V
  #print(est_V)
  dsm<-GlobalAttempt$dsm
  #print(dsm)
  b_est<-GlobalAttempt$b_est
  #print(b_est)
  
  AICc_value<-AICc(n=n,k=3,MLEs$value)  
  r_sq_value<-r_sq(b_est=b_est,est_V=est_V,dsm=dsm,predictor=predictor,response=response)
  b_est[2]<-phyloCorrectFactor.OUBMBM(MLEs$par,predictor=predictor,x.tre=x.tre)*b_est[2]
  #abline(a=b_est[1],b=b_est[2],col="blue")
  
  result<-list(MLEs=MLEs,b_est=b_est,est_V=est_V,r_sq_value=r_sq_value,AICc_value=AICc_value,dsm=dsm,OUBMBMmean=OUBMBMmean,OUBMBMsigma=OUBMBMsigma,OUBMBMvcv=OUBMBMvcv)
  return(result)
}#oubmbm
#-----------------------------------------------------------
#OUOUBM model
#Input Tree, Trait:response, precitor
#output MLEs, Likelihood, regression coefficient
#-----------------------------------------------------------
phyloCorrectFactor.OUOUBM<-function(x,predictor=predictor,x.tre=x.tre){
  alp<-x[1]
  bet<-x[2]
  xi<-x[3]
  b1<-x[4]
  
  rho_y_th<-0
  rho_y_sig<-0
  rho_th_sig<-0
  
  y_0<-0
  th_0<-0
  sig_0<-0
  
  G<-vcv(x.tre)
  G<-G/max(G)
  distG<-2*(1-G)
  n<-dim(G)[1]
  one<-array(1,c(n,1))
  solveG<-solve(G)
  
  inv_A_bet<-pseudoinverse(OUvcv(bet,G=G,distG=distG))
  mean_predictor<-c(t(one)%*%inv_A_bet%*%predictor)
  mean_predictor<-c(mean_predictor/t(one)%*%inv_A_bet%*%one)
  var_predictor<-t(predictor-mean_predictor*one)%*%inv_A_bet%*%(predictor-mean_predictor*one)/n
  sig1<-abs(b1*sqrt(abs(var_predictor)))
  
  cov_y_th<-
    -((alp - bet)*y_0*exp(bet*1) + (alp*exp(alp*1) - alp*exp(bet*1))*th_0)*th_0*exp(-alp*1 - 2*bet*1)/(alp - bet) + 1/2*(2*(alp^2*bet - bet^3)*th_0*y_0*exp(alp*1 + bet*1) + (alp^2 - alp*bet)*sig1^2 - 2*(alp^2*bet - alp*bet^2 - (alp^2*bet - alp*bet^2)*exp(alp*1 + bet*1))*th_0^2 + 2*(alp*bet*sig1^2 + (alp*bet + bet^2)*sig1*rho_y_th*sig_0)*exp(alp*1 + bet*1) - (2*(alp*bet + bet^2)*sig1*rho_y_th*sig_0 + (alp^2 + alp*bet)*sig1^2)*exp(2*bet*1))*exp(-2*bet*1)/(alp^2*bet - bet^3)
  var_th<-
    -1/2*bet*sig1^2*(exp(-2*bet*1) - 1)*exp(-bet*0)
  return(c(cov_y_th/var_th))
  }

ououbm<-function(x.tre,response=response, predictor=predictor){
  CovRes<-function(x,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor){
    alp<-x[1]
    bet<-x[2]
    #sig1<-x[3]
    xi<-x[3]
    b1<-x[4]
    
    #print(c(alp,bet,xi,b1))
    
    rho_y_th<-0
    rho_y_sig<-0
    rho_th_sig<-0
    
    y_0<-0
    th_0<-0
    sig_0<-0
    
    
    inv_A_bet<-pseudoinverse(OUvcv(bet,G=G,distG=distG))
    mean_predictor<-c(t(one)%*%inv_A_bet%*%predictor)
    mean_predictor<-c(mean_predictor/t(one)%*%inv_A_bet%*%one)
    var_predictor<-t(predictor-mean_predictor*one)%*%inv_A_bet%*%(predictor-mean_predictor*one)/n
    #print(inv_A_alp)
    sig1<-abs(b1*sqrt(abs(var_predictor)))
    
    #print(sig1)
    
    cov_yi_yj<-array(0,c(n,n))
    cov_y_th<-array(0,c(n,n))
    cov_thi_thj<-array(0,c(n,n))
    b1_hat<-array(0,c(n,n))
    cov_ri_rj<-array(0,c(n,n))  
    
    
    for(i in 1:n){
      for(j in 1:n){
        
        cov_yi_yj[i,j]<-
          -(2*((alp - bet)*y_0*exp(bet*G[i,j]) + (alp*exp(alp*G[i,j]) - alp*exp(bet*G[i,j]))*th_0)*th_0*exp(-alp*G[i,j] - 2*bet*G[i,j])/(alp - bet) - (2*(alp^2*bet - bet^3)*th_0*y_0*exp(alp*G[i,j] + bet*G[i,j]) + (alp^2 - alp*bet)*sig1^2 - 2*(alp^2*bet - alp*bet^2 - (alp^2*bet - alp*bet^2)*exp(alp*G[i,j] + bet*G[i,j]))*th_0^2 + 2*(alp*bet*sig1^2 + (alp*bet + bet^2)*sig1*rho_y_th*sig_0)*exp(alp*G[i,j] + bet*G[i,j]) - (2*(alp*bet + bet^2)*sig1*rho_y_th*sig_0 + (alp^2 + alp*bet)*sig1^2)*exp(2*bet*G[i,j]))*exp(-2*bet*G[i,j])/(alp^2*bet - bet^3))*(alp*exp(1/2*alp*distG[i,j]) - alp*exp(1/2*bet*distG[i,j]))*exp(-alp*distG[i,j] - 1/2*bet*distG[i,j])/(alp - bet) - 1/2*(2*th_0^2*exp(-2*bet*G[i,j]) - (2*bet*th_0^2 + sig1^2*exp(2*bet*G[i,j]) - sig1^2)*exp(-2*bet*G[i,j])/bet)*(alp*exp(1/2*alp*distG[i,j]) - alp*exp(1/2*bet*distG[i,j]))^2*exp(-alp*distG[i,j] - bet*distG[i,j])/(alp - bet)^2 + 1/4*((2*(3*alp^5 - alp^4*bet)*sig1^2*exp(2*alp*G[i,j]) + (3*alp^3*bet - alp^2*bet^2 - 3*alp*bet^3 + bet^4 - (3*alp^3*bet - alp^2*bet^2 - 3*alp*bet^3 + bet^4 - 2*(3*alp^4*bet - alp^3*bet^2 - 3*alp^2*bet^3 + alp*bet^4)*G[i,j])*exp(2*alp*G[i,j]))*xi^2*exp(2*bet*G[i,j]) + 4*(3*alp^5*bet - alp^4*bet^2 - 3*alp^3*bet^3 + alp^2*bet^4)*y_0^2*exp(2*bet*G[i,j]) - 4*((3*alp^5*bet - alp^4*bet^2)*exp(2*alp*G[i,j]) - 2*(alp^5*bet - alp^4*bet^2)*exp(3*alp*G[i,j] + bet*G[i,j]) - (alp^5*bet + alp^4*bet^2)*exp(2*bet*G[i,j]))*th_0^2 + 8*((alp^5*bet - alp^3*bet^3)*exp(3*alp*G[i,j] + bet*G[i,j]) - (alp^5*bet - alp^3*bet^3)*exp(2*bet*G[i,j]))*th_0*y_0 + 8*(alp^4*bet*sig1^2 + (alp^4*bet + alp^3*bet^2)*sig1*rho_y_th*sig_0)*exp(3*alp*G[i,j] + bet*G[i,j]) + 2*(2*(alp^4*bet - alp^2*bet^3)*sig1*rho_y_th*sig_0 - (alp^4*bet + alp^3*bet^2)*sig1^2 - (3*alp^4*bet - alp^3*bet^2 - 3*alp^2*bet^3 + alp*bet^4)*sig_0^2 - (2*(3*alp^4*bet + 2*alp^3*bet^2 - alp^2*bet^3)*sig1*rho_y_th*sig_0 + (3*alp^5 + 2*alp^4*bet - alp^3*bet^2)*sig1^2 - (3*alp^4*bet - alp^3*bet^2 - 3*alp^2*bet^3 + alp*bet^4)*sig_0^2)*exp(2*alp*G[i,j]))*exp(2*bet*G[i,j]))*exp(-2*alp*G[i,j] - 2*bet*G[i,j])/(3*alp^5*bet - alp^4*bet^2 - 3*alp^3*bet^3 + alp^2*bet^4) - 4*((alp - bet)*y_0*exp(bet*G[i,j]) + (alp*exp(alp*G[i,j]) - alp*exp(bet*G[i,j]))*th_0)^2*exp(-2*alp*G[i,j] - 2*bet*G[i,j])/(alp - bet)^2)*exp(-alp*distG[i,j])
        
        cov_y_th[i,j]<-
          -((alp - bet)*y_0*exp(bet*G[i,j]) + (alp*exp(alp*G[i,j]) - alp*exp(bet*G[i,j]))*th_0)*th_0*exp(-alp*G[i,j] - 2*bet*G[i,j])/(alp - bet) + 1/2*(2*(alp^2*bet - bet^3)*th_0*y_0*exp(alp*G[i,j] + bet*G[i,j]) + (alp^2 - alp*bet)*sig1^2 - 2*(alp^2*bet - alp*bet^2 - (alp^2*bet - alp*bet^2)*exp(alp*G[i,j] + bet*G[i,j]))*th_0^2 + 2*(alp*bet*sig1^2 + (alp*bet + bet^2)*sig1*rho_y_th*sig_0)*exp(alp*G[i,j] + bet*G[i,j]) - (2*(alp*bet + bet^2)*sig1*rho_y_th*sig_0 + (alp^2 + alp*bet)*sig1^2)*exp(2*bet*G[i,j]))*exp(-2*bet*G[i,j])/(alp^2*bet - bet^3)
        
        cov_thi_thj[i,j]<-
          -1/2*bet*sig1^2*(exp(-2*bet*G[i,j]) - 1)*exp(-bet*distG[i,j])
        
        if(distG[i,j]<10^(-10)){
          b1_hat[i,j]<- (2*alp*(-alp/2 +3/2*bet) +alp^2+alp*bet )  / (alp^2-bet^2)
        }else{
          b1_hat[i,j]<-
            (2*((alp - bet)*y_0*exp(bet*distG[i,j]/2) + (alp*exp(alp*distG[i,j]/2) - alp*exp(bet*distG[i,j]/2))*th_0)*th_0*exp(-alp*distG[i,j]/2 - 2*bet*distG[i,j]/2)/(alp - bet) - (2*alp*bet*sig1^2*exp(alp*distG[i,j]/2 + bet*distG[i,j]/2) + 2*(alp^2*bet - bet^3)*th_0*y_0*exp(alp*distG[i,j]/2 + bet*distG[i,j]/2) - (alp^2 + alp*bet)*sig1^2*exp(2*bet*distG[i,j]/2) + (alp^2 - alp*bet)*sig1^2 + 2*((alp*bet + bet^2)*sig1*exp(alp*distG[i,j]/2 + bet*distG[i,j]/2) - (alp*bet + bet^2)*sig1*exp(2*bet*distG[i,j]/2))*rho_y_th*sig_0 - 2*(alp^2*bet - alp*bet^2 - (alp^2*bet - alp*bet^2)*exp(alp*distG[i,j]/2 + bet*distG[i,j]/2))*th_0^2)*exp(-2*bet*distG[i,j]/2)/(alp^2*bet - bet^3))/(2*th_0^2*exp(-2*bet*distG[i,j]/2) - (2*bet*th_0^2 + sig1^2*exp(2*bet*distG[i,j]/2) - sig1^2)*exp(-2*bet*distG[i,j]/2)/bet)
        }
        
        cov_ri_rj[i,j]<-
          cov_yi_yj[i,j]-2*b1_hat[i,j]*cov_y_th[i,j]+(b1_hat[i,j])^2*cov_thi_thj[i,j]
        
      }#end of j loop
    }#end of i loop
    #print(c(cov_yi_yj[i,j],-2*b1_hat[i,j]*cov_y_th[i,j],+(b1_hat[i,j])^2*cov_thi_thj[i,j]))
    #print(is.positive.definite (cov_ri_rj))
    return(cov_ri_rj)
  }#end of CovRes fcn
  
  
  NegLogLike<-function(x,regb=regb,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor,dsm=dsm){
    badval<-(0.5)*.Machine$double.xmax
    alp<-x[1]
    bet<-x[2]
    #sig1<-x[3]
    xi<-x[3]
    b1<-x[4]
    
    p_est<-c(alp,bet,xi,b1)
    V<-CovRes(p_est,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor)
    #V<-as.matrix(nearPD(V)$mat)
    
    
    regb[2]<-b1
    negloglike<-n/2*log(2*pi)+1/2*log(abs(det(V)))
    #print(V)
    negloglike<-negloglike+1/2*t(response-dsm%*%regb)%*%pseudoinverse(V)%*%(response-dsm%*%regb)
    #print("alp,bet,xi,negloglike")
    #print( c(alp,bet,xi,negloglike)  )
    if(alp<0||alp>100||bet<0||bet>100 ||xi<0|| !is.finite(negloglike) || negloglike <= -1000) {
      negloglike<-badval 
    }
    #print(negloglike)
    matrix.condition <- kappa(V, exact=TRUE)
    #print(paste("ououbm matrix condition: ", log(matrix.condition) ))
    if(log(matrix.condition)>2){ #2 is the precision
      proportions <- seq(from=1, to=0, length.out=101) 
      lnl.vector <- rep(NA, length(proportions))
      max.diff <- 0
      
      
      proportions <- seq(from=1, to=0, length.out=101) 
      lnl.vector <- rep(NA, length(proportions))
      max.diff <- 0
      for(i in sequence(length(proportions))) {
        V.modified.by.proportions<-(1-proportions[i]) * V + proportions[i] * diag(dim(V)[1]) * diag(V)
        
        
        try(local.lnl <- (n/2)*log(2*pi)+(1/2)*t(response-dsm%*%regb)%*%pseudoinverse(V.modified.by.proportions)%*%(response-dsm%*%regb) + (1/2)*log(abs(det(V.modified.by.proportions)))) 
        if(i>6) {
          very.local.lnl <- lnl.vector[(i-6):(i-1)]
          max.diff <- NA
          try(max.diff <- max(abs(very.local.lnl[-1] - very.local.lnl[-length(very.local.lnl)]))) #looking locally for jumps in the likelihood)
          current.diff <- NA
          try(current.diff <- abs(local.lnl - lnl.vector[i-1]))
          if(!is.na(max.diff) && !is.na(current.diff)) {
            if(current.diff > 2 * max.diff) {
              #print(paste("breaking after ", i))
              break() #the modified matrix is still poorly conditioned, so stop here  
            }
          }	
        }
        lnl.vector[i] <- local.lnl
      }
      
      proportions<-proportions[which(!is.na(lnl.vector))]
      lnl.vector<-lnl.vector[which(!is.na(lnl.vector))]
      
      proportions<-proportions[which(is.finite(lnl.vector))]
      lnl.vector<-lnl.vector[which(is.finite(lnl.vector))]
      
      
      NegLogML <- predict(smooth.spline(proportions, lnl.vector), data.frame(proportions =0.000))$y
    }#end of precision
    
    
    return(negloglike[1])
  }
  
  #tree in newick format with branch length and topology
  G<-vcv(x.tre)
  G<-G/max(G)
  distG<-2*(1-G)
  n<-dim(G)[1]
  one<-array(1,c(n,1))
  
  dsm<-cbind(one,predictor)#design matrix for initial estimate of b
  #print(dsm)
  #print(response)
  b_ini<- solve(t(dsm)%*%dsm)%*%t(dsm)%*%response
  #print(b_ini)
  #p0<-c(rexp(1),rexp(1),sd(predictor),b_ini[2])
  p0<-c(0.05,0.12,0.3,0.72)
  GlobalAttempt<-GlobalOptim(NegLogLike,CovRes,precision=1,p0=p0,regb=b_ini,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor,dsm=dsm,model="OUOUBM")   
  #NegLogLike(p0,regb=b_ini,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor,dsm=dsm)
  MLEs<-GlobalAttempt$MLEs
  est_V<-GlobalAttempt$est_V
  dsm<-GlobalAttempt$dsm
  b_est<-GlobalAttempt$b_est
  
  #print(b_est)
  AICc_value<-AICc(n=n,k=4,MLEs$value) 
  r_sq_value<-r_sq(b_est=b_est,est_V=est_V,dsm=dsm,predictor=predictor,response=response)
  
  #print(r_sq_value)
  b_est[2]<-phyloCorrectFactor.OUOUBM(MLEs$par,predictor=predictor,x.tre=x.tre)*b_est[2]
  #abline(a=b_est[1],b=b_est[2],col="purple")
  
  oumean<-OUmean(MLEs$par[2],predictor=predictor,one=one,G=G,distG=distG)
  ousigma<-OUvcv(MLEs$par[2],G=G,distG=distG)
  
  result<-list(MLEs=MLEs,b_est=b_est,est_V=est_V,r_sq_value=r_sq_value,AICc_value=AICc_value,dsm=dsm,oumean=oumean,ousigma=ousigma)
  return(result)
}#ououbm
#---------------------------------------------------------------------
#OUBMCIR model
#Input Tree, Trait:response, precitor1 , precitor2
#output MLEs, Likelihood, regression coefficient

phyloCorrectFactor.OUBMCIR <- function(x , predictor_1 = predictor_1 , predictor_2 = predictor_2 , x.tre = x.tre){
  alp<-x[1]
  bet<-x[2]
  th_1<-x[3]
  
  xi<-x[4]
  b1<-x[5]
  b2<-x[6]
  
  
  y_0<-0
  th_0<-0
  sig_0<-0
  
  G <- vcv(x.tre)
  G <- G/max(G)
  distG <- 2*(1-G)
  n <- dim(G)[1]
  one <- array(1,c(n,1))
  solveG <- solve(G)
  
  inv_BMG<-pseudoinverse(G)
  mean_predictor_1<-c(t(one)%*%inv_BMG%*%predictor_1)
  mean_predictor_1<-c(mean_predictor_1/(t(one)%*%inv_BMG%*%one))
  var_predictor_1<-t(predictor_1-mean_predictor_1*one)%*%inv_BMG%*%(predictor_1-mean_predictor_1*one)/n
  
  
  mean_predictor_2<-c(t(one)%*%inv_BMG%*%predictor_2)
  mean_predictor_2<-c(mean_predictor_2/(t(one)%*%inv_BMG%*%one))
  var_predictor_2<-t(predictor_2-mean_predictor_2*one)%*%inv_BMG%*%(predictor_2-mean_predictor_2*one)/n
  
  sig_1 <- sqrt(b1^2*var_predictor_1 + b2^2*var_predictor_2)
  
  cov_y_th <- -(th_0*(exp(alp*1) - 1) + y_0)*th_0*exp(-alp*1) +
    ((alp*exp(alp*1) - alp)*th_0^2 + alp*th_0*y_0 + sig_1^2 +
       (alp*sig_1^2*1 - sig_1^2)*exp(alp*1))*exp(-alp*1)/alp
  var_th <- sig_1^2*1
  return(c(cov_y_th/var_th))
  
}

oubmcir<-function(x.tre,response=response, predictor_1=predictor_1, predictor_2=predictor_2){
  CovRes <- function(x,n=n,G=G,distG=distG,one=one,response=response,predictor_1=predictor_1,predictor_2=predictor_2)
    {

     alp<-x[1]
     bet<-x[2]
     th_1<-x[3]
     
     xi<-x[4]
     b1<-x[5]
     b2<-x[6]
  
   
   
    
    
   
    
    y_0<-0
    th_0<-0
    sig_0<-0

    inv_BMG<-pseudoinverse(G)
    mean_predictor_1<-c(t(one)%*%inv_BMG%*%predictor_1)
    mean_predictor_1<-c(mean_predictor_1/(t(one)%*%inv_BMG%*%one))
    var_predictor_1<-t(predictor_1-mean_predictor_1*one)%*%inv_BMG%*%(predictor_1-mean_predictor_1*one)/n
  
  
    mean_predictor_2<-c(t(one)%*%inv_BMG%*%predictor_2)
    mean_predictor_2<-c(mean_predictor_2/(t(one)%*%inv_BMG%*%one))
    var_predictor_2<-t(predictor_2-mean_predictor_2*one)%*%inv_BMG%*%(predictor_2-mean_predictor_2*one)/n
  
    sig_1 <- sqrt(b1^2*var_predictor_1 + b2^2*var_predictor_2)
    cov_yi_yj<-array(0,c(n,n))
    cov_y_th<-array(0,c(n,n))
    cov_thi_thj<-array(0,c(n,n))
    b1_hat<-array(0,c(n,n))
    cov_ri_rj<-array(0,c(n,n))

  for(i in 1:n){
    for(j in 1:n){
  
  
  
    if(bet!=2*alp){
    cov_yi_yj[i,j] <-
    sig_1^2*G[i,j]*(exp(1/2*alp*distG[i,j]) - 1)^2*exp(-alp*distG[i,j]) -
    2*((th_0*(exp(alp*G[i,j]) - 1) + y_0)*th_0*exp(-alp*G[i,j]) -
    ((alp*exp(alp*G[i,j]) - alp)*th_0^2 + alp*th_0*y_0 + sig_1^2 +
    (alp*sig_1^2*G[i,j] - sig_1^2)*exp(alp*G[i,j]))*exp(-alp*G[i,j])/alp)*(exp(1/2*alp*distG[i,j]) - 1) -
    1/4*(4*(th_0*(exp(alp*G[i,j]) - 1) + y_0)^2*exp(-2*alp*G[i,j]) -
    (4*(2*(2*alp^2*bet^2 - alp*bet^3)*th_1^2*exp(bet*G[i,j]) -
    (4*alp^3*bet - alp*bet^3)*th_1*exp(bet*G[i,j]) + (2*alp^3*bet -
    alp^2*bet^2)*exp(bet*G[i,j]))*y_0^2*exp(2*bet*G[i,j]*th_1) -
    8*(2*(2*alp^2*bet^2 - alp*bet^3 - (2*alp^2*bet^2 -
    alp*bet^3)*exp(alp*G[i,j]))*th_0*th_1^2*exp(bet*G[i,j]) - (4*alp^3*bet
    - alp*bet^3 - (4*alp^3*bet -
    alp*bet^3)*exp(alp*G[i,j]))*th_0*th_1*exp(bet*G[i,j]) + (2*alp^3*bet -
    alp^2*bet^2 - (2*alp^3*bet -
    alp^2*bet^2)*exp(alp*G[i,j]))*th_0*exp(bet*G[i,j]))*y_0*exp(2*bet*G[i,j]*th_1
    ) - 2*(2*alp^2*bet - alp*bet^2 + 2*(2*alp^2*bet -
    alp*bet^2)*sig_0^2 - 2*(2*alp^2*bet -
    alp*bet^2)*sig_0)*th_1*exp(2*alp*G[i,j] + bet*G[i,j]) -
    (4*alp^2*sig_0*xi^2*exp(2*alp*G[i,j]) +
    4*(2*alp*bet^2*exp(2*alp*G[i,j]) - (bet^3 + (2*alp*bet^2 -
    bet^3)*exp(2*alp*G[i,j]))*exp(bet*G[i,j]))*th_1^3 - 4*(2*alp^3*bet -
    alp^2*bet^2 + (2*alp^3*bet - alp^2*bet^2)*exp(2*alp*G[i,j]) -
    2*(2*alp^3*bet - alp^2*bet^2)*exp(alp*G[i,j]))*th_0^2*exp(bet*G[i,j])
    - 2*(4*(2*alp^2*bet^2 - alp*bet^3 + (2*alp^2*bet^2 -
    alp*bet^3)*exp(2*alp*G[i,j]) - 2*(2*alp^2*bet^2 -
    alp*bet^3)*exp(alp*G[i,j]))*th_0^2*exp(bet*G[i,j]) +
    2*(2*alp*bet^2*sig_0 - alp*bet*xi^2 +
    2*alp^2*bet)*exp(2*alp*G[i,j]) - (4*alp*bet^2*sig_0 -
    bet^2*xi^2 - 8*(2*alp*bet^2 - bet^3)*sig_1^2*exp(alp*G[i,j]) +
    bet^3 + 2*(2*alp*bet^2 - bet^3)*sig_1^2 - (4*(2*alp^2*bet^2
    - alp*bet^3)*sig_1^2*G[i,j] - 4*alp^2*bet + bet^3 -
    6*(2*alp*bet^2 - bet^3)*sig_1^2 + (2*alp*bet -
    bet^2)*xi^2)*exp(2*alp*G[i,j]))*exp(bet*G[i,j]))*th_1^2 +
    (4*(4*alp^3*bet - alp*bet^3 + (4*alp^3*bet -
    alp*bet^3)*exp(2*alp*G[i,j]) - 2*(4*alp^3*bet -
    alp*bet^3)*exp(alp*G[i,j]))*th_0^2*exp(bet*G[i,j]) +
    4*(2*alp^2*bet*sig_0 - (alp*bet*sig_0 +
    alp^2)*xi^2)*exp(2*alp*G[i,j]) - (4*alp*bet^2*sig_0 -
    8*(4*alp^2*bet - bet^3)*sig_1^2*exp(alp*G[i,j]) + 4*(2*alp^2*bet
    - alp*bet^2)*sig_0^2 + 2*(4*alp^2*bet - bet^3)*sig_1^2 -
    (4*alp*bet*sig_0 + bet^2)*xi^2 - (4*(4*alp^3*bet -
    alp*bet^3)*sig_1^2*G[i,j] - 4*alp^2*bet + 2*alp*bet^2 -
    6*(4*alp^2*bet - bet^3)*sig_1^2 + (4*alp^2 -
    bet^2)*xi^2)*exp(2*alp*G[i,j]))*exp(bet*G[i,j]))*th_1 -
    (2*alp*bet*sig_0*xi^2 + 8*(2*alp^2*bet -
    alp*bet^2)*sig_1^2*exp(alp*G[i,j]) - 2*(2*alp^2*bet -
    alp*bet^2)*sig_0^2 - 2*(2*alp^2*bet - alp*bet^2)*sig_1^2 +
    (4*(2*alp^3*bet - alp^2*bet^2)*sig_1^2*G[i,j] - 6*(2*alp^2*bet -
    alp*bet^2)*sig_1^2 + (2*alp^2 -
    alp*bet)*xi^2)*exp(2*alp*G[i,j]))*exp(bet*G[i,j]))*exp(2*bet*G[i,j]*th_1)
    + (2*(2*alp^2*bet - alp*bet^2)*sig_0^2 - (2*alp^2 -
    alp*bet - 2*(2*alp^2 - alp*bet)*sig_0)*xi^2)*exp(2*alp*G[i,j]
    + bet*G[i,j]))*exp(-2*bet*G[i,j]*th_1)/(2*(2*alp^2*bet^2 -
    alp*bet^3)*th_1^2*exp(2*alp*G[i,j] + bet*G[i,j]) - (4*alp^3*bet -
    alp*bet^3)*th_1*exp(2*alp*G[i,j] + bet*G[i,j]) + (2*alp^3*bet -
    alp^2*bet^2)*exp(2*alp*G[i,j] + bet*G[i,j]))) *exp(alp*distG[i,j])
}
    
    
    cov_y_th[i,j] <-
    -(th_0*(exp(alp*G[i,j]) - 1) + y_0)*th_0*exp(-alp*G[i,j]) +
    ((alp*exp(alp*G[i,j]) - alp)*th_0^2 + alp*th_0*y_0 + sig_1^2 +
    (alp*sig_1^2*G[i,j] - sig_1^2)*exp(alp*G[i,j]))*exp(-alp*G[i,j])/alp

    cov_thi_thj[i,j] <-
    sig_1^2*G[i,j]
if(distG[i,j]<10^(-10)){b1_hat[i,j]<-0}else{
    b1_hat[i,j] <-
    -((th_0*(exp(alp*distG[i,j]/2) - 1) + y_0)*th_0*exp(-alp*distG[i,j]/2) -
    ((alp*exp(alp*distG[i,j]/2) - alp)*th_0^2 + alp*th_0*y_0 + sig_1^2 +
    (alp*distG[i,j]/2*sig_1^2 -
    sig_1^2)*exp(alp*distG[i,j]/2))*exp(-alp*distG[i,j]/2)/alp)/(distG[i,j]/2*sig_1^2)
}
    cov_ri_rj[i,j]<-
      cov_yi_yj[i,j]-2*b1_hat[i,j]*cov_y_th[i,j]+(b1_hat[i,j])^2*cov_thi_thj[i,j]
      }#end of j loop
  }#end of i loop
    return(cov_ri_rj)
}# end of CovRes function

  NegLogLike <- function(x,regb=regb,n=n,G=G,distG=distG,one=one,response=response,predictor_1=predictor_1,predictor_2=predictor_2,dsm=dsm){
    badval<-(0.5)*.Machine$double.xmax
    alp<-x[1]
    bet<-x[2]
    th_1<-x[3]
    
    xi<-x[4]
    b1<-x[5]
    b2<-x[6]
    p_est <- c(alp,bet,th_1,xi,b1,b2)
    V<-CovRes(p_est,n=n,G=G,distG=distG,one=one,response=response,predictor_1=predictor_1,predictor_2=predictor_2)
    regb[2]<-b1
    regb[3]<-b2
    negloglike<- n/2*log(2*pi)+1/2*log(abs(det(V)))
    negloglike<- negloglike + 1/2*t(response-dsm%*%regb)%*%pseudoinverse(V)%*%(response-dsm%*%regb)
    if(alp<0 || alp>100 || xi<0 || !is.finite(negloglike) || bet<0 || bet>100 || negloglike <= -1000){
      return(badval)
    }  
  
  # matrix.condition <- kappa(V , exact = TRUE)
  # if(log(matrix.condition)>2){
  #   proportions <- seq(from=1, to=0, length.out=101) 
  #   lnl.vector <- rep(NA, length(proportions))
  #   max.diff <- 0
  #   
  #   
  #   proportions <- seq(from=1, to=0, length.out=101) 
  #   lnl.vector <- rep(NA, length(proportions))
  #   max.diff <- 0
  #   for(i in sequence(length(proportions))){
  #     V.modified.by.proportions <- (1-proportions[i])*V+proportions[i]*diag(dim(V)[1])*diag(V)
  #     try(local.lnl <- (n/2)*log(2*pi)+(1/2)*t(response-dsm%*%regb)%*%pseudoinverse(V.modified.by.proportions)%*%(response-dsm%*%regb)+
  #           (1/2)*log(abs(det(V.modified.by.proportions))))
  #     if(i>6){
  #       very.local.lnl <- lnl.vector[(i-6):(i-1)]
  #       max.diff<-NA
  #       try(max.diff<-max(abs(very.local.lnl[-1] - very.local.lnl[-length(very.local.lnl)])))
  #       current.diff <- NA
  #       
  #     }
  #   }
  # }
  # 
    return(negloglike[1])
      }

  #tree in newick format with branch length and topology
  G<-vcv(x.tre)
  G<-G/max(G)
  distG<-2*(1-G)
  n<-dim(G)[1]
  one<-array(1,c(n,1))
  solveG<-solve(G)
  
  BMmean_1<-t(one)%*%solveG%*%predictor_1
  BMmean_1<-BMmean_1/t(one)%*%solveG%*%one
  BMmean_1<-c(BMmean_1)*one
  BMsigma_1<-t(predictor_1-BMmean_1)%*%solveG%*%(predictor_1-BMmean_1)/n
  BMvcv_1<-c(BMsigma_1)*G
  
  BMmean_2<-t(one)%*%solveG%*%predictor_2
  BMmean_2<-BMmean_2/t(one)%*%solveG%*%one
  BMmean_2<-c(BMmean_2)*one
  BMsigma_2<-t(predictor_2-BMmean_2)%*%solveG%*%(predictor_2-BMmean_2)/n
  BMvcv_2<-c(BMsigma_2)*G
  # put BMOUCIR est here function with return
  dsm<-cbind(one,predictor_1,predictor_2)#design matrix for initial estimate of b
  b_ini<- solve(t(dsm)%*%dsm)%*%t(dsm)%*%response
  p0 <- c(2,3,5,10,20,50)


  GlobalAttempt<-
    GlobalOptim(NegLogLike,CovRes,precision=1,p0=p0,regb=b_ini,n=n,G=G,distG=distG,one=one,response=response,predictor_1=predictor_1,predictor_2=predictor_2,dsm=dsm,model="OUBMCIR")
    MLEs <- GlobalAttempt$MLEs
    est_V<-GlobalAttempt$est_V
    dsm<-GlobalAttempt$dsm
    b_est<-GlobalAttempt$b_est
    r_sq_value<-r_sq(b_est=b_est,est_V=est_V,dsm=dsm,predictor_1=predictor_1,predictor_2=predictor_2,response = response)
    b_est[2]<-b_est[2]*phyloCorrectFactor.OUBMCIR(MLEs$par,predictor_1=predictor_1,predictor_2=predictor_2,x.tre=x.tre)
    AICc_value <- AICc(n=n,k=3,MLEs$value)
  result<-
    list(MLEs=MLEs,b_est=b_est,est_V=est_V,dsm=dsm,r_sq_value=r_sq_value,AICc_value=AICc_value,BMmean_1=BMmean_1,BMmean_2=BMmean_2,BMsigma_1=BMsigma_1,BMsigma_2=BMsigma_2,BMvcv_1=BMvcv_1,BMvcv_2=BMvcv_2)
  return(result)


    }# end of oubmcir function
#-----------------------------------------------------
#OUOUCIR model
#Input Tree, Trait:response, precitor1 , precitor2
#output MLEs, Likelihood, regression coefficient

#OU model
OUvcv <- function(alp , G=G , distG=distG){ #vcv for OU , not inculde sigma_sq
  A <- (1 - exp(-2*alp*G))/(2*alp)
  A <- A*exp(-alp*distG)
  return(A)
}

#---------------------------------------------------

phyloCorrectFactor.OUOUCIR <- function(x , predictor_1 = predictor_1 , predictor_2 = predictor_2 , x.tre = x.tre){
  
  alp<-x[1]
  bet<-x[2]
  delta<-x[3]
  th_1<-x[4]
  
  xi<-x[5]
  b1<-x[6]
  b2<-x[7]
  
  y_0<-0
  th_0<-0
  sig_0<-0
  
  G<-vcv(x.tre)
  G<-G/max(G)
  distG<-2*(1-G)
  n<-dim(G)[1]
  one<-array(1,c(n,1))
  
  inv_A_bet<-pseudoinverse(OUvcv(bet,G=G,distG=distG))
  mean_predictor_1<-c(t(one)%*%inv_A_bet%*%predictor_1)
  mean_predictor_1<-c(mean_predictor_1/(t(one)%*%inv_A_bet%*%one))
  var_predictor_1<-t(predictor_1-mean_predictor_1*one)%*%inv_A_bet%*%(predictor_1-mean_predictor_1*one)/n
  
  
  mean_predictor_2<-c(t(one)%*%inv_A_bet%*%predictor_2)
  mean_predictor_2<-c(mean_predictor_2/(t(one)%*%inv_A_bet%*%one))
  var_predictor_2<-t(predictor_2-mean_predictor_2*one)%*%inv_A_bet%*%(predictor_2-mean_predictor_2*one)/n
  
  sig_1 <- sqrt(b1^2*var_predictor_1 + b2^2*var_predictor_2)

  cov_y_th <- -((alp - bet)*y_0*exp(bet*1) + (alp*exp(alp*1) -
            alp*exp(bet*1))*th_0)*th_0*exp(-alp*1 - 2*bet*1)/(alp - bet) + 1/2*(2*alp*bet*sig_1^2*exp(bet*1) + (alp^2 -
            alp*bet)*sig_1^2*exp(alp*1 + 2*bet*1) - (alp^2 + alp*bet)*sig_1^2*exp(alp*1) + 2*(alp^2*bet -bet^3)*th_0*y_0*exp(bet*1) + 2*((alp^2*bet + alp*bet^2)*exp(alp*1) - (alp^2*bet +alp*bet^2)*exp(bet*1))*th_0^2)*exp(-alp*1 -
            2*bet*1)/(alp^2*bet - bet^3)
  var_th <- -th_0^2*exp(-2*bet*1) + 1/2*(2*bet*th_0^2 + sig_1^2*exp(2*bet*1) - sig_1^2)*exp(-2*bet*1)/bet
  
  return(c(cov_y_th/var_th))
  
}


#---------------------------------------------------

OUmean <- function(x , predictor_1=predictor_1 , predictor_2=predictor_2 , one = one , G=G , distG=distG){
  alp <- x[1]
  V <- OUvcv(alp,G=G,distG=distG)
  InvOUvcv <- pseudoinverse(V)
  mean_predictor_1 <- t(one)%*%V%*%predictor_1
  mean_predictor_1 <- c(mean_predictor_1/t(one)%*%V%*%one)
  
  mean_predictor_2 <- t(one)%*%V%*%predictor_2
  mean_predictor_2 <- c(mean_predictor_2/t(one)%*%V%*%one)
  return(c(mean_predictor_1*one , mean_predictor_2*one))
  # return(mean_predictor_2*one)
  
}


ououcir <- function(x.tre , response = response , predictor_1 = predictor_1 , predictor_2 = predictor_2){
  
  CovRes <- function(x,n=n,G=G,distG=distG,one=one,response=response,predictor_1=predictor_1,predictor_2=predictor_2)
  {
   # x<-c(2,3,4,1,1,2,5)
    alp<-x[1]
    bet<-x[2]
    delta<-x[3]
    th_1<-x[4]
  
    xi<-x[5]
    b1<-x[6]
    b2<-x[7]
  
    y_0<-0
    th_0<-0
    sig_0<-0
  #print("hello")
  inv_A_bet<-pseudoinverse(OUvcv(bet,G=G,distG=distG))
  #print(inv_A_bet)
  mean_predictor_1<-c(t(one)%*%inv_A_bet%*%predictor_1)
  mean_predictor_1<-c(mean_predictor_1/(t(one)%*%inv_A_bet%*%one))
  var_predictor_1<-t(predictor_1-mean_predictor_1*one)%*%inv_A_bet%*%(predictor_1-mean_predictor_1*one)/n
  
  
  mean_predictor_2<-c(t(one)%*%inv_A_bet%*%predictor_2)
  mean_predictor_2<-c(mean_predictor_2/(t(one)%*%inv_A_bet%*%one))
  var_predictor_2<-t(predictor_2-mean_predictor_2*one)%*%inv_A_bet%*%(predictor_2-mean_predictor_2*one)/n
  
  print((b1^2*var_predictor_1 + b2^2*var_predictor_2) )
  sig_1 <- sqrt(b1^2*var_predictor_1 + b2^2*var_predictor_2)
  cov_yi_yj<-array(0,c(n,n))
  cov_y_th<-array(0,c(n,n))
  cov_thi_thj<-array(0,c(n,n))
  b1_hat<-array(0,c(n,n))
  cov_ri_rj<-array(0,c(n,n))

  for(i in 1:n){
    for(j in 1:n){
 
      if(2*alp != delta){
      A_1 <- exp(2*alp*G[i,j] - delta*G[i,j])/(2*alp - delta) -1/(2*alp - delta)
      }else{
      A_1 <- G[i,j]
      }
       if(alp != bet ){
         A_2 <-2*(alp*exp(alp*G[i,j]) -alp*exp(bet*G[i,j]))*th_0*th_0*exp(-alp*G[i,j] - 2*bet*G[i,j])/(alp -bet) 
         A_3 <-(alp^2 -alp*bet)*sig_1^2*exp(alp*G[i,j] + 2*bet*G[i,j])/ (alp^2*bet - bet^3)* exp(-alp*G[i,j] - 2*bet*G[i,j])
         A_4 <-  (alp*exp(1/2*alp*distG[i,j]) - alp*exp(1/2*bet*distG[i,j]))^2*exp(-alp*distG[i,j] - bet*distG[i,j])/(alp - bet)^2
         #A_5 <-  2*((alp^2*bet +alp*bet^2)*exp(alp*G[i,j]) - (alp^2*bet + alp*bet^2)*exp(bet*G[i,j]))*th_0^2*exp(-alp*G[i,j] - 2*bet*G[i,j])/(alp^2*bet - bet^3)
        # A_6 <- (2*alp*bet*exp(bet*G[i,j]) - (alp^2 +alp*bet)*exp(alp*G[i,j]))/(alp^2*bet - bet^3)
         }else{
         A_2 <- G[i,j]*exp(bet*G[i,j])
         A_3 <- alp/(bet*(alp+bet))*sig_1^2*exp(alp*G[i,j] + 2*bet*G[i,j])* exp(-alp*G[i,j] - 2*bet*G[i,j])
         A_4 <- -(alp/2*G[i,j])^2*exp(alp*G[i,j])
         #A_5 <- 2*bet*exp(bet*G[i,j])*th_0^2*exp(-alp*G[i,j] - 2*bet*G[i,j])
        # A_6 <- exp(bet*G[i,j])
       }

    cov_yi_yj[i,j] <-
    1/2*(2*y_0^2*exp(-2*alp*G[i,j]) + (alp*xi^2*th_1*(exp(2*alp*G[i,j])/alp
    - 1/alp)/delta + 2*alp*th_1^2*(exp(2*alp*G[i,j])/alp - 1/alp) +
    4*alp*sig_0*xi^2*(A_1)/delta 
  - 8*alp*th_1^2*(A_1) -
    4*alp*sig_0*th_1*(exp(2*alp*G[i,j] - 2*delta*G[i,j])/(alp - delta) -
    1/(alp - delta)) + alp*xi^2*th_1*(exp(2*alp*G[i,j] -
    2*delta*G[i,j])/(alp - delta) - 1/(alp - delta))/delta +
    2*alp*th_1^2*(exp(2*alp*G[i,j] - 2*delta*G[i,j])/(alp - delta) - 1/(alp -
    delta)) + xi^2*th_1*(exp(2*alp*G[i,j])/alp - 1/alp)/delta +
    2*th_1^2*(exp(2*alp*G[i,j])/alp - 1/alp) +
    4*sig_0*xi^2*(A_1)/delta + 4*(2*delta*sig_0 -
    xi^2)*alp*th_1*(A_1)/delta - 8*th_1^2*(A_1) +
    2*sig_0^2*(exp(2*alp*G[i,j] - 2*delta*G[i,j])/(alp - delta) - 1/(alp -
    delta)) - 2*sig_0*xi^2*(exp(2*alp*G[i,j] - 2*delta*G[i,j])/(alp - delta) -
    1/(alp - delta))/delta - 4*sig_0*th_1*(exp(2*alp*G[i,j] -
    2*delta*G[i,j])/(alp - delta) - 1/(alp - delta)) +
    xi^2*th_1*(exp(2*alp*G[i,j] - 2*delta*G[i,j])/(alp - delta) - 1/(alp -
    delta))/delta + 2*th_1^2*(exp(2*alp*G[i,j] - 2*delta*G[i,j])/(alp - delta) -
    1/(alp - delta)) + 4*(2*delta*sig_0 - xi^2)*th_1*(A_1)/delta +
    2*(delta*sig_0^2 - sig_0*xi^2)*alp*(exp(2*alp*G[i,j] -
    2*delta*G[i,j])/(alp - delta) - 1/(alp - delta))/delta)*exp(-2*alp*G[i,j])
    - 2*((alp - bet)*y_0*exp(bet*G[i,j]) + (alp*exp(alp*G[i,j]) -
    alp*exp(bet*G[i,j]))*th_0)^2*exp((-2*alp*G[i,j] - 2*bet*G[i,j])/(alp -bet)^2))*(alp - bet)^2*exp(bet*distG[i,j]) - 
    (
      2*(y_0*exp(bet*G[i,j]))+A_2 - A_3-2*th_0*y_0*exp(bet*G[i,j])*exp(-alp*G[i,j] - 2*bet*G[i,j]) 
    -(2*alp*bet*exp(bet*G[i,j]) - (alp^2 +alp*bet)*exp(alp*G[i,j]))/(alp^2*bet - bet^3)*sig_1^2*exp(-alp*G[i,j] - 2*bet*G[i,j])+ 2*((alp^2*bet +alp*bet^2)*exp(alp*G[i,j]) - (alp^2*bet + alp*bet^2)*exp(bet*G[i,j]))*th_0^2*exp(-alp*G[i,j] - 2*bet*G[i,j])/(alp^2*bet - bet^3)
  )*(alp*exp(1/2*alp*distG[i,j])-alp*exp(1/2*bet*distG[i,j]))*exp(-1/2*alp*distG[i,j]) -
    1/2*(2*th_0^2*exp(-2*bet*G[i,j]) - (2*bet*th_0^2 + sig_1^2*exp(2*bet*G[i,j])- sig_1^2)*exp(-2*bet*G[i,j])/bet)*A_4
      
  
    if(alp==bet){
      cov_y_th[i,j]<-(sig_1^2/(2*bet^3))*(2*bet^3*G[i,j]*exp(bet*G[i,j])+bet^2*exp(bet*G[i,j])-bet^2*exp(3*bet*G[i,j]))
    }else{
    cov_y_th[i,j] <-
    -((alp - bet)*y_0*exp(bet*G[i,j]) + (alp*exp(alp*G[i,j]) -
    alp*exp(bet*G[i,j]))*th_0)*th_0*exp(-alp*G[i,j] - 2*bet*G[i,j])/(alp -
    bet) + 1/2*(2*alp*bet*sig_1^2*exp(bet*G[i,j]) + (alp^2 -
    alp*bet)*sig_1^2*exp(alp*G[i,j] + 2*bet*G[i,j]) - (alp^2 +
    alp*bet)*sig_1^2*exp(alp*G[i,j]) + 2*(alp^2*bet -
    bet^3)*th_0*y_0*exp(bet*G[i,j]) + 2*((alp^2*bet +
    alp*bet^2)*exp(alp*G[i,j]) - (alp^2*bet +
    alp*bet^2)*exp(bet*G[i,j]))*th_0^2)*exp(-alp*G[i,j] -
    2*bet*G[i,j])/(alp^2*bet - bet^3)
}
    cov_thi_thj[i,j] <-
     -1/2*bet*sig_1^2*(exp(-2*bet*G[i,j]) - 1)*exp(-bet*distG[i,j])
    if(distG[i,j]<10^(-10)){
      b1_hat[i,j]<- (2*alp*(-alp/2 +3/2*bet) +alp^2+bet*bet )  / (alp^2-bet^2)
    }else{
        if( alp!=bet ){
          B_1<- ((alp - bet)*y_0*exp(bet*distG[i,j]/2) + (alp*exp(alp*distG[i,j]/2) -alp*exp(bet*distG[i,j]/2))*th_0)/(alp-bet)*th_0*exp(-alp*distG[i,j]/2 -2*bet*distG[i,j]/2)
        }else{
          B_1<- (y_0*exp(bet*distG[i,j]/2) + exp(bet*distG[i,j]/2))*th_0*exp(-alp*distG[i,j]/2 -2*bet*distG[i,j]/2)
        }
      
      b1_hat[i,j] <-
    (
      2*B_1 -
    (2*alp*bet*sig_1^2*exp(bet*distG[i,j]/2) + (alp^2 -alp*bet)*sig_1^2*exp(alp*distG[i,j]/2 + 2*bet*distG[i,j]/2)
      - (alp^2 +alp*bet)*sig_1^2*exp(alp*distG[i,j]/2) + 
        2*(alp^2*bet -bet^3)*th_0*y_0*exp(bet*distG[i,j]/2) + 
        2*((alp^2*bet +alp*bet^2)*exp(alp*distG[i,j]/2)-
    (alp^2*bet +alp*bet^2)*exp(bet*distG[i,j]/2))*th_0^2 )*exp(-alp*distG[i,j]/2 - 2*bet*distG[i,j]/2)/(alp^2*bet - bet^3)
    )/(2*th_0^2*exp(-2*bet*distG[i,j]/2) - (2*bet*th_0^2 + sig_1^2*exp(2*bet*distG[i,j]/2) -sig_1^2)*exp(-2*bet*distG[i,j]/2)/bet)
    }
    
    cov_ri_rj[i,j]<-
      cov_yi_yj[i,j]-2*b1_hat[i,j]*cov_y_th[i,j]+(b1_hat[i,j])^2*cov_thi_thj[i,j]
   }# end of j loop
 }#end of i loop
return(cov_ri_rj)
}#end of CovRes  function

NegLogLike <- function(x,regb=regb,n=n,G=G,distG=distG,one=one,response=response,predictor_1=predictor_1,predictor_2=predictor_2,dsm=dsm){
  badval<-(0.5)*.Machine$double.xmax
 
  alp<-x[1]
  bet<-x[2]
  delta<-x[3]
  th_1<-x[4]
  
  xi<-x[5]
  b1<-x[6]
  b2<-x[7]
  p_est<-c(alp,bet,delta,th_1,xi,b1,b2)

  V<-CovRes(p_est,n=n,G=G,distG=distG,one=one,response=response,predictor_1=predictor_1,predictor_2=predictor_2)
  regb[2] <- b1
  regb[3] <- b2
  negloglike<-n/2*log(2*pi)+1/2*log(abs(det(V)))
  negloglike<-negloglike+1/2*t(response-dsm%*%regb)%*%pseudoinverse(V)%*%(response-dsm%*%regb)
  
  if(alp<0 || alp>100 || bet<0 || bet>100 || delta<0 || delta>100 || xi<0 || !is.finite(negloglike)|| negloglike <= -1000){
    return(badval)
  }# if 
return(negloglike[1])
  } # end of NegLogLike function

#tree in newick format with branch length and topology
G<-vcv(x.tre)
G<-G/max(G)
distG<-2*(1-G)
n<-dim(G)[1]
one<-array(1,c(n,1))

dsm<-cbind(one,predictor_1,predictor_2)#design matrix for initial estimate of b
b_ini<- solve(t(dsm)%*%dsm)%*%t(dsm)%*%response
p0 <- c(2,3,4,1,1,2,5)


GlobalAttempt<-
  GlobalOptim(NegLogLike,CovRes,precision = 1,p0=p0,regb=b_ini,n=n,G=G,distG=distG,one=one,response=response,predictor_1=predictor_1
              ,predictor_2=predictor_2,dsm = dsm,model = "OUOUCIR" )
MLEs <- GlobalAttempt$MLEs
est_V <- GlobalAttempt$est_V
dsm <- GlobalAttempt$dsm
b_est <- GlobalAttempt$b_est
AICc_value <- AICc(n=n,k=7,MLEs$value)
r_sq_value <- r_sq(b_est=b_est,est_V=est_V,dsm=dsm,predictor_1=predictor_1,predictor_2=predictor_2,response=response)
b_est[2] <- phyloCorrectFactor.OUOUCIR(MLEs$par,predictor_1 = predictor_1 , predictor_2 = predictor_2 , x.tre = x.tre)*b_est[2]
oumean <- OUmean(MLEs$par[2] , predictor_1 = predictor_1 , predictor_2=predictor_2 , one = one , G=G , distG=distG)
ousigma <- OUvcv(MLEs$par[2] , G=G , distG=distG)

result <- 
  list(MLEs=MLEs,b_est=b_est,est_V=est_V,r_sq_value=r_sq_value,AICc_value=AICc_value,dsm=dsm,oumean=oumean,ousigma=ousigma)
return(result)

  } # end of ououcir function

####function####



r_sq<-function(b_est=b_est,est_V=est_V,dsm=dsm,predictor=predictor,response=response){
  one<-array(1,c(length(response),1))
  inv_V_est<-pseudoinverse(est_V)
  mean_response<-(t(one)%*%inv_V_est%*%response)/sum(inv_V_est)
  SST<-t(response-c(mean_response))%*%inv_V_est%*%(response-c(mean_response))
  SSE<-t(response-dsm%*%b_est)%*%inv_V_est%*%(response-dsm%*%b_est)
  return((SST-SSE)/SST)
}

AICc<-function(n,k,negloglike){
  AIC<-2*k+negloglike
  AICc<-AIC+ 2*k*(k+1)/(n-k-1)
  return(AICc)
}


GenerateValues <- function(par, lower, upper, max.tries=100, expand.prob=0, examined.max, examined.min) {
  pass=FALSE
  tries=0
  while(!pass && tries<=max.tries) {
    tries <- tries+1
    pass=TRUE
    new.vals <- rep(NA, length(par))
    for(i in sequence(length(par))) {
      examined.max[i]<-max(0.001, examined.max[i])
      left_end<-min(max(lower[i], 0.9*examined.min[i]), min(upper[i], 1.1*examined.max[i]))
      right_end<-max(max(lower[i], 0.9*examined.min[i]), min(upper[i], 1.1*examined.max[i]))
      new.vals[i]<-runif(1,left_end,right_end)
      if(new.vals[i]<lower[i]) {
        pass=FALSE
      }
      if(new.vals[i]>upper[i]) {
        pass=FALSE
      }
    }
  }
  if(tries>max.tries) {
    return(NA)
  }
  return(new.vals)
}
#### change not yet 
GlobalOptim<-function(NegLogLike,CovRes,precision=precision,p0=p0,regb=regb,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor,dsm=dsm,model=model){
  free.parameters<-rep(TRUE, 5)
  names(free.parameters) <- c("alp", "bet", "sig", "xi", "b1")
  if(model=="OUBM") {
    free.parameters[which(names(free.parameters) == "bet")]<-FALSE
    free.parameters[which(names(free.parameters) == "xi")]<-FALSE
  }
  if(model=="OUOU") {
    free.parameters[which(names(free.parameters) == "xi")]<-FALSE
    
  }
  if(model== "OUBMBM") {
    free.parameters[which(names(free.parameters) == "sig")]<-FALSE
    free.parameters[which(names(free.parameters) == "bet")]<-FALSE
  }
  if(model== "OUOUBM") {
    free.parameters[which(names(free.parameters)=="sig")]<-FALSE
  }
  
  tmp_p0<-p0
  tmp_regb<-regb
  best.run<-list(value=NegLogLike(p0,regb=regb,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor,dsm=dsm),par=p0)
  print("starting likelihood")
  print(best.run$value)
  trial<-0
  impr_count<-0
  while(trial <=1){
    trial<-trial+1
    #print(paste("improvement search ",trial))
    tryout<-0
    while(tryout<=3){
      tryout<-tryout+1
      #print(paste("This is the search for", trial, "th improvement where the ",tryout, "th search were trying to get the convergence estimates"))
      #print("starting points")
      #print(p0)
      if(model=="OUBM"){dsm<-sim.dsm.oubm(p0,predictor=predictor,x.tre=x.tre)}
      if(model=="OUOU"){dsm<-sim.dsm.ouou(p0,predictor=predictor,x.tre=x.tre)}
      if(model=="OUBMBM"){dsm<-sim.dsm.oubmbm(p0,predictor=predictor,x.tre=x.tre)}
      if(model=="OUOUBM"){dsm<-sim.dsm.ououbm(p0,predictor=predictor,x.tre=x.tre)}
        
      new.run<-optim(p0,NegLogLike,method="Nelder-Mead",regb=regb,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor,dsm=dsm)       
     
      if(new.run$convergence==0){
        #print("get convergence estimates")
        #print(new.run$par)
        best.run<-new.run
        break}else{#print("not find convergence estimates yet, start from neighbor point")
                   p0<-GenerateValues(best.run$par, lower=c(0, 0, 0, 0, -10)[which(free.parameters)], upper=c(100,100,(range(response)[2]-range(response)[1]),(range(response)[2]-range(response)[1]),10 ), examined.max=10*best.run$par, examined.min=0.1*best.run$par)                 
                   #print(p0)
                   }
        }#end of tryout loop
    
    #if(abs(new.run$value)+0.0001<abs(best.run$value) && new.run$value>0 ){
    #if(abs(new.run$value)<abs(best.run$value) && new.run$value>0 ){
    #  impr_count<-impr_count+1
    #  print("find an improvement")
    #  print(paste("old value is", best.run$value  ,"new improvement is" ,new.run$value))
    #  best.run<-new.run
    #  print("try to search a better improvement")
    #  p0<-GenerateValues(best.run$par, lower=c(0, 0, 0, 0, -10)[which(free.parameters)], upper=c(100,100,(range(response)[2]-range(response)[1]),(range(response)[2]-range(response)[1]),10 ), examined.max=10*best.run$par, examined.min=0.1*best.run$par)
    #  print("start point")
    #  print(p0)
      #regb<-b_est+runif(1,-1,1)
      #p0[length(best.run$par)]<-best.run$par[length(best.run$par)]
      #regb[2]<-best.run$par[length(best.run$par)]
      #print(paste("retry generate values",p0))         
    #  if(impr_count==2){print("find two improvement shall be good")
    #                    print(paste("old value is", best.run$value  ,"new improvement is" ,new.run$value))                        
    #                    ;break}
    #  }else{
    #  print("not find improvement, keep searching")
    #    print( paste("best vs new ", best.run$value,new.run$value,sep=" ") )       
    #     if(trial %%3==0){
    #       p0<-best.run$par
            #regb<-b_est+runif(1,-1,1)
    #       p0[length(best.run$par)]<-best.run$par[length(best.run$par)]
    #       #regb[2]<-best.run$par[length(best.run$par)]
    #       print("research with mls values")
    #       print(p0)
    #       }else{
    #       p0<-GenerateValues(best.run$par, lower=c(0, 0, 0, 0, -10)[which(free.parameters)], upper=c(100,100,(range(response)[2]-range(response)[1]),(range(response)[2]-range(response)[1]),10 ), examined.max=10*best.run$par, examined.min=0.1*best.run$par)   
    #       print("research with generated values")
    #       print(p0)
    #       }
    #  }
    #if(trial==2){print("search does not find any improvement")}
    }#end of trial loop 
    # minLikeIndex<-  which(Likevalue==min(Likevalue))
  
  
  
  MLEs<-best.run #so either convergence estimate or the true parameter
  print("searched likelihood value")
  print(best.run$value)
  #print(model)
  #print("the estimator values are")
  print(MLEs$par)
  
  est_V<- CovRes(best.run$par,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor)
  #print(est_V)
  #we need to put differnt phyloCorrFact here
  if(model=="OUBM"){dsm<-sim.dsm.oubm(best.run$par,predictor=predictor,x.tre=x.tre)}
  if(model=="OUOU"){dsm<-sim.dsm.ouou(best.run$par,predictor=predictor,x.tre=x.tre)}
  if(model=="OUBMBM"){dsm<-sim.dsm.oubmbm(best.run$par,predictor=predictor,x.tre=x.tre)}
  if(model=="OUOUBM"){dsm<-sim.dsm.ououbm(best.run$par,predictor=predictor,x.tre=x.tre)}
  
  est_V_inv<-pseudoinverse(est_V) 
  b_est<-pseudoinverse(t(dsm)%*%est_V_inv%*%dsm)%*%t(dsm)%*%est_V_inv%*%response
  
  result<-list(MLEs=MLEs,est_V=est_V,dsm=dsm,b_est=b_est)
  return(result)
  }
#now we can do the package for this file

GetModels<-function(x.tre,predictor=predictor, response=response,model=model){
  if(model=="BM"){output<-bm(x.tre,response=response,predictor=predictor)}
  if(model=="OU"){output<-ou(x.tre,response=response,predictor=predictor)}
  if(model=="OUBM"){output<-oubm(x.tre,response=response,predictor=predictor)}
  if(model=="OUOU"){output<-ouou(x.tre,response=response,predictor=predictor)}
  if(model=="OUBMBM"){output<-oubmbm(x.tre,response=response,predictor=predictor)}
  if(model=="OUOUBM"){output<-ououbm(x.tre,response=response,predictor=predictor)}
  return(output)
  }

RegressionLine<-function(x.tre,predictor=predictor,response=response){
  numModel<-4
  Models<-array(c("OUBM","OUOU","OUBMBM","OUOUBM"),c(numModel,1))
  negloglike_array<-array(0,c(numModel,1))
  regCoef<-array(0,c(numModel,2))
  AICcArray<-array(0,c(numModel,1))
  r_sq_Array<-array(0,c(numModel,1))
  output.table<-NULL
  for(modelIndex in 1:length(Models)){ 
    model<-Models[modelIndex]
    #print(model)
    assign(paste("output",Models[modelIndex],sep=""),GetModels(x.tre,response=response,predictor=predictor,model=model))
    AICcArray[modelIndex,1]<-get(paste("output",Models[modelIndex],sep=""))$AICc_value
    r_sq_Array[modelIndex,1]<-get(paste("output",Models[modelIndex],sep=""))$r_sq_value
    regCoef[modelIndex,]<-get(paste("output",Models[modelIndex],sep=""))$b_est
  }
  OutTable<-cbind(Models,regCoef,AICcArray,r_sq_Array)
  colnames(OutTable)<-c("Models","b0","b1","AICc","r_sq")
  #print(OutTable)
  result<-list(outputOUBM=outputOUBM,outputOUOU=outputOUOU,outputOUBMBM=outputOUBMBM,outputOUOUBM=outputOUOUBM,OutTable=OutTable)
  return(result)
  }

sim.dsm.oubm<-function(model.params,predictor=predictor,x.tre=x.tre){
  alp<-model.params[1]
  sig<-model.params[2]  
  #b1<-x[3]

  n<-length(predictor)
  one<-array(1,c(n,1))
  #print(one)
  second<-array(phyloCorrectFactor.OUBM(c(alp,sig),predictor=predictor,x.tre=x.tre)*predictor)
  return(cbind(one,second))
  }

sim.dsm.ouou<-function(x,predictor=predictor,x.tre=x.tre){
  alp<-x[1]
  bet<-x[2]
  sig<-x[3]
  b1<-x[4]
  
  one<-array(1,c(length(predictor),1))
  second<-array(phyloCorrectFactor.OUOU(c(alp,bet,sig,b1),predictor=predictor,x.tre=x.tre)*predictor)
  return(cbind(one,second))
  }

sim.dsm.oubmbm<-function(x,predictor=predictor,x.tre=x.tre){
  alp<-x[1]
  xi <-x[2]
  b1<-x[3]
  
  one<-array(1,c(length(predictor),1))
  second<-array(phyloCorrectFactor.OUBMBM(c(alp,xi,b1),predictor=predictor,x.tre=x.tre)*predictor)
  return(cbind(one,second))
}

sim.dsm.ououbm<-function(x,predictor=predictor,x.tre=x.tre){
  alp<-x[1]
  bet<-x[2]
  xi<-x[3]
  b1<-x[4]
  
  one<-array(1,c(length(predictor),1))
  second<-array(phyloCorrectFactor.OUOUBM(c(alp,bet,xi,b1),predictor=predictor,x.tre=x.tre)*predictor)
  return(cbind(one,second))
  }

sim.dsm.oubmcir <- function(x , predictor_1 = predictor_1 , predictor_2 = predictor_2 , x.tre = x.tre){
  
  alp_1 <- x[1]
  alp_2 <- x[2]
  th_1 <- x[3]
  
  sig_2 <- x[4]
  b1 <- x[5]
  b2 <- x[6]
  n <- length(predictor_1)
  one <- array(1,c(n,1))
  second<-array(phyloCorrectFactor.OUBMCIR(c(alp_1,alp_2,th_1,sig_2,b1,b2),predictor_1 = predictor_1,predictor_2 = predictor_2,x.tre = x.tre)*predictor_1)
  third<-array(phyloCorrectFactor.OUBMCIR(c(alp_1,alp_2,th_1,sig_2,b1,b2),predictor_1 = predictor_1,predictor_2 = predictor_2,x.tre = x.tre)*predictor_2)
  return(cbind(one,second,third))
  }

sim.dsm.ououcir <- function(x , predictor_1 = predictor_1 , predictor_2 = predictor_2 , x.tre = x.tre){
  alp_1<-x[1]
  alp_2<-x[2]
  alp_3<-x[3]
  th_2<-x[4]
  
  sig_2<-x[5]
  b1<-x[6]
  b2<-x[7]
  p_est<-c(alp_1,alp_2,alp_3,th_2,sig_2,b1,b2)
  n <- length(predictor_1)
  one <- array(1,c(n,1))
  second <- array(phyloCorrectFactor.OUOUCIR(p_est , predictor_1 = predictor_1 , predictor_2 = predictor_2,x.tre = x.tre)*predictor_1)
  third <- array(phyloCorrectFactor.OUOUCIR(p_est , predictor_1 = predictor_1 , predictor_2 = predictor_2,x.tre = x.tre)*predictor_2)
  return(cbind(one,second,third))
}
