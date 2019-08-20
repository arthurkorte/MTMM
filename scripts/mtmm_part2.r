
mtmm_part2<-function(X,incl.singleGWAS=FALSE,use.SNP_INFO=FALSE,mac=5) {

X_<-X[rownames(X)%in%ecot_id,]
rm(X)

AC <- data.frame(colnames(X_),apply(X_,2,function(x){length(which(x==0))}),apply(X_,2,function(x){length(which(x==0.5))}),apply(X_,2,function(x){length(which(x==1))}))
colnames(AC)<-c('SNP','AC_0','AC_0.5','AC_1')

MAC<-data.frame(AC,MAC=apply(AC[,c(2,4)],1,min)+0.5*AC$AC_0.5)

MAF<-data.frame(MAC,MAF=(MAC$MAC/nrow(X_)))


MAF_ok<-subset(MAF,MAC>mac)

X_ok<-X_[,colnames(X_) %in% MAF_ok[,1]]

rm(AC,MAC,MAF,X_)


# scaling of the SNPs to interpret GxE results
X_ok1<-X_ok*sqrt(varcov[1,1])
X_ok2<-X_ok*sqrt(varcov[2,2])

m<-ncol(X_ok)

### now the calculation with summary(lm) to get effect sizes for glob and gxe pval. not working for the full p-value

X_<-rbind(X_ok1,X_ok2)

X_t<-array(dim=c(nrow(X_),ncol(X_),2))
  X_t[,,1]<-crossprod(M,X_)
  X_t[,,2]<-crossprod(M,(Env*X_))
  
  glob <-apply(X_t,2,function(x){summary(lm(Y_t~0+cof_t+x[,1]))$coeff[3,]})
  ge <-apply(X_t,2,function(x){summary(lm(Y_t~0+cof_t+x[,1]+x[,2]))$coeff[4,]})
  
  ## idea for getting effect size of the full model, but this is not working at the moment 
  
  #full<-apply(X_t,2,function(x){a<-summary(lm(Y_t~0+cof_t+x[,1]+x[,2]))$fstatistic
  #1-pf(a[1],a[2],a[3])})
  
  
  Res_<-as.data.frame(cbind(glob[4,],glob[1,],glob[2,],ge[4,],ge[1,],ge[2,]))
 ##
  Res<-cbind(colnames(X_),Res_)
 colnames(Res)<-c('SNP','pval_glob','beta_glob','se_glob','pval_ge','beta_ge','se_ge')
 
 RSS_full<-apply(X_t,2,function(x){sum(lsfit(cbind(cof_t,x),Y_t,intercept=FALSE)$residuals^2)})
 #RSS_glob<-apply(X_t[,,1],2,function(x){sum(lsfit(cbind(cof_t,x),Y_t,intercept=FALSE)$residuals^2)})
 RSS_env<-sum(lsfit(cof_t,Y_t,intercept = FALSE)$residuals^2)


 par_env<-ncol(cof_t)
 #par_glob<-par_env+1
 par_full<-par_env+2
 
 ##FTESTS

 F_full<-(rep(RSS_env,m)/RSS_full-1)*(2*n-par_full)/(par_full-par_env)
 #F_ge<-(RSS_glob/RSS_full-1)*(2*n-par_full)/(par_full-par_glob)
 #F_glob<-(rep(RSS_env,m)/RSS_glob-1)*(2*n-par_glob)/(par_glob-par_env)
 
 pval_full<-pf(F_full,par_full-par_env,2*n-par_full,lower.tail=FALSE)
 #pval_ge<-pf(F_ge,par_full-par_glob,2*n-par_full,lower.tail=FALSE)
 #pval_glob<-pf(F_glob,par_glob-par_env,2*n-par_glob,lower.tail=FALSE)
 
results<-cbind(Res,pval_full)

## generate SNP_INFO for Arabidopsis data, alternativley provide your own file linking SNP names to chr and pos.
if(use.SNP_INFO == FALSE) { 
  options(stringsAsFactors = FALSE)
  SNP_INFO<-data.frame(cbind(colnames(X_ok),matrix(nrow=ncol(X_ok),ncol=2,data=unlist(strsplit(colnames(X_ok),split='- ')),byrow=T)))
  colnames(SNP_INFO)<-c('SNP','Chr','Pos')
  SNP_INFO[,2]<-as.numeric(SNP_INFO[,2])
  SNP_INFO[,3]<-as.numeric(SNP_INFO[,3])
}

data.1<-merge(SNP_INFO,MAF_ok,by='SNP')
data.out<-merge(data.1,results,by='SNP')

## und jetzt single GWAS 
if (incl.singleGWAS==T) {
  out1<-amm_gwas(Y,m=2,X=X_ok,K=K_stand,report=FALSE)
  colnames(out1)[8]<-'pval_Y1'
  out2<-amm_gwas(Y,m=3,X=X_ok,K=K_stand,report=FALSE)
  colnames(out2)[8]<-'pval_Y2'
  data.out<-merge(data.out,out1[,c(1,8)],by='SNP')
  data.out<-merge(data.out,out2[,c(1,8)],by='SNP')
}

#outputs

data.out_<-data.out[order(data.out[,3]),]
out<-data.out_[order(data.out_[,2]),]
return(out)
}
