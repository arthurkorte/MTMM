

### function for estimating the MTMM parameters using ASREML, to get passed on to a non-ASREML programm !!
#library(lattice)
#library(asreml)
#library(msm)
#library(nadiv)
#source('emma.r')

mtmm_estimates<-function(Y,k=2,l=3,K,method='default',only.vca=FALSE)  {

name<-paste(colnames(Y)[k],colnames(Y)[l],'mtmm_estimates.rda',sep='_')
Y<-Y[,c(1,k,l)]
rm(k,l)


if(anyDuplicated(Y[,1])>0) { cat('duplicated ecotypes detected','\n')
stopifnot(anyDuplicated(Y[,1])==0)}



Y_<-na.omit(Y)


Y<-Y_[which(Y_[,1]%in%rownames(K)),]
cat('analysis performed on',nrow(Y),'ecotypes',sep=' ')


Y1<-(Y[,2])
Y2<-(Y[,3])
names(Y1)<-Y[,1]
names(Y2)<-Y[,1]

# ecot_id<-as.integer(names(Y1))
ecot_id<-names(Y1)

K<-as.matrix(K[which(rownames(K)%in%names(Y1)),which(colnames(K)%in%names(Y1))])


n<-nrow(Y)




# combining the traits

Y_ok<-c(Y1,Y2)

#environment

Env<-c(rep(0,n),rep(1,n))

#standardize kinship 

K_stand<-(n-1)/sum((diag(n)-matrix(1,n,n)/n)*K)*K

K_inv<<-solve(K_stand)

#############
#ANALYSIS
###
#####Y1



Xo<-rep(1,n)
ex<-as.matrix(Xo)


null1<<-emma.REMLE(Y1,ex,K_stand)
herit1<<-null1$vg/(null1$vg+null1$ve)
null2<<-emma.REMLE(Y2,ex,K_stand)
herit2<<-null2$vg/(null2$vg+null2$ve)


#### if single analysis is to include !:. ##



### attention if exclude ==FALSE ncol(out_models1)<> ncol(out_models2) !


##### START  MTMM   #####

#intercept
Xo<-rep(1,2*n)

# preparing the data
data_asr<-data.frame(ecot_id=c(ecot_id,ecot_id),Env=Env,Y_ok=Y_ok)
data_asr$ecot_id<-as.factor(data_asr$ecot_id)
data_asr$Env<-as.factor(data_asr$Env)


## setting starting values

sigma1<<-null1$vg
sigma2<<-null2$vg
rho<<-cor(Y1,Y2)/(sqrt(herit1)*sqrt(herit2))
###set borders for rho between -1 and 1
if (rho>0.98) {rho<-0.98
}else { if (rho< -0.98) {rho<--0.98}}


#alpha_e<-rho*sqrt((null1$vg+null1$ve)*(null2$vg+null2$ve))
alpha<<-rho*sqrt(null1$vg*null2$vg)
beta2<<-cor(Y1,Y2)*sqrt(null1$ve*null2$ve)
alpha2<<-cor(Y1,Y2)*sqrt(null1$vg*null2$vg)



if ( method == "default" ) {


 
mixmod.asr<-asreml(fixed=Y_ok~Env,random=~us(Env,init=c(sigma1,alpha,sigma2)):giv(ecot_id),ginverse=list(ecot_id=K_inv),rcov=~at(Env,init=c(null1$ve,null2$ve)):units,data=data_asr,maxit=50)
varcov<-matrix(data=c(summary(mixmod.asr)$varcomp$component[1:2],summary(mixmod.asr)$varcomp$component[2:3]),nrow=2,ncol=2)

ve<-matrix(data=c(summary(mixmod.asr)$varcomp$component[4],rep(0,2),summary(mixmod.asr)$varcomp$component[5]),nrow=2,ncol=2)




} else if ( method == "errorcorrelation" ) {
   
### including autocorrelation 


mixmod.asr<-asreml(fixed=Y_ok~Env,random=~us(Env,init=c(sigma1,alpha2,sigma2)):giv(ecot_id),ginverse=list(ecot_id=K_inv),rcov=~us(Env,init=c(null1$ve,beta2,null2$ve)):id(ecot_id),data=data_asr,maxit=100)

varcov<-matrix(data=c(summary(mixmod.asr)$varcomp$component[1:2],summary(mixmod.asr)$varcomp$component[2:3]),nrow=2,ncol=2)

ve<-matrix(data=c(summary(mixmod.asr)$varcomp$component[5:6],summary(mixmod.asr)$varcomp$component[6:7]),nrow=2,ncol=2)
}


rho_g<-varcov[1,2]/sqrt(varcov[1,1]*varcov[2,2])
rho_e<-ve[1,2]/sqrt(ve[1,1]*ve[2,2])
rho_p<-(varcov[1,2]+ve[1,2])/sqrt((varcov[1,1]+ve[1,1])*(varcov[2,2]+ve[2,2]))
pearson<-cor(Y1,Y2)
heriti1<-varcov[1,1]/(varcov[1,1]+ve[1,1])
heriti2<-varcov[2,2]/(varcov[2,2]+ve[2,2])
G<-varcov[1,2]
G1<-varcov[1,1]-varcov[1,2]
G2<-varcov[2,2]-varcov[1,2]
E<-ve[1,1]+ve[2,2]
EE<-ve[1,2]
tot_var<-G+G1+G2+E+EE
var_G<-G/tot_var
var_GE<-(G1+G2)/tot_var
var_E<-E/tot_var
var_EE<-EE/tot_var
hallo<-mixmod.asr
if ( method == "default" ) {


werte<-summary(hallo)$varcomp$component [1:3]
werte2<-summary(hallo)$varcomp$component [4:5]

covG<-aiFun(hallo)[1:3,1:3]

covH1<-matrix(nrow=2,ncol=2,data=c(aiFun(hallo)[1,1],aiFun(hallo)[4,1],aiFun(hallo)[4,1],aiFun(hallo)[4,4]))
covH2<-matrix(nrow=2,ncol=2,data=c(aiFun(hallo)[3,3],aiFun(hallo)[5,3],aiFun(hallo)[5,3],aiFun(hallo)[5,5]))
covP<-aiFun(hallo)
covG[upper.tri(covG)==T]<-covG[lower.tri(covG)==T]

covP[upper.tri(covP)==T]<-covP[lower.tri(covP)==T]
null.g<-asreml(fixed=Y_ok~Env,random=~idh(Env,init=c(sigma1,sigma2)):giv(ecot_id),ginverse=list(ecot_id=K_inv),rcov=~at(Env,init=c(null1$ve,null2$ve)):units,data=data_asr,maxit=50)


ap<-hallo$loglik

cp<-null.g$loglik
se_g<-deltamethod( ~ x2/sqrt(x1*x3),werte,covG)
p_g<-1-pchisq(2*(ap-cp),1)
se_e<-NA
p_e<-NA
se_p<-deltamethod( ~ (x2)/sqrt((x1+x4)*(x3+x5)), c(werte,werte2),covP)
se_h1<-deltamethod( ~ x1/(x1+x2),c(werte[1],werte2[1]),covH1)
se_h2<-deltamethod( ~ x1/(x1+x2),c(werte[3],werte2[2]),covH2)



} else if ( method == "errorcorrelation" ) {

werte<-summary(hallo)$varcomp$component [1:3]
werte2<-summary(hallo)$varcomp$component [5:7]

covG<-aiFun(hallo)[1:3,1:3]
covE<-aiFun(hallo)[5:7,5:7]
covH1<-matrix(nrow=2,ncol=2,data=c(aiFun(hallo)[1,1],aiFun(hallo)[5,1],aiFun(hallo)[5,1],aiFun(hallo)[5,5]))
covH2<-matrix(nrow=2,ncol=2,data=c(aiFun(hallo)[3,3],aiFun(hallo)[7,3],aiFun(hallo)[7,3],aiFun(hallo)[7,7]))
covP<-aiFun(hallo)[-4,-4]
covG[upper.tri(covG)==T]<-covG[lower.tri(covG)==T]
covE[upper.tri(covE)==T]<-covE[lower.tri(covE)==T]
covP[upper.tri(covP)==T]<-covP[lower.tri(covP)==T]
null.e<-asreml(fixed=Y_ok~Env,random=~us(Env,init=c(sigma1,alpha2,sigma2)):giv(ecot_id),ginverse=list(ecot_id=K_inv),rcov=~idh(Env,init=c(null1$ve,null2$ve)):id(ecot_id),data=data_asr,maxit=50)
null.g<-asreml(fixed=Y_ok~Env,random=~idh(Env,init=c(sigma1,sigma2)):giv(ecot_id),ginverse=list(ecot_id=K_inv),rcov=~us(Env,init=c(null1$ve,beta2,null2$ve)):id(ecot_id),data=data_asr,maxit=50)


ap<-hallo$loglik
bp<-null.e$loglik
cp<-null.g$loglik
se_g<-deltamethod( ~ x2/sqrt(x1*x3),werte,covG)
p_g<-1-pchisq(2*(ap-cp),1)
se_e<-deltamethod( ~ x2/sqrt(x1*x3),werte2,covE)
p_e<-1-pchisq(2*(ap-bp),1)
se_p<-deltamethod( ~ (x2+x5)/sqrt((x1+x4)*(x3+x6)), c(werte,werte2),covP)
se_h1<-deltamethod( ~ x1/(x1+x2),c(werte[1],werte2[1]),covH1)
se_h2<-deltamethod( ~ x1/(x1+x2),c(werte[3],werte2[3]),covH2)

}
correlation<-list(converge=hallo$converge,pearson=pearson,gen_cor=rho_g,SE_gen_cor=se_g,p_gen_cor=p_g,env_cor=rho_e,SE_env_cor=se_e,p_env_cor=p_e,phen_cor=rho_p,SE_phen_cor=se_p,h1_joint=heriti1,SE_h1_joint=se_h1,h2_joint=heriti2,SE_h2_joint=se_h2,var_G=var_G,var_GE=var_GE,var_E=var_E,var_EE=var_EE)






#extract variances from asreml output

K_comb<-kronecker(varcov,K_stand)
rownames(K_comb)<-c(ecot_id,ecot_id)
colnames(K_comb)<-c(ecot_id,ecot_id)

I_comb<-kronecker(ve,diag(n))
rownames(I_comb)<-rownames(K_comb)
colnames(I_comb)<-colnames(K_comb)


### what is this step doing ???   
bigK<-K_comb+I_comb


M<-solve(chol(bigK))


# scaling of the SNPs to interpret GxE results


Y_t<-crossprod(M,Y_ok)
cof_t<-crossprod(M,cbind(Xo,Env))



if(only.vca==T) {return(unlist(correlation))
} else {

save(n,cof_t,M,Env,correlation,varcov,Y_t,ecot_id,Xo,K_stand,ecot_id,file=name)

}}



