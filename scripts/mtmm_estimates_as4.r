### function for estimating the MTMM parameters using ASREML, to get passed on to a non-ASREML programm !!
#library(lattice)
#library(asreml)
#library(msm)
#library(nadiv)
#source('/storage/evolgen/R/scripts/emma.r')

mtmm_estimates<-function(Y,k=2,l=3,K,method='default',only.vca=FALSE)  {

name<-paste(colnames(Y)[k],colnames(Y)[l],'mtmm_estimates.rda',sep='_')
Y<-Y[,c(1,k,l)]
rm(k,l)

if(anyDuplicated(Y[,1])>0) { cat('duplicated ecotypes detected','\n')
stopifnot(anyDuplicated(Y[,1])==0)}
if((is.numeric(Y[,1])|is.integer(Y[,1]))==FALSE) stop('ecot_id is not numeric or integer','\n')




Y<-Y[order(Y[,1]),]
Y_<-na.omit(Y)
Y<-Y_[which(Y_[,1]%in%rownames(K)),]
cat('analysis performed on',nrow(Y),'ecotypes',sep=' ')

Y1<-(Y[,2])
Y2<-(Y[,3])
names(Y1)<-Y[,1]
names(Y2)<-Y[,1]
K<-as.matrix(K[which(rownames(K)%in%names(Y1)),which(colnames(K)%in%names(Y1))])
ecot_id<-as.integer(Y[,1])
Y[,1]<-ecot_id
n<-nrow(Y)
colnames(K)=rownames(K)=ecot_id
# combining the traits

Y_ok<-c(Y1,Y2)

#environment

Env<-c(rep(0,n),rep(1,n))

#standardize kinship 

K_stand<-(n-1)/sum((diag(n)-matrix(1,n,n)/n)*K)*K
Xo<-rep(1,n)
ex<-as.matrix(Xo)

null1<-emma.REMLE(Y1,ex,K_stand)
herit1<-null1$vg/(null1$vg+null1$ve)
null2<-emma.REMLE(Y2,ex,K_stand)
herit2<-null2$vg/(null2$vg+null2$ve)

##### START  MTMM   #####

#intercept
Xo<-rep(1,2*n)

# preparing the data
data_asr<-data.frame(ecot_id=c(ecot_id,ecot_id),Env=Env,Y_ok=Y_ok)
data_asr$ecot_id<-as.factor(data_asr$ecot_id)
data_asr$Env<-as.factor(data_asr$Env)
# bind K_stand to global env, that asreml can see it 
assign('K_stand',K_stand,'.GlobalEnv')
## setting starting values

sigma1<-null1$vg
sigma2<-null2$vg
rho<-cor(Y1,Y2)/(sqrt(herit1)*sqrt(herit2))
###set borders for rho between -1 and 1
if (rho>0.98) {rho<-0.98
}else { if (rho< -0.98) {rho<--0.98}}


alpha<-rho*sqrt(null1$vg*null2$vg)
beta2<-cor(Y1,Y2)*sqrt(null1$ve*null2$ve)
alpha2<-cor(Y1,Y2)*sqrt(null1$vg*null2$vg)

if ( method == "default" ) {

mixmod.asr<-asreml(fixed=Y_ok~Env,random=~us(Env):vm(ecot_id,K_stand),residual=~dsum(~id(ecot_id)|Env),data=data_asr,maxit=30)
varcov<-matrix(data=c(summary(mixmod.asr)$varcomp$component[1:2],summary(mixmod.asr)$varcomp$component[2:3]),nrow=2,ncol=2)
ve<-matrix(data=c(summary(mixmod.asr)$varcomp$component[4],rep(0,2),summary(mixmod.asr)$varcomp$component[5]),nrow=2,ncol=2)

} else if ( method == "errorcorrelation" ) {
   
mixmod.asr<-asreml(fixed=Y_ok~Env,random=~us(Env):vm(ecot_id,K_stand),residual=~us(Env):id(ecot_id),data=data_asr,maxit=30)
varcov<-matrix(data=c(summary(mixmod.asr)$varcomp$component[1:2],summary(mixmod.asr)$varcomp$component[2:3]),nrow=2,ncol=2)
ve<-matrix(data=c(summary(mixmod.asr)$varcomp$component[5:6],summary(mixmod.asr)$varcomp$component[6:7]),nrow=2,ncol=2)
}
# scale environmental amnd genetic correlation 

rho_g<-(varcov[1,2]/sqrt(varcov[1,1]*varcov[2,2]))/2

rho_e<-(ve[1,2]/sqrt(ve[1,1]*ve[2,2]))/2

#rho_p<-(varcov[1,2]+ve[1,2])/sqrt((varcov[1,1]+ve[1,1])*(varcov[2,2]+ve[2,2]))
rho_p<-rho_g+rho_e
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
### aiFun() is currently not working use mixmod.asr$ai instead   aiFun(hallo) = hallo$ai
if ( method == "default" ) {

werte<-summary(hallo)$varcomp$component [1:3]
werte2<-summary(hallo)$varcomp$component [4:5]

covG<-hallo$ai[1:3,1:3]

covH1<-matrix(nrow=2,ncol=2,data=c(hallo$ai[1,1],hallo$ai[4,1],hallo$ai[4,1],hallo$ai[4,4]))
covH2<-matrix(nrow=2,ncol=2,data=c(hallo$ai[3,3],hallo$ai[5,3],hallo$ai[5,3],hallo$ai[5,5]))
covP<-hallo$ai
#covG[upper.tri(covG)==T]<-covG[lower.tri(covG)==T]

#covP[upper.tri(covP)==T]<-covP[lower.tri(covP)==T]

null.g<-asreml(fixed=Y_ok~Env,random=~idh(Env):vm(ecot_id,K_stand),residual=~dsum(~id(ecot_id) | Env),data=data_asr,maxit=30)
ap<-hallo$loglik
cp<-null.g$loglik

se_g<-deltamethod( ~ (x2/sqrt(x1*x3))/2,werte,covG)
p_g<-1-pchisq(2*(ap-cp),1)
se_e<-NA
p_e<-NA
#se_p<-deltamethod( ~ (x2)/sqrt((x1+x4)*(x3+x5)), c(werte,werte2),covP)
se_p<-se_g
se_h1<-deltamethod( ~ x1/(x1+x2),c(werte[1],werte2[1]),covH1)
se_h2<-deltamethod( ~ x1/(x1+x2),c(werte[3],werte2[2]),covH2)

} else if ( method == "errorcorrelation" ) {

werte<-summary(hallo)$varcomp$component [1:3]
werte2<-summary(hallo)$varcomp$component [5:7]

covG<-hallo$ai[1:3,1:3]
covE<-hallo$ai[5:7,5:7]
covH1<-matrix(nrow=2,ncol=2,data=c(hallo$ai[1,1],hallo$ai[5,1],hallo$ai[5,1],hallo$ai[5,5]))
covH2<-matrix(nrow=2,ncol=2,data=c(hallo$ai[3,3],hallo$ai[7,3],hallo$ai[7,3],hallo$ai[7,7]))
covP<-hallo$ai[-4,-4]
#covG[upper.tri(covG)==T]<-covG[lower.tri(covG)==T]
#covE[upper.tri(covE)==T]<-covE[lower.tri(covE)==T]
#covP[upper.tri(covP)==T]<-covP[lower.tri(covP)==T]

#both calls are identical, number of vcs differs ! 
#null.e1<-asreml(fixed=Y_ok~Env,random=~us(Env):vm(ecot_id,K_stand),residual=~idh(Env):id(ecot_id),data=data_asr,maxit=30)
null.e<-asreml(fixed=Y_ok~Env,random=~us(Env):vm(ecot_id,K_stand),residual=~dsum(~id(ecot_id)|Env),data=data_asr,maxit=30)

null.g<-asreml(fixed=Y_ok~Env,random=~idh(Env):vm(ecot_id,K_stand),residual=~us(Env):id(ecot_id),data=data_asr,maxit=30)
ap<-hallo$loglik
bp<-null.e$loglik
cp<-null.g$loglik
se_g<-deltamethod( ~ (x2/sqrt(x1*x3)+x2)/2,werte,covG)
p_g<-1-pchisq(2*(ap-cp),1)
se_e<-deltamethod( ~ (x2/sqrt(x1*x3)+x2)/2,werte2,covE)
p_e<-1-pchisq(2*(ap-bp),1)
#se_p<-deltamethod( ~ (x2+x5)/sqrt((x1+x4)*(x3+x6)), c(werte,werte2),covP)
se_p<-se_g+se_e
se_h1<-deltamethod( ~ x1/(x1+x2),c(werte[1],werte2[1]),covH1)
se_h2<-deltamethod( ~ x1/(x1+x2),c(werte[3],werte2[3]),covH2)

}

# need to update what to output 
correlation<-list(converge=hallo$converge,pearson=pearson,gen_cor=rho_g,se_gen_cor=se_g,p_gen_cor=p_g,env_cor=rho_e,se_env_cor=se_e,p_env_cor=p_e,phen_cor=rho_p,se_phen_cor=se_p,h1=herit1,h1_j=heriti1,se_h1_j=se_h1,h2=herit2,h2_j=heriti2,se_h2_j=se_h2,var_G=var_G,var_GE=var_GE,var_E=var_E,var_EE=var_EE)

#extract variances from asreml output

K_comb<-kronecker(varcov,K_stand)
rownames(K_comb)<-c(ecot_id,ecot_id)
colnames(K_comb)<-c(ecot_id,ecot_id)

I_comb<-kronecker(ve,diag(n))
rownames(I_comb)<-rownames(K_comb)
colnames(I_comb)<-colnames(K_comb)

bigK<-K_comb+I_comb

M<-solve(chol(bigK))

Y_t<-crossprod(M,Y_ok)
cof_t<-crossprod(M,cbind(Xo,Env))

if(only.vca==T) {return(unlist(correlation))
} else {

save(correlation,Y,n,cof_t,M,Env,varcov,Y_t,ecot_id,Xo,K_stand,ecot_id,file=name)

}}



