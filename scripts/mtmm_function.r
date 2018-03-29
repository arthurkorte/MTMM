
##############################################################################################################################################
###MTMM - Multi-Trait Mixed Model
###SET OF FUNCTIONS FOR GWAS ON CORRELATED TRAITS WITH CORRECTING FOR POPULATION STRUCTURE
#######
## 
##note: requires the R-packages  EMMA and AsREML 
## The AsREML package needs a valid license ( see http://www.vsni.co.uk/software/asreml)
##		
#library(emma)
#source('emma.r')
#library(lattice)
#library(asreml)
#Sys.setenv(ASREML_LICENSE_FILE='...')
#
##REQUIRED DATA & FORMAT
#
#PHENOTYPE - Y: a n by 3 matrix, where n=number of individuals and the first column contains the individual names, the second the phenotypic values for the first trait and the third column the phenotypic values of the second trait
#                  note, that the default script will only work, if in each row phenotypic measurements for both traits exists.
#GENOTYPE - X: a n by m matrix, where n=number of individuals, m=number of SNPs, with rownames(X)=individual names, and colnames(X)=SNP names in the format ('Chr- Pos'), e.g. '1- 567'..)
#KINSHIP - K: a n by n matrix, with rownames(K)=colnames(K)=individual names, you can calculate K as IBS matrix using the emma package K<-emma.kinship(t(X)) 
#each of these data being sorted in the same way, according to the individual name
#
##FOR PLOTING THE GWAS RESULTS

#
###### analysis will include single GWAS on both traits including the calculation of the respective effect sizes
###### might be slow for big datasets, user specified faster scans are possible
#
##GWAS MANHATTAN PLOTS
#
##QQPLOTS of pvalues

#

###

###  use default method for traits measured on different, but genetically identical individuals (e.g. A.thaliana data) or errorcorrelation, if both traits are measured on the same individual (e.g. for different human blood lipids) 
#method<-'default'
#method<-'errorcorrelation'
# run=FALSE will only report the correlation estimates and not perform the GWAS, the estimates are stored in a list called correlation
###
###   Methods for handling reps and estimates of non-overlapping sets are in progress and will be updated accordingly  

mtmm<-function(Y,X,K,method='default',include.single.analysis=T,calculate.effect.size=T,gen.data='binary',exclude=T,run=T)  {
## removing unique measurments

if(anyDuplicated(Y[,1])>0) { cat('duplicated ecotypes detected','\n')
stopifnot(anyDuplicated(Y[,1])==0)}


options(stringsAsFactors = FALSE)
SNP_INFO<-data.frame(cbind(colnames(X),matrix(nrow=ncol(X),ncol=2,data=unlist(strsplit(colnames(X),split='- ')),byrow=T)))
colnames(SNP_INFO)<-c('SNP','Chr','Pos')
SNP_INFO[,2]<-as.numeric(SNP_INFO[,2])
SNP_INFO[,3]<-as.numeric(SNP_INFO[,3])

if (exclude==T) {
#exclude non overlapping datasets
Y_<-na.omit(Y)
 
cat(nrow(Y)-nrow(Y_),'values excluded, leaving',nrow(Y_),'values','\n')
}else { Y_<-Y
cat('single anlysis performed on',length(na.omit(Y[,2])),'values for',colnames(Y)[2],'and ', length(na.omit(Y[,3])),'values for',colnames(Y)[3],'\n')

}





#stopifnot(nrow(Y_)==nrow(Y)) 

## preparing genotype and Kinship data 
X<-X[which(rownames(X)%in%Y_[,1]),]
Y<-Y_[which(Y_[,1]%in%rownames(X)),]
### enter here filtering step for 2 different kinships for different subsets .

Y1<-(Y[,2])
Y2<-(Y[,3])
names(Y1)<-Y[,1]
names(Y2)<-Y[,1]
Y1_<-na.omit(Y1)
Y2_<-na.omit(Y2)

ecot_id1<-as.integer(names(Y1_))
ecot_id2<-as.integer(names(Y2_))

K<-as.matrix(K[which(rownames(K)%in%names(Y1)),which(colnames(K)%in%names(Y1))])
K1<-as.matrix(K[which(rownames(K)%in%ecot_id1),which(colnames(K)%in%ecot_id1)])
K2<-as.matrix(K[which(rownames(K)%in%ecot_id2),which(colnames(K)%in%ecot_id2)])

n<-nrow(Y)
n1<-length(ecot_id1)
n2<-length(ecot_id2)


# combining the traits

Y_ok<-c(Y1,Y2)

#environment

Env<-c(rep(0,n),rep(1,n))

#standardize kinship 

K_stand<-(n-1)/sum((diag(n)-matrix(1,n,n)/n)*K)*K
K_inv<<-solve(K_stand)
K_stand1<-(n1-1)/sum((diag(n1)-matrix(1,n1,n1)/n1)*K1)*K1
K_stand2<-(n2-1)/sum((diag(n2)-matrix(1,n2,n2)/n2)*K2)*K2


#### need to test this section, if its still working !

if (gen.data=='binary') {
#calculate MAF&MAC Arabidopsis

AC_1 <- data.frame(colnames(X),apply(X,2,sum))
colnames(AC_1)<-c('SNP','AC_1')

MAF_1<-data.frame(AC_1,AC_0=nrow(X)-AC_1$AC_1)

MAF_2<-data.frame(MAF_1,MAC=apply(MAF_1[,2:3],1,min))

MAF_3<-data.frame(MAF_2,MAF=(MAF_2$MAC/nrow(X)))

MAF_ok<<-merge(SNP_INFO,MAF_3,by='SNP')

rm(AC_1,MAF_1,MAF_2,MAF_3)

#Filter for MAF

MAF<-subset(MAF_ok,MAF==0)[,1]

X_ok<-X[,!colnames(X) %in% MAF]

rm(MAF)
}else { if (gen.data=='heterozygot') {

count2<-function(x) {length(which(x==2))}
count1<-function(x) {length(which(x==1))}

AC_2 <- data.frame(colnames(X),apply(X,2,count2))
colnames(AC_2)<-c('SNP','AC_2')
AC_1 <- data.frame(colnames(X),apply(X,2,count1))
colnames(AC_1)<-c('SNP','AC_1')

MAF_1<-data.frame(AC_2,AC_1$AC_1,AC_0=nrow(X)-AC_1$AC_1-AC_2$AC_2)

MAF_2<-data.frame(MAF_1,MAC=apply(MAF_1[,c(2,4)],1,min))

MAF_3<-data.frame(MAF_2,MAF=(MAF_2$MAC/nrow(X)))

MAF_ok<<-merge(SNP_INFO,MAF_3,by='SNP')

rm(AC_1,MAF_1,MAF_2,MAF_3)

#Filter for MAF

MAF<-subset(MAF_ok,MAF==0)[,1]

X_ok<-X[,!colnames(X) %in% MAF]



} else {X_ok<-X}}
#############
#ANALYSIS
###
#####Y1



Xo1<-rep(1,n1)
ex1<-as.matrix(Xo1)
Xo2<-rep(1,n2)
ex2<-as.matrix(Xo2)

null1<<-emma.REMLE(Y1_,ex1,K_stand1)
herit1<-null1$vg/(null1$vg+null1$ve)
null2<<-emma.REMLE(Y2_,ex2,K_stand2)
herit2<-null2$vg/(null2$vg+null2$ve)


#### if single analysis is to include !:. ##
if (include.single.analysis==T & run==T) {

X1<-X_ok[rownames(X_ok)%in%names(Y1_),]

AC_1 <- data.frame(colnames(X1),apply(X1,2,sum))
colnames(AC_1)<-c('SNP','AC_1')

MAF_1<-data.frame(AC_1,AC_0=nrow(X)-AC_1$AC_1)

MAF_2<-data.frame(MAF_1,MAC=apply(MAF_1[,2:3],1,min))

MAF_3<-data.frame(MAF_2,MAF=(MAF_2$MAC/nrow(X)))

MAF_ok<<-merge(SNP_INFO,MAF_3,by='SNP')

rm(AC_1,MAF_1,MAF_2,MAF_3)

#Filter for MAF

MAF<-subset(MAF_ok,MAF==0)[,1]

X_ok1<-X1[,!colnames(X1) %in% MAF]

rm(MAF)


M1<-solve(chol(null1$vg*K_stand1+null1$ve*diag(dim(K_stand1)[1])))
Y_t1<-crossprod(M1,Y1_)
int_t1<-crossprod(M1,(rep(1,n1)))
X_t1<-crossprod(M1,X_ok1)
## two options : a with calculating the effect sizes of each SNP (slow)

if (calculate.effect.size==T) {
models1<-apply(X_t1,2,function(x){summary(lm(Y_t1~0+int_t1+x))$coeff[2,]})
out_models1<<-data.frame(SNP=colnames(models1),Pval_Y1=models1[4,],beta_Y1=models1[1,])
} else {
# or b : Fast scan without effect sizes

RSS_env<-rep(sum(lsfit(int_t1,Y_t1,intercept = FALSE)$residuals^2),ncol(X_ok1))
R1_full<-apply(X_ok1,2,function(x){sum(lsfit(cbind(int_t1,crossprod(M1,x)),Y_t1,intercept = FALSE)$residuals^2)})
m<-length(Y1)
F_1<-((RSS_env-R1_full)/1)/(R1_full/(m-3))
pval_Y1<-pf(F_1,1,(m-3),lower.tail=FALSE)
out_models1<-data.frame(SNP=names(pval_Y1),Pval_Y1=pval_Y1)}

X2<-X_ok[rownames(X_ok)%in%names(Y2_),]



AC_1 <- data.frame(colnames(X2),apply(X2,2,sum))
colnames(AC_1)<-c('SNP','AC_1')

MAF_1<-data.frame(AC_1,AC_0=nrow(X)-AC_1$AC_1)

MAF_2<-data.frame(MAF_1,MAC=apply(MAF_1[,2:3],1,min))

MAF_3<-data.frame(MAF_2,MAF=(MAF_2$MAC/nrow(X)))

MAF_ok<<-merge(SNP_INFO,MAF_3,by='SNP')

rm(AC_1,MAF_1,MAF_2,MAF_3)

#Filter for MAF

MAF<-subset(MAF_ok,MAF==0)[,1]

X_ok2<-X2[,!colnames(X2) %in% MAF]

rm(MAF)



M2<-solve(chol(null2$vg*K_stand2+null2$ve*diag(dim(K_stand2)[1])))
Y_t2<-crossprod(M2,Y2_)
int_t2<-crossprod(M2,(rep(1,n2)))
X_t2<-crossprod(M2,X_ok2)

if (calculate.effect.size==T) {
models2<-apply(X_t2,2,function(x){summary(lm(Y_t2~0+int_t2+x))$coeff[2,]})

out_models2<<-data.frame(SNP=colnames(models2),Pval_Y2=models2[4,],beta_Y2=models2[1,])

}else {
RSS_env<-rep(sum(lsfit(int_t2,Y_t2,intercept = FALSE)$residuals^2),ncol(X_ok2))
R2_full<-apply(X_ok2,2,function(x){sum(lsfit(cbind(int_t2,crossprod(M2,x)),Y_t2,intercept = FALSE)$residuals^2)})
m<-length(Y2)
F_2<-((RSS_env-R2_full)/1)/(R2_full/(m-3))
pval_Y2<-pf(F_2,1,(m-3),lower.tail=FALSE)
out_models2<-data.frame(SNP=names(pval_Y2),Pval_Y2=pval_Y2)}

} else { out_models1<-NA
out_models2<-NA}


### attention if exclude ==FALSE ncol(out_models1)<> ncol(out_models2) !


##### START  MTMM   #####

#intercept
Xo<-rep(1,2*n)

# preparing the data
data_asr<-data.frame(ecot_id=c(ecot_id1,ecot_id2),Env=Env,Y_ok=Y_ok)
data_asr$ecot_id<-as.factor(data_asr$ecot_id)
data_asr$Env<-as.factor(data_asr$Env)


## setting starting values

sigma1<<-null1$vg
sigma2<<-null2$vg
rho<-cor(Y1,Y2)/(sqrt(herit1)*sqrt(herit2))
###set borders for rho between -1 and 1
if (rho>0.98) {rho<-0.98
}else { if (rho< -0.98) {rho<--0.98}}


#alpha_e<-rho*sqrt((null1$vg+null1$ve)*(null2$vg+null2$ve))
alpha<<-rho*sqrt(null1$vg*null2$vg)
beta2<<-cor(Y1,Y2)*sqrt(null1$ve*null2$ve)
alpha2<<-cor(Y1,Y2)*sqrt(null1$vg*null2$vg)



if ( method == "default" ) {


 
mixmod.asr<-asreml(fixed=Y_ok~Env,random=~us(Env,init=c(sigma1,alpha,sigma2)):giv(ecot_id),ginverse=list(ecot_id=K_inv),rcov=~at(Env,init=c(null1$ve,null2$ve)):units,data=data_asr,maxit=30)
varcov<-matrix(data=c(summary(mixmod.asr)$varcomp$component[1:2],summary(mixmod.asr)$varcomp$component[2:3]),nrow=2,ncol=2)

ve<-matrix(data=c(summary(mixmod.asr)$varcomp$component[4],rep(0,2),summary(mixmod.asr)$varcomp$component[5]),nrow=2,ncol=2)




} else if ( method == "errorcorrelation" ) {
   
### including autocorrelation 



mixmod.asr<-asreml(fixed=Y_ok~Env,random=~us(Env,init=c(sigma1,alpha2,sigma2)):giv(ecot_id),ginverse=list(ecot_id=K_inv),rcov=~us(Env,init=c(null1$ve,beta2,null2$ve)):id(ecot_id),data=data_asr,maxit=30)

varcov<-matrix(data=c(summary(mixmod.asr)$varcomp$component[1:2],summary(mixmod.asr)$varcomp$component[2:3]),nrow=2,ncol=2)

ve<-matrix(data=c(summary(mixmod.asr)$varcomp$component[5:6],summary(mixmod.asr)$varcomp$component[6:7]),nrow=2,ncol=2)
}


rho_g<-varcov[1,2]/sqrt(varcov[1,1]*varcov[2,2])
rho_e<-ve[1,2]/sqrt(ve[1,1]*ve[2,2])
rho_p<-(varcov[1,2]+ve[1,2])/sqrt((varcov[1,1]+ve[1,1])*(varcov[2,2]+ve[2,2]))
pearson<-cor(Y1,Y2)
heriti1<-varcov[1,1]/(varcov[1,1]+ve[1,1])
heriti2<-varcov[2,2]/(varcov[2,2]+ve[2,2])
correlation<-list(pearson=pearson,gen_cor=rho_g,env_cor=rho_e,phen_cor=rho_p,h1_joint=heriti1,h2_joint=heriti2,converge=mixmod.asr$converge)

if ( run == FALSE) {
correlation<<-list(pearson=pearson,gen_cor=rho_g,env_cor=rho_e,phen_cor=rho_p,h1_joint=heriti1,h2_joint=heriti2,converge=mixmod.asr$converge)

print(unlist(correlation))
stop('No GWAS performed ! \n')}

#extract variances from asreml output

K_comb<-kronecker(varcov,K_stand)
rownames(K_comb)<-c(ecot_id1,ecot_id2)
colnames(K_comb)<-c(ecot_id1,ecot_id2)

I_comb<-kronecker(ve,diag(n))
rownames(I_comb)<-rownames(K_comb)
colnames(I_comb)<-colnames(K_comb)

bigK<-K_comb+I_comb


M<-solve(chol(bigK))


# scaling of the SNPs to interpret GxE results
X_ok1<-X_ok*sqrt(varcov[1,1])
X_ok2<-X_ok*sqrt(varcov[2,2])

Y_t<-crossprod(M,Y_ok)
cof_t<-crossprod(M,cbind(Xo,Env))

RSS_env<-sum(lsfit(cof_t,Y_t,intercept = FALSE)$residuals^2)

nbchunks<-5
m<-ncol(X_ok)

RSS_full<-list()
RSS_glob<-list()

for (j in 1:(nbchunks-1)){
X_<-rbind(X_ok1[,((j-1)*round(m/nbchunks)+1):(j*round(m/nbchunks))],X_ok2[,((j-1)*round(m/nbchunks)+1):(j*round(m/nbchunks))])
X_t<-array(dim=c(nrow(X_),ncol(X_),2))
X_t[,,1]<-crossprod(M,X_)
X_t[,,2]<-crossprod(M,(Env*X_))
RSS_full[[j]]<-apply(X_t,2,function(x){sum(lsfit(cbind(cof_t,x),Y_t,intercept=FALSE)$residuals^2)})
RSS_glob[[j]]<-apply(X_t[,,1],2,function(x){sum(lsfit(cbind(cof_t,x),Y_t,intercept=FALSE)$residuals^2)})
rm(X_,X_t)}
X_<-rbind(X_ok1[,((j)*round(m/nbchunks)+1):m],X_ok2[,((j)*round(m/nbchunks)+1):m])
X_t<-array(dim=c(nrow(X_),ncol(X_),2))
X_t[,,1]<-crossprod(M,X_)
X_t[,,2]<-crossprod(M,(Env*X_))
RSS_full[[nbchunks]]<-apply(X_t,2,function(x){sum(lsfit(cbind(cof_t,x),Y_t,intercept=FALSE)$residuals^2)})
RSS_glob[[nbchunks]]<-apply(X_t[,,1],2,function(x){sum(lsfit(cbind(cof_t,x),Y_t,intercept=FALSE)$residuals^2)})
rm(X_,X_t,j)

RSS_full<-unlist(RSS_full)
RSS_glob<-unlist(RSS_glob)

###############test###############
################################
###################################


#nb parameters in each models
#env
par_env<-ncol(cof_t)
par_glob<-par_env+1
par_full<-par_glob+1

##FTESTS
#FULL vs NULL

F_full<-(rep(RSS_env,m)/RSS_full-1)*(2*n-par_full)/(par_full-par_env)
F_ge<-(RSS_glob/RSS_full-1)*(2*n-par_full)/(par_full-par_glob)
F_glob<-(rep(RSS_env,m)/RSS_glob-1)*(2*n-par_glob)/(par_glob-par_env)

pval_full<-pf(F_full,par_full-par_env,2*n-par_full,lower.tail=FALSE)
pval_ge<-pf(F_ge,par_full-par_glob,2*n-par_full,lower.tail=FALSE)
pval_glob<-pf(F_glob,par_glob-par_env,2*n-par_glob,lower.tail=FALSE)

#outputs
outfull<<-data.frame('SNP'=colnames(X_ok),pval=pval_full)
outge<<-data.frame('SNP'=colnames(X_ok),pval=pval_ge)
outglob<<-data.frame('SNP'=colnames(X_ok),pval=pval_glob)

if (include.single.analysis==T & run==T) {
if (calculate.effect.size==T) {
## note to remove betas from matrix, if not calculated  
data.out<-merge(MAF_ok,cbind(out_models1[1:2],out_models2[,2],outfull[,2],outge[,2],outglob[,2],out_models1[,3],out_models2[,3]),by='SNP')
colnames(data.out)[8:14]<-c('pval_Y1','pval_Y2','pval_full','pval_trait_specific','pval_trait_common','beta_Y1','beta_Y2')
data.out_<-data.out[order(data.out[,3]),]
out<-data.out_[order(data.out_[,2]),]
} else { 
data.out<-merge(MAF_ok,cbind(out_models1[1:2],out_models2[,2],outfull[,2],outge[,2],outglob[,2]),by='SNP')
colnames(data.out)[8:12]<-c('pval_Y1','pval_Y2','pval_full','pval_trait_specific','pval_trait_common')
data.out_<-data.out[order(data.out[,3]),]
out<-data.out_[order(data.out_[,2]),]
}}else {

data.out<-merge(MAF_ok,cbind(outfull,outge[,2],outglob[,2]),by='SNP')
colnames(data.out)[8:10]<-c('pval_full','pval_trait_specific','pval_trait_common')
data.out_<-data.out[order(data.out[,3]),]
out<-data.out_[order(data.out_[,2]),]}

results<<-list(phenotype=Y,pvals=out,statistics=correlation,kinship=K)


}




#####      plots


 plot_gwas<-function(output,m=250000,h=8,maf=0.05,black=T) {
 colnames(output)[h]<-'Pval'

new_<-subset(output,output$MAF>maf&output$Pval<0.01)

vu<-ceiling(-log10(min(new_[,8]))+1)

output_<-new_[order(new_$Pos),]
output_ok<-output_[order(output_$Chr),]

maxpos<-c(0,cumsum(aggregate(output_ok$Pos,list(output_ok$Chr),max)$x+max(cumsum(aggregate(output_ok$Pos,list(output_ok$Chr),max)$x))/100))
if(black==T) {
plot_col<-rep(c('gray10','gray60'),ceiling(max(unique(output_ok$Chr))/2))
} else { plot_col<-c('blue','darkgreen','red','cyan','purple')}
size<-aggregate(output_ok$Pos,list(output_ok$Chr),length)$x

a<-rep(maxpos[1],size[1])
b<-rep(plot_col[1],size[1])
for (i in 2:max(unique(output_ok$Chr))){
a<-c(a,rep(maxpos[i],size[i]))
b<-c(b,rep(plot_col[i],size[i]))}

output_ok$xpos<-output_ok$Pos+a
output_ok$col<-b

d<-(aggregate(output_ok$xpos,list(output_ok$Chr),min)$x+aggregate(output_ok$xpos,list(output_ok$Chr),max)$x)/2

plot(output_ok$xpos,-log10(output_ok$Pval),col=output_ok$col,pch=16,ylab='-log10(pval)',xaxt='n',xlab='chromosome',axes=F,cex=1.2)
axis(1,tick=FALSE,at=d,labels=c(1:max(unique(output_ok$Chr))))
axis(2,lwd=2)

abline(h=-log10(0.05/m),lty=3,col='black',lwd=2)
}





#qqplots

qq_plot<-function(output,h=8,maf=0.05,farbe='red') {
colnames(output)[h]<-'Pval'

out<-subset(output,output$MAF>maf&output$Pval<0.05)
e<--log10(ppoints(nrow(output)))[1:nrow(out)]
o<--log10(sort(out$Pval))
maxp<-ceiling(max(unlist(o)))+1


plot(e,o,type='l',cex=0.8,col=farbe,xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))),xlim=c(0,max(e)+1),ylim=c(0,maxp),lwd=2,axes=FALSE)
abline(0,1,col="dark grey")
axis(1,lwd=2)
axis(2,lwd=2)
}

qq_plot_all<-function(output,h=c(8,9,10,11,12),maf=0.05,farbe=c('blue','cornflowerblue','goldenrod','chartreuse','brown2')) {
maxp<-ceiling(max(-log10(output[,h])))

for (z in 1: length(h)) {
colnames(output)[h[[z]]]<-'Pval'

out<-subset(output,output$MAF>maf&output[,h[[z]]]<0.05)
e<--log10(ppoints(nrow(output)))[1:nrow(out)]
o<--log10(sort(out[,h[[z]]]))



plot(e,o,type='l',cex=0.8,col=farbe[[z]],xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))),xlim=c(0,max(e)+1),ylim=c(0,maxp),lwd=2,axes=FALSE)
abline(0,1,col="dark grey")
axis(1,lwd=2)
axis(2,lwd=2)
par(new=T)}
}







