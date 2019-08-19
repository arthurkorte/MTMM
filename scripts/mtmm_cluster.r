mtmm_cluster<-function(out.name,save.output=T,generate.plot=T,X.folder='/storage/evolgen/data/2029_data',incl.singleGWAS=FALSE) {

for ( u in 1:22) {

filename<-paste(X.folder,'/X_2029_',u,'.rda',sep='')


load(filename)
if (u==1) {
results<-mtmm_part2(X,incl.singleGWAS=incl.singleGWAS)

rm(X)
} else {

out<-mtmm_part2(X,incl.singleGWAS=incl.singleGWAS)
results<-rbind(results,out)
rm(X,out)
}

}



if (save.output==T) {

name1<-paste(out.name,'.rda',sep='')
save(results,correlation,file=name1)
}
if (generate.plot==T) {
name2<-paste(out.name,'.pdf',sep='')
plot_mtmm(results,correlation,name2,incl.singleGWAS=incl.singleGWAS)
}

}
