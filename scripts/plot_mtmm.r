plot_mtmm<-function(results,correlation,name2,incl.singleGWAS=T) {

pie<-data.frame(
  class=c('G','GxE','E'),
  perc=unlist(correlation)[17:19])

pie<- pie %>%
  arrange(desc(class)) %>%
  mutate(lab.ypos = cumsum(perc) - 0.5*perc)
pie[,2]<-round(pie[,2],digits=2)
mycols <- c("#696969", "#C0C0C0","#2E8B57")

vca<-ggplot(pie, aes(x = 2, y = perc, fill = class)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0)+
  geom_text(aes(y = lab.ypos, label = perc), color = "white")+
  scale_fill_manual(values = mycols) +
  theme_void()+
  xlim(0.5, 2.5)+
  ggtitle(paste('VCA',colnames(Y)[2],colnames(Y)[3],sep='_'))

pdf(file=name2)
print(vca)
plot_gwas(results,h=9,maf=0.05,name='MTMM_common_effect')
plot_gwas(results,h=12,maf=0.05,name='MTMM_GxE_effect')
plot_gwas(results,h=15,maf=0.05,name='MTMM_full_model')
if (incl.singleGWAS==T) {
  plot_gwas(results,h=16,maf=0.05,name=paste('single GWAS',colnames(Y)[2]))
  plot_gwas(results,h=17,maf=0.05,name=paste('single GWAS',colnames(Y)[3]))
}
dev.off()
}