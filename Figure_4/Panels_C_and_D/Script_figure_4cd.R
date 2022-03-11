library(plotrix)

##Figure 4D##

## Import data
plate <- read.csv("20210411_dox_HCT116_miR34a_miR21_for_analysis.csv", header = T, sep = ';',dec=",",stringsAsFactors = F)
# The control replicate n°3 collected 72h after doxorubicin treatment (ctrl-72h-r3) was not amplified correctly for both miR-34a and miR-21.
# We started over the reverse-transcription and ddPCR for this sample and applied new measurements to complete previous data.
plate_bis <- read.csv("20210416_dox_HCT116_miR34a_miR21_for_analysis.csv", header = T, sep = ';',dec=",",stringsAsFactors = F)
sample_name=c(rep(c("ctrl-24h-r1","ctrl-24h-r2","ctrl-24h-r3","ctrl-72h-r1","ctrl-72h-r2","ctrl-72h-r3",
                    "dox-24h-r1","dox-24h-r2","dox-24h-r3","dox-72h-r1","dox-72h-r2","dox-72h-r3"),2))
target_miRNA=c(rep("miR-34a",12),rep("miR-21",12))
coeff_dilution=c(rep(1,12),rep(10,12)) #Reverse-transcription reactions were diluted to 1/10 for miR-21 quantification.
ng_RT=rep(10,24) #Reverse-transcription reactions were performed with 10 ng total RNA in final volume of 15 µL.
volume_RT=15
volume_RT_for_ddPCR=1.33 #ddPCR quantifications were performed from 1.33 µL reverse-transcription reactions.

table=as.data.frame(cbind(plate$Well,sample_name,target_miRNA,coeff_dilution,ng_RT,plate$Copies.20µLWell),stringsAsFactors = F)
colnames(table)=c("well","sample","target_miRNA","dilution","ng_RT","copies_per_well")
#As said before, to replace NA values in the first plate (control replicate n°3 for miR-34a and miR-21), we applied the mean of 4 technical replicates from the second plate.
table$copies_per_well[which(table$sample == "ctrl-72h-r3" & table$target_miRNA == "miR-34a")]=mean(plate_bis$Copies.20µLWell[1:4])
table$copies_per_well[which(table$sample == "ctrl-72h-r3" & table$target_miRNA == "miR-21")]=mean(plate_bis$Copies.20µLWell[5:8])
table$copies_per_ng = as.numeric(table$copies_per_well)/(((volume_RT_for_ddPCR*as.numeric(table$ng_RT))/volume_RT)/as.numeric(table$dilution))

#miR-34a and miR-21 levels in biological replicates
#par(mfrow = c(1, 2))
#barplot(table$copies_per_ng[which(table$target_miRNA == "miR-34a")], main ='miR-34a abundance\nper sample',names.arg =table$sample[which(table$target_miRNA == "miR-34a")] ,las =2 , ylab = "Copies of miRNA / ng of total RNA", cex.names = 0.75)
#barplot(table$copies_per_ng[which(table$target_miRNA == "miR-21")], main ='miR-21 abundance\nper sample',names.arg =table$sample[which(table$target_miRNA == "miR-21")] ,las =2, ylab = "Copies of miRNA / ng of total RNA", cex.names = 0.75)

mean_table=c()
pas=seq(1,22,3)
rep=rep(2,length(pas))
for (i in 1:length(pas)) {
  mean_table=rbind(mean_table,c(paste(table$sample[pas[i]],table$target_miRNA[pas[i]],sep="_"), mean(table$copies_per_ng[pas[i]:(pas[i]+rep[i])],na.rm=T),sd(table$copies_per_ng[pas[i]:(pas[i]+rep[i])],na.rm=T)))
}
mean_table=as.data.frame(mean_table)
mean_table$V2=as.numeric(as.character(mean_table$V2))
mean_table$V3=as.numeric(as.character(mean_table$V3))
#pdf(file="RTddPCR_dox_miRNAs_per_ng_split.pdf", width= 5.5, height = 4.5)
par(mfrow = c(1, 2),mar= c(5.1, 4.1, 4.1, 2.1))
plotCI(barplot(mean_table$V2[5:8],ylim=c(0,max(mean_table$V2[5:8]+mean_table$V3[5:8])),main ='miR-21 abundance per group', ylab="miR-21 copies per ng of total RNA",col=rep("dimgrey",4), names.arg = c("ctrl_24h","ctrl_72h","dox_24h","dox_72h"),las=2, cex.axis = 0.8 , space=c(0.2,0.2,0.5,0.2)), mean_table$V2[5:8],mean_table$V3[5:8],add=TRUE)
plotCI(barplot(mean_table$V2[1:4],ylim=c(0,max(mean_table$V2[1:4]+mean_table$V3[1:4])),main ='miR-34a abundance per group', ylab="miR-34a copies per ng of total RNA",col=rep("#f94144",4) ,names.arg = c("ctrl_24h","ctrl_72h","dox_24h","dox_72h"),las=2, cex.axis = 0.8, space=c(0.2,0.2,0.5,0.2)) , mean_table$V2[1:4],mean_table$V3[1:4],add=TRUE)
#dev.off()
#pdf(file="RTddPCR_dox_miRNAs_per_ng.pdf", width=5.5 , height =4.5 )
par(mfrow = c(1, 1),mar= c(5.1, 4.1, 4.1, 2.1))
plotCI(barplot(mean_table$V2,ylim=c(0,max(mean_table$V2+mean_table$V3)),main ='miR-34a and miR-21 abundance per group', ylab="Copies of miRNA / ng of total RNA", col=c(rep("#f94144",4), rep("dimgrey",4)) ,names.arg = rep(c("ctrl_24h","ctrl_72h","dox_24h","dox_72h"),2),las=2, cex.axis = 0.8, space = c(0.2,0.2,0.5,0.2,0.5,0.2,0.5,0.2) ), mean_table$V2,mean_table$V3,add=TRUE)
#dev.off()
pdf(file="RTddPCR_dox_miRNAs_per_ng_log_scale.pdf", width=5.5 , height =4.5 )
plotCI(barplot(mean_table$V2,ylim=c(100,100000),log="y",main ='miR-34a and miR-21 abundance per group', ylab="Copies of miRNA / ng of total RNA", col=c(rep("#f94144",4), rep("dimgrey",4))  ,names.arg = rep(c("ctrl_24h","ctrl_72h","dox_24h","dox_72h"),2),las=2, cex.axis = 0.8, space = c(0.2,0.2,0.5,0.2,0.5,0.2,0.5,0.2) ), mean_table$V2,mean_table$V3,add=TRUE)
axis(2, c(2e2,3e2,4e2,5e2,6e2,7e2,8e2,9e2,2e3,3e3,4e3,5e3,6e3,7e3,8e3,9e3,2e4,3e4,4e4,5e4,6e4,7e4,8e4,9e4), labels= F,col.ticks = "dimgrey",las=2)
dev.off()

#ANOVA test for miR-34a
stat_table=data.frame(var = table$copies_per_ng[which(table$target_miRNA == "miR-34a")],
                      treatment = factor(c(rep("ctrl",6),rep("dox",6))),
                      time_point= factor(c(rep("24h",3),rep("72h",3),rep("24h",3),rep("72h",3))))
summary(aov(var ~ treatment + time_point, stat_table))
t.test(var ~ treatment, stat_table[which(stat_table$time_point == "24h"),]) #No significant difference between ctrl_24h and dox_24h for miR-34
t.test(var ~ treatment, stat_table[which(stat_table$time_point == "72h"),]) #Significant difference between ctrl_72h and dox_72h for miR-34
FC_dox_24h=mean_table[which(mean_table$V1=="dox-24h-r1_miR-34a"),2]/mean_table[which(mean_table$V1=="ctrl-24h-r1_miR-34a"),2] #Fold-change between not treated and treated with doxorubicin after 24 hours
FC_dox_72h=mean_table[which(mean_table$V1=="dox-72h-r1_miR-34a"),2]/mean_table[which(mean_table$V1=="ctrl-72h-r1_miR-34a"),2] #Fold-change between not treated and treated with doxorubicin after 72 hours
#ANOVA test for miR-21
stat_table=data.frame(var = table$copies_per_ng[which(table$target_miRNA == "miR-21")],
                      treatment = factor(c(rep("ctrl",6),rep("dox",6))),
                      time_point= factor(c(rep("24h",3),rep("72h",3),rep("24h",3),rep("72h",3))))
summary(aov(var ~ treatment + time_point, stat_table))
#No significant effect of treatment or treatment time for miR-21



##Figure 4C##

## Import data
plate <- read.csv("20210412_mimic_HCT116_miR34a_for_analysis.csv", header = T, sep = ';',dec=",",stringsAsFactors = F)
sample_name=c(rep(c("ctrl-1h-r1","ctrl-1h-r2","ctrl-1h-r3","ctrl-24h-r1","ctrl-24h-r2","ctrl-24h-r3",
                    "1nM-1h-r1","1nM-1h-r2","1nM-1h-r3","1nM-24h-r1","1nM-24h-r2","1nM-24h-r3",
                    "10nM-1h-r1","10nM-1h-r2","10nM-1h-r3","10nM-24h-r1","10nM-24h-r2","10nM-24h-r3"),2))
target_miRNA=c(rep("miR-34a",18),rep("miR-21",18))
coeff_dilution=c(rep(1,6),rep(100,6),rep(1000,6),rep(10,18)) 
#Reverse-transcription reactions were diluted to 1/100 for miR-34a quantification of samples transfected with 1 nM miR-34a duplex,
#to 1/1000 for miR-34a quantification of samples transfected with 10 nM miR-34a duplex,
#and to 1/10 for miR-21 quantification of all samples.
ng_RT=rep(10,36) #Reverse-transcription reactions were performed with 10 ng total RNA in final volume of 15 µL.
volume_RT=15
volume_RT_for_ddPCR=1.33 #ddPCR quantifications were performed from 1.33 µL reverse-transcription reactions.

table=as.data.frame(cbind(plate$Well,sample_name,target_miRNA,coeff_dilution,ng_RT,plate$Copies.20µLWell),stringsAsFactors = F)
colnames(table)=c("well","sample","target_miRNA","dilution","ng_RT","copies_per_well")
table$copies_per_ng = as.numeric(table$copies_per_well)/(((volume_RT_for_ddPCR*as.numeric(table$ng_RT))/volume_RT)/as.numeric(table$dilution))

#miR-34a and miR-21 levels in biological replicates
#par(mfrow = c(1, 2))
#barplot(table$copies_per_ng[which(table$target_miRNA == "miR-34a")], main ='miR-34a abundance\nper sample',names.arg =table$sample[which(table$target_miRNA == "miR-34a")] ,las =2 , ylab = "Copies of miRNA / ng of total RNA", cex.names = 0.75)
#barplot(table$copies_per_ng[which(table$target_miRNA == "miR-21")], main ='miR-21 abundance\nper sample',names.arg =table$sample[which(table$target_miRNA == "miR-21")] ,las =2, ylab = "Copies of miRNA / ng of total RNA", cex.names = 0.75)

mean_table=c()
pas=seq(1,34,3)
rep=rep(2,length(pas))
for (i in 1:length(pas)) {
  mean_table=rbind(mean_table,c(paste(table$sample[pas[i]],table$target_miRNA[pas[i]],sep="_"), mean(table$copies_per_ng[pas[i]:(pas[i]+rep[i])]),sd(table$copies_per_ng[pas[i]:(pas[i]+rep[i])])))
}
mean_table=as.data.frame(mean_table)
mean_table$V2=as.numeric(as.character(mean_table$V2))
mean_table$V3=as.numeric(as.character(mean_table$V3))
#pdf(file="RTddPCR_mimic_miRNAs_per_ng_split.pdf", width= 5.5, height = 4.5)
par(mfrow = c(1, 2),mar= c(5.1, 4.1, 4.1, 2.1))
plotCI(barplot(mean_table$V2[7:12],ylim=c(0,max(mean_table$V2[7:12]+mean_table$V3[7:12])),main ='miR-21 abundance per group', ylab="miR-21 copies per ng of total RNA",col=rep("dimgrey",6), names.arg = c("ctrl_1h","ctrl_24h","1nM_1h","1nM_24h","10nM_1h","10nM_24h"),las=2, cex.axis = 0.8 , space=c(0.2,0.2,0.5,0.2,0.5,0.2)), mean_table$V2[7:12],mean_table$V3[7:12],add=TRUE)
plotCI(barplot(mean_table$V2[1:6],ylim=c(0,max(mean_table$V2[1:6]+mean_table$V3[1:6])),main ='miR-34a abundance per group', ylab="miR-34a copies per ng of total RNA",col=rep("#f94144",6) ,names.arg = c("ctrl_1h","ctrl_24h","1nM_1h","1nM_24h","10nM_1h","10nM_24h"),las=2, cex.axis = 0.8, space=c(0.2,0.2,0.5,0.2,0.5,0.2)) , mean_table$V2[1:6],mean_table$V3[1:6],add=TRUE)
#dev.off()
#pdf(file="RTddPCR_mimic_miRNAs_per_ng.pdf", width=5.5 , height =4.5 )
par(mfrow = c(1, 1),mar= c(5.1, 4.1, 4.1, 2.1))
plotCI(barplot(mean_table$V2,ylim=c(0,max(mean_table$V2+mean_table$V3)),main ='miR-34a and miR-21 abundance per group', ylab="Copies of miRNA / ng of total RNA", col=c(rep("#f94144",6), rep("dimgrey",6)) ,names.arg = rep(c("ctrl_1h","ctrl_24h","1nM_1h","1nM_24h","10nM_1h","10nM_24h"),2),las=2, cex.axis = 0.8, space = c(0.2,0.2,0.5,0.2,0.5,0.2,0.2,0.2,0.5,0.2,0.5,0.2) ), mean_table$V2,mean_table$V3,add=TRUE)
#Zoom in
plotCI(barplot(mean_table$V2,ylim=c(0,4e4),main ='miR-34a and miR-21 abundance per group', ylab="Copies of miRNA / ng of total RNA", col=c(rep("#f94144",6), rep("dimgrey",6)) ,names.arg = rep(c("ctrl_1h","ctrl_24h","1nM_1h","1nM_24h","10nM_1h","10nM_24h"),2),las=2, cex.axis = 0.8, space = c(0.2,0.2,0.5,0.2,0.5,0.2,0.2,0.2,0.5,0.2,0.5,0.2) ), mean_table$V2,mean_table$V3,add=TRUE)
#dev.off()
pdf(file="RTddPCR_mimic_miRNAs_per_ng_log_scale.pdf", width=5.5 , height =4.5 )
par(mfrow = c(1, 1),mar= c(5.1, 4.1, 4.1, 2.1))
#  plotCI(barplot(mean_table$V2,ylim=c(0,30000),main ='miR-34a and miR-21 RT-ddPCR', ylab="miRNA copies per ng", col=c(rep("#f94144",6), rep("dimgrey",6))  ,names.arg = rep(c("ctrl_1h","ctrl_24h","ctrl_48h","mimic_1h","mimic_24h","mimic_48h"),2),las=2, cex.axis = 0.8, space = c(0.2,0.2,0.2,0.5,0.2,0.2,0.5,0.2,0.2,0.5,0.2,0.2) ), mean_table$V2,mean_table$V3,add=TRUE)
plotCI(barplot(mean_table$V2,ylim=c(100,10000000),log="y",main ='miR-34a and miR-21 abundance per group', ylab="Copies of miRNA / ng of total RNA", col=c(rep("#f94144",6), rep("dimgrey",6))  ,names.arg = rep(c("ctrl_1h","ctrl_24h","1nM_1h","1nM_24h","10nM_1h","10nM_24h"),2),las=2, cex.axis = 0.8, space = c(0.2,0.2,0.5,0.2,0.5,0.2,0.5,0.2,0.5,0.2,0.5,0.2) ), mean_table$V2,mean_table$V3,add=TRUE)
axis(2, c(1e3,1e5,1e7), labels= T,cex.axis = 0.8,las=2)
axis(2, c(2e2,3e2,4e2,5e2,6e2,7e2,8e2,9e2,2e3,3e3,4e3,5e3,6e3,7e3,8e3,9e3,2e4,3e4,4e4,5e4,6e4,7e4,8e4,9e4,2e5,3e5,4e5,5e5,6e5,7e5,8e5,9e5,2e6,3e6,4e6,5e6,6e6,7e6,8e6,9e6), labels= F,col.ticks = "dimgrey",las=2)
dev.off()

#Fold-change calculations (transfected with miR-34a duplex compared with transfected with control siRNA)
FC_miR34a_1nM_1h=mean_table[which(mean_table$V1=="1nM-1h-r1_miR-34a"),2]/mean_table[which(mean_table$V1=="ctrl-1h-r1_miR-34a"),2] # 1 hour of treatment with 1 nM duplex
FC_miR34a_1nM_24h=mean_table[which(mean_table$V1=="1nM-24h-r1_miR-34a"),2]/mean_table[which(mean_table$V1=="ctrl-24h-r1_miR-34a"),2] # 24 hours of treatment with 1 nM duplex
FC_miR34a_10nM_1h=mean_table[which(mean_table$V1=="10nM-1h-r1_miR-34a"),2]/mean_table[which(mean_table$V1=="ctrl-1h-r1_miR-34a"),2] # 1 hour of treatment with 10 nM duplex
FC_miR34a_10nM_24h=mean_table[which(mean_table$V1=="10nM-24h-r1_miR-34a"),2]/mean_table[which(mean_table$V1=="ctrl-24h-r1_miR-34a"),2] # 24 hours of treatment with 10 nM duplex
