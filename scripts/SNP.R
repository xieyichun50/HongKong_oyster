library(tidyr)

SNP.raw<-read.delim("SNP_table_hox.txt", header = FALSE)
#SNP.raw<-read.delim("SNP_table_hsp.txt", header = FALSE)
names(SNP.raw)<-SNP.raw[1,]
SNP.raw<-SNP.raw[-1,]
names(SNP.raw)[1]="CHROM:POS"

SNPlist<-as.data.frame(matrix(NA, nrow = 0, ncol = 3))
names(SNPlist)<-c("CHROM:POS","count","sample")

for (i in 2:ncol(SNP.raw)){
  SNP.sub<-SNP.raw[,c(1,i)]
  SNP.sub$sample<-names(SNP.sub)[2]
  names(SNP.sub)[2]<-"count"
  SNPlist<-rbind(SNPlist, SNP.sub)
}

#sampleID<-read.delim("sample_info_hsp.txt", header = TRUE)
sampleID<-read.delim("sample_info_hox.txt", header = TRUE)
SNPlist<-merge(SNPlist,sampleID, by = "sample", all.x = TRUE)

SNPlist<-separate(SNPlist, count, into = c("Ref","A1","A2"), sep = ",")
SNPlist$A1<-as.numeric(SNPlist$A1)
SNPlist$A2<-as.numeric(SNPlist$A2)
SNPlist$Ref<-as.numeric(SNPlist$Ref)
SNPlist<-SNPlist[is.na(SNPlist$Ref)==F,]

SNPlist$A2[is.na(SNPlist$A2)]=0
SNPlist$ratio<-(SNPlist$A1+SNPlist$A2)/(SNPlist$A1+SNPlist$A2+SNPlist$Ref)

SNP.matrix<-as.data.frame(SNP.raw$`CHROM:POS`)
names(SNP.matrix)[1]<-"CHROM:POS"
for (i in 1:nrow(sampleID)) {
  SNP.matrix.sub<-SNPlist[SNPlist$sample==sampleID$sample[i],c("CHROM:POS", "ratio")]
  #SNP.matrix.sub$ratio[SNP.matrix.sub$ratio=="NaN"]<-0
  #SNP.matrix.sub$ratio[SNP.matrix.sub$ratio==0]<-0.001
  names(SNP.matrix.sub)[2]<-sampleID$sample[i]
  
  SNP.matrix<-merge(SNP.matrix, SNP.matrix.sub, by = "CHROM:POS", all.x = TRUE)
}

row.names(SNP.matrix)<-SNP.matrix[,1]
SNP.matrix<-SNP.matrix[,-1]
#write.table(SNP.matrix, "SNP.matrix.hsp.txt", quote = F, sep = "\t", row.names = T)
#write.table(SNP.matrix, "SNP.matrix.hox.txt", quote = F, sep = "\t", row.names = T)

SNPlist<-SNPlist[SNPlist$ratio !="NaN",]
SNPlist<-SNPlist[SNPlist$ratio >=0.05,]

SNPlist.summary<-as.data.frame(xtabs(~`CHROM:POS`+batch, SNPlist))
SNPlist.summary<-SNPlist.summary[SNPlist.summary$Freq>0,]

batchID<-read.delim("batch_info.txt", header = TRUE)
SNPlist.summary<-merge(SNPlist.summary, batchID, by = "batch", all.x = TRUE)

SNPlist.summary$sample_ratio<-SNPlist.summary$Freq/SNPlist.summary$sample_n
sitelist<-unique(SNPlist.summary$CHROM.POS[duplicated(SNPlist.summary$CHROM.POS)==T])
SNPlist.summary.test<-SNPlist.summary[-which(SNPlist.summary$CHROM.POS %in% sitelist),]
SNPlist.summary.test<-SNPlist.summary.test[SNPlist.summary.test$Freq>1 & SNPlist.summary.test$sample_ratio>0,]
write.table(SNPlist.summary.test, 
#            "SNPlist_hsp_single.txt", sep = "\t", quote = F, row.names = F)
            "SNPlist_hox_single.txt", sep = "\t", quote = F, row.names = F)

##merge Shenzhen Bay
SNPlist1<-SNPlist
SNPlist1$batch[SNPlist1$batch != "HKU"]<-"SZB"

SNPlist.summary1<-as.data.frame(xtabs(~`CHROM:POS`+batch, SNPlist1))
SNPlist.summary1<-SNPlist.summary1[SNPlist.summary1$Freq>1,]

batchID1<-read.delim("batch_info1.txt", header = TRUE)
SNPlist.summary1<-merge(SNPlist.summary1, batchID1, by = "batch", all.x = TRUE)

SNPlist.summary1$sample_ratio<-SNPlist.summary1$Freq/SNPlist.summary1$sample_n
#SNPlist.summary1<-SNPlist.summary1[SNPlist.summary1$sample_ratio>0.5,]
sitelist<-unique(SNPlist.summary1$CHROM.POS[duplicated(SNPlist.summary1$CHROM.POS)])
SNPlist.summary1<-SNPlist.summary1[-which(SNPlist.summary1$CHROM.POS %in% sitelist),]
SNPlist.summary1<-SNPlist.summary1[SNPlist.summary1$sample_ratio>0,]
write.table(SNPlist.summary1, 
#            "SNPlist_hsp.txt", sep = "\t", quote = F, row.names = F)
             "SNPlist_hox.txt", sep = "\t", quote = F, row.names = F)
