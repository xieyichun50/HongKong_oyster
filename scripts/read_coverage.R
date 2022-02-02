library(dplyr)
library(tidyr)
library(ggplot2)
library(eoffice)

coverage.tab<-as.data.frame(matrix(NA, ncol=5, nrow=0))
names(coverage.tab)<-c("Chr","Pos","Depth","Sample","Gene")
samples<-c("R25Ca","R25Cb","R32Ca","R32Cb")

for (i in 1:4) {
  coverage.tab.sub<-read.delim(paste0("al.", samples[i], ".coverage"), header = F)
  names(coverage.tab.sub)<-c("Chr","Pos","Depth")
  for (j in 1:nrow(coverage.tab.sub)) {
    if (coverage.tab.sub$Chr[j] == "ScQ1DuM_12767") {
      coverage.tab.sub$Gene[j]<-"Mhon_031090"
      coverage.tab.sub$Pos[j]<-(0-coverage.tab.sub$Pos[j])
    } else if (coverage.tab.sub$Chr[j] == "ScQ1DuM_6034" & 
               coverage.tab.sub$Pos[j] >= 7925334 &
               coverage.tab.sub$Pos[j] <= 7928380) {
      coverage.tab.sub$Gene[j]<-"Mhon_006787"
    } else if (coverage.tab.sub$Chr[j] == "ScQ1DuM_6034" & 
               coverage.tab.sub$Pos[j] >= 60708770 &
               coverage.tab.sub$Pos[j] <= 60709816) {
      coverage.tab.sub$Gene[j]<-"Mhon_010093"
      coverage.tab.sub$Pos[j]<-(0-coverage.tab.sub$Pos[j])
    }
  }
  coverage.tab.sub$Sample=samples[i]
  coverage.tab<-rbind(coverage.tab, coverage.tab.sub)
}

geneid<-"Mhon_006787"
plotdata<-coverage.tab[coverage.tab$Gene==geneid,]
plotdata<-plotdata[plotdata$Sample == "R25Ca" | plotdata$Sample == "R32Cb",]
geneid<-"Mhon_010093"
plotdata<-coverage.tab[coverage.tab$Gene==geneid,]
plotdata<-plotdata[plotdata$Sample == "R25Cb" | plotdata$Sample == "R32Cb",]
geneid<-"Mhon_031090"
plotdata<-coverage.tab[coverage.tab$Gene==geneid,]
plotdata<-plotdata[plotdata$Sample == "R25Ca" | plotdata$Sample == "R32Cb",]
a<-plotdata %>% 
  #mutate(Groups = factor(Groups, levels = group.order.sub)) %>% 
  ggplot(aes(x = Pos, y = Depth, 
             fill = Depth))+
  labs(title = geneid, 
       fill = "Coverage",
       x = "", 
       y = "")+
  geom_bar(stat = "identity", position = "dodge")+
  scale_fill_gradient(high = "#FF0000", low = "#0000FF", 
                      #limits = c(0.5, max(log10(plotdata$Depth+1)))
                      limits = c(0, max(plotdata$Depth)))+
  #geom_vline(xintercept = 0, color = "black", linetype = "solid")+
  guides(size = guide_legend(order = 1))+
  theme(panel.background = element_rect(colour = "black", fill = "white", size = 0.5), 
        axis.title.x = element_text(colour = "black", size = 12), 
        axis.title.y = element_text(colour = "white"), 
        axis.text = element_text(colour = "black", size = 12), 
        axis.text.x = element_text(size = 12, colour = "white"), 
        axis.text.y = element_text(size = 12),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.text = element_blank(),
        legend.title = element_text(size = 12), 
        plot.title = element_text(size = 12),
        legend.position = "none")+
  facet_grid(Sample~., scales = "fixed", space = "fixed")

a

f = paste0("coverage.",geneid,".pptx")
topptx(a,f, width = 4, height = 4, units = "cm")
ggsave(paste0("coverage.",geneid,".tiff"), width = 4, height = 4, units = "in", dpi = 300)
