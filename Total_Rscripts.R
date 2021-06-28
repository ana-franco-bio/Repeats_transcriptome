###### R code Setaria transcriptome Repeats######
#Calculate the GC ratio from files of python Script
#create the function to apply on tables
l <- function(my_data) {colnames(my_data) <- c("1","2","3","G", "C", "A", "T")
Total_G <- sum(my_data$G)
Total_C <- sum(my_data$C)
Total_A <- sum(my_data$A)
Total_T <- sum(my_data$T)
Total <- sum(Total_G, Total_C, Total_A, Total_T)
PercentageGC <- ((Total_G + Total_C)/Total)
return(PercentageGC)
}

l(my_data)

#list of txt files
fs <- list.files(path = "C:\\Users\\Ana Luiza\\Desktop\\GC_propP2", recursive = TRUE,
                 pattern = ".*.txt", 
                 full.names = TRUE)

lst <- vector("list", length(fs))

#Read files in to list
for(i in 1:length(fs)) {
  lst[[i]] <- read.csv(fs[i], header = TRUE, check.names = FALSE, sep = "\t")
}

names(lst) <- fs
fs

#Apply a function to the list
results <- sapply(lst, l)
results

resultsOK <- cbind(results)
resultsOK

#Save as a .xlsx table
library(openxlsx)
write.xlsx(results4, "C:\\Users\\Ana Luiza\\Desktop\\GC_ratio_Outputtt.xlsx")

###Calculate the mode of repeats similarity
library(stats)
filelist = list.files(path = "C:\\Users\\Ana Luiza\\Desktop\\Files_similarity",
                      pattern = ".*.txt" ,
                      full.names = TRUE)
str(filelist)
Peaksimilarity=c()
num = 0
for (file in filelist) {
  df <- read.table(file = file)
  b <- density(df$V5)
  e = mean(df$V5)
  cv = CV(mean(df$V5), sd(df$V5))
  Peaksimilarity=rbind(Peaksimilarity, c(file, b$x[which.max(b$y)], e, cv))
}
Peaksimilarity

#Save in a spreadsheet
library(openxlsx)
write.xlsx(Peaksimilarity, "C:\\Users\\Ana Luiza\\Desktop\\Peak_Similarity.xlsx")


#Normality tests among genimic data
shapiro.test(Genomic_data)
Genomic_data$Similarity <- log10(Genomic_data$Similarity)
shapiro.test(Genomic_data$Similarity)
ks.test(Genomic_data$Similarity, "pnorm", mean(Genomic_data$Similarity), sd(Genomic_data$Similarity))
ks.test(Genomic_data$Log_GP, "pnorm", mean(Genomic_data$Log_GP), sd(Genomic_data$Log_GP))
ks.test(Genomic_data$GC_ratio, "pnorm", mean(Genomic_data$GC_ratio), sd(Genomic_data$GC_ratio))

#Correlation tests
library(corrplot)
library(RColorBrewer)
library(ggpubr)
#Pearson's coefficient
M<- cor(Genomic_data, method = "pearson")
M
Cor_Graph<- corrplot(M, method = "color", addCoef.col="black", col = brewer.pal(n = 8, name = "PuOr")) 
ggsave(M, file="Correlation.png", height=5, width=6)
#Test the significance
res <- cor.test(Genomic_data$GC_ratio, Genomic_data$Log_GP, 
                method = "pearson")
res

####Linear models####

m1 <- lm(Log_TP~ Tissue + Condition + Repeat_Type, data= new_data)
summary (m1)

m2 <- lm(Log_TP~ Tissue + Condition + Annotation, data= new_data)
summary (m2)

m3 <- lm(Log_TP~ (Similarity + GC_ratio)*Repeat_Type + Tissue + Condition, data= Copia_Gypsy_)
summary (m3)

m4 <- lm(Log_TP~ (Log_GP)*Repeat_Type + Tissue + Condition, data= Copia_Gypsy_)
summary (m4)

m5 <- lm(Log_TP~ (Similarity + GC_ratio)*Repeat_Type, data= Copia_Gypsy_Leaf)
summary (m5)

m6 <- lm(Log_TP~ (Log_GP)*Repeat_Type, data= Copia_Gypsy_Leaf)
summary (m5)

m7 <- lm(Log_TP~ (Similarity + GC_ratio)*Repeat_Type, data= Copia_Gypsy_Flor)
summary (m7)

m7 <- lm(Log_TP~ (Log_GP)*Repeat_Type, data= Copia_Gypsy_Flor)
summary (m8)

#Linear model graphs
d1 <- ggplot(Copia_Gypsy_G, aes(x=GC_ratio, y=Log_GP, color=Repeat_Type, shape=Repeat_Type, alpha=1)) +
  geom_point(aes(size=Repeat_Type, alpha=1)) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  scale_color_manual(values=c('purple', 'orange')) +
  scale_shape_manual(values=c(17,19))  +
  scale_size_manual(values=c(3.5,3.5)) +
  theme_bw()
d1
d2 <- ggplot(Copia_Gypsy_G, aes(x=Similarity, y=Log_GP, color=Repeat_Type, shape=Repeat_Type, alpha=1)) +
  geom_point(aes(size=Repeat_Type, alpha=1)) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  scale_color_manual(values=c('purple', 'orange')) +
  scale_shape_manual(values=c(17,19))  +
  scale_size_manual(values=c(3.5,3.5)) +
  theme_bw()
d2
#Plot all graphs together
library(ggpmisc)
install.packages("ggpmisc")
figure_d <- ggarrange(d1, d2, 
                      labels = c("a", "b", size=1),
                      ncol = 2,
                      common.legend = TRUE,legend="bottom",
                      heights = c(1,1)) 


####Graphs####
#Bar graphs and boxplots
#Bar graphs genome proportion
#Genome Proportion
GP<-read.table("clipboard", header=TRUE) 
GP
names(GP)
str(GP)
attach(GP)
GP1 <- as.data.frame(GP)

GP1$Repeat_Type <-factor(GP1$Repeat_Type,
                         levels = c("Gypsy", "LTR", "Unknown", "Copia", "Satellite", "ClassII", "LINE", "Pararetrovirus"))
GP_plot <- ggplot(GP1, aes(x = reorder(Repeat_Type, GP),fill= Repeat_Type, y=GP)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Genome Proportion") +
  coord_flip() +
  guides(fill = guide_legend(reverse = FALSE))

GP_plot
GP_plot + scale_fill_manual(values = c("darkblue", "grey", "darkmagenta", "red", "yellow", "lightsalmon", "darkgreen", "black"))
ggsave("GPFINAL_OK2.png", dpi=200, height=4, width=5, units="in")

#Boxplots range similarity and gc ratio
rangedata <- read.table("clipboard", header = TRUE, check.names = TRUE)
rangedatasim <- as.data.frame(rangedata)

rangedatasim$Repeat_Type <-factor(rangedatasim$Repeat_Type,
                                  levels = c("Gypsy", "LTR", "Unknown", "Copia", "Satellite", "ClassII", "LINE", "Pararetrovirus"))
rangedatasim$Repeat_Type

#plot similarity
p<-ggplot(rangedatasim, aes(x=Repeat_Type, y=Similarity, fill = Repeat_Type)) +
  geom_boxplot(lwd=1)
p + theme_classic() +
  
  scale_fill_manual(values = c("darkblue", "grey", "darkmagenta", "red", "yellow", "lightsalmon", "darkgreen", "black")) 

ggsave("rangeOKsimilarity fill.png", dpi=300, height=4, width=5, units="in")

#Plot GC ratio
g<-ggplot(rangedatasim, aes(x=Repeat_Type, y=GC_ratio, fill=Repeat_Type)) +
  geom_boxplot(lwd=1)
g + theme_classic() + scale_fill_manual(values = c("darkblue", "grey", "darkmagenta", "red", "yellow", "lightsalmon", "darkgreen", "black"))
ggsave("range GC ratio fill.png", dpi=300, height=4, width=5, units="in")


scale_fill_manual(values = c("darkblue", "grey", "darkmagenta", "red", "yellow", "lightsalmon", "darkgreen", "black"))

GP_plot + scale_color_manual(values = c("darkblue", "grey", "darkmagenta", "red", "yellow", "lightsalmon", "darkgreen", "black"))
ggsave("range similarity.png", dpi=300, height=4, width=5, units="in")

#Bar graphs TP
#load data
Data_<-read.table("clipboard", header=TRUE) 
Data_
names(Data_)
str(Data_)
attach(Data_)
Data1 <- as.data.frame(Data_)
Data1$Condition <- as.character(Data1$Condition)
str(Data1)
#select specific data
subStem <- subset(Data1, Tissue=="Stem")
subCrown <- subset(Data1, Tissue=="Crown")
subLeaf <- subset(Data1, Tissue=="Leaf")
subFlor <- subset(Data1, Tissue =="Inflorescence")

#Graphs - Stem and Crown bars
theme_set(theme_bw())
bar1 <- ggplot(subCrown, aes(x = reorder(Repeat_Type, -TP), fill=Condition, y=TP)) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle("Transcriptome Proportion - Crown") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(breaks = seq(0, 0.8, by=0.1), limits = c(0, 0.8))

theme_set(theme_bw())
bar2 <- ggplot(subStem, aes(x = reorder(Repeat_Type, -TP), fill=Condition, y=TP)) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle("Transcriptome Proportion - Stem") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(breaks = seq(0, 0.8, by=0.1), limits = c(0, 0.8))

#colors
bar + scale_fill_manual(values=c("red", "darkblue"))
ggsave("June_Crown.png", dpi=200, height=4, width=5, units="in")

#Subfamily types
theme_set(theme_bw())

bar3 <- ggplot(subCrown, aes(x = reorder(Annotation, TP), fill=Condition, y=TP)) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle("Families Proportion - Crown") +
  theme(axis.text.x = element_text(angle = 90)) +
  coord_flip()

bar3 + scale_fill_manual(values=c("red", "darkblue")) + scale_y_continuous(limits = c(0, 30))
ggsave("June_Crown_FamiliespaperOK.png", dpi=300, height=4, width=5, units="in")

bar4 <- ggplot(subStem, aes(x = reorder(Annotation, TP), fill=Condition, y=TP)) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle("Families Proportion - Stem") +
  theme(axis.text.x = element_text(angle = 90)) +
  coord_flip()

bar4 + scale_fill_manual(values=c("red", "darkblue")) + scale_y_continuous(limits = c(0, 30))
ggsave("June_Stem_FamiliespaperOK.png", dpi=300, height=4, width=5, units="in")

