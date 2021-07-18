###### R code Setaria transcriptome repeats######
##Calculate the GC ratio
#Create the function to apply on outputs of script of "get_cgcontent.pl" from Meneguin, 2009
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

#list of txt files outputs
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

#Apply the "l" function to the list of files
results_GC <- sapply(lst, l)
results_GC_final <- cbind(results_GC)

###Calculate the mode of repeats similarity
#Output of scripts in Get_reads_similarity
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


#Plot simple histogram
Similarity_ <- Sim_dataOK$similarity
a1 <- hist(Similarity_, col = NULL, breaks = 40, xlim=c(85,100))

#Graphs Figure1
#Bar graphs and boxplots
#Genome Proportion (GP)
GP <-read.table("Genomic_features.txt",header = TRUE, check.names = FALSE, sep = "\t") 
attach(GP)
GP1 <- as.data.frame(GP)
#GP bar graph
GP1$Repeat_Type <-factor(GP1$Repeat_Type,
                         levels = c("Gypsy", "LTR", "Unknown", "Copia", "Satellite", "ClassII", "LINE", "Pararetrovirus"))
GP_plot <- ggplot(GP1, aes(x = reorder(Repeat_Type, GP),fill= Repeat_Type, y=GP)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Genome Proportion") +
   guides(fill = guide_legend(reverse = FALSE))

GP_plot + scale_fill_manual(values = c("darkblue", "grey", "darkmagenta", "red", "yellow", "lightsalmon", "darkgreen", "black"))
ggsave("GPFINAL_OK2.png", dpi=200, height=4, width=5, units="in")

#Boxplots range similarity and gc ratio
GP$Repeat_Type <-factor(GP1$Repeat_Type,
                                  levels = c("Gypsy", "LTR", "Unknown", "Copia", "Satellite", "ClassII", "LINE", "Pararetrovirus"))
GP$Repeat_Type

#plot similarity
p<-ggplot(GP1, aes(x=Repeat_Type, y=Similarity, fill = Repeat_Type)) +
  geom_boxplot(lwd=1)
p + theme_classic() +
  scale_fill_manual(values = c("darkblue", "grey", "darkmagenta", "red", "yellow", "lightsalmon", "darkgreen", "black")) 
p
#Plot GC ratio
g<-ggplot(GP1, aes(x=Repeat_Type, y=GC_ratio, fill=Repeat_Type)) +
  geom_boxplot(lwd=1)
g + theme_classic() + scale_fill_manual(values = c("darkblue", "grey", "darkmagenta", "red", "yellow", "lightsalmon", "darkgreen", "black"))
g
#Bar graphs TP
#load data
Data_<-read.table("Mean_TP.txt", header = TRUE, check.names = FALSE, sep = "\t")
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
bar1 + scale_fill_manual(values=c("red", "darkblue"))

bar2 <- ggplot(subStem, aes(x = reorder(Repeat_Type, -TP), fill=Condition, y=TP)) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle("Transcriptome Proportion - Stem") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(breaks = seq(0, 0.8, by=0.1), limits = c(0, 0.8))
bar2 + scale_fill_manual(values=c("red", "darkblue"))

#Subfamily types
theme_set(theme_bw())
bar3 <- ggplot(subCrown, aes(x = reorder(Annotation, TP), fill=Condition, y=TP)) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle("Families Proportion - Crown") +
  theme(axis.text.x = element_text(angle = 90)) +
  coord_flip()
bar3 + scale_fill_manual(values=c("red", "darkblue")) + scale_y_continuous(limits = c(0, 30))


bar4 <- ggplot(subStem, aes(x = reorder(Annotation, TP), fill=Condition, y=TP)) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle("Families Proportion - Stem") +
  theme(axis.text.x = element_text(angle = 90)) +
  coord_flip()

bar4 + scale_fill_manual(values=c("red", "darkblue")) + scale_y_continuous(limits = c(0, 30))

#Leaf and Inflorescence
#Transcriptome proportion Leaf
#Leaf TP
subLeaf$Repeat_Type <-factor(subLeaf$Repeat_Type,
                             levels = c("Unknown", "Copia", "ClassII", "LTR", "Gypsy", "LINE", "Satellite", "Pararetrovirus"))
#Plot
TPok<- ggplot(subLeaf, aes(x = reorder(Repeat_Type, - TP),fill= Repeat_Type, y=TP)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Leaf")
#Set colors
TPok + scale_fill_manual(values=c("grey", "red", "lightsalmon", "darkmagenta", "darkblue", "darkgreen", "yellow", "black"))

#Inflorescence TP
subFlor$Repeat_Type <-factor(subFlor$Repeat_Type,
                             levels = c("Unknown", "ClassII", "Copia", "Gypsy", "LTR", "LINE", "Satellite", "Pararetrovirus"))
#Plot
TP_inflor<- ggplot(subFlor, aes(x = reorder(Repeat_Type, - TP),fill= Repeat_Type, y=TP)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Inflorescence")
#Set colors
TP_inflor + scale_fill_manual(values=c("grey", "lightsalmon", "red", "darkblue", "darkmagenta", "darkgreen", "yellow", "black"))

#Normality tests among genomic data
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
library(ggplot2)
#Pearson's coefficient
M<- cor(Genomic_data, method = "pearson")
M
Cor_Graph<- corrplot(M, method = "color", addCoef.col="black", col = brewer.pal(n = 8, name = "PuOr")) 
ggsave(M, file="Correlation.png", height=5, width=6)
#Test the significance
res <- cor.test(Genomic_data$GC_ratio, Genomic_data$Log_GP, 
                method = "pearson")
res


###Linear models###

#Load data
GP <-read.table("Genomic_features.txt",header = TRUE, check.names = FALSE, sep = "\t")
Ribo_dep <- read.table("All_repeats_ribodepleted.txt",header = TRUE, check.names = FALSE, sep = "\t")
Ribo_dep_no_outlier <- read.table("Ribo_dep_no_outlier.txt", header = TRUE, check.names = FALSE, sep = "\t")
Leaf_polyA <- read.table("Leaf_polyA.txt",header = TRUE, check.names = FALSE, sep = "\t")
Inflor_polyA <- read.table("Inflor_polyA.txt",header = TRUE, check.names = FALSE, sep = "\t")

#Preparing Data (do it for each library)
Ribo_dep$Cluster <- as.factor(Ribo_dep$Cluster)
Ribo_dep$Tissue <- as.factor(Ribo_dep$Tissue)
Ribo_dep$Condition <- as.factor(Ribo_dep$Condition)
Ribo_dep$Repeat_Type <- as.factor(Ribo_dep$Repeat_Type)
Ribo_dep$Annotation <- as.factor(Ribo_dep$Annotation)
str(Ribo_dep)

#Subset Ty1/copia and Ty3/gypsy elements
Copia_Gypsy_G <- subset(GP, GP$Repeat_Type == "Copia"|GP$Repeat_Type =="Gypsy")
  
Copia_Gypsy_ <-subset(Ribo_dep, Ribo_dep$Repeat_Type == "Copia"|Ribo_dep$Repeat_Type =="Gypsy")

Copia_Gypsy_Leaf <-subset(Leaf_polyA, Leaf_polyA$Repeat_Type == "Copia"|Leaf_polyA$Repeat_Type =="Gypsy")

Copia_Gypsy_Inflor <-subset(Inflor_polyA, Inflor_polyA$Repeat_Type == "Copia"|Inflor_polyA$Repeat_Type =="Gypsy")

#Genomic features analysis
G1 <- lm(Log_GP ~ GC_ratio*Repeat_Type, data= Copia_Gypsy_G)
summary (G1)

G2 <- lm(Log_GP~ Similarity*Repeat_Type, data= Copia_Gypsy_G)
summary (G2)

G3 <- lm(GC_ratio~ Similarity*Repeat_Type, data= Copia_Gypsy_G)
summary (G3)

#TP analysis - applied for both ribo depleted data (with and without the outlier replicate)
m1 <- lm(Log_TP~ Tissue + Condition + Repeat_Type, data= Ribo_dep)
summary (m1)

m2 <- lm(Log_TP~ Tissue + Condition + Annotation, data= Ribo_dep)
summary (m2)

m3 <- lm(Log_TP~ (Similarity + GC_ratio)*Repeat_Type + Tissue + Condition, data= Copia_Gypsy_)
summary (m3)

m4 <- lm(Log_TP~ (Log_GP)*Repeat_Type + Tissue + Condition, data= Copia_Gypsy_)
summary (m4)

m5 <- lm(Log_TP~ Log_GP + Repeat_Type*Condition, data= Copia_Gypsy_)
summary (m5)

#TP analysis in mRNA libraries
m5 <- lm(Log_TP~ (Similarity + GC_ratio)*Repeat_Type, data= Copia_Gypsy_Leaf)
summary (m5)

m6 <- lm(Log_TP~ (Log_GP)*Repeat_Type, data= Copia_Gypsy_Leaf)
summary (m5)

m7 <- lm(Log_TP~ (Similarity + GC_ratio)*Repeat_Type, data= Copia_Gypsy_Inflor)
summary (m7)

m8 <- lm(Log_TP~ (Log_GP)*Repeat_Type, data= Copia_Gypsy_Inflor)
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

#Plot linear model graphs by Repeat_Type
library(ggplot2)
library(gridExtra)

#Ribodepleted libraries
p1 <- ggplot(data = Copia_Gypsy_, aes(Log_GP, Log_TP)) +
  geom_point(aes(shape= Repeat_Type, color = Repeat_Type, alpha= 1), size=2.5) +
  geom_smooth(aes(color = Repeat_Type), method = lm, 
              se = FALSE, fullrange = TRUE) +
  theme_bw() +
  ggtitle("Stem and Crown") +
    theme(strip.background = element_blank(), strip.placement = "outside") +
  facet_wrap(~Tissue + Condition) +
  scale_color_manual(values=c('purple', 'orange')) +
  scale_shape_manual(values=c(17,19))


#mRNA libraries
k1 <- ggplot(Copia_Gypsy_Leaf, aes(x=Log_GP, y=Log_TP, color=Repeat_Type, shape=Repeat_Type, alpha=1)) +
  geom_point(aes(size=Repeat_Type, alpha=1)) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  scale_color_manual(values=c('purple', 'orange')) +
  scale_shape_manual(values=c(17,19)) + 
  scale_size_manual(values=c(3.5,3.5)) +
  theme_bw()

k2 <- ggplot(Copia_Gypsy_Flor, aes(x=Log_GP, y=Log_TP, color=Repeat_Type, shape=Repeat_Type, alpha=1)) +
  geom_point(aes(size=Repeat_Type, alpha=1)) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  scale_color_manual(values=c('purple', 'orange')) +
  scale_shape_manual(values=c(17,19))  +
  scale_size_manual(values=c(3.5,3.5)) +
  theme_bw() 
