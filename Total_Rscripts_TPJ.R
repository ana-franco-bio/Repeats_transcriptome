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


#Plot simple histogram - Figure S1
Similarity_ <- Sim_dataOK$similarity
a1 <- hist(Similarity_, col = NULL, breaks = 40, xlim=c(85,100))

#Graphs Figure2
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


#Bar graphs TP - Figures 5 and 6
#load data
Data_<-read.table("Mean_TP.txt", header = TRUE, check.names = FALSE, sep = "\t")
attach(Data_)
Data1 <- as.data.frame(Data_)
str(Data1)
#select specific data
subStem <- subset(Data1, Tissue=="Stem")
subCrown <- subset(Data1, Tissue=="Crown")
subLeaf <- subset(Data1, Tissue=="Leaf")
subFlor <- subset(Data1, Tissue =="Inflorescence")

#Crown and Stem TP - ribo-depleted
#Crown
subCrown$Repeat_Type <-factor(subCrown$Repeat_Type,
                                  levels = c("Unknown", "LTR", "Copia", "Gypsy", "ClassII", "LINE", "Satellite", "Pararetrovirus"))
#Plot
TPok<- ggplot(subcrown_bar, aes(x = reorder(Repeat_Type, - TP),fill= Repeat_Type, y=TP)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Crown") +
  theme_bw()
#Set colors
TPok + scale_fill_manual(values=c("grey", "darkmagenta", "red", "darkblue", "lightsalmon", "darkgreen", "yellow", "black"))

#Stem
#Plot
TP_Stem<- ggplot(subStem, aes(x = reorder(Repeat_Type, - TP),fill= Repeat_Type, y=TP)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Stem")+
  theme_bw()
#Set colors
TP_Stem + scale_fill_manual(values=c("grey","red", "darkmagenta", "darkblue", "lightsalmon", "darkgreen" ,"yellow", "black"))

#Leaf and Inflorescence TP - poly-A
#Leaf 
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
#Test the significance
res <- cor.test(Genomic_data$GC_ratio, Genomic_data$Log_GP, 
                method = "pearson")
res


###Linear models###

#Load data
GP <-read.table("Genomic_features.txt",header = TRUE, check.names = FALSE, sep = "\t")
Ribo_dep <- read.table("All_repeats_ribodepleted.txt",header = TRUE, check.names = FALSE, sep = "\t")
Leaf_polyA <- read.table("Leaf_polyA.txt",header = TRUE, check.names = FALSE, sep = "\t")
Inflor_polyA <- read.table("Inflor_polyA.txt",header = TRUE, check.names = FALSE, sep = "\t")

#Preparing Data (do it for each library)
Ribo_dep$Cluster <- as.factor(Ribo_dep$Cluster)
Ribo_dep$Tissue <- as.factor(Ribo_dep$Tissue)
Ribo_dep$Repeat_Type <- as.factor(Ribo_dep$Repeat_Type)
Ribo_dep$Annotation <- as.factor(Ribo_dep$Annotation)
str(Ribo_dep)

#Subset Ty1/copia and Ty3/gypsy elements
Copia_Gypsy_Genomic <- subset(GP, GP$Repeat_Type == "Copia"|GP$Repeat_Type =="Gypsy")
Transposons_Genomic <- subset(GP, GP$Repeat_Type == "Class_II")
Copia_Gypsy_Ribo <-subset(Ribo_dep, Ribo_dep$Repeat_Type == "Copia"|Ribo_dep$Repeat_Type =="Gypsy")
Transposons_Ribo <- subset(Ribo_dep, Ribo_dep$Repeat_Type == "Class_II")
Copia_Gypsy_Leaf <-subset(Leaf_polyA, Leaf_polyA$Repeat_Type == "Copia"|Leaf_polyA$Repeat_Type =="Gypsy")
Transposons_Leaf <-subset(Leaf_polyA, Leaf_polyA$Repeat_Type == "Class_II")
Copia_Gypsy_Inflor <-subset(Inflor_polyA, Inflor_polyA$Repeat_Type == "Copia"|Inflor_polyA$Repeat_Type =="Gypsy")
Transposons_Inflor <-subset(Inflor_polyA, Inflor_polyA$Repeat_Type == "Class_II")

#Genomic features analysis
G1 <- lm(Log_GP ~ GC_ratio) #applied to all subsets separately
summary (G1)

G2 <- lm(Log_GP ~ Similarity, data= Copia_G1) #applied to all subsets separately
summary (G2)

#TP analysis - applied for both ribo depleted and Poly-A data
m1 <- lm(Log_TP~ Log_GP + Tissue) #applied to all subsets separately
summary (m1)

m2 <- lm(Log_TP~ Similarity + GC_ratio + Tissue)#applied to all subsets separately
summary (m2)


##Linear model graphs- Genomic features
d1 <- ggplot(Copia_Gypsy_Genomic, aes(x=GC_ratio, y=Log_GP, color=Repeat_Type, shape=Repeat_Type, alpha=1)) +
  geom_point(aes(size=Repeat_Type, alpha=1)) + 
  geom_smooth(data= subset(Copia_Gypsy_G, Repeat_Type == "Gypsy"),method=lm, se=FALSE, fullrange=TRUE)+
  scale_color_manual(values=c('purple', 'orange')) +
  scale_shape_manual(values=c(17,19))  +
  scale_size_manual(values=c(3.5,3.5)) +
  theme_bw()
d1
d2 <- ggplot(Copia_Gypsy_Genomic, aes(x=Similarity, y=Log_GP, color=Repeat_Type, shape=Repeat_Type, alpha=1)) +
  geom_point(aes(size=Repeat_Type, alpha=1)) + 
  geom_smooth(data= subset(Copia_Gypsy_G, Repeat_Type == "Copia"),method=lm, se=FALSE, fullrange=TRUE)+
  scale_color_manual(values=c('purple', 'orange')) +
  scale_shape_manual(values=c(17,19))  +
  scale_size_manual(values=c(3.5,3.5)) +
  theme_bw()
d2
#Plot all graphs together
library(ggpmisc)
install.packages("Rcpp")
library(Rcpp)
library(ggpubr)
install.packages("ggpmisc")
figure_d <- ggarrange(d1, d2, 
                      labels = c("a", "b", size=1),
                      ncol = 2,
                      common.legend = TRUE,legend="bottom",
                      heights = c(1,1)) 
figure_d

##Linear models graphs - genomic and transcriptome features
All_tissues <- read.table("clipboard", header=TRUE)
All_tissues$Cluster <- as.factor(All_tissues$Cluster)
All_tissuesr$Tissue <- as.factor(All_tissues$Tissue)
All_tissues$Repeat_Type <- as.factor(All_tissues$Repeat_Type)

Copia_Gypsy<-subset(All_tissues, All_tissues$Repeat_Type == "Copia"|All_tissues$Repeat_Type =="Gypsy")

Scatter1 <- ggplot(data = Copia_Gypsy, aes(Log_GP, Log_TP)) +
  geom_point(aes(shape= Repeat_Type, color = Repeat_Type, alpha= 1), size=2.5) +
  geom_smooth(aes(color = Repeat_Type), method = lm, 
              se = FALSE, fullrange = TRUE) +
  theme_bw() +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  facet_wrap(~Tissue) +
  scale_color_manual(values=c('purple', 'orange')) +
  scale_shape_manual(values=c(17,19))

Scatter2 <- ggplot(data = Copia_Gypsy, aes(GC_ratio, Log_TP)) +
  geom_point(aes(shape= Repeat_Type, color = Repeat_Type, alpha= 1), size=2.5) +
  geom_smooth(aes(color = Repeat_Type), method = lm, 
              se = FALSE, fullrange = TRUE) +
  theme_bw() +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  facet_wrap(~Tissue) +
  scale_color_manual(values=c('purple', 'orange')) +
  scale_shape_manual(values=c(17,19))

Scatter3 <- ggplot(data = Copia_Gypsy, aes(Similarity, Log_TP)) +
  geom_point(aes(shape= Repeat_Type, color = Repeat_Type, alpha= 1), size=2.5) +
  geom_smooth(aes(color = Repeat_Type), method = lm, 
              se = FALSE, fullrange = TRUE) +
  theme_bw() +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  facet_wrap(~Tissue) +
  scale_color_manual(values=c('purple', 'orange')) +
  scale_shape_manual(values=c(17,19))



 