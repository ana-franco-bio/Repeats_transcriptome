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

Sim_dataOK
#Plot simple histogram - Figure S2
Similarity_ <- Sim_dataOK$similarity
a1 <- hist(Similarity_, col = NULL, breaks = 40, xlim=c(85,100))
moda <- names(sort(-table(Similarity_)))[1]

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


#Bar graphs TP - Figure 6
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


###Linear models###

#Load data
GP <-read.table("Genomic_features_final.txt",header = TRUE, check.names = FALSE, sep = "\t")
Ribo_dep <- read.table("All_repeats_ribodepleted.txt",header = TRUE, check.names = FALSE, sep = "\t")
Leaf_polyA <- read.table("Leaf_polyA.txt",header = TRUE, check.names = FALSE, sep = "\t")
Inflor_polyA <- read.table("Inflor_polyA.txt",header = TRUE, check.names = FALSE, sep = "\t")

#Preparing Data (do it for each library)
Ribo_dep$Cluster <- as.factor(Ribo_dep$Cluster)
Ribo_dep$Tissue <- as.factor(Ribo_dep$Tissue)
Ribo_dep$Repeat_Type <- as.factor(Ribo_dep$Repeat_Type)
Ribo_dep$Annotation <- as.factor(Ribo_dep$Annotation)
str(Ribo_dep)

GP$Repeat_Type <- as.factor(GP$Repeat_Type)
GPok <- as.data.frame(GP)
#Subset Ty1/copia and Ty3/gypsy elements
Copia_Gypsy_Genomic <- subset(GPok, GPok$Repeat_Type == "Copia"|GPok$Repeat_Type =="Gypsy")
Transposons_Genomic <- subset(GPok, GPok$Repeat_Type == "ClassII")
Gypsy_genomic <- subset(GPok, GPok$Repeat_Type=="Gypsy")
Copia_genomic <- subset(GPok, GPok$Repeat_Type=="Copia")
Copia_Gypsy_Ribo <-subset(Ribo_dep, Ribo_dep$Repeat_Type == "Copia"|Ribo_dep$Repeat_Type =="Gypsy")
Transposons_Ribo <- subset(Ribo_dep, Ribo_dep$Repeat_Type == "Class_II")
Copia_Gypsy_Leaf <-subset(Leaf_polyA, Leaf_polyA$Repeat_Type == "Copia"|Leaf_polyA$Repeat_Type =="Gypsy")
Transposons_Leaf <-subset(Leaf_polyA, Leaf_polyA$Repeat_Type == "Class_II")
Copia_Gypsy_Inflor <-subset(Inflor_polyA, Inflor_polyA$Repeat_Type == "Copia"|Inflor_polyA$Repeat_Type =="Gypsy")
Transposons_Inflor <-subset(Inflor_polyA, Inflor_polyA$Repeat_Type == "Class_II")

#Genomic features analysis
Copia2G <- read.table("clipboard", header =TRUE)
Genomic_1 <- lm(Log_GP ~ GC_ratio, data = Transposons_Genomic) #applied to all subsets separately
summary(Genomic_1)

Genomic_2 <- lm(Log_GP ~ GC_ratio, data = Copia_genomic) #applied to all subsets separately
summary(Genomic_2)

Genomic_2 <- lm(Log_GP ~ Similarity, data = Gypsy_genomic) #applied to all subsets separately
summary(Genomic_2)

Genomic_2 <- lm(Log_GP ~ Similarity, data = Copia_genomic) #applied to all subsets separately
summary(Genomic_2)

Genomic3 <- lm(GC_ratio ~ Similarity*Repeat_Type, data= Copia_Gypsy_Genomic)
summary(Genomic3)
#TP analysis - applied for both ribo depleted and Poly-A data
m1 <- lm(Log_TP~ Log_GP + Tissue) #applied to all subsets separately
summary (m1)

m2 <- lm(Log_TP~ Similarity + GC_ratio + Tissue)#applied to all subsets separately
summary (m2)


##Linear model graphs- Genomic features
library(ggplot2)
d1 <- ggplot(Copia_Gypsy_Genomic, aes(x=GC_ratio, y=Log_GP, color=Repeat_Type, shape=Repeat_Type, alpha=1)) +
  geom_point(aes(size=Repeat_Type), alpha =1) + 
  geom_smooth(data= subset(Copia_Gypsy_Genomic, Repeat_Type == "Gypsy"),method=lm, se=FALSE, fullrange=TRUE)+
  scale_color_manual(values=c('orange', 'purple')) +
  scale_shape_manual(values=c(21,24))  +
  scale_size_manual(values=c(3.5,3.5)) +
  theme_bw()
d1
d1 <- ggplot(Copia_Gypsy_Genomic, aes(x=GC_ratio, y=Log_GP, shape=Repeat_Type, group =Repeat_Type)) +
  geom_point(aes(size=Repeat_Type, fill= Repeat_Type)) + 
  geom_smooth(data= subset(Copia_Gypsy_Genomic, Repeat_Type == "Gypsy"),method=lm, se=FALSE, aes(color=Repeat_Type), fullrange=TRUE)+
  scale_color_manual(values=c('orange', 'purple')) +
  scale_fill_manual(values=c('orange', 'purple'))+
  scale_shape_manual(values=c(21,24))  +
  scale_size_manual(values=c(3.5,3.5)) +
  theme_bw()
d1

d2 <- ggplot(Copia_Gypsy_Genomic, aes(x=Similarity, y=Log_GP, color=Repeat_Type, shape=Repeat_Type, alpha=1)) +
  geom_point(aes(size=Repeat_Type, alpha=1)) + 
  geom_smooth(data= subset(Copia_Gypsy_Genomic, Repeat_Type == "Copia"),method=lm, se=FALSE, fullrange=TRUE)+
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

#Subset Class II, Ty1/copia and Ty3/gypsy elements
Copia_Setaria_F <-subset(Setaria_Flor2023, Setaria_Flor2023$Repeat_Type == "Copia")
Gypsy_Setaria_F <- subset (Setaria_Flor2023, Setaria_Flor2023$Repeat_Type =="Agypsy")
Transposons_Setaria_F <-subset(Setaria_Flor2023, Setaria_Flor2023$Repeat_Type == "ClassII")
Retransposons_Setaria_F <-subset(Setaria_Flor2023, Setaria_Flor2023$Repeat_Type == "Copia"|Setaria_Flor2023$Repeat_Type =="Agypsy")

#TP analysis - applied for both ribo depleted and Poly-A data
F1 <- lm(Log_TP~ Log_GP, data = Copia_Setaria_F) #Copia Elements
summary (F1)


#For Gypsy elements
F2 <- lm(Log_TP~ Log_GP, data = Gypsy_Setaria_F) #Gypsy Elements
summary (F2)

# For transposons Class II
F3 <- lm(Log_TP~ Log_GP, data = Transposons_Setaria_F) #Class II Elements
summary (F3)

#TP with similarity and GC ratio
#For Copia elements
F4Sim <- lm(Log_TP~ Similarity + GC_ratio, data = Copia_Setaria_F)#applied to all subsets separately
summary (F4Sim)

#For gypsy
F5 <- lm(Log_TP~ Similarity + GC_ratio, data = Gypsy_Setaria_F)#applied to all subsets separately
summary (F5)

#For Transposons
F6 <- lm(Log_TP~ Similarity + GC_ratio, data = Transposons_Setaria_F)#applied to all subsets separately
summary (F6)

####################
##Linear models graphs - genomic and transcriptome features
All_tissues <- read.table("clipboard", header=TRUE)
All_tissues$Cluster <- as.factor(All_tissues$Cluster)
All_tissuesr$Tissue <- as.factor(All_tissues$Tissue)
All_tissues$Repeat_Type <- as.factor(All_tissues$Repeat_Type)

#Crown and Stem
Scatter_newSetaria4 <- ggplot(Gypsy_Setaria_SC, aes(x= GC_ratio, y=  Log_TP, group = Tissue)) +
  geom_point(aes( shape = Tissue, fill = Tissue), size=3.0) +
  geom_smooth(aes(color = Tissue), method = lm, 
              se = FALSE, fullrange = TRUE)+
  theme_bw() +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  scale_fill_manual(values=c('salmon', 'darkgrey'))+
  scale_color_manual(values=c('salmon', 'darkgrey')) +
  scale_shape_manual(values=c(21,21))
Scatter_newSetaria4

#Inflorescence graphs
Scatter_flor <- ggplot(Transposons_Setaria_F, aes(x= GC_ratio, y=  Log_TP)) +
  geom_point(color= 'black', fill = 'red', shape = 21,size=3.5) +
  geom_smooth(method = lm, color = 'red', 
              se = FALSE, fullrange = TRUE)+
  theme_bw() +
  theme(strip.background = element_blank(), strip.placement = "outside") 
Scatter_flor
