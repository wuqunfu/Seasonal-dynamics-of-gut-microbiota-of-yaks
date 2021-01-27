#Procrustes analysis
library("vegan")
library("ggplot2")
library("cowplot")
#set working environment
setwd("C:/Users/Fuyu Shi/Desktop/bar/man/procrustes_analysis")
getwd()
#bray-curtis similarity of diet and microbiome, Mantel test 
metadata = read.csv(file = "C:/Users/Fuyu Shi/Desktop/bar/man/EPIPHYTE_METADATA1.csv", header = TRUE, check.names = FALSE, row.names = 1)

df1 = read.csv(file = "C:/Users/Fuyu Shi/Desktop/bar/man/procrustes_analysis/grass_diet1.csv", header= TRUE,check.names = FALSE, row.names = 1)

df2 = read.csv(file = "C:/Users/Fuyu Shi/Desktop/bar/man/procrustes_analysis/grass_microbiome1.csv", header= TRUE,check.names = FALSE, row.names = 1)

#abundance data frame - bray curtis dissimilarity
dist.abund <- vegdist(df1, method = "bray")
dist.abund <- as.dist(dist.abund)
mdist.abund = vegdist(df2, method = "bray")
mdist.abund <- as.dist(mdist.abund)
#make pcoas
dpcoa <- as.data.frame(cmdscale(dist.abund)) 

mpcoa <- as.data.frame(cmdscale(mdist.abund))

#procrustes analysis
pro <- procrustes(X = dpcoa, Y = mpcoa, scale = TRUE,symmetric = TRUE)

pro_test <- protest(dpcoa,mpcoa,perm=9999)

eigen <- sqrt(pro$svd$d)
percent_var <- signif(eigen/sum(eigen), 4)*100

beta_pro <- data.frame(pro$X)
trans_pro <- data.frame(pro$Yrot)
beta_pro$UserName <- rownames(beta_pro)
beta_pro$type <- "Diet(Bray_curtis)"

seasons=metadata[,3]

beta_pro=cbind(beta_pro,seasons)

trans_pro$UserName <- rownames(trans_pro)
trans_pro$type <- "Microbiome"
seasons=metadata[,3]

trans_pro=cbind(trans_pro,seasons)
colnames(trans_pro) <- colnames(beta_pro)
pval <- signif(pro_test$signif, 1)
plot <- rbind(beta_pro, trans_pro)
col1 = rgb(250/255,60/255,60/255)
col4 = rgb(0/255,200/255,200/255)
col8 = rgb(160/255,0/255,200/255)
col10 = rgb(0/255,160/255,255/255)


grass_food_micro <- ggplot(plot) +
            geom_point(size = 4, alpha=0.75, aes(x = V1, y = V2, color = seasons,shape=type))+ 
            scale_color_manual(values = c(col8,col1,col10,col4)) +
            theme_classic() +
            scale_x_continuous(limits = c(-0.13,0.13))+
            scale_y_continuous(limits = c(-0.13,0.1))+
            geom_line(aes(x= V1, y=V2, group=UserName), col = "darkgrey", alpha = 0.6,size=0.2) +
            theme(panel.background = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             legend.title = element_blank(),
             legend.text = element_text(size=10,colour="black"),
             legend.position = 'bottom',
             axis.text = element_text(size=10,colour="black"),
             axis.title = element_text(size=13,colour="black"),
             aspect.ratio = 1) +
  guides(color = guide_legend(ncol = 1)) +
  annotate("text", x = 0.07, y = -0.13, label = paste0("p-value=",pval), size = 4) +
  xlab(paste0("PC 1 [",percent_var[1],"%]")) +
  ylab(paste0("PC 2 [",percent_var[2],"%]")) 

grass_food_micro_leg <- get_legend(grass_food_micro) 

grass_food_micro + theme(legend.position = "right")


