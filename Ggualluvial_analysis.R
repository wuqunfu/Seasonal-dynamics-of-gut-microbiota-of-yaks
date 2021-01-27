#ggualluvial
setwd('C:/Users/Fuyu Shi/Desktop/bar/ggalluvial/season/1')
library("reshape2", quietly=T, warn.conflicts=F)
library("ggalluvial")
library("ggplot2")
library("RColorBrewer")
main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   axis.line.x=element_line(size=1, colour="black"),
                   axis.line.y=element_line(size=1, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=18),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=13),
                   text=element_text(family="sans", size=18))
design = read.table("design1.txt", header=T, row.names= 1, sep="\t")
otu_table = read.delim("grass_otu_table.txt", row.names= 1,  header=T, sep="\t")
taxonomy = read.delim("rep_seqs_tax.txt", row.names= 1,header=F, sep="\t")
colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus","species")
idx = taxonomy$family == "f__"
taxonomy$full=as.character(taxonomy$family) 
taxonomy[idx,]$full=as.character(taxonomy[idx,]$family)
tax_count = merge(taxonomy, otu_table, by="row.names")
tax_count_sum = aggregate(tax_count[,-(1:9)], by=tax_count[9], FUN=sum) # mean
rownames(tax_count_sum) = tax_count_sum$full
tax_count_sum = tax_count_sum[,-1]
per = t(t(tax_count_sum)/colSums(tax_count_sum,na=T)) * 100 # normalization to total 100

mean_sort = per[(order(-rowSums(per))), ] # decrease sort
colSums(mean_sort)

mean_sort=as.data.frame(mean_sort)
mean_sort=mean_sort[rownames(mean_sort)!=" c__unknown",]
other = colSums(mean_sort[14:dim(mean_sort)[1], ])
mean_sort = mean_sort[1:(14-1), ]
mean_sort = rbind(mean_sort,other)
rownames(mean_sort)[14] = c("Others")

b=c("Others"," f__unknown")
mean_sort1=mean_sort[rownames(mean_sort)[!rownames(mean_sort)%in%b],]
100-apply(mean_sort1,2,sum)
others=100-apply(mean_sort1,2,sum)
others=t(as.data.frame(others))
mean_sort=rbind(mean_sort1,others)




write.table(mean_sort, file="Top14family_grass_RosFamily_14.txt", append = F, sep="\t", quote=F, row.names=T, col.names=T)

topN=rownames(mean_sort)


sub_design = subset(design,Season %in% c("Spring","Summer","Autumn","Winter") )


sub_design$group=sub_design$Season


sub_design$group  = factor(sub_design$group, levels=c("Spring","Summer","Autumn","Winter"))


print(paste("Number of group: ",length(unique(sub_design$group)),sep="")) # show group numbers


idx = rownames(sub_design) %in% colnames(mean_sort) 
sub_design = sub_design[idx,]
mean_sort = mean_sort[, rownames(sub_design)] # reorder according to design


mean_sort$family = rownames(mean_sort)

data_all = as.data.frame(melt(mean_sort, id.vars=c("family")))

data_all = merge(data_all, sub_design[c("group")], by.x="variable", by.y = "row.names")



p = ggplot(data_all, aes(x=variable, y = value, fill = family )) + 
  geom_bar(stat = "identity",position="fill", width=1)+ 
  scale_y_continuous(labels = scales::percent) + 
  facet_grid( ~ group, scales = "free_x", switch = "x") +  main_theme +
  theme(axis.ticks.x = element_blank(), legend.position="top", axis.text.x = element_blank(), strip.background = element_blank())+
  xlab("Groups")+ylab("Percentage (%)")           

p





mat = mean_sort[,1:(dim(mean_sort)[2]-1)]

mat_t = t(mat)

mat_t2 = merge(sub_design[c("group")], mat_t, by="row.names")

mat_t2 = mat_t2[,-1]

mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
mat_mean_final = do.call(rbind, mat_mean)[-1,]


geno = mat_mean$group
colnames(mat_mean_final) = geno
mat_mean_final = as.data.frame(mat_mean_final)
mat_mean_final$family = rownames(mat_mean_final)

# Table transform
data_all = as.data.frame(melt(mat_mean_final, id.vars=c("family")))
d1=data_all[data_all[,2]=="Spring",]
d11=d1[,2:3]
rownames(d11)=d1[,1]
d111=d11[order(-d11$value),]
d1111=data.frame(family=rownames(d111),d111)
d2=data_all[data_all[,2]=="Summer",]
d22=d2[,2:3]
rownames(d22)=d2[,1]
d222=d22[rownames(d111),]
d2222=data.frame(family=rownames(d111),d222)
d3=data_all[data_all[,2]=="Autumn",]
d33=d3[,2:3]
rownames(d33)=d3[,1]
d333=d33[rownames(d111),]
d3333=data.frame(family=rownames(d111),d333)
d4=data_all[data_all[,2]=="Winter",]
d44=d4[,2:3]
rownames(d44)=d4[,1]
d444=d44[rownames(d111),]
d4444=data.frame(family=rownames(d111),d444)

#data_all=data.frame(
#                    family=rep(d1111$family,4),
#                    variable=factor(
#                                   c(
#                                      as.character(d1111$variable),
#                                      as.character(d2222$variable),
#                                      as.character(d3333$variable),
#                                      as.character(d4444$variable)
#                                      )
#                                   ),
#                    value=c(d1111$value,d2222$value,d3333$value,d4444$value)
#                    )


data_all=rbind(d1111,d2222,d3333,d4444)
p = ggplot(data_all, aes(x=variable, y = value, fill = family )) + 
  geom_bar(stat = "identity",position="fill", width=0.7)+ 
  scale_y_continuous(labels = scales::percent) + 
  xlab("Groups")+ylab("Percentage (%)")+main_theme+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
p


Family=factor(data_all$family,level=rownames(d111))
#Variable=factor(as.character(data_all$variable))#,level=c("spring","summer","autumn","winter")
#ff=data_all$variable
#levels(ff)=c("spring","summer","autumn","winter")

colourCount = length(unique(mtcars$hp))
getPalette = colorRampPalette(brewer.pal(9, "Set3"))

col1 = rgb(250/255,60/255,60/255)
col2 = rgb(0/255,220/255,0/255)
col3 = rgb(30/255,60/255,255/255)
col4 = rgb(0/255,200/255,200/255)
col5 = rgb(240/255,0/255,130/255)
col6 = rgb(230/255,220/255,50/255)
col7 = rgb(240/255,130/255,40/255)
col8 = rgb(160/255,0/255,200/255)
col9 = rgb(160/255,230/255,50/255)
col10 = rgb(0/255,160/255,255/255)
col11 = rgb(230/255,175/255,45/255)
col12 = rgb(0/255,210/255,140/255)
col13 = rgb(247/255,153/255,209/255)

color <- c(col1,col2,col4,col6,col5,col3,col9,col7,col8,col10,col12,col11,col13)

p = ggplot(data = data_all, aes(x = data_all$variable, y = value,alluvium = Family), order(Polygonaceae,Rosaceae,Cyperaceae,Gramineae,Compositae,Ranunculaceae,Leguminosae,Salicaceae,Scrophulariaceae,Others)) +
  geom_alluvium(aes(fill = Family, colour = Family), alpha = .7,size=0) +
  main_theme + theme(axis.text.x = element_text( hjust = 0.5)) +
  ylab("Relative abundance (%)")+xlab("Season")+
  theme(axis.title=element_text(face="bold",size="10")) +
  ggtitle("Family changes among seasons under trans-humance management")+
  scale_fill_manual(values =color)
p