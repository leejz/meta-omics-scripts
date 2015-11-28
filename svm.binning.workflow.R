### SVM binning workflow R script v.1.0 ############################################################
#based on the metagenomics workflow of Albertsen et al. (2013) and updated with binning and visualizaiton code

# Install the required packages
install.packages("vegan")           #Used for PCA/CA analysis
install.packages("RColorBrewer")    #Easy color selection
install.packages("alphahull")       #Generation of convex spaces for extraction of points
install.packages("fpc")             #library for DBSCAN
install.packages("e1071")           #library for SVM
install.packages("fields")          #for drawing splines
install.packages("rgl")             #Plot3D functionality

# Load the required packages
library(vegan)
library(RColorBrewer)
library(alphahull)
library(fpc)
library(e1071)
library(fields)
library(rgl)

### Read and prepare data ##########################################################################
# Read Data - all data have been generated using previous steps - R is used to combine them

Coassembly_coverage <- read.csv("k29.Scaffolds.All.Counts.txt", header = F, sep = "\t", quote = "")               
colnames(Coassembly_coverage) = c("Contig.Name","Reference.length","Coverage_1","Coverage_2","Coverage_3","Coverage_4")
gc <- read.delim("scaffolds.gc.tab", header = T)
kmer <- read.delim("scaffolds.kmer.tab", header = T)
colnames(kmer)[1] = "name"
ess <- read.table("scaffolds.assembly.800.orfs.hmm.id.txt", header = F)
colnames(ess) = c("name","orf","hmm.id")
ess.tax <- read.delim("scaffolds.assembly.800.orfs.hmm.blat.tax.tab", header = F) 
colnames(ess.tax) = c("name","orf","phylum")
cons.tax <- read.delim("scaffolds.assembly.800.tax.consensus.txt", header = T)
colnames(cons.tax) = c("name","phylum","tax.color","all.assignments")

# reduced filtered df, if needed
filter_cutoff = 1500
filtered_df <- Coassembly_coverage[Coassembly_coverage$Reference.length > filter_cutoff,]
filtered_df <- merge(filtered_df, gc, by.x ="Contig.Name", by.y="contig")

#PCA based analysis
#log coverage PCA
data <- cbind(log10(filtered_df[c("Coverage_1","Coverage_2","Coverage_3","Coverage_4")]+1))
pca_mat <- cbind(t(data))

pca_results <- prcomp(pca_mat)
print(pca_results)
summary(pca_results)
plot(pca_results)
pca_results

PC1 <- data.frame(filtered_df$Contig.Name, pca_results$rotation[,1]+1,filtered_df$Reference.length)
PC2  <- data.frame(filtered_df$Contig.Name, pca_results$rotation[,2]+1,filtered_df$Reference.length)
PC3  <- data.frame(filtered_df$Contig.Name, pca_results$rotation[,3]+1,filtered_df$Reference.length)
PC4  <- data.frame(filtered_df$Contig.Name, pca_results$rotation[,4]+1,filtered_df$Reference.length)
colnames(PC1) = c("Name","Loading","Reference.length")
colnames(PC2) = c("Name","Loading","Reference.length")
colnames(PC3) = c("Name","Loading","Reference.length")
colnames(PC4) = c("Name","Loading","Reference.length")

# Combine all data to one data matrix "d"
d <- data.frame(PC1$Name,     
                PC1$Reference.length,
                filtered_df$gc,
                PC1$Loading*100,
                PC2$Loading*100,
                PC3$Loading*100
                )

colnames(d) = c("name",
               "length",
               "gc",
               "HPminus",
               "HPplus",
               "HP3"
                )

d <- merge(d,cons.tax, by = "name", all.x = T)
# Combine data on essential genes 
ess<- merge(ess, d, by = "name", all.x = T)
ess<- merge(ess, ess.tax, by = c("name","orf"), all.x = T)

### Preliminary coverage v. length analysis ##########################################################################

pre_d <- merge(d[c("name","tax.color","phylum","length")], Coassembly_coverage[c("Contig.Name", "Coverage_1","Coverage_2","Coverage_3","Coverage_4")], by.x="name",by.y="Contig.Name")
pre_d$avg <- rowMeans(pre_d[c("Coverage_1","Coverage_2","Coverage_3","Coverage_4")])

byor<-colorRampPalette(c("blue","orange","red"))
phspan <- length(unique(pre_d$tax.color))-1
byor_key <- arrange(pre_d[!duplicated(factor(pre_d$tax.color)),],tax.color)[1:phspan,"tax.color"]
palette(adjustcolor(byor(phspan), alpha.f = 1))

plot(x = pre_d$length,
     y = pre_d$avg,
     ylim = c(2, 250),
     xlim = c(1500,200000),
     log = "yx",
     cex = .8,
     pch=pre_d$tax.color,
     col=match(pre_d$tax.color,byor_key),
     xlab = "Contig Length",
     ylab = "Avg Contig Coverage"     
)

points(x = pre_d$length[is.na(pre_d$tax.color)],
       y = pre_d$avg[is.na(pre_d$tax.color)],
       cex = .4,
       pch=20,
       col=rgb(0,0,0,0.1),
)

# Add phylum affiliation legend
legend(x = 2e5,
       y = 200,       
       arrange(pre_d[!duplicated(factor(pre_d$tax.color)),],tax.color)[1:phspan,"phylum"],
       bty="n",
       pch=arrange(pre_d[!duplicated(factor(pre_d$tax.color)),],tax.color)[1:phspan,"tax.color"],
       col=1:length(arrange(pre_d[!duplicated(factor(pre_d$tax.color)),],tax.color)[1:phspan,"tax.color"]),
       pt.cex=2,
       xjust=1,
       y.intersp=.5
)

### split off training dataset based on length cutoff ##########################################################################

cutoff = 10000
d2 <- subset(d,length > cutoff)
ess2 = subset(ess,length > cutoff)

### TN PCA of cutoff testing  ##########################################################################

#TN_data <- kmer[match(d4$name, kmer$name),]

#rownames(TN_data) = TN_data$name
#TN_data$name = {}
#TN_pca_mat <- cbind(t(TN_data))
#TN_pca_results <- prcomp(TN_pca_mat)
#print(TN_pca_results)
#plot(TN_pca_results)
#summary(TN_pca_results)
#plot(TN_pca_results$rotation[,1:2],
#     pch=20,
#     cex=.7,
#     #xlim = c(-10, 10), 
#     #ylim = c(-20, 20),
#     #xlim = c(-4, 4), 
#     #ylim = c(-5, 5),
#     #col = d2$tax.color[d2$db_cluster==9]
#     col = d2$gc[d2$db_cluster==9]-min(d2$gc)
#)

#plot(TN_pca_results$rotation[,1:2],
#     pch=20,
#     cex=.7,
#     #xlim = c(-10, 10), 
#     #ylim = c(-20, 20),
#     #xlim = c(-4, 4), 
#     #ylim = c(-5, 5),
#     col = d4$tax.color
#     #col = d4$gc-min(d2$gc)
#)
#TN_cca <- cca(TN_data)
#plot(scores(TN_cca)$sites,
#     pch=20,
#     cex=.7,
#     #xlim = c(-10, 10), 
#     #ylim = c(-20, 20),
#     xlim = c(-4, 4), 
#     ylim = c(-5, 5),
#     #col = d2$tax.color[d2$db_cluster==9]
#     col = d2$gc[d2$db_cluster==9]-min(d2$gc)
#)

#TN PCA of cutoff 
#TN_data <- kmer[match(d2$name, kmer$name),]
#rownames(TN_data) = TN_data$name
#TN_data$name = {}
#TN_pca_mat <- cbind(t(TN_data))
#TN_pca_results <- prcomp(TN_pca_mat)
#print(TN_pca_results)
#plot(TN_pca_results)

#TN_cca <- cca(d2[d2$db_cluster == 9,c("HPplus", "HPminus", "HP3")],TN_pca_results$rotation[,1])
#plot(TN_cca)
# plot TN results with phylum
#palette(brewer.pal(11,"Paired"))
#plot(scores(TN_cca)$sites,
#     pch=20,
#     cex=.7,
#     #xlim = c(-10, 10), 
#     #ylim = c(-20, 20),
#     xlim = c(-4, 4), 
#     ylim = c(-5, 5),
#     col = d2$tax.color[d2$db_cluster==9]
#     #col = d2$gc[d2$db_cluster==9]-min(d2$gc)
#)

#points(scores(TN_cca)$sites,
#       pch=20,
#       cex=1,
#       col = d2$tax.color
#)

#legend(x = -4,
#       y = 10,
#       arrange(d2[!duplicated(factor(d2$tax.color)),],tax.color)[1:11,"phylum"],
#       bty="n",
#       pch=20,
#       col=1:11,
#       pt.cex=4,
#       xjust=1,
#       y.intersp=1.2
#)

#gbr<-colorRampPalette(c("green","blue","orange","red"))
#gcspan <- round(max(d2$gc)-min(d2$gc)+1)
#palette(adjustcolor(gbr(gcspan), alpha.f = 0.5))
#points(scores(TN_cca)$sites,
#       pch=20,
#       cex=1,
#       col = d2$gc-min(d2$gc)
#)

### Clustering by DBSCAN #######################################################################################

#DBSCAN analysis
filtered_data <- data.frame(PC1$Loading*1000, PC2$Loading*1000, PC3$Loading*1000, d$gc*.1, d$length)
colnames(filtered_data) <- c("PC1","PC2","PC3","gc","length")
filtered_data = subset(filtered_data, length > cutoff)
filtered_data$length = NULL
  
db_results <- dbscan(filtered_data, eps=.35, MinPts=4)
#db_results <- dbscan(filtered_data, eps=.45, MinPts=4)
db_results
plot(db_results, filtered_data)

d2$db_cluster = db_results$cluster
ess2 <- merge(ess2, d2[c("db_cluster","name")], by="name",All.x=T)
unique_cluster = unique(d2$db_cluster)

dev.off()
# Define a good color palette and make it transparrent
gbr<-colorRampPalette(c("green","blue","orange","red"))
gcspan <- round(max(d$gc)-min(d$gc)+1)
palette(adjustcolor(gbr(70), alpha.f = 0.1))

dbspan <- length(unique_cluster)
palette(adjustcolor(rainbow(dbspan), alpha.f = 1))

#dbspan <- round(max(d2$db_cluster)-min(d2$db_cluster)+1)
#palette(adjustcolor(rainbow(dbspan), alpha.f = 1))
#palette(brewer.pal(11,"Paired"))
#gbr<-colorRampPalette(c("green","blue","orange","red"))
#gcspan <- round(max(d2$gc)-min(d2$gc)+1)
#palette(adjustcolor(gbr(gcspan), alpha.f = 0.5))

#Visualize PC1-3 in 3D

plot3d(d2$HPminus[d2$db_cluster > 0], 
       d2$HPplus[d2$db_cluster > 0], 
       d2$HP3[d2$db_cluster > 0], 
       col=d2$db_cluster[d2$db_cluster > 0],        
       size=3, 
       type='p',
       xlab="PC1",
       ylab="PC2",
       zlab="PC3")

points3d(d2$HPminus[d2$db_cluster < 1], 
       d2$HPplus[d2$db_cluster < 1], 
       d2$HP3[d2$db_cluster < 1])

bintext <-{}
# Add a label in the center of each extracted bin
for (i in 1:length(unique_cluster)) {
  bintext <- rbind(bintext, c(as.character(unique_cluster[i]),
                              mean(d2$HPminus[d2$db_cluster ==unique_cluster[i]]), 
                              mean(d2$HPplus[d2$db_cluster == unique_cluster[i]]), 
                              mean(d2$HP3[d2$db_cluster == unique_cluster[i]])))
}
text3d(x=bintext[,2], y=bintext[,3], z=bintext[,4],text = bintext[,1])

#record a movie, if needed
M <- par3d("userMatrix") 
play3d( spin3d(c(0,1,0), rpm=2 ), duration=30) 
movie3d( spin3d(c(1,0,0), rpm=3 ), duration=20, dir="./raw\ movie", convert=FALSE, top= TRUE, clean=FALSE, type="gif", startTime=0)

#spot check point
scaffold_name = "scaffold-5988"
points3d(x=d3$HPplus[d3$name == scaffold_name], y=d3$HPminus[d3$name == scaffold_name], z=d3$HP3[d3$name == scaffold_name], size=15,col=9)
text3d(x=d3$HPplus[d3$name == scaffold_name], y=d3$HPminus[d3$name == scaffold_name], z=d3$HP3[d3$name == scaffold_name],text=scaffold_name)

#View ess only
#plot3d(ess2$HPminus[!is.na(ess2$tax.color)], 
#       ess2$HPplus[!is.na(ess2$tax.color)], 
#       ess2$HP3[!is.na(ess2$tax.color)], 
#d2$HPminus[!is.na(d2$tax.color)], 
#d2$HPplus[!is.na(d2$tax.color)], 
#d2$HP3[!is.na(d2$tax.color)], 
#col=d2$db_cluster[!is.na(d2$tax.color)]+1, 
#       col=ess2$db_cluster[!is.na(ess2$tax.color)]+1, 
#       size=5, 
#       type='p',
#       xlab="PC1",
#       ylab="PC2",
#       zlab="PC3")
#plot3d(x=data$Coverage_1,y=data$Coverage_2,z=data$Coverage_3,size=3,col=d$gc-min(d$gc)+1,type='p')

#Visualize dbscan clusters in 2D
plot(x = d2$HPminus,
     y = d2$HPplus,
     xlim = c(99.3, 101.2),
     ylim = c(98.4, 101.8),    
     cex = .8,
     pch=20,
     col=d2$db_cluster,     
     xlab = "PC1",
     ylab = "PC2"     
)

  borderpts <- subset(d2, db_cluster == 0)
  points(x = borderpts$HPminus,
         y = borderpts$HPplus,
         pch=20,
         cex = sqrt(d_prot$length)/100,
         col=rgb(0,0,0,1)
         lwd=.8
  )

# Add a label in the center of each extracted bin

for (i in 1:length(unique_cluster)) {
       text(mean(d2$HPminus[d2$db_cluster ==unique_cluster[i]]),
              mean(d2$HPplus[d2$db_cluster == unique_cluster[i]]),
              labels = as.character(unique_cluster[i]),
              cex=1,
              font=1
       )
}

#optional visualization code

#plot the %gc of a selected cluster
#d_prot <- subset(d2, db_cluster == 9)
#palette(adjustcolor(rainbow(max(d_prot$gc)-min(d_prot$gc)), alpha.f = 0.5))
#points(x = d_prot$HPminus,
#       y = d_prot$HPplus,
#       pch=20,
#       cex = 2,
#       col=d_prot$gc-min(d_prot$gc)
#)

#spot check clusters
#d_check_clust <- subset(d2, db_cluster == 35)
#hist(d_check_clust$gc,breaks=40,xlim=c(25,75))
#d_prot1 <- subset(d_check_clust, gc < 68.5)
#d_prot2 <- subset(d_check_clust, gc > 68.5)

#cross reference phylum with bins
table(d2$phylum,d2$db_cluster)

#optional bin manipulation code

# Manually join bins
#view %gc and join 15-501
#hist(d2$gc[d2$db_cluster == 501],breaks=20,xlim=c(25,70))
#hist(d2$gc[d2$db_cluster == 15],breaks=20,xlim=c(25,70))
#d2$db_cluster[d2$db_cluster == 15 | d2$db_cluster == 501] = 15501

#Manually split bins

#split bin 5 into 500 and 501
# Interactively chose 6 points on the plot that include the subset of scaffolds targeted
#z.def <- ahull(locator(6, type="p", pch=20), alpha=100000)                 
# Mark the defined subset on the plot
#plot(z.def,add=T, col="black")
#sourcebin = 5
# Extract the scaffolds (gcent.out) and essential genes (ecent.out) in the defined subset
#gcent.out <- {}
#for (i in 1:nrow(d2)) { if (inahull(z.def, c(d2$HPminus[i],d2$HPplus[i]))) gcent.out <- rbind(gcent.out,d2[i,])}
#gcent.out <- subset(gcent.out, db_cluster == sourcebin)
#g_exclude <- subset(d2, db_cluster == sourcebin & !(name %in% gcent.out$name))
#hist(gcent.out$gc,breaks=40,xlim=c(30,40))
#hist(g_exclude$gc,breaks=40,xlim=c(30,40))
#d2$db_cluster[d2$name %in% gcent.out$name] = 500
#d2$db_cluster[d2$name %in% g_exclude$name] = 501
#manual gc filtering
#d2$db_cluster[d2$name %in% g_exclude$name & d2$gc < 31] = 0

#manual gc filtering direct to dataframe, keeps the same binnum
#filter bin 9 for sequences > 68.5 %GC
#sourcebin = 9
#remove members with %GC threshold
#d2$db_cluster[d2$db_cluster == sourcebin & d2$gc > 68.5] = 0
#view result
#hist(d2$gc[d2$db_cluster == sourcebin],breaks=20,xlim=c(65,75))
#d2[d2$db_cluster == sourcebin,c("gc","phylum")] 

#manually create bin from graph
# Interactively chose 6 points on the plot that include the subset of scaffolds targeted
#z.def <- ahull(locator(6, type="p", pch=20), alpha=100000)                 
# Mark the defined subset on the plot
#plot(z.def,add=T, col="black")
#sourcebin = 7420
# Extract the scaffolds (gcent.out) and essential genes (ecent.out) in the defined subset
#gcent.out <- {}
#for (i in 1:nrow(d2)) { if (inahull(z.def, c(d2$HPminus[i],d2$HPplus[i]))) gcent.out <- rbind(gcent.out,d2[i,])}
#hist(gcent.out$gc,breaks=40,xlim=c(40,50))
#d2$db_cluster[d2$name %in% gcent.out$name] = sourcebin


### # Display galaxy chart with training contigs up front #######################################################################################

#get top 11 phyla/class  BROKEN
freqtab2 <- as.data.frame(margin.table(table(d3[c("tax.color","phylum")]),1))
colorlist <- freqtab[with(freqtab,order(-Freq)),]
topcolors <- colorlist$tax.color[1:11]
freqtab <- as.data.frame(margin.table(table(d3[c("tax.color","phylum")]),2))
phylumlist <- freqtab[with(freqtab,order(-Freq)),]
topphylum <- phylumlist$phylum[1:11]

plot(x = d2$HPminus,
     y = d2$HPplus,
     xlim = c(99.3, 101.2),
     ylim = c(98.4, 101.8),
     cex = .3,
     pch=20,
     xlab = "PC1",
     ylab = "PC2"     
)

palette(brewer.pal(11,"Paired"))
points(x = d$HPminus[d$tax.color %in% topcolors],
       y = d$HPplus[d$tax.color %in% topcolors],
       cex = sqrt(d$length[d$tax.color %in% topcolors])/100,
       col=d$tax.color[d$tax.color %in% topcolors],
       lwd=.8
)
points(x = d$HPminus,
       y = d$HPplus,
       cex = .2,
       col=rgb(0,0,0,0.05),
       lwd=1,
)

points(x = d$HPminus[d$name %in% esfc],
       y = d$HPplus[d$name %in% esfc],
       cex = 2,
       col=rgb(1,1,0,1),
       pch=20,
       lwd=1,
)

# Add scaffold size legend
t<-as.character(c(100,50,10,1))
legend(x = 99.35,
       y = 101.7,
       t,
       bty="n",
       pch=20,
       col=rgb(0,0,0,0.2),
       pt.cex=sqrt(as.integer(t)*1000)/100,
       x.intersp=4,
       y.intersp=1,
       adj=c(0.5,0.5),
       title="Length (kbp)"
)

# Add phylum affiliation legend
legend(x = 101.3,
       y = 101.7,
       topphylum,
       bty="n",
       pch=20,
       col=topcolors,
       pt.cex=2,
       xjust=1,
       y.intersp=.5
)

#PC based coloring to visualize long vs all contigs
plot(x = d$HPminus,
     y = d$HPplus,
     xlim = c(99, 102),
     ylim = c(99.5, 101),
     cex = sqrt(d$length)/200,
     pch=20,
     col=rgb(0,0,1,0.05),
     xlab = "PC1",
     ylab = "PC2"     
)

points(x = d2$HPminus,
       y = d2$HPplus,
       cex = .8,
       pch=20,
)

palette(brewer.pal(11,"Paired"))
points(x = d$HPminus[d$tax.color<12],
       y = d$HPplus[d$tax.color<12],
       cex = sqrt(d$length[d$tax.color<12])/100*0.7,
       col=d$tax.color[d$tax.color<12]      
)

#GC coloring at same scale (optional)
#gbr<-colorRampPalette(c("green","blue","orange","red"))
#gcspan <- round(max(d$gc)-min(d$gc)+1)
#palette(adjustcolor(gbr(gcspan), alpha.f = 0.1))
#plot(x = d$HPminus,
#     y = d$HPplus,
#     xlim = c(99, 102),
#     ylim = c(99.5, 101),
#     cex = sqrt(d$length)/200,
#     pch=20,
#     col=d$gc-min(d$gc),
#     xlab = "PC1",
#     ylab = "PC2"     
#)

### SVM classification of remaining data #####################################################

#use this as the training set for an SVM to train and classify the filtered_df dataset

trainset <- d2[c("name","HPplus","HPminus","HP3","gc","db_cluster")]
trainset <- rename(trainset, replace=c("db_cluster" = "Category"))
trainset <- subset(trainset, Category >= 0)
index <- 1:nrow(trainset)
testindex <- sample(index, trunc(length(index)/3))
testset <- trainset[testindex,]

## svm testing, for determination of SVM parameters only
svm.model <- svm(Category~., data = trainset[2:6], gamma = 1, cost = 10, type = "C-classification")
svm.pred <- predict(svm.model, testset[2:5])
#plot(svm.model,trainset[2:7])
tab <- table(pred = svm.pred, true = testset[,7])
tab
classAgreement(tab)
drops <- c("row.names","Contig.Name")
tuned <- tune.svm(Category~., data = trainset[2:6], gamma = 10^(0:1), cost = 10^(3:4))
summary(tuned)

# full svm binning run
runset <- d[c("name","HPplus","HPminus","HP3","gc")]
rownames(runset) <- runset$name
svm.model <- svm(Category~., data = trainset[2:6], gamma = 100, cost = 10, type = "C-classification")
svm.pred <- predict(svm.model, runset[2:5])
svmpred <- as.data.frame(svm.pred)
d3 <- merge(d,svmpred, by.x="name",by.y="row.names",All.x=T)
d3$svm.pred <- as.numeric(as.character(d3$svm.pred))
e3 <- merge(ess,svmpred, by.x="name",by.y="row.names",All.x=T)
e3$svm.pred <- as.numeric(as.character(e3$svm.pred))

#spot check individual bins (optional code)
#spot_bin = 10
#d4 <- subset(d3, svm.pred == spot_bin)
#hist(d4$gc,breaks=30)
#d4[is.na(d4$tax.color),"tax.color"] = 1
#e4 <- subset(e3, svm.pred == spot_bin)
#d4_cov <- subset(filtered_df, Contig.Name %in% d4$name,select=c("Contig.Name","Coverage_1","Coverage_2","Coverage_3","Coverage_4"))
#rownames(d4_cov) <- d4_cov$Contig.Name

#completeness stats
# Find the duplicated single copy essential genes TIGR00436, PF01795 and PF00750 is not single copy

summary.out <- {}
for (i in 1:length(unique_cluster)) { running.hmmid <- e3[e3$svm.pred == unique_cluster[i],]
                                      running.out <- running.hmmid[which(duplicated(running.hmmid$hmm.id) | duplicated(running.hmmid$hmm.id, fromLast=TRUE)),] 
                                      summary.out <- rbind(summary.out, c(unique_cluster[i], 
                                                                          length(d3$length[d3$svm.pred == unique_cluster[i]]),
                                                                          sum(d3$length[d3$svm.pred == unique_cluster[i]]),
                                                                          length(unique(running.hmmid$hmm.id)), 
                                                                          nrow(running.out),
                                                                          length(unique(running.out$hmm.id)),
                                                                          mean(d3$length[d3$svm.pred == unique_cluster[i]]),
                                                                          max(d3$length[d3$svm.pred == unique_cluster[i]]),
                                                                          mean(colMeans(subset(Coassembly_coverage, Contig.Name %in% d3$name[d3$svm.pred == unique_cluster[i]])[,3:6])),
                                                                          length(unique(d3$phylum[d3$svm.pred == unique_cluster[i]]))-1))
}
summary.out <- as.data.frame(summary.out)
colnames(summary.out) = c("bin","total_contigs", "total_bp", "completeness", "total_duplicates", "unique_duplicates", "avg_len","max_len","mean_cov","est_phyla")
summary.out

#show bins
dev.off()
palette(adjustcolor(gbr(dbspan+1), alpha.f = 1))
plot(x = d3$HPminus,
     y = d3$HPplus,
     xlim = c(99.3, 101.2),
     ylim = c(98.4, 101.8),
     cex = sqrt(d$length)/100,
     pch=20,
     col=rgb(0,0,0,0.05),
     xlab = "PC1",
     ylab = "PC2"     
)

points(x = d3$HPminus[d3$svm.pred > 0],
     y = d3$HPplus[d3$svm.pred > 0],
     cex = .8,
     pch=20,
     col=d3$svm.pred[d3$svm.pred > 0]
)

# Add a label in the center of each extracted bin

for (i in 1:length(unique_cluster)) {
  text(mean(d3$HPminus[d3$svm.pred ==unique_cluster[i]]),
       mean(d3$HPplus[d3$svm.pred == unique_cluster[i]]),
       labels = as.character(unique_cluster[i]),
       cex=1,
       font=1
  )
}

#label a specific point (optional)
#scaffold_name = "scaffold-9416"
#points(x = d3$HPminus[d3$name == scaffold_name],
#       y = d3$HPplus[d3$name == scaffold_name],
#       cex = 3,
#       pch=20,
#       col=rgb(0,0,0,1)
#)

#manually create svm bin
# Interactively chose 6 points on the plot that include the subset of scaffolds targeted
#z.def <- ahull(locator(6, type="p", pch=20), alpha=100000)                 
# Mark the defined subset on the plot
#plot(z.def,add=T, col="black")
#sourcebin = 7420
# Extract the scaffolds (gcent.out) and essential genes (ecent.out) in the defined subset
#gcent.out <- {}
#for (i in 1:nrow(d3)) { if (inahull(z.def, c(d3$HPminus[i],d3$HPplus[i]))) gcent.out <- rbind(gcent.out,d3[i,])}
#ecent.out<-{}
#for (i in 1:nrow(e3)) { if (inahull(z.def, c(e3$HPminus[i],e3$HPplus[i]))) ecent.out <- rbind(ecent.out,e3[i,])}
#hist(gcent.out$gc,breaks=40,xlim=c(30,60))
#gcent.out <- gcent.out[gcent.out$gc < 52  & gcent.out$gc > 35,]
#ecent.out <- ecent.out[ecent.out$gc < 52  & ecent.out$gc > 35,]
#d3$svm.pred[d3$name %in% gcent.out$name] = sourcebin
#e3$svm.pred[e3$name %in% ecent.out$name] = sourcebin
#unique_cluster <- c(unique_cluster, sourcebin)

#spot check bins and split bins
#hist(d3$gc[d3$svm.pred == 1])

#cross reference phylum with svm bins
table(d3$phylum,d3$svm.pred)

#view in 3d

plot3d(d3$HPminus[d3$svm.pred == 0], 
       d3$HPplus[d3$svm.pred == 0], 
       d3$HP3[d3$svm.pred == 0],  
       col=rgb(0,0,0,0.1),         
       size=.3, 
       type='p',
       xlab="PC1",
       ylab="PC2",
       zlab="PC3")

points3d(d3$HPminus[d3$svm.pred > 0], 
         d3$HPplus[d3$svm.pred > 0], 
         d3$HP3[d3$svm.pred > 0],
         col=d3$svm.pred[d3$svm.pred > 0],
         size=3)

bintext <-{}

# Add a label in the center of each extracted bin
for (i in 1:length(unique_cluster)) {
  bintext <- rbind(bintext, c(as.character(unique_cluster[i]),
                              mean(d3$HPminus[d3$svm.pred ==unique_cluster[i]]), 
                              mean(d3$HPplus[d3$svm.pred == unique_cluster[i]]), 
                              mean(d3$HP3[d3$svm.pred == unique_cluster[i]])))
}
text3d(x=bintext[,2], y=bintext[,3], z=bintext[,4],text = bintext[,1])

#record a movie, if needed
M <- par3d("userMatrix") 
play3d( spin3d(c(0,1,0), rpm=2 ), duration=30) 
movie3d( spin3d(c(1,0,0), rpm=3 ), duration=20, dir="./raw\ movie", convert=FALSE, top= TRUE, clean=FALSE, type="gif", startTime=0)

### Write SVM List to file #####################################################

#manual list of bins to export
#writelist <- c(2,4,8,9,11,21,29,30,32,34,610,15501)
#for (i in 1:length(writelist)) {
#  write.table(d3$name[d3$svm.pred == writelist[i]],file=paste(writelist[i],"_svm_bin.txt",sep=""),quote=F,row.names=F,col.names=F)
#}

#export all bins
write.table(d3[,c("name","svm.pred")],file="10k_EPS_0_35_all_svm_bin.txt",quote=F,row.names=F,col.names=F)

### Load in saved bin names for viewing purposes #####################################################
filelist <- list.files("./svm_binning/",pattern="^[0-9]*_svm_bin.txt")

bin_df <- {}
mean_df <- {}
for (i in 1:(length(filelist))) { 
  
  read_df <- read.csv(paste("./svm_binning/", filelist[i], sep=""), header = F, sep = "\t", quote = "")
  binnum <- strsplit(filelist[i],"_")[[1]][1]
  colnames(read_df)[1] = "contig"
  read_df["gc"] <- gc[match(read_df$contig, gc$contig),"gc"]
  read_df["length"] <- Coassembly_coverage[match(read_df$contig, Coassembly_coverage$Contig.Name),"Reference.length"]
  mean_df <- rbind(mean_df, c(binnum, mean(Coassembly_coverage[match(read_df$contig, Coassembly_coverage$Contig.Name),"Reference.length"])))
  mean(Coassembly_coverage[match(read_df$contig, Coassembly_coverage$Contig.Name),"Reference.length"])
  
  read_df["bin"] = binnum
  bin_df <- rbind(bin_df, read_df)
}

hist(bin_df[bin_df$bin == 4,"gc"],breaks=15)

file_svm_pred <- read.csv("all_svm_bin.txt",header=F, sep=" ")
colnames(file_svm_pred) <- c("name","file_svmpred")
d3_file <- merge(d3, file_svm_pred, by.x ="name", by.y="name")
e3_file <- merge(e3, file_svm_pred, by.x ="name", by.y="name")

file_summary.out <- {}
file_unique_cluster <- unique(d3_file$file_svmpred)
for (i in 1:length(file_unique_cluster)) { running.hmmid <- e3[e3_file$file_svmpred == file_unique_cluster[i],]
                                      running.out <- running.hmmid[which(duplicated(running.hmmid$hmm.id) | duplicated(running.hmmid$hmm.id, fromLast=TRUE)),] 
                                      file_summary.out <- rbind(file_summary.out, c(file_unique_cluster[i], 
                                                                          length(d3_file$length[d3_file$file_svmpred == file_unique_cluster[i]]),
                                                                          sum(d3_file$length[d3_file$file_svmpred == file_unique_cluster[i]]),
                                                                          length(unique(running.hmmid$hmm.id)), 
                                                                          nrow(running.out),
                                                                          length(unique(running.out$hmm.id)),
                                                                          mean(d3_file$length[d3_file$file_svmpred == file_unique_cluster[i]]),
                                                                          max(d3_file$length[d3_file$file_svmpred == file_unique_cluster[i]]),
                                                                          mean(colMeans(subset(Coassembly_coverage, Contig.Name %in% d3_file$name[d3_file$file_svmpred == file_unique_cluster[i]])[,2:5]))))
}
file_summary.out <- as.data.frame(file_summary.out)
colnames(file_summary.out) = c("bin","total_contigs", "total_bp", "completeness", "total_duplicates", "unique_duplicates", "avg_len","max_len","mean_cov")
file_summary.out

