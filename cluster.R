## clustering of sequences
library(dplyr)
library(ggplot2)
library(cowplot)

#### read unique S-protein sequences
#### users need to prepare the file by downloading the data from GISAID and getting 
#### them into the required format as shown in the supplementary document
protein <- read.csv("sprotein_unique.csv")
names(sprotein_unique)<- "Sequence"
sprotein_unique <- data.frame(Number = 1:length(sprotein_unique$Sequence), 
                              Sequence = sprotein_unique$Sequence)

#### create distance matrix that calculates pairwise mismatched letters
string.diff <- function (a, b){
  diff.a <-unlist(strsplit(a,split=""))
  diff.b <-unlist(strsplit(b,split=""))
  result <- sum(diff.a != diff.b)
  return (result)
}
n <- length(sprotein_unique$Sequence)
dist.matrix <- sapply(1:n, 
                      function(k) sapply(1:n, function(i) 
                        string.diff(sprotein_unique$Sequence[k], sprotein_unique$Sequence[i])))
dist <- as.dist(dist.matrix)

#### clustering based on ward D method with size k=5
hc <- hclust(dist, method = "ward.D")
plot(hc, cex = 0.03)
sub_grp <- cutree(hc, k = 5)
sub_grp[sub_grp == 2] <- "temp" # set cluster 2 as reference
sub_grp[sub_grp == 4] <- 2
sub_grp[sub_grp == "temp"] <- 4
table(sub_grp)

#### summarize mutation info in each cluster
##### split sequences into 1273 letters
sprotein_unique_unlist <- sprotein_unique
index <- 3:1275
sprotein_unique_unlist[index] <- list(character(0))
names(sprotein_unique_unlist)[index] <- 1:1273
for (i in 1:n){
  sprotein_unique_unlist[i, index] <-  unlist(strsplit(sprotein_unique$Sequence[i],split=""))
}

##### read and split reference sequence
#### users need to prepare the file by downloading the data from GISAID and getting 
#### them into the required format as shown in the supplementary document
ref <- read.csv("sprotein_reference.csv")
names(ref)<- "Sequence"
ref <- data.frame(Number = 1, Sequence = ref$Sequence)
ref[index] <- list(character(0))
names(ref)[index] <- 1:1273
ref[1,index] <- unlist(strsplit(ref[1,2],split=""))
refnum <- sprotein_unique$Number[sprotein_unique[,2]== ref[1,2]]

##### create dataframe that records mutations in each sequence
diff <- data.frame()
diff[1:50] <- list(character(0))
for (k in 1:n){
  if (k == refnum){
    diff[k, ] <- NA
  } else {
    func <- sprotein_unique_unlist[k,index] != sprotein_unique_unlist[refnum,index]
    diff.num <- sum(func)
    diff[k,1:diff.num] <- sapply(1:diff.num, function(i)
      paste(ref[1,which(func)+2][i], which(func)[i], 
            sprotein_unique_unlist[k,which(func)+2][i], sep = ""))
  }
}

##### mutations for each sequence in each cluster
diff.g5 <- diff[sub_grp == 5, ]
diff.g4 <- diff[sub_grp == 4, ]
diff.g3 <- diff[sub_grp == 3, ]
diff.g2 <- diff[sub_grp == 2, ]
diff.g1 <- diff[sub_grp == 1, ]

##### sort mutations by frequencies in each cluster

tbl1 <- table(stack(diff.g1)[,1]) %>%
  as.data.frame() %>% 
  arrange(desc(Freq))
tbl1$Percent <- tbl1$Freq/sum(sub_grp == 1)

tbl2 <- table(stack(diff.g2)[,1]) %>%
  as.data.frame() %>% 
  arrange(desc(Freq))
tbl2$Percent <- tbl2$Freq/sum(sub_grp == 2)

tbl3 <- table(stack(diff.g3)[,1]) %>%
  as.data.frame() %>% 
  arrange(desc(Freq))
tbl3$Percent <- tbl3$Freq/sum(sub_grp == 3)

tbl4 <- table(stack(diff.g4)[,1]) %>%
  as.data.frame() %>% 
  arrange(desc(Freq))
tbl4$Percent <- tbl4$Freq/sum(sub_grp == 4)

tbl5 <- table(stack(diff.g5)[,1]) %>%
  as.data.frame() %>% 
  arrange(desc(Freq))
tbl5$Percent <- tbl5$Freq/sum(sub_grp == 5)

list.df <- list(tbl1, tbl2, tbl3, tbl4, tbl5)
max.rows <- max(unlist(lapply(list.df, nrow), use.names = F))
list.df <- lapply(list.df, function(x) {
  na.count <- max.rows - nrow(x)
  if (na.count > 0L) {
    na.dm <- matrix(NA, na.count, ncol(x))
    colnames(na.dm) <- colnames(x)
    rbind(x, na.dm)
  } else {
    x
  }
})

tbl <- do.call(cbind, list.df)
write.csv(tbl, file="top_mutations.csv")

#### Top unique sequences in each cluster, ranked by frequency
#### users need to prepare the file by downloading the data from GISAID and getting 
#### them into the required format as shown in the supplementary document
protein <- read.csv("cluster_top_seq.csv")
string.diff2 <- function (a, b){
  diff.a <-unlist(strsplit(a,split=""))
  diff.b <-unlist(strsplit(b,split=""))
  result1 <- sum(diff.a != diff.b)
  result2 <- which(diff.a != diff.b)
  return (paste(diff.b[result2], result2, diff.a[result2], sep=""))
}

protein$mutation <- sapply(1:15, function(i) 
  paste0(string.diff2(protein$Seq[i], protein$Seq[4]), collapse="", sep=" "))

write.csv(protein, file = "top_sequences.csv")


#### plot the piechart
#### users need to prepare the file by downloading the data from GISAID and getting 
#### them into the required format as shown in the supplementary document
df <- read.csv("piechart.csv", fileEncoding="UTF-8-BOM")

df_US <- data.frame(Cluster = c("Cluster I", "Cluster II", "Cluster III", "Cluster IV","Cluster V"),
                    Count = df$US)
df_CA <- data.frame(Cluster = c("Cluster I", "Cluster II", "Cluster III", "Cluster IV","Cluster V"),
                    Count = df$CA)
df_UK <- data.frame(Cluster = c("Cluster I", "Cluster II", "Cluster III", "Cluster IV","Cluster V"),
                    Count = df$UK)
df_NL <- data.frame(Cluster = c("Cluster I", "Cluster II", "Cluster III", "Cluster IV","Cluster V"),
                    Count = df$NL)
df_FR <- data.frame(Cluster = c("Cluster I", "Cluster II", "Cluster III", "Cluster IV","Cluster V"),
                    Count = df$FR)
df_SP <- data.frame(Cluster = c("Cluster I", "Cluster II", "Cluster III", "Cluster IV","Cluster V"),
                    Count = df$SP)
df_CN <- data.frame(Cluster = c("Cluster I", "Cluster II", "Cluster III", "Cluster IV","Cluster V"),
                    Count = df$CN)
df_IN <- data.frame(Cluster = c("Cluster I", "Cluster II", "Cluster III", "Cluster IV","Cluster V"),
                    Count = df$IN)
df_AU <- data.frame(Cluster = c("Cluster I", "Cluster II", "Cluster III", "Cluster IV","Cluster V"),
                    Count = df$AU)

p1 <- ggplot(df_US,aes(factor(1),y=Count, group=Cluster, fill = Cluster)) + 
  geom_bar(width = 1, stat = "identity", color = "white") +  
  coord_polar("y", start=0)  + ylab("US") + xlab("") + labs(fill="")+
  theme(legend.position="bottom",
        panel.background = element_blank(),
        axis.ticks = element_blank(), 
        panel.grid = element_blank(), 
        axis.text  = element_blank()) +
  theme(legend.title = element_blank())

p2 <- ggplot(df_CA,aes(factor(1),y=Count, group=Cluster, fill = Cluster)) + 
  geom_bar(width = 1, stat = "identity", color = "white") +  
  coord_polar("y", start=0)  + ylab("CA") + xlab("") + labs(fill="")+
  theme(legend.position="bottom",
        panel.background = element_blank(),
        axis.ticks = element_blank(), 
        panel.grid = element_blank(), 
        axis.text  = element_blank()) +
  theme(legend.title = element_blank())

p3 <- ggplot(df_UK,aes(factor(1),y=Count, group=Cluster, fill = Cluster)) + 
  geom_bar(width = 1, stat = "identity", color = "white") +  
  coord_polar("y", start=0)  + ylab("UK") + xlab("") + labs(fill="")+
  theme(legend.position="bottom",
        panel.background = element_blank(),
        axis.ticks = element_blank(), 
        panel.grid = element_blank(), 
        axis.text  = element_blank()) +
  theme(legend.title = element_blank())

p4 <- ggplot(df_NL,aes(factor(1),y=Count, group=Cluster, fill = Cluster)) + 
  geom_bar(width = 1, stat = "identity", color = "white") +  
  coord_polar("y", start=0)  + ylab("NL") + xlab("") + labs(fill="")+
  theme(legend.position="bottom",
        panel.background = element_blank(),
        axis.ticks = element_blank(), 
        panel.grid = element_blank(), 
        axis.text  = element_blank()) +
  theme(legend.title = element_blank())

p5 <- ggplot(df_FR,aes(factor(1),y=Count, group=Cluster, fill = Cluster)) + 
  geom_bar(width = 1, stat = "identity", color = "white") +  
  coord_polar("y", start=0)  + ylab("FR") + xlab("") + labs(fill="")+
  theme(legend.position="bottom",
        panel.background = element_blank(),
        axis.ticks = element_blank(), 
        panel.grid = element_blank(), 
        axis.text  = element_blank()) +
  theme(legend.title = element_blank())

p6 <- ggplot(df_SP,aes(factor(1),y=Count, group=Cluster, fill = Cluster)) + 
  geom_bar(width = 1, stat = "identity", color = "white") +  
  coord_polar("y", start=0)  + ylab("SP") + xlab("") + labs(fill="")+
  theme(legend.position="bottom",
        panel.background = element_blank(),
        axis.ticks = element_blank(), 
        panel.grid = element_blank(), 
        axis.text  = element_blank()) +
  theme(legend.title = element_blank())

p7 <- ggplot(df_CN,aes(factor(1),y=Count, group=Cluster, fill = Cluster)) + 
  geom_bar(width = 1, stat = "identity", color = "white") +  
  coord_polar("y", start=0)  + ylab("CN") + xlab("") + labs(fill="")+
  theme(legend.position="bottom",
        panel.background = element_blank(),
        axis.ticks = element_blank(), 
        panel.grid = element_blank(), 
        axis.text  = element_blank()) +
  theme(legend.title = element_blank())

p8 <- ggplot(df_IN,aes(factor(1),y=Count, group=Cluster, fill = Cluster)) + 
  geom_bar(width = 1, stat = "identity", color = "white") +  
  coord_polar("y", start=0)  + ylab("IN") + xlab("") + labs(fill="")+
  theme(legend.position="bottom",
        panel.background = element_blank(),
        axis.ticks = element_blank(), 
        panel.grid = element_blank(), 
        axis.text  = element_blank()) +
  theme(legend.title = element_blank())

p9 <- ggplot(df_AU,aes(factor(1),y=Count, group=Cluster, fill = Cluster)) + 
  geom_bar(width = 1, stat = "identity", color = "white") +  
  coord_polar("y", start=0)  + ylab("AU") + xlab("") + labs(fill="")+
  theme(legend.position="bottom",
        panel.background = element_blank(),
        axis.ticks = element_blank(), 
        panel.grid = element_blank(), 
        axis.text  = element_blank()) +
  theme(legend.title = element_blank())

legend1 <- get_legend(p1)
prow <- plot_grid(p1+ theme(legend.position = "none"), 
                  p2+ theme(legend.position = "none"),
                  p3+ theme(legend.position = "none"), 
                  p4+ theme(legend.position = "none"),
                  p5+ theme(legend.position = "none"), 
                  p6+ theme(legend.position = "none"),
                  p7+ theme(legend.position = "none"), 
                  p8+ theme(legend.position = "none"),
                  p9+ theme(legend.position = "none"),
                  align = 'none',labels = c("a", "b","c","d","e","f","g","h","i"),
                  label_size = 10,
                  hjust = -1, nrow = 2,axis="1")

prow <- plot_grid(p1+ theme(legend.position = "none"), 
                  p2+ theme(legend.position = "none"),
                  p3+ theme(legend.position = "none"), 
                  p4+ theme(legend.position = "none"),
                  p5+ theme(legend.position = "none"), 
                  align = 'v',
                  label_size = 10,
                  hjust = -1, nrow = 1,axis="1")

prow2 <- plot_grid(NULL,p6+ theme(legend.position = "none"),
                   p7+ theme(legend.position = "none"), 
                   p8+ theme(legend.position = "none"),
                   p9+ theme(legend.position = "none"), NULL,
                   align = 'vh', 
                   label_size = 10, scale = c(0.001,1,1,1,1,0.001),
                   hjust = -1, nrow = 1,axis="1")

p <- plot_grid(prow, prow2, legend1, rel_heights = c(1,1,0.1), nrow = 3, scale = c(1,1.2,0.1))
p

plot(legend1)
pdf(file="piechart2.pdf", width = 9, height = 4)  
p
dev.off()


