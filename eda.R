## exploratory analysis
library(dplyr)
library(ggplot2)
library(cowplot)

#### read sequences and their respective clusters from each country
#### users need to prepare these files by downloading the data from GISAID and getting 
#### them into the required format as shown in the supplementary document
US_cluster <- read.csv("US_cluster.csv")
CA_cluster <- read.csv("CA_cluster.csv")
AU_cluster <- read.csv("AU_cluster.csv")
UK_cluster <- read.csv("UK_cluster.csv")
NL_cluster <- read.csv("NL_cluster.csv")
FR_cluster <- read.csv("FR_cluster.csv")
SP_cluster <- read.csv("SP_cluster.csv")
CN_cluster <- read.csv("CN_cluster.csv")
IN_cluster <- read.csv("IN_cluster.csv")

#### group sequences by dates
cluster_data <- function(df){
  
  df[["Cluster1"]] <- ifelse(df[["Cluster"]]== 1,1,0)
  df[["Cluster2"]] <- ifelse(df[["Cluster"]]== 2,1,0)
  df[["Cluster3"]] <- ifelse(df[["Cluster"]]== 3,1,0)
  df[["Cluster4"]] <- ifelse(df[["Cluster"]]== 4,1,0)
  df[["Cluster5"]] <- ifelse(df[["Cluster"]]== 5,1,0)
  df <- df %>% group_by(.data[["Date"]]) %>%
    summarise(Cluster1=sum(.data[["Cluster1"]]), 
              Cluster2=sum(.data[["Cluster2"]]), 
              Cluster3=sum(.data[["Cluster3"]]),
              Cluster4=sum(.data[["Cluster4"]]),
              Cluster5=sum(.data[["Cluster5"]]))
  return(df)
}

US_data <- cluster_data(US_cluster)
CA_data <- cluster_data(CA_cluster)
UK_data <- cluster_data(UK_cluster)
NL_data <- cluster_data(NL_cluster)
FR_data <- cluster_data(FR_cluster)
SP_data <- cluster_data(SP_cluster)
CN_data <- cluster_data(CN_cluster)
IN_data <- cluster_data(IN_cluster)
AU_data <- cluster_data(AU_cluster)

#### summarize the number of sequences by dates for each country
ts.data <- data.frame(Date = seq(as.Date("2019-12-24"), as.Date("2020-10-08"), by="days"))

cont.fun <- function(df){
  df$Date <- as.Date(df$Date, format = "%Y-%m-%d") 
  df.tmp <-left_join(ts.data, df, by = "Date")
  df.tmp <- cbind(Index = 1:length(df.tmp$Date), df.tmp)
  df.tmp <- df.tmp[,1:7]
  df.tmp[is.na(df.tmp)] <- 0
  return(df.tmp)
}

US_daily <- cont.fun(US_data)
CA_daily <- cont.fun(CA_data)
UK_daily <- cont.fun(UK_data)
NL_daily <- cont.fun(NL_data)
FR_daily <- cont.fun(FR_data)
SP_daily <- cont.fun(SP_data)
CN_daily <- cont.fun(CN_data)
IN_daily <- cont.fun(IN_data)
AU_daily <- cont.fun(AU_data)

#### plot EDA for US and UK as in Figure 2
data_prop <- function(df){
  
  df[["Sum"]] <- rowSums(df[,c("Cluster1", "Cluster2", "Cluster3",
                              "Cluster4", "Cluster5")])
  
  df[["Prop1"]] <- df[["Cluster1"]]/df[["Sum"]]
  df[["Prop2"]] <- df[["Cluster2"]]/df[["Sum"]]
  df[["Prop3"]] <- df[["Cluster3"]]/df[["Sum"]]
  df[["Prop4"]] <- df[["Cluster4"]]/df[["Sum"]]
  df[["Prop5"]] <- df[["Cluster5"]]/df[["Sum"]]
  return(df)
}

UK_prop <- data_prop(UK_daily)
US_prop <- data_prop(US_daily)

clist <- c("US", "UK")
clist <- factor(clist, levels=c("US","UK"))
count <- c(US_prop$Cluster1, US_prop$Cluster2,
           US_prop$Cluster3,US_prop$Cluster4, US_prop$Cluster5,
           UK_prop$Cluster1, UK_prop$Cluster2,
           UK_prop$Cluster3,UK_prop$Cluster4, UK_prop$Cluster5)
prop <- c(US_prop$Prop1, US_prop$Prop2,
         US_prop$Prop3,US_prop$Prop4, US_prop$Prop5,
         UK_prop$Prop1, UK_prop$Prop2,
         UK_prop$Prop3,UK_prop$Prop4, UK_prop$Prop5)

df <- data.frame(Date = rep(seq(as.Date("2019-12-24"), 
                                as.Date("2020-10-08"), by="days"), 5),
                 count = count,
                 prop = prop,
                 category = rep(rep(c("Cluster I","Cluster II", "Cluster III", 
                                  "Cluster IV", "Cluster V"), each = 290),2),
                 country = rep(clist, each = 290*5))



p1 <- ggplot(df, aes(x=Date, y = count, color = category)) + geom_point() + 
  facet_wrap(~country)+
  geom_smooth(method = "loess", se=FALSE, span =0.3) +
  xlab("Dates") + ylab("Count") + ylim(c(0,400))+
  theme(legend.position="bottom")+
  theme(legend.title = element_blank())

p2 <- ggplot(df, aes(x=Date, y = prop, color = category)) + geom_point() + 
  facet_wrap(~country)+
  geom_smooth(method = loess, se=FALSE) +
  xlab("Dates") + ylab("Proportion") + ylim(c(0,1)) + 
  theme(legend.position="bottom")+
  theme(legend.title = element_blank())


legend1 <- get_legend(p1)
prow <- plot_grid(p1+ theme(legend.position = "none"), 
                  p2+ theme(legend.position = "none"),
                  align = 'vh',
                  hjust = -1, nrow = 2,axis="1")
p <- plot_grid(prow, legend1, rel_heights = c(1, 0.1), nrow = 2, scale = c(1,0.1))

pdf(file="plot_cluster_USUK.pdf", width = 6, height = 5.5)  
p
dev.off()


#### save daily cluster proportions in arrays used for sampling

array_allcont <- array(dim=c(290,9,5))
for (n in 1:290){
  array_allcont[n,1,] <- as.numeric(US_daily[n,c(3:7)])
  array_allcont[n,2,] <- as.numeric(CA_daily[n,c(3:7)])
  array_allcont[n,3,] <- as.numeric(UK_daily[n,c(3:7)])
  array_allcont[n,4,] <- as.numeric(NL_daily[n,c(3:7)])
  array_allcont[n,5,] <- as.numeric(FR_daily[n,c(3:7)])
  array_allcont[n,6,] <- as.numeric(SP_daily[n,c(3:7)])
  array_allcont[n,7,] <- as.numeric(CN_daily[n,c(3:7)])
  array_allcont[n,8,] <- as.numeric(IN_daily[n,c(3:7)])
  array_allcont[n,9,] <- as.numeric(AU_daily[n,c(3:7)])
}

save(array_allcont, file = "array_allcont.rda")



