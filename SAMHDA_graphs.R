####### graph the clusters 
# 1) line chart
# x = year; y = p(yes); line/colour = cluster
# 2) bar chart
# x = cluster; y = total freq across year
# cluster = sum(reasons in cluster); cluster_bin = assignment
# red = avoidance; blue = embarrassment; green = cost

# need to add cluster, cluster_bin to df
# 1 = endorsed; 0 = not endorsed
cluster.df <- reasons_df

# lists of reasons by cluster
avoid.clust <- notmt_clust$cNames$cluster1
embar.clust <- notmt_clust$cNames$cluster2
cost.clust  <- notmt_clust$cNames$cluster3

# cluster -- sum(reasons in cluster)
cluster.df$clust1 <- rowSums(cluster.df[,names(cluster.df) %in% avoid.clust]) # max = 7
cluster.df$clust2 <- rowSums(cluster.df[,names(cluster.df) %in% embar.clust]) # max = 4
cluster.df$clust3 <- rowSums(cluster.df[,names(cluster.df) %in% cost.clust]) # max = 4

# cluster_bin
cluster.df$bin1 <- ifelse(cluster.df$clust1==0,0,1) # endorsed = 2619
cluster.df$bin2 <- ifelse(cluster.df$clust2==0,0,1) # endorsed = 1775
cluster.df$bin3 <- ifelse(cluster.df$clust3==0,0,1) # endorsed = 3730

# reshape cluster df
cluster.long1 <- melt(cluster.df[,c(1,19:21)],id.vars="year")
cluster.long2 <- melt(cluster.df[,c(1,22:24)],id.vars="year")



############## graph 1
# table of endorsed/not by year by bin
freq.df        <- table(cluster.long2$year,cluster.long2$value,cluster.long2$variable) %>% as.data.frame
names(freq.df) <- c("year","endorsed","bin","count")

# just endorsed + calculate percent
freq.df.end <- freq.df[freq.df$endorsed=="1",]
totals      <- table(cluster.long2$year,cluster.long2$variable) %>% as.data.frame()
freq.df.end <- cbind(freq.df.end,totals$Freq)

names(freq.df.end)[5]   <- "total"
freq.df.end$per         <- freq.df.end$count / freq.df.end$total
levels(freq.df.end$bin) <- c("avoidance","embarrassment","cost")


# graph 1
# percent of patients per bin by year
# (bin = patients who endorsed any items within a cluster)
cluster.year <- ggplot(freq.df.end, aes(x=year,y=per,group=bin,colour=bin)) + 
  geom_line() +
  geom_point() +
  scale_colour_manual(values=c("red","blue","green")) +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=13, color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=13, color="black"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  xlab("Year") +
  ylab("Percent Endorsed")
cluster.year



############## graph 2
# table of endorsed/not by bin (across year)
freq.df2        <- table(cluster.long2$value,cluster.long2$variable) %>% as.data.frame
names(freq.df2) <- c("endorsed","bin","count")

# just endorsed + calculate percent
freq.df.end2 <- freq.df2[freq.df2$endorsed=="1",]
totals2      <- table(cluster.long2$variable) %>% as.data.frame()
freq.df.end2 <- cbind(freq.df.end2,totals2$Freq)

names(freq.df.end2)[4]   <- "total"
freq.df.end2$per         <- freq.df.end2$count / freq.df.end2$total
levels(freq.df.end2$bin) <- c("avoidance","embarrassment","cost")

# graph 2
# percent of patients in each bin across years
cluster.coll <- ggplot(freq.df.end2, aes(x=bin,y=per,fill=bin)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("red","blue","green")) +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=13, color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  xlab("") +
  ylab("Percent Endorsed (2008-2014)")
cluster.coll
