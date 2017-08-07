### for the paper
# importance df
fig1           <- caret::varImp(out.model.enet, scale=TRUE)[[1]] %>% as.data.frame
fig1           <- cbind(rownames(fig1), fig1) %>% arrange(-Overall)
names(fig1)    <- c("Predictor", "Importance")
fig1$Predictor <- factor(fig1$Predictor,
                         levels = fig1$Predictor[order(fig1$Importance)])

# code by var type
fig1$var_type <- ifelse(fig1$Predictor %in% mdd_vars,"clinical",
                        ifelse(fig1$Predictor %in% ses_vars,"demographic",
                               ifelse(fig1$Predictor %in% med_vars,"med_hx",
                                      ifelse(fig1$Predictor %in% k6_vars,"clinical","FIX"))))
fig1[c(4,7,8,11,13,20),3] <- "demographic"     # i'm lazy and only fixing the top 20

# rename vars for readability
longer <- lapply(fig1$Predictor, defineMe) %>%
  unlist %>%
  stringr::str_replace(., pattern = "NO MH TMT ", replacement = "")

pred.dict   <- data.frame(fig1$Predictor,longer)
fig1$labels <- pred.dict[match(fig1$Predictor,pred.dict$fig1.Predictor),2] %>% as.character

fig1[4,4]  <- "MULTIRACIAL"           # only fixing top 20
fig1[7,4]  <- "NATIVE AMERICAN"
fig1[8,4]  <- "FEMALE"
fig1[11,4] <- "WIDOWED"
fig1[13,4] <- "PACIFIC ISLANDER"
fig1[20,4] <- "BLACK"

# importance has to be ordered by x-axis
fig1.20          <- fig1[1:20,]
fig1.20$labels   <- factor(fig1.20$labels,
                           levels = fig1.20$labels[order(fig1.20$Importance)])
fig1.20$var_type <- as.factor(fig1.20$var_type)
levels(fig1.20$var_type) <- c("symptoms","demographic","medical history")

# figure 
ggplot(fig1.20, aes(x=labels,y=Importance)) + 
  geom_bar(aes(fill=var_type),stat="identity",width=.85) +
  geom_text(stat="identity",y=18, hjust=0, aes(label=labels), colour="black", size=4) +
  scale_fill_manual(values=c("#3399FF","#33CC00","#FF3300")) +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=13,color="black"),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=13),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust=-.06,vjust=1.8,size=18,face="bold")) +
  ylab("\nVariable Importance") +
  xlab("") +
  coord_flip(ylim = c(20,100))