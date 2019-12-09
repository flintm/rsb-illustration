# Bar chart for alternative configuration ranking
# RSB Overview Illustration
# Madeleine Flint, 2019-10-29

# plotting setup
library(ggplot2)
PPcols = c(base = 'gray50', # black
           F_y = '#f5c842', # mustard, F_y
           delta_p  = '#12a2db', # blue, delta_p
           delta_pc = '#42d483', # delta_pc, green
           combined = '#25a88e') # teal, combined

# read results and organize dataframe
df.rank = read.table('rank.txt', header = TRUE, sep = "\t", stringsAsFactors = FALSE)
df.rank$Struct = factor(grepl('conc',df.rank$alt_rank), 
                        levels = c(TRUE, FALSE), 
                        labels = c("concrete~MRF", "steel~MRF"))
df.rank$Rank = factor(1:nrow(df.rank), levels = nrow(df.rank):1, labels = as.character(nrow(df.rank):1))

PPs = c("base", "F_y", "delta_p", "delta_pc", "combined")
PPlabels = c("--", "F_y", "delta_p", "delta_pc", "F_y, delta_p, delta_pc")
df.rank$PP[grepl("Base",df.rank$alt_rank)] = "base"
df.rank$PP[grepl("combined",df.rank$alt_rank)] = "combined"
df.rank$PP[grepl("theta_pc",df.rank$alt_rank)] = "delta_pc"
df.rank$PP[grepl("theta_p\\>",df.rank$alt_rank)] = "delta_p"
df.rank$PP[grepl("F_y",df.rank$alt_rank)] = "F_y"
df.rank$PP <- factor(df.rank$PP, levels = PPs)#, labels = PPlabels)

# performance parameter values
df.rank$thetaL[df.rank$PP=="base"] = 0
df.rank$thetaL[grepl("1,",df.rank$alt_rank)] = 1
df.rank$thetaL[grepl("0.5,",df.rank$alt_rank)] = 0.5
df.rank$thetaL = factor(df.rank$thetaL, levels = c(0,0.5,1), labels = c("0","0.5","1"))

df.rank$thetaE[grepl("(: 0)$",df.rank$alt_rank)] = 0
df.rank$thetaE[grepl("(: 1)$",df.rank$alt_rank)] = 1
df.rank$thetaE[grepl("(: 0.5)$",df.rank$alt_rank)] = 0.5
df.rank$thetaE = factor(df.rank$thetaE, levels = c(0,0.5,1), labels = c("0","0.5","1"))

# names
df.rank$names = paste0("italic(l)[italic(k)]==",df.rank$Struct,"~bold(theta)[list(k,kappa)]==~phantom(1)~list(0, 0, ",df.rank$thetaL,", ",df.rank$thetaE,")")
write.table(df.rank, file = "rank_org.txt", sep = "\t")

# plot------
#as.expression(parse(text=levels(factor(limits, levels = limits, labels = limits))))
p <- ggplot(subset(df.rank,as.integer(as.character(Rank))<=10), aes(x=Rank, y = pF_rank, fill = PP)) +
  geom_col(width = 0.5) + coord_flip() + scale_x_discrete(breaks = as.character(10:1), 
                        labels = as.expression(parse(text=df.rank$names[10:1])), 
                        name = 'Alternative Configuration', expand = c(0,0)) +
  scale_fill_manual(values=PPcols,breaks=PPs,
                    labels=c(base = PPlabels[1],
                             F_y = expression(italic(F)[italic(y)]),
                             delta_p = expression(italic(delta)[italic(p)]),
                             delta_pc = expression(italic(delta)[italic(pc)]),
                             combined = expression(list(italic(F)[italic(y)], italic(delta)[italic(p)], italic(delta)[italic(pc)]))), 
                    name = "Lateral PP", drop=FALSE) +
  scale_y_continuous(limits = c(0,0.5),expand = c(0,0.01)
                     ,name='Probability of Failing Preference System') +
  guides(fill=guide_legend(direction = "horizontal"))+
  theme(panel.background   = element_rect(fill = "white"),
        legend.key         = element_rect(fill = "white") ,
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        #panel.spacing      = unit(c(0,0,0,0),"cm"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        axis.title.x = element_text(color = "black", size = 10),
        axis.title.y = element_text(color = "black", size = 10),
        legend.text  = element_text(color = "black", size = 10,hjust=0),
        legend.title = element_text(color = "black", size = 10),
        legend.key.height=unit(0.7,"line"),
        legend.position = "bottom")#,
        # plot.margin     = unit(c(0,0,0,0), "cm"))
p
ggsave(p,filename = "bar.pdf", width = 6.5, height = 3,useDingbats = F)
ggsave(p,filename = "bar.eps", width = 6.5, height = 3)
