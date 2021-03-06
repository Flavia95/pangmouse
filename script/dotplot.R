library(ggplot2)
library(tidyverse)
library(argparse)
args = commandArgs(trailingOnly=TRUE)
input=as.numeric(args[1])


myd <- read.table(args[1], header=T, sep="\t", row.names=NULL)

x <- myd[,1, drop=FALSE]
x$homo  <- "homo"
y <- myd[,2, drop=FALSE]
y$mouse  <- "mouse"
df = data.frame (x, y)
colnames(df) = c("value_homo", "id_homo", "value_mouse", "id_mouse")
myd = df %>% gather(species,values,starts_with("value_"))
p = ggplot(myd, aes(values, values)) + geom_point(aes(colour = as.factor(species)))
r =  p  +  theme(panel.background = element_rect(fill = 'white', colour = 'white')) + labs(color='species') + scale_color_manual(labels = c("homosapiens", "mouse"), values = c("blue", "red"))
#q = ggplot(myd, aes(species, values)) + geom_point(aes(colour = as.factor(species)))
ggsave("./aln.png", plot= r, device="png", width = 20, height = 15, units = "cm", dpi = 300)
