library(dplyr)
library(ggplot2)
library(ggrepel)
library(forcats)
library(scales)

mydf <- structure(list(Group = structure(c(1L, 2L, 3L, 4L, 5L), .Label = c("snv", 
"del", "ins", "mnp", "clumped"), class = "factor"), value = c( 241844L, 246950L, 
234064L,30336L, 25262L)), .Names = c("Group", "value"), class = "data.frame", row.names = c("1", 
"2", "3", "4", "5"))

mydf %>%
arrange(desc(value)) %>%
mutate(prop = percent(value / sum(value))) -> mydf 

pie <- ggplot(mydf, aes(x = "", y = value, fill = fct_inorder(Group))) +
       geom_bar(width = 1, stat = "identity") +
       coord_polar("y", start = 0) +
       geom_label_repel(aes(label = prop), size=5, show.legend = F, nudge_x = 1) +
       guides(fill = guide_legend(title = "Group")) + theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank())
       
bxp <- pie + labs(title = "Distribution of the variants from C57BL6J and DBA2J",
              x = " ", y = " ") + theme(panel.background = element_rect(fill = 'white', colour = 'white'))
