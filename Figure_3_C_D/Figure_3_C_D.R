#load libraries 

library(tidyverse)

##############################################################################################################
#setup the directories
Basedirectory=file.path("~","Opto-CLIP","Figure_3_C_D")
list.files(Basedirectory)

Datadirectory=file.path(Basedirectory,"Data")
list.files(Datadirectory)

Outdirectory=file.path(Basedirectory,"Graphs")
dir.create(Outdirectory,showWarnings = TRUE)

##############################################################################################################
#Bring in current clamp data
cc200pA= read_delim(file.path(Datadirectory,"CurrentClamp_CA1_mV.csv"), col_names = TRUE)

cc200pA_correct <-  cc200pA %>% 
  mutate_at(vars(L:V), 
            .funs = funs(. * 1000000000000))

df.mV <- cc200pA_correct %>% 
  select(Time,A:K) %>%
  gather(key = "trace", value = "mV", -Time)

df.pA <- cc200pA_correct  %>% 
  select(Time,L:V) %>%
  gather(key = "trace", value = "pA", -Time)

#Figure 3C- label the rheostat current green
mV_plot=ggplot(df.mV, aes(x = Time, y = mV)) + 
  xlim(5,6) +
  geom_line(aes(color = trace), show.legend=FALSE, linewidth=1) +
  theme_void() +
  scale_color_manual(values = c(gray.colors(4, start = 0.1, end = 0.4), 'green', gray.colors(6, start = 0.4, end = 0.7))) + 
  annotate("segment", x = 5.9, y= -80, xend = 6.0, yend= -80, color="black", linewidth=1) + 
  annotate("segment", x = 6.0, y= -80, xend = 6.0, yend= -70, color="black", linewidth=1)
mV_plot
ggsave(filename = file.path(Outdirectory,"Figure3C_scale_10mV_100ms.pdf"), plot=mV_plot, device='pdf') 

#Figure 3D- label the rheostat current green
pA_plot=ggplot(df.pA, aes(x = Time, y = pA)) + 
  geom_line(aes(color = trace), show.legend=FALSE, linewidth=1) +
  theme_void() +
  scale_color_manual(values = c(gray.colors(4, start = 0.1, end = 0.4), 'green', gray.colors(6, start = 0.4, end = 0.7))) + 
  xlim(5,6) +
  annotate("segment", x = 5.9, y= -50, xend = 6.0, yend= -50, color="black", linewidth=1) + 
  annotate("segment", x = 6.0, y= -50, xend = 6.0, yend= -25, color="black", linewidth=1)
pA_plot
ggsave(filename = file.path(Outdirectory,"Figure3D_scale_25pA_100ms.pdf"), plot=pA_plot, device='pdf') 
