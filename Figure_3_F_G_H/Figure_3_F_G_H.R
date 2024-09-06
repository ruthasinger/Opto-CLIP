library(tidyverse)

##############################################################################################################
#setup the directories
Basedirectory=file.path("~","Opto-CLIP","Figure_3_F_G_H")
list.files(Basedirectory)

Datadirectory=file.path(Basedirectory,"Data")
list.files(Datadirectory)

Outdirectory=file.path(Basedirectory,"Graphs")
dir.create(Outdirectory,showWarnings = TRUE)

##############################################################################################################
#Bring in AAV-Control patch clamping data
LED_1Hz_Control <- read_delim(file.path(Datadirectory,"Patch_Control_1Hz.csv"), col_names = TRUE)
colnames(LED_1Hz_Control)=c("Time","LED_1Hz_Control")

LED_5Hz_Control <- read_delim(file.path(Datadirectory,"Patch_Control_5Hz.csv"), col_names = TRUE)
colnames(LED_5Hz_Control)=c("Time","LED_5Hz_Control")

LED_10Hz_Control <- read_delim(file.path(Datadirectory,"Patch_Control_10Hz.csv"), col_names = TRUE)
colnames(LED_10Hz_Control)=c("Time","LED_10Hz_Control")

#ChR2 files
LED_1Hz_ChR2 <- read_delim(file.path(Datadirectory,"Patch_ChR2_1Hz.csv"), col_names = TRUE)
colnames(LED_1Hz_ChR2)=c("Time","LED_1Hz_ChR2")
LED_1Hz_ChR2=LED_1Hz_ChR2[1:210000,]

LED_5Hz_ChR2 <- read_delim(file.path(Datadirectory,"Patch_ChR2_5Hz.csv"), col_names = TRUE)
colnames(LED_5Hz_ChR2)=c("Time","LED_5Hz_ChR2")

LED_10Hz_ChR2 <- read_delim(file.path(Datadirectory,"Patch_ChR2_10Hz.csv"), col_names = TRUE)
colnames(LED_10Hz_ChR2)=c("Time","LED_10Hz_ChR2")

#1Hz
segment_data1hz = data.frame(
  x = seq(1.25,20,2),
  xend=seq(1.25,20,2),
  y = -75,
  yend = -72)

OneHz_plot=ggplot() +
  geom_line(data = LED_1Hz_ChR2, aes(x = Time, y = LED_1Hz_ChR2, color = "ChR2"), linewidth=1) +
  geom_line(data = LED_1Hz_Control, aes(x = Time, y = LED_1Hz_Control, color = "Control"), linewidth=1) +
  scale_color_manual(name = "AAV",
                     values = c("ChR2" = "red3", "Control" = "black")) +
  geom_smooth(method = "loess") +
  xlim(0,20) +
  geom_segment(data = segment_data1hz, aes(x = x, y = y, xend = xend, yend = yend), color="steelblue2",linewidth=1) +
  annotate("segment",x = 15, y =40, xend = 17, yend =40, color="black",linewidth=1) + 
  annotate("segment", x = 17, y =40, xend =17, yend =50, color="black",linewidth=1) +
  theme_void() +
  theme(legend.position="none")
OneHz_plot
ggsave(filename = file.path(Outdirectory,"Figure3F_1Hz_scale_1s_10mV.pdf"), plot= OneHz_plot, device='pdf')

OneHz_plot_zoom=ggplot() +
  geom_line(data = LED_1Hz_ChR2, aes(x = Time, y = LED_1Hz_ChR2, color = "ChR2"), linewidth=1) +
  geom_line(data = LED_1Hz_Control, aes(x = Time, y = LED_1Hz_Control, color = "Control"), linewidth=1) +
  scale_color_manual(name = "AAV",
                     values = c("ChR2" = "red3", "Control" = "black")) +
  geom_smooth(method = "loess") +
  xlim(4,6) +
  annotate("segment",x = 5.7, y =40, xend = 5.9, yend =40, color="black",linewidth=1) + 
  annotate("segment", x = 5.9, y =40, xend =5.9, yend =50, color="black",linewidth=1) +
  geom_segment(data = segment_data1hz, aes(x = x, y = y, xend = xend, yend = yend), color="steelblue2",linewidth=1) +
  theme_void() +
  theme(legend.position="none")
OneHz_plot_zoom
ggsave(filename = file.path(Outdirectory,"Figure3F_1Hz_scale_100ms_10mV_zoom_4to6.pdf"), plot= OneHz_plot_zoom, device='pdf')

#5hz
segment_data5hz = data.frame(
  x = seq(1.25,20,.4),
  xend=seq(1.25,20,.4),
  y = -75,
  yend = -72)

FiveHz_plot=ggplot() +
  geom_line(data = LED_5Hz_ChR2, aes(x = Time, y = LED_5Hz_ChR2, color = "ChR2"), linewidth=1) +
  geom_line(data = LED_5Hz_Control, aes(x = Time, y = LED_5Hz_Control, color = "Control"), linewidth=1) +
  scale_color_manual(name = "AAV",
                     values = c("ChR2" = "red3", "Control" = "black")) +
  geom_smooth(method = "loess") +
  xlim(0,20) +
  geom_segment(data = segment_data5hz, aes(x = x, y = y, xend = xend, yend = yend), color="steelblue2",linewidth=1) +
  annotate("segment",x = 15, y =40, xend = 17, yend =40, color="black",linewidth=1) + 
  annotate("segment", x = 17, y =40, xend =17, yend =50, color="black",linewidth=1) +
  theme_void() +
  theme(legend.position="none")
FiveHz_plot
ggsave(filename = file.path(Outdirectory,"Figure3G_5Hz_scale_1s_10mV.pdf"), plot= FiveHz_plot, device='pdf')

FiveHz_plot_zoom=ggplot() +
  geom_line(data = LED_5Hz_ChR2, aes(x = Time, y = LED_5Hz_ChR2, color = "ChR2"), linewidth=1) +
  geom_line(data = LED_5Hz_Control, aes(x = Time, y = LED_5Hz_Control, color = "Control"), linewidth=1) +
  scale_color_manual(name = "AAV",
                     values = c("ChR2" = "red3", "Control" = "black")) +
  geom_smooth(method = "loess") +
  xlim(4,6) +
  annotate("segment",x = 5.7, y =40, xend = 5.9, yend =40, color="black",linewidth=1) + 
  annotate("segment", x = 5.9, y =40, xend =5.9, yend =50, color="black",linewidth=1) +
  geom_segment(data = segment_data5hz, aes(x = x, y = y, xend = xend, yend = yend), color="steelblue2",linewidth=1) +
  theme_void() +
  theme(legend.position="none")
FiveHz_plot_zoom
ggsave(filename = file.path(Outdirectory,"Figure3G_5Hz_scale_100ms_10mV_zoom_4to6.pdf"), plot= FiveHz_plot_zoom, device='pdf')

#10Hz
segment_data10hz = data.frame(
  x = seq(1.25,20,.2),
  xend=seq(1.25,20,.2),
  y = -75,
  yend = -72)

TenHz_plot=ggplot() +
  geom_line(data = LED_10Hz_ChR2, aes(x = Time, y = LED_10Hz_ChR2, color = "ChR2"), linewidth=1) +
  geom_line(data = LED_10Hz_Control, aes(x = Time, y = LED_10Hz_Control, color = "Control"), linewidth=1) +
  scale_color_manual(name = "AAV",
                     values = c("ChR2" = "red3", "Control" = "black")) +
  geom_smooth(method = "loess") +
  xlim(0,20) +
  geom_segment(data = segment_data10hz, aes(x = x, y = y, xend = xend, yend = yend), color="steelblue2",linewidth=1) +
  annotate("segment",x = 15, y =40, xend = 17, yend =40, color="black",linewidth=1) + 
  annotate("segment", x = 17, y =40, xend =17, yend =50, color="black",linewidth=1) +
  theme_void() +
  theme(legend.position="none")
TenHz_plot
ggsave(filename = file.path(Outdirectory,"Figure3H_10Hz_scale_1s_10mV.pdf"), plot= TenHz_plot, device='pdf')

TenHz_plot_zoom=ggplot() +
  geom_line(data = LED_10Hz_ChR2, aes(x = Time, y = LED_10Hz_ChR2, color = "ChR2"), linewidth=1) +
  geom_line(data = LED_10Hz_Control, aes(x = Time, y = LED_10Hz_Control, color = "Control"), linewidth=1) +
  scale_color_manual(name = "AAV",
                     values = c("ChR2" = "red3", "Control" = "black")) +
  geom_smooth(method = "loess") +
  xlim(4,6) +
  annotate("segment",x = 5.7, y =40, xend = 5.9, yend =40, color="black",linewidth=1) + 
  annotate("segment", x = 5.9, y =40, xend =5.9, yend =50, color="black",linewidth=1) +
  geom_segment(data = segment_data10hz, aes(x = x, y = y, xend = xend, yend = yend), color="steelblue2",linewidth=1) +
  theme_void() +
  theme(legend.position="none")
TenHz_plot_zoom
ggsave(filename = file.path(Outdirectory,"Figure3H_10Hz_scale_100ms_10mV_zoom_4to6.pdf"), plot= TenHz_plot_zoom, device='pdf')
