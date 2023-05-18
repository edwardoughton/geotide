library(tidyverse)
library(ggpubr)

iso3 = 'BFA'

folder <- dirname(rstudioapi::getSourceEditorContext()$path)
filename = paste('acled_',iso3,'.csv', sep="")
path = file.path(folder, '..', 'results', iso3,  filename)
data = read.csv(path)
data$count = 1

filename = 'mobile_codes.csv'
path = file.path(folder, '..', 'data', 'raw',  filename)
mobile_codes = read.csv(path)
mobile_codes = mobile_codes[(mobile_codes$iso3 == iso3),]
mobile_codes = select(mobile_codes, mnc, network)

data = merge(data, mobile_codes, by.x='net', by.y='mnc')

data$radio = factor(data$radio,
                              levels=c("GSM","UMTS","LTE"),
                              labels=c("2G GSM","3G UMTS","4G LTE")
)

data$infra = factor(data$infra,
    levels=c("antenna","installation","installation (mention)","radio station", "telecommunications station"),
    labels=c("Antenna","Installation","Installation","Radio Station", "Telecom\nStation")
)

data$damage_inf = factor(data$damage_inf,
    levels=c("arson","battle","destruction","ied","looting", "ransacking", "sabotage", "sabotage_removed"),
    labels=c("Arson","Battle","Destruction","IED","Looting", "Ransacking", "Sabotage", "Sabotage (Removed)")
)

data$actor1 = factor(data$actor1,
     levels=c("JNIM: Group for Support of Islam and Muslims", 
              "Rioters (Burkina Faso)",
              "Unidentified Armed Group (Burkina Faso)"),
     labels=c("JNIM: Group\nfor Support of\nIslam and Muslims", 
              "Rioters\n(Burkina Faso)",
              "Unidentified\nArmed Group\n(Burkina Faso)")
)


data$severity_i = factor(data$severity_i,
                         levels=c("burned", "destroyed","looting","partially burned","unknown"),
                         labels=c("Burned", "Destroyed","Looting","Partially Burned","Unknown")
)

data = select(data, year, country_1, infra, sub_event_, damage_inf, 
              severity_i, actor1, radio, network, count)

subset = data %>%
  group_by(year, country_1, infra, sub_event_, damage_inf, severity_i, actor1, radio, network) %>%
  summarise(
    count = sum(count),
  )

subset2 = subset[(subset$count > 60),] 
subset3 = subset[(subset$count <= 60),] 
subset2$count = subset2$count / 8
subset = rbind(subset2, subset3)

remove(subset2, subset3)

ggplot(subset, aes(x=year, y=count, fill=radio)) + 
  geom_bar(stat="identity") +
  theme(legend.position = 'bottom',
      axis.text.x = element_text(angle=45, hjust=1)) +
  labs(colour=NULL,
       title = "Damaged Cell Count Estimate: Burkina Faso",
       subtitle = "Reported by Year and Generation.", 
       x = "Year", y = "Estimated Cell Count", fill="Generation") +
  theme(panel.spacing = unit(0.6, "lines")) + 
  expand_limits(y=0) +
  guides(fill=guide_legend(ncol=4, title='Scenario')) +
  scale_fill_viridis_d(direction=1) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(2011, 2023, by=1)) +
  scale_y_continuous(expand = c(0, 0), limits=c(0, 620))

path = file.path(folder, 'figures', paste(iso3,'plot1.png',sep=''))
ggsave(path, units="in", width=6, height=4, dpi=300)

plot2 = ggplot(subset, aes(x=sub_event_, y=count, fill=infra)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  theme(legend.position = 'bottom',
      axis.text.x = element_text(angle=45, hjust=1)) +
  labs(colour=NULL,
       title = "Damaged Cell Count Estimate: Burkina Faso",
       subtitle = "Reported by Sub-Event and Generation.", 
       x = "Sub-Event", y = "Estimated Cell Count", fill="Generation") +
  theme(panel.spacing = unit(0.6, "lines")) + 
  expand_limits(y=0) +
  guides(fill=guide_legend(ncol=6, title='Infra Type')) +
  scale_fill_viridis_d(direction=1) +
  # scale_x_discrete(expand = c(0, 0.15)) +
  scale_y_continuous(expand = c(0, 0), limits=c(0, 620))

path = file.path(folder, 'figures', paste(iso3,'plot2.png',sep=''))
ggsave(path, units="in", width=6, height=4, dpi=300)

plot3 = ggplot(subset, aes(x=damage_inf, y=count, fill=infra)) + 
  geom_bar(stat="identity") +
  coord_flip() + 
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle=45, hjust=1)) +
  labs(colour=NULL,
       title = "Damaged Cell Count Estimate: Burkina Faso",
       subtitle = "Reported by Damage Type and Infrastructure Type", 
       x = "Damage Type", y = "Estimated Cell Count", fill="Generation") +
  theme(panel.spacing = unit(0.6, "lines")) + 
  expand_limits(y=0) +
  guides(fill=guide_legend(ncol=4, title='Infra Type')) +
  scale_fill_viridis_d(direction=1) +
  # scale_x_discrete(expand = c(0, 0.15)) +
  scale_y_continuous(expand = c(0, 0), limits=c(0, 620))

path = file.path(folder, 'figures', paste(iso3,'plot3.png',sep=''))
ggsave(path, units="in", width=6, height=4, dpi=300)

plot4 = ggplot(subset, aes(x=severity_i, y=count, fill=radio)) + 
  geom_bar(stat="identity") +
coord_flip() + 
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle=45, hjust=1)) +
  labs(colour=NULL,
       title = "Damaged Cell Count Estimate: Burkina Faso",
       subtitle = "Reported by Damage Severity and Generation.", 
       x = "Damage Severity", y = "Estimated Cell Count", fill="Generation") +
  theme(panel.spacing = unit(0.6, "lines")) + 
  expand_limits(y=0) +
  guides(fill=guide_legend(ncol=4, title='Generation')) +
  scale_fill_viridis_d(direction=1) +
  # scale_x_discrete(expand = c(0, 0.15)) +
  scale_y_continuous(expand = c(0, 0), limits=c(0, 620))

path = file.path(folder, 'figures', paste(iso3,'plot4.png',sep=''))
ggsave(path, units="in", width=6, height=4, dpi=300)

plot5 = ggplot(subset, aes(x=actor1, y=count, fill=radio)) + 
  geom_bar(stat="identity") +
  coord_flip() + 
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle=45, hjust=1)) +
  labs(colour=NULL,
       title = "Damaged Cell Count Estimate: Burkina Faso",
       subtitle = "Reported by Damage Actor and Generation.", 
       x = "Damage Actor", y = "Estimated Cell Count", fill="Generation") +
  theme(panel.spacing = unit(0.6, "lines")) + 
  expand_limits(y=0) +
  guides(fill=guide_legend(ncol=4, title='Generation')) +
  scale_fill_viridis_d(direction=1) +
  # scale_x_discrete(expand = c(0, 0.15)) +
  scale_y_continuous(expand = c(0, 0), limits=c(0, 620))

path = file.path(folder, 'figures', paste(iso3,'plot5.png',sep=''))
ggsave(path, units="in", width=6, height=4, dpi=300)

plot6 = ggplot(subset, aes(x=network, y=count, fill=radio)) + 
  geom_bar(stat="identity") +
  coord_flip() + 
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle=45, hjust=1)) +
  labs(colour=NULL,
       title = "Damaged Cell Count Estimate: Burkina Faso",
       subtitle = "Reported by Mobile Network Operator and Generation.", 
       x = "Mobile Network Operator", y = "Estimated Cell Count", fill="Generation") +
  theme(panel.spacing = unit(0.6, "lines")) + 
  expand_limits(y=0) +
  guides(fill=guide_legend(ncol=4, title='Generation')) +
  scale_fill_viridis_d(direction=1) +
  # scale_x_discrete(expand = c(0, 0.15)) +
  scale_y_continuous(expand = c(0, 0), limits=c(0, 620))

path = file.path(folder, 'figures', paste(iso3,'plot6.png',sep=''))
ggsave(path, units="in", width=6, height=4, dpi=300)

