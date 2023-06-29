library(tidyverse)
library(ggpubr)

folder <- dirname(rstudioapi::getSourceEditorContext()$path)

get_data <- function(iso3) { # create a function with the name my_function
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
}

data1 = get_data('BFA')
data2 = get_data('MLI')
data3 = get_data('NER')

data = rbind(data1, data2, data3)

remove(data1, data2, data3)

data$radio = factor(data$radio,
                              levels=c("GSM","UMTS","LTE"),
                              labels=c("2G GSM","3G UMTS","4G LTE")
)

data$infra = factor(data$infra,
    levels=c("antenna","installation","installation (mention)","radio station", "telecommunications station", "other"),
    labels=c("Antenna","Installation","Installation","Radio Station", "Telecom\nStation", "Other")
)

data$damage_inf = factor(data$damage_inf,
    levels=c("arson","battle","destruction","gunfire","ied","looting", "ransacking", "sabotage", "sabotage_removed", "vac"),
    labels=c("Arson","Battle","Destruction","Gunfire","IED","Looting", "Ransacking", "Sabotage", "Sabotage", "Vac")
)

data$actor1 = factor(data$actor1,
     levels=c("Islamic State and/or JNIM",
              "Islamic State (West Africa) - Greater Sahara Faction",
              "JNIM: Group for Support of Islam and Muslims",
              "Katiba Macina",
              "MUJAO: Movement for Unity and Jihad in West Africa",
              "Mutiny of Military Forces of Mali (2002-2012)",
              "Rioters (Burkina Faso)",
              "Unidentified Armed Group (Burkina Faso)",
              "Unidentified Armed Group (Mali)",
              "Unidentified Communal Militia (Mali)"
              ), 
     labels=c("Islamic State and/or JNIM",
              "Islamic State (West Africa) - Greater Sahara Faction",
              "JNIM: Group for Support of Islam and Muslims",
              "Katiba Macina",
              "MUJAO: Movement for Unity and Jihad in West Africa",
              "Mutiny of Military Forces of Mali (2002-2012)",
              "Rioters (Burkina Faso)",
              "Unidentified Armed Group (Burkina Faso)",
              "Unidentified Armed Group (Mali)",
              "Unidentified Communal Militia (Mali)")
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

subset$year_period = ""
subset$year_period[subset$year < 2022] = 'Pre-2022'
subset$year_period[subset$year >= 2022] = 'Post-2022'

subset$handle = paste(subset$country_1, " (",subset$year_period,")", sep = "")

subset$handle = factor(subset$handle,
       levels=c(
         "Burkina Faso (Pre-2022)", 
         "Burkina Faso (Post-2022)",
         "Mali (Pre-2022)",
         "Mali (Post-2022)",
         "Niger (Pre-2022)",
         "Niger (Post-2022)" 
       )
)

totals <- subset %>%
  group_by(handle, year) %>%
  summarize(count = sum(count))
max_y_value = max(totals$count, na.rm = TRUE)

#######by year
ggplot(subset, aes(x=year, y=count, fill=radio)) +
  geom_col() +
  geom_text(aes(year, count, label = count,  #y=0, 
                fill = NULL), show.legend = FALSE, ##FF0000FF
            size = 2, data = totals, vjust=-.9, hjust=.5) +
  theme(legend.position = 'bottom',
      axis.text.x = element_text(angle=45, hjust=1)) +
  labs(colour=NULL,
       title = "Damaged Cell Count Estimate",
       subtitle = "Reported by Year, Generation, Country and Time Period.",
       x = "Year", y = "Estimated Cell Count", fill="Generation") +
  theme(panel.spacing = unit(0.6, "lines")) +
  expand_limits(y=0) +
  guides(fill=guide_legend(ncol=4, title='Generation')) +
  scale_fill_viridis_d(direction=1) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(2011, 2023, by=1)) +
  scale_y_continuous(expand = c(0, 0), limits=c(0, max_y_value + 20)) +
  facet_wrap(~handle, ncol=2)

path = file.path(folder, 'figures', 'year_generation.png')
ggsave(path, units="in", width=8, height=8, dpi=300)

ggplot(subset, aes(x=year, y=count, fill=sub_event_)) +
  geom_bar(stat="identity") +
  geom_text(aes(year, count, label = count,  #y=0, 
                fill = NULL), show.legend = FALSE, ##FF0000FF
            size = 2, data = totals, vjust=-.9, hjust=.5) +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle=45, hjust=1)) +
  labs(colour=NULL,
       title = "Damaged Cell Count Estimate",
       subtitle = "Reported by Year, Sub-Event, Country and Time Period.",
       x = "Year", y = "Estimated Cell Count", fill="Sub-Event") +
  theme(panel.spacing = unit(0.6, "lines")) +
  expand_limits(y=0) +
  guides(fill=guide_legend(ncol=4, title='Sub-Event')) +
  scale_fill_viridis_d(direction=1) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(2011, 2023, by=1)) +
  scale_y_continuous(expand = c(0, 0), limits=c(0, max_y_value+20)) +
  facet_wrap(~handle, ncol=2)

path = file.path(folder, 'figures', 'year_sub-event.png')
ggsave(path, units="in", width=8, height=8, dpi=300)

ggplot(subset, aes(x=year, y=count, fill=damage_inf)) +
  geom_bar(stat="identity") +
  geom_text(aes(year, count, label = count,  #y=0, 
                fill = NULL), show.legend = FALSE, ##FF0000FF
            size = 2, data = totals, vjust=-.9, hjust=.5) +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle=45, hjust=1)) +
  labs(colour=NULL,
       title = "Damaged Cell Count Estimate",
       subtitle = "Reported by Year, Damage, Country and Time Period.",
       x = "Year", y = "Estimated Cell Count", fill="Damage_Inf") +
  theme(panel.spacing = unit(0.6, "lines")) +
  expand_limits(y=0) +
  guides(fill=guide_legend(ncol=6, title='Damage')) +
  scale_fill_viridis_d(direction=1) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(2011, 2023, by=1)) +
  scale_y_continuous(expand = c(0, 0), limits=c(0, max_y_value+20)) +
  facet_wrap(~handle, ncol=2)

path = file.path(folder, 'figures', 'year_damage.png')
ggsave(path, units="in", width=8, height=8, dpi=300)

ggplot(subset, aes(x=year, y=count, fill=actor1)) +
  geom_bar(stat="identity") +
  geom_text(aes(year, count, label = count,  #y=0, 
                fill = NULL), show.legend = FALSE, ##FF0000FF
            size = 2, data = totals, vjust=-.9, hjust=.5) +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle=45, hjust=1)) +
  labs(colour=NULL,
       title = "Damaged Cell Count Estimate",
       subtitle = "Reported by Year, Actor, Country and Time Period.",
       x = "Year", y = "Estimated Cell Count", fill="Damage_Inf") +
  theme(panel.spacing = unit(0.6, "lines")) +
  expand_limits(y=0) +
  guides(fill=guide_legend(ncol=2, title='Actor')) +
  scale_fill_viridis_d(direction=1) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(2011, 2023, by=1)) +
  scale_y_continuous(expand = c(0, 0), limits=c(0, max_y_value+20)) +
  facet_wrap(~handle, ncol=2)

path = file.path(folder, 'figures', 'year_actor.png')
ggsave(path, units="in", width=8, height=8, dpi=300)

#######by sub-event
totals <- subset %>%
  group_by(handle, sub_event_) %>%
  summarize(count = sum(count))
max_y_value = max(totals$count, na.rm = TRUE)

ggplot(subset, aes(x=sub_event_, y=count, fill=infra)) +
  geom_bar(stat="identity") +
  geom_text(aes(sub_event_, count, label = count,  #y=0,
                fill = NULL), show.legend = FALSE, ##FF0000FF
            size = 2, data = totals, vjust=.3, hjust=-.5) +
  coord_flip() +
  theme(legend.position = 'bottom',
      axis.text.x = element_text(angle=45, hjust=1)) +
  labs(colour=NULL,
       title = "Damaged Cell Count Estimate",
       subtitle = "Reported by Sub-Event, Infrastructure Type, Country and Time Period.",
       x = "Sub-Event", y = "Estimated Cell Count", fill="Infra Type") +
  theme(panel.spacing = unit(0.6, "lines")) +
  expand_limits(y=0) +
  guides(fill=guide_legend(ncol=6, title='Infra Type')) +
  scale_fill_viridis_d(direction=1) +
  # scale_x_discrete(expand = c(0, 0.15)) +
  scale_y_continuous(expand = c(0, 0), limits=c(0, max_y_value+10)) +
  facet_wrap(~handle, ncol=2)

path = file.path(folder, 'figures', 'sub-event_infra-type.png')
ggsave(path, units="in", width=8, height=8, dpi=300)

ggplot(subset, aes(x=sub_event_, y=count, fill=actor1)) +
  geom_bar(stat="identity") +
  geom_text(aes(sub_event_, count, label = count,  #y=0,
                fill = NULL), show.legend = FALSE, ##FF0000FF
            size = 2, data = totals, vjust=.3, hjust=-.5) +
  coord_flip() +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle=45, hjust=1)) +
  labs(colour=NULL,
       title = "Damaged Cell Count Estimate",
       subtitle = "Reported by Sub-Event, Actor, Country and Time Period.",
       x = "Sub-Event", y = "Estimated Cell Count", fill="Infra Type") +
  theme(panel.spacing = unit(0.6, "lines")) +
  expand_limits(y=0) +
  guides(fill=guide_legend(ncol=2, title='Actor')) +
  scale_fill_viridis_d(direction=1) +
  # scale_x_discrete(expand = c(0, 0.15)) +
  scale_y_continuous(expand = c(0, 0), limits=c(0, max_y_value+10)) +
  facet_wrap(~handle, ncol=2)

path = file.path(folder, 'figures', 'sub-event_actor.png')
ggsave(path, units="in", width=8, height=8, dpi=300)

ggplot(subset, aes(x=sub_event_, y=count, fill=damage_inf)) +
  geom_bar(stat="identity") +
  geom_text(aes(sub_event_, count, label = count,  #y=0,
                fill = NULL), show.legend = FALSE, ##FF0000FF
            size = 2, data = totals, vjust=.3, hjust=-.5) +
  coord_flip() +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle=45, hjust=1)) +
  labs(colour=NULL,
       title = "Damaged Cell Count Estimate",
       subtitle = "Reported by Sub-Event, Damage, Country and Time Period.",
       x = "Sub-Event", y = "Estimated Cell Count", fill="Infra Type") +
  theme(panel.spacing = unit(0.6, "lines")) +
  expand_limits(y=0) +
  guides(fill=guide_legend(ncol=6, title='Damage')) +
  scale_fill_viridis_d(direction=1) +
  # scale_x_discrete(expand = c(0, 0.15)) +
  scale_y_continuous(expand = c(0, 0), limits=c(0, max_y_value+10)) +
  facet_wrap(~handle, ncol=2)

path = file.path(folder, 'figures', 'sub-event_damage-inf.png')
ggsave(path, units="in", width=8, height=8, dpi=300)


#####damage type
totals <- subset %>%
  group_by(handle, damage_inf) %>%
  summarize(count = sum(count))
max_y_value = max(totals$count, na.rm = TRUE)

ggplot(subset, aes(x=damage_inf, y=count, fill=infra)) +
  geom_bar(stat="identity") +
  geom_text(aes(damage_inf, count, label = count,  #y=0,
                fill = NULL), show.legend = FALSE, ##FF0000FF
            size = 2, data = totals, vjust=.4, hjust=-.7) +
  coord_flip() +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle=45, hjust=1)) +
  labs(colour=NULL,
       title = "Damaged Cell Count Estimate",
       subtitle = "Reported by Damage, Infrastructure Type, Country and Time Period",
       x = "Damage Type", y = "Estimated Cell Count", fill="Infra Type") +
  theme(panel.spacing = unit(0.6, "lines")) +
  expand_limits(y=0) +
  guides(fill=guide_legend(ncol=4, title='Infra Type')) +
  scale_fill_viridis_d(direction=1) +
  # scale_x_discrete(expand = c(0, 0.15)) +
  scale_y_continuous(expand = c(0, 0), limits=c(0, max_y_value+7)) +
  facet_wrap(~handle, ncol=2)

path = file.path(folder, 'figures', 'damage_infra-type.png')
ggsave(path, units="in", width=8, height=8, dpi=300)


ggplot(subset, aes(x=damage_inf, y=count, fill=actor1)) +
  geom_bar(stat="identity") +
  geom_text(aes(damage_inf, count, label = count,  #y=0,
                fill = NULL), show.legend = FALSE, ##FF0000FF
            size = 2, data = totals, vjust=.4, hjust=-.7) +
  coord_flip() +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle=45, hjust=1)) +
  labs(colour=NULL,
       title = "Damaged Cell Count Estimate",
       subtitle = "Reported by Damage, Actor, Country and Time Period",
       x = "Damage Type", y = "Estimated Cell Count", fill="Actor") +
  theme(panel.spacing = unit(0.6, "lines")) +
  expand_limits(y=0) +
  guides(fill=guide_legend(ncol=2, title='Actor')) +
  scale_fill_viridis_d(direction=1) +
  # scale_x_discrete(expand = c(0, 0.15)) +
  scale_y_continuous(expand = c(0, 0), limits=c(0, max_y_value+7)) +
  facet_wrap(~handle, ncol=2)

path = file.path(folder, 'figures', 'damage_actor.png')
ggsave(path, units="in", width=8, height=8, dpi=300)

ggplot(subset, aes(x=damage_inf, y=count, fill=sub_event_)) +
  geom_bar(stat="identity") +
  geom_text(aes(damage_inf, count, label = count,  #y=0,
                fill = NULL), show.legend = FALSE, ##FF0000FF
            size = 2, data = totals, vjust=.4, hjust=-.7) +
  coord_flip() +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle=45, hjust=1)) +
  labs(colour=NULL,
       title = "Damaged Cell Count Estimate",
       subtitle = "Reported by Damage, Sub-Event, Country and Time Period",
       x = "Damage Type", y = "Estimated Cell Count", fill="Sub-Event") +
  theme(panel.spacing = unit(0.6, "lines")) +
  expand_limits(y=0) +
  guides(fill=guide_legend(ncol=4, title='Sub-Event')) +
  scale_fill_viridis_d(direction=1) +
  # scale_x_discrete(expand = c(0, 0.15)) +
  scale_y_continuous(expand = c(0, 0), limits=c(0, max_y_value+7)) +
  facet_wrap(~handle, ncol=2)

path = file.path(folder, 'figures', 'damage_sub-event.png')
ggsave(path, units="in", width=8, height=8, dpi=300)


# plot4 = ggplot(subset, aes(x=severity_i, y=count, fill=radio)) +
#   geom_bar(stat="identity") +
# coord_flip() +
#   theme(legend.position = 'bottom',
#         axis.text.x = element_text(angle=45, hjust=1)) +
#   labs(colour=NULL,
#        title = "Damaged Cell Count Estimate: Burkina Faso",
#        subtitle = "Reported by Damage Severity and Generation.",
#        x = "Damage Severity", y = "Estimated Cell Count", fill="Generation") +
#   theme(panel.spacing = unit(0.6, "lines")) +
#   expand_limits(y=0) +
#   guides(fill=guide_legend(ncol=4, title='Generation')) +
#   scale_fill_viridis_d(direction=1) +
#   # scale_x_discrete(expand = c(0, 0.15)) +
#   scale_y_continuous(expand = c(0, 0), limits=c(0, 620))
# 
# path = file.path(folder, 'figures', paste(iso3,'plot4.png',sep=''))
# ggsave(path, units="in", width=6, height=4, dpi=300)
# 
# plot5 = ggplot(subset, aes(x=actor1, y=count, fill=radio)) +
#   geom_bar(stat="identity") +
#   coord_flip() +
#   theme(legend.position = 'bottom',
#         axis.text.x = element_text(angle=45, hjust=1)) +
#   labs(colour=NULL,
#        title = "Damaged Cell Count Estimate: Burkina Faso",
#        subtitle = "Reported by Damage Actor and Generation.",
#        x = "Damage Actor", y = "Estimated Cell Count", fill="Generation") +
#   theme(panel.spacing = unit(0.6, "lines")) +
#   expand_limits(y=0) +
#   guides(fill=guide_legend(ncol=4, title='Generation')) +
#   scale_fill_viridis_d(direction=1) +
#   # scale_x_discrete(expand = c(0, 0.15)) +
#   scale_y_continuous(expand = c(0, 0), limits=c(0, 620))
# 
# path = file.path(folder, 'figures', paste(iso3,'plot5.png',sep=''))
# ggsave(path, units="in", width=6, height=4, dpi=300)
# 
# plot6 = ggplot(subset, aes(x=network, y=count, fill=radio)) +
#   geom_bar(stat="identity") +
#   coord_flip() +
#   theme(legend.position = 'bottom',
#         axis.text.x = element_text(angle=45, hjust=1)) +
#   labs(colour=NULL,
#        title = "Damaged Cell Count Estimate: Burkina Faso",
#        subtitle = "Reported by Mobile Network Operator and Generation.",
#        x = "Mobile Network Operator", y = "Estimated Cell Count", fill="Generation") +
#   theme(panel.spacing = unit(0.6, "lines")) +
#   expand_limits(y=0) +
#   guides(fill=guide_legend(ncol=4, title='Generation')) +
#   scale_fill_viridis_d(direction=1) +
#   # scale_x_discrete(expand = c(0, 0.15)) +
#   scale_y_continuous(expand = c(0, 0), limits=c(0, 620))
# 
# path = file.path(folder, 'figures', paste(iso3,'plot6.png',sep=''))
# ggsave(path, units="in", width=6, height=4, dpi=300)
# 
