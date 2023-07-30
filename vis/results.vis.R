library(tidyverse)
library(ggpubr)

folder <- dirname(rstudioapi::getSourceEditorContext()$path)

get_data <- function(iso3) { # create a function with the name my_function
  filename = paste('cba_site_results_',iso3,'.csv', sep="")
  path = file.path(folder, '..', 'results', iso3,  filename)
  data = read.csv(path)
}

data1 = get_data('BFA')
data1$iso3 = 'BFA'
data2 = get_data('MLI')
data2$iso3 = 'MLI'
data3 = get_data('NER')
data3$iso3 = 'NER'

data = rbind(data1, data2, data3)

remove(data1, data2, data3)

data$bca = data$value_at_risk / data$protection_costs

data$iso3 = factor(data$iso3,
                    levels=c("BFA","MLI","NER"),
                    labels=c("BFA","MLI","NER")
)
data$value_at_risk = data$value_at_risk /1e6

# Density plots
plot1 = ggplot(data, aes(x=iso3, y=value_at_risk, colour=iso3)) + geom_boxplot() +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle=0, hjust=.5)) +
  labs(colour=NULL,
       title = "(A) Distribution of Value at Risk Estimates",
       # subtitle = "Reported by Year, Sub-Event, Country and Time Period.",
       x = "Country ISO3 Code", y = "Value@Risk (US$)") +
  theme(panel.spacing = unit(0.6, "lines")) +
  expand_limits(y=0) +
  guides(fill=guide_legend(ncol=4, title='Sub-Event')) +
  scale_fill_viridis_d(direction=1) +
  scale_y_continuous(expand = c(0, 0), limits=c(0, .7)) 

# Density plots
plot2 = ggplot(data, aes(x=iso3, y=bca, colour=iso3)) + geom_boxplot() +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle=0, hjust=.5)) +
  labs(colour=NULL,
       title = "(B) Distribution of CBR Estimates",
       # subtitle = "Reported by Year, Sub-Event, Country and Time Period.",
       x = "Country ISO3 Code", y = "Benefit-Cost Ratio") +
  theme(panel.spacing = unit(0.6, "lines")) +
  expand_limits(y=0) +
  guides(fill=guide_legend(ncol=4, title='Sub-Event')) +
  scale_fill_viridis_d(direction=1) +
  scale_y_continuous(expand = c(0, 0), limits=c(0, 100)) 

panel <- ggarrange(plot1, plot2,
                   ncol = 2, nrow = 1, align = c("hv"),
                   common.legend = FALSE#,
                   # legend='bottom',
                   # heights=c(1,1,.85)
                   )

path = file.path(folder, 'figures', 'test_panel.png')
dir.create(file.path(folder, 'figures'), showWarnings = FALSE)

ggsave(
  'panel.png',
  plot = last_plot(),
  device = "png",
  path=file.path(folder, 'figures'),
  units = c("in"),
  width = 9,
  height = 5,
  bg="white"
)



