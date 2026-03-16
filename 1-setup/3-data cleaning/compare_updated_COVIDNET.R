rm(list = ls())

library(readr)

og <- read_csv("data/clean-data/weekly-incidence-estimates-US-validationPeriod-modelInput.csv")

updated <- read_csv("data/clean-data/weekly-incidence-estimates-US-validationPeriod-modelInput-updated.csv")

#Plot
plot_data <- og %>%
  pivot_longer(
    cols = -c(age_group, week, weeks_since),  # all columns except these
    names_to = "variable",
    values_to = "value"
  )

# Prepare your second dataframe
plot_data2 <- updated %>%
  pivot_longer(
    cols = -c(age_group, week, weeks_since),
    names_to = "variable",
    values_to = "value"
  ) %>%
  filter(week <= max(plot_data$week))

# Create the plot with both datasets
plot_data %>% 
  filter(variable %in% c("healthy_inc")) %>%
  ggplot(aes(week, value, color = age_group)) + 
  geom_line(linetype = "solid") + 
  geom_point() +
  # Add second dataframe with dashed lines
  geom_line(data = plot_data2 %>% filter(variable %in% c("healthy_inc")),
            linetype = "dashed") +
  geom_point(data = plot_data2 %>% filter(variable %in% c("healthy_inc"))) +
  #ylim(0, 30) +
  xlab("Time") +
  ylab("Weekly incidence (per 100,000 persons)") +
  scale_x_date(date_labels = ("%b-%Y"),
               limits = c(as.Date("2023-07-01"), as.Date("2025-05-01")),
               breaks = "1 months") +
  labs(color = "Age Group") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = 'none') +
  ggtitle(#"Weekly Incidence of COVID Hospitalizations in US\n
    "Risk Group: Healthy") + 
  facet_wrap(~age_group, scales = "free_y")

plot_data %>% 
  filter(variable %in% c("higher_risk_inc")) %>%
  ggplot(aes(week, value, color = age_group, linetype = variable)) + 
  geom_line(linetype = "solid") + 
  geom_point() +
  # Add second dataframe with dashed lines
  geom_line(data = plot_data2 %>% filter(variable %in% c("higher_risk_inc")),
            linetype = "dashed") +
  geom_point(data = plot_data2 %>% filter(variable %in% c("higher_risk_inc"))) +
  #ylim(0, 60) +
  xlab("Time") +
  ylab("Weekly incidence (per 100,000 persons)") +
  scale_x_date(date_labels = ("%b-%Y"),
               limits = c(as.Date("2023-07-01"), as.Date("2025-05-01")),
               breaks = "1 months") +
  labs(color = "Age Group", linetype = "Risk Group") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = 'none')+
  ggtitle(#"Weekly Incidence of COVID Hospitalizations in US\n
          "Risk Group: Higher risk") + 
  facet_wrap(~age_group, scales = "free_y")

plot_data %>% 
  filter(variable %in% c("immunocompromised_inc")) %>%
  ggplot(aes(week, value, color = age_group, linetype = variable)) + 
  geom_line(linetype = "solid") + 
  geom_point() +
  # Add second dataframe with dashed lines
  geom_line(data = plot_data2 %>% filter(variable %in% c("immunocompromised_inc")),
            linetype = "dashed") +
  geom_point(data = plot_data2 %>% filter(variable %in% c("immunocompromised_inc"))) +
  #ylim(0, 75) +
  xlab("Time") +
  ylab("Weekly incidence (per 100,000 persons)") +
  scale_x_date(date_labels = ("%b-%Y"),
               limits = c(as.Date("2023-07-01"), as.Date("2025-05-01")),
               breaks = "1 months") +
  labs(color = "Age Group", linetype = "Risk Group") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = 'none') +
  ggtitle(#"Weekly Incidence of COVID Hospitalizations in US\n
    "Risk Group: Immunocompromised") + 
  facet_wrap(~age_group, scales = "free_y")

grid.arrange(p1, p2, p3, nrow = 1)

ggplot() + 
  # First dataset - solid lines
  geom_line(data = plot_data %>% filter(variable %in% c("adj_inc")),
            aes(x = week, y = value, color = age_group),
            linetype = "solid") + 
  geom_point(data = plot_data %>% filter(variable %in% c("adj_inc")),
             aes(x = week, y = value, color = age_group)) +
  # Second dataset - dashed lines
  geom_line(data = plot_data2 %>% filter(variable %in% c("adj_inc")),
            aes(x = week, y = value, color = age_group),
            linetype = "dashed") +
  geom_point(data = plot_data2 %>% filter(variable %in% c("adj_inc")),
             aes(x = week, y = value, color = age_group)) +
  #ylim(0, 100) +
  xlab("Time") +
  ylab("Weekly incidence (per 100,000 persons)") +
  scale_x_date(date_labels = ("%b-%Y"),
               limits = c(as.Date("2023-07-01"), as.Date("2025-05-01")),
               breaks = "1 months") +
  labs(color = "Age Group") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ggtitle("Weekly Incidence of COVID Hospitalizations in US")

##################

cumulative_original <- plot_data %>%
  group_by(age_group) %>%
  filter(variable == "adj_inc") %>%
  summarize(cumulative_inc = sum(value, na.rm = TRUE), .groups = 'drop')

cumulative_updated <- plot_data2 %>%
  group_by(age_group) %>%
  filter(variable == "adj_inc") %>%
  summarize(cumulative_inc = sum(value, na.rm = TRUE), .groups = 'drop')
