# Bind sensitivity and exposure
# for a first bivariate map

rm(list=ls())


library(tidyverse)


# select degree resolution
deg = "01" # can be "1" or "01"

# select the type of normalization
norm = "log" # can be "rank", "log", or "lin"

# load exposure 
expo <- readRDS(paste0("Output/Exposure/RISK_24_expo_norm_", norm,"_r", deg))
# load sensitivity
sensi <- readRDS(paste0("Output/Sensitivity/RISK_33_sensitivity_norm_", norm,"_r", deg))



####### Normalized exposure ########

# final expo = sum of three exposure components
plot(expo %>% dplyr::ungroup() %>% dplyr::select(-c(x,y)))

expo_fin <- expo %>% ungroup() %>%
  mutate(expo_tot = SR_tot_ias + range_med + med_iasa_tot) %>%
  mutate(expo_tot_norm = (expo_tot/max(expo_tot))) %>%
  dplyr::select(x, y, expo_tot_norm)

quantile(expo_fin$expo_tot_norm, c(0.33, 0.66))

hist(expo_fin$expo_tot_norm)

# plot final expo
ggplot(data = expo_fin) +
  geom_raster(aes(x = x, y = y, fill = expo_tot_norm)) +
  scale_fill_gradient(
    low = "green", 
    high = "red")+
  theme_classic()

# plot each component
# richness
ggplot(data = expo) +
  geom_raster(aes(x = x, y = y, fill = SR_tot_ias)) +
  scale_fill_gradient(
    low = "green", 
    high = "red")+
  theme_classic()

# range
ggplot(data = expo) +
  geom_raster(aes(x = x, y = y, fill = range_med)) +
  scale_fill_gradient(
    low = "green", 
    high = "red")+
  theme_classic()

# strength (nb of native associated)
ggplot(data = expo) +
  geom_raster(aes(x = x, y = y, fill = med_iasa_tot)) +
  scale_fill_gradient(
    low = "green", 
    high = "red")+
  theme_classic()



####### Normalized sensitivity ########

mcor = cor(sensi)
ggcorrplot(mcor)


sensi_fin <- sensi %>% dplyr::select(x, y, SR_tot_bmr, SR_iasa_bmr)

nrow(sensi_fin %>% filter(SR_iasa_bmr==0))

summary(sensi_fin)

ggplot(data = sensi_fin) +
  geom_raster(aes(x = x, y = y, fill = SR_iasa_bmr)) +
  scale_fill_gradient(
    low = "green", 
    high = "red")+
  theme_classic()



#### Combine exposure and sensitivity #####

# create a theme map
theme_map <- function(...) {
  theme_minimal() +
    theme(
      text = element_text(family = "Ubuntu Regular", color = "#22211d"),
      # remove all axes
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      # add a subtle grid
      panel.grid.major = element_line(color = "#ebebe5", size = 0.2),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "#f5f5f2", color = NA), 
      plot.margin = unit(c(.5, .5, .2, .5), "cm"),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "#f5f5f2", color = NA), 
      panel.spacing = unit(c(-.1, 0.2, .2, 0.2), "cm"),
      legend.background = element_rect(fill = "#f5f5f2", color = NA),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 11, hjust = 0, color = "#4e4d47"),
      plot.title = element_text(size = 16, hjust = 0.5, color = "#4e4d47"),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "#4e4d47", 
                                   margin = margin(b = -0.1, 
                                                   t = -0.1, 
                                                   l = 2, 
                                                   unit = "cm"), 
                                   debug = F),
      plot.caption = element_text(size = 9, 
                                  hjust = .5, 
                                  margin = margin(t = 0.2, 
                                                  b = 0, 
                                                  unit = "cm"), 
                                  color = "#939184"),
      ...
    )
}


# data 
expo_sensi <- full_join(expo_fin, sensi_fin, by = c("x","y"))
expo_sensi[is.na(expo_sensi)] <- 0

# create 3 buckets for expo

#when deg == 1
if(deg=="1"){
  probs_expo = c(0, 0.5, 0.75, 1)
}
if(deg=="01"){
  probs_expo = c(0, 0.88, 0.95, 1)
}

quantiles_expo <- expo_sensi %>%
  pull(expo_tot_norm) %>%
  quantile(probs = probs_expo)

# create 3 buckets for mean income
quantiles_sensi <- expo_sensi %>%
  pull(SR_iasa_bmr) %>%
  quantile(probs = seq(0, 1, length.out = 4))

# create color scale that encodes two variables
# red for gini and blue for mean income
# the special notation with gather is due to readibility reasons
bivariate_color_scale <- tibble(
  "3 - 3" = "#3F2949", # high inequality, high income
  "2 - 3" = "#435786",
  "1 - 3" = "#4885C1", # low inequality, high income
  "3 - 2" = "#77324C",
  "2 - 2" = "#806A8A", # medium inequality, medium income
  "1 - 2" = "#89A1C8",
  "3 - 1" = "#AE3A4E", # high inequality, low income
  "2 - 1" = "#BC7C8F",
  "1 - 1" = "#CABED0" # low inequality, low income
) %>%
  gather("group", "fill")

# cut into groups defined above and join fill
expo_sensi %<>%
  mutate(
    expo_quantiles = cut(
      expo_tot_norm,
      breaks = quantiles_expo,
      include.lowest = TRUE
    ),
    sensi_quantiles = cut(
      SR_iasa_bmr,
      breaks = quantiles_sensi,
      include.lowest = TRUE
    ),
    # by pasting the factors together as numbers we match the groups defined
    # in the tibble bivariate_color_scale
    group = paste(
      as.numeric(expo_quantiles), "-",
      as.numeric(sensi_quantiles)
    )
  ) %>%
  # we now join the actual hex values per "group"
  # so each municipality knows its hex value based on the his gini and avg
  # income value
  left_join(bivariate_color_scale, by = "group")



# Map

map <- ggplot(data = expo_sensi)+
  # color cells according to their expo/sensi combination
  geom_raster(aes( x = x, y = y, fill = fill)) +
  scale_fill_identity() +
  # add titles
  labs(x = NULL,
       y = NULL,
       title = "Global risk of amniotes facing biological invasions")+
  # add the theme
  theme_map()

map


# legend
bivariate_color_scale %<>%
  separate(group, into = c("expo", "sensi"), sep = " - ") %>%
  mutate(expo = as.integer(expo),
         sensi = as.integer(sensi))

legend <- ggplot() +
  geom_tile(data = bivariate_color_scale,
    mapping = aes(x = expo, y = sensi, fill = fill)) +
  scale_fill_identity() +
  labs(x = "Exposure ⟶️",
       y = "Sensitivity ⟶️") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank())+
  # make font small enough
  theme(axis.title = element_text(size = 8)) +
  # quadratic tiles
  coord_fixed()

legend
library(cowplot)
ggdraw() +
  draw_plot(map, 0, 0, 1, 1) +
  draw_plot(legend, 0.04, 0.85, 0.15, 0.15)
