
# I always start a script with this line, just to make sure the workspace is empty
rm(list = ls())

# If you don’t have INLA you can install it with this
# install.packages(“INLA”,repos=c(getOption(“repos”),
# INLA=”https://inla.r-inla-download.org/R/stable”), dep=TRUE)

# load in the packages
require(tidyverse)
require(INLA)

# load in the data
df <- readr::read_csv("~/PitcherList/INLA/Model_1_data_xwOBAcon.csv") %>% 
  dplyr::filter(!is.na(launch_angle))
  

# calculate the number of combinations of EV and LA in the data
df %>% 
    dplyr::filter(!is.na(launch_speed) | !is.na(launch_angle)) %>% 
    group_by(launch_speed, launch_angle) %>% 
    tally %>%
    arrange(desc(n))


# randomly sample 35 points from each of these combinations of launch_speed and launch_angle
temp <- df %>% 
  dplyr::filter(launch_speed == 80 & launch_angle == 69 |
                  launch_speed == 82.9 & launch_angle == -20.7 |
                  launch_speed == 41 & launch_angle == -39 |
                  launch_speed == 90.3 & launch_angle == -17.3 |
                  launch_speed == 89.2 & launch_angle == 39.3) %>% 
  group_by(launch_speed, launch_angle) %>% 
  sample_n(35) 


# remove all combinations where the values are likely more than expected
df <- df %>% 
  dplyr::filter(!(launch_speed == 80 & launch_angle == 69 |
                  launch_speed == 82.9 & launch_angle == -20.7 |
                  launch_speed == 41 & launch_angle == -39 |
                  launch_speed == 90.3 & launch_angle == -17.3 |
                  launch_speed == 89.2 & launch_angle == 39.3))


# bind the orifinal data and the sub-sampled data back together
df <- bind_rows(df, temp)


# now those 5 grouping that had a lot of observations are all dropeed to 35 points in each bin
df %>% 
  dplyr::filter(!is.na(launch_speed) | !is.na(launch_angle)) %>% 
  group_by(launch_speed, launch_angle) %>% 
  tally %>% 
  arrange(desc(n))
  
  
  
# Let's do a quick visualization of the data
df %>% 
  sample_n(5000) %>% 
  ggplot(., aes(launch_speed, launch_angle, col = woba_value)) + geom_point() + 
  scale_colour_gradient() + xlab("Launch Speed") + ylab("Launch Angle") + 
  labs(colour = "wOBA") + theme_bw() + theme(legend.position = "bottom")
  
  
######################

# Step 1. Create the mesh
max.edge <- 3
locs <- cbind(df$launch_speed, df$launch_angle)

mesh <- inla.mesh.2d(loc = locs,
                      max.edge = c(4, 12),
                      cutoff = 3)
                      
# FYI this might take a while depending on your CPU. If you want it to run more quickly
# update your mesh to something like this:
# mesh <- inla.mesh.2d(loc = locs,
#                      max.edge = c(6, 12),
#                      cutoff = 5)
                      

# let's see how many nodes our mesh has
mesh$n
# the one I've vreated here has ~1900. For model testing purposes you should aim for 800-900

# So what is a mesh? Let's plot it out and see
plot(mesh)
# if you'd like you can also add points to the mesh that show the launch speed and launch angle to show what the mesh is actually doing
#points(df$launch_speed, df$launch_angle, pch = 19, col = 2)


# create a convex hull around the points
hpts <- chull(df$launch_speed, df$launch_angle)
hpts <- c(hpts, hpts[1])

poly <- df[hpts,1:2] %>%
  dplyr::rename(x = launch_speed, y = launch_angle)

mesh_locs <- data.frame(X = mesh$loc[,1], Y = mesh$loc[,2])
mesh_locs$intercept <- GEOmap::inpoly(mesh_locs$X, mesh_locs$Y, POK = list(x = poly$x, y = poly$y))
mesh_locs$intercept <- ifelse(mesh_locs$intercept == 0, NA, 1)


# Step 2. Define the weighting factors a_ik (also called the projector matrix).
A <- inla.spde.make.A(mesh = mesh,
                      loc = locs)


# Step 3: Define the SPDE.
# Define spatial model
spde <- inla.spde2.pcmatern(mesh,
                            alpha = 2, ### give it some flexibility # https://groups.google.com/g/r-inla-discussion-group/c/ZhZVu8YPI8I/m/pUUWL1UZBAAJ
                            prior.range = c(50, 0.5), ### P(practic.range < 100) = 0.9 
                            prior.sigma = c(1, 0.1),   ### P(sigma > 5) = 0.1
                            constr = FALSE)

# Step 4. Define the spatial field
s.index <- inla.spde.make.index(
  name = 'spatial.field',
  n.spde = spde$n.spde)

# First: Stack for model inputs
stk.fit <- inla.stack(
  tag = 'fit',
  data = list(Y = df$woba_value, link = 1),
  A = list(1, 1, A), 
  effects = list(intercept = rep(1, nrow(df)),
                 team = df$home_team,
                 s.index))

# Second: Stack for predictions
stk.pred <- inla.stack(
  tag = 'pred',
  data = list(Y = matrix(NA, nrow = ncol(A)), link = 1),
  A = list(1, 1), 
  effects = list(intercept = mesh_locs$intercept,
                 s.index))


stk.all <- inla.stack(stk.fit, stk.pred)

# create the formula for our model. home_team is our random effect and spatial.field if what we used to define
# the spatial field
form <- Y ~ -1 + 
  intercept + f(team, model = "iid") + f(spatial.field, model = spde) 

# run the model. 
mod <- inla(form,
            family = "gaussian",
            data = inla.stack.data(stk.all), 
            control.predictor = list(A = inla.stack.A(stk.all), 
                                     compute = T, link = link),    
            control.compute = list(cpo = TRUE, waic = TRUE))

# summary of the finished model
summary(mod)




# Plot the fitted values from our model
# 

i.x = inla.stack.index(stk.all, tag = "pred")$data
prop = matrix(mod$summary.fitted.values$`0.5quant`[i.x],
              spde$n.spde)
proj = inla.mesh.projector(mesh, dims = c(300, 300))
field.proj <- inla.mesh.project(proj, prop)

loc = expand.grid(proj$x, proj$y)
colnames(loc) = c("x","y")
loc$z = as.vector(field.proj)


# clip the data so that we only plot data inside a convex polygon (plotting close to real world observed values)
hpts <- chull(df$launch_speed, df$launch_angle)
hpts <- c(hpts, hpts[1])

poly <- df[hpts,1:2] %>%
  dplyr::rename(x = launch_speed, y = launch_angle)

mesh_locs <- data.frame(X = mesh$loc[,1], Y = mesh$loc[,2])
xx <- GEOmap::inpoly(loc$x, loc$y, POK = list(x = poly$x, y = poly$y))

# only include points inside the polygon
loc <- loc[xx == 1,]

# finally make the plot
ggplot() + 
  geom_tile(data = loc, aes(x = x, y = y, fill = z)) + 
  xlab("Launch Speed (mph)") + ylab("Launch Angle (degrees)") + 
  scale_fill_gradient2(midpoint = 0.7) + theme_bw() + 
  ggtitle("xwOBAcon Model Prediction") +
  labs(fill = "wOBA") + theme(legend.position = "bottom",
                                 plot.title = element_text(hjust = 0.5))





# Remember how we included the home team as a random effect? Well here is how we can see the results
# We should see Colorado at the top (showing that Coors field favours hitter), and
# Seattle and St. Louis at or near the bottom (showing that those are pitchers parks)

mod$summary.random$team %>% 
  arrange(`0.5quant`) %>% 
  mutate(ID = factor(ID, levels = ID)) %>% 
  ggplot(., aes(`0.5quant`, ID)) + geom_point() + 
  geom_errorbarh(aes(xmin = `0.025quant`, xmax = `0.975quant`)) + 
  ylab("Home Team Park") + xlab("Relative Park Effect") + theme_bw()
  


# print out the values for our random effects

mod$summary.random$team %>% 
  arrange(`0.5quant`) %>% 
  mutate(ID = factor(ID, levels = ID)) %>% 
  rename(Team = ID) %>% 
  dplyr::select(Team, mean, sd, `0.025quant`, `0.975quant`) %>%
  arrange(desc(mean))
