# Project "Effect of allelic dominance on stable polymorphism"
# Hannah Heidelmeyer, Sofie Rutenbeck, Benjamin Winkler
# 22.08.24 - 16.09.24
# Start ----

# loading packages:
library(tidyverse)                # data formating
library(ggplot2)                  # data plotting
library(patchwork)                # arrangement of plots
library(gridExtra)                # arrangement of plots
library(stickylabeller)           # labelling inside facet_wrap() command

# options(scipen = 666)   # turn scientific notation off (optionally)


# defining the base parameters: ----
Sf <- 0.1      # Sf = fitness effect of F allele in females
Sm <- 0        # Sm = fitness effect of M allele in males                      
Hf <- 1        # Hf = degree of allelic dominance of F allele         
Hm <- 0        # Hm = degree of allelic dominance of M allele          
p <- 0.01      # p = frequency of F allele       
q <- (1 - p)   # q = frequency of M allele          


# creating contingency square of fitness effects per sex in a single-locus model: ----
## using the parameters above, excluding frequencies: 
square <- as.data.frame(rbind(c((1 + Sf), (1 + Sf * Hf), 1),
                              c(1, (1 + Sm * Hm), (1 + Sm))))  # 3x2 df

rownames(square) <- c("Female", "Male")    # renaming rows

colnames(square) <- c("FF", "FM", "MM")    # renaming columns 

square # shows allelic fitness effects in each sex for a specific set of parameters


# equation for a single generation change, no iterations: ----
## equations based on the Hardy-Weinberg equilibrium:
## for further explanation, see Methods:
P_f_1 <- ((p^2 * (1 + Sf) + (2 * p * q * (1 + Sf * Hf)) * 0.5) / 
            (p^2 * (1 + Sf) + 2 * p *q *(1 + Sf * Hf) + q^2))

P_m_1 <- ((p^2 + (2 * p * q *(1 + Sm * Hm)) * 0.5) / 
            (p^2 + (2 * p * q * (1 + Sm * Hm)) + q^2 * (1 + Sm)))

P_both_1 <- (P_f_1 + P_m_1) / 2  # _1 as for first offspring generation


# iterating equation; calculating output for one specific generation: ----
population <- tibble(p = 0.01)   # starting at 0.01% of the population having the F allele
current_gen <- 1                 # starting in 1st generation

while(current_gen <= 10) {           # example for calculation after 10th generation
  population <- population %>% 
    mutate(q = 1 - p,
           P_f = ((p^2 * (1 + Sf)) + ((2 * p * q * (1 + Sf * Hf)) * (0.5))) / 
                    ((p^2 * (1 + Sf)) + ((2 * p *q *(1 + Sf * Hf)) + (q^2))),
           P_m = ((p^2) + ((2 * p * q *(1 + Sm * Hm))) * (0.5)) / 
                    ((p^2) + ((2 * p * q) * (1 + Sm * Hm)) + (q^2 * (1 + Sm))),
           p = ((P_f + P_m) / 2))
  
  current_gen <- (current_gen + 1)  # iterating function unless current_gen becomes <= 10
  
  population  # data for values stored in population df
}
                  

# creating an iterating function to track per generation change: ----
trajectory_function <- function(Sf, Sm, Hf, Hm, generations){
  p <- 0.01
  
  current_gen <- 1
  
  pop <- tibble(p = 0.01,
                current_gen = 1,
                q = 1 - p)
  
  while(current_gen <= generations & 
        p >= 0.01 & p <= 0.99) {        # stops the calculation after hitting p tresholds
    q = (1 - p)
    
    p_freq_f <- ((p^2 * (1 + Sf)) + ((2 * p * q * (1 + Sf * Hf)) * (0.5))) / 
      ((p^2 * (1 + Sf)) + ((2 * p * q *(1 + Sf * Hf)) + (q^2)))
    
    p_freq_m <- ((p^2) + ((2 * p * q *(1 + Sm * Hm))) * (0.5)) / 
      ((p^2) + ((2 * p * q) * (1 + Sm * Hm)) + (q^2 * (1 + Sm))) 
    
    p <- ((p_freq_f + p_freq_m) / 2)
    
    pop <- pop %>% 
      add_row(p = p,
              current_gen = current_gen,
              q = q)
    
    current_gen <- (current_gen + 1) 
  }
pop
}


# calculation of the trajectory of frequencies after the invasion of the F allele, with specific parameters:
## high fitness effect for females, no fitness effect for males, same dominance for both:
inv_sf_only <- trajectory_function(Sf = 0.5, Sm = 0.0, Hf = 0.5, Hm = 0.5, generations = 1000)

## similar high fitness effect for both, same dominance for both:
inv_both <- trajectory_function(Sf = 0.55, Sm = 0.45, Hf = 0.5, Hm = 0.5, generations = 1000)


# plotting the trajectories: ----
## Figure 1, a):
(line_plot_1 <- ggplot(data = inv_sf_only,           # dataset used
                       aes(x = current_gen,          # x-axis
                           y = p)) +                 # y-axis
    geom_line(alpha = 2,                             # type of plot; opaqueness of line
              colour = "black",                      # colour of lines
              size = 1.5) +                          # thickness of lines
    geom_line(aes(x = current_gen,                   # new line plot added
                  y = q),
              colour = "black",
              alpha = 2,
              size = 1.5,
              linetype = 2) +                        # changing linetype to dashed
    scale_x_continuous(expand = c(0, 0),                        # starting plot at x = 0
                       breaks = c(0, 20, 41, 60, 80, 100)) +    # labelling of x-axis; breaks
    scale_y_continuous(expand = c(0, 0),
                       limits = c(0, 1)) +           # not defining breaks but limits on y-axis
    geom_vline(xintercept = 41,                      # adding vertical line at x = 41
               linetype = "dotted",                  # dotted line
               colour = "red",
               size = 1.5) + 
    geom_label(label = "Frequency F allele",         # adding label to differentiate both lines plotted
               x = 70,                               # label at x = 70
               y = 0.80,                             # + y = 0.80
               colour = "black",                     # colour of text inside label
               fill = "white",                       # colour of box surrounding label
               size = 5) +                           # total size of label
    geom_label(label = "Frequency M allele",
               x = 70,
               y = 0.20,
               colour = "black",
               fill = "white",
               size = 5) + 
    ggtitle("a) Sf = 0.5, Sm = 0.0") +               # adding title to plot
    xlab("Generations") +                            # adding x-axis label
    ylab("Frequency") +                              # adding y-axis label
    theme_bw() +                                     # adding white background
    theme(text = element_text(size = 20)))           # increasing all text sizes in plot

## Figure 1, b):
(line_plot_2 <- ggplot(data = inv_both,           
                       aes(x = current_gen,
                           y = p)) +
    geom_line(alpha = 2, 
              colour = "black", 
              size = 1.5) +
    geom_line(aes(x = current_gen,
                  y = q),
              colour = "black",
              alpha = 2,
              size = 1.5,
              linetype = 2) +
    scale_x_continuous(expand = c(0, 0),
                       breaks = c(0, 113, 250, 500, 750, 1000),
                       limits = c(0, 1030)) +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(0, 1)) +
    geom_vline(xintercept = 113,
               linetype = "dotted",
               colour = "red",
               size = 1.5) + 
    geom_label(label = "Frequency F allele",
               x = 800,
               y = 0.80,
               colour = "black",
               fill = "white",
               size = 5) +
    geom_label(label = "Frequency M allele",
               x = 800,
               y = 0.20,
               colour = "black",
               fill = "white",
               size = 5) + 
    ggtitle("b) Sf = 0.55, Sm = 0.45") +
    xlab("") +
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = 20)))


grid.arrange(line_plot_1, line_plot_2,  # defining the 2 plots that should be arranged
             ncol = 2,                  # plotted in 2 columns next to each other
             nrow = 1)                  # arranged in 1 row


# iterating calculation for a combination of multiple parameters: ----
param <- expand_grid(Sf = seq(0, 1, length = 50),   # 50 values of Sf between 0 and 1
                     Sm = seq(0, 1, length = 50),   # 50 values of Sm between 0 and 1
                     Hm = c(0, 0.5, 1)) %>%         # 3 values of Hm
  mutate(Hf = (1 - Hm)) %>%                         # 3 values of H
  select("Sf", "Sm", "Hf", "Hm")                    # reordering columns

## 50 * 50 * 3 = 7500 combinations of values in total

run_sim <- function(generations, row, parameters){
  
library(tidyverse)

  p <- 0.01
  
  population = tibble(p = 0.01,
                      q = 1 - p)
  
  current_gen = 1
  
  Sf = parameters$Sf[row]
  
  Sm = parameters$Sm[row]
  
  Hf = parameters$Hf[row]
  
  Hm = parameters$Hm[row]
  
  while(current_gen <= generations & 
        p >= 0.01 & p <= 0.99){
    population <- population %>% 
      mutate(q = 1 - p,
             P_f = ((p^2 * (1 + Sf)) + ((2 * p * q * (1 + Sf * Hf)) * (0.5))) / 
               ((p^2 * (1 + Sf)) + ((2 * p *q *(1 + Sf * Hf)) + (q^2))),
             P_m = ((p^2) + ((2 * p * q *(1 + Sm * Hm))) * (0.5)) / 
               ((p^2) + ((2 * p * q) * (1 + Sm * Hm)) + (q^2 * (1 + Sm))),
             p = ((P_f + P_m) / 2),
             gen= current_gen)
    
    current_gen <- (current_gen + 1)
    
    p <- population$p
  }
  bind_cols(parameters[row, ], population)
}

## calculation for a specific row in param df:
# test <- run_sim(row = 151, 1000, param) # for 151st row

## calculations for the first 10 rows of param:
# simulation_output <- map_dfr(1:10, run_sim,
#                              generations = 1000,
#                              parameters = param)  


# using parallel cores to speed up the run_sim function: ----
n_cores <- parallel::detectCores()            # find cores

cluster <- parallel::makeCluster(n_cores)     # create cluster

parallel::clusterExport(cluster, ("param"))   # split parameters across cluster

output <- parallel::parLapply(cluster,
                              1:7500,              # all 7500 rows
                              run_sim,             # using run_sim function
                              generations = 5000,  # run for 1000/ 5000 generations
                              parameters = param)  # using param df

df_param <- do.call(rbind.data.frame, output) # formatting the list into a df

df_param_5000 <- do.call(rbind.data.frame, output) # when run with 5000 generations

## saving dfs:
# save.image("C:/Uni Mainz/7. Semester/Theor_Evo&Eco/Project_Stable_Polymorphism/df_param.RData")
# save.image("C:/Uni Mainz/7. Semester/Theor_Evo&Eco/Project_Stable_Polymorphism/df_param_5000.RData")

## importing dfs:
# load("C:/Uni Mainz/7. Semester/Theor_Evo&Eco/Project_Stable_Polymorphism/df_param.RData")
# load("C:/Uni Mainz/7. Semester/Theor_Evo&Eco/Project_Stable_Polymorphism/df_param_5000.RData")


# creating new columns in dfs to highlight specific tiles or regions in heatmap: ----
## df with 1000 generations: 
df_param_types <- df_param %>%                     # using df_param to create new df
  mutate(case_when("State of alleles" =            # adding a new column
                   p > 0.99 ~ "Only F allele",   # specify threshold and define name of values within
                   p < 0.01 ~ "Only M allele",
                   p < 0.99 & p > 0.55 ~ "Putative stable polymorphism",
                   p < 0.6 & p > 0.4 ~ "Putative equalised polymorphism",
                   p < 0.45 & p > 0.01 ~ "Putative stable polymorphism"))

colnames(df_param_types)[10] <- "state"  # renaming "State of alleles" column for simplicity

table(df_param_types$state, useNA = "always") # show counts of states; check for NAs

df_param_types <- df_param_types %>% 
  mutate(row_ID = as.factor(1:n()))  # necessary for finding and highlighting a single tiles

new_order <- c("Only M allele",   # defining new order of state levels; just for aesthetical reasons later 
               "Putative stable polymorphism",
               "Putative equalised polymorphism",
               "Only F allele")                      

df_param_types$state <- factor(df_param_types$state,
                               levels = new_order) # applying new order of state levels


## df with 5000 generations:
df_param_5000_types <- df_param_5000 %>%    # same things as above just with changed df
  mutate(case_when("State of alleles" =
                     p > 0.99 ~ "Only F allele",
                   p < 0.01 ~ "Only M allele",
                   p < 0.99 & p > 0.55 ~ "Putative stable polymorphism",
                   p < 0.6 & p > 0.4 ~ "Putative equalised polymorphism",
                   p < 0.45 & p > 0.01 ~ "Putative stable polymorphism"))

colnames(df_param_5000_types)[10] <- "state"

df_param_5000_types <- df_param_5000_types %>% 
  mutate(row_ID = as.factor(1:n()))  

new_order <- c("Only M allele",
               "Putative stable polymorphism",
               "Putative equalised polymorphism",
               "Only F allele")                      

df_param_5000_types$state <- factor(df_param_5000_types$state, levels = new_order) 


# plotting the outcome of all combinations as a heatmap: ----
## Figure 2:
(heat_plot <- ggplot(data = na.omit(df_param_types), # na.omit to remove 3 NAs where Sf = Sm = 0
                     aes(x = Sf,
                         y = Sm,
                         fill = state)) +              # colouring tiles per state
    scale_fill_manual(name = "Frequency of F allele:",   # defining header of label
                      labels = c("< 0.01",               # defining texts of label
                                 "> 0.01 & < 0.40 or\n> 0.60 & < 0.99",
                                 "> 0.40 & < 0.60",
                                 "> 0.99"),
                      values = c("bisque4", "darkseagreen1", "coral2", "darkgrey")) + # defining colour of states
    geom_tile(colour = "black") +       # defining heatmap as plot
    facet_wrap(~ Hf + Hm,               # splitting plots per Hf/ Hm values                  
               labeller = label_glue("hf = {Hf}, hm = {Hm}")) +  # defining headers of split plots
    geom_rect(data = filter(df_param_types,            # highlighting one specific tile
                           row_ID == "3752"),          # row_ID correlates to fig. 1, a)
             size=1, fill=NA, colour="cyan",           # defining aesthetics of this tile
             aes(xmin=Sf - 0.01, xmax=Sf + 0.01, ymin=Sm - 0.01, ymax=Sm + 0.01)) + # size of one tile
    geom_rect(data = filter(df_param_types, 
                           row_ID == "4121"),          # row_ID correlates to fig. 1, b)
             size=1, fill=NA, colour="magenta",        
             aes(xmin=Sf - 0.01, xmax=Sf + 0.01, ymin=Sm - 0.01, ymax=Sm + 0.01)) +
    scale_x_continuous(expand = c(0, 0)) +             # removes empty space caused by facet_wrap() on x-axis
    scale_y_continuous(expand = c(0, 0)) +             # same as above for y-axis 
    coord_fixed() +                                    # forces plots to be squares with same length/ height
    xlab("Fitness effect of F allele in females") + 
    ylab("Fitness effect of F allele in males") +
    theme_bw() +
    theme(text = element_text(size = 30)) +         # increase general text size
    theme(panel.spacing = unit(4, "lines")) +       # add space between facet_wrap() panels
    theme(legend.title = element_text(size = 15),   # decreases legend header size
          legend.text = element_text(size = 15)))   # decreases legend text size


## Figure 3:
(heat_plot_2 <- ggplot(data = na.omit(df_param_5000_types), 
                     aes(x = Sf,
                         y = Sm,
                         fill = state)) +
    scale_fill_manual(name = "Frequency of F allele:",
                      labels = c("< 0.01",
                                 "> 0.01 & < 0.40 or\n> 0.60 & < 0.99",
                                 "> 0.40 & < 0.60",
                                 "> 0.99"),
                      values = c("bisque4", "darkseagreen1", "coral2", "darkgrey")) +
    geom_tile(colour = "black") +
    facet_wrap(~ Hf + Hm,                                
               labeller = label_glue("hf = {Hf}, hm = {Hm}")) + 
    scale_x_continuous(expand = c(0, 0)) +             
    scale_y_continuous(expand = c(0, 0)) +             
    coord_fixed() +   
    xlab("Fitness effect of F allele in females") +
    ylab("Fitness effect of F allele in males") +
    theme_bw() +
    theme(text = element_text(size = 30)) +             
    theme(panel.spacing = unit(4, "lines")) +          
    theme(legend.title = element_text(size = 15), 
          legend.text = element_text(size = 15)))        
