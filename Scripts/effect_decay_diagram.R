rm(list = ls())

library(tidyr)
library(ggplot2)

# start age
start_age  <- 18
end_age    <- 30
mean_METS  <- 3000
perc_habit <- 0.25
exp_decay  <- 0.5


# intervention effect decay diagram
v_time <- seq(start_age, 
              end_age, 
              1)

# create baseline PA decay
v_base <- (1 - pweibull(q = seq(0, 1, length.out = length(v_time)),
                           shape = 0.7, 
                           scale = 1)) * mean_METS
# intervention effect remains forever
v_int <- (1 - pweibull(q = seq(0, 1, length.out = length(v_time)),
                           shape = 0.8, 
                           scale = 1)) * 
  seq(1.1, 
      1.5, 
      length.out = length(v_time)) * mean_METS

# decay after set number of periods (e.g. 2 years)
v_imm <- c(v_int[1:3], v_base[4:length(v_time)])
# with proportion (assume 25% habit)
v_imm_prop <- v_imm * (1 - perc_habit) + v_int * perc_habit

# linear decay
decay_lin_mult <- c(seq(1, 0, length.out = 10), 
                    rep(0, length.out = length(v_time) - 10))
v_lin <- v_int * decay_lin_mult + v_base * (1-decay_lin_mult)
# with proportion (assume 25% habit)
v_lin_prop <- v_lin * (1 - perc_habit) + v_int * perc_habit

# exponential decay
decay_exp_mult <- (1-exp_decay)^(0:(length(v_time)-1))
v_exp <- v_int * decay_exp_mult + v_base * (1-decay_exp_mult)
# with proportion (assume 25% habit)
v_exp_prop <- v_exp * (1 - perc_habit) + v_int * perc_habit

# data.frame of intervention effect decays...
df <- data.frame(v_time, v_base, v_int, 
           v_imm, v_lin, v_exp, 
           v_imm_prop, v_lin_prop, v_exp_prop) 

# convert for plotting function
df_plot <- tidyr::pivot_longer(data = df, 
                               cols = c("v_imm", "v_lin", "v_exp",
                                        "v_imm_prop", "v_lin_prop", "v_exp_prop"))

# create ggplot
plot_of_decays <- ggplot() +
  geom_line(data = df_plot,
            mapping = aes(v_time, col = name, y = value)) +
  geom_line(data = df_plot,
            mapping = aes(v_time, y = v_base),
            col = "black") +
  geom_line(data = df_plot,
            mapping = aes(v_time, y = v_int),
            col = "black") +
  facet_wrap( ~ name,
              ncol = 2,
              labeller = as_labeller(
                c(
                  v_imm = "Fixed Duration, 2 year",
                  v_lin = "Linear, 10 year",
                  v_exp = "Exponential, 50% per year",
                  v_imm_prop = "Fixed Duration, 2 year, 25% develop Habit",
                  v_lin_prop = "Linear, 10 year, 25% develop Habit",
                  v_exp_prop = "Exponential, 50% per year, 25% develop Habit"
                )
              )) +
  theme_minimal() +
  scale_x_continuous(name = "Age", 
                     limits = c(start_age, 30),
                     breaks = seq(start_age, 30, 2)) +
  scale_y_continuous(name = "MET minutes/wk", limits = c(1000, 3500))+
  scale_color_discrete(guide = "none")

ggsave(filename = "Figures/plot_of_decays.png",
       plot = plot_of_decays,
       width = 12, 
       height = 15,
       scale = 0.5)


