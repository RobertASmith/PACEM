# diagram showing baseline trajectories...
rm(list = ls())

library(ggplot2)



# assign distributions of physical activity
v_age_mids <- c(11, 21, 30, 40, 50, 60, 70, 80)
v_ages     <- sample(x = v_age_mids, size = 1000, replace = T)
v_mets     <- 4000 - 40 * v_ages + rnorm(n = length(v_ages), mean = 0, sd = 1000)
v_mets[v_mets < 1] <- 1
v_mets[v_ages == 60] <- v_mets[v_ages == 60] * 1.4
v_mets[v_ages == 40] <- v_mets[v_ages == 40] * 0.9

# physical activity ggplot
df <- data.frame(age = v_ages,
           mets = v_mets)

# linear regression model fit...
lm_fit <- lm(mets ~ age, data = df)
k <- 5
sigma <- sigma(lm_fit)
ab <- coef(lm_fit); a <- ab[1]; b <- ab[2]

# create x and y values
x <- seq(-k*sigma, k*sigma, length.out = 50)
y <- dnorm(x, 0, sigma)/dnorm(0, 0, sigma)*6

# create a segments for each age.
getSegment <- function(mid_point_age, a, b) {
  x0 <- mid_point_age
  y0 <- a + b * x0
  path <- data.frame(x = y + x0, y = x + y0)
  return(path)
}

# create a list of segments to be input into the model...
l_segs <- lapply(X = as.list(v_age_mids),
                 getSegment,
                 a = lm_fit$coefficients['(Intercept)'],
                 b = lm_fit$coefficients['age'])
names(l_segs) <- paste0("age_", v_age_mids)

# make the median trajectory...
df_medians <- df %>% group_by(age) %>% summarise(med_mets = median(mets))
df_medTraj <- as.data.frame(approx(x = df_medians$age, y = df_medians$med_mets, 
                                    xout = 10:80, method = "linear"))
df_medTraj[df_medTraj$x < 18, "y"] <- NA

# create the plot
p2 <- ggplot(df, 
       mapping = aes(x = age, y = mets)) + 
  geom_jitter(color="black", width = 1, size = 1, shape = 17 ) + 
  geom_violin(mapping = aes(x = age, y = mets, group = age), alpha = 0.2)+
  geom_line(data = df_medTraj, 
            mapping = aes(x = x, y = y), 
            size = 1.5)+
  scale_x_continuous(name = "Age in years", limits = c(5, 85)) +
  scale_y_continuous(name = "MET minutes per week", limits = c(0, 7000))+
  theme_minimal()+
annotate("text",
         x = 14,
         y = 6500,
         label = "Intervention") +
  annotate(
    "rect",
    xmin = 10,
    xmax = 18,
    ymin = 0,
    ymax = 7000,
    alpha = .2,
    fill = "grey",
    col = "grey"
  )+
  annotate("text",
           x = 50,
           y = 6500,
           label = "Projection") +
  annotate(
    "rect",
    xmin = 18,
    xmax = 85,
    ymin = 0,
    ymax = 7000,
    alpha = .2,
    fill = "lightblue",
    col = "lightblue"
  )

ggsave(filename = "Figures/Final/trajectoryPA.png",
       plot = p2, 
       width = 15, 
       height = 12,
       scale = 0.5)


#for(a in names(l_segs)){
#  
#  p2 <- p2 + geom_path(data = l_segs[[a]], 
#                       mapping = aes(x,y), 
#                       color = "black", 
#                       size = 1, 
#                       alpha = 0.5)
#
#  }
#
#p2
#
#
#df <- data.frame(age = c(11, 21, 30, 40, 50, 60, 70, 80),
#                 METmins = c(3000, 2000, 1700, 1500, 1400, 1600, 1200, 1000))
#
#df_METsmooth = as.data.frame(approx(x = df$age, y = df$METmins, 
#                                    xout = 10:80, method = "linear"))
#
#
#p1 <- ggplot() +
#  geom_line(data = df_METsmooth,
#            mapping = aes(x = x, y = y)) +
#  geom_point(data = df,
#             mapping = aes(x = age, y = METmins), size = 3) +
#  scale_x_continuous(name = "Age in years", limits = c(10, 80)) +
#  scale_y_continuous(name = "MET minutes per week", limits = c(0, 6000)) +
#  theme_minimal() +
#  annotate("text",
#           x = 14,
#           y = 1000,
#           label = "Intervention \n Period") +
#  annotate(
#    "rect",
#    xmin = 10,
#    xmax = 18,
#    ymin = 0,
#    ymax = 4000,
#    alpha = .2,
#    fill = "grey",
#    col = "grey"
#  )
## add in the blue sections showing the periods for which PA is estimated.
#v_ages <- c(18, 25, 35, 45, 55, 65, 75, 80)
#for(x in 1:length(v_ages)) {
#  p1 <- p1 + annotate(
#    "rect",
#    xmin = v_ages[x],
#    xmax = v_ages[x+1],
#    ymin = 0,
#    ymax = 4000,
#    alpha = .2,
#    fill = x+1,
#    col = "grey40"
#  )
#}
#  
#  
#ggsave(filename = "Figures/Final/trajectoryPA.png", 
#       plot = p1,
#       width = 15,
#       height = 12)
#
#
#
#
#
#
#library(ggplot2)
#
#x <- runif(n = 100, min =  0,max =  15)
#y <- 1000 + 200*x + rnorm(100, 0, 300)
#df <- data.frame(x, y)
#lm_fit <- lm(y ~ x, data = df)
#
#k <- 2.5
#sigma <- sigma(lm_fit)
#ab <- coef(lm_fit); a <- ab[1]; b <- ab[2]
#
#x <- seq(-k*sigma, k*sigma, length.out = 50)
#y <- dnorm(x, 0, sigma)/dnorm(0, 0, sigma) * 3
#
#x0 <- 0
#y0 <- a+b*x0
#path1 <- data.frame(x = y + x0, y = x + y0)
#segment1 <- data.frame(x = x0, y = y0 - k*sigma, xend = x0, yend = y0 + k*sigma)
#
#x0 <- 5
#y0 <- a+b*x0
#path2 <- data.frame(x = y + x0, y = x + y0)
#segment2 <- data.frame(x = x0, y = y0 - k*sigma, xend = x0, yend = y0 + k*sigma)
#
#x0 <- 10
#y0 <- a+b*x0
#path3 <- data.frame(x = y + x0, y = x + y0)
#segment3 <- data.frame(x = x0, y = y0 - k*sigma, xend = x0, yend = y0 + k*sigma)
#
#ggplot(df, 
#       mapping = aes(x=x, y=y)) + 
#  geom_point(color="blue") + 
#  geom_smooth(method='lm', se=FALSE, color="red") + 
#  geom_path(aes(x,y), data = path1, color = "green") + 
#  geom_segment(aes(x=x,y=y,xend=xend,yend=yend), data = segment1) +
#  geom_path(aes(x,y), data = path2, color = "green") + 
#  geom_segment(aes(x=x,y=y,xend=xend,yend=yend), data = segment2) +
#  geom_path(aes(x,y), data = path3, color = "green") + 
#  geom_segment(aes(x=x,y=y,xend=xend,yend=yend), data = segment3)  
#  #