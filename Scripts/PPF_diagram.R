#devtools::install_github("R-CoderDotCom/econocharts")

library(econocharts)
library(ggplot2)

# Figure 1, PPF tradeoff 
p <- ppf(
  x = c(5, 4),
  # Intersections
  geom = "text",
  # Intersection labels as text
  generic = TRUE,
  # Generic axis tick labels
  labels = c("A", "B"),
  # Custom labels
  xlab = "External Validity",
  # X-axis label
  ylab = "Usability"
) 

final_plot <- p$p+
  geom_point(data = data.frame(x = 5, y = 5), 
             size = 3) +
  annotate(geom = "text",
           x = 5.25,
           y = 5.25,
           label = "C") +
  annotate(
    "segment",
    x = 3.25,
    xend = 4.25,
    y = 5,
    yend = 5,
    arrow = arrow(length = unit(0.5, "lines")),
    colour = 1,
    lwd = 1
  ) +
  annotate(
    "segment",
    x = 4.25,
    xend = 4.25,
    y = 5,
    yend = 4.25,
    arrow = arrow(length = unit(0.5, "lines")),
    colour = 1,
    lwd = 1
  )


ggsave(filename = "Figures/PPF.png",
       plot = final_plot,
       height = 7, width = 7,
       scale = 0.8)
