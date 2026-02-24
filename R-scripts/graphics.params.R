# graphics parameters ----

## Palettes to use ----
tscolor <- "red4"
aggtscolor <-"deepskyblue4"
aggtscolor2 <- "deepskyblue"
axistitle <- 20
textsize <- 16
linesize <- 3
nPanel <- 12

# Graph of variables
PCA_color_gradient <- c("#36648B", "#FFA500", "#8B2500")
PCA_color_gradient <- brewer.pal(9, name="YlOrBr")[c(3,4,5,6,7,8,9)]

# for ETI
colors <- data.frame(
  val = 1:11,
  col = brewer.pal(11,'RdYlGn'),
  labels = c("Collapse",
             "Shocks Likely",
             "Shocks Possible",
             "Low Integrity",
             "Med Integrity",
             "Med-High Integrity",
             "High Integrity",
             "Robust Integrity",
             "Close to Pristine",
             "Pristine",
             "Pristine"))
#display.brewer.pal(10, name='RdYlGn')
