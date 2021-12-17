require(ggplot2)

curveDf <- data.frame("batchSize" = c(16,8,4,2,1), "s" = c(1,1.7,3,5.3,7.8))
plot(curveDf$s ~ rev(curveDf$batchSize), type = "b")
lines(rev(curveDf$batchSize), 1/curveDf$s)
plotDf <- data.frame("batchSize" = rev(curveDf$batchSize), "speedUp" = curveDf$s, "relativeEfficiency"= 1/curveDf$s)

temperatureColor <- "#69b3a2"
priceColor <- rgb(0.2, 0.6, 0.9, 1)
pdf("./imgResults/fig1_efficiency.pdf", width = 5,height = 3)
ggplot(plotDf, aes(x=batchSize)) +
    geom_line( aes(y=relativeEfficiency*7, color = temperatureColor)) + 
    geom_line( aes(y=speedUp)) + 
    geom_point( aes(y=relativeEfficiency*7, color = temperatureColor)) + 
    geom_point( aes(y=speedUp)) + 
    scale_y_continuous(
        name = "Speed up",
        sec.axis = sec_axis(~./7, name="Relative Efficiency"),
        breaks = c(0,2,4,6,8,10)) + 
    scale_x_continuous(name = "Number of Nodes",
                       breaks = c(0,2,4,6,8,10,12,14,16)) +
    theme(
        axis.title.y = element_text(color = 1, size=13),
        axis.title.y.right = element_text(color = 2, size=13, angle = 90),
        legend.position = "none"
    )
dev.off()
