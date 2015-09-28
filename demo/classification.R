
xp <- NULL
xp <- rbind(xp, c( 1,  5))
xp <- rbind(xp, c( 5,  4))
xp <- rbind(xp, c(-1, -1))
xp <- rbind(xp, c( 1,  0))
xp <- rbind(xp, c(-5,  1))
xp <- rbind(xp, c(-1, -5))
yp <- NULL
yp <- rbind(yp, c( 1,  0))
yp <- rbind(yp, c( 1,  0))
yp <- rbind(yp, c( 1,  0))
yp <- rbind(yp, c( 0,  1))
yp <- rbind(yp, c( 0,  1))
yp <- rbind(yp, c( 0,  1))

x  <- as.matrix(expand.grid(x = seq(-5, 5, length.out = 100),
                            y = seq(-5, 5, length.out = 100)))

gp <- new.gp(0.5, kernel.squared.exponential(8, 10),
             likelihood=NULL,
             link=new.link("probit"),
             dim=2)
gp <- posterior(gp, xp, yp)

p <- plot(gp, x, verbose=T)
p$p1 <- p$p1 + scale_shape_identity()
p$p1 <- p$p1 + geom_point(data=data.frame(x=xp[,1], y=xp[,2],
                                    z = sapply(yp[,1], function(x) if(x == 1) 43 else 45)),
                    aes(x = x, y = y, shape = z, size = 100)) +
    scale_size_area()
grid.arrange(p$p1, p$p2, ncol=2)
