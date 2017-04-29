install.packages("BayesFactor", repos = "https://cloud.r-project.org/", type="source", dependencies = T)
library(BayesFactor)
options(scipen = 1)


N = c(20, 40, 80) # sample size, index with n
ES = c(0, .2, .5) # effect size, index with e (is mean, sd is always 1)
Trails = c(2, 3, 5) # number of Trails, index with t 
D = c(-.2, 0, .2) # inferiority-margin delta, index with d
Runs = 500 # index with r
IntervalPoint = c(1, 2) # 1 for interval null 2 for point null

LowHighBF01 = array(NA, dim =c(2, length(IntervalPoint), length(N), length(ES), length(Trails), length(D), Runs))
# first dimension = lowest highest BFs, lowest is 1 


#Optim Point func------
RPointVary <- function(data, par){# No need for D[d] for point BFs because you only compute them in D[d] = 0 cases
  BF <- 1/as.vector(ttestBF(x = data[,1], y= data[,2], nullInterval = c(0, Inf), rscale = par))[1]
  BF
}

## Simulation 1 -----

for(d in 3:length(D)){
  #Optim Interval func ------
  RInterVary <- function(data, par){ # CORRECTION: should have been: BF <- ttestBF(x = (data[,1]-D[d]), y = (data[,2]), nullInterval = c(0, Inf), rscale = par)
    BF <- ttestBF(x = (data[,1]-D[d]), y = (data[,2]-D[d]), nullInterval = c(0, Inf), rscale = par)
    as.vector(BF[2]/BF[1])[1]
  }
  #----
  for(e in 1:length(ES)){
    for(t in 1:length(Trails)){
      for(n in 1:length(N)){
        cat("N = ", N[n], ", ", "Trails = ", Trails[t], ", ", "ES = ", ES[e], ", ", "D = ", D[d],  "\n ")
        for(r in 1:Runs){
          p = rep(0.5, Trails[t])
          while(!(sum(p < .025) == 2 & sum(p > .975) == 0)){
            Control = matrix(rnorm(N[n]*Trails[t], 0, 1), nrow = Trails[t], ncol = N[n])
            Treat = matrix(rnorm(N[n]*Trails[t], ES[e], 1), nrow = Trails[t], ncol = N[n])
            for(i in 1:Trails[t]){
              p[i] = t.test (Treat[i,], Control[i,], mu = D[d], alternative = "greater", var.eq=TRUE)$p.value
            }
          }
          ### optim()----------------
          Data= matrix(c(Treat, Control), ncol=2)
          OptimMin <- optim(par = sqrt(2)/2, fn =RInterVary, data= Data, method="Brent",
                            lower=sqrt(2)/2, upper=sqrt(2))
          OptimMax <- optim(par = 0, fn = RInterVary, data= Data,  method="Brent",
                            lower=sqrt(2)/2, upper=sqrt(2), control=list(fnscale= -1))
          LowHighBF01[ , IntervalPoint[1], n, e, t, d, r] <- c(OptimMin$value, OptimMax$value)
          if(!D[d]){# Point Nulls ----
            OptimMinP <- optim(par = sqrt(2)/2, fn = RPointVary, data= Data, method="Brent", 
                               lower=sqrt(2)/2, upper=sqrt(2))
            OptimMaxP <- optim(par = sqrt(2)/2, fn = RPointVary, data= Data,  method="Brent",
                               lower=sqrt(2)/2, upper=sqrt(2), control=list(fnscale= -1))
            LowHighBF01[ , IntervalPoint[2], n, e, t, d, r] <- c(OptimMinP$value, OptimMaxP$value)
            #----
          }
          #-----
        }
        save(LowHighBF01, file = paste("/home/anja/Desktop/ReMa/FDA/C_R/data/","Sim1BF", ".RData", sep=""))
      }
    }
  }
} #-------


## Simulation 2 ----- 
N = c(20, 40, 80) # sample size, index with n
ES = c(0, .2, .3) # effect size, index with e (is mean, sd is always 1)
Trails = 2 # number of Trails, not varied
D <- matrix(c(-.2, -.3, .2, .3), 2) # inferiority-margin delta, not varied
Runs = 500 # index with r


EquivBF01 = array(NA, dim =c(2, length(D[,1]), length(N), length(ES), Runs))
# first element = lowest highest, lowest is 1 



# Sim 2----
for(d in 1:length(D[,1])){
  # Equivalence Optim Interval func ------
  REquivalenceVary <- function(data, par){ 
    BF <- ttestBF(x = data[,1], y= data[,2], nullInterval = c(D[d, 1], D[d, 2]), rscale = par)
    as.vector(BF[1]/BF[2])[1]#BF01: H0 is the equivalence region youdo not want to reject
  }
  for(e in 1:length(ES)){
    for(n in 1:length(N)){
      cat("D = ", D[d, 2], "N = ", N[n], ", ", "ES = ", ES[e],  "\n ")
      for(r in 1:Runs){
        p = matrix(rep(0, 2*Trails), ncol = 2)
        while(!(sum(p < .025) == 0 & sum(p > .975) == 0)){ 
          Control = matrix(rnorm(N[n]*Trails, 0, 1), nrow = Trails, ncol = N[n])
          Treat = matrix(rnorm(N[n]*Trails, ES[e], 1), nrow = Trails, ncol = N[n])
          for(i in 1:Trails){
            p[i, 1] = t.test (Treat[i,], Control[i,], mu = D[d, 1], alternative = "less", var.eq=TRUE)$p.value
            p[i, 2] = t.test (Treat[i,], Control[i,], mu = D[d, 2], alternative = "greater", var.eq=TRUE)$p.value
          }
        }
        ### optim()----------------
        Data = matrix(c(Treat, Control), ncol = 2)
        OptimMin <- optim(par = sqrt(2)/2, fn = REquivalenceVary, data= Data, method="Brent",
                          lower = sqrt(2)/2, upper = sqrt(2))
        OptimMax <- optim(par = 0, fn = REquivalenceVary, data= Data,  method="Brent",
                          lower = sqrt(2)/2, upper = sqrt(2), control = list(fnscale= -1))
        EquivBF01[ , d, n, e, r] <- c(OptimMin$value, OptimMax$value)
      }#-----
    }
  }
}
save(EquivBF01, file = paste("/home/anja/Desktop/ReMa/FDA/C_R/data/","EquivBF01", ".RData", sep=""))


## Result PlotSimulation 1------------------
N = c(20, 40, 80) # sample size, index with n
ES = c(0, .2, .5) # effect size, index with e (is mean, sd is always 1)
Trails = c(2, 3, 5) # number of Trails, index with t 
D = c(-.2, 0, .2) # inferiority-margin delta, index with d
Runs = 500 # index with r
IntervalPoint = c(1, 2)  # 1 for interval null 2 for point null

options(scipen = 1)

Ylim <- c(-12, 2)
Cex1 <- 1.3
Cex2 <- .8 #axis ticks
Cex3  <- 1.3 #ES heading


for(d in 1:length(D)){
  png(filename = paste("/home/anja/Desktop/ReMa/FDA/C_R/plots/","boxplot",d,".png",sep=""), 
      width = 960, height = 960)
  layout(matrix (1:9, 3, 3, byrow = T))
  par(mar = c(4.5, 6.5, 8, 0.5), cex.main = 2, cex.axis = 2)
  for(e in 1:length(ES)){
    for(t in 1:length(Trails)){
      if(!D[d]){
        boxplot(col = adjustcolor("white", alpha=0.2), border= "red", at = c(2, 4, 6), x = t(log(LowHighBF01[1, 1, , e, t, d, ], base = 10)), axes = F, ylim = Ylim, xlim= c(0,7), xlab = "", ylab = "", main = paste ("Trials:", Trails[t]), pch = 16)
        boxplot(col = adjustcolor("red", alpha=0.2), add= TRUE, at = c(1, 3, 5), x = t(log(LowHighBF01[1, 2, , e, t, d, ], base = 10)), axes = F, ylim = Ylim, xlab = "", ylab = "", main = paste ("Trials:", Trails[t]), pch = 16)
        abline(h = 0, lty = 2, lwd = 2)
        axis(1, at = c(0, 1.7 , 3.7, 5.7, 7), labels = c( "", N, ""))
        axis(2, at = Ylim[1]:Ylim[2], labels = 10^(Ylim[1]:Ylim[2]), las = 1, cex= Cex2)
        mtext(text = "N", side = 1, line = 3, cex =Cex1)
        boxplot(add= TRUE, col = adjustcolor("white", alpha=0.2), border="blue", at = c(2.2, 4.2, 6.2), x = t(log(LowHighBF01[2, 1, , e, t, d, ], base = 10)), axes = F, ylim = Ylim, xlab = "", ylab = "", main = paste ("Trials:", Trails[t]), pch = 16)
        boxplot(add= TRUE, col=adjustcolor("blue", .2), at = c(1.2, 3.2, 5.2), x = t(log(LowHighBF01[2, 2, , e, t, d, ], base = 10)), axes = F, ylim = Ylim, xlab = "", ylab = "", main = paste ("Trials:", Trails[t]), pch = 16)
      }else{
        boxplot(col = adjustcolor("white", alpha=0.2), border= "red",  at = c(1, 2, 3), x = t(log(LowHighBF01[1, 1, , e, t, d, ], base = 10)), axes = F, ylim = Ylim, xlim= c(0,4), xlab = "", ylab = "", main = paste ("Trials:", Trails[t]), pch = 16)
        abline(h = 0, lty = 2, lwd = 2)
        axis(1, at = c(0, 1.1:3.1, 4), labels = c( "", N, ""))
        axis(2, at = Ylim[1]:Ylim[2], labels = 10^(Ylim[1]:Ylim[2]), las = 1)
        mtext(text = "N", side = 1, line = 3, cex =Cex1)
        boxplot(add= TRUE, col=adjustcolor("white", .2), border= "blue", at = c(1.2, 2.2, 3.2), x = t(log(LowHighBF01[2, 1, , e, t, d, ], base = 10)), axes = F, ylim = Ylim, xlim= c(0,4), xlab = "", ylab = "", main = paste ("Trials:", Trails[t]), pch = 16)
      }
      if(t%%2 == 0){
        mtext(text = paste("ES =", ES[e]), side = 3, line = 6, cex = Cex3, font = 2)
      }
      if(t*e == 1){
        mtext(text = paste("Non-Inferiority:", D[d]), side = 3, line = 5.4, cex = 2.2, font = 2)
      }
    }
  }
  dev.off ()
}


# Histogram -------------------
# Add the box.cex parameter to the legend function by altering legend source code
# source code from: http://code.metager.de/source/xref/gnu/R/src/library/graphics/R/legend.R
Legend <- function( #--------------------
                    x, y = NULL, legend, fill = NULL, col = par("col"), border="black",
                    lty, lwd, pch, angle = 45, density = NULL, bty = "o", bg = par("bg"),
                    box.lwd = par("lwd"), box.lty = par("lty"), box.col = par("fg"),
                    pt.bg = NA, cex = 1, pt.cex = cex, pt.lwd = lwd,
                    xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1, adj = c(0, 0.5),
                    text.width = NULL, text.col = par("col"), text.font = NULL,
                    merge = do.lines && has.pch, trace = FALSE,
                    plot = TRUE, ncol = 1, horiz = FALSE, title = NULL,
                    inset = 0, xpd, title.col = text.col, title.adj = 0.5,
                    seg.len = 2, box.cex =c(0.8, 0.5))
{
  ## the 2nd arg may really be `legend'
  if(missing(legend) && !missing(y) &&
     (is.character(y) || is.expression(y))) {
    legend <- y
    y <- NULL
  }
  mfill <- !missing(fill) || !missing(density)
  
  if(!missing(xpd)) {
    op <- par("xpd")
    on.exit(par(xpd=op))
    par(xpd=xpd)
  }
  title <- as.graphicsAnnot(title)
  if(length(title) > 1) stop("invalid 'title'")
  legend <- as.graphicsAnnot(legend)
  n.leg <- if(is.call(legend)) 1 else length(legend)
  if(n.leg == 0) stop("'legend' is of length 0")
  auto <-
    if (is.character(x))
      match.arg(x, c("bottomright", "bottom", "bottomleft", "left",
                     "topleft", "top", "topright", "right", "center"))
  else NA
  
  if (is.na(auto)) {
    xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
    nx <- length(x)
    if (nx < 1 || nx > 2) stop("invalid coordinate lengths")
  } else nx <- 0
  
  xlog <- par("xlog")
  ylog <- par("ylog")
  
  rect2 <- function(left, top, dx, dy, density = NULL, angle, ...) {
    r <- left + dx; if(xlog) { left <- 10^left; r <- 10^r }
    b <- top  - dy; if(ylog) {  top <- 10^top;  b <- 10^b }
    rect(left, top, r, b, angle = angle, density = density, ...)
  }
  segments2 <- function(x1, y1, dx, dy, ...) {
    x2 <- x1 + dx; if(xlog) { x1 <- 10^x1; x2 <- 10^x2 }
    y2 <- y1 + dy; if(ylog) { y1 <- 10^y1; y2 <- 10^y2 }
    segments(x1, y1, x2, y2, ...)
  }
  points2 <- function(x, y, ...) {
    if(xlog) x <- 10^x
    if(ylog) y <- 10^y
    points(x, y, ...)
  }
  text2 <- function(x, y, ...) {
    ##--- need to adjust  adj == c(xadj, yadj) ?? --
    if(xlog) x <- 10^x
    if(ylog) y <- 10^y
    text(x, y, ...)
  }
  if(trace)
    catn <- function(...)
      do.call("cat", c(lapply(list(...),formatC), list("\n")))
  
  cin <- par("cin")
  Cex <- cex * par("cex")		# = the `effective' cex for text
  
  ## at this point we want positive width even for reversed x axis.
  if(is.null(text.width))
    text.width <- max(abs(strwidth(legend, units="user",
                                   cex=cex, font = text.font)))
  else if(!is.numeric(text.width) || text.width < 0)
    stop("'text.width' must be numeric, >= 0")
  
  xc <- Cex * xinch(cin[1L], warn.log=FALSE) # [uses par("usr") and "pin"]
  yc <- Cex * yinch(cin[2L], warn.log=FALSE)
  if(xc < 0) text.width <- -text.width
  
  xchar  <- xc
  xextra <- 0
  yextra <- yc * (y.intersp - 1)
  ## watch out for reversed axis here: heights can be negative
  ymax   <- yc * max(1, strheight(legend, units="user", cex=cex)/yc)
  ychar <- yextra + ymax
  if(trace) catn("  xchar=", xchar, "; (yextra,ychar)=", c(yextra,ychar))
  
  if(mfill) {
    ##= sizes of filled boxes.
    xbox <- xc * box.cex[1]
    ybox <- yc * box.cex[2]
    dx.fill <- xbox ## + x.intersp*xchar
  }
  do.lines <- (!missing(lty) && (is.character(lty) || any(lty > 0))
  ) || !missing(lwd)
  
  ## legends per column:
  n.legpercol <-
    if(horiz) {
      if(ncol != 1)
        warning(gettextf("horizontal specification overrides: Number of columns := %d",
                         n.leg), domain = NA)
      ncol <- n.leg
      1
    } else ceiling(n.leg / ncol)
  
  has.pch <- !missing(pch) && length(pch) > 0 # -> default 'merge' is available
  if(do.lines) {
    x.off <- if(merge) -0.7 else 0
  } else if(merge)
    warning("'merge = TRUE' has no effect when no line segments are drawn")
  
  if(has.pch) {
    if(is.character(pch) && !is.na(pch[1L]) &&
       nchar(pch[1L], type = "c") > 1) {
      if(length(pch) > 1)
        warning("not using pch[2..] since pch[1L] has multiple chars")
      np <- nchar(pch[1L], type = "c")
      pch <- substr(rep.int(pch[1L], np), 1L:np, 1L:np)
    }
    ## this coercion was documented but not done in R < 3.0.0
    if(!is.character(pch)) pch <- as.integer(pch)
  }
  
  if (is.na(auto)) {
    ##- Adjust (x,y) :
    if (xlog) x <- log10(x)
    if (ylog) y <- log10(y)
  }
  if(nx == 2) {
    ## (x,y) are specifiying OPPOSITE corners of the box
    x <- sort(x)
    y <- sort(y)
    left <- x[1L]
    top  <- y[2L]
    w <- diff(x)# width
    h <- diff(y)# height
    w0 <- w/ncol # column width
    
    x <- mean(x)
    y <- mean(y)
    if(missing(xjust)) xjust <- 0.5
    if(missing(yjust)) yjust <- 0.5
    
  }
  else {## nx == 1  or  auto
    ## -- (w,h) := (width,height) of the box to draw -- computed in steps
    h <- (n.legpercol + !is.null(title)) * ychar + yc
    w0 <- text.width + (x.intersp + 1) * xchar
    if(mfill)	w0 <- w0 + dx.fill
    if(do.lines)		w0 <- w0 + (seg.len + x.off)*xchar
    w <- ncol*w0 + .5* xchar
    if (!is.null(title)
        && (abs(tw <- strwidth(title, units="user", cex=cex) + 0.5*xchar)) > abs(w)) {
      xextra <- (tw - w)/2
      w <- tw
    }
    
    ##-- (w,h) are now the final box width/height.
    
    if (is.na(auto)) {
      left <- x - xjust * w
      top	 <- y + (1 - yjust) * h
    } else {
      usr <- par("usr")
      inset <- rep_len(inset, 2)
      insetx <- inset[1L]*(usr[2L] - usr[1L])
      left <- switch(auto, "bottomright" =,
                     "topright" =, "right" = usr[2L] - w - insetx,
                     "bottomleft" =, "left" =, "topleft" = usr[1L] + insetx,
                     "bottom" =, "top" =, "center" = (usr[1L] + usr[2L] - w)/2)
      insety <- inset[2L]*(usr[4L] - usr[3L])
      top <- switch(auto, "bottomright" =,
                    "bottom" =, "bottomleft" = usr[3L] + h + insety,
                    "topleft" =, "top" =, "topright" = usr[4L] - insety,
                    "left" =, "right" =, "center" = (usr[3L] + usr[4L] + h)/2)
    }
  }
  
  if (plot && bty != "n") { ## The legend box :
    if(trace)
      catn("  rect2(", left, ",", top,", w=", w, ", h=", h, ", ...)",
           sep = "")
    rect2(left, top, dx = w, dy = h, col = bg, density = NULL,
          lwd = box.lwd, lty = box.lty, border = box.col)
  }
  
  ## (xt[],yt[]) := `current' vectors of (x/y) legend text
  xt <- left + xchar + xextra +
    (w0 * rep.int(0:(ncol-1), rep.int(n.legpercol,ncol)))[1L:n.leg]
  yt <- top -	0.5 * yextra - ymax -
    (rep.int(1L:n.legpercol,ncol)[1L:n.leg] - 1 + !is.null(title)) * ychar
  
  if (mfill) {		#- draw filled boxes -------------
    if(plot) {
      if(!is.null(fill)) fill <- rep_len(fill, n.leg)
      rect2(left = xt, top=yt+ybox/2, dx = xbox, dy = ybox,
            col = fill,
            density = density, angle = angle, border = border)
    }
    xt <- xt + dx.fill
  }
  if(plot && (has.pch || do.lines))
    col <- rep_len(col, n.leg)
  
  ## NULL is not documented but people use it.
  if(missing(lwd) || is.null(lwd))
    lwd <- par("lwd") # = default for pt.lwd
  if (do.lines) {			#- draw lines ---------------------
    ## NULL is not documented
    if(missing(lty) || is.null(lty)) lty <- 1
    lty <- rep_len(lty, n.leg)
    lwd <- rep_len(lwd, n.leg)
    ok.l <- !is.na(lty) & (is.character(lty) | lty > 0) & !is.na(lwd)
    if(trace)
      catn("  segments2(",xt[ok.l] + x.off*xchar, ",", yt[ok.l],
           ", dx=", seg.len*xchar, ", dy=0, ...)")
    if(plot)
      segments2(xt[ok.l] + x.off*xchar, yt[ok.l],
                dx = seg.len*xchar, dy = 0,
                lty = lty[ok.l], lwd = lwd[ok.l], col = col[ok.l])
    # if (!merge)
    xt <- xt + (seg.len+x.off) * xchar
  }
  if (has.pch) {			#- draw points -------------------
    pch <- rep_len(pch, n.leg)
    pt.bg <- rep_len(pt.bg, n.leg)
    pt.cex <- rep_len(pt.cex, n.leg)
    pt.lwd <- rep_len(pt.lwd, n.leg)
    ok <- !is.na(pch)
    if (!is.character(pch)) {
      ## R 2.x.y omitted pch < 0
      ok <- ok & (pch >= 0 | pch <= -32)
    } else {
      ## like points
      ok <- ok & nzchar(pch)
    }
    x1 <- (if(merge && do.lines) xt-(seg.len/2)*xchar else xt)[ok]
    y1 <- yt[ok]
    if(trace)
      catn("  points2(", x1,",", y1,", pch=", pch[ok],", ...)")
    if(plot)
      points2(x1, y1, pch = pch[ok], col = col[ok],
              cex = pt.cex[ok], bg = pt.bg[ok], lwd = pt.lwd[ok])
    ##D	if (!merge) xt <- xt + dx.pch
  }
  
  xt <- xt + x.intersp * xchar
  if(plot) {
    if (!is.null(title))
      text2(left + w*title.adj, top - ymax, labels = title,
            adj = c(title.adj, 0), cex = cex, col = title.col)
    
    text2(xt, yt, labels = legend, adj = adj, cex = cex,
          col = text.col, font = text.font)
  }
  invisible(list(rect = list(w = w, h = h, left = left, top = top),
                 text = list(x = xt, y = yt)))
} #----------------------------

palf <- colorRampPalette(c("yellow", "blue"))
Ylim <- c(-13, 1)
Cex2 = 2

png(filename = paste("/home/anja/Desktop/ReMa/FDA/C_R/plots/","histogram",".png",sep=""), 
    width = 960, height = 960)
layout(matrix (1:9, 3, 3, byrow = T))
par(mar = c(4.5, 6.5, 8, 0.5), cex.main = 2, cex.axis = 2)

d = 2
for(e in 1:length(ES)){
  for(t in 1:length(Trails)){
    p <- hist(t(log(LowHighBF01[1, 1, 1, e, t, d, ], base = 10)), breaks = (Ylim[2]+.4):(Ylim[1]-.4), plot = F)
    p2 <- hist(t(log(LowHighBF01[1, 1, 2, e, t, d, ], base = 10)), breaks = (Ylim[2]+.2):(Ylim[1]-.2), plot = F)
    p3 <- hist(t(log(LowHighBF01[1, 1, 3, e, t, d, ], base = 10)), breaks = Ylim[2]:Ylim[1], plot = F)
    plot(p, col = adjustcolor(palf(100)[15], alpha=0.3), xlab = "", xaxt='n', ann=FALSE, cex.axis = Cex2, main = paste ("Trials:", Trails[t]), pch = 16, xlim =Ylim, ylim = c(0, 250))
    plot(p2, col= adjustcolor(palf(100)[20], .15), axes = FALSE, xlab = "", ylab = "", main = "", xlim = Ylim, ylim = c(0, 250), add = T)
    plot(p3, col= adjustcolor(palf(100)[100], .08), axes = FALSE, xlab = "", ylab = "", main = "", xlim = Ylim, ylim = c(0, 250), add = T)
    axis(1, at = Ylim[1]:Ylim[2], labels = 10^(Ylim[1]:Ylim[2]), las = 1, cex= Cex2)
    abline(v =t(log(.05, base = 10)), lty = 2, lwd = 2)
    Legend(x= -13, y= 270, legend = c(paste("N =", N[1]), paste("N =", N[2]), paste("N =", N[3]))
           , fill= c(adjustcolor(palf(100)[15], 0.3), 
                     adjustcolor(palf(100)[25], .15),
                     adjustcolor(palf(100)[100], .08)
           ),
           cex =Cex2, bty = "n",
           x.intersp = .7, y.intersp = 1.4, box.cex=c(1.3 , 1.3)
    )
    mtext(text = expression(paste("Lowest Interval ", BF[0][1])), side = 1, line = 3, cex =Cex1)
    if(t%%2 == 0){
      mtext(text = paste("ES =", ES[e]), side = 3, line = 6, cex = Cex3, font = 2)
    }
  }
}

dev.off ()




#Equivalence boxplots----
N = c(20, 40, 80) # sample size, index with n
ES = c(0, .2, .3) # effect size, index with e (is mean, sd is always 1)
Trails = 2# number of Trails, not varied
D <- matrix(c(-.2, -.3, .2, .3), 2) # inferiority-margin delta, not varied
Runs = 500 # index with r

options(scipen = 1)
Ylim <- c(-1, 3)
Cex1 <- 1.3
Cex2 <- .8 #axis ticks
Cex3  <- 1.3 #ES heading

png(filename = paste("/home/anja/Desktop/ReMa/FDA/C_R/plots/","boxplotEquivalence",".png",sep=""), 
    width = 960, height = 960)
layout(matrix (1:6, 2, 3, byrow = T))
par(mar = c(4.5, 6.5, 8, 0.5), cex.main = 2, cex.axis = 2)


for(d in 1:length(D[,1])){
  for(e in 1:length(ES)){
    boxplot(col = adjustcolor("white", alpha=0.2), border= "red",  at = c(1, 2, 3), x = t(log(EquivBF01[1, d, , e, ], base = 10)), 
            axes = F, ylim = Ylim, xlim= c(0,4), xlab = "", ylab = "", main = paste ("Effect Size:", ES[e]), pch = 16)
    abline(h = 0, lty = 2, lwd = 2)
    axis(1, at = c(0, 1.1:3.1, 4), labels = c( "", N, ""))
    axis(2, at = Ylim[1]:Ylim[2], labels = 10^(Ylim[1]:Ylim[2]), las = 1)
    mtext(text = "N", side = 1, line = 3, cex =Cex1)
    boxplot(add= TRUE, col=adjustcolor("white", .2), border= "blue",at = c(1.2, 2.2, 3.2), x = t(log(EquivBF01[2, d, , e, ], base = 10)),
            axes = F, ylim = Ylim, xlim= c(0,4), xlab = "", ylab = "", main = paste ("Effect Size:", ES[e]), pch = 16)
    if(e == 1){
      mtext(text = paste("Equivalence Region: (", D[d, 1], " - ", D[d, 2], ")"), side = 3, line = 5.4, cex = 1.5, font = 2)
    }
  }
}
dev.off ()




