### The log link has a strictly smaller deviance than the inverse link, suggesting better fit.
### Using the log link tightens the confidence intervals.
### Results are overall similar, and the best two models choose the same predictors.
### Furthermore, it identifies a single outlier, and has better diagnostic plots for the best model.
### That single outlier is Cercopithecus campbelli, which has a network with a large (b/c)*. 
### Dropping that outlier, which we otherwise have no particular reason to do, gives good coefficients across the board.
### summary(glm(bc ~ Brain + kmean + logwclust, data = modelframe, family = Gamma("log"), subset = -3))


                                        # For a brute force search through possible models
library(MuMIn)

                                        # Load data
modelframe <- read.csv("./data/modelframe.csv")

                                        # settings
save_plots <- FALSE # TRUE
pal <- "Tableau 10"
palette(pal)
pal.colors <- palette.colors(palette = pal)

                                        # print data table
cbind(
    modelframe[, c("genus", "species")],
    round(modelframe[, c("Brain", "ncr", "N", "kmean", "logsmean", "clust", "logwclust")], 2)
)

                                        # Collect all predictors to make a full model from which to
                                        # dredge
lhs <- "bc ~ "
xterms <- colnames(modelframe)[4:ncol(modelframe)]
rhs <- paste(xterms, collapse = " + ")
form <- formula(paste0(lhs, rhs))
                                        # A positive, right-skewed, continuous variable may be modeled
                                        # with gamma-distributed errors
fullmodel <- glm(form, Gamma, modelframe, na.action = "na.fail")
                                        # Make versions of the full model with other reasonable
                                        # choices
modeloptions <- list(
    fullmodel,
    update(fullmodel, family = gaussian(link = "log")),
    update(fullmodel, family = gaussian(link = "inverse")),
    update(fullmodel, family = Gamma(link = "log")),
    update(fullmodel, family = quasipoisson)
)
                                        # Check model fit. Small deviance reflects better optimization
                                        # results.
dev <- sapply(modeloptions, deviance) # deviance is a measure of model fit
                                        # Calculate estimates of model fit, Faraway ~p. 157
print(dev)
1 - pchisq(dev, sapply(modeloptions, df.residual))

names(modeloptions) <- c(
    "Gamma, Inverse",
    "Gaussian, Log",
    "Gaussian, Inverse",
    "Gamma, Log",
    "Quasipoisson"
)
                                        # Check residuals for all chosen models.
                                        # Gaussian models are clearly poor from visual inspection,
                                        # confirming deviance analysis. 
## dev.new(height = 10, width = 15); par(mfrow = c(2, 3))
## for(i in 1:length(modeloptions)) {
##     model <- modeloptions[[i]]
##     plot(
##         residuals(model) ~ predict(model, type = "link"),
##         xlab = expression(hat(eta)), ylab = "Deviance residuals", main = names(modeloptions)[i]
##     )
## }

## Dredge
                                        # A single function call from the MuMIn library searches all
                                        # models which are subsets of the full model specified above
                                        # and have from 0 to five predictor variables. There are
                                        # several warnings having to do with missing values produced
                                        # in model fit. These would be poor models anyway.
                                        # ---> VERIFY

                                        # to forbid Brain and ncr from being in the same model
## termmatrix <- matrix(
##     NA, nrow = length(xterms), ncol = length(xterms), dimnames = list(sort(xterms), sort(xterms))
## )
## termmatrix[lower.tri(termmatrix)] <- TRUE
## termmatrix["ncr", "Brain"] <- FALSE
## dredged <- dredge(fullmodel, m.lim = c(0, 5))#, subset = termmatrix)
                                        # Choose the minimum deviance model
dredged <- dredge(modeloptions[[which.min(dev)]], m.lim = c(0, 5))
## dredged <- dredge(modeloptions[["Gamma, Inverse"]], m.lim = c(0, 5))
                                        # Focus on models that have an AICc value no greater than 3
                                        # more than the best model.
bestmodels <- get.models(dredged, subset = delta < 3)


                                        # Best model including the NCR term (Model 2)
best <- 
                                        # demonstrate effect of removing influential data point
## best <- update(bestmodels[[1]], subset = species != "campbelli")
summary(best)
                                        # Diagnostics
## dev.new()
## plot(residuals(best) ~ predict(best, type = "link"),
##      xlab = expression(hat(eta)), ylab = "Deviance Residuals")

## dev.new(height = 11, width = 11)
## plot(
##     data.frame(
##         Response = 1/modelframe$bc,
##         kmean = modelframe$kmean,
##         NCR = modelframe$ncr)
## )

## dev.new()
## par(mfrow = c(2, 2))
## plot(best)


                                        # Plot the data and model fit w/ 2*SE for the best model which
                                        # includes NCR as a predictor
add_model <- function(model, xvar, pdata, color1, color2) {
    y <- predict(model, pdata, type = "response", se = TRUE)
    lines(pdata[, xvar], y$fit, lwd = 3, col = color1)
    lines(pdata[, xvar], y$fit + 2*y$se.fit, lwd = 1.5, col = color2, lty = 3)
    lines(pdata[, xvar], y$fit - 2*y$se.fit, lwd = 1.5, col = color2, lty = 3)
}

pdata <- data.frame(
    logwclust = median(modelframe$logwclust),
    kmean = median(modelframe$kmean),
    ncr = seq(min(modelframe$ncr), max(modelframe$ncr), length.out = 50),
    Brain = seq(min(modelframe$Brain), max(modelframe$Brain), length.out = 50)
)
ht <- 7; wd <- 7
if(save_plots) {
    cairo_pdf("./img/modelfig.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
par(mar = c(4, 4, 1, 1) + 0.5)
plot(
    bc ~ Brain, data = modelframe, type = "p", log = "y",
    pch = 19, cex = 2, col = pal.colors["lightgray"], #"gray50",
    #ylim = c(0, 1.05*max(modelframe$bc)), xlim = range(modelframe$ncr),#c(2, 3.5),
    xlab = "Neocortex ratio (NCR)", ylab = "Cooperation threshold (b/c)*",
    cex.lab = 1.5, cex.axis = 1.5
)
##add_model(best, "ncr", pdata, pal.colors["blue"], pal.colors["blue"])
add_model(best, "Brain", pdata, 1, 1)
points(modelframe$Brain[3], modelframe$bc[3], col = 2, pch = 19, cex = 3)
if(save_plots) dev.off()

                                        # Plot the coefficients of all the best models, and the 95%
                                        # confidence intervals of those coefficients
palette("Tableau 10")
yticklabels <- c(
    "Neocortex ratio", #"NCR",
    "Brain mass", # "ln(Brain mass)",
    "Body mass", #"ln(Body mass)",
    "Average degree", # expression(group(langle, italic(k), rangle)),#"〈k〉",
    "Weighted\nclustering\ncoefficient" # expression(italic(tilde(C))[w])
)
matcher <- c(
    "ncr",
    "Brain",
    "Body",
    "kmean",
    "logwclust"
)
ypos <- data.frame(
    name = matcher,
    pos = rev(1:length(yticklabels))
)
plot_coefs <- function(model, ptsize, color, pch, adjust) {
    coefs <- as.data.frame(coefficients(model))
    cis <- confint(model)
    if(nrow(coefs) == 1) {
        coefs[, 2] <- cis[1]
        coefs[, 3] <- cis[2]
    } else {
        coefs <- cbind(coefs, cis)
    }
    coefs$ypos <- ypos$pos[match(rownames(coefs), ypos$name)]
    pch <- pch; lwd <- ptsize*1.5; cex <- ptsize; lty <- 1
    y <- coefs[-1, 4] + adjust
    points(x = coefs[-1, 1], y = y, pch = pch, lwd = lwd, cex = cex, col = color)
    segments(x0 = coefs[-1, 2], x1 = coefs[-1, 3], y0 = y, lwd = lwd, col = color, lty = 1)
}
##xlim <- c(-.3, .3)
xlim <- c(-2, 2)
ylim <- range(ypos[-9, "pos"]) + c(-.5, .5)
ht <- 7
wd <- 14
pchs <- c(0, 1, 2, 5, 6)
if(save_plots) {
    cairo_pdf("./img/searchfig.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
par(mar = c(4, 11, 2, 2))
plot(NULL, xlim = xlim, ylim = ylim, axes = FALSE, xlab = "", ylab = "")
box()
axis(1, cex.lab = 1.5, cex.axis = 1.5)
axis(2, at = ypos$pos, labels = yticklabels, las = 2, cex.lab = 1.5, cex.axis = 1.5) # ypos$label
title(xlab = "Coefficient value", cex.lab = 1.5)
abline(v = 0, col = "black", lty = 1, lwd = .75)
howmany <- length(bestmodels)
adjusts <- seq(.33, -.33, length.out = howmany)
abline(h = seq(min(ypos$pos) + .5, max(ypos$pos) - .5, by = 1),
       lwd = 1, lty = 2, col = pal.colors["lightgray"])
for(i in 1:length(bestmodels)) {
    plot_coefs(
        model = bestmodels[[i]],
        ptsize = 2, color = i + 2, pch = pchs[i], adjust = adjusts[i]
    )
}
text(-.3, .5, "Makes cooperation less likely", adj = 0, cex = 1.25)
text(.3, .5, "Makes cooperation more likely", adj = 1, cex = 1.25)
legend(
    "topright", bty = "n", ncol = howmany,
    legend = paste("Model", 1:howmany, "      "),
    pch = pchs, col = 1:howmany + 2, pt.cex = 1.5, lwd = 2, pt.lwd = 2, cex = 1.25
)
if(save_plots) dev.off()

                                        # Show the AIC results
ht <- 7; wd <- 7
if(save_plots) {
    cairo_pdf("./img/AICcomp.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
par(mar = c(4, 4, 1, 1) + 0.5)
plot(
    seq(nrow(dredged)), dredged$AICc, type = "o", lwd = 2, lty = 1, col = 1,
    ##NULL,
    xlab = "Model rank", ylab = "AICc", cex.axis = 1.5, cex.lab = 1.5,
    xlim = range(seq(nrow(dredged))), ylim = range(dredged$AICc)
)
## lines(seq(nrow(dredged)), dredged$AICc, col = pal.colors["lightgray"], lwd = 2, lty = 1)
## points(
##     6:nrow(dredged), dredged$AICc[6:nrow(dredged)],
##     pch = 1, lwd = 2, cex = 2, col = pal.colors["lightgray"]
## )
## points(1:5, dredged$AICc[1:5], col = 1:5, pch = 1, cex = 2, lwd = 4)
if(save_plots) dev.off()


                                        # Print summaries of the models with ΔAICc < 3 (top five)
summarylocal <- function(model) {
    print(round(summary(model)$coefficients[, c(1, 2, 4)], 3))
    print(round(AICc(model), 3))
}

lapply(get.models(dredged, delta < 3), summarylocal)
