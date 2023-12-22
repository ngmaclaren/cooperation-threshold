                                        # For a brute force search through possible models
library(MuMIn)

                                        # Load data
modelframe <- read.csv("./data/modelframe.csv")
                                        # Settings
save_plots <- FALSE # TRUE
pal <- "Tableau 10"
palette(pal)
pal.colors <- palette.colors(palette = pal)

                                        # print data table
with(list(x = modelframe), {
    x$Brain <- exp(x$Brain)
    x$smean <- exp(x$logsmean)
    x$wclust <- exp(x$logwclust)

    cbind(
        x[, c("genus", "species")],
        round(x[, c("bc", "ncr", "Brain", "N", "kmean", "smean", "clust", "wclust")], 2)
    )
})

                                        # Collect all predictors to make a full model for dredge
lhs <- "bc ~ "
xterms <- colnames(modelframe)[-which(colnames(modelframe) %in% c("genus", "species", "bc"))]
rhs <- paste(xterms, collapse = " + ")
form <- formula(paste0(lhs, rhs))
                                        # A positive, right-skewed, continuous variable may be modeled
                                        # with gamma-distributed errors
initialmodel <- glm(form, Gamma, modelframe, na.action = "na.fail")
                                        # Make versions of the full model with other reasonable
                                        # choices
modeloptions <- list(
    initialmodel,
    update(initialmodel, family = Gamma(link = "log")),
    update(initialmodel, family = gaussian(link = "log")), # does not converge
    update(initialmodel, family = gaussian(link = "inverse")), # does not converge
    update(initialmodel, family = quasipoisson)
)
names(modeloptions) <- c(
    "Gamma, Inverse",
    "Gamma, Log",
    "Gaussian, Log",
    "Gaussian, Inverse",
    "Quasipoisson"
)                                        # Calculate estimates of model fit, Faraway ~p. 157
                                        # Small deviance reflects better fit.
dev <- sapply(modeloptions, deviance)
print(dev)
1 - pchisq(dev, sapply(modeloptions, df.residual))

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
                                        # and have from 0 to five predictor variables.
                                        # Choose the minimum deviance model as the base for dredging.

                                        # Changes Friday, September 1, 2023
                                        # Remove limitation on number of model terms.
                                        # , m.lim = c(0, 5)
                                        # Add limitation that the two brain size variables cannot
                                        # appear in the same model
chosen_model <- modeloptions[[which.min(dev)]]
termsmatrix <- matrix(
    TRUE, ncol = length(xterms), nrow = length(xterms), dimnames = list(xterms, xterms)
)
termsmatrix[upper.tri(termsmatrix)] <- NA
diag(termsmatrix) <- NA
termsmatrix["Brain", "ncr"] <- FALSE
dredged <- dredge(chosen_model, subset = termsmatrix)
                                        # Alternately, choose the canonical link, with slightly poorer
                                        # deviance
## dredged <- dredge(modeloptions[["Gamma, Inverse"]], subset = termsmatrix)

                                        # Make the AIC comp figure here
ht <- 7; wd <- 7
if(save_plots) {
    pdf("./img/AICcomp.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
par(mar = c(4, 4, 1, 1) + 0.5)
plot(
    seq(nrow(dredged)), dredged$AICc, type = "o", lwd = 2, lty = 1, col = 1,
    xlab = "Model rank", ylab = "AICc", cex.axis = 1.5, cex.lab = 1.5,
    xlim = range(seq(nrow(dredged))), ylim = range(dredged$AICc)
)
plotrix::axis.break(2, min(dredged$AICc)-0.5, style = "zigzag")
if(save_plots) dev.off()

                                        # Focus on models that have an AICc value no greater than 3
                                        # more than the best model.
bestmodels <- get.models(dredged, subset = delta < 3)
                                        # Inspect basic summary
lapply(bestmodels, summary)
                                        # and profile likelihood confidence intervals
lapply(bestmodels, confint)
                                        # check variance inflation factors
lapply(bestmodels, car::vif)
car::vif(chosen_model)
                                        # pseudo-R^2
sapply(bestmodels, function(m) 1 - (m$deviance/m$null.deviance)) # McFadden's
sapply(bestmodels, function(m) (m$null.deviance - m$deviance)/m$null.deviance) # Cohen's

                                        # Best model including the NCR term (Model 2)
best <- bestmodels[[2]]

                                        # Diagnostics
## plot(best)

## plot(
##     data.frame(
##         Response = log(modelframe$bc),
##         kmean = modelframe$kmean,
##         NCR = modelframe$ncr)
## )

##                                         # demonstrate effect of removing influential data point
## best_alt <- update(best, subset = species != "campbelli")
## summary(best_alt)
## plot(best_alt)

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
    pdf("./img/modelfig.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
par(mar = c(4, 4, 1, 1) + 0.5)
plot(
    bc ~ ncr, data = modelframe, type = "p", # log = "y",
    pch = 19, cex = 2, col = pal.colors["lightgray"], #"gray50",
    xlab = "Neocortex ratio", ylab = "Cooperation threshold (b/c)*",
    cex.lab = 1.5, cex.axis = 1.5
)
add_model(best, "ncr", pdata, pal.colors["blue"], pal.colors["blue"])
## add_model(best, "Brain", pdata, 1, 1)
                                        # Influential data point
## points(modelframe$Brain[3], modelframe$bc[3], col = 2, pch = 19, cex = 3)
if(save_plots) dev.off()

                                        # Plot the coefficients of all the best models, and the 95%
                                        # confidence intervals of those coefficients
palette("Tableau 10")
yticklabels <- c(
    "Neocortex ratio", #"NCR",
    "Brain mass", # "ln(Brain mass)",
    "Average degree", # expression(group(langle, italic(k), rangle)),#"〈k〉",
    "Weighted\nclustering\ncoefficient" # expression(italic(tilde(C))[w])
)
matcher <- c(
    "ncr",
    "Brain",
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
    pdf("./img/searchfig.pdf", height = ht, width = wd)
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
adjusts <- seq(.25, -.25, length.out = howmany)
abline(h = seq(min(ypos$pos) + .5, max(ypos$pos) - .5, by = 1),
       lwd = 1, lty = 2, col = pal.colors["lightgray"])
for(i in 1:length(bestmodels)) {
    plot_coefs(
        model = bestmodels[[i]],
        ptsize = 2, color = i + 2, pch = pchs[i], adjust = adjusts[i]
    )
}
text(2, .5, "Makes cooperation less likely", adj = 1, cex = 1.25)
text(-2, .5, "Makes cooperation more likely", adj = 0, cex = 1.25)
legend(
    "topright", bty = "n", ncol = howmany,
    legend = paste("Model", 1:howmany, "      "),
    pch = pchs, col = 1:howmany + 2, pt.cex = 1.5, lwd = 2, pt.lwd = 2, cex = 1.25
)
if(save_plots) dev.off()



                                        # Print summaries of the models with ΔAICc < 3 (top five)
summarylocal <- function(model) {
    print(round(summary(model)$coefficients[, c(1, 2, 4)], 3))
    print(round(AICc(model), 3))
}

lapply(get.models(dredged, delta < 3), summarylocal)
