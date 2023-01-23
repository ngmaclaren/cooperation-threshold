                                        # For a brute force search through possible models
library(MuMIn)
add_model <- function(model, xvar, pdata, color1, color2) {
    y <- predict(model, pdata, type = "response", se = TRUE)
    lines(pdata[, xvar], y$fit, lwd = 3, col = color1)
    lines(pdata[, xvar], y$fit + 2*y$se.fit, lwd = 1.5, col = color2, lty = 3)
    lines(pdata[, xvar], y$fit - 2*y$se.fit, lwd = 1.5, col = color2, lty = 3)
}

## file settings
save_plots <- FALSE # TRUE
pal <- "Tableau 10"
palette(pal)
pal.colors <- palette.colors(palette = pal)
##
##
## Data Prep
                                        # This data is a processed version of the output from
                                        # "evaluate-threshold.py"
df <- read.csv("./data/network-data.csv")
                                        # Focus on positive social interactions. These networks
                                        # exclude connections formed by proximity only and by sexual
                                        # and antagonistic interactions. 
interaction_types <- c("grooming", "physical_contact", "overall_mix")
                                        # This network should also be excluded (wrong type).
dropthis <- "macaques_pugagonzalez_4"
                                        # Focus on networks that fit the requirements listed so far
df <- df[
    df$interaction_type %in% interaction_types &
    df$order == "Primates" &
    !(df$filename %in% dropthis) # but negative b/c anyway
   ,
    ]
                                        # Load neocortex ratio data
ncr_genus <- read.csv("./data/ncr-genus.csv")
ncr_genus_species <- read.csv("./data/ncr-genus_species.csv")
                                        # Make a reference table of NCR values that prefers species
                                        # data if available, otherwise takes the genus-level data
ncr <- ncr_genus_species
for(i in 1:nrow(ncr)) {
    if(is.na(ncr[i, "ncr"])) {
        ncr[i, "ncr"] <- ncr_genus$ncr[ncr_genus$genus == ncr[i, "genus"]]
    }
}
                                        # Merge
df <- merge(df, ncr, by = c("genus", "species"), all.x = TRUE)
                                        # Same process with the brain/body size data
brainsize <- read.csv("./data/smaers-etal-2021-brainsize.csv")
                                        # make P. papio a copy of P. cynocephalus for brain/body size
Papiopapio <- data.frame( 
    Status = "Extant", Genus = "Papio", Species = "papio",
    Body = brainsize[brainsize$Species == "cynocephalus", "Body"],
    Brain = brainsize[brainsize$Species == "cynocephalus", "Brain"]
)
brainsize <- rbind(brainsize, Papiopapio)
                                        # Sapajus is a synonym for Cebus
brainsize[brainsize$Genus == "Cebus" & brainsize$Species == "apella", "Genus"] <- "Sapajus"
                                        # Merge
df <- merge(df, brainsize, by.x = c("genus", "species"), by.y = c("Genus", "Species"), all.x = TRUE)
                                        # Drop any NA values
                                        # one P. cynocephalus network with NA in b/c
df <- df[, -which(colnames(df) %in% "population_type")]
                                        # The following line drops baboon_franz_grooming_group_12, a
                                        # disconnected network
df <- df[complete.cases(df), ] 

                                        # Remove negative b/c ratios as irrelevant
dfr <- df[
    df$bc > 0, # this is fragile, as comparisons with NA always result in NA
    ]
                                        # Take the median to aggregate to one value per species. 
newdf <- aggregate(
    cbind(bc, N, kmean, clust, smean, wclust, ncr, Body, Brain) ~ species, data = dfr, median
)
                                        # Generate log-transform versions of the weighted variables
newdf$logsmean <- log(newdf$smean)
newdf$logwclust <- log(newdf$wclust)
                                        # and drop the un-log'ed versions
newdf <- newdf[, -which(colnames(newdf) %in% c("smean", "wclust"))]
                                        # For convenience in interactive use
newdf$genus <- dfr$genus[match(newdf$species, dfr$species)]
newdf <- newdf[, c(ncol(newdf), 1:(ncol(newdf) - 1))]

                                        # What is the correlation between N and b/c?
cor(newdf$bc, newdf$N)

## Full Model

modelframe <- newdf
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
dev.new(height = 10, width = 15); par(mfrow = c(2, 3))
for(i in 1:length(modeloptions)) {
    model <- modeloptions[[i]]
    plot(
        residuals(model) ~ predict(model, type = "link"),
        xlab = expression(hat(eta)), ylab = "Deviance residuals", main = names(modeloptions)[i]
    )
}

## Dredge
                                        # A single function call from the MuMIn library searches all
                                        # models which are subsets of the full model specified above
                                        # and have from 0 to four predictor variables. There are
                                        # several warnings having to do with missing values produced
                                        # in model fit. These would be poor models anyway.
                                        # ---> VERIFY
dredged <- dredge(fullmodel, m.lim = c(0, 5))
                                        # Focus on models that have an AICc value no greater than 3
                                        # more than the best model.
bestmodels <- get.models(dredged, subset = delta < 3)

                                        # Show the AIC results
ht <- 7; wd <- 7
if(save_plots) {
    cairo_pdf("./img/AICcomp.pdf", height = ht, width = wd, family = "Bitstream Vera Sans")
} else {
    dev.new(height = ht, width = wd)
}
par(mar = c(4, 4, 1, 1) + 0.5)
plot(1:nrow(dredged), dredged$AICc, type = "o",
     pch = 1, lwd = 2, col = pal.colors["blue"],
     xlab = "Model rank", ylab = "AICc", cex.axis = 1.5, cex.lab = 1.5)
if(save_plots) dev.off()

                                        # Best model including the ncr term
best <- bestmodels[[3]]
summary(best)
                                        # Diagnostics
dev.new()
plot(residuals(best) ~ predict(best, type = "link"),
     xlab = expression(hat(eta)), ylab = "Deviance Residuals")

dev.new(height = 11, width = 11)
plot(
    data.frame(
        Response = 1/modelframe$bc,
        kmean = modelframe$kmean,
        NCR = modelframe$ncr)
)

dev.new()
par(mfrow = c(2, 2))
plot(best)


                                        # Plot the data and model fit w/ 2*SE for the best model which
                                        # includes NCR as a predictor
pdata <- data.frame(
    logwclust = median(modelframe$logwclust),
    kmean = median(modelframe$kmean),
    ncr = seq(min(modelframe$ncr), max(modelframe$ncr), length.out = 50)
)
ht <- 7; wd <- 7
if(save_plots) {
    cairo_pdf("./img/modelfig.pdf", height = ht, width = wd, family = "Bitsream Vera Sans")
} else {
    dev.new(height = ht, width = wd)
}
par(mar = c(4, 4, 1, 1) + 0.5)
plot(
    bc ~ ncr, data = modelframe, type = "p",
    pch = 19, cex = 2, col = pal.colors["lightgray"], #"gray50",
    ylim = c(0, 1.05*max(modelframe$bc)), xlim = range(modelframe$ncr),#c(2, 3.5),
    xlab = "Neocortex ratio (NCR)", ylab = "Cooperation threshold (b/c)*",
    cex.lab = 1.5, cex.axis = 1.5
)
add_model(best, "ncr", pdata, pal.colors["blue"], pal.colors["blue"])
if(save_plots) dev.off()

                                        # Plot the coefficients of all the best models, and the 95%
                                        # confidence intervals of those coefficients
palette("Tableau 10")
yticklabels <- c(
    "NCR",
    "ln(Brain mass)",
    "ln(Body mass)",
    "〈k〉",
    expression(italic(tilde(C))[w])
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
plot_coefs <- function(model, ptsize, color, adjust) {
    coefs <- as.data.frame(coefficients(model))
    cis <- confint(model)
    if(nrow(coefs) == 1) {
        coefs[, 2] <- cis[1]
        coefs[, 3] <- cis[2]
    } else {
        coefs <- cbind(coefs, cis)
    }
    coefs$ypos <- ypos$pos[match(rownames(coefs), ypos$name)]
    pch <- 1; lwd <- ptsize*1.5; cex <- ptsize; lty <- 1
    y <- coefs[-1, 4] + adjust
    points(x = coefs[-1, 1], y = y, pch = pch, lwd = lwd, cex = cex, col = color)
    segments(x0 = coefs[-1, 2], x1 = coefs[-1, 3], y0 = y, lwd = lwd, col = color, lty = 1)
}
xlim <- c(-.3, .3)
ylim <- range(ypos[-9, "pos"]) + c(-.5, .5)
ht <- 7
wd <- 14
if(save_plots) {
    cairo_pdf("./img/searchfig.pdf", height = ht, width = wd, family = "Bitsream Vera Sans")
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
        ptsize = 2, color = i, adjust = adjusts[i]
    )
}
text(-.3, .5, "Makes cooperation less likely", adj = 0, cex = 1.25)
text(.3, .5, "Makes cooperation more likely", adj = 1, cex = 1.25)
legend(
    "topright", bty = "n", ncol = howmany,
    legend = paste("Model", 1:howmany, "      "),
    pch = 1, col = 1:howmany, pt.cex = 1.5, lwd = 1.5, pt.lwd = 1.5, cex = 1.25
)
if(save_plots) dev.off()

                                        # Print summaries of the models with ΔAICc < 3 (top five)
summarylocal <- function(model) {
    print(round(summary(model)$coefficients[, c(1, 2, 4)], 3))
    print(round(AICc(model), 3))
}

lapply(get.models(dredged, delta < 3), summarylocal)
