                                        # Data Prep
                                        #
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

                                        # Zero-order correlations of interest
cor(newdf$bc, newdf$N)
cor(newdf$Brain, newdf$Body)
cor(newdf$ncr, newdf$N)

                                        # Full Model
                                        # Drop body size due to high correlation with brain size
modelframe <- newdf[, -grep("Body", colnames(newdf))]

write.csv(modelframe, file = "./data/modelframe.csv", row.names = FALSE)
