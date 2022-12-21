####ALLISON DENNERT, Pollinator Visit Analyses
####"Marine-derived nutrients from salmon indirectly affect pollinator visits 
####in a coastal wildflower meadow"
####by A. Dennert, E. Elle, J. Reynolds

#This R code:
#1) creates a PCA of plant floral/inflorescence availability
#splices PC scores to pollinator data
#2) manipulates/cleans pollinator visit data into analyzable format, 
#3) creates various stacked bar charts of cumulative visits by insects
#4) makes summary calculations of those cumulative visits
#5) uses structural equation modelling to model floral availability and total 
#insect visits to a wildflower meadow field experiment with 4 treatments: salmon, 
#algae, both, control. The experiment had 25 randomized blocks, with 4 treatments 
#within a block. Number of visits were measured in 10 minute intervals per 
#observation period, and floral availabilitycounted before each 10 minute interval. 

####0.0 PACKAGE LOADING####
set.seed(1)
library(glmmTMB)
library(taRifx)
library(car)
library(MASS)
library(DHARMa)
library(emmeans)
library(plotrix)
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(ggfortify)
library(PNWColors)
library(patchwork)
#packages for SEM/path analysis
library(semPlot)
library(OpenMx)
library(knitr)
library(kableExtra)
library(GGally)
library(piecewiseSEM)

####1.0 PLANT PCA AND DATA SPLICING####
#This code takes a matrix of available floral resources (open flowers and/or 
#inflorescences at the time of the pollinator surverys) and reduces the number of
#axes so that they may be modelled as a predictor in later pollinator visit 
#models. This involves constraining the number of axes into fewer dimensions, 
#through a PCA and appending these data to the pollinator visit data. 

####1.1 PREPARING DATA TO RUN PCA####
ph <- read.csv("Raw Data/plantphenology.csv", stringsAsFactors = TRUE, 
               strip.white = TRUE)

#Filter and select only the floral observations done during pollinator surveys
ph <- ph %>% filter(pollinator != "n")
str(ph)
ph$year <- as.factor(ph$year)

#to run the PCA, convert integer column types to numeric
ph <- ph %>% mutate_at(c('acmi', 'anlu', 'assu', 'caly', 'cami', 'capl', 'clsi', 
                         'copa', 'frca', 'gatr', 'glma', 'grass', 'juba', 'madi', 
                         'oesa', 'plmc', 'plmr', 'poan', 'sili', 'trar', 'trma', 
                         'trwo'), as.numeric)

#PCA can only run on columns with values it can calculate variance on. Must
#find columns with only 0-1 different values with minimal variance, and remove. 
lapply(ph, function(x) unique(x))
#remove the following plant species from the ordination due to minimal variance
ph <- ph %>% select(-madi, -oesa, -trar,) 

#The PCA will remove entire rows if some values are NA; There are 18 such rows due 
#to field logistical errors, where there were no inflorescence values. Thus, this
#segment of code calculates the column median and imputes the NA values to the 
#median in cases where the value was missing. 
f = function(x){
  x[is.na(x)] = median(x, na.rm=TRUE) #convert the NA to column median
  x #display the column
}
ph.na <- data.frame(apply(ph[,7:25],2,f))

#Since skew and the magnitude of the variables influence the resulting PCs, 
#it is good practice to apply skew transformation to variables prior to the 
#application of PCA. 
log.ph <- log(ph.na + 1)

####1.2 RUNNING THE PCA####
#This uses prior log transformed values to scale abundant plants, therefore do not
#scale the data. This fits a covariance matrix, and not a correlation matrix
log.ph.pca <- prcomp(na.omit(log.ph), center = TRUE, scale = FALSE)
biplot(log.ph.pca)
summary(log.ph.pca)
screeplot(log.ph.pca)
abline(a = 1, b = 0)
log.ph.pca$rotation

#PC1 and PC2 represent 68.73% of the variance in plant species, therefore I can 
#extract the values for PC1 and PC2, and append them as predictor variables
#for the number of insect visits to each plot. PC1 varies positively by yarrow 
#(ACMI) and negatively by common rush (JUBA) and silverweed (POAN). PC2 varies 
#positively by Douglas' aster (ASSU) and negatively by grasses. This aligns with
#field surveys of percent cover and anecdotal observations.

#ggplot2 will be used to customize the output of the PCA, so PCA scores will be 
#extracted and appended them to the data frame. Scores will then be used as 
#predictors to model the pollinator visit data. 

ph$pc1 <- log.ph.pca$x[, 1] #indexing the first column
ph$pc2 <- log.ph.pca$x[, 2] #indexing the second column
ph$pc3 <- log.ph.pca$x[, 3] #indexing the third column
ph$pc4 <- log.ph.pca$x[, 4] #indexing the fourth column

####1.3 PLOTTING THE PCA WITH GGPLOT2####
autoplot(log.ph.pca, data = ph,
         loadings = TRUE, loadings.colour = 'black', loadings.label = FALSE,
         loadings.label.size = 5, alpha = 0.15) +
  theme_classic(30)
ggsave("Figures/PCA_loadings.png",  height=9, width=16)
#note that this excludes the species labels that correspond to each loading for 
#ease of display in the publication (where the latin names were added). These 
#labels can be added back in with loadings.label = TRUE

####1.4 SPLICING EXTRACTED PC SCORES TO POLLINATOR DATA####
#Now that the PC1/PC2/PC3/PC4 scores for each plot have been extracted, and 
#these scores must be appended to the pollinator data set

#add unique identifiers to each data set row w/ so that data can be spliced
ph <- ph %>% mutate(ref.id = paste(date, year, pollinator, block,
                                   plot, sep = "_"))

#read in pollinator data set, summarize all visits by plot and, and add ref.id
poll <- read.csv("Raw Data/pollinator.csv", stringsAsFactors=TRUE, 
                 strip.white=TRUE)
poll_sum <- poll %>% group_by(date, year, pollinator, block, plot) %>% 
  summarize(sum.visits = sum(total.visits), bumblebees = sum(bumblebees),
            hoverflies = sum(hoverflies), other.flies = sum(other.flies)) 
str(poll_sum)

poll_sum <- poll_sum %>%  mutate(ref.id = paste(date, year, pollinator, block,
                                                plot, sep = "_"))

#unite plant and pollinator data by the unique identifiers, keeping plant 
#flowering stem numbers from phenology data set to use in SEM models 
joined <- left_join(poll_sum, ph[,7:30],  
                    by = c("ref.id"="ref.id"))

write.csv(joined, "Manipulated Data/summarized_pollinator.csv")


####2.0 FINAL POLLINATOR DATA LOADING & DATA CLEANING####
poll <- read.csv("Manipulated Data/summarized_pollinator.csv", 
                 stringsAsFactors=TRUE, 
                 strip.white=TRUE)

#Ensure all variables are numeric or factor when needed, adjust order of factors
#for ease of coeff interpretation
levels(poll[,'plot']) 
poll$plot <- factor(poll$plot,levels(poll$plot)[c(4,1,2,3)])
levels(poll[,'plot']) 
str(poll)
poll <- japply(poll, which(sapply(poll, class) == "integer"), as.numeric)
poll$year <- as.factor(poll$year)
poll$pollinator <- as.ordered(poll$pollinator)
poll$block <- as.factor(poll$block)

str(poll)

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-(1.96*std.error(x))
  ymax <- m+(1.96*std.error(x))
  return(c(y=m,ymin=ymin,ymax=ymax)) }

####3.0 CUMULATIVE VISIT BAR PLOTS####

#Code to creat stack bar charts of total visits
pollraw <- read.csv("Raw Data/pollinator.csv", stringsAsFactors=TRUE, 
                 strip.white=TRUE)

poll_sum_by_plant <- pollraw %>% 
  group_by(date, year, pollinator, block, plot, plant.species) %>% 
  gather("visitor.id", "visits", 14:21)

poll_sum_by_plant$visitor.id <- as.factor(poll_sum_by_plant$visitor.id)
levels(poll_sum_by_plant$visitor.id)

#rename factor levels for legend
poll_sum_by_plant$visitor.id <- recode_factor(poll_sum_by_plant$visitor.id,
                                              blowflies = "Blow Flies",
                                              bumblebees = "Bumblebees", 
                                              hoverflies = "Hover Flies",
                                              other.flies = "Other Flies",
                                              wasps = "Wasps",
                                              mason.bees = "Mason Bees",
                                              sweat.mining.bees = 
                                                "Sweat & Mining Bees",
                                              other = "Other")

#select only the most commonly visited plants and visitors to display in figure;
#remove rare visitors and plants (5 or fewer visits by insect, 17 or fewer visits
#to plant species). Note that there are 991 total visits by pollinators, yet only
#982 visits to plants. This is due to field observation recording error, which failed
#to assign 9 floral visits to a particular plant species. 

poll_sum_by_plant <- poll_sum_by_plant %>% filter(visitor.id %in% 
                                                    c("Bumblebees", "Hover Flies",
                                                      "Other Flies", "Wasps", 
                                                      "Other"))
levels(poll_sum_by_plant$visitor.id)

#create month labels
month_labels <- c(
  "1" = "June",
  "2" = "July",
  "3" = "August")

#create a colour pallette
pal <- pnw_palette("Cascades", 5)

#plot a bar chart of cumulative visits by plant species and year
ggplot(poll_sum_by_plant, aes(fill = visitor.id, y = visits, 
                              x = plant.species)) + 
  geom_bar(position = "stack", stat = "identity") +
  facet_grid(rows = vars(year), cols = NULL) +
  scale_x_discrete(limits = c("acmi", "assu", 
                              "cami", "grass", 
                              "juba", "poan"),
                   labels = c("A. millefolium", 
                              " S. subspicatum", "C. miniata",
                              "Grasses", "J. balticus", 
                              "P. anserina")) +
  scale_fill_manual(values = pal) +
  theme_classic(30) +
  theme(axis.text.x = element_text(angle = -45, vjust = -0.1, hjust = 0.1,
                                   face = c("italic","italic","italic",
                                            "italic", NULL, "italic")), 
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", size = 1, fill = NA),
        panel.spacing.x = unit(0.3,"line"),
        panel.spacing.y = unit(0.3,"line"),
        strip.background = element_rect(colour = "white",size = 1, fill = NA),
        line = element_line(size = 0.25),
        strip.text.y.right = element_text(angle = 0)) +
  labs(x = "Plant Species", y = "# Cumulative Visits", fill = "Floral Visitor")
ggsave("Figures/cumulativevisits_byplant_year.png",  height=9, width=16)

##plot a bar chart of cumulative visits by plant species and year
ggplot(poll_sum_by_plant, aes(fill = visitor.id, y = visits, 
                              x = plant.species)) + 
  geom_bar(position = "stack", stat = "identity") +
  facet_grid(rows = vars(year), cols = vars(pollinator), labeller = 
               labeller(pollinator = month_labels)) +
  scale_x_discrete(limits = c("acmi", "assu", 
                              "cami", "grass", 
                              "juba", "poan"),
                   labels = c("A. millefolium", 
                              " S. subspicatum", "C. miniata",
                              "Grasses", "J. balticus", 
                              "P. anserina")) +
  scale_fill_manual(values = pal) +
  theme_classic(30) +
  theme(axis.text.x = element_text(angle = -45, vjust = -0.1, hjust = 0.1,
                                   face = c("italic","italic","italic",
                                            "italic", NULL, "italic")), 
        axis.ticks = element_blank(), 
        panel.border = element_rect(colour = "black", size = 1, fill = NA),
        panel.spacing.x = unit(0.3,"line"),
        panel.spacing.y = unit(0.3,"line"),
        strip.background = element_rect(colour = "white",size = 1, fill = NA),
        line = element_line(size = 0.25),
        strip.text.y.right = element_text(angle = 0)) +
  labs(x = "Plant Species", y = "# Cumulative Visits", fill = "Floral Visitor")
ggsave("Figures/cumulativevisits_byplant_monthyear.png",  height=9, width=16)

#plot a bar chart of cumulative visits by treatment plot
ggplot(poll_sum_by_plant, aes(fill = visitor.id, y = visits, 
                              x = plot)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = pal) +
  theme_classic(30) +
  theme(axis.text.x = element_text(hjust = 0,
                                   vjust = -2.5),
        axis.title.x = element_text(margin = margin(t = 35))) +
  labs(x = "Treatment", y = "# Cumulative Visits", fill = "Floral Visitor")
ggsave("Figures/cumulativevisits_bytreatment.png",  height=9, width=16)

####4.0 TOTAL VISIT CALCULATIONS####
#Data summaries for calculations of proportion of visits by different groups to 
#different plants species and treatments to report in the results

#total visits to each plant species
p1 <- poll_sum_by_plant %>% group_by(plant.species) %>% summarise(sum(visits))
View(p1)

#total visits by insect and year
p2 <- poll_sum_by_plant %>% group_by(year, visitor.id) %>% summarise(sum(visits))
View(p2)

#total visits by treatment
p3 <- poll_sum_by_plant %>% group_by(plot) %>% summarise(sum(visits))
View(p3)



####5.0 STRUCTURAL EQUATION MODELLING DATA CLEANING####
#Recall poll data, which contains visits by plot, block, year, observation 
#period, PC1/PC2/PC3/PC4, and number of flowering stems by plant species

poll_SEM <- read.csv("Manipulated Data/summarized_pollinator.csv", 
                 stringsAsFactors = TRUE, 
                 strip.white = TRUE)
poll_SEM <- poll_SEM[,-1]

#Ensure all variables are numeric or factor when needed, adjust order of factors
#for ease of coefficient interpretation with "control" (plot D) as the intercept 
str(poll_SEM)
poll_SEM$year <- as.factor(poll_SEM$year)
poll_SEM$pollinator <- as.factor(poll_SEM$pollinator)
poll_SEM$block <- as.factor(poll_SEM$block)
levels(poll_SEM[,'plot']) 
poll_SEM$plot <- factor(poll_SEM$plot,levels(poll_SEM$plot)[c(4,1,2,3)])
levels(poll_SEM[,'plot']) 
poll_SEM$grass <- as.numeric(poll_SEM$grass)
poll_SEM$juba <- as.numeric(poll_SEM$juba)
poll_SEM$poan <- as.numeric(poll_SEM$poan)
poll_SEM$cami <- as.numeric(poll_SEM$cami)
poll_SEM$acmi <- as.numeric(poll_SEM$acmi)
poll_SEM$assu <- as.numeric(poll_SEM$assu)
poll_SEM$anlu <- as.numeric(poll_SEM$anlu)
poll_SEM$copa <- as.numeric(poll_SEM$copa)
poll_SEM$trma <- as.numeric(poll_SEM$trma)

##SEM will not run if NA values are present; There are 18 such rows due 
#to field logistical errors, where there were no inflorescence values. Thus, this
#segment of code calculates the column median and imputes the NA values to the 
#median in cases where the value was missing. This code is identical to that which
#was fed into the PCA for the same reason. 
poll_SEM <- poll_SEM %>% 
  mutate(grass = ifelse(is.na(grass), median(grass, na.rm = T), grass)) %>%
  mutate(juba = ifelse(is.na(juba), median(juba, na.rm = T), juba)) %>% 
  mutate(poan = ifelse(is.na(poan), median(poan, na.rm = T), poan)) %>% 
  mutate(cami = ifelse(is.na(cami), median(cami, na.rm = T), cami)) %>% 
  mutate(acmi = ifelse(is.na(acmi), median(acmi, na.rm = T), acmi)) %>% 
  mutate(assu = ifelse(is.na(assu), median(assu, na.rm = T), assu))


####5.1 SEM COMPONENT MODELS####
# Create component models of the SEM and store in list. There are 7 models as a
#part of the SEM: number of floral visits, grass stems, rush stems, silverweed
#stems, paintbrush steams, yarrow stems, and Douglas' aster stems

#visit model specification, diagnostics, and least-squares means 
visit_mod <- glmmTMB(sum.visits ~ plot + pollinator * year
                     + grass + juba + poan + cami + acmi + assu
                     + (1|block), data = poll_SEM, 
                     family = nbinom2(), ziformula = ~pollinator,
                     control = glmmTMBControl(optimizer = nlminb), na.action = "na.fail")
summary(visit_mod)
MuMIn::r.squaredGLMM(visit_mod)
res <- simulateResiduals(fittedModel = visit_mod, n = 500)
plot(res)
plotResiduals(res$scaledResiduals, poll_SEM$plot)
plotResiduals(res$scaledResiduals, poll_SEM$year)
plotResiduals(res$scaledResiduals, poll_SEM$pollinator)
plotResiduals(res$scaledResiduals, poll_SEM$grass)
plotResiduals(res$scaledResiduals, poll_SEM$juba)
plotResiduals(res$scaledResiduals, poll_SEM$poan)
plotResiduals(res$scaledResiduals, poll_SEM$cami)
plotResiduals(res$scaledResiduals, poll_SEM$acmi)
plotResiduals(res$scaledResiduals, poll_SEM$assu)
testOutliers(simulationOutput = res) 
testDispersion(res) 
testZeroInflation(res)
emsum <- emmeans(visit_mod, ~plot)
emsum
contrast(emsum, method = "dunnett", adjust = "dunnett")


#grass model specification, diagnostics, and least-squares means 
grass_mod <- glmmTMB(grass ~ plot + pollinator * year
                     + (1|block), data = poll_SEM, 
        family = nbinom2(), dispformula = ~year+pollinator, ziformula = ~year,
        control = glmmTMBControl(optimizer = nlminb), na.action = "na.fail")
summary(grass_mod)
MuMIn::r.squaredGLMM(grass_mod)
res <- simulateResiduals(fittedModel = grass_mod, n = 500)
plot(res)
plotResiduals(res$scaledResiduals, poll_SEM$plot)
plotResiduals(res$scaledResiduals, poll_SEM$year)
plotResiduals(res$scaledResiduals, poll_SEM$pollinator)
testOutliers(simulationOutput = res) 
testDispersion(res) 
testZeroInflation(res)
emsum <- emmeans(grass_mod, ~plot)
emsum
contrast(emsum, method = "dunnett", adjust = "dunnett")

#rush model specification, diagnostics, and least-squares means 
juba_mod <- glmmTMB(juba ~ plot + pollinator * year + (1|block), data = poll_SEM, 
                    family = nbinom2(), dispformula = ~year, ziformula = ~pollinator,
                    control = glmmTMBControl(optimizer = nlminb), na.action = "na.fail")
summary(juba_mod)
MuMIn::r.squaredGLMM(juba_mod)
res <- simulateResiduals(fittedModel = juba_mod, n = 500)
plot(res)
plotResiduals(res$scaledResiduals, poll_SEM$plot)
plotResiduals(res$scaledResiduals, poll_SEM$year)
plotResiduals(res$scaledResiduals, poll_SEM$pollinator)
testOutliers(simulationOutput = res) 
testDispersion(res) 
testZeroInflation(res)
emsum <- emmeans(juba_mod, ~plot)
emsum
contrast(emsum, method = "dunnett", adjust = "dunnett")

#silverweed model specification, diagnostics, and least-squares means 
poan_mod <- glmmTMB(poan ~ plot + pollinator * year + (1|block), data = poll_SEM, 
                    family = nbinom2(), dispformula = ~plot,
                    control = glmmTMBControl(optimizer = nlminb), na.action = "na.fail")
summary(poan_mod)
MuMIn::r.squaredGLMM(poan_mod)
res <- simulateResiduals(fittedModel = poan_mod, n = 500)
plot(res)
plotResiduals(res$scaledResiduals, poll_SEM$plot)
plotResiduals(res$scaledResiduals, poll_SEM$year)
plotResiduals(res$scaledResiduals, poll_SEM$pollinator)
testOutliers(simulationOutput = res) 
testDispersion(res) 
testZeroInflation(res)
emsum <- emmeans(poan_mod, ~plot)
emsum
contrast(emsum, method = "dunnett", adjust = "dunnett")

#paintbrush model specification, diagnostics, and least-squares means 
cami_mod <- glmmTMB(cami ~ plot + pollinator * year + (1|block), data = poll_SEM, 
                    family = nbinom2(),
                    control = glmmTMBControl(optimizer = nlminb), na.action = "na.fail")
summary(cami_mod)
MuMIn::r.squaredGLMM(cami_mod)
res <- simulateResiduals(fittedModel = cami_mod, n = 500)
plot(res)
plotResiduals(res$scaledResiduals, poll_SEM$plot)
plotResiduals(res$scaledResiduals, poll_SEM$year)
plotResiduals(res$scaledResiduals, poll_SEM$pollinator)
testOutliers(simulationOutput = res) 
testDispersion(res) 
testZeroInflation(res)
emsum <- emmeans(cami_mod, ~plot)
emsum
contrast(emsum, method = "dunnett", adjust = "dunnett")

#yarrow model specification, diagnostics, and least-squares means 
acmi_mod <- glmmTMB(acmi ~ plot + pollinator * year + (1|block), data = poll_SEM, 
                    family = nbinom2(),
                    control = glmmTMBControl(optimizer = nlminb), na.action = "na.fail")
summary(acmi_mod)
MuMIn::r.squaredGLMM(acmi_mod)
res <- simulateResiduals(fittedModel = acmi_mod, n = 500)
plot(res)
plotResiduals(res$scaledResiduals, poll_SEM$plot)
plotResiduals(res$scaledResiduals, poll_SEM$year)
plotResiduals(res$scaledResiduals, poll_SEM$pollinator)
testOutliers(simulationOutput = res) 
testDispersion(res) 
testZeroInflation(res)
emsum <- emmeans(acmi_mod, ~plot)
emsum
contrast(emsum, method = "dunnett", adjust = "dunnett")

#Douglas' aster model specification, diagnostics, and least-squares means 
assu_mod <- glmmTMB(assu ~ plot + pollinator * year + (1|block), data = poll_SEM, 
                    family = nbinom2(),
                    control = glmmTMBControl(optimizer = nlminb), na.action = "na.fail")
summary(assu_mod)
MuMIn::r.squaredGLMM(assu_mod)
res <- simulateResiduals(fittedModel = assu_mod, n = 500)
plot(res)
plotResiduals(res$scaledResiduals, poll_SEM$plot)
plotResiduals(res$scaledResiduals, poll_SEM$year)
plotResiduals(res$scaledResiduals, poll_SEM$pollinator)
testOutliers(simulationOutput = res) 
testDispersion(res) 
testZeroInflation(res)
emsum <- emmeans(assu_mod, ~plot)
emsum
contrast(emsum, method = "dunnett", adjust = "dunnett")

#create SEM list of component models
SEM_list <- list(visit_mod, grass_mod, juba_mod, poan_mod, cami_mod, acmi_mod, 
                 assu_mod)

####5.3 EVALUATE SEM####
#install older version of piecewiseSEM with support for glmmTMB
packageurl <- "https://cran.r-project.org/src/contrib/Archive/piecewiseSEM/piecewiseSEM_1.2.1.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
library(piecewiseSEM)

# Evaluate path significance using unstandardized coefficients
results <- piecewiseSEM::sem.coefs(SEM_list, poll_SEM, standardize = "none")
View(results)
