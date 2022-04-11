## ------------------------------------------------------------------------
## 'Biodiversity misuse'
## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
# 'R script to reproduce the full analysis'
## ------------------------------------------------------------------------

# Analysis performed with R (v. R 4.1.0) and R studio (v. 1.4.1103)
# Authors: Stefano Mammola

####################
# Data preparation #
####################

# Clean the workspace -----------------------------------------------------

rm(list = ls())

# Loading R package -------------------------------------------------------

library("Amelia")
library("dplyr")
library("gam")
library("ggplot2")
library("MuMIn")
library("parameters")
library("performance")
library("sjPlot")
library("PupillometryR")
library("tidyr")
library("tidytext")
library("tidyverse")
library("wordcloud")
library("RColorBrewer")

# Source functions and plot parameters ------------------------------------

source("Script/Functions.R")

# Loading database --------------------------------------------------------

db <- read.csv(file = "Data/Database_Full_study.csv", sep = '\t', dec='.', header = TRUE, as.is = FALSE)

str(db)
dim(db)

# Selecting paper to analyse
db <- db %>% filter(Analysis == "yes") ; db <- droplevels(db)

###########################
# Data preparation part 1 #
###########################

# Checking factor levels
levels(db$Geography)
levels(db$Domain)
levels(db$Method_data_collection)

# Calculating proportion of biodiversity
db$Biodiversity_prop <- rowSums(db[,36:91]) / length(36:91)
db$Animals_prop      <- rowSums(db[,c(36:70,89)]) / length(c(36:70,89))
db$Plants_prop       <- rowSums(db[,c(71:75,90)]) / length(c(71:75,90))
db$Fungi_prop        <- rowSums(db[,c(76:84,91)]) / length(c(76:84,91))
db$Micro_prop        <- rowSums(db[,c(85:88)])    / length(c(85:88))

#Title fanciness
db$Title_fanciness <- rowSums(db[,c(23,24)])
table(db$Title_fanciness) #too few obs

# Calculating total number of specifics to the title
db$Title_adjecties <- as.factor(rowSums(db[,27:29]))
table(db$Title_adjecties)

# Number of biodiversity facets
db$Facets_biodiversity <- as.factor(rowSums(db[,31:34]))
table(db$Facets_biodiversity)

# Missing data 
Amelia::missmap(db)

# Chcking N° of coutnris in the reference list

db$author_country <- as.character(db$author_country)

country_diversity <- c()

for(i in 1:nrow(db)){
  
  diversity_i <- strsplit(db[i,]$author_country, split = ";")[[1]]
  diversity_i <- trimws(diversity_i, which = c("both"))
  diversity_i <- unique(diversity_i)
  
  country_diversity <- append(country_diversity, length(diversity_i))}

db <- cbind(db,country_diversity) 

###########################
# Data preparation part 2 #
###########################

# Correcting citations by year of publication ----------------------

# What is the trend of citation over time?
db %>% ggplot(aes(x = Publication_year, y = tot_cites)) + 
  geom_point(size = 1, alpha = 0.7, color = "grey40")+
  labs(x = "year of publication", y = "citation")+
  theme_custom()

# Modeling temporal trend
M0 <- gam::gam(tot_cites ~ s(Publication_year), family = poisson, data = db) #constraining dof to 2 so that you don't get negative predicted citations
summary(M0)

par(mar= c(rep(2,4)))
plot(M0, se = TRUE)

# Taking the residual of citations
db <- data.frame(db, citation_residuals = resid(M0, type="pearson"))

# Checking how good it predict by comparing with a simple normalization by year of publication
db %>%  ggplot(aes(x = tot_cites/Publication_year, y = citation_residuals)) + #normalize by year
  geom_point(size = 1, alpha = 0.7, color = "grey40")+
  labs(x = "citation / year", y = "citation residuals")+
  theme_custom()

# Correcting altmetrics by year of publication ----------------------------

db_alt <- db %>% select(year = Publication_year,
                        Altmetrics = Altmetrics,
                        row_ID) %>% na.omit
                 
db_alt %>% ggplot(aes(x = year, y = Altmetrics)) + 
  geom_point(size = 1, alpha = 0.7, color = "grey40")+
  labs(x = "year of publication", y = "altmetric score")+
  theme_custom()

# Modeling temporal trend
M1 <- gam::gam(Altmetrics ~ s(year), family = poisson, data = db_alt) #constraining dof to 2 so that you don't get negative predicted citations
summary(M1)

dev.off()
par(mar= c(rep(2,4)))
plot(M1, se = TRUE)

# Taking the residual of citations & merge with the original db
db_alt <- data.frame(db_alt, Altmetrics_residuals = resid(M1,type="pearson")) # Taking the residual of citations

db_alt <- db_alt %>% select(row_ID,Altmetrics_residuals)

db <- db %>% dplyr::left_join(db_alt, by = "row_ID") ; rm(db_alt, M0, M1)

#####################################
# Summary stats  & Data exploration #
#####################################

# How many studies have no biodiversity proportion?
nrow(db[db$Biodiversity_prop == 0,])/nrow(db) * 100 #22.1% do not consider any biodiversity group

# Removing study with no Biodiversity 
db2 <- db %>% filter(Biodiversity_prop > 0)

# Biodiversity facets (%)
round((table(db2$Facets_biodiversity)/sum(table(db2$Facets_biodiversity)))*100,2)
table(db2$Phylogenetic_div)
table(db2$Functional_div)
table(db2$Taxonomic_div)
table(db2$Other_div)

#What proportion of biodiversity across studies?
mean(db2$Biodiversity_prop, na.rm = TRUE) #mean 
std(db2$Biodiversity_prop) #sd 
range(db2$Biodiversity_prop, na.rm = TRUE) #range

# Checking temporal distribution
plot1a <- db2 %>% ggplot(aes(x = Publication_year, y = Biodiversity_prop)) + 
  geom_point(col = "grey10", fill = "grey30", size = 5, shape = 21, alpha = 0.3) +
  geom_smooth(method = "glm", formula = y ~ x, 
              method.args = list(family = quasibinomial(link = "logit")),
              col="blue", fill = "blue") +
  labs(x = NULL, y = Y.label)+
  theme_custom()

#Split multiple regions separated by ";"
box1 <- semi_colon_splitter(input1 = db2$Geography,
                            input2 = db2$Biodiversity_prop, 
                            names  = c("Geography","Biodiversity_prop"))

box1$Biodiversity_prop <- as.numeric(as.character(box1$Biodiversity_prop))

#Sort levels
box1$Geography <- factor(box1$Geography,
                         c(levels(box1$Geography)[4],levels(box1$Geography)[-4]))

plot1b <- box1 %>% 
    drop_na(Geography,Biodiversity_prop) %>% 
    filter(Biodiversity_prop < 0.4) %>%  #removing 1 outlier
    ggplot(aes(x = Geography, y = Biodiversity_prop)) +
    geom_flat_violin(position = position_nudge(x = 0.2, y = 0), 
                     alpha = 0.4, fill= "blue", col = "white", adjust = 1.5) +
    geom_point(position = position_jitter(width = 0.15), size = 1, alpha = 0.7, color = "grey40") +
    geom_boxplot(width = 0.2,  col = "blue", outlier.shape = NA, alpha = 0) +
    labs(y = Y.label, x = NULL) +
    theme_custom()

#Split multiple domains separated by ";"
box2 <- semi_colon_splitter(input1 = db2$Domain,
                            input2 = db2$Biodiversity_prop, 
                            names = c("Domain","Biodiversity_prop"))

box2$Biodiversity_prop <- as.numeric(as.character(box2$Biodiversity_prop))

plot1c <- box2 %>% 
    drop_na(Domain,Biodiversity_prop) %>% 
    ggplot(aes(x = Domain, y = Biodiversity_prop)) +
    geom_flat_violin(position = position_nudge(x = 0.2, y = 0), 
                     alpha = 0.4, fill= "blue", col = "white", adjust = 1.5) +
    geom_point(position = position_jitter(width = 0.15), size = 1, alpha = 0.7, color = "grey40") +
    geom_boxplot(width = 0.2, col = "blue", outlier.shape = NA, alpha = 0) +
    labs(y = Y.label, x = NULL) +
    theme_custom()

#Split multiple methds separated by ";"
box3 <- semi_colon_splitter(input1 = db2$Method_data_collection,
                            input2 = db2$Biodiversity_prop, 
                            names = c("Method","Biodiversity_prop"))

box3$Biodiversity_prop <- as.numeric(as.character(box3$Biodiversity_prop))

#Sort levels
box3$Method <- factor(box3$Method,
                         c(levels(box3$Method)[-4],levels(box3$Method)[4]))

levels(box3$Method)[c(2,5)] <- "Other"

plot1d <- box3 %>% 
    drop_na(Method,Biodiversity_prop) %>% 
    ggplot(aes(x = Method, y = Biodiversity_prop)) +
    geom_flat_violin(position = position_nudge(x = 0.2, y = 0), 
                     alpha = 0.4, fill= "blue", col = "white", adjust = 1.5) +
    geom_point(position = position_jitter(width = 0.15), size = 1, alpha = 0.7, color = "grey40") +
    geom_boxplot(width = 0.2, col = "blue", outlier.shape = NA, alpha = 0) +
    labs(y = Y.label, x = NULL) +
    theme_custom()

pdf(file = "Figure/Figure_2.pdf", width = 19, height = 14)

ggpubr::ggarrange(plot1a, plot1b, plot1c, plot1d,
                  common.legend = FALSE,
                  hjust = -5,
                  align = "hv",
                  labels = c("A", "B", "C", "D"),
                  ncol = 2, nrow = 2) #warnings() are due to NA removal

dev.off()

#clean
rm(box1,box2,box3)
rm(plot1a, plot1b, plot1c, plot1d)

# Figure S1 (by group) --------------------------------------------------

#Loading silhouettes
animal_png <- png::readPNG("Phylopics/Animal.png")
fungi_png  <- png::readPNG("Phylopics/Fungi.png")
micro_png  <- png::readPNG("Phylopics/Micro.png")
plant_png  <- png::readPNG("Phylopics/Plant.png")

# Plotting
plot2a <- ggplot(data = db[db$Animals_prop>0,], aes(x = Publication_year, y = Animals_prop)) + 
    geom_point(col = "grey10", fill = "grey30", size = 5, shape = 21, alpha = 0.3) +
    geom_smooth(method = "glm", formula = y ~ x, 
                method.args = list(family = quasibinomial(link = "logit")), 
                col = "blue", fill = "blue") +
    labs(title = "Only animals", x = NULL , y = Y.label)+
    xlim(1992,2020)+
    ylim(0,1)+
    annotation_custom(grid::rasterGrob(animal_png), xmin = 1990, xmax = 1998, ymin = 0.75, ymax = 1)+ 
    theme_custom()

plot2b <- ggplot(db[db$Plants_prop>0,], aes(x = Publication_year, y = Plants_prop)) + 
    geom_point(col = "grey10", fill = "grey30", size = 5, shape = 21, alpha = 0.3) +
    geom_smooth(method = "glm", formula = y ~ x, 
                method.args = list(family = quasibinomial(link = "logit")), 
                col = "blue", fill = "blue") +
    labs(title = "Only plants", x = NULL , y = NULL)+
    xlim(1992,2020)+
    ylim(0,1)+
    annotation_custom(grid::rasterGrob(plant_png), xmin = 1992, xmax = 1995, ymin = 0.75, ymax = 1)+ 
    theme_custom()

plot2c <- ggplot(data = db[db$Fungi_prop>0,], aes(x = Publication_year, y = Fungi_prop)) + 
    geom_point(col = "grey10", fill = "grey30", size = 5, shape = 21, alpha = 0.3) +
    geom_smooth(method = "glm", formula = y ~ x, 
                method.args = list(family = quasibinomial(link = "logit")), 
                col = "blue", fill = "blue") +
    labs(title = "Only fungi", x = "Publication year" , y = Y.label)+
    xlim(1992,2020)+
    ylim(0,1)+
    annotation_custom(grid::rasterGrob(fungi_png), xmin = 1992, xmax = 1997, ymin = 0.75, ymax =1)+ 
    theme_custom()

plot2d <- ggplot(data = db[db$Micro_prop>0,], aes(x = Publication_year, y = Micro_prop)) + 
    geom_point(col = "grey10", fill = "grey30", size = 5, shape = 21, alpha = 0.3) +
    geom_smooth(method = "glm", formula = y ~ x, 
                method.args = list(family = quasibinomial(link = "logit")), 
                col = "blue", fill = "blue") +
    labs(title = "Only microorganisms", x = "Publication year" , y = NULL)+
    xlim(1992,2020)+
    ylim(0,1)+
    annotation_custom(grid::rasterGrob(micro_png), xmin = 1992, xmax = 1997, ymin = 0.75, ymax = 1)+
    theme_custom()

pdf(file = "Figure/Figure_S1.pdf", width = 16, height =12)
ggpubr::ggarrange(plot2a,plot2b,plot2c,plot2d,
                  common.legend = FALSE,
                  hjust = -5,
                  align = "hv",
                  labels = c("A", "B", "C", "D"),
                  ncol=2, nrow=2) #warnings() are due to NA removal
dev.off()

#clean
rm(plot2a,plot2b,plot2c,plot2d)
rm(animal_png, fungi_png, micro_png, plant_png)

# End of data exploration

#####################################################
# Data preparation: part 2 (after data exploration) #
#####################################################

# Converting multiples levels in factors
method_split <- strsplit(as.character(db2$Method_data_collection), ";")

method <- c()
for(i in 1:length(method_split))
  method <- c(method, ifelse(length(method_split[[i]]) > 1, "Multiple", method_split[[i]]) )

geography_split <- strsplit(as.character(db2$Geography), ";")

geography <- c()
for(i in 1:length(geography_split))
  geography <- c(geography, ifelse(length(geography_split[[i]]) > 1, "Global", geography_split[[i]]) )

domain_split <- strsplit(as.character(db2$Domain), ";")

domain <- c()
for(i in 1:length(domain_split))
  domain <- c(domain, ifelse(length(domain_split[[i]]) > 1, "Multiple", domain_split[[i]]) )

db2$Method    <- method
db2$Geography <- geography
db2$Domain    <- domain

db2 <- db2 %>% mutate_at(vars("Method","Geography","Domain","Title_geo","Title_hab","Title_taxon"), as_factor)

# Set baseline
db2 <- within(db2, Geography <- relevel(Geography, ref = "Global"))
db2 <- within(db2, Domain    <- relevel(Domain,    ref = "Multiple"))
db2 <- within(db2, Method    <- relevel(Method,    ref = "Multiple"))

#clean
rm(method_split, method, geography, geography_split, domain, domain_split, i)

########################################################################
# Regression model (all data) ------------------------------------------
########################################################################

# Data exploration --------------------------------------------------------

# Subset
db_glm <- db2 %>% select(year = Publication_year,
                        n_aut,
                        citation_residuals,
                        Method,
                        Phylogenetic_div,      
                        Functional_div,
                        Other_div,
                        Title_geo,
                        Title_hab,
                        Title_taxon,
                        Title_adjecties,
                        Biodiversity_prop,
                        Geography,
                        Domain,
                        country_diversity)

# Checking balancing of factors
table(db_glm$Method) #Citizen science/simulation too few records
levels(db_glm$Method)[c(6,7)] <- "Other"

table(db_glm$Geography) #Antartica a bit unbalanced

table(db_glm$Domain) #OK

db_glm <- droplevels(db_glm)

# Dependent var
db_glm$prop  <- rowSums(db2[,36:91])
db_glm$total <- length(36:91)

# Checking outliers
# par(mar= c(rep(2,4)))
# dotchart(db_glm$citation_residuals) # 1 outlier
# dotchart(db_glm$prop) #2 outliers

db_glm <- db_glm[db_glm$citation_residuals < 200,]
db_glm <- db_glm[db_glm$prop < 20,]

# Scale continuous variables
db_glm <- db_glm %>% mutate_at(vars(year,citation_residuals,country_diversity),scale)

# Modelling  --------------------------------------------------------------

# Set formula
model_1 <- as.formula("prop ~ 
                      year + 
                      Domain + Geography + Method +
                      Phylogenetic_div + Functional_div + Other_div +
                      Title_geo + Title_hab + Title_taxon + country_diversity")

#Initial model
m1  <- glm(model_1, data = db_glm, family = "poisson")
performance::check_overdispersion(m1) #overdispersed

#Refit with quasibinomial due to overdispersion 
m1b  <- MASS::glm.nb(model_1, data = db_glm)

parameters::parameters(m1b, df_method="wald")
performance::check_collinearity(m1b)
performance::r2(m1b)

# PostHoc
pairs(emmeans::emmeans(m1b, ~ Domain), simple=c("Domain"))
pairs(emmeans::emmeans(m1b, ~ Geography), simple="Geography")
pairs(emmeans::emmeans(m1b, ~ Method), simple="Method")

# Estract estimates
Estimates_m1 <- 
  m1b %>% 
  summary %>% 
  magrittr::extract2("coefficients") %>% # extract estimates
  as.data.frame %>% rownames_to_column("Variable") %>% 
  dplyr::filter(!row_number() %in% 1) %>%  #remove intercept
  dplyr::rename(SE = 3, z = 4, p = 5) #rename

# Set variable order and rename
order_var1 <- c("Year of publication",
               "Domain [Terrestrial]",
               "Domain [Saltwater]",
               "Domain [Freshwater]",
               "Geographic [Palearctic]",
               "Geographic [Afrotropical]",
               "Geographic [Indomalayan]",
               "Geographic [Neartic]",
               "Geographic [Australasian]",
               "Geographic [Neotropical]",
               "Geographic [Antartic]",
               "Method [Review/Opinion]",
               "Method [Field sampling]",
               "Method [Big data]",
               "Method [Other]",
               "Phylogenetic diversity [yes]",
               "Functional diversity [yes]",
               "Other diversity [yes]",
               "Mention of location in title [yes]",
               "Mention of habitat in title [yes]",
               "Mention of taxon/a in title [yes]",
               "Number of coauthors' countries")

Estimates_m1$Variable <- order_var1 #Rename
Estimates_m1$Variable <- factor(Estimates_m1$Variable, rev(order_var1)) #Sort

sign <- ifelse(Estimates_m1$p > 0.05, "", ifelse(Estimates_m1$Estimate>0.01," *", " **")) #Significance
col_p <- ifelse(Estimates_m1$p > 0.05, "grey5", ifelse(Estimates_m1$Estimate>0,"orange","blue")) #Significance

# Plot
plot_model1 <- ggplot2::ggplot(data = Estimates_m1) +

    geom_pointrange(aes(x = Variable, 
                        y = Estimate,
                        ymin = Estimate-SE, 
                        ymax = Estimate+SE), col = col_p, size = 0.5) + 
    
    geom_hline(lty = 3, size = 0.7, col = "grey50", yintercept = 0) +
    
    geom_text(aes(Variable, Estimate),
              label = paste0(round(Estimates_m1$Estimate,2),sign), 
              vjust = -1, size = 3, col = col_p) +
   
    labs(title = paste0("Articles with biodiversity proportion > 0 [N = ",nrow(db_glm),"]"),
         y = expression(paste("Odds ratio" %+-% "Standard Error")),
         x = NULL)+
    theme_custom() + theme(axis.text.y  = element_text(colour = rev(col_p))) + coord_flip()+
  annotate(geom = 'text', x = 2, y = -0.9, size =5,
           label = paste0("R^2 ==",round(as.numeric(performance::r2(m1b)[1]),2)), parse = TRUE)

# Clean
rm(m1, m1b, sign, order_var1, model_1, col_p)

########################################################################
# Regression model (partial data) --------------------------------------
########################################################################

# Repeating the analysis only with generic titles (no moderators)
db_glm2 <- db2 %>% drop_na(Title_adjecties) %>% 
           filter(Title_adjecties == 0) %>% 
           select(year = Publication_year,
           n_aut,
           citation_residuals,
           Method,
           Phylogenetic_div,      
           Functional_div,
           Other_div,
           Biodiversity_prop,
           Geography,
           Domain,
           country_diversity) 

db_glm2 <- na.omit(db_glm2)

# Balance of factors

table(db_glm2$Method) #Problems!
levels(db_glm2$Method)[c(6,7)] <- "Other"

table(db_glm2$Geography) #problems!

db_glm2 <- db_glm2 %>% filter(!Geography == "Antarctic")
db_glm2$Geography <- droplevels(db_glm2$Geography)

table(db_glm2$Geography) #problems!

levels(db_glm2$Geography)[4] <-  "Afrotropical"
db_glm2$Geography <- droplevels(db_glm2$Geography)

# db_glm2$Geography2 <- db_glm2$Geography
# levels(db_glm2$Geography2)[c(2:7)] <-  c("Global North", "Global South", "Global South", "Global North", "Global North", "Global South")
# db_glm2$Geography2 <- droplevels(db_glm2$Geography2)

table(db_glm2$Domain)  #problems!

levels(db_glm2$Domain)[c(3,4)] <- "Aquatic"
db_glm2$Domain <- droplevels(db_glm2$Domain)

# Dependent var
db_glm2$prop  <- db_glm2$Biodiversity_prop*56
db_glm2$total <- length(36:91)

# Checking outliers
# par(mar= c(rep(2,4)))
# dotchart(db_glm2$citation_residuals) # OK
# dotchart(db_glm2$prop) #1 outliers

# Removing outliers
db_glm2 <- db_glm2 %>% filter(prop < 11)

# Scale continuous variables
db_glm2 <- db_glm2 %>% mutate_at(vars(year,citation_residuals,country_diversity),scale)

# Modelling  --------------------------------------------------------------

# Set formula
model_2 <- as.formula("prop ~ 
                      year + 
                      Domain + 
                      Geography + 
                      Method +
                      Phylogenetic_div + 
                      Functional_div +
                      Other_div +
                      country_diversity")

# Initial model
m2  <- glm(model_2, data = db_glm2, family = "poisson")

performance::check_overdispersion(m2) 
parameters::parameters(m2, df_method="wald")
performance::check_collinearity(m2)
performance::r2(m2)

# PostHoc
pairs(emmeans::emmeans(m2, ~ Domain), simple=c("Domain"))
pairs(emmeans::emmeans(m2, ~ Geography), simple="Geography")
pairs(emmeans::emmeans(m2, ~ Method), simple="Method")

# Estract estimates
Estimates_m2 <- 
  m2 %>% 
  summary %>% 
  magrittr::extract2("coefficients") %>% # extract estimates
  as.data.frame %>% rownames_to_column("Variable") %>% 
  dplyr::filter(!row_number() %in% 1) %>%  #remove intercept
  dplyr::rename(SE = 3, z = 4, p = 5) #rename

# Set variable order and rename
order_var2 <- c("Year of publication",
                "Domain [Terrestrial]",
                "Domain [Aquatic]",
                "Geographic [Palearctic]",
                "Geographic [Afrotropical]",
                "Geographic [Neartic]",
                "Geographic [Australasian]",
                "Geographic [Neotropical]",
                "Method [Review/Opinion]",
                "Method [Field sampling]",
                "Method [Big data]",
                "Method [Other]",
                "Phylogenetic diversity [yes]",
                "Functional diversity [yes]",
                "Other diversity [yes]",
                "Number of coauthors' countries")

Estimates_m2$Variable <- order_var2 #Rename
Estimates_m2$Variable <- factor(Estimates_m2$Variable, rev(order_var2)) #sort

sign <- ifelse(Estimates_m2$p > 0.05, "", ifelse(Estimates_m2$Estimate>0.01," *", " **")) #Significance
col_p <- ifelse(Estimates_m2$p > 0.05, "grey5", ifelse(Estimates_m2$Estimate>0,"orange","blue") )

plot_model2 <- ggplot2::ggplot(data = Estimates_m2) +
    
    geom_pointrange(aes(x = Variable, 
                        y = Estimate,
                        ymin = Estimate-SE, 
                        ymax = Estimate+SE), col = col_p, size = 0.5) + 
    
    geom_hline(lty = 3, size = 0.7, col = "grey50", yintercept = 0) +
    
    geom_text(aes(Variable, Estimate),
              label = paste0(round(Estimates_m2$Estimate,2),sign), 
              vjust = -1, size = 3, col = col_p) +
    
    labs(title = paste0("Articles with no moderators in the title [N = ",nrow(db_glm2),"]"),
         y = expression(paste("Odds ratio" %+-% "Standard Error")),
         x = NULL)+
    theme_custom() + theme(axis.text.y  = element_text(colour = rev(col_p))) + coord_flip()+
  annotate(geom = 'text', x = 1.2, y = -0.7, size =5,
           label = paste0("R^2 ==",round(as.numeric(performance::r2(m2)[1]),2)), parse = TRUE)


rm(m2, sign, order_var2, model_2, col_p)

# Figure 4 ----------------------------------------------------------------

vector <- sort(apply(db2[,36:86],2, sum, na.rm = TRUE), decreasing = TRUE)

bar1 <- data.frame(Phyla = names(vector), N = vector) 

threshold <- 10 #cut of to merge in multiple category
Taxa_to_rename <- ifelse(bar1$N < threshold, bar1$N, NA)
N_to_rename    <- sum(Taxa_to_rename,na.rm = TRUE)
Taxa_to_rename <- length(na.omit(Taxa_to_rename))

bar1 <- bar1[ 1 : (nrow(bar1)-Taxa_to_rename), ]

bar1 <- rbind(bar1, data.frame(Phyla = paste0("Others (n = ",Taxa_to_rename,")"), N = N_to_rename))

bar1 <- cbind(bar1, Type = c(rep("Animal",2),
                           "Plant",
                           rep("Animal",2),
                           "Plant",
                           "Microorganism",
                           "Plant",
                           rep("Animal",3),
                           rep("Plant",2),
                           "Animal",
                           "Fungi", #ascoo
                           rep("Animal",2),
                           "Multiple",
                           "Fungi",
                           "Multiple"))

bar1$Phyla <- factor(bar1$Phyla,levels = bar1$Phyla)
levels(bar1$Phyla)[18] <- "Protista *"

col_bar <- rev(RColorBrewer::brewer.pal(5, "Blues"))

barchart1 <- ggplot(bar1, aes(x = Phyla, y = N))+
    geom_bar(aes(fill= Type),stat="identity", alpha=1, colour = "black")  +
    labs(x = NULL, y = "Count")+
    scale_fill_manual(values = col_bar)+
    theme_custom()+ 
    theme(legend.position = "none")+ 
    coord_flip()+
    annotate(geom = "text", y = bar1[nrow(bar1),2] + 2, 
             x = nrow(bar1), 
             label = paste0("Phyla/division with sample size < ", threshold),
             hjust = 0)

rm(Taxa_to_rename, N_to_rename, vector, bar1, threshold)

## Only with titles containong no moderators

bar2 <- db2[db2$Title_adjecties == 0,]

vector <- sort(apply(bar2[,36:86],2, sum, na.rm = TRUE), decreasing = TRUE)

bar2 <- data.frame(Phyla = names(vector), N = vector)

threshold2 <- 3 #cut of to merge in multiple category
Taxa_to_rename <- ifelse(bar2$N < threshold2, bar2$N, NA)
N_to_rename    <- sum(Taxa_to_rename,na.rm=TRUE)
Taxa_to_rename <- length(na.omit(Taxa_to_rename))

bar2 <- bar2[ 1 : (nrow(bar2)-Taxa_to_rename), ]

bar2 <- rbind(bar2, data.frame(Phyla = paste0("Others (n = ",Taxa_to_rename,")"), N = N_to_rename))

bar2 <- cbind(bar2, Type = c(rep("Animal",2),
                             rep("Plant",2),
                             "Microorganism",
                             rep("Animal",3),
                             rep("Plant",2),
                             rep("Animal",2),
                             "Multiple"))

bar2$Phyla <- factor(bar2$Phyla, levels = bar2$Phyla)

col_bar2 <- col_bar[c(1,3:5)]

barchart2 <- ggplot(bar2, aes(x = Phyla, y = N))+
    geom_bar(aes(fill = Type),stat="identity", alpha=1, colour = "black")  +
    labs(x = NULL, y = "Count")+
    scale_fill_manual(values = col_bar2) +
    theme_custom() + 
    theme(legend.position = c(0.75,0.45)) + 
    coord_flip() +
    annotate(geom = "text", y = bar2[nrow(bar2),2] + 1, 
             x = nrow(bar2), 
             label = paste0("Phyla/division with sample size < ",threshold2),
             hjust = 0)

# Save
pdf(file = "Figure/Figure_3.pdf", width = 15, height = 16)

ggpubr::ggarrange(plot_model1,barchart1,plot_model2,barchart2,
                  common.legend = FALSE,
                  #hjust = -5,
                  align = "h",
                  labels = c("A", "B", "C", "D"),
                  ncol=2, nrow=2) #warnings() are due to NA removal
dev.off()

rm(Taxa_to_rename, N_to_rename, vector, col_bar, col_bar2, bar2, threshold2)

#############################
# Modelling article impact ##
#############################

#missing impact factor

# Preparing the data ------------------------------------------------------

db2$bio <- db2$Biodiversity_prop*56

db_impact <- db2 %>% select(year = Publication_year,
                        Altmetrics,
                        cit = tot_cites,
                        Altmetrics_residuals,
                        citation_residuals,
                        Biodiversity = bio,
                        Moderators = Title_adjecties,
                        Title_geo,
                        Title_taxon,              
                        Title_hab,
                        n_aut,
                        journal,
                        country_diversity) %>% 
                        mutate_at(vars(starts_with("Title_"), journal), as_factor) %>% 
                        mutate_at(vars(n_aut, year, country_diversity), scale)

# Setting baseline
levels(db_impact$Title_taxon) <- c("Mention","No mention")
levels(db_impact$Title_hab)   <- c("Mention","No mention")
levels(db_impact$Title_geo)   <- c("Mention","No mention")
levels(db_impact$Moderators)  <- c("None","One moderator", "Two moderators", "Three moderators")

db_impact <- within(db_impact, Title_geo   <- relevel(Title_geo, ref = "Mention"))
db_impact <- within(db_impact, Title_taxon <- relevel(Title_taxon, ref = "Mention"))
db_impact <- within(db_impact, Title_hab   <- relevel(Title_hab, ref = "Mention"))

#Random factor levels
nlevels(db_impact$journal) 

# Checking outliers
#dotchart(db_impact$Biodiversity) #1 outliers
db_impact$Biodiversity <- scale(log(db_impact$Biodiversity+1)) #log transform

# N° of authors and country diversity are collinear
cor(db_impact$n_aut, db_impact$country_diversity) 

# Testing Citations -----------------------------------------------------------

#Check outliers
par(mar= c(rep(2,4)))
#dotchart(db_impact$citation_residuals)  #1 outlier
db_cit <- db_impact %>% filter(citation_residuals < 150) #removing outlier

# First model (all moderators sum,ed)

model_3 <- as.formula("citation_residuals ~ Title_geo + 
                                            Title_taxon + 
                                            Title_hab + 
                                            country_diversity +
                                            Biodiversity : Moderators +
                                            (1|journal)")
#Biodiversity + 

m3  <- lme4::lmer(model_3, data = db_impact)
parameters::parameters(m3)
check_collinearity(m3)
performance::r2(m3)[1]

# Estract estimates
Estimates_m3 <- 
  m3 %>% 
  summary %>% 
  magrittr::extract2("coefficients") %>% # extract estimates
  as.data.frame %>% rownames_to_column("Variable") %>% 
  dplyr::filter(!row_number() %in% 1) %>%  #remove intercept
  dplyr::rename(SE = 3, t = 4) #rename

# Set variable order and rename
order_var3 <- c("Mention of location in title [no]",
                "Mention of taxon/a in title [no]",
                "Mention of habitat in title [no]",
                "Number of coauthors' countries",
                "Sampled biodiversity * No moderators",
                "Sampled biodiversity * 1 moderator",
                "Sampled biodiversity * 2 moderators", 
                "Sampled biodiversity * 3 moderators")

Estimates_m3$Variable <- order_var3 #Rename
Estimates_m3$Variable <- factor(Estimates_m3$Variable, rev(order_var3)) #Sort

par <- parameters::parameters(m3) %>%
  dplyr::filter(!row_number() %in% 1)

par <- na.omit(par$p)

sign <- ifelse(par > 0.05, "", ifelse(par > 0.01," *", " **")) #Significance
col_p <- ifelse(par > 0.05, "grey5", "blue")

#Set variable order and rename

(plot_model3 <- ggplot2::ggplot(data = Estimates_m3) +
    
    geom_pointrange(aes(x = Variable, 
                        y = Estimate,
                        ymin = Estimate-SE, 
                        ymax = Estimate+SE), col = col_p, size = 0.5) + 
    
    geom_hline(lty = 3, size = 0.7, col = "grey50", yintercept = 0) +
    geom_text(aes(Variable, Estimate), label = paste0(round(Estimates_m3$Estimate,2),sign), 
              vjust = -1, size = 3, col = col_p) +
    
    labs(title = paste0("Citations [N = ", nrow(db_impact),"]"),
         y = expression(paste("Estimate beta" %+-% "Standard Error")),
         x = NULL)+
    theme_custom() + theme(axis.text.y  = element_text(colour = rev(col_p))) + coord_flip() +
    annotate(geom = 'text', x = 1, y = 3.5, size = 5,
             label = paste0("R^2 ==",round(as.numeric(performance::r2(m3)[1]),2)), parse = TRUE)
)

# Check the interaction
(plot_interaction1 <- db_cit %>% ggplot2::ggplot(aes(x = Biodiversity, y = citation_residuals)) + 
    facet_wrap( ~ Moderators, nrow = 2, ncol = 2) +
    geom_point(col = "grey10", fill = "grey30", size = 5, shape = 21, alpha = 0.3)+
    geom_smooth(method = "lm",  se = TRUE, col = "blue", fill = "blue",
                formula = y ~ x) +
    labs(x = "", 
         y = "Citations [residuals]",
         title = "Interaction sampled biodiversity * N° of moderators")+ theme_custom())

# Testing altmetric -------------------------------------------------------

db_alt <- na.omit(db_impact)

# Checking outliers
# par(mar= c(rep(2,4)))
# dotchart(db_alt$Altmetrics_residuals) # OK

# Set formula
model_4 <- as.formula("Altmetrics_residuals ~ Title_geo + 
                                            Title_taxon + 
                                            Title_hab + 
                                            country_diversity +
                                            Biodiversity : Moderators +
                                            (1|journal)")

# Initial model
m4 <- lme4::lmer(model_4, data = db_alt)
performance::r2(m4)
parameters::parameters(m4)

# Estract estimates
Estimates_m4 <- 
  m4 %>% 
  summary %>% 
  magrittr::extract2("coefficients") %>% # extract estimates
  as.data.frame %>% rownames_to_column("Variable") %>% 
  dplyr::filter(!row_number() %in% 1) %>%  #remove intercept
  dplyr::rename(SE = 3, t = 4) #rename

# Set variable order and rename
order_var4 <- c("Mention of location in title [no]",
                "Mention of taxon/a in title [no]",
                "Mention of habitat in title [no]",
                "Number of coauthors' countries",
                "Sampled biodiversity * No moderators",
                "Sampled biodiversity * 1 moderator",
                "Sampled biodiversity * 2 moderators", 
                "Sampled biodiversity * 3 moderators")

Estimates_m4$Variable <- order_var4 #Rename
Estimates_m4$Variable <- factor(Estimates_m4$Variable, rev(order_var4)) #Sort

par <- parameters::parameters(m4) %>%
  dplyr::filter(!row_number() %in% 1)

par <- na.omit(par$p)

sign <- ifelse(par > 0.05, "", ifelse(par > 0.01," *", " **")) #Significance
col_p <- ifelse(par > 0.05, "grey5", "blue")

(plot_model4 <- ggplot2::ggplot(data = Estimates_m4) +
    
    geom_pointrange(aes(x = Variable, 
                        y = Estimate,
                        ymin = Estimate-SE, 
                        ymax = Estimate+SE), col = col_p, size = 0.5) + 
    
    geom_hline(lty = 3, size = 0.7, col = "grey50", yintercept = 0) +
    
    geom_text(aes(Variable, Estimate), label = paste0(round(Estimates_m4$Estimate,2),sign), 
              vjust = -1, size = 3, col = col_p) +
    
    labs(title = paste0("Altmetric score [N = ", nrow(db_alt),"]"),
         y = expression(paste("Estimate beta" %+-% "Standard Error")),
         x = NULL)+
    theme_custom() + theme(axis.text.y  = element_text(colour = rev(col_p))) + coord_flip()+
    annotate(geom = 'text', x = 1, y = 4.5, size = 5,
             label = paste0("R^2 ==",round(as.numeric(performance::r2(m4)[1]),2)), parse = TRUE)
)

# Check the interaction
(plot_interaction2 <- db_alt %>% ggplot2::ggplot(aes(x = Biodiversity, y = Altmetrics_residuals)) + 
    facet_wrap( ~ Moderators, nrow = 2, ncol = 2) +
    geom_point(col = "grey10", fill = "grey30", size = 5, shape = 21, alpha = 0.3)+
    geom_smooth(method = "lm",  se = TRUE, col = "blue", fill = "blue",
                formula = y ~ x) +
    labs(x = "Sampled biodiversity [log-transformed]", 
         y = "Altmetric score [residuals]",
         title = NULL)+ theme_custom())

pdf(file = "Figure/Figure_4.pdf", width = 16, height = 14)

ggpubr::ggarrange(plot_model3,plot_interaction1,plot_model4,plot_interaction2,
                  common.legend = FALSE,
                  #hjust = -5,
                  align = "h",
                  labels = c("A", "B", "C", "D"),
                  ncol=2, nrow=2) #warnings() are due to NA removal
dev.off()

rm(col_p, sign, par)

# Map ---------------------------------------------------------------------

# Loading data
world <- map_data("world")
biog_regions <- raster::shapefile("Shapefiles/biogeographic_regions.shp")

(map1 <- ggplot() +
  geom_map(map = world, data = world,
           aes(map_id = region), 
           color = "gray10", fill = "gray10", size = 0.3) +
  
  labs(title = NULL) +
  
  #Add bioregion
  geom_path(data = fortify(biog_regions),
            aes(x = long, y = lat, group = group),
            color = 'blue', size = .2) +
  # 
  # #Antartica
  # annotate(geom="text", x=120, y=-80, label="Antartic",
  #          color="grey10")+
  # 
  # #Antartica
  # annotate(geom="text", x=50, y=-60, label="Antartic",
  #          color="blue")+
  # 
  # #Palearctic
  # annotate(geom="text", x=24, y=94, label="Palearctic",
  #          color="blue")+
  # 
  # #Nearctic
  # annotate(geom="text", x=-118, y=94, label="Nearctic",
  #          color="blue")+
  # 
  # #Neotropical
  # annotate(geom="text", x=-105, y=-28, label="Neotropical",
  #          color="blue")+
  # 
  # #Afrotropical
  # annotate(geom="text", x=-13, y=-3, label="Afrotropical",
  #          color="blue")+
  # 
  # #Indomalaysian
  # annotate(geom="text", x=86, y=-12, label="Indomalaysian",
  #          color="blue")+
  # 
  # #Australasian
  # annotate(geom="text", x=189, y=-19, label="Oceanian",
  #          color="blue")+
  
  theme_bw() +
  theme(
    axis.text.x  = element_blank(), 
    axis.text.y  = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(), 
    axis.line.x = element_blank(), 
    axis.line.y = element_blank(),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),                                          
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(), 
    axis.ticks = element_blank(),
    plot.margin = unit(c(1,1,1,1), 'cm'),
    plot.title = element_text(size = 18, vjust = 1, hjust = 0),
    legend.text = element_text(size = 12),          
    legend.title = element_blank(),                              
    legend.position = c(0.1, 0.2), 
    legend.key = element_blank(),
    legend.background = element_rect(color = "black", 
                                     fill = "white", 
                                     size = 2, linetype = "blank")))

map1


# Wordcloud ---------------------------------------------------------------

db_full <- read.csv(file = "Data/Biodiversity_WOS_V1.csv",sep='\t', dec='.',header=TRUE,as.is=FALSE)

db_full$title <- as.character(db_full$title)

# A list of boring and non-useful words, bundled with `tidytext`
data(stop_words)

# remove all numbers from titles
db_full$title <- gsub('[0-9]+', '', db_full$title)

# Count words
Title <- db_full  %>%
  unnest_tokens(output = title_word,
                input = title) %>%
  anti_join(stop_words, by = c("title_word" = "word")) %>%
  count(title_word, sort = TRUE) 

# Plot
dev.off()
Title[3:300,] %>% with(wordcloud(words = title_word, 
                                 freq = n, 
                                 max.words = 200,
                                 scale=c(4,.2),
                                 random.color=TRUE, color = c("orange","aquamarine3","aquamarine4","darkblue","black")))


#### APPUNTI:::

# Trends over time in methods and geography -------------------------------

t1 <- semi_colon_splitter(input1 = db$Method_data_collection,
                          input2 = db$Publication_year, 
                          names = c("Method","year"))

t1_prop <- t1  %>% group_by(year) %>% count(Method)
t1_tot <- data.frame(table(t1$year)) ; colnames(t1_tot) <- c("year", "Tot")

t1 <- data.frame(dplyr::left_join(t1_prop,t1_tot, by = "year"))

rm(t1_prop,t1_tot) #cleanù

t1$year <- as.numeric(as.character(t1$year))

# Modelling the temporal trends
model   <- list()
par     <- list()

for (i in levels(factor(t1$Method))) {
  
  db_i <- t1[t1$Method==i, ]
  model[[i]]   <- glm(cbind(n,Tot) ~ year, data = db_i, family = "binomial")
  par[[i]]     <- parameters::model_parameters(model[[i]])
  
}  

# m_t1 <- glm(cbind(n,Tot) ~ year * Geography, data = t1, 
#             family = "binomial")
# summary(m_t1)

# Summary
# for (i in 1:nlevels(factor(t1$Method))) {
#   
#   message(paste("::::::  ",levels(factor(t1$Method))[i],"  :::::"))
#   print(parameters::model_parameters(model[[i]]))   
#   message(paste(":::::::::::"))
#   
# }  

y2 <- seq(from = min(t1$year), to = max(t1$year), 1) #temporal series of interest

(Plot_trend1 <- ggplot() +
    ylab("Relative proportion of studies") + xlab(NULL) +
    #trend lines
    geom_line(aes(y = logisticline(y2,model[[1]]), x = y2), colour = COL[1],linetype="solid",size=1.1,alpha=1)+
    geom_line(aes(y = logisticline(y2,model[[2]]), x = y2), colour = COL[2],linetype="solid",size=1.1,alpha=1)+
    geom_line(aes(y = logisticline(y2,model[[3]]), x = y2), colour = COL[3],linetype="solid",size=1.1,alpha=1)+
    geom_line(aes(y = logisticline(y2,model[[4]]), x = y2), colour = COL[4],linetype="solid",size=1.1,alpha=1)+
    geom_line(aes(y = logisticline(y2,model[[5]]), x = y2), colour = COL[5],linetype="solid",size=1.1,alpha=1)+
    #confidence intervals
    geom_ribbon(aes(ymax = logisticline_max(y2, model[[01]]),
                    ymin = logisticline_min(y2, model[[01]]),x = y2),alpha = 0.5,fill=COL[1])+
    geom_ribbon(aes(ymax = logisticline_max(y2, model[[02]]),
                    ymin = logisticline_min(y2, model[[02]]),x = y2),alpha = 0.5,fill=COL[2])+
    geom_ribbon(aes(ymax = logisticline_max(y2, model[[03]]),
                    ymin = logisticline_min(y2, model[[03]]),x = y2),alpha = 0.5,fill=COL[3])+
    geom_ribbon(aes(ymax = logisticline_max(y2, model[[04]]),
                    ymin = logisticline_min(y2, model[[04]]),x = y2),alpha = 0.5,fill=COL[4])+
    geom_ribbon(aes(ymax = logisticline_max(y2, model[[05]]),
                    ymin = logisticline_min(y2, model[[05]]),x = y2),alpha = 0.5,fill=COL[5])+
   
    #Text
    annotate(geom="text", hjust = 0,vjust = 0.3,
             x= 2020.5, y= logisticline_max(y2, model[[01]])[21], 
             label = levels(factor(t1$Method))[1],
             color=COL[1],alpha=1)+
    
    annotate(geom="text", hjust = 0,vjust = 0,
             x= 2020.5, y= logisticline_max(y2, model[[02]])[21], 
             label = levels(factor(t1$Method))[2],
             color=COL[2],alpha=1)+
    
    annotate(geom="text", hjust = 0,vjust = 0,
             x= 2020.5, y= logisticline_max(y2, model[[03]])[21], 
             label = levels(factor(t1$Method))[3],
             color=COL[3],alpha=1)+
    
    annotate(geom="text", hjust = 0,vjust = 0,
             x= 2020.5, y= logisticline_max(y2, model[[04]])[21], 
             label = levels(factor(t1$Method))[4],
             color=COL[4],alpha=1)+
    
    annotate(geom="text", hjust = 0,vjust = 0,
             x= 2020.5, y= logisticline_max(y2, model[[05]])[21], 
             label = levels(factor(t1$Method))[5],
             color=COL[5])+
     
    coord_cartesian(xlim = c(1992, 2020), # This focuses the x-axis on the range of interest
                    clip = 'off') +   # This keeps the labels from disappearing
    
    theme_classic() + theme(plot.margin = unit(c(0.5,4,0.5,0.5), 'cm'))
)

################

t2 <- semi_colon_splitter(input1 = db$Geography,
                          input2 = db$Publication_year, 
                          names = c("Geography","year"))

t2_prop <- t2  %>% group_by(year) %>% count(Geography)
t2_tot <- data.frame(table(t2$year)) ; colnames(t2_tot) <- c("year", "Tot")

t2 <- data.frame(dplyr::left_join(t2_prop,t2_tot, by = "year"))

rm(t2_prop,t2_tot) #cleanù

t2$year <- as.numeric(as.character(t2$year))

# Modelling the temporal trends
model   <- list()
par     <- list()

for (i in levels(factor(t2$Geography))) {
  
  db_i <- t2[t2$Geography==i, ]
  model[[i]]   <- glm(cbind(n,Tot) ~ year, data = db_i, family = "binomial")
  par[[i]]     <- parameters::model_parameters(model[[i]])
  
}  

# m_t2 <- glm(cbind(n,Tot) ~ year * Geography, data = t2, 
#             family = "binomial")
# summary(m_t2)

# Summary
# for (i in 1:nlevels(factor(t2$Geography))) {
# 
#   message(paste("::::::  ",levels(factor(t2$Geography))[i],"  :::::"))
#   print(parameters::model_parameters(model[[i]]))
#   message(paste(":::::::::::"))
# 
# }

y2 <- seq(from = min(t2$year), to = max(t2$year), 1) #temporal series of interest

(Plot_trend2 <- ggplot() +
    ylab("Relative proportion of studies") + xlab(NULL)+ #ylim(0,0.3)+
    #trend lines
    geom_line(aes(y = logisticline(y2,model[[1]]), x = y2), colour = COL[1],linetype="solid",size=1.1,alpha=1)+
    geom_line(aes(y = logisticline(y2,model[[2]]), x = y2), colour = COL[2],linetype="solid",size=1.1,alpha=1)+
    geom_line(aes(y = logisticline(y2,model[[3]]), x = y2), colour = COL[3],linetype="solid",size=1.1,alpha=1)+
    geom_line(aes(y = logisticline(y2,model[[4]]), x = y2), colour = COL[4],linetype="solid",size=1.1,alpha=1)+
    geom_line(aes(y = logisticline(y2,model[[5]]), x = y2), colour = COL[5],linetype="solid",size=1.1,alpha=1)+
    geom_line(aes(y = logisticline(y2,model[[6]]), x = y2), colour = COL[6],linetype="solid",size=1.1,alpha=1)+
    geom_line(aes(y = logisticline(y2,model[[7]]), x = y2), colour = COL[7],linetype="solid",size=1.1,alpha=1)+
    geom_line(aes(y = logisticline(y2,model[[8]]), x = y2), colour = COL[8],linetype="solid",size=1.1,alpha=1)+
    #confidence intervals
    geom_ribbon(aes(ymax = logisticline_max(y2, model[[01]]),
                    ymin = logisticline_min(y2, model[[01]]),x = y2),alpha = 0.5,fill=COL[1])+
    geom_ribbon(aes(ymax = logisticline_max(y2, model[[02]]),
                    ymin = logisticline_min(y2, model[[02]]),x = y2),alpha = 0.5,fill=COL[2])+
    geom_ribbon(aes(ymax = logisticline_max(y2, model[[03]]),
                    ymin = logisticline_min(y2, model[[03]]),x = y2),alpha = 0.5,fill=COL[3])+
    geom_ribbon(aes(ymax = logisticline_max(y2, model[[04]]),
                    ymin = logisticline_min(y2, model[[04]]),x = y2),alpha = 0.5,fill=COL[4])+
    geom_ribbon(aes(ymax = logisticline_max(y2, model[[05]]),
                    ymin = logisticline_min(y2, model[[05]]),x = y2),alpha = 0.5, fill=COL[5])+
    geom_ribbon(aes(ymax = logisticline_max(y2, model[[06]]),
                    ymin = logisticline_min(y2, model[[06]]),x = y2),alpha = 0.5, fill=COL[6])+
    geom_ribbon(aes(ymax = logisticline_max(y2, model[[07]]),
                    ymin = logisticline_min(y2, model[[07]]),x = y2),alpha = 0.5, fill=COL[7])+
    geom_ribbon(aes(ymax = logisticline_max(y2, model[[08]]),
                    ymin = logisticline_min(y2, model[[08]]),x = y2),alpha = 0.5, fill=COL[8])+
    #Text
    annotate(geom="text", hjust = 0,vjust = 0.3,
             x= 2020.5, y= logisticline_max(y2, model[[01]])[21], 
             label = levels(factor(t2$Geography))[1],
             color=COL[1],alpha=1)+
    
    annotate(geom="text", hjust = 0,vjust = 0,
             x= 2020.5, y= logisticline_max(y2, model[[02]])[21], 
             label = levels(factor(t2$Geography))[2],
             color=COL[2],alpha=1)+
    
    annotate(geom="text", hjust = 0,vjust = 0,
             x= 2020.5, y= logisticline_max(y2, model[[03]])[21], 
             label = levels(factor(t2$Geography))[3],
             color=COL[3],alpha=1)+
    
    annotate(geom="text", hjust = 0,vjust = 0,
             x= 2020.5, y= logisticline_max(y2, model[[04]])[21], 
             label = levels(factor(t2$Geography))[4],
             color=COL[4],alpha=1)+
    
    annotate(geom="text", hjust = 0,vjust = 0,
             x= 2020.5, y= logisticline_max(y2, model[[05]])[21], 
             label = levels(factor(t2$Geography))[5],
             color=COL[5])+
    
    annotate(geom="text", hjust = 0,vjust = 0,
             x= 2020.5, y= logisticline_max(y2, model[[06]])[21], 
             label = levels(factor(t2$Geography))[6],
             color=COL[6])+
    
    annotate(geom="text", hjust = 0,vjust = 0,
             x= 2020.5, y= logisticline_max(y2, model[[07]])[21], 
             label = levels(factor(t2$Geography))[7],
             color=COL[7])+
    
    annotate(geom="text", hjust = 0,vjust = 0,
             x= 2020.5, y= logisticline_max(y2, model[[08]])[21], 
             label = levels(factor(t2$Geography))[8],
             color=COL[8])+
    
    coord_cartesian(xlim = c(1992, 2020), # This focuses the x-axis on the range of interest
                    clip = 'off') +   # This keeps the labels from disappearing
    
    theme_classic() + theme(plot.margin = unit(c(0.5,4,0.5,0.5), 'cm'))
)

################

t3_prop <- db  %>% group_by(Publication_year) %>% 
           summarize(Taxonomic_div = sum(Taxonomic_div,na.rm = T ),
                     Phylogenetic_div = sum(Phylogenetic_div,na.rm = T),
                     Functional_div = sum(Functional_div,na.rm = T),
                     Other_div = sum(Other_div,na.rm = T))
                     
t3 <- data.frame(t3_prop, Tot = rowSums(t3_prop[,2:5])) 
colnames(t3)[1] <- "year"
  
# Modelling the temporal trends
model   <- list()

model[[1]]   <- glm(cbind(Taxonomic_div,Tot) ~ year, data = t3, family = "binomial")
model[[2]]   <- glm(cbind(Phylogenetic_div,Tot) ~ year, data = t3, family = "binomial")
model[[3]]   <- glm(cbind(Functional_div,Tot) ~ year, data = t3, family = "binomial")
model[[4]]   <- glm(cbind(Other_div,Tot) ~ year, data = t3, family = "binomial")

y2 <- seq(from = min(t3$year), to = max(t3$year), 1) #temporal series of interest

(Plot_trend3 <- ggplot() +
    ylab("Relative proportion of studies") + xlab(NULL)+ #ylim(0,0.3)+
    #trend lines
    geom_line(aes(y = logisticline(y2,model[[1]]), x = y2), colour = COL[1],linetype="solid",size=1.1,alpha=1)+
    geom_line(aes(y = logisticline(y2,model[[2]]), x = y2), colour = COL[2],linetype="solid",size=1.1,alpha=1)+
    geom_line(aes(y = logisticline(y2,model[[3]]), x = y2), colour = COL[3],linetype="solid",size=1.1,alpha=1)+
    geom_line(aes(y = logisticline(y2,model[[4]]), x = y2), colour = COL[4],linetype="solid",size=1.1,alpha=1)+
    #confidence intervals
    geom_ribbon(aes(ymax = logisticline_max(y2, model[[01]]),
                    ymin = logisticline_min(y2, model[[01]]),x = y2),alpha = 0.5,fill=COL[1])+
    geom_ribbon(aes(ymax = logisticline_max(y2, model[[02]]),
                    ymin = logisticline_min(y2, model[[02]]),x = y2),alpha = 0.5,fill=COL[2])+
    geom_ribbon(aes(ymax = logisticline_max(y2, model[[03]]),
                    ymin = logisticline_min(y2, model[[03]]),x = y2),alpha = 0.5,fill=COL[3])+
    geom_ribbon(aes(ymax = logisticline_max(y2, model[[04]]),
                    ymin = logisticline_min(y2, model[[04]]),x = y2),alpha = 0.5,fill=COL[4])+
    #Text
    annotate(geom="text", hjust = 0,vjust = 0.3,
             x= 2020.5, y= logisticline_max(y2, model[[01]])[21], 
             label = "Taxonomic diversity",
             color=COL[1],alpha=1)+
    
    annotate(geom="text", hjust = 0,vjust = 0,
             x= 2020.5, y= logisticline_max(y2, model[[02]])[21], 
             label = "Phylogenetic diversity",
             color=COL[2],alpha=1)+
    
    annotate(geom="text", hjust = 0,vjust = 0,
             x= 2020.5, y= logisticline_max(y2, model[[03]])[21], 
             label = "Functional diversity",
             color=COL[3],alpha=1)+
    
    annotate(geom="text", hjust = 0,vjust = 0,
             x= 2020.5, y= logisticline_max(y2, model[[04]])[21], 
             label = "Other diversity",
             color=COL[4],alpha=1)+
    
    coord_cartesian(xlim = c(1992, 2020), # This focuses the x-axis on the range of interest
                    clip = 'off') +   # This keeps the labels from disappearing
    
    theme_classic() + theme(plot.margin = unit(c(0.5,4,0.5,0.5), 'cm'))
)
