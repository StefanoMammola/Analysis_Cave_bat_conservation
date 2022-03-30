## ------------------------------------------------------------------------
## 'Biodiversity misuse'
## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
# 'R script to reproduce the full analysis'
## ------------------------------------------------------------------------

# Analysis performed with R (v. R 4.1.0) and R studio (v. 1.4.1103)
# Authors: Stefano Mammola

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

# Data cleaning -----------------------------------------------------------

# Selecting paper to analyse
db <- db %>% filter(Analysis == "yes") ; db <- droplevels(db)
dim(db)

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

#biodiversity facets (%)
round((table(db$Facets_biodiversity)/sum(table(db$Facets_biodiversity)))*100,2)

table(db$Phylogenetic_div)
table(db$Functional_div)
table(db$Taxonomic_div)


# Sorting factors ---------------------------------------------------------

# Data exploration -------------------------------------------------------

###### temporal trends ######

# all
range(db$Biodiversity_prop, na.rm = TRUE) ; mean(db$Biodiversity_prop, na.rm = TRUE) ; std(db$Biodiversity_prop)
nrow(db)

(plot1a <- ggplot(data = db, aes(x = Publication_year, y = Biodiversity_prop)) + 
  geom_point(col = "grey10", fill = "grey30", size = 5, shape = 21, alpha = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, col="blue", fill = "blue") +
  labs(x = NULL, y = Y.label)+
  annotate(geom="text", x=2000, y=0.35, 
           label=expression(paste("Average biodiversity proportion" %+-% "SE: 0.03" %+-% "0.001")),
           color="grey10", size = 5)+
  theme_custom()) #warnings() are due to NA removal

#Split multiple regions separated by ";"
box1 <- semi_colon_splitter(input1 = db$Geography,
                            input2 = db$Biodiversity_prop, 
                            names  = c("Geography","Biodiversity_prop"))

box1$Biodiversity_prop <- as.numeric(as.character(box1$Biodiversity_prop))

#Sort levels
box1$Geography <- factor(box1$Geography,
                         c(levels(box1$Geography)[4],levels(box1$Geography)[-4]))

(plot1b <- ggplot(data = box1 %>% drop_na(Geography,Biodiversity_prop), 
                  aes(x = Geography, y = Biodiversity_prop)) +
    geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.4, fill= "blue", col = "white") +
    geom_point(position = position_jitter(width = 0.15), size = 1, alpha = 0.7, color = "grey40") +
    geom_boxplot(width = 0.2,  col = "blue", outlier.shape = NA, alpha = 0) +
    labs(y = Y.label, x = NULL) +
    theme_custom())

#Split multiple domains separated by ";"
box2 <- semi_colon_splitter(input1 = db$Domain,
                            input2 = db$Biodiversity_prop, 
                            names = c("Domain","Biodiversity_prop"))

box2$Biodiversity_prop <- as.numeric(as.character(box2$Biodiversity_prop))

(plot1c <- ggplot(data = box2 %>% drop_na(Domain,Biodiversity_prop), 
                  aes(x = Domain, y = Biodiversity_prop)) +
    geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.4, fill= "blue", col = "white") +
    geom_point(position = position_jitter(width = 0.15), size = 1, alpha = 0.7, color = "grey40") +
    geom_boxplot(width = 0.2, col = "blue", outlier.shape = NA, alpha = 0) +
    labs(y = Y.label, x = NULL) +
    theme_custom())

#Split multiple methds separated by ";"
box3 <- semi_colon_splitter(input1 = db$Method_data_collection,
                            input2 = db$Biodiversity_prop, 
                            names = c("Method","Biodiversity_prop"))

box3$Biodiversity_prop <- as.numeric(as.character(box3$Biodiversity_prop))

#Sort levels
box3$Method <- factor(box3$Method,
                         c(levels(box3$Method)[-4],levels(box3$Method)[4]))


levels(box3$Method)[c(2,5)] <- "Other"

(plot1d <- ggplot(data = box3 %>% drop_na(Method,Biodiversity_prop), 
                  aes(x = Method, y = Biodiversity_prop)) +
    geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.4, fill= "blue", col = "white") +
    geom_point(position = position_jitter(width = 0.15), size = 1, alpha = 0.7, color = "grey40") +
    geom_boxplot(width = 0.2, col = "blue", outlier.shape = NA, alpha = 0) +
    labs(y = Y.label, x = NULL) +
    theme_custom())

pdf(file = "Figure/Figure_1.pdf", width = 19, height = 14)

ggpubr::ggarrange(plot1a, plot1b, plot1c, plot1d,
                  common.legend = FALSE,
                  hjust = -5,
                  align = "hv",
                  labels = c("A", "B", "C", "D"),
                  ncol = 2, nrow = 2) #warnings() are due to NA removal

dev.off()

# Figure 2 ----------------------------------------------------------------

#Loading silouhettes
animal_png <- png::readPNG("Phylopics/Animal.png")
fungi_png  <- png::readPNG("Phylopics/Fungi.png")
micro_png  <- png::readPNG("Phylopics/Micro.png")
plant_png  <- png::readPNG("Phylopics/Plant.png")

# Plotting
(plot2a <- ggplot(data = db[db$Animals_prop>0,], aes(x = Publication_year, y = Animals_prop)) + 
    geom_point(col = "grey10", fill = "grey30", size = 5, shape = 21, alpha = 0.3) +
    geom_smooth(method = "glm", formula = y ~ x, 
                method.args = list(family = quasibinomial(link = "logit")), 
                col = "blue", fill = "blue") +
    labs(title = "Only animals", x = NULL , y = Y.label)+
    xlim(1992,2020)+
    ylim(0,1)+
    annotation_custom(grid::rasterGrob(animal_png), xmin = 1992, xmax = 2000, ymin = 0.75, ymax = 1)+ 
    theme_custom()) #warnings() are due to NA removal

(plot2b <- ggplot(db[db$Plants_prop>0,], aes(x = Publication_year, y = Plants_prop)) + 
    geom_point(col = "grey10", fill = "grey30", size = 5, shape = 21, alpha = 0.3) +
    geom_smooth(method = "glm", formula = y ~ x, 
                method.args = list(family = quasibinomial(link = "logit")), 
                col = "blue", fill = "blue") +
    labs(title = "Only plants", x = NULL , y = NULL)+
    xlim(1992,2020)+
    ylim(0,1)+
    annotation_custom(grid::rasterGrob(plant_png), xmin = 1992, xmax = 1995, ymin = 0.75, ymax = 1)+ 
    theme_custom()) #warnings() are due to NA removal

(plot2c <- ggplot(data = db[db$Fungi_prop>0,], aes(x = Publication_year, y = Fungi_prop)) + 
    geom_point(col = "grey10", fill = "grey30", size = 5, shape = 21, alpha = 0.3) +
    geom_smooth(method = "glm", formula = y ~ x, 
                method.args = list(family = quasibinomial(link = "logit")), 
                col = "blue", fill = "blue") +
    labs(title = "Only fungi", x = "Publication year" , y = Y.label)+
    xlim(1992,2020)+
    ylim(0,1)+
    annotation_custom(grid::rasterGrob(fungi_png), xmin = 1992, xmax = 1997, ymin = 0.75, ymax =1)+ 
    theme_custom()) #warnings() are due to NA removal

(plot2d <- ggplot(data = db[db$Micro_prop>0,], aes(x = Publication_year, y = Micro_prop)) + 
    geom_point(col = "grey10", fill = "grey30", size = 5, shape = 21, alpha = 0.3) +
    geom_smooth(method = "glm", formula = y ~ x, 
                method.args = list(family = quasibinomial(link = "logit")), 
                col = "blue", fill = "blue") +
    labs(title = "Only microorganisms", x = "Publication year" , y = NULL)+
    xlim(1992,2020)+
    ylim(0,1)+
    annotation_custom(grid::rasterGrob(micro_png), xmin = 1992, xmax = 1997, ymin = 0.75, ymax = 1)+
    theme_custom()) #warnings() are due to NA removal

pdf(file = "Figure/Figure_2.pdf", width = 16, height =12)
ggpubr::ggarrange(plot2a,plot2b,plot2c,plot2d,
                  common.legend = FALSE,
                  hjust = -5,
                  align = "hv",
                  labels = c("A", "B", "C", "D"),
                  ncol=2, nrow=2) #warnings() are due to NA removal
dev.off()

# Figure 3 ----------------------------------------------------------------

# Regression model -------------------------------------------------------

colnames(db)

# Subset
db_glm <- db %>% select(year = Publication_year,
                        n_aut,
                        cit = tot_cites,
                        Method = Method_data_collection,
                        Phylogenetic_div,      
                        Functional_div,
                        Title_geo,
                        Title_hab,
                        Title_taxon,
                        Title_adjecties,
                        Biodiversity_prop,
                        Geography,
                        Domain)

# Converting multiples
method_split <- strsplit(as.character(db_glm$Method), ";")

method <- c()
for(i in 1:length(method_split))
  method <- c(method, ifelse(length(method_split[[i]]) > 1, "Multiple", method_split[[i]]) )

geography_split <- strsplit(as.character(db_glm$Geography), ";")

geography <- c()
for(i in 1:length(geography_split))
  geography <- c(geography, ifelse(length(geography_split[[i]]) > 1, "Global", geography_split[[i]]) )

domain_split <- strsplit(as.character(db_glm$Domain), ";")

domain <- c()
for(i in 1:length(domain_split))
  domain <- c(domain, ifelse(length(domain_split[[i]]) > 1, "Multiple", domain_split[[i]]) )

db_glm$Method    <- method
db_glm$Geography <- geography
db_glm$Domain    <- domain

db_glm <- db_glm %>% mutate_at(vars("Method","Geography","Domain",
                                    "Title_geo","Title_hab","Title_taxon"), as_factor)

rm(method_split,method,geography,geography_split,domain,domain_split) #clean

# Data exploration --------------------------------------------------------

# Checking balancing of factors
table(db_glm$Method) #Citizen science/simulation too few records
levels(db_glm$Method)[c(6,7)] <- "Other"

table(db_glm$Geography) #Antartica a bit unbalanced

table(db_glm$Domain) #OK
db_glm <- droplevels(db_glm)

# Correcting citation by year

## What is the trend of citation over time?
ggplot(data = db_glm, aes(x = year, y = cit)) + 
  geom_point(size = 1, alpha = 0.7, color = "grey40")+
  labs(x = "year of publication", y = "citation")+
  theme_custom()

## Correcting citation by article age
M0 <- gam::gam(cit ~ s(year), family = poisson, data = db_glm) #constraining dof to 2 so that you don't get negative predicted citations
summary(M0)

par(mar= c(rep(2,4)))
plot(M0, se = TRUE)

## Taking the residual of citations
db_glm <- data.frame(db_glm, citation_residuals = resid(M0,type="pearson"))

#Checking how good it predict by comparing with a simple normalization by year of publication

ggplot(data = db_glm, aes(x = citation_residuals, 
                          y = cit/year)) + #normalize by year
         geom_point(size = 1, alpha = 0.7, color = "grey40")+
         labs(x = "citation residuals", y = "citation normalized")+
         theme_custom()

# Dependent var
db_glm$prop  <- rowSums(db[,36:91])
db_glm$total <- length(36:91)

# Checking outliers
par(mar= c(rep(2,4)))
dotchart(db_glm$citation_residuals) # 1 outlier
dotchart(db_glm$prop) #2 outliers

# Removing outliers
db_glm <- db_glm[db_glm$citation_residuals < 200,]
db_glm <- db_glm[db_glm$prop < 20,]

#Set baseline
db_glm <- within(db_glm, Geography <- relevel(Geography, ref = "Global"))
db_glm <- within(db_glm, Domain    <- relevel(Domain,    ref = "Multiple"))
db_glm <- within(db_glm, Method    <- relevel(Method,    ref = "Multiple"))

# Modelling  --------------------------------------------------------------

#scale contionuos variables
db_glm$year <- scale(db_glm$year)
db_glm$citation_residuals <- scale(db_glm$citation_residuals)

#Initial model
m1  <- glm(cbind(prop,total) ~ 
             year + 
             citation_residuals + 
             Domain + 
             Geography + 
             Method +
             Title_geo + 
             Title_hab + 
             Title_taxon, 
           data = db_glm, family = "binomial")

parameters(m1)
performance::check_overdispersion(m1) #overdispersed
performance::check_collinearity(m1)

#Refit with quasibinomial due to overdispersion 
m1b  <- glm(cbind(prop,total) ~ 
              year + 
              citation_residuals + 
              Domain + 
              Geography + 
              Method +
              Title_geo + 
              Title_hab + 
              Title_taxon, 
            data = db_glm, family = "quasibinomial")

parameters(m1b)
performance::check_collinearity(m1b)

# posthoc
pairs(emmeans::emmeans(m1b, ~ Domain), simple=c("Domain"))
pairs(emmeans::emmeans(m1b, ~ Geography), simple="Geography")
pairs(emmeans::emmeans(m1b, ~ Method), simple="Method")

# sjPlot::plot_model(m1, title ="Factors correlating with biodiversity proportion",
#                    sort.est = FALSE,  vline.color = "grey80",
#                    show.values = TRUE, value.offset = .3, se = TRUE, show.p = TRUE) + theme_classic()
# 

# Plot...

# Estract estimates
Estimates_m1 <- 
  m1b %>% 
  summary %>% 
  magrittr::extract2("coefficients") %>% # extract estimates
  as.data.frame %>% rownames_to_column("Variable") %>% 
  dplyr::filter(!row_number() %in% 1) %>%  #remove intercept
  dplyr::rename(SE = 3, z = 4, p = 5) #rename

#Set variable order and rename
order_var <- c("Year of publication",
               "Citations (corrected by year)",
               "Domain [Terrestrial]",
               "Domain [Saltwater]",
               "Domain [Freshwater]",
               "Geographic [Palearctic]",
               "Geographic [Afrotropical]",
               "Geographic [Indomalayan]",
               "Geographic [Neartic]",
               "Geographic [Australasian]",
               "Geographic [Antartic]",
               "Geographic [Neotropical]",
               "Method [Review/Opinion]",
               "Method [Field sampling]",
               "Method [Big data]",
               "Method [Other]",
               "Mention of location in title",
               "Mention of habitat in title",
               "Mention of taxon/a in title")

Estimates_m1$Variable <- order_var #Rename
Estimates_m1$Variable <- factor(Estimates_m1$Variable, rev(order_var)) #sort

sign <- ifelse(Estimates_m1$p > 0.05, "", " *")
col_p <- ifelse(Estimates_m1$p > 0.05, "grey5", ifelse(Estimates_m1$Estimate>0,"orange","blue") )

(plot_model3 <- ggplot2::ggplot(data = Estimates_m1) +

    geom_pointrange(aes(x = Variable, 
                        y = Estimate,
                        ymin = Estimate-SE, 
                        ymax = Estimate+SE), col = col_p, size = 0.5) + 
    
    geom_hline(lty = 3, size = 0.7, col = "grey50", yintercept = 0) +
    
    geom_text(aes(Variable, Estimate),
              label = paste0(round(Estimates_m1$Estimate,2),sign), 
              vjust = -1, size = 3, col = col_p) +
   
    labs(y = expression(paste("Odds ratio" %+-% "Standard Error")),
         x = NULL)+
    theme_custom() + theme(axis.text.y  = element_text(colour = rev(col_p))) + coord_flip()
    # geom_vline(lty = 1, size = 0.2, col = "blue", xintercept = 3.5)+
    # geom_vline(lty = 1, size = 0.2, col = "blue", xintercept = 7.5)+
    # geom_vline(lty = 1, size = 0.2, col = "blue", xintercept = 14.5)+
    # geom_vline(lty = 1, size = 0.2, col = "blue", xintercept = 17.5)+
    # geom_vline(lty = 1, size = 0.2, col = "blue", xintercept = 18.5)
)

pdf(file = "Figure/Figure_3.pdf", width = 12, height =8)
plot_model3
dev.off()


# Modelling #2  ------------------------------------------------------------

# Repeating the analysis only with generic titles

db_glm2 <- db_glm[db_glm$Title_taxon == 0,]
db_glm2 <- db_glm2[db_glm2$Title_hab == 0,]
db_glm2 <- db_glm2[db_glm2$Title_geo == 0,]
db_glm2 <- db_glm2[db_glm2$Geography != "Antarctic",] ; db_glm2$Geography <- droplevels(db_glm2$Geography)

#Balance of factors
table(db_glm2$Method)
table(db_glm2$Geography)
table(db_glm2$Domain) 

#Initial model
m2  <- glm(cbind(prop,total) ~ 
             year + 
             citation_residuals + 
             Domain + 
             Geography + 
             Method, 
           data = db_glm2, family = "binomial")

performance::check_overdispersion(m2) #overdispersed

#Refit with quasibinomial due to overdispersion 
m2b  <- glm(cbind(prop,total) ~ 
              year + 
              citation_residuals + 
              Domain + 
              Geography + 
              Method, 
           data = db_glm2, family = "quasibinomial")

parameters(m2b)
performance::check_collinearity(m2b)

# posthoc
pairs(emmeans::emmeans(m2b, ~ Domain), simple=c("Domain"))
pairs(emmeans::emmeans(m2b, ~ Geography), simple="Geography")
pairs(emmeans::emmeans(m2b, ~ Method), simple="Method")

# sjPlot::plot_model(m1, title ="Factors correlating with biodiversity proportion",
#                    sort.est = FALSE,  vline.color = "grey80",
#                    show.values = TRUE, value.offset = .3, se = TRUE, show.p = TRUE) + theme_classic()
# 

# Plot...

# Estract estimates
Estimates_m2 <- 
  m2b %>% 
  summary %>% 
  magrittr::extract2("coefficients") %>% # extract estimates
  as.data.frame %>% rownames_to_column("Variable") %>% 
  dplyr::filter(!row_number() %in% 1) %>%  #remove intercept
  dplyr::rename(SE = 3, z = 4, p = 5) #rename

#Set variable order and rename
order_var <- c("Year of publication",
               "Citations (corrected by year)",
               "Domain [Terrestrial]",
               "Domain [Saltwater]",
               "Domain [Freshwater]",
               "Geographic [Palearctic]",
               "Geographic [Afrotropical]",
               "Geographic [Indomalayan]",
               "Geographic [Neartic]",
               "Geographic [Australasian]",
               "Geographic [Neotropical]",
               "Method [Review/Opinion]",
               "Method [Field sampling]",
               "Method [Big data]",
               "Method [Other]")

Estimates_m2$Variable <- order_var #Rename
Estimates_m2$Variable <- factor(Estimates_m2$Variable, rev(order_var)) #sort

sign <- ifelse(Estimates_m2$p > 0.05, "", " *")
col_p <- ifelse(Estimates_m2$p > 0.05, "grey5", ifelse(Estimates_m2$Estimate>0,"orange","blue") )

(plot_model3bis <- ggplot2::ggplot(data = Estimates_m2) +
    
    geom_pointrange(aes(x = Variable, 
                        y = Estimate,
                        ymin = Estimate-SE, 
                        ymax = Estimate+SE), col = col_p, size = 0.5) + 
    
    geom_hline(lty = 3, size = 0.7, col = "grey50", yintercept = 0) +
    
    geom_text(aes(Variable, Estimate),
              label = paste0(round(Estimates_m2$Estimate,2),sign), 
              vjust = -1, size = 3, col = col_p) +
    
    labs(y = expression(paste("Odds ratio" %+-% "Standard Error")),
         x = NULL)+
    theme_custom() + theme(axis.text.y  = element_text(colour = rev(col_p))) + coord_flip()
  # geom_vline(lty = 1, size = 0.2, col = "blue", xintercept = 3.5)+
  # geom_vline(lty = 1, size = 0.2, col = "blue", xintercept = 7.5)+
  # geom_vline(lty = 1, size = 0.2, col = "blue", xintercept = 14.5)+
  # geom_vline(lty = 1, size = 0.2, col = "blue", xintercept = 17.5)+
  # geom_vline(lty = 1, size = 0.2, col = "blue", xintercept = 18.5)
)

pdf(file = "Figure/Figure_3bis.pdf", width = 12, height =8)
plot_model3bis
dev.off()

# Figure 4 ----------------------------------------------------------------

vector <- sort(apply(db[,36:86],2, sum, na.rm = TRUE), decreasing = TRUE)

db2 <- data.frame(Phyla = names(vector), N = vector)

Taxa_to_rename <- ifelse(db2$N<10,db2$N,NA)
N_to_rename    <- sum(Taxa_to_rename,na.rm=TRUE)
Taxa_to_rename <- length(na.omit(Taxa_to_rename))

db2 <- db2[ 1 : (nrow(db2)-Taxa_to_rename), ]

db2 <- rbind(db2, data.frame(Phyla = paste0("Others (n = ",Taxa_to_rename,")"), N = N_to_rename))

db2 <- cbind(db2, Type = c(rep("Animal",2),
                           "Plant",
                           rep("Animal",2),
                           "Plant",
                           "Microorganism",
                           "Plant",
                           rep("Animal",3),
                           "Plant",
                           "Animal",
                           "Plant",
                           "Animal",
                           "Fungi",
                           "Microorganism",
                           rep("Animal",2),
                           "Fungi",
                           "Multiple"
                           ))

db2$Phyla <- factor(db2$Phyla,levels = db2$Phyla)

(figure_4 <- ggplot(db2, aes(x= Phyla, y=N))+
    geom_bar(aes(fill= Type),stat="identity", alpha=1, colour = "black")  +
    
    labs(x = NULL, y = "Count")+
    
    scale_fill_manual(values = rev(RColorBrewer::brewer.pal(5, "Blues")))+
    theme_custom() + theme(axis.text.x = element_text(angle = 70, hjust=1 ))
)

pdf(file = "Figure/Figure_4.pdf", width = 12, height =8)
figure_4
dev.off()

## Only with titles containong no moderators

db_no_moderators <- db[c(db$Title_taxon | db$Title_hab | db$Title_geo) != 1,]

vector <- sort(apply(db_no_moderators[,36:86],2, sum, na.rm = TRUE), decreasing = TRUE)

db2 <- data.frame(Phyla = names(vector), N = vector)

Taxa_to_rename <- ifelse(db2$N<2,db2$N,NA)
N_to_rename    <- sum(Taxa_to_rename,na.rm=TRUE)
Taxa_to_rename <- length(na.omit(Taxa_to_rename))

db2 <- db2[ 1 : (nrow(db2)-Taxa_to_rename), ]

db2 <- rbind(db2, data.frame(Phyla = paste0("Others (n = ",Taxa_to_rename,")"), N = N_to_rename))

db2 <- cbind(db2, Type = c(rep("Animal",2),
                           "Plant",
                           rep("Animal",2),
                           "Plant",
                           "Microorganism",
                           "Plant",
                           rep("Animal",3),
                           "Plant",
                           "Animal",
                           "Plant",
                           "Animal",
                           "Fungi",
                           "Microorganism",
                           rep("Animal",2),
                           "Fungi",
                           "Multiple"
))

db2$Phyla <- factor(db2$Phyla,levels = db2$Phyla)

(figure_4bis <- ggplot(db2, aes(x= Phyla, y=N))+
    geom_bar(stat="identity", alpha=1, colour = "black")  +
    
    labs(x = NULL, y = "Count")+
    
    scale_fill_manual(values = rev(RColorBrewer::brewer.pal(5, "Blues")))+
    theme_custom() + theme(axis.text.x = element_text(angle = 70, hjust=1 ))
)

pdf(file = "Figure/Figure_4.pdf", width = 12, height =8)
figure_4bis
dev.off()

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

# Adjectives --------------------------------------------------------------

ggplot(data = db %>% drop_na(Title_adjecties,Biodiversity_prop), 
       aes(x = Title_adjecties, y = Biodiversity_prop)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.4, fill= "turquoise3") +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.4,fill= "orangered") +
  geom_point(position = position_jitter(width = 0.15), size = 1, alpha = 0.7, fill= "grey20") +
  labs(y = "Biodiversity (Proportion)", x = "N° of classifiers for biodiversity") +
  theme_classic()

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

# Regression model -------------------------------------------------------

colnames(db)

# Subset
db_glm <- db %>% select(year = Publication_year,
                        n_aut,
                        cit = tot_cites,
                        Method = Method_data_collection,
                        Facets_biodiversity,
                        Title_geo,
                        Title_hab,
                        Title_taxon,
                        Title_adjecties,
                        Biodiversity_prop,
                        Geography,
                        Domain)
str(db_glm)
db_glm$Facets_biodiversity <- as.numeric(db_glm$Facets_biodiversity)
# Converting multiples
method_split <- strsplit(as.character(db_glm$Method), ";")

method <- c()
for(i in 1:length(method_split))
  method <- c(method, ifelse(length(method_split[[i]]) > 1, "Multiple", method_split[[i]]) )

geography_split <- strsplit(as.character(db_glm$Geography), ";")

geography <- c()
for(i in 1:length(geography_split))
  geography <- c(geography, ifelse(length(geography_split[[i]]) > 1, "Global", geography_split[[i]]) )

domain_split <- strsplit(as.character(db_glm$Domain), ";")

domain <- c()
for(i in 1:length(domain_split))
  domain <- c(domain, ifelse(length(domain_split[[i]]) > 1, "Multiple", domain_split[[i]]) )

db_glm$Method    <- method
db_glm$Geography <- geography
db_glm$Domain    <- domain

db_glm <- db_glm %>% mutate_at(vars("Method","Geography","Domain",
                                    "Title_geo","Title_hab","Title_taxon"), as_factor)

rm(method_split,method,geography,geography_split,domain,domain_split) #clean

# Correcting citation by year (NOTE:: to do better with gam)
db_glm$cit <- db_glm$cit / db_glm$year

# Dependent var
db_glm$prop  <- rowSums(db[,36:91])
db_glm$total <- length(36:91)

#Set baseline
db_glm <- within(db_glm, Geography <- relevel(Geography, ref = "Global"))
db_glm <- within(db_glm, Domain    <- relevel(Domain, ref = "Multiple"))
db_glm <- within(db_glm, Method    <- relevel(Method, ref = "Multiple"))

# Fitting the model
m1  <- glm(cbind(prop,total) ~ year + n_aut + cit + Domain + Geography + 
             Title_geo + Title_hab  + Title_taxon,
           data = db_glm, family = "binomial")

performance::check_model(m1)
performance::check_overdispersion(m1)
performance::check_collinearity(m1)

cor(db_glm$n_aut,db_glm$cit)

str(db_glm)

db_glm %>% ggplot(aes(x=Geography , y = n_aut)) + geom_boxplot() + theme_classic()

# drop n_aut
m2  <- glm(cbind(prop,total) ~ year + 
                                cit + 
                              Domain + 
                           Geography + 
             Method +
                              Title_geo + Title_hab  + Title_taxon, 
           data = db_glm, family = "binomial")

performance::check_collinearity(m2)
performance::check_model(m2)
performance::check_overdispersion(m2)

sjPlot::plot_model(m2, title ="Factors correlating with biodiversity proportion",
                   sort.est = FALSE,  vline.color = "grey80",
                   show.values = TRUE, value.offset = .3, se = TRUE, show.p = TRUE) + theme_classic()


parameters::model_parameters(m2)

# quasibinomial?
table(db$Biodiversity_prop)

m2  <- glm(cbind(prop,total) ~ year + 
             cit + 
             Domain + 
             Geography + 
             Method +
             Title_geo + 
             Title_hab + 
             Title_taxon, data = db_glm, family = "binomial")

summary(m2)

sjPlot::plot_model(m2, title ="Factors correlating with biodiversity proportion",
                   sort.est = FALSE,  vline.color = "grey80",
                   show.values = TRUE, value.offset = .3, se = TRUE, show.p = TRUE) + theme_classic()

MuMIn::r.squaredLR(m2)[[1]]


m2_bis  <- glm(cbind(prop,total) ~ year + 
             cit + 
             Domain + 
             Geography + 
             Method +
             Title_geo + 
             Title_hab + 
             Title_taxon, data = db_glm[db_glm$prop>0,], family = "binomial")

sjPlot::plot_model(m2_bis, title ="Factors correlating with biodiversity proportion",
                   sort.est = FALSE,  vline.color = "grey80",
                   show.values = TRUE, value.offset = .3, se = TRUE, show.p = TRUE) + theme_classic()


# Wordcloud ---------------------------------------------------------------

db_full <- read.csv(file = "Data/Biodiversity_WOS_V2.csv",sep='\t', dec='.',header=TRUE,as.is=FALSE)

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

#end





