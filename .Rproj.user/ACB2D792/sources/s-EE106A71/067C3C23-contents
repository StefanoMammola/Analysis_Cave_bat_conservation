###############################################################

## Practical conservation of bats in subterranean habitats
## Meierhofer M., ... Mammola S.

## ------------------------------------------------------------------------
# 'R script to reproduce the analyses'
## ------------------------------------------------------------------------

# Analysis performed with R (v. R 4.1.0) and R studio (v. 1.4.1103)
# Authors (code): Stefano Mammola & Melissa B. Meierhofer

###############################################################

# clean the workspace -----------------------------------------------------

rm(list = ls())

# Loading R package -------------------------------------------------------

library("metafor")   # Meta-Analysis Package for R
library("ggplot2")       
library("tidyverse")


library("circlize")      
library("cowplot")
library("dplyr")         
library("ggplot2")       
library("ggpubr")
library("grid")
library("gridExtra")
library("maps")
library("rsq")
library("parameters")
library("performance")
library("scatterpie")
library("tidyr")
library("tidygraph")
library("ggraph")

# Sourcing useful functions ------------------------------------------------

source("Functions/Functions_bat.R")

###############################################################

## Data preparation:

###############################################################

# Loading the Database ----------------------------------------------------

db_full <-
  read.csv(
    file = "Data/Master_Database_Cave_Conservation.csv",
    sep = '\t',
    dec = ',',
    header = TRUE,
    as.is = FALSE
  )

# Removing duplicates

db_full <- db_full[db_full$Remove != "yes",] ; db_full <- droplevels(db_full)

#Subselecting BAT studies

db <- db_full[db_full$Bat_analysis == "yes",] ; db <- droplevels(db)

db$N <- as.numeric(db$N)

# How many study consider bats and other groups?
a <- ifelse(table(db$ID, db$Taxon_Group)>0,1,0)
table(rowSums(a))

sum(table(rowSums(a))[c(2,3)])/sum(table(rowSums(a)))*100 #3.87 % of studies consider multiple organisms

# Selecting only bats
db <- db[db$Taxon_Group == "Bats",] ; db <- droplevels(db)

#Database only with distinct paper
db_unique <- distinct(db, ID, .keep_all = TRUE) 

#Checking levels of factors
levels(db$Taxon_Group)
levels(db$Tested_statistically)
levels(db$Higher_Geography)
levels(db$System)
levels(db$Domain)
levels(db$Impact)
levels(db$Conservation_Group)
levels(db$Conservation_Action)

#Type of actions
table(db$Publication_type)

#Summary statistics (Literature)
table(db_unique$Source) ; sum(table(db_unique$Source)) # N° of unique sources

mean(table(db$ID)) ; SE(table(db$ID)) # mean number of actions/paper

#Summary statistics (Testing)
table(db$Tested_statistically)[2] / sum(table(db$Tested_statistically)) #N° and % testing

#How many estimates would be usable for meta analysis?
n_studies    <- c() 
n_estimates  <- c()
perc_testing <- c()
usable       <- c()
unusable     <- c()
perc_usable  <- c()

for(i in 1:nlevels(db$Conservation_Action)){

  db_i_tot <- db[db$Conservation_Action == levels(db$Conservation_Action)[i],]
  db_i     <- db_i_tot[db_i_tot$Tested_statistically == "yes",]
  
  table_i        <- table(db_i$Pearson_r_conversion) #% of usable statistics
  n_studies      <- c(n_studies, nrow(distinct(db_i, ID, .keep_all = TRUE)) ) #unique studies
  n_estimates    <- c(n_estimates, nrow(db_i) ) #unique estimates
  perc_testing   <- c(perc_testing, round(nrow(db_i)/nrow(db_i_tot),2)*100 )
  usable         <- c(usable, sum(table_i[1],table_i[3]))
  unusable       <- c(unusable, sum(table_i[2]))
  perc_usable    <- c(perc_usable, round((usable[i]/sum(table_i)),2)*100)
  
}

Table_1 <- data.frame(Intervention = levels(db$Conservation_Action), n_studies, n_estimates, perc_testing, usable, unusable, perc_usable)
Table_1[is.na(Table_1)] <- 0
colnames(Table_1) <- c("Intervention", "N° studies", "N° interventions", "% testing", "N° usable", "N° unusable", "% usable")

write.csv(Table_1,"Tables/Table_1.csv")

#Redefining levels
levels(db$Impact)[2] <- "Multiple"

levels(db$Conservation_Action)[5] <- "Gating"

###############################################################

## Meta-Analysis

###############################################################

db_metafor <- db[db$Tested_statistically == "yes",]
db_metafor <- db_metafor[db_metafor$Pearson_r_conversion == "converted",]
db_metafor <- droplevels(db_metafor)

dim(db_metafor)
nlevels(db_metafor$ID) #250 references 

db_metafor <- db_metafor %>% select(ID, 
                             N,
                             Domain,
                             System,
                             Family,
                             Response_Group,
                             Predictor_Group,
                             r = Pearson.s_r)

# Derive Fischer's Z and its variance

db_metafor <- metafor::escalc(measure = "COR", ri = r, ni = N, data = db_metafor)

# Gate
db_metafor <- db_metafor[db_metafor$Predictor_Group == "Gate" | 
                         db_metafor$Predictor_Group == "Disturbance reduction" |
                         db_metafor$Predictor_Group == "Restoration", ] ; db_metafor <- droplevels(db_metafor)

table(db_metafor$Predictor_Group,db_metafor$Response_Group) # Disturbance reduction & Gate

db_metafor <- db_metafor[!c(db_metafor$Predictor_Group == "Disturbance reduction" & db_metafor$Response_Group == "Population"),]

#Check sample size for each predictors
table_n <- data.frame(predictor = NULL, n = NULL, n_papers = NULL)

for(i in 1:length(unique(levels(db_metafor$Response_Group))))
  table_n <- rbind(table_n, 
                          data.frame(predictor = levels(db_metafor$Response_Group)[i],
                                     n = nrow(db_metafor[db_metafor$Response_Group == levels(db_metafor$Response_Group)[i], ]),
                                     n_papers = length(unique(db_metafor[db_metafor$Response_Group == levels(db_metafor$Response_Group)[i], ]$ID)))  
                          
  )

db_metafor <- db_metafor[db_metafor$Response_Group != "Survival",]

actions_to_analyse    <- c("Disturbance reduction", "Gate", "Restoration")

SUBSET    <- list()
MODEL     <- list()

result_for_plot <- data.frame(label_action = NULL,
                              label_pred = NULL,
                              size = NULL,
                              b     = NULL,
                              ci.lb = NULL,
                              ci.ub = NULL,
                              ES    = NULL,
                              L     = NULL,
                              U     = NULL)

# Modelling
for (j in 1:length(actions_to_analyse)){ 
  
  data_j  <- db_metafor[db_metafor$Predictor_Group == actions_to_analyse[j], ] ; data_j <- droplevels(data_j)
  
  predictors_to_analyse <- levels(data_j$Response_Group)
  
  
  for (i in 1:length(predictors_to_analyse)){  
    
    #subset the predictor
    data_i  <- data_j[data_j$Response_Group == predictors_to_analyse[i], ]
    
    if(nrow(data_i)<2) { NULL } else {
    
    #fitting the model
    model_i <- rma.mv(yi, vi, random =  ~ 1 | ID, data = na.omit(data_i)) 
    
    #extracting coefficients
    result_for_plot_i <- data.frame(label_action =  actions_to_analyse[j],
                                  label_pred = paste(predictors_to_analyse[i]," (" ,
                                                  nrow(data_i),", ",
                                                  length(unique(data_i$ID)),")",sep=''),
                                    size = length(unique(data_i$ID)),
                                    b     = model_i$b,
                                    ci.lb = model_i$ci.lb,
                                    ci.ub = model_i$ci.ub,
                                    ES    = ((exp(model_i$b)-1))/((exp(model_i$b)+1)),
                                    L     = ((exp(model_i$ci.lb)-1)/(exp(model_i$ci.lb)+1)),
                                    U     = ((exp(model_i$ci.ub)-1)/(exp(model_i$ci.ub)+1)))
    
    
    
    
    #store the data 
    SUBSET[[i]]     <- data_i
    MODEL[[i]]      <- model_i
    result_for_plot <- rbind(result_for_plot,result_for_plot_i)
    
    }}
}

rownames(result_for_plot) <- NULL

ORDER <- as.character(result_for_plot$label_pred)

result_for_plot$label_pred <- factor(result_for_plot$label_pred, ORDER) #sort

#Converting multiple families as multiple
family_split <- strsplit(as.character(db_metafor$Family), ";")

family <- c()
for(i in 1:length(family_split))
  family <- c(family, ifelse(length(family_split[[i]]) > 1, "Multiple", family_split[[i]]) )

db_metafor$Family <- family

# renaming Response group as in the result_for_plot

new_name <- factor(paste(db_metafor$Predictor_Group, db_metafor$Response_Group))
levels(new_name) <- ORDER

db_metafor <- data.frame(db_metafor,new_name)
colnames(db_metafor)

(meta_analysis <- ggplot(data= result_for_plot) +
     geom_hline(yintercept = 0, lty = 2, col = "grey50") +  # add a dotted line at x=1 after flip
     xlab("")+
     ylab("Effect size [r]")+
     geom_jitter(data = db_metafor, aes(x = new_name, y = r, shape = Family, col = Predictor_Group), 
                 size = 1.5, width = 0.2)+
     geom_pointrange(aes(x=label_pred, y=ES, ymin=L, ymax=U, col= label_action), size = 1) + 
     scale_color_manual("Conservation action", values = c("darkmagenta","grey10","darkcyan"))+
     scale_shape_manual("Taxon", values = c(1,2,3))+
     coord_flip() + 
     theme_custom() + theme(legend.position = "right", 
                            legend.direction = "vertical",
        legend.title = element_text(size = 12, face = "bold"),
                            axis.text.y = element_text(face= c("plain","bold","plain","bold","plain")))) # flip coordinates (puts labels on y axis)


#Save figure
pdf(file = "Figure/Meta_analysis.pdf", width = 7, height =5)
meta_analysis
dev.off()

plot_r <- semi_colon_splitter(input1 = db_metafor$Family,
                              input2 = db_metafor$r,
                             names = c("Species","r"))

plot_r$r <- as.numeric(as.character(plot_r$r))

ggplot(data = plot_r, aes(x=1, y = r))+
  geom_boxplot()+
  geom_point(aes(fill = Species), posiion = "jitter",size = 2, alpha = 0.8, pch = 21)+
  scale_fill_manual(values = c("darkmagenta","grey10","darkcyan")) +
    theme_custom() 

# Action by region

COL =  c("grey20","darkorange")

geo_action <- semi_colon_splitter3(input1 = db$Higher_Geography,
                                  input2 = db$Conservation_Action,
                                  input3 = db$Tested_statistically,
                                   names = c("Geography","Action","Test"))

geo_action <- na.omit(geo_action)

bar_1 <- data.frame(table(geo_action$Geography,geo_action$Action,geo_action$Test))

colnames(bar_1) <- c("geo","action","test","N")

bar_1 <- bar_1[bar_1$geo != "Global",] ; bar_1 <- droplevels(bar_1)

bar_1$action <- factor(bar_1$action,levels = 
                       c("Assessment", "Education","Legislation","Monitoring","Prioritization",
                         "Risk assessment",
                         "Decontamination","Eradication","Gating",
                         "Habitat creation","Habitat restoration",
                         "Protected area",  "Regulate access"
                         ))

(bar_p1 <-  ggplot(bar_1, aes(x = action,y = N, fill = test)) +
    facet_wrap( ~ geo, nrow = 2, ncol = 3) +
    geom_bar(stat="identity",position= "stack", color = "grey10")+
    scale_fill_manual("",labels=c("Not tested", "Tested"),values = COL)+
    labs(title=NULL, subtitle = NULL,x=NULL, y = "Frequency")+
    geom_vline(xintercept = 6.5, linetype="dotted",
               color = "grey40", size=0.4)+
    theme_custom()+
    theme(legend.position =  "bottom",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          plot.margin = unit(c(0.2,0.2,0.2,0.2), 'cm'),
          strip.text.x=element_text(color = "grey10", face = "bold", size=12),
          strip.background = element_rect(colour=NA, fill=NA)
          ) 
)

pdf(file = "Figure/Conservation_measures.pdf", width = 12, height = 7)
bar_p1
dev.off()

# Action by region
geo_impact <- semi_colon_splitter(input1 = db$Higher_Geography,
                                  input2 = db$Impact,
                                  names = c("Geography","Impact"))

geo_impact <- na.omit(geo_impact)

bar_2 <- data.frame(table(geo_impact$Geography,geo_impact$Impact))

colnames(bar_2) <- c("geo","Impact","N")

bar_2 <- bar_2[bar_2$geo != "Global",] ; bar_2 <- droplevels(bar_2)

bar_2$Impact <- factor(bar_2$Impact,levels = 
                         c("None identified",
                           "Multiple",
                           "Alien species",
                           "Climate change",
                           "Pathogen",                    
                           "Poaching",
                           "Pollution",
                           "Subterranean habitat change",
                           "Surface habitat change",     
                           "Visitors"))

(bar_p2 <-  ggplot(bar_2, aes(x=Impact,y=N)) +
    facet_wrap( ~ geo, nrow = 2, ncol = 3) +
    #ylim(0,70) +
    geom_bar(stat="identity",position=position_dodge(), color = "grey20")+
    labs(title=NULL, subtitle = NULL,x=NULL, y = "Frequency")+
    theme_custom()+
    theme(legend.position =  "bottom",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          plot.margin = unit(c(0.2,0.2,0.2,0.2), 'cm'),
          strip.text.x=element_text(color = "grey10", face = "bold", size=12),
          strip.background = element_rect(colour=NA, fill=NA)) 
)

pdf(file = "Figure/Threats.pdf", width = 12, height = 7)
bar_p2
dev.off()

ci
# Network -----------------------------------------------------------------

library("ggraph")
library("igraph")        
library("tidygraph")

db_full <- db_full[db_full$System != "Anchialine/Marine",] ; db_full <- droplevels(db_full)
#db_full <- db_full[db_full$Direction_of_effect != "Negative",] ; db_full <- droplevels(db_full)

db_network <- semi_colon_splitter(input1 = db_full$Taxon_Group,
                                  input2 = db_full$Conservation_Action, 
                                   names = c("Taxon","Action"))

db_network$Taxon <- factor(db_network$Taxon, levels = 
                         c("Bats",
                           "Arthropoda",
                           "Microorganisms",
                           "Not specific",       
                           "Other invertebrates",
                           "Other vertebrates","Plants"))

Graph_bipartite <- db_network %>% 
  dplyr::select(Taxon,Action) %>% 
  table() %>% 
  igraph::graph_from_incidence_matrix(directed = TRUE) %>% 
  tidygraph::as_tbl_graph(directed = TRUE)

bipartite_matrix <- as_incidence_matrix(Graph_bipartite)  # Extract the matrix

animal_jaccard       <- ade4::dist.binary(bipartite_matrix, method=1, upper=TRUE, diag = FALSE) # Method #1 is "Jaccard Index"
conservation_jaccard <- ade4::dist.binary(t(bipartite_matrix), method=1, upper=TRUE, diag = FALSE) 

animal_jaccard <- as.matrix(animal_jaccard)   
diag(animal_jaccard)<-0

#animal_jaccard <- ifelse(animal_jaccard>0.5, animal_jaccard, 0)     # Binarize

animal_jaccard_net <- graph_from_adjacency_matrix(animal_jaccard, 
                                         mode = "undirected",
                                         weighted = TRUE) %>% as_tbl_graph()

#Select bat nodes
animal_jaccard_net <- animal_jaccard_net %>% activate(edges) %>% mutate(Bat = edge_is_from(1)) %>% filter(Bat == TRUE)

#Add sample size
N <- data.frame(table(db_network$Taxon)) ; colnames(N) <- c("name","N")
animal_jaccard_net <- animal_jaccard_net %>% activate(nodes) %>%  dplyr::left_join(N,by = "name")

plot <- ggraph::ggraph(animal_jaccard_net,  layout_with_kk(animal_jaccard_net)) +
  geom_edge_arc(aes(width=weight), color = "orange", alpha = 0.8, strength = .1) +
  
  geom_node_point(aes(size = N), fill="blue", alpha = 1, 
                  shape = 21) + 
  
  scale_edge_width_continuous("Jaccard similarity",range=c(0.2,3))+
  scale_size_continuous("Sample size")+
  
  geom_node_text(aes(label = name), size=8, color="gray10", repel=TRUE, force = 13) +
   theme_void() + theme(legend.position = "top",legend.direction = "horizontal") + coord_fixed()

pdf(file = "Figure/Network.pdf", width = 12, height = 7)
plot
dev.off()

## APPNTI:::





# Conservation actions by Family ------------------------------------------

levels(db$Higher_Geography)
family_action <- semi_colon_splitter(input1 = db$Family,
                            input2 = db$Conservation_Action, 
                            names = c("Family","Action"))


family_action <- na.omit(family_action)

Graph_bipartite <- family_action %>% 
  dplyr::select(Family,Action) %>% 
  table() %>% 
  igraph::graph_from_incidence_matrix(directed = FALSE) %>% 
  tidygraph::as_tbl_graph(directed = FALSE)

# Graph_bipartite <- Graph_bipartite %>% tidygraph::activate(nodes) %>% 
#   left_join(rbind(data.frame(table(db$Taxon)),
#                   data.frame(table(db$Conservation_Action))), by = c("name" = "Var1"))


# Collapse it into an unipartite 
Graph_unipartite_full <- igraph::bipartite_projection(Graph_bipartite)

# Takes the unipartite project graph
Graph_unipartite <- Graph_unipartite_full$proj1  %>% as_tbl_graph(directed = FALSE) %>% 
  activate(edges) %>% #%>% mutate(weight = 1) 
  igraph::simplify(edge.attr.comb = "sum") %>% 
  as_tbl_graph


########################
# Plotting the network #
########################

#SpatialLayout
Graph_plot <- Graph_unipartite %>% 
  igraph::simplify(edge.attr.comb = "sum") 

# Graph_plot_bat <- to_subgraph(Graph_unipartite, to %in% c(9) | from %in% c(9), subset_by = "edges")$subgraph %>% 
#   igraph::simplify(edge.attr.comb = "sum") 

library("igraph")

ggraph::ggraph(Graph_plot,  layout_with_kk(Graph_plot)) +
   geom_edge_arc(aes(width=weight) , strength = .1,
                alpha = 0.1) +
  
  geom_node_point(fill="grey30", alpha = .8, 
                  shape = 21) + 
  geom_node_text(aes(label = name), size=4, color="gray10", repel=TRUE, force = 10) +
  scale_color_gradient2("Connection strength",
                        low= "#f7fcfd",
                        mid = "#8c96c6",
                        high = "#4d004b") +
   theme_void() + theme(legend.position = "bottom",legend.direction = "horizontal") + coord_fixed()

Graph_bipartite %>%  igraph::simplify(edge.attr.comb = "sum") %>% 
  ggraph::ggraph(.,  layout = "bipartite") +
  geom_edge_link0(edge_colour = "grey66")+
  # geom_node_point(alpha = .8, 
  #                 aes(size=Freq, fill=type), shape = 21) +
  scale_fill_manual(values=c("blue","red"))+
  scale_colour_manual(values=c("blue","red"))+
  geom_node_text(aes(label = name, color=type), size=4, repel=TRUE, force = 10) +
  
  theme_void() + theme(legend.position = "left",legend.direction = "horizontal") 

bipartite_matrix <- as_incidence_matrix(Graph_bipartite)  # Extract the matrix

animal_jaccard <- dist.binary(bipartite_matrix, method=1, upper=TRUE, diag = FALSE) # Method #1 is "Jaccard Index"
conservation_jaccard <- dist.binary(t(bipartite_matrix), method=1, upper=TRUE, diag = FALSE) 

animal_jaccard <- as.matrix(animal_jaccard)   
diag(animal_jaccard)<-0

# women_jaccard          # Look at the matrix before you binarize
animal_jaccard <- ifelse(animal_jaccard>0.7, 1, 0)     # Binarize

# jaccard_women      # Take a look at the matrix if you like.

animal_jaccard <- graph_from_adjacency_matrix(animal_jaccard,    # Create an igraph network
                                              mode = "undirected")
plot(animal_jaccard)

# Network test ------------------------------------------------------------

db_full <- db_full[db_full$System != "Anchialine/Marine",] ; db_full <- droplevels(db_full)


levels(db_full$Taxon_Group)[c(2,3,4)] <- "Arthropoda"

db <- db[db$Direction_of_effect != "Negative",] ; db <- droplevels(db)

Graph_bipartite <- db_full %>% 
  dplyr::select(Taxon_Group,Conservation_Action) %>% 
  table() %>% 
  igraph::graph_from_incidence_matrix(directed = FALSE) %>% 
  tidygraph::as_tbl_graph(directed = FALSE) 

Graph_bipartite <- Graph_bipartite %>% tidygraph::activate(nodes) %>% 
  left_join(rbind(data.frame(table(db$Taxon)),
                  data.frame(table(db$Conservation_Action))), by = c("name" = "Var1"))

# Collapse it into an unipartite 
Graph_unipartite_full <- igraph::bipartite_projection(Graph_bipartite)

# Takes the unipartite project graph
Graph_unipartite <- Graph_unipartite_full$proj1  %>% as_tbl_graph(directed = FALSE) %>% 
  activate(edges) %>% #%>% mutate(weight = 1) 
  igraph::simplify(edge.attr.comb = "sum") %>% 
  as_tbl_graph

########################
# Plotting the network #
########################

#SpatialLayout
Graph_plot <- Graph_unipartite %>% 
    igraph::simplify(edge.attr.comb = "sum") 

Graph_plot_bat <- to_subgraph(Graph_unipartite, to %in% c(2) | from %in% c(2), subset_by = "edges")$subgraph %>% 
  igraph::simplify(edge.attr.comb = "sum") 

ggraph::ggraph(Graph_plot_bat,  layout_with_kk(Graph_plot)) +
    #geom_edge_density(fill="orange", alpha=1) +
    geom_edge_arc(aes(width=weight) , strength = .1,
                  alpha = 0.1) +
    
    geom_node_point(fill="grey30", alpha = .8, 
                    aes(size=Freq), shape = 21) + 
    geom_node_text(aes(label = name), size=4, color="gray10", repel=TRUE, force = 10) +
    scale_color_gradient2("Connection strength",
                          low= "#f7fcfd",
                          mid = "#8c96c6",
                          high = "#4d004b") +
    #scale_fill_manual(values = c("blue", "orange", "pink","purple", "grey15")) +
    theme_void() + theme(legend.position = "bottom",legend.direction = "horizontal") + coord_fixed()

Graph_bipartite %>%  igraph::simplify(edge.attr.comb = "sum") %>% 
  ggraph::ggraph(.,  layout_with_kk(.)) +
   geom_edge_arc(
                 strength = .1,
                alpha = 0.3) +
   geom_node_point(alpha = .8, 
                  aes(size=Freq, fill=type), shape = 21) +
  scale_fill_manual(values=c("blue","red"))+
  scale_colour_manual(values=c("blue","red"))+
  geom_node_text(aes(label = name, color=type), size=4, repel=TRUE, force = 10) +
  
  theme_void() + theme(legend.position = "bottom",legend.direction = "horizontal") + coord_fixed()

library(ade4) # If you have not already done so

bipartite_matrix <- as_incidence_matrix(Graph_bipartite)  # Extract the matrix

animal_jaccard <- dist.binary(bipartite_matrix, method=1, upper=TRUE, diag = FALSE) # Method #1 is "Jaccard Index"
conservation_jaccard <- dist.binary(t(bipartite_matrix), method=1, upper=TRUE, diag = FALSE) 

animal_jaccard <- as.matrix(animal_jaccard)   
diag(animal_jaccard)<-0

# women_jaccard          # Look at the matrix before you binarize
animal_jaccard <- ifelse(animal_jaccard>0.8, 1, 0)     # Binarize

# jaccard_women      # Take a look at the matrix if you like.

animal_jaccard <- graph_from_adjacency_matrix(animal_jaccard,    # Create an igraph network
                                          mode = "undirected")
plot(animal_jaccard)


