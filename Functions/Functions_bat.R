###############################################################

## Let's not wing it: Practical conservation of subterranean-roosting bats
## Meierhofer M., et al.

## ------------------------------------------------------------------------
# 'R script to reproduce the analyses'
## ------------------------------------------------------------------------

# Analysis performed with R (v. R 4.1.0) and R studio (v. 1.4.1103)
# Authors (code): Stefano Mammola & Melissa B. Meierhofer

###############################################################

# Loading useful functions ------------------------------------------------

# Custom function to split columns having semicolon as a separator
semi_colon_splitter <- function(input1, input2, names = c("input1","input2")){
  
  df        <- data.frame(input1,input2)  
  df$input1 <- as.factor(df$input1)
  df$input2 <- as.factor(df$input2)
  
  to_separate <- levels(df$input1)[grepl(";", levels(df$input1))]
  
  df_all <- df[df$input1 %in% to_separate ,]
  df     <- df[!df$input1 %in% to_separate,]
  df$input1 <- droplevels(df$input1)
  
  df_all$input1 <- as.character(df_all$input1)
  
  for(i in nrow(df_all)) {
    
    df_i <- df_all[i,]
    split   <- strsplit(df_all$input1, ";")[[1]]
    split   <- trimws(split, which = c("both"))
    
    df <- rbind(df,data.frame(input1  = split,
                              input2  = rep(df_i$input2, length(split))))
    
  }
  
  colnames(df) <- names
  return(df)
}

# Custom function to split columns having semicolon as a separator (3 columns version)
semi_colon_splitter3 <- function(input1, input2,input3, names = c("input1","input2","input3")){
  
  df        <- data.frame(input1,input2,input3)  
  df$input1 <- as.factor(df$input1)
  df$input2 <- as.factor(df$input2)
  df$input3 <- as.factor(df$input3)
  
  to_separate <- levels(df$input1)[grepl(";", levels(df$input1))]
  
  df_all <- df[df$input1 %in% to_separate ,]
  df     <- df[!df$input1 %in% to_separate,]
  df$input1 <- droplevels(df$input1)
  
  df_all$input1 <- as.character(df_all$input1)
  
  for(i in nrow(df_all)) {
    
    df_i <- df_all[i,]
    split   <- strsplit(df_all$input1, ";")[[1]]
    split   <- trimws(split, which = c("both"))
    
    df <- rbind(df,data.frame(input1  = split,
                              input2  = rep(df_i$input2, length(split)),
                              input3  = rep(df_i$input3, length(split))))
    
  }
  
  colnames(df) <- names
  return(df)
}

# Custom function to get standard error
SE <- function(x) sd(x) / sqrt( length(x) )

# Parameters for plots ----------------------------------------------------

# Plot style (ggplot2)

theme_custom <- function(){
  theme_bw() +
    theme(#text = element_text(family = "Arial"),
      axis.text = element_text(size = 10), 
      axis.title = element_text(size = 12),
      axis.line.x = element_line(color="black"), 
      axis.line.y = element_line(color="black"),
      #panel.border = element_blank(),
      panel.grid.major.x = element_blank(),                                          
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),  
      plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
      plot.title = element_text(size = 18, vjust = 1, hjust = 0),
      legend.text = element_text(size = 12),          
      legend.title = element_blank(),                              
      legend.position = "top", 
      legend.direction = "horizontal",
      legend.key = element_blank(),
      legend.background = element_rect(color = "black", 
                                       fill = "transparent", 
                                       size = 2, linetype = "blank"))
}

