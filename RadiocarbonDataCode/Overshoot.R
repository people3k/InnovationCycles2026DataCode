##Overhoot plotting for fit models
library(tidyverse)
library(cowplot)
library(viridis)
library(ggplot2)

##Plot theme
my_theme <- theme_bw() +
  theme(
    axis.text.x  = element_text(size=20, colour="black"),
    axis.text.y  = element_text(size=20, colour="black"),
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size=18),
    plot.title   = element_text(size=18, face="bold"),
    legend.title = element_text(size=14),
    legend.text  = element_text(size=12),
    strip.text.x = element_text(size = 18, face = "bold"),
    strip.text.y = element_text(size = 18, face = "bold"),
    
    # optional: facet background cleaner
    strip.background = element_rect(fill = "grey90", colour = "black")
  )



##S. Colorado Plat.===================
###Load data
na_rec <- read.csv("Inno_Pop_Fitting/Plotting/NARecessions.csv")

fits <- read.csv("Inno_Pop_Fitting/Plotting/NullModel_Envelopes_AllObjectives.csv")

box <- read.csv("Inno_Pop_Fitting/Plotting/Stablity_Robust_Region2.csv")

kde <- read.csv("Inno_Pop_Fitting/Plotting/AllKDE_Hybrid_Corr_Features2.csv")


# ============================================================
# SUBSET DATA
# ============================================================

# Stability table
box <- subset(box, fit_type== "hybrid")

# Create classification variable
box$Classification <- apply(
  box[, c("eig1_real","eig1_imag",
          "eig2_real","eig2_imag",
          "eig3_real","eig3_imag")],
  1,
  function(x) {
    
    real_vals <- c(x[1], x[3], x[5])
    imag_vals <- c(x[2], x[4], x[6])
    
    # Any complex unstable eigenvalues
    if (any(real_vals > 0 & imag_vals != 0)) {
      return("complex unstable")
      
      # Any imaginary component present
    } else if (any(imag_vals != 0)) {
      return("complex stable")
      
      # Otherwise purely real
    } else {
      return("real")
    }
  }
)

# Calculate overshoot values for each eigen pair
ov1 <- log10((abs(box$eig1_imag) / 
                (abs(box$eig1_real) + 0.001)) + 0.001)

ov2 <- log10((abs(box$eig2_imag) / 
                (abs(box$eig2_real) + 0.001)) + 0.001)

ov3 <- log10((abs(box$eig3_imag) / 
                (abs(box$eig3_real) + 0.001)) + 0.001)

# Select highest overshoot across eig1–eig3
box$Overshoot <- pmax(ov1, ov2, ov3)



# Fits table
fits <-subset(fits, objective=="shape")

uplandsw1<-subset(na_rec, region_id=="S. Colorado Plat.")
uplandswNull<-subset(fits, region_id=="S. Colorado Plat.")
uplandswFit<-subset(kde, region_id=="S. Colorado Plat." )
uplandswEq<-subset(box, region_id=="S. Colorado Plat." )

###Add calBP to the data sets for plotting
calBP<-uplandswNull$time*-1
uplandswNull$calBP <- uplandswNull$time * -1
uplandswFit<-cbind(uplandswFit, calBP)

# Create lookup vectors
class_lookup <- setNames(
  uplandswEq$Classification,
  uplandswEq$segment_index
)

overshoot_lookup <- setNames(
  uplandswEq$Overshoot,
  uplandswEq$segment_index
)

# Transfer values into uplandswFit
uplandswFit$Classification <- class_lookup[
  as.character(uplandswFit$segment_index)
]

uplandswFit$Overshoot <- overshoot_lookup[
  as.character(uplandswFit$segment_index)
]

###Plot------------------------------------------
library(ggplot2)
library(dplyr)
library(viridis)

segshade <- uplandswFit %>%
  group_by(segment_index, Overshoot) %>%
  summarise(
    xmin = min(calBP, na.rm = TRUE),
    xmax = max(calBP, na.rm = TRUE),
    .groups = "drop"
  )
# Segment boundary positions
seg_bounds <- sort(unique(c(segshade$xmax, segshade$xmax)))

vir_cols <- viridis(10)

Cpt <- ggplot() +
  
  # ==========================================================
# Segment shading
# ==========================================================
geom_rect(
  data = segshade,
  aes(
    xmin = xmin,
    xmax = xmax,
    ymin = -Inf,
    ymax = Inf,
    fill = Overshoot
  ),
  alpha = 0.18,
  inherit.aes = FALSE
) +
  
  # ----------------------------------------------------------
# Segment boundary lines
# ----------------------------------------------------------
geom_vline(
  xintercept = seg_bounds,
  linetype = "dashed",
  color = "grey40",
  linewidth = 1,
  alpha = 0.8
) +
  
  # ==========================================================
# KDE envelope
# ==========================================================
geom_ribbon(
  data = uplandsw1,
  aes(
    x = calBP,
    ymin = lo * 100,
    ymax = hi * 100
  ),
  fill = vir_cols[1],
  alpha = 0.5
) +
  
  # ==========================================================
# Null envelope
# ==========================================================
geom_ribbon(
  data = uplandswNull,
  aes(
    x = calBP,
    ymin = Y_min,
    ymax = Y_max
  ),
  fill = vir_cols[9],
  alpha = 0.35
) +
  
  # ==========================================================
# KDE curve
# ==========================================================
geom_path(
  data = uplandsw1,
  aes(
    x = calBP,
    y = MKDE * 100,
    color = "KDE"
  ),
  linewidth = 1.5
) +
  
  # ==========================================================
# Fitted curve
# ==========================================================
geom_path(
  data = uplandswFit,
  aes(
    x = calBP,
    y = YPredicted,
    color = "Fitted"
  ),
  linewidth = 1.75
) +
  
  # ==========================================================
# Null mean
# ==========================================================
geom_path(
  data = uplandswNull,
  aes(
    x = calBP,
    y = Y_mean,
    color = "Null mean"
  ),
  linewidth = 1
) +
  
  # ==========================================================
# Axes
# ==========================================================
scale_x_reverse(
  limits = c(4000, 200)
) +
  
  # ==========================================================
# Line colors
# ==========================================================
scale_color_manual(
  name = "Curves",
  values = c(
    "KDE" = vir_cols[1],
    "Fitted" = vir_cols[5],
    "Null mean" = vir_cols[9]
  )
) +
  
  # ==========================================================
# Segment colors
# ==========================================================
scale_fill_viridis(
  name = "Overshoot",
  option = "D",
  direction=1
) +
  
  theme_bw() +
  my_theme +
  
  labs(
    x = "Years cal BP",
    y = "Mean KDE density",
    title = "E. S. Colorado Plat. Output (Y) vs. Time"
  )

Cpt


###Plot for Chihuahua Desert===================================
##Chi Desert
#===============================================================
na_rec <- read.csv("Inno_Pop_Fitting/Plotting/NARecessions.csv")

fits <- read.csv("Inno_Pop_Fitting/Plotting/NullModel_Envelopes_AllObjectives.csv")

box <- read.csv("Inno_Pop_Fitting/Plotting/Stablity_Robust_Region2.csv")

kde <- read.csv("Inno_Pop_Fitting/Plotting/AllKDE_Hybrid_Corr_Features2.csv")


# ============================================================
# SUBSET DATA
# ============================================================

# Stability table
box <- subset(box, fit_type== "hybrid")

# Create classification variable
box$Classification <- apply(
  box[, c("eig1_real","eig1_imag",
          "eig2_real","eig2_imag",
          "eig3_real","eig3_imag")],
  1,
  function(x) {
    
    real_vals <- c(x[1], x[3], x[5])
    imag_vals <- c(x[2], x[4], x[6])
    
    # Any complex unstable eigenvalues
    if (any(real_vals > 0 & imag_vals != 0)) {
      return("complex unstable")
      
      # Any imaginary component present
    } else if (any(imag_vals != 0)) {
      return("complex stable")
      
      # Otherwise purely real
    } else {
      return("real")
    }
  }
)

# Calculate overshoot values for each eigen pair
ov1 <- log10((abs(box$eig1_imag) / 
                (abs(box$eig1_real) + 0.001)) + 0.001)

ov2 <- log10((abs(box$eig2_imag) / 
                (abs(box$eig2_real) + 0.001)) + 0.001)

ov3 <- log10((abs(box$eig3_imag) / 
                (abs(box$eig3_real) + 0.001)) + 0.001)

# Select highest overshoot across eig1–eig3
box$Overshoot <- pmax(ov1, ov2, ov3)


# Fits table
fits <-subset(fits, objective=="shape")

chides1<-subset(na_rec, region_id=="Chihuahua Desert")
chidesNull<-subset(fits, region_id=="Chihuahua Desert")
chidesFit<-subset(kde, region_id=="Chihuahua Desert" )
chidesEq<-subset(box, region_id=="Chihuahua Desert" )

###Add calBP to the data sets for plotting
calBP<-chidesNull$time*-1
chidesNull$calBP <- chidesNull$time * -1
chidesFit<-cbind(chidesFit, calBP)

# Create lookup vectors
class_lookup <- setNames(
  chidesEq$Classification,
  chidesEq$segment_index
)

overshoot_lookup <- setNames(
  chidesEq$Overshoot,
  chidesEq$segment_index
)

# Transfer values into uplandswFit
chidesFit$Classification <- class_lookup[
  as.character(chidesFit$segment_index)
]

chidesFit$Overshoot <- overshoot_lookup[
  as.character(chidesFit$segment_index)
]

###Plot------------------------------------------
library(ggplot2)
library(dplyr)
library(viridis)

segshade <- chidesFit %>%
  group_by(segment_index, Overshoot) %>%
  summarise(
    xmin = min(calBP, na.rm = TRUE),
    xmax = max(calBP, na.rm = TRUE),
    .groups = "drop"
  )
# Segment boundary positions
seg_bounds <- sort(unique(c(segshade$xmax, segshade$xmax)))

vir_cols <- viridis(10)

chidesert <- ggplot() +
  
  # ==========================================================
# Segment shading
# ==========================================================
geom_rect(
  data = segshade,
  aes(
    xmin = xmin,
    xmax = xmax,
    ymin = -Inf,
    ymax = Inf,
    fill = Overshoot
  ),
  alpha = 0.18,
  inherit.aes = FALSE
) +
  
  # ----------------------------------------------------------
# Segment boundary lines
# ----------------------------------------------------------
geom_vline(
  xintercept = seg_bounds,
  linetype = "dashed",
  color = "grey40",
  linewidth = 1,
  alpha = 0.8
) +
  
  # ==========================================================
# KDE envelope
# ==========================================================
geom_ribbon(
  data = chides1,
  aes(
    x = calBP,
    ymin = lo * 100,
    ymax = hi * 100
  ),
  fill = vir_cols[1],
  alpha = 0.5
) +
  
  # ==========================================================
# Null envelope
# ==========================================================
geom_ribbon(
  data = chidesNull, aes(x = calBP, ymin = Y_min,
                         ymax = Y_max ), fill = vir_cols[9],  alpha = 0.35) +
  
  # ==========================================================
# KDE curve
# ==========================================================
geom_path(
  data = chides1,
  aes(
    x = calBP,
    y = MKDE * 100,
    color = "KDE"
  ),
  linewidth = 1.5
) +
  
  # ==========================================================
# Fitted curve
# ==========================================================
geom_path(
  data = chidesFit,
  aes(
    x = calBP,
    y = YPredicted,
    color = "Fitted"
  ),
  linewidth = 1.75
) +
  
  # ==========================================================
# Null mean
# ==========================================================
geom_path(
  data = chidesNull,
  aes(
    x = calBP,
    y = Y_mean,
    color = "Null mean"
  ),
  linewidth = 1
) +
  
  # ==========================================================
# Axes
# ==========================================================
scale_x_reverse(
  limits = c(4000, 200)
) +
  
  # ==========================================================
# Line colors
# ==========================================================
scale_color_manual(
  name = "Curves",
  values = c(
    "KDE" = vir_cols[1],
    "Fitted" = vir_cols[5],
    "Null mean" = vir_cols[9]
  )
) +
  
  # ==========================================================
# Segment colors
# ==========================================================
scale_fill_viridis(
  name = "Overshoot",
  option = "D",
  direction=1
) +
  
  theme_bw() +
  my_theme +
  
  labs(
    x = "Years cal BP",
    y = "Mean KDE density",
    title = "C. Chihuahua Desert Output (Y) vs. Time"
  )

chidesert


###Plot for Sonoran Desert===================================
##Sonoran Desert
#===============================================================
na_rec <- read.csv("Inno_Pop_Fitting/Plotting/NARecessions.csv")

fits <- read.csv("Inno_Pop_Fitting/Plotting/NullModel_Envelopes_AllObjectives.csv")

box <- read.csv("Inno_Pop_Fitting/Plotting/Stablity_Robust_Region2.csv")

kde <- read.csv("Inno_Pop_Fitting/Plotting/AllKDE_Hybrid_Corr_Features2.csv")


# ============================================================
# SUBSET DATA
# ============================================================

# Stability table
box <- subset(box, fit_type== "hybrid")

# Create classification variable
box$Classification <- apply(
  box[, c("eig1_real","eig1_imag",
          "eig2_real","eig2_imag",
          "eig3_real","eig3_imag")],
  1,
  function(x) {
    
    real_vals <- c(x[1], x[3], x[5])
    imag_vals <- c(x[2], x[4], x[6])
    
    # Any complex unstable eigenvalues
    if (any(real_vals > 0 & imag_vals != 0)) {
      return("complex unstable")
      
      # Any imaginary component present
    } else if (any(imag_vals != 0)) {
      return("complex stable")
      
      # Otherwise purely real
    } else {
      return("real")
    }
  }
)

# Calculate overshoot values for each eigen pair
ov1 <- log10((abs(box$eig1_imag) / 
                (abs(box$eig1_real) + 0.001)) + 0.001)

ov2 <- log10((abs(box$eig2_imag) / 
                (abs(box$eig2_real) + 0.001)) + 0.001)

ov3 <- log10((abs(box$eig3_imag) / 
                (abs(box$eig3_real) + 0.001)) + 0.001)

# Select highest overshoot across eig1–eig3
box$Overshoot <- pmax(ov1, ov2, ov3)


# Fits table
fits <-subset(fits, objective=="shape")

sondes1<-subset(na_rec, region_id=="Sonoran Desert")
sondesNull<-subset(fits, region_id=="Sonoran Desert")
sondesFit<-subset(kde, region_id=="Sonoran Desert" )
sondesEq<-subset(box, region_id=="Sonoran Desert" )

###Add calBP to the data sets for plotting
calBP<-sondesNull$time*-1
sondesNull$calBP <- sondesNull$time * -1
sondesFit<-cbind(sondesFit, calBP)

# Create lookup vectors
class_lookup <- setNames(
  sondesEq$Classification,
  sondesEq$segment_index
)

overshoot_lookup <- setNames(
  sondesEq$Overshoot,
  sondesEq$segment_index
)

# Transfer values into uplandswFit
sondesFit$Classification <- class_lookup[
  as.character(sondesFit$segment_index)
]

sondesFit$Overshoot <- overshoot_lookup[
  as.character(sondesFit$segment_index)
]

###Plot------------------------------------------
library(ggplot2)
library(dplyr)
library(viridis)

segshade <- sondesFit %>%
  group_by(segment_index, Overshoot) %>%
  summarise(
    xmin = min(calBP, na.rm = TRUE),
    xmax = max(calBP, na.rm = TRUE),
    .groups = "drop"
  )
# Segment boundary positions
seg_bounds <- sort(unique(c(segshade$xmax, segshade$xmax)))

vir_cols <- viridis(10)

sondesert <- ggplot() +
  
  # ==========================================================
# Segment shading
# ==========================================================
geom_rect(
  data = segshade,
  aes(
    xmin = xmin,
    xmax = xmax,
    ymin = -Inf,
    ymax = Inf,
    fill = Overshoot
  ),
  alpha = 0.18,
  inherit.aes = FALSE
) +
  
  # ----------------------------------------------------------
# Segment boundary lines
# ----------------------------------------------------------
geom_vline(
  xintercept = seg_bounds,
  linetype = "dashed",
  color = "grey40",
  linewidth = 1,
  alpha = 0.8
) +
  
  # ==========================================================
# KDE envelope
# ==========================================================
geom_ribbon(
  data = sondes1,
  aes(
    x = calBP,
    ymin = lo * 100,
    ymax = hi * 100
  ),
  fill = vir_cols[1],
  alpha = 0.5
) +
  
  # ==========================================================
# Null envelope
# ==========================================================
geom_ribbon(
  data = sondesNull, aes(x = calBP, ymin = Y_min,
                         ymax = Y_max ), fill = vir_cols[9],  alpha = 0.35) +
  
  # ==========================================================
# KDE curve
# ==========================================================
geom_path(
  data = sondes1,
  aes(
    x = calBP,
    y = MKDE * 100,
    color = "KDE"
  ),
  linewidth = 1.5
) +
  
  # ==========================================================
# Fitted curve
# ==========================================================
geom_path(
  data = sondesFit,
  aes(
    x = calBP,
    y = YPredicted,
    color = "Fitted"
  ),
  linewidth = 1.75
) +
  
  # ==========================================================
# Null mean
# ==========================================================
geom_path(
  data = sondesNull,
  aes(
    x = calBP,
    y = Y_mean,
    color = "Null mean"
  ),
  linewidth = 1
) +
  
  # ==========================================================
# Axes
# ==========================================================
scale_x_reverse(
  limits = c(4000, 200)
) +
  
  # ==========================================================
# Line colors
# ==========================================================
scale_color_manual(
  name = "Curves",
  values = c(
    "KDE" = vir_cols[1],
    "Fitted" = vir_cols[5],
    "Null mean" = vir_cols[9]
  )
) +
  
  # ==========================================================
# Segment colors
# ==========================================================
scale_fill_viridis(
  name = "Overshoot",
  option = "D",
  direction=1
) +
  
  theme_bw() +
  my_theme +
  
  labs(
    x = "Years cal BP",
    y = "Mean KDE density",
    title = "D. Sonoran Desert Output (Y) vs. Time"
  )

sondesert

###Plot for N. Colorado Plat.===================================
##N. Colorado Plat.
#===============================================================
na_rec <- read.csv("Inno_Pop_Fitting/Plotting/NARecessions.csv")

fits <- read.csv("Inno_Pop_Fitting/Plotting/NullModel_Envelopes_AllObjectives.csv")

box <- read.csv("Inno_Pop_Fitting/Plotting/Stablity_Robust_Region2.csv")

kde <- read.csv("Inno_Pop_Fitting/Plotting/AllKDE_Hybrid_Corr_Features2.csv")


# ============================================================
# SUBSET DATA
# ============================================================

# Stability table
box <- subset(box, fit_type== "hybrid")

# Create classification variable
box$Classification <- apply(
  box[, c("eig1_real","eig1_imag",
          "eig2_real","eig2_imag",
          "eig3_real","eig3_imag")],
  1,
  function(x) {
    
    real_vals <- c(x[1], x[3], x[5])
    imag_vals <- c(x[2], x[4], x[6])
    
    # Any complex unstable eigenvalues
    if (any(real_vals > 0 & imag_vals != 0)) {
      return("complex unstable")
      
      # Any imaginary component present
    } else if (any(imag_vals != 0)) {
      return("complex stable")
      
      # Otherwise purely real
    } else {
      return("real")
    }
  }
)

# Calculate overshoot values for each eigen pair
ov1 <- log10((abs(box$eig1_imag) / 
                (abs(box$eig1_real) + 0.001)) + 0.001)

ov2 <- log10((abs(box$eig2_imag) / 
                (abs(box$eig2_real) + 0.001)) + 0.001)

ov3 <- log10((abs(box$eig3_imag) / 
                (abs(box$eig3_real) + 0.001)) + 0.001)

# Select highest overshoot across eig1–eig3
box$Overshoot <- pmax(ov1, ov2, ov3)


# Fits table
fits <-subset(fits, objective=="shape")

ncol1<-subset(na_rec, region_id=="N. Colorado Plat.")
ncolNull<-subset(fits, region_id=="N. Colorado Plat.")
ncolFit<-subset(kde, region_id=="N. Colorado Plat." )
ncolEq<-subset(box, region_id=="N. Colorado Plat." )

###Add calBP to the data sets for plotting
calBP<-ncolNull$time*-1
ncolNull$calBP <- ncolNull$time * -1
ncolFit<-cbind(ncolFit, calBP)

# Create lookup vectors
class_lookup <- setNames(
  ncolEq$Classification,
  ncolEq$segment_index
)

overshoot_lookup <- setNames(
  ncolEq$Overshoot,
  ncolEq$segment_index
)

# Transfer values into uplandswFit
ncolFit$Classification <- class_lookup[
  as.character(ncolFit$segment_index)
]

ncolFit$Overshoot <- overshoot_lookup[
  as.character(ncolFit$segment_index)
]

###Plot------------------------------------------
library(ggplot2)
library(dplyr)
library(viridis)

segshade <- ncolFit %>%
  group_by(segment_index, Overshoot) %>%
  summarise(
    xmin = min(calBP, na.rm = TRUE),
    xmax = max(calBP, na.rm = TRUE),
    .groups = "drop"
  )
# Segment boundary positions
seg_bounds <- sort(unique(c(segshade$xmax, segshade$xmax)))

vir_cols <- viridis(10)

ncolplat <- ggplot() +
  
  # ==========================================================
# Segment shading
# ==========================================================
geom_rect(
  data = segshade,
  aes(
    xmin = xmin,
    xmax = xmax,
    ymin = -Inf,
    ymax = Inf,
    fill = Overshoot
  ),
  alpha = 0.18,
  inherit.aes = FALSE
) +
  
  # ----------------------------------------------------------
# Segment boundary lines
# ----------------------------------------------------------
geom_vline(
  xintercept = seg_bounds,
  linetype = "dashed",
  color = "grey40",
  linewidth = 1,
  alpha = 0.8
) +
  
  # ==========================================================
# KDE envelope
# ==========================================================
geom_ribbon(
  data = ncol1,
  aes(
    x = calBP,
    ymin = lo * 100,
    ymax = hi * 100
  ),
  fill = vir_cols[1],
  alpha = 0.5
) +
  
  # ==========================================================
# Null envelope
# ==========================================================
geom_ribbon(
  data = ncolNull, aes(x = calBP, ymin = Y_min,
                       ymax = Y_max ), fill = vir_cols[9],  alpha = 0.35) +
  
  # ==========================================================
# KDE curve
# ==========================================================
geom_path(
  data = ncol1,
  aes(
    x = calBP,
    y = MKDE * 100,
    color = "KDE"
  ),
  linewidth = 1.5
) +
  
  # ==========================================================
# Fitted curve
# ==========================================================
geom_path(
  data = ncolFit,
  aes(
    x = calBP,
    y = YPredicted,
    color = "Fitted"
  ),
  linewidth = 1.75
) +
  
  # ==========================================================
# Null mean
# ==========================================================
geom_path(
  data = ncolNull,
  aes(
    x = calBP,
    y = Y_mean,
    color = "Null mean"
  ),
  linewidth = 1
) +
  
  # ==========================================================
# Axes
# ==========================================================
scale_x_reverse(
  limits = c(2700, 200)
) +
  
  # ==========================================================
# Line colors
# ==========================================================
scale_color_manual(
  name = "Curves",
  values = c(
    "KDE" = vir_cols[1],
    "Fitted" = vir_cols[5],
    "Null mean" = vir_cols[9]
  )
) +
  
  # ==========================================================
# Segment colors
# ==========================================================
scale_fill_viridis(
  name = "Overshoot",
  option = "D",
  direction=1
) +
  
  theme_bw() +
  my_theme +
  
  labs(
    x = "Years cal BP",
    y = "Mean KDE density",
    title = "F. N. Colorado Plat. Output (Y) vs. Time"
  )

ncolplat


###Plot for US Southwest.===================================
##US Southwest Region
#===============================================================
na_rec <- read.csv("Inno_Pop_Fitting/Plotting/NARecessions.csv")

fits <- read.csv("Inno_Pop_Fitting/Plotting/NullModel_Envelopes_AllObjectives.csv")

box <- read.csv("Inno_Pop_Fitting/Plotting/Stablity_Robust_Region2.csv")

kde <- read.csv("Inno_Pop_Fitting/Plotting/AllKDE_Hybrid_Corr_Features2.csv")


# ============================================================
# SUBSET DATA
# ============================================================

# Stability table
box <- subset(box, fit_type== "hybrid")

# Create classification variable
box$Classification <- apply(
  box[, c("eig1_real","eig1_imag",
          "eig2_real","eig2_imag",
          "eig3_real","eig3_imag")],
  1,
  function(x) {
    
    real_vals <- c(x[1], x[3], x[5])
    imag_vals <- c(x[2], x[4], x[6])
    
    # Any complex unstable eigenvalues
    if (any(real_vals > 0 & imag_vals != 0)) {
      return("complex unstable")
      
      # Any imaginary component present
    } else if (any(imag_vals != 0)) {
      return("complex stable")
      
      # Otherwise purely real
    } else {
      return("real")
    }
  }
)

# Calculate overshoot values for each eigen pair
ov1 <- log10((abs(box$eig1_imag) / 
                (abs(box$eig1_real) + 0.001)) + 0.001)

ov2 <- log10((abs(box$eig2_imag) / 
                (abs(box$eig2_real) + 0.001)) + 0.001)

ov3 <- log10((abs(box$eig3_imag) / 
                (abs(box$eig3_real) + 0.001)) + 0.001)

# Select highest overshoot across eig1–eig3
box$Overshoot <- pmax(ov1, ov2, ov3)


# Fits table
fits <-subset(fits, objective=="shape")

swr1<-subset(na_rec, region_id=="Southwest US")
swrNull<-subset(fits, region_id=="Southwest US")
swrFit<-subset(kde, region_id=="Southwest US" )
swrEq<-subset(box, region_id=="Southwest US" )

###Add calBP to the data sets for plotting
calBP<-swrNull$time*-1
swrNull$calBP <- swrNull$time * -1
swrFit<-cbind(swrFit, calBP)

# Create lookup vectors
class_lookup <- setNames(
  swrEq$Classification,
  swrEq$segment_index
)

overshoot_lookup <- setNames(
  swrEq$Overshoot,
  swrEq$segment_index
)

# Transfer values into uplandswFit
swrFit$Classification <- class_lookup[
  as.character(swrFit$segment_index)
]

swrFit$Overshoot <- overshoot_lookup[
  as.character(swrFit$segment_index)
]

###Plot------------------------------------------
library(ggplot2)
library(dplyr)
library(viridis)

segshade <- swrFit %>%
  group_by(segment_index, Overshoot) %>%
  summarise(
    xmin = min(calBP, na.rm = TRUE),
    xmax = max(calBP, na.rm = TRUE),
    .groups = "drop"
  )
# Segment boundary positions
seg_bounds <- sort(unique(c(segshade$xmax, segshade$xmax)))

vir_cols <- viridis(10)

swrplat <- ggplot() +
  
  # ==========================================================
# Segment shading
# ==========================================================
geom_rect(
  data = segshade,
  aes(
    xmin = xmin,
    xmax = xmax,
    ymin = -Inf,
    ymax = Inf,
    fill = Overshoot
  ),
  alpha = 0.18,
  inherit.aes = FALSE
) +
  
  # ----------------------------------------------------------
# Segment boundary lines
# ----------------------------------------------------------
geom_vline(
  xintercept = seg_bounds,
  linetype = "dashed",
  color = "grey40",
  linewidth = 1,
  alpha = 0.8
) +
  
  # ==========================================================
# KDE envelope
# ==========================================================
geom_ribbon(
  data = swr1,
  aes(
    x = calBP,
    ymin = lo * 100,
    ymax = hi * 100
  ),
  fill = vir_cols[1],
  alpha = 0.5
) +
  
  # ==========================================================
# Null envelope
# ==========================================================
geom_ribbon(
  data = swrNull, aes(x = calBP, ymin = Y_min,
                      ymax = Y_max ), fill = vir_cols[9],  alpha = 0.35) +
  
  # ==========================================================
# KDE curve
# ==========================================================
geom_path(
  data = swr1,
  aes(
    x = calBP,
    y = MKDE * 100,
    color = "KDE"
  ),
  linewidth = 1.5
) +
  
  # ==========================================================
# Fitted curve
# ==========================================================
geom_path(
  data = swrFit,
  aes(
    x = calBP,
    y = YPredicted,
    color = "Fitted"
  ),
  linewidth = 1.75
) +
  
  # ==========================================================
# Null mean
# ==========================================================
geom_path(
  data = swrNull,
  aes(
    x = calBP,
    y = Y_mean,
    color = "Null mean"
  ),
  linewidth = 1
) +
  
  # ==========================================================
# Axes
# ==========================================================
scale_x_reverse(
  limits = c(4000, 200)
) +
  
  # ==========================================================
# Line colors
# ==========================================================
scale_color_manual(
  name = "Curves",
  values = c(
    "KDE" = vir_cols[1],
    "Fitted" = vir_cols[5],
    "Null mean" = vir_cols[9]
  )
) +
  
  # ==========================================================
# Segment colors
# ==========================================================
scale_fill_viridis(
  name = "Overshoot",
  option = "D",
  direction=1
) +
  
  theme_bw() +
  my_theme +
  
  labs(
    x = "Years cal BP",
    y = "Mean KDE density",
    title = "B. US Southwest Output (Y) vs. Time"
  )

swrplat


###Plot for US Continent===================================
##Continguous US
#===============================================================
na_rec <- read.csv("Inno_Pop_Fitting/Plotting/NARecessions.csv")

fits <- read.csv("Inno_Pop_Fitting/Plotting/NullModel_Envelopes_AllObjectives.csv")

box <- read.csv("Inno_Pop_Fitting/Plotting/Stablity_Robust_Region2.csv")

kde <- read.csv("Inno_Pop_Fitting/Plotting/AllKDE_Hybrid_Corr_Features2.csv")


# ============================================================
# SUBSET DATA
# ============================================================

# Stability table
box <- subset(box, fit_type== "hybrid")

# Create classification variable
box$Classification <- apply(
  box[, c("eig1_real","eig1_imag",
          "eig2_real","eig2_imag",
          "eig3_real","eig3_imag")],
  1,
  function(x) {
    
    real_vals <- c(x[1], x[3], x[5])
    imag_vals <- c(x[2], x[4], x[6])
    
    # Any complex unstable eigenvalues
    if (any(real_vals > 0 & imag_vals != 0)) {
      return("complex unstable")
      
      # Any imaginary component present
    } else if (any(imag_vals != 0)) {
      return("complex stable")
      
      # Otherwise purely real
    } else {
      return("real")
    }
  }
)

# Calculate overshoot values for each eigen pair
ov1 <- log10((abs(box$eig1_imag) / 
                (abs(box$eig1_real) + 0.001)) + 0.001)

ov2 <- log10((abs(box$eig2_imag) / 
                (abs(box$eig2_real) + 0.001)) + 0.001)

ov3 <- log10((abs(box$eig3_imag) / 
                (abs(box$eig3_real) + 0.001)) + 0.001)

# Select highest overshoot across eig1–eig3
box$Overshoot <- pmax(ov1, ov2, ov3)


# Fits table
fits <-subset(fits, objective=="shape")

usc1<-subset(na_rec, region_id=="Continental US")
uscNull<-subset(fits, region_id=="Continental US")
uscFit<-subset(kde, region_id=="Continental US" )
uscEq<-subset(box, region_id=="Continental US" )

###Add calBP to the data sets for plotting
calBP<-uscNull$time*-1
uscNull$calBP <- uscNull$time * -1
uscFit<-cbind(uscFit, calBP)

# Create lookup vectors
class_lookup <- setNames(
  uscEq$Classification,
  uscEq$segment_index
)

overshoot_lookup <- setNames(
  uscEq$Overshoot,
  uscEq$segment_index
)

# Transfer values into uplandswFit
uscFit$Classification <- class_lookup[
  as.character(uscFit$segment_index)
]

uscFit$Overshoot <- overshoot_lookup[
  as.character(uscFit$segment_index)
]

###Plot------------------------------------------
library(ggplot2)
library(dplyr)
library(viridis)

segshade <- uscFit %>%
  group_by(segment_index, Overshoot) %>%
  summarise(
    xmin = min(calBP, na.rm = TRUE),
    xmax = max(calBP, na.rm = TRUE),
    .groups = "drop"
  )
# Segment boundary positions
seg_bounds <- sort(unique(c(segshade$xmax, segshade$xmax)))

vir_cols <- viridis(10)

uscplat <- ggplot() +
  
  # ==========================================================
# Segment shading
# ==========================================================
geom_rect(
  data = segshade,
  aes(
    xmin = xmin,
    xmax = xmax,
    ymin = -Inf,
    ymax = Inf,
    fill = Overshoot
  ),
  alpha = 0.18,
  inherit.aes = FALSE
) +
  
  # ----------------------------------------------------------
# Segment boundary lines
# ----------------------------------------------------------
geom_vline(
  xintercept = seg_bounds,
  linetype = "dashed",
  color = "grey40",
  linewidth = 1,
  alpha = 0.8
) +
  
  # ==========================================================
# KDE envelope
# ==========================================================
geom_ribbon(
  data = usc1,
  aes(
    x = calBP,
    ymin = lo * 100,
    ymax = hi * 100
  ),
  fill = vir_cols[1],
  alpha = 0.5
) +
  
  # ==========================================================
# Null envelope
# ==========================================================
geom_ribbon(
  data = uscNull, aes(x = calBP, ymin = Y_min,
                      ymax = Y_max ), fill = vir_cols[9],  alpha = 0.35) +
  
  # ==========================================================
# KDE curve
# ==========================================================
geom_path(
  data = usc1,
  aes(
    x = calBP,
    y = MKDE * 100,
    color = "KDE"
  ),
  linewidth = 1.5
) +
  
  # ==========================================================
# Fitted curve
# ==========================================================
geom_path(
  data = uscFit,
  aes(
    x = calBP,
    y = YPredicted,
    color = "Fitted"
  ),
  linewidth = 1.75
) +
  
  # ==========================================================
# Null mean
# ==========================================================
geom_path(
  data = uscNull,
  aes(
    x = calBP,
    y = Y_mean,
    color = "Null mean"
  ),
  linewidth = 1
) +
  
  # ==========================================================
# Axes
# ==========================================================
scale_x_reverse(
  limits = c(4000, 200)
) +
  
  # ==========================================================
# Line colors
# ==========================================================
scale_color_manual(
  name = "Curves",
  values = c(
    "KDE" = vir_cols[1],
    "Fitted" = vir_cols[5],
    "Null mean" = vir_cols[9]
  )
) +
  
  # ==========================================================
# Segment colors
# ==========================================================
scale_fill_viridis(
  name = "Overshoot",
  option = "D",
  direction=1
) +
  
  theme_bw() +
  my_theme +
  
  labs(
    x = "Years cal BP",
    y = "Mean KDE density",
    title = "A. US Continent Output (Y) vs. Time"
  )

uscplat

#### Plot figure for overhoot severity and fit
FigSWOvershoot<-plot_grid(uscplat,swrplat, chidesert, sondesert, 
                      Cpt, ncolplat, ncol=1, align="hv", axis = "rl")
FigSWOvershoot

pdf("figures/OvershootRev.pdf", width=16, height=18.55)
FigSWOvershoot
dev.off()


###Stability Classification Plotting ===========================

##================================================================
na_rec <- read.csv("Inno_Pop_Fitting/Plotting/NARecessions.csv")

fits <- read.csv("Inno_Pop_Fitting/Plotting/NullModel_Envelopes_AllObjectives.csv")

box <- read.csv("Inno_Pop_Fitting/Plotting/Stablity_Robust_Region2.csv")

kde <- read.csv("Inno_Pop_Fitting/Plotting/AllKDE_Hybrid_Corr_Features2.csv")


# ============================================================
# SUBSET DATA
# ============================================================

# Stability table
box <- subset(box, fit_type== "hybrid")

# Create classification variable
box$Classification <- apply(
  box[, c("eig1_real","eig1_imag",
          "eig2_real","eig2_imag",
          "eig3_real","eig3_imag")],
  1,
  function(x) {
    
    real_vals <- c(x[1], x[3], x[5])
    imag_vals <- c(x[2], x[4], x[6])
    
    # Any complex unstable eigenvalues
    if (any(real_vals > 0 & imag_vals != 0)) {
      return("complex unstable")
      
      # Any imaginary component present
    } else if (any(imag_vals != 0)) {
      return("complex stable")
      
      # Otherwise purely real
    } else {
      return("real")
    }
  }
)

# Calculate overshoot values for each eigen pair
ov1 <- log10((abs(box$eig1_imag) / 
                (abs(box$eig1_real) + 0.001)) + 0.001)

ov2 <- log10((abs(box$eig2_imag) / 
                (abs(box$eig2_real) + 0.001)) + 0.001)

ov3 <- log10((abs(box$eig3_imag) / 
                (abs(box$eig3_real) + 0.001)) + 0.001)

# Select highest overshoot across eig1–eig3
box$Overshoot <- pmax(ov1, ov2, ov3)



# Fits table
fits <-subset(fits, objective=="shape")

uplandsw1<-subset(na_rec, region_id=="S. Colorado Plat.")
uplandswNull<-subset(fits, region_id=="S. Colorado Plat.")
uplandswFit<-subset(kde, region_id=="S. Colorado Plat." )
uplandswEq<-subset(box, region_id=="S. Colorado Plat." )

###Add calBP to the data sets for plotting
calBP<-uplandswNull$time*-1
uplandswNull$calBP <- uplandswNull$time * -1
uplandswFit<-cbind(uplandswFit, calBP)

# Create lookup vectors
class_lookup <- setNames(
  uplandswEq$Classification,
  uplandswEq$segment_index
)

overshoot_lookup <- setNames(
  uplandswEq$Overshoot,
  uplandswEq$segment_index
)

# Transfer values into uplandswFit
uplandswFit$Classification <- class_lookup[
  as.character(uplandswFit$segment_index)
]

uplandswFit$Overshoot <- overshoot_lookup[
  as.character(uplandswFit$segment_index)
]

###Plot------------------------------------------
library(ggplot2)
library(dplyr)
library(viridis)

segshade <- uplandswFit %>%
  group_by(segment_index, Classification) %>%
  summarise(
    xmin = min(calBP, na.rm = TRUE),
    xmax = max(calBP, na.rm = TRUE),
    .groups = "drop"
  )
# Segment boundary positions
seg_bounds <- sort(unique(c(segshade$xmax, segshade$xmax)))

vir_cols <- viridis(10)

Cpt <- ggplot() +
  
  # ==========================================================
# Segment shading
# ==========================================================
geom_rect(
  data = segshade,
  aes(
    xmin = xmin,
    xmax = xmax,
    ymin = -Inf,
    ymax = Inf,
    fill = Classification
  ),
  alpha = 0.18,
  inherit.aes = FALSE
) +
  
  # ----------------------------------------------------------
# Segment boundary lines
# ----------------------------------------------------------
geom_vline(
  xintercept = seg_bounds,
  linetype = "dashed",
  color = "grey40",
  linewidth = 1,
  alpha = 0.8
) +
  
  # ==========================================================
# KDE envelope
# ==========================================================
geom_ribbon(
  data = uplandsw1,
  aes(
    x = calBP,
    ymin = lo * 100,
    ymax = hi * 100
  ),
  fill = vir_cols[1],
  alpha = 0.5
) +
  
  # ==========================================================
# Null envelope
# ==========================================================
geom_ribbon(
  data = uplandswNull,
  aes(
    x = calBP,
    ymin = Y_min,
    ymax = Y_max
  ),
  fill = vir_cols[9],
  alpha = 0.35
) +
  
  # ==========================================================
# KDE curve
# ==========================================================
geom_path(
  data = uplandsw1,
  aes(
    x = calBP,
    y = MKDE * 100,
    color = "KDE"
  ),
  linewidth = 1.5
) +
  
  # ==========================================================
# Fitted curve
# ==========================================================
geom_path(
  data = uplandswFit,
  aes(
    x = calBP,
    y = YPredicted,
    color = "Fitted"
  ),
  linewidth = 1.75
) +
  
  # ==========================================================
# Null mean
# ==========================================================
geom_path(
  data = uplandswNull,
  aes(
    x = calBP,
    y = Y_mean,
    color = "Null mean"
  ),
  linewidth = 1
) +
  
  # ==========================================================
# Axes
# ==========================================================
scale_x_reverse(
  limits = c(4000, 200)
) +
  
  # ==========================================================
# Line colors
# ==========================================================
scale_color_manual(
  name = "Curves",
  values = c(
    "KDE" = vir_cols[1],
    "Fitted" = vir_cols[5],
    "Null mean" = vir_cols[9]
  )
) +
  
  # ==========================================================
# Segment colors
# ==========================================================
scale_fill_viridis_d(
  name = "Classification",
  option = "D",
  direction=-1
) +
  
  theme_bw() +
  my_theme +
  
  labs(
    x = "Years cal BP",
    y = "Mean KDE density",
    title = "E. S. Colorado Plat. Output (Y) vs. Time"
  )

Cpt


###Plot for Chihuahua Desert===================================
##Chi Desert
#===============================================================
na_rec <- read.csv("Inno_Pop_Fitting/Plotting/NARecessions.csv")

fits <- read.csv("Inno_Pop_Fitting/Plotting/NullModel_Envelopes_AllObjectives.csv")

box <- read.csv("Inno_Pop_Fitting/Plotting/Stablity_Robust_Region2.csv")

kde <- read.csv("Inno_Pop_Fitting/Plotting/AllKDE_Hybrid_Corr_Features2.csv")


# ============================================================
# SUBSET DATA
# ============================================================

# Stability table
box <- subset(box, fit_type== "hybrid")

# Create classification variable
box$Classification <- apply(
  box[, c("eig1_real","eig1_imag",
          "eig2_real","eig2_imag",
          "eig3_real","eig3_imag")],
  1,
  function(x) {
    
    real_vals <- c(x[1], x[3], x[5])
    imag_vals <- c(x[2], x[4], x[6])
    
    # Any complex unstable eigenvalues
    if (any(real_vals > 0 & imag_vals != 0)) {
      return("complex unstable")
      
      # Any imaginary component present
    } else if (any(imag_vals != 0)) {
      return("complex stable")
      
      # Otherwise purely real
    } else {
      return("real")
    }
  }
)

# Calculate overshoot values for each eigen pair
ov1 <- log10((abs(box$eig1_imag) / 
                (abs(box$eig1_real) + 0.001)) + 0.001)

ov2 <- log10((abs(box$eig2_imag) / 
                (abs(box$eig2_real) + 0.001)) + 0.001)

ov3 <- log10((abs(box$eig3_imag) / 
                (abs(box$eig3_real) + 0.001)) + 0.001)

# Select highest overshoot across eig1–eig3
box$Overshoot <- pmax(ov1, ov2, ov3)


# Fits table
fits <-subset(fits, objective=="shape")

chides1<-subset(na_rec, region_id=="Chihuahua Desert")
chidesNull<-subset(fits, region_id=="Chihuahua Desert")
chidesFit<-subset(kde, region_id=="Chihuahua Desert" )
chidesEq<-subset(box, region_id=="Chihuahua Desert" )

###Add calBP to the data sets for plotting
calBP<-chidesNull$time*-1
chidesNull$calBP <- chidesNull$time * -1
chidesFit<-cbind(chidesFit, calBP)

# Create lookup vectors
class_lookup <- setNames(
  chidesEq$Classification,
  chidesEq$segment_index
)

overshoot_lookup <- setNames(
  chidesEq$Overshoot,
  chidesEq$segment_index
)

# Transfer values into uplandswFit
chidesFit$Classification <- class_lookup[
  as.character(chidesFit$segment_index)
]

chidesFit$Overshoot <- overshoot_lookup[
  as.character(chidesFit$segment_index)
]

###Plot------------------------------------------
library(ggplot2)
library(dplyr)
library(viridis)

segshade <- chidesFit %>%
  group_by(segment_index, Classification) %>%
  summarise(
    xmin = min(calBP, na.rm = TRUE),
    xmax = max(calBP, na.rm = TRUE),
    .groups = "drop"
  )
# Segment boundary positions
seg_bounds <- sort(unique(c(segshade$xmax, segshade$xmax)))

vir_cols <- viridis(10)

chidesert <- ggplot() +
  
  # ==========================================================
# Segment shading
# ==========================================================
geom_rect(
  data = segshade,
  aes(
    xmin = xmin,
    xmax = xmax,
    ymin = -Inf,
    ymax = Inf,
    fill = Classification
  ),
  alpha = 0.18,
  inherit.aes = FALSE
) +
  
  # ----------------------------------------------------------
# Segment boundary lines
# ----------------------------------------------------------
geom_vline(
  xintercept = seg_bounds,
  linetype = "dashed",
  color = "grey40",
  linewidth = 1,
  alpha = 0.8
) +
  
  # ==========================================================
# KDE envelope
# ==========================================================
geom_ribbon(
  data = chides1,
  aes(
    x = calBP,
    ymin = lo * 100,
    ymax = hi * 100
  ),
  fill = vir_cols[1],
  alpha = 0.5
) +
  
  # ==========================================================
# Null envelope
# ==========================================================
geom_ribbon(
  data = chidesNull, aes(x = calBP, ymin = Y_min,
                         ymax = Y_max ), fill = vir_cols[9],  alpha = 0.35) +
  
  # ==========================================================
# KDE curve
# ==========================================================
geom_path(
  data = chides1,
  aes(
    x = calBP,
    y = MKDE * 100,
    color = "KDE"
  ),
  linewidth = 1.5
) +
  
  # ==========================================================
# Fitted curve
# ==========================================================
geom_path(
  data = chidesFit,
  aes(
    x = calBP,
    y = YPredicted,
    color = "Fitted"
  ),
  linewidth = 1.75
) +
  
  # ==========================================================
# Null mean
# ==========================================================
geom_path(
  data = chidesNull,
  aes(
    x = calBP,
    y = Y_mean,
    color = "Null mean"
  ),
  linewidth = 1
) +
  
  # ==========================================================
# Axes
# ==========================================================
scale_x_reverse(
  limits = c(4000, 200)
) +
  
  # ==========================================================
# Line colors
# ==========================================================
scale_color_manual(
  name = "Curves",
  values = c(
    "KDE" = vir_cols[1],
    "Fitted" = vir_cols[5],
    "Null mean" = vir_cols[9]
  )
) +
  
  # ==========================================================
# Segment colors
# ==========================================================
scale_fill_viridis_d(
  name = "Classification",
  option = "D",
  direction=-1
) +
  
  theme_bw() +
  my_theme +
  
  labs(
    x = "Years cal BP",
    y = "Mean KDE density",
    title = "C. Chihuahua Desert Output (Y) vs. Time"
  )

chidesert


###Plot for Sonoran Desert===================================
##Sonoran Desert
#===============================================================
na_rec <- read.csv("Inno_Pop_Fitting/Plotting/NARecessions.csv")

fits <- read.csv("Inno_Pop_Fitting/Plotting/NullModel_Envelopes_AllObjectives.csv")

box <- read.csv("Inno_Pop_Fitting/Plotting/Stablity_Robust_Region2.csv")

kde <- read.csv("Inno_Pop_Fitting/Plotting/AllKDE_Hybrid_Corr_Features2.csv")


# ============================================================
# SUBSET DATA
# ============================================================

# Stability table
box <- subset(box, fit_type== "hybrid")

# Create classification variable
box$Classification <- apply(
  box[, c("eig1_real","eig1_imag",
          "eig2_real","eig2_imag",
          "eig3_real","eig3_imag")],
  1,
  function(x) {
    
    real_vals <- c(x[1], x[3], x[5])
    imag_vals <- c(x[2], x[4], x[6])
    
    # Any complex unstable eigenvalues
    if (any(real_vals > 0 & imag_vals != 0)) {
      return("complex unstable")
      
      # Any imaginary component present
    } else if (any(imag_vals != 0)) {
      return("complex stable")
      
      # Otherwise purely real
    } else {
      return("real")
    }
  }
)

# Calculate overshoot values for each eigen pair
ov1 <- log10((abs(box$eig1_imag) / 
                (abs(box$eig1_real) + 0.001)) + 0.001)

ov2 <- log10((abs(box$eig2_imag) / 
                (abs(box$eig2_real) + 0.001)) + 0.001)

ov3 <- log10((abs(box$eig3_imag) / 
                (abs(box$eig3_real) + 0.001)) + 0.001)

# Select highest overshoot across eig1–eig3
box$Overshoot <- pmax(ov1, ov2, ov3)


# Fits table
fits <-subset(fits, objective=="shape")

sondes1<-subset(na_rec, region_id=="Sonoran Desert")
sondesNull<-subset(fits, region_id=="Sonoran Desert")
sondesFit<-subset(kde, region_id=="Sonoran Desert" )
sondesEq<-subset(box, region_id=="Sonoran Desert" )

###Add calBP to the data sets for plotting
calBP<-sondesNull$time*-1
sondesNull$calBP <- sondesNull$time * -1
sondesFit<-cbind(sondesFit, calBP)

# Create lookup vectors
class_lookup <- setNames(
  sondesEq$Classification,
  sondesEq$segment_index
)

overshoot_lookup <- setNames(
  sondesEq$Overshoot,
  sondesEq$segment_index
)

# Transfer values into uplandswFit
sondesFit$Classification <- class_lookup[
  as.character(sondesFit$segment_index)
]

sondesFit$Overshoot <- overshoot_lookup[
  as.character(sondesFit$segment_index)
]

###Plot------------------------------------------
library(ggplot2)
library(dplyr)
library(viridis)

segshade <- sondesFit %>%
  group_by(segment_index, Classification) %>%
  summarise(
    xmin = min(calBP, na.rm = TRUE),
    xmax = max(calBP, na.rm = TRUE),
    .groups = "drop"
  )
# Segment boundary positions
seg_bounds <- sort(unique(c(segshade$xmax, segshade$xmax)))

vir_cols <- viridis(10)

sondesert <- ggplot() +
  
  # ==========================================================
# Segment shading
# ==========================================================
geom_rect(
  data = segshade,
  aes(
    xmin = xmin,
    xmax = xmax,
    ymin = -Inf,
    ymax = Inf,
    fill = Classification
  ),
  alpha = 0.18,
  inherit.aes = FALSE
) +
  
  # ----------------------------------------------------------
# Segment boundary lines
# ----------------------------------------------------------
geom_vline(
  xintercept = seg_bounds,
  linetype = "dashed",
  color = "grey40",
  linewidth = 1,
  alpha = 0.8
) +
  
  # ==========================================================
# KDE envelope
# ==========================================================
geom_ribbon(
  data = sondes1,
  aes(
    x = calBP,
    ymin = lo * 100,
    ymax = hi * 100
  ),
  fill = vir_cols[1],
  alpha = 0.5
) +
  
  # ==========================================================
# Null envelope
# ==========================================================
geom_ribbon(
  data = sondesNull, aes(x = calBP, ymin = Y_min,
                         ymax = Y_max ), fill = vir_cols[9],  alpha = 0.35) +
  
  # ==========================================================
# KDE curve
# ==========================================================
geom_path(
  data = sondes1,
  aes(
    x = calBP,
    y = MKDE * 100,
    color = "KDE"
  ),
  linewidth = 1.5
) +
  
  # ==========================================================
# Fitted curve
# ==========================================================
geom_path(
  data = sondesFit,
  aes(
    x = calBP,
    y = YPredicted,
    color = "Fitted"
  ),
  linewidth = 1.75
) +
  
  # ==========================================================
# Null mean
# ==========================================================
geom_path(
  data = sondesNull,
  aes(
    x = calBP,
    y = Y_mean,
    color = "Null mean"
  ),
  linewidth = 1
) +
  
  # ==========================================================
# Axes
# ==========================================================
scale_x_reverse(
  limits = c(4000, 200)
) +
  
  # ==========================================================
# Line colors
# ==========================================================
scale_color_manual(
  name = "Curves",
  values = c(
    "KDE" = vir_cols[1],
    "Fitted" = vir_cols[5],
    "Null mean" = vir_cols[9]
  )
) +
  
  # ==========================================================
# Segment colors
# ==========================================================
scale_fill_viridis_d(
  name = "Classification",
  option = "D",
  direction=-1
) +
  
  theme_bw() +
  my_theme +
  
  labs(
    x = "Years cal BP",
    y = "Mean KDE density",
    title = "D. Sonoran Desert Output (Y) vs. Time"
  )

sondesert

###Plot for N. Colorado Plat.===================================
##N. Colorado Plat.
#===============================================================
na_rec <- read.csv("Inno_Pop_Fitting/Plotting/NARecessions.csv")

fits <- read.csv("Inno_Pop_Fitting/Plotting/NullModel_Envelopes_AllObjectives.csv")

box <- read.csv("Inno_Pop_Fitting/Plotting/Stablity_Robust_Region2.csv")

kde <- read.csv("Inno_Pop_Fitting/Plotting/AllKDE_Hybrid_Corr_Features2.csv")


# ============================================================
# SUBSET DATA
# ============================================================

# Stability table
box <- subset(box, fit_type== "hybrid")

# Create classification variable
box$Classification <- apply(
  box[, c("eig1_real","eig1_imag",
          "eig2_real","eig2_imag",
          "eig3_real","eig3_imag")],
  1,
  function(x) {
    
    real_vals <- c(x[1], x[3], x[5])
    imag_vals <- c(x[2], x[4], x[6])
    
    # Any complex unstable eigenvalues
    if (any(real_vals > 0 & imag_vals != 0)) {
      return("complex unstable")
      
      # Any imaginary component present
    } else if (any(imag_vals != 0)) {
      return("complex stable")
      
      # Otherwise purely real
    } else {
      return("real")
    }
  }
)

# Calculate overshoot values for each eigen pair
ov1 <- log10((abs(box$eig1_imag) / 
                (abs(box$eig1_real) + 0.001)) + 0.001)

ov2 <- log10((abs(box$eig2_imag) / 
                (abs(box$eig2_real) + 0.001)) + 0.001)

ov3 <- log10((abs(box$eig3_imag) / 
                (abs(box$eig3_real) + 0.001)) + 0.001)

# Select highest overshoot across eig1–eig3
box$Overshoot <- pmax(ov1, ov2, ov3)


# Fits table
fits <-subset(fits, objective=="shape")

ncol1<-subset(na_rec, region_id=="N. Colorado Plat.")
ncolNull<-subset(fits, region_id=="N. Colorado Plat.")
ncolFit<-subset(kde, region_id=="N. Colorado Plat." )
ncolEq<-subset(box, region_id=="N. Colorado Plat." )

###Add calBP to the data sets for plotting
calBP<-ncolNull$time*-1
ncolNull$calBP <- ncolNull$time * -1
ncolFit<-cbind(ncolFit, calBP)

# Create lookup vectors
class_lookup <- setNames(
  ncolEq$Classification,
  ncolEq$segment_index
)

overshoot_lookup <- setNames(
  ncolEq$Overshoot,
  ncolEq$segment_index
)

# Transfer values into uplandswFit
ncolFit$Classification <- class_lookup[
  as.character(ncolFit$segment_index)
]

ncolFit$Overshoot <- overshoot_lookup[
  as.character(ncolFit$segment_index)
]

###Plot------------------------------------------
library(ggplot2)
library(dplyr)
library(viridis)

segshade <- ncolFit %>%
  group_by(segment_index, Classification) %>%
  summarise(
    xmin = min(calBP, na.rm = TRUE),
    xmax = max(calBP, na.rm = TRUE),
    .groups = "drop"
  )
# Segment boundary positions
seg_bounds <- sort(unique(c(segshade$xmax, segshade$xmax)))

vir_cols <- viridis(10)

ncolplat <- ggplot() +
  
  # ==========================================================
# Segment shading
# ==========================================================
geom_rect(
  data = segshade,
  aes(
    xmin = xmin,
    xmax = xmax,
    ymin = -Inf,
    ymax = Inf,
    fill = Classification
  ),
  alpha = 0.18,
  inherit.aes = FALSE
) +
  
  # ----------------------------------------------------------
# Segment boundary lines
# ----------------------------------------------------------
geom_vline(
  xintercept = seg_bounds,
  linetype = "dashed",
  color = "grey40",
  linewidth = 1,
  alpha = 0.8
) +
  
  # ==========================================================
# KDE envelope
# ==========================================================
geom_ribbon(
  data = ncol1,
  aes(
    x = calBP,
    ymin = lo * 100,
    ymax = hi * 100
  ),
  fill = vir_cols[1],
  alpha = 0.5
) +
  
  # ==========================================================
# Null envelope
# ==========================================================
geom_ribbon(
  data = ncolNull, aes(x = calBP, ymin = Y_min,
                       ymax = Y_max ), fill = vir_cols[9],  alpha = 0.35) +
  
  # ==========================================================
# KDE curve
# ==========================================================
geom_path(
  data = ncol1,
  aes(
    x = calBP,
    y = MKDE * 100,
    color = "KDE"
  ),
  linewidth = 1.5
) +
  
  # ==========================================================
# Fitted curve
# ==========================================================
geom_path(
  data = ncolFit,
  aes(
    x = calBP,
    y = YPredicted,
    color = "Fitted"
  ),
  linewidth = 1.75
) +
  
  # ==========================================================
# Null mean
# ==========================================================
geom_path(
  data = ncolNull,
  aes(
    x = calBP,
    y = Y_mean,
    color = "Null mean"
  ),
  linewidth = 1
) +
  
  # ==========================================================
# Axes
# ==========================================================
scale_x_reverse(
  limits = c(2700, 200)
) +
  
  # ==========================================================
# Line colors
# ==========================================================
scale_color_manual(
  name = "Curves",
  values = c(
    "KDE" = vir_cols[1],
    "Fitted" = vir_cols[5],
    "Null mean" = vir_cols[9]
  )
) +
  
  # ==========================================================
# Segment colors
# ==========================================================
scale_fill_viridis_d(
  name = "Classification",
  option = "D",
  direction=-1
) +
  
  theme_bw() +
  my_theme +
  
  labs(
    x = "Years cal BP",
    y = "Mean KDE density",
    title = "F. N. Colorado Plat. Output (Y) vs. Time"
  )

ncolplat


###Plot for US Southwest.===================================
##US Southwest Region
#===============================================================
na_rec <- read.csv("Inno_Pop_Fitting/Plotting/NARecessions.csv")

fits <- read.csv("Inno_Pop_Fitting/Plotting/NullModel_Envelopes_AllObjectives.csv")

box <- read.csv("Inno_Pop_Fitting/Plotting/Stablity_Robust_Region2.csv")

kde <- read.csv("Inno_Pop_Fitting/Plotting/AllKDE_Hybrid_Corr_Features2.csv")


# ============================================================
# SUBSET DATA
# ============================================================

# Stability table
box <- subset(box, fit_type== "hybrid")

# Create classification variable
box$Classification <- apply(
  box[, c("eig1_real","eig1_imag",
          "eig2_real","eig2_imag",
          "eig3_real","eig3_imag")],
  1,
  function(x) {
    
    real_vals <- c(x[1], x[3], x[5])
    imag_vals <- c(x[2], x[4], x[6])
    
    # Any complex unstable eigenvalues
    if (any(real_vals > 0 & imag_vals != 0)) {
      return("complex unstable")
      
      # Any imaginary component present
    } else if (any(imag_vals != 0)) {
      return("complex stable")
      
      # Otherwise purely real
    } else {
      return("real")
    }
  }
)

# Calculate overshoot values for each eigen pair
ov1 <- log10((abs(box$eig1_imag) / 
                (abs(box$eig1_real) + 0.001)) + 0.001)

ov2 <- log10((abs(box$eig2_imag) / 
                (abs(box$eig2_real) + 0.001)) + 0.001)

ov3 <- log10((abs(box$eig3_imag) / 
                (abs(box$eig3_real) + 0.001)) + 0.001)

# Select highest overshoot across eig1–eig3
box$Overshoot <- pmax(ov1, ov2, ov3)


# Fits table
fits <-subset(fits, objective=="shape")

swr1<-subset(na_rec, region_id=="Southwest US")
swrNull<-subset(fits, region_id=="Southwest US")
swrFit<-subset(kde, region_id=="Southwest US" )
swrEq<-subset(box, region_id=="Southwest US" )

###Add calBP to the data sets for plotting
calBP<-swrNull$time*-1
swrNull$calBP <- swrNull$time * -1
swrFit<-cbind(swrFit, calBP)

# Create lookup vectors
class_lookup <- setNames(
  swrEq$Classification,
  swrEq$segment_index
)

overshoot_lookup <- setNames(
  swrEq$Overshoot,
  swrEq$segment_index
)

# Transfer values into uplandswFit
swrFit$Classification <- class_lookup[
  as.character(swrFit$segment_index)
]

swrFit$Overshoot <- overshoot_lookup[
  as.character(swrFit$segment_index)
]

###Plot------------------------------------------
library(ggplot2)
library(dplyr)
library(viridis)

segshade <- swrFit %>%
  group_by(segment_index, Classification) %>%
  summarise(
    xmin = min(calBP, na.rm = TRUE),
    xmax = max(calBP, na.rm = TRUE),
    .groups = "drop"
  )
# Segment boundary positions
seg_bounds <- sort(unique(c(segshade$xmax, segshade$xmax)))

vir_cols <- viridis(10)

swrplat <- ggplot() +
  
  # ==========================================================
# Segment shading
# ==========================================================
geom_rect(
  data = segshade,
  aes(
    xmin = xmin,
    xmax = xmax,
    ymin = -Inf,
    ymax = Inf,
    fill = Classification
  ),
  alpha = 0.18,
  inherit.aes = FALSE
) +
  
  # ----------------------------------------------------------
# Segment boundary lines
# ----------------------------------------------------------
geom_vline(
  xintercept = seg_bounds,
  linetype = "dashed",
  color = "grey40",
  linewidth = 1,
  alpha = 0.8
) +
  
  # ==========================================================
# KDE envelope
# ==========================================================
geom_ribbon(
  data = swr1,
  aes(
    x = calBP,
    ymin = lo * 100,
    ymax = hi * 100
  ),
  fill = vir_cols[1],
  alpha = 0.5
) +
  
  # ==========================================================
# Null envelope
# ==========================================================
geom_ribbon(
  data = swrNull, aes(x = calBP, ymin = Y_min,
                      ymax = Y_max ), fill = vir_cols[9],  alpha = 0.35) +
  
  # ==========================================================
# KDE curve
# ==========================================================
geom_path(
  data = swr1,
  aes(
    x = calBP,
    y = MKDE * 100,
    color = "KDE"
  ),
  linewidth = 1.5
) +
  
  # ==========================================================
# Fitted curve
# ==========================================================
geom_path(
  data = swrFit,
  aes(
    x = calBP,
    y = YPredicted,
    color = "Fitted"
  ),
  linewidth = 1.75
) +
  
  # ==========================================================
# Null mean
# ==========================================================
geom_path(
  data = swrNull,
  aes(
    x = calBP,
    y = Y_mean,
    color = "Null mean"
  ),
  linewidth = 1
) +
  
  # ==========================================================
# Axes
# ==========================================================
scale_x_reverse(
  limits = c(4000, 200)
) +
  
  # ==========================================================
# Line colors
# ==========================================================
scale_color_manual(
  name = "Curves",
  values = c(
    "KDE" = vir_cols[1],
    "Fitted" = vir_cols[5],
    "Null mean" = vir_cols[9]
  )
) +
  
  # ==========================================================
# Segment colors
# ==========================================================
scale_fill_viridis_d(
  name = "Classification",
  option = "D",
  direction=-1
) +
  
  theme_bw() +
  my_theme +
  
  labs(
    x = "Years cal BP",
    y = "Mean KDE density",
    title = "B. US Southwest Output (Y) vs. Time"
  )

swrplat


###Plot for US Continent===================================
##Continguous US
#===============================================================
na_rec <- read.csv("Inno_Pop_Fitting/Plotting/NARecessions.csv")

fits <- read.csv("Inno_Pop_Fitting/Plotting/NullModel_Envelopes_AllObjectives.csv")

box <- read.csv("Inno_Pop_Fitting/Plotting/Stablity_Robust_Region2.csv")

kde <- read.csv("Inno_Pop_Fitting/Plotting/AllKDE_Hybrid_Corr_Features2.csv")


# ============================================================
# SUBSET DATA
# ============================================================

# Stability table
box <- subset(box, fit_type== "hybrid")

# Create classification variable
box$Classification <- apply(
  box[, c("eig1_real","eig1_imag",
          "eig2_real","eig2_imag",
          "eig3_real","eig3_imag")],
  1,
  function(x) {
    
    real_vals <- c(x[1], x[3], x[5])
    imag_vals <- c(x[2], x[4], x[6])
    
    # Any complex unstable eigenvalues
    if (any(real_vals > 0 & imag_vals != 0)) {
      return("complex unstable")
      
      # Any imaginary component present
    } else if (any(imag_vals != 0)) {
      return("complex stable")
      
      # Otherwise purely real
    } else {
      return("real")
    }
  }
)

# Calculate overshoot values for each eigen pair
ov1 <- log10((abs(box$eig1_imag) / 
                (abs(box$eig1_real) + 0.001)) + 0.001)

ov2 <- log10((abs(box$eig2_imag) / 
                (abs(box$eig2_real) + 0.001)) + 0.001)

ov3 <- log10((abs(box$eig3_imag) / 
                (abs(box$eig3_real) + 0.001)) + 0.001)

# Select highest overshoot across eig1–eig3
box$Overshoot <- pmax(ov1, ov2, ov3)


# Fits table
fits <-subset(fits, objective=="shape")

usc1<-subset(na_rec, region_id=="Continental US")
uscNull<-subset(fits, region_id=="Continental US")
uscFit<-subset(kde, region_id=="Continental US" )
uscEq<-subset(box, region_id=="Continental US" )

###Add calBP to the data sets for plotting
calBP<-uscNull$time*-1
uscNull$calBP <- uscNull$time * -1
uscFit<-cbind(uscFit, calBP)

# Create lookup vectors
class_lookup <- setNames(
  uscEq$Classification,
  uscEq$segment_index
)

overshoot_lookup <- setNames(
  uscEq$Overshoot,
  uscEq$segment_index
)

# Transfer values into uplandswFit
uscFit$Classification <- class_lookup[
  as.character(uscFit$segment_index)
]

uscFit$Overshoot <- overshoot_lookup[
  as.character(uscFit$segment_index)
]

###Plot------------------------------------------
library(ggplot2)
library(dplyr)
library(viridis)

segshade <- uscFit %>%
  group_by(segment_index, Classification) %>%
  summarise(
    xmin = min(calBP, na.rm = TRUE),
    xmax = max(calBP, na.rm = TRUE),
    .groups = "drop"
  )
# Segment boundary positions
seg_bounds <- sort(unique(c(segshade$xmax, segshade$xmax)))

vir_cols <- viridis(10)

uscplat <- ggplot() +
  
  # ==========================================================
# Segment shading
# ==========================================================
geom_rect(
  data = segshade,
  aes(
    xmin = xmin,
    xmax = xmax,
    ymin = -Inf,
    ymax = Inf,
    fill = Classification
  ),
  alpha = 0.18,
  inherit.aes = FALSE
) +
  
  # ----------------------------------------------------------
# Segment boundary lines
# ----------------------------------------------------------
geom_vline(
  xintercept = seg_bounds,
  linetype = "dashed",
  color = "grey40",
  linewidth = 1,
  alpha = 0.8
) +
  
  # ==========================================================
# KDE envelope
# ==========================================================
geom_ribbon(
  data = usc1,
  aes(
    x = calBP,
    ymin = lo * 100,
    ymax = hi * 100
  ),
  fill = vir_cols[1],
  alpha = 0.5
) +
  
  # ==========================================================
# Null envelope
# ==========================================================
geom_ribbon(
  data = uscNull, aes(x = calBP, ymin = Y_min,
                      ymax = Y_max ), fill = vir_cols[9],  alpha = 0.35) +
  
  # ==========================================================
# KDE curve
# ==========================================================
geom_path(
  data = usc1,
  aes(
    x = calBP,
    y = MKDE * 100,
    color = "KDE"
  ),
  linewidth = 1.5
) +
  
  # ==========================================================
# Fitted curve
# ==========================================================
geom_path(
  data = uscFit,
  aes(
    x = calBP,
    y = YPredicted,
    color = "Fitted"
  ),
  linewidth = 1.75
) +
  
  # ==========================================================
# Null mean
# ==========================================================
geom_path(
  data = uscNull,
  aes(
    x = calBP,
    y = Y_mean,
    color = "Null mean"
  ),
  linewidth = 1
) +
  
  # ==========================================================
# Axes
# ==========================================================
scale_x_reverse(
  limits = c(4000, 200)
) +
  
  # ==========================================================
# Line colors
# ==========================================================
scale_color_manual(
  name = "Curves",
  values = c(
    "KDE" = vir_cols[1],
    "Fitted" = vir_cols[5],
    "Null mean" = vir_cols[9]
  )
) +
  
  # ==========================================================
# Segment colors
# ==========================================================
scale_fill_viridis_d(
  name = "Classification",
  option = "D",
  direction=-1
) +
  
  theme_bw() +
  my_theme +
  
  labs(
    x = "Years cal BP",
    y = "Mean KDE density",
    title = "A. US Continent Output (Y) vs. Time"
  )

uscplat

#### Plot figure for overhoot severity and fit
FigSWClass<-plot_grid(uscplat,swrplat, chidesert, sondesert, 
                          Cpt, ncolplat, ncol=1, align="hv", axis = "rl")
FigSWClass

pdf("figures/ClassRev.pdf", width=16, height=18.55)
FigSWClass
dev.off()


###Model ``validation" exercises=======================================

Inq<-read.csv("gini_database.csv", header=TRUE)
Inqna<-subset(Inq, Subregion=="Upland SW" & TotalAreaHouse>0 & EndHouse>-2501 & EndHouse<1801)


####Graph fit data for the South Colorado Plat.
fituplandSW<-subset(kde, region_id=="S. Colorado Plat.")

#Transform house dates into cal BP
Inqna$EndHouseBP<-abs(1950-Inqna$EndHouse)


Inqna$PhaseTime <- as.numeric(as.character(cut(
  Inqna$EndHouseBP,
  breaks = c(150, 550, 1450, 2400, 3000, 3500, 4500),
  labels = c(300, 1000, 1950, 3000, 3950, 4500),
  include.lowest = TRUE
)))

Inqna$PeriodID <- cut(Inqna$EndHouseBP,
                      breaks=c(200, 320, 670, 820, 1050, 1210, 1455, 3430, 4000),
                      labels=c('PV','PIV','PIII','PII','PI','BMIII','BMII','Early Agricultural Period'))


#Inqna$PhaseTime

Meanseq <- Inqna%>% group_by(PhaseTime) %>%
  summarize(Avg = median(TotalAreaHouse))
Meanseq

library(ineq)

Meanseq2 <- Inqna %>%
  group_by(PhaseTime) %>%
  summarize(
    gini = ineq(TotalAreaHouse, type = "Gini", na.rm = TRUE)
  )
Meanseq2

Inqna2<-subset(Inqna, EndHouseBP<4001 & EndHouseBP>199)

phase_bands <- data.frame(
  xmin = c(4000, 3430, 1455, 1210, 1050, 820, 670, 320),
  xmax = c(3430, 1455, 1210, 1050, 820, 670, 320, 200),
  phase = c("EAP", "BMII", "BMIII", "PI", "PII", "PIII", "PIV", "PV")
)


library(viridis)

# consistent viridis picks
vir_cols <- viridis(4) # [1]=purple, [2]=blue, [3]=green, [4]=yellow

hs1 <- ggplot(Inqna2, aes(x = EndHouseBP, y = log(TotalAreaHouse))) +
  
  # 🔹 Phase boundaries (vertical lines)
  geom_vline(
    data = phase_bands,
    aes(xintercept = xmin),
    inherit.aes = FALSE,
    linetype = "dashed",
    color = "grey30",
    linewidth = 0.8
  ) +

  geom_vline(
    data = phase_bands,
    aes(xintercept = xmax),
    inherit.aes = FALSE,
    linetype = "dashed",
    color = "grey30",
    linewidth = 0.8
  ) +
  
  # Points (keep neutral)
  geom_point(size = 1, color = "black") +
  
  # Phase labels
  geom_text(
    data = phase_bands,
    aes(x = (xmin + xmax)/2,
        y = 5.7,
        label = phase),
    inherit.aes = FALSE,
    angle = 90,
    vjust = 0.5,
    size = 5
  ) +
  
  # Fitted curve (green, viridis)
  geom_path(
    data = fituplandSW,
    aes(x = calBP, y = log(KPredicted * 60), color = "predicted K"),
    size = 2
  ) +
  
  # Smoothed trend
 # geom_smooth(
  #  aes(color = "Trend"),
   # se = FALSE
  #) +
  
  theme_bw() +
  scale_x_reverse(limits = c(4000, 200)) +
  
  # Viridis mapping (consistent with earlier figures)
  scale_color_manual(
    name = "Curves",
    values = c(
      "Trend" = vir_cols[4],   # purple
      "predicted K" = vir_cols[3]   # green
    )
  ) +
  
  labs(
    x = "Years cal BP",
    y = "ln Structure area and predicted K",
    title = "A. South Colorado Plat. Structure Area (K) vs. Time"
  ) +
  my_theme
hs1



###Graph C Alpha, beta and Phi=================================


Params <- ggplot(fituplandSW, aes(x = calBP)) +
  
  # 🔹 Phase boundaries (vertical lines)
  geom_vline(
    data = phase_bands,
    aes(xintercept = xmin),
    inherit.aes = FALSE,
    linetype = "dashed",
    color = "grey30",
    linewidth = 0.8
  ) +
  
  geom_vline(
    data = phase_bands,
    aes(xintercept = xmax),
    inherit.aes = FALSE,
    linetype = "dashed",
    color = "grey30",
    linewidth = 0.8
  ) +
  
  # Phase labels
  geom_text(
    data = phase_bands,
    aes(x = (xmin + xmax)/2,
        y = 0.8,
        label = phase),
    inherit.aes = FALSE,
    angle = 90,
    vjust = 0.5,
    size = 5
  ) +
  
  # be (blue)
  geom_path(
    aes(y = be, color = "be"),
    size = 2.5
  ) +
  
  # A (purple)
  geom_path(
    aes(y = A, color = "A"),
    size = 1.5
  ) +
  
  # phi (green)
  geom_path(
    aes(y = phi, color = "phi"),
    size = 1
  ) +
  
  theme_bw() +
  scale_x_reverse(limit = c(4000, 200)) +
  scale_y_continuous(limits = c(0, 0.85)) +
  
  #Same viridis mapping logic as previous plots
  scale_color_manual(
    name = "Parameters",
    values = c(
      "A" = vir_cols[1],   # purple
      "be" = vir_cols[4],  # yellow
      "phi" = vir_cols[3]  # green
    )
  ) +
  
  labs(
    x = "Years cal BP",
    y = "Fit parameter values (A, be, phi)",
    title = "B. South Colorado Plat. Parameters vs. Time"
  ) +
  
  my_theme +
  
  # Labels (kept as-is)
  annotate("text", x = 3500, y = .65, label = "be", size = 6) +
  annotate("text", x = 3500, y = .314, label = "phi", size = 6) +
  annotate("text", x = 3500, y = .05, label = "A", size = 6)
Params

FigUplandSW<-plot_grid(hs1, Params, ncol=1, align="hv", axis = "rl")
FigUplandSW

pdf("figures/UplandSW.pdf", width=12.55, height=10)
FigUplandSW
dev.off()


