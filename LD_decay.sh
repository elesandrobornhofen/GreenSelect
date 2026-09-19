#!/bin/bash
#SBATCH --job-name=LD_decay
#SBATCH --output=LD_decay_%j.out
#SBATCH --error=LD_decay_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

# Activate your conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate !!!!!!!!!!!!!YOUR_ENVIRONMENT!!!!!!!!!!!!!

# Run R code directly from this Bash script
R --vanilla <<'EOF'

options(scipen = 6, digits = 4)

## Load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  ggplot2,
  dplyr,
  data.table,
  glue,
  tidyr,
  scales
)

## Run PLINK
system(glue(
  "plink --bfile ../greenSelectParents ",
  "--allow-extra-chr ",
  "--ld-window 999999 ",
  "--ld-window-kb 1000 ",
  "--ld-window-r2 0 ",
  "--r2 ",
  "--maf 0.05 ",
  "--out pairwise_sqr"
))

## Read PLINK output
ld_out <- fread("pairwise_sqr.ld", data.table = FALSE)

## Prepare file
ld_out_sorted <- ld_out |>
  mutate(dist = BP_B - BP_A) |>
  arrange(dist)

## Create bins
bin_size <- 100

ld_out_sorted$distc <- cut(
  ld_out_sorted$dist,
  breaks = seq(
    from = min(ld_out_sorted$dist),
    to = max(ld_out_sorted$dist),
    by = bin_size
  ),
  right = FALSE
)

ld_averages <- ld_out_sorted %>%
  group_by(distc) %>%
  summarise(
    Avg_distance = mean(dist),
    Avg_R2 = mean(R2)
  )

## Compute half-decay using LOESS
fit <- loess(
  Avg_R2 ~ Avg_distance,
  data = ld_averages,
  span = 0.1
)

new_df <- data.frame(
  Avg_distance = seq(
    min(ld_averages$Avg_distance),
    max(ld_averages$Avg_distance),
    length.out = 1000
  )
)

new_df$pred <- predict(
  fit,
  newdata = new_df
)

half_max <- max(new_df$pred, na.rm = TRUE) * 0.5

half_decay_bp <- new_df[
  min(which(new_df$pred <= half_max)),
  "Avg_distance"
]

cat(
  "Half decay distance:",
  half_decay_bp,
  "bp\n"
)

## Plot LD decay
ld_averages |>
  drop_na() |>
  ggplot() +
  aes(
    x = Avg_distance,
    y = Avg_R2
  ) +
  geom_point(
    shape = "circle",
    size = 1.5,
    colour = "grey50",
    alpha = 0.5
  ) +
  stat_smooth(
    method = "loess",
    se = FALSE,
    colour = "darkred",
    alpha = 1,
    size = 1,
    span = 0.1
  ) +
  labs(
    x = "Distance (kb)",
    y = "R2"
  ) +
  geom_hline(
    yintercept = half_max,
    linetype = "dashed",
    color = "black",
    size = 1
  ) +
  xlab(
    paste(
      "Distance in bp.",
      "LD decay =",
      half_decay_bp
    )
  ) +
  scale_y_continuous(
    limits = c(0.0, 0.2),
    n.breaks = 10
  ) +
  scale_x_continuous(
    labels = scales::label_number(
      scale = 1e-3,
      suffix = "kb"
    ),
    n.breaks = 10
  ) +
  theme_classic() -> p1

ggsave(
  plot = p1,
  filename = "ld_decay_plot.pdf"
)

EOF