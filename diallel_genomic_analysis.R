##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##
## Script name: Diallel analysis of perennial ryegrass
##
## ::::::::::::::::::::::::::::
##
## Notes: This script uses an external package (DMU) for some of the analysis
##
## ::::::::::::::::::::::::::::

options(scipen = 6, digits = 4) # Non-scientific notation

## ::::::::::::::::::::::::::::

## load up the packages we will need

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, magrittr, data.table, ppcor, glue, viridis, lme4, doParallel, foreach,
               stringr, RColorBrewer, ggh4x, cowplot) # List all necessary packages here

##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Custom Functions ----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Compute the sparse inverse of a genomic relationship matrix
#' @param mat A square genomic relationship matrix (GRM)
#' @return A sparse inverse formatted for DMU software
sparse_inverse <- function(mat) {
  eigG <- eigen(mat)
  values <- eigG$values
  ii <- which(values<=0)
  if(length(ii) <= 1){
    values[length(values)] <- values[length(values) - 1] * 0.8
  } else {
    for(i in 1:length(ii)){
      values[length(values) - length(ii) + i] <- values[length(values) - length(ii) + i - 1] * 0.8
    }
  }
  #  tail(values, (length(ii) + 5))
  Geigcor_inv <- eigG$vectors %*% diag(1 / values) %*% t(eigG$vectors)
  
  ## Sparse inverse -
  colnames(Geigcor_inv) = rownames(Geigcor_inv) = rownames(mat)
  Ginv_sparse <- Geigcor_inv
  Ginv_sparse[lower.tri(Ginv_sparse)] <- NA
  Ginv_sparse <- reshape2::melt(Ginv_sparse)
  Ginv_sparse[, 1:2] <- lapply(Ginv_sparse[, 1:2], as.character)
  Ginv_sparse <- Ginv_sparse[!is.na(Ginv_sparse[,3]),]
  G.det <- as.numeric(determinant(mat,logarithm = T)[1]) # log determinant
  if(is.na(G.det) == TRUE | G.det == Inf | G.det == -Inf) {
    G.det <- 1
  }
  G.det.line1 <- c(0, 0, G.det)
  Ginv_sparse_inverse <- rbind(G.det.line1, Ginv_sparse)
  
  return(Ginv_sparse_inverse)
}

#' Wrapper to execute DMU for mixed models
#' @param y Dataframe containing phenotype records
#' @param trait Name of the trait being analyzed
#' @param dir Text string containing the DMU directive file parameters
#' @param G1 Sparse inverse of Parental GRM (default NULL)
#' @param G2 Sparse inverse of Crosses GRM (default NULL)
#' @param dmu4_prior Priors for cross-validation evaluation (default NULL)
#' @param save_lst Logical, if TRUE saves the DMU .lst file to results
#' @param path_to_dmu Path to DMU executables
run_dmu <- function(y, trait, dir, G1 = NULL, G2 = NULL, dmu4_prior = NULL, save_lst = FALSE, path_to_dmu = "../softwares/dmu/") {
  
  original_dir <- getwd()
  
  # Create a temporary directory to work in
  temp_dir <- tempfile(pattern = "dmu_", tmpdir = "")
  temp_dir <- gsub("\\\\", "", temp_dir)                      # Remove backslash
  dir.create(temp_dir)
  on.exit({
    setwd(original_dir)                                       # Return to original working directory
    #unlink(temp_dir, recursive = TRUE)                        # Delete temp dir on exit
  }, add = TRUE)
  
  # Copy the DMU to the temporary directory
  dmu_files <- list.files(path_to_dmu, recursive = TRUE)
  for (file in dmu_files) {
    file.copy(paste0(path_to_dmu, "/", file), file.path(temp_dir))
  }
  
  # Run DMU
  setwd(temp_dir)
  
  # Write files to disk
  fwrite(y, "data_frame.dat", sep = "\t", col.names = FALSE)  # Write phenotype file
  writeLines(dir, con = paste0(trait, ".DIR"))                # Write the directive file
  if (!is.null(G1)) {
    fwrite(G1, "GRM1.txt", sep = "\t", col.names = F)         # Write GRM 1
  }
  if (!is.null(G2)) {
    fwrite(G2, "GRM2.txt", sep = "\t", col.names = F)         # Write GRM 2
  }
  
  if (!is.null(dmu4_prior)) {
    
    fwrite(dmu4_prior[[trait]][["par"]], "prior.txt", sep = "\t", col.names = F)
    
    system(glue("run_dmu4 {trait}"), intern = T)               # Run DMU
    
    # Read outputs
    sol <- data.table::fread(paste0(trait, ".SOL"), data.table = F)
    
    # Return outputs
    list(sol = sol)
    
  } else {
    
    system(glue("run_dmuai {trait}"), intern = T)               # Run DMU
    
    status <- readLines(paste0(trait, ".LLIK"))
    delta <- gsub("[^0-9.E-]+", "", status[3]) |> as.numeric()
    gradient <- gsub("[^0-9.E-]+", "", status[4]) |> as.numeric()
    
    if (delta < 1.0e-7 | gradient < 1.0e-6) {
      print("Converged! Life is good =)")
    } else {
      stop("Not converged!")
    }
    
    # Read outputs
    var <- as.numeric(read.table(paste0(trait, ".PAROUT_STD"), nrows = 1))
    var <- read.table(paste0(trait, ".PAROUT_STD"), skip = 1, nrows = var)
    par <- read.table(paste0(trait, ".PAROUT"))
    sol <- data.table::fread(paste0(trait, ".SOL"), data.table = F)
    res <- data.table::fread(paste0(trait, ".RESIDUAL"), data.table = F)
    
    if (save_lst) {
      lst_path <- paste0(trait, ".lst")
      new_path <- paste0("../../results/", trait, ".lst")
      file.rename(from = lst_path, to = new_path)
    }
    
    # Return outputs
    list(convergency = list(delta = delta, gradient = gradient), var = var, par = par, sol = sol, res = res)
    
  }
  
}

#' Construct an additive Genomic Relationship Matrix (Method 1 VanRaden)
#' @param X Marker matrix (individuals in rows, markers in columns)
#' @param use_clump Optional vector of marker names to subset
#' @param maf Minor allele frequency threshold (0 to 0.5)
#' @param n.core Number of cores for parallel processing
myGmat = function(X, use_clump = NULL, maf = NULL, n.core = 1) {
  if (!is.matrix(X)) {
    if (is.data.frame(X)) {
      X <- as.matrix(X)
    } else {
      stop("Input X must be a matrix or data frame.")
    }
  }
  if (!is.null(maf)) {
    if (!is.numeric(maf) || length(maf) != 1 || maf < 0 || maf > 0.5) {
      stop("maf must be a single number between 0 and 0.5.")
    }
  }
  if (!is.null(use_clump)) {
    X <- X[, which(colnames(X) %in% use_clump), drop = FALSE]
  }
  
  n = nrow(X)
  m = ncol(X)
  
  # Allele frequencies
  # Only for Linux machine
  if (n.core > 1) {
    it = split(1:m, factor(cut(1:m, n.core, labels = FALSE)))
    resit = parallel::mclapply(it, function(markers) {
      apply(X[, markers, drop = FALSE], 2, mean) / 2
    }, mc.cores = n.core)
    p = unlist(resit, use.names = FALSE)
  } else {
    p = apply(X, 2, mean) / 2
  }
  
  # Minor allele frequency filter
  if (!is.null(maf)) {
    keep = pmin(p, 1 - p) >= maf
    if (!any(keep)) stop("No SNPs remain after MAF filtering at maf = ", maf, ".")
    message(sum(!keep), " of ", length(keep), " SNPs removed by MAF filter (maf < ", maf, ").")
    X = X[, keep, drop = FALSE]
    p = p[keep]
    m = ncol(X)
  }
  
  v1 = matrix(1, n, 1)
  q = 1 - p
  var.A = 2 * mean(p * q)
  Mp = tcrossprod(v1, matrix(p, m, 1))
  W = X - 2 * Mp
  A = tcrossprod(W) / var.A / m
  return(A)
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Data Loading & Preparation ----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Phenotypes
pheno <- read_csv("../data/2023.02.23_phenotypes.csv", show_col_types = F)

# Genotypes
GD <- fread("../data/greenSelectParents_QUAL20_MQ10_maf005_miss50_imputed_GD.csv", data.table = F)
GD %<>% column_to_rownames("Taxa")
GM <- fread("../data/greenSelectParents_QUAL20_MQ10_maf005_miss50_imputed_GM.csv", data.table = F)

# Prepare data for dmu
ids <- data.frame(parent = unique(c(pheno$gp1, pheno$gp2)))
numerical_dictionary <- ids |> 
  mutate(parent_numeric = as.numeric(as.factor(parent)))

pheno_dmu <- pheno |> 
  filter(!is.na(gp1) & !is.na(gp2)) |> # Remove checks
  left_join(numerical_dictionary, by = c("gp1" = "parent")) |> 
  rename(P1 = parent_numeric) |> 
  left_join(numerical_dictionary, by = c("gp2" = "parent")) |> 
  rename(P2 = parent_numeric)

pheno_dmu <- pheno |> 
  filter(!is.na(gp1) & !is.na(gp2)) |> # Remove checks
  left_join(numerical_dictionary, by = c("gp1" = "parent")) |> 
  rename(P1 = parent_numeric) |> 
  left_join(numerical_dictionary, by = c("gp2" = "parent")) |> 
  rename(P2 = parent_numeric) |> 
  mutate(P1P2 = paste0(P1, P2),
         E = paste0(year, loc), # Environment is equal to the year and location combinations
         E_N = paste0(E, manag),
         E_N_B = paste0(E_N, block),
         P1_E = paste0(P1, E),
         P2_E = paste0(P2, E),
         P1P2_E = paste0(P1P2, E),
         P1_N = paste0(P1, manag),
         P2_N = paste0(P2, manag),
         P1P2_N = paste0(P1P2, manag)) |> 
  dplyr::select(E, E_N, cov, E_N_B, P1, P2, P1P2, P1_E, P2_E, P1P2_E, P1_N, P2_N, P1P2_N, row, col, sumDMY:WSC) |> 
  mutate(across(sumDMY:WSC, ~replace_na(., -999)))

## Make GRMs ----
# Change parental ids to numeric in the marker matrix
GD_tmp <- GD |> 
  rownames_to_column("parent") |> 
  left_join(numerical_dictionary, by = "parent") |> 
  column_to_rownames("parent_numeric") |> 
  dplyr::select(-parent) |> 
  as.matrix()

# Make relationship matrices
G_parents <- myGmat(GD_tmp, maf = 0.05)
G_crosses <- kronecker(G_parents, G_parents, make.dimnames = T)
colnames(G_crosses) = rownames(G_crosses) = sub(":", "", colnames(G_crosses))
cross_ids <- unique(paste0(pheno_dmu$P1, pheno_dmu$P2))
G_crosses = G_crosses[cross_ids, cross_ids]; rm(cross_ids)

# Make sparse inverse
G_inv_parents <- sparse_inverse(G_parents)
G_inv_crosses <- sparse_inverse(G_crosses)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Model Fitting (Full LMM: GCA + SCA) ----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Full model ----
trait_list <- c("sumDMY", "ADF","ADL", "DMDig", "NDF", "NDFD", "Prot", "WSC" )
out_full <- list()
for (i in 1:length(trait_list)) {
  
  trait = trait_list[i]
  
  glue_directive <- glue("
  $ANALYSE 1 1 0 0
  
  $DATA  ASCII (15, 8, -999) data_frame.dat
   
  $VARIABLE
  {glue_collapse(colnames(pheno_dmu),  sep = ' ')}


  $MODEL
  1
  0
  {i} 0 13 0 2 3 4 5 6 7 8 9 10 11 12 13
  10             1 2+2 3 4+4 5   6+6  7 
  0
  0

  $VAR_STR 2 GREL ASCII GRM1.txt
  $VAR_STR 3 GREL ASCII GRM2.txt
  $RESIDUALS ASCII

  $DMUAI
  10
  1.0d-7
  1.0d-6
  1
  0
  0
  ")
  
  pheno_in <- pheno_dmu |> 
    filter(!!sym(trait) != -999)
  
  fm <- run_dmu(y = pheno_in, 
                trait = trait, 
                G1 = G_inv_parents, 
                G2 = G_inv_crosses, 
                save_lst = F, 
                dir = glue_directive)
  out_full[[trait]] <- fm
  
  df <- cbind(pheno_in, resid = out_full[[trait]]$res[,4])
  df |> 
    ggplot(aes(row, col, fill = resid)) +
    geom_tile() +
    ggtitle(trait) +
    scale_fill_viridis(direction = -1, option = "turbo") +
    facet_wrap(~E, ncol = 2) -> p
  ggsave(plot = p, filename = paste0("../figures/", trait, ".png"))
  
}
saveRDS(out_full, "../results/full_model_outputs.RDS")
if (!exists("out_full") || is.null(out_full) || identical(out_full, list())) {
  out_full <- readRDS("../results/full_model_outputs.RDS")
}

### Variance Components & Heritability ----
# Variance components
combined_df <- lapply(names(out_full), function(sublist_name) {
  df <- out_full[[sublist_name]][["var"]]
  df$Sublist <- sublist_name
  df
}) |> 
  bind_rows() |> 
  dplyr::select(trait = Sublist, random_mat = V1, var = V2) |> 
  mutate(trait = ifelse(trait == "sumDMY", "DMY", trait))

var_comp <- combined_df |> 
  mutate(var = ifelse(random_mat == 2, var * 2 * mean(diag(G_parents)), var),
         var = ifelse(random_mat == 3, var * mean(diag(G_crosses)), var),
         var = ifelse(random_mat == 4, var * 2, var),
         var = ifelse(random_mat == 6, var * 2, var)) |> 
  group_by(trait) |> 
  mutate(total_var = sum(var),
         var_percentage = var / total_var * 100) |> 
  ungroup() |> 
  mutate(random_mat = case_when(
    random_mat == 1 ~ "Intra block",
    random_mat == 2 ~ "GCA",
    random_mat == 3 ~ "SCA",
    random_mat == 4 ~ "GCAxE",
    random_mat == 5 ~ "SCAxE",
    random_mat == 6 ~ "GCAxN",
    random_mat == 7 ~ "SCAxN",
    random_mat == 8 ~ "Residual"
  )) |> 
  mutate(random_mat = as.factor(random_mat)) |> 
  mutate(random_mat = factor(random_mat, levels = c("GCA", "SCA", "GCAxE", "SCAxE", "GCAxN", "SCAxN", "Intra block", "Residual")))

var_comp |> 
  ggplot(aes(x = fct_inorder(trait), y = var_percentage, fill = random_mat)) +
  geom_col(color = "black", linewidth = .3, alpha = .6) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(n.breaks = 10) +
  xlab("Trait") +
  ylab("Percentage of phenotypic variance") +
  labs(fill = "Variance component") +
  theme_classic()+
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15),
    axis.title = element_text(size = 15),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12)
  ) -> p1
p1

# Heritability
herit <- var_comp |> 
  group_by(trait) |> 
  mutate(
    total_var = sum(var),
    add = 2 * var,
    genetic = ifelse(
      random_mat == "GCA", var * 2,
      ifelse(random_mat == "SCA", var * 4, 0)
    ),
    genetic = sum(genetic)
  ) |> 
  filter(random_mat == "GCA") |> 
  mutate(
    H2 = genetic / total_var,
    h2 = add / total_var
  )

p2 <- herit |> 
  tidyr::pivot_longer(
    cols = c(H2, h2),
    names_to = "heritability",
    values_to = "value"
  ) |> 
  ggplot(aes(
    x = forcats::fct_inorder(trait),
    y = value,
    fill = heritability
  )) +
  
  geom_col(
    position = "identity",
    color = "black",
    linewidth = .3,
    alpha = .6
  ) +
  
  # Broad-sense H2: above the bar
  geom_text(
    data = \(x) dplyr::filter(x, heritability == "H2"),
    aes(label = sprintf("%.2f", value)),
    vjust = -0.5,
    size = 4
  ) +
  
  # Narrow-sense h2: inside the green bar
  geom_text(
    data = \(x) dplyr::filter(x, heritability == "h2"),
    aes(label = sprintf("%.2f", value)),
    vjust = 1.5,
    size = 4
  ) +
  
  ylab("Heritability") +
  labs(fill = "Heritability") +
  scale_fill_manual(
    values = c("#5aae61", "#9970ab"),
    labels = c("Narrow-sense", "Broad-sense")
  ) +
  scale_y_continuous(
    n.breaks = 6,
    limits = c(0, 1),
    expand = expansion(mult = c(0, .08))
  ) +
  theme_classic() +
  theme(
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 15),
    axis.title.x = element_blank(),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 11),
    axis.ticks.x = element_blank()
  )

p2

# Combine plots
p3 <- plot_grid(p2, p1, align = 'v', nrow = 2, rel_heights = c(4,7))
p3
ggsave("../figures/variance_components.png", p3, width = 20, height = 15, units = "cm")

### Baseline Cross-Validation (Full Model) ----
# Get BLUES for all traits
adj.means <- NULL
ctr = 0
for (trait in trait_list) {
  ctr = ctr + 1
  # Adjusted means
  lmm <- lmer(trait ~ P1P2 - 1 + cov + (1|E_N_B), data = pheno_dmu |> 
                dplyr::select(1:15, trait = {{trait}}) |> 
                mutate(P1P2 = as.factor(P1P2),
                       cov = as.factor(cov)) |> 
                filter(!(trait == -999)))
  adj.tmp <- data.frame(tmp_name = fixef(lmm)) |> 
    rownames_to_column("P1P2") |> 
    filter(!str_detect(P1P2, "^cov")) |> 
    mutate(P1P2 = gsub("P1P2", "", P1P2)) |> 
    dplyr::rename(!!trait := tmp_name)
  
  if (ctr == 1) {
    adj.means <- adj.tmp
  } else {
    adj.means <- adj.means |> left_join(adj.tmp, by = "P1P2")
  }
}
adj.means.final <- adj.means |> 
  left_join(pheno_dmu |> dplyr::select(P1, P2, P1P2) |> distinct(), by = "P1P2") |> 
  dplyr::select(P1, P2, P1P2, everything())

# Make repeated k-folds
n_folds = 5
n_reps = 10
set.seed(1234)
fold_lists <- replicate(n_reps, {
  pheno_dmu %>%
    distinct(P1P2) %>%
    sample_n(size = n(), replace = FALSE) %>%
    mutate(fold = cut(row_number(), breaks = n_folds, labels = FALSE)) %>%
    group_split(fold) %>%
    map(~ .x$P1P2)
}, simplify = FALSE)

cl <- makeCluster(4)
registerDoParallel(cl)
out_all <- NULL
for (i in 1:n_reps) {
  for (j in 1:n_folds) {
    
    test_set <- fold_lists[[i]][[j]]
    
    
    out <- foreach(k = 1:8, .packages = c("dplyr", "data.table", "glue", "lme4", "tibble", "stringr"), .combine = "rbind") %dopar% {
      
      trait <- trait_list[k]
      
      train_df <- pheno_dmu |>
        dplyr::select(1:15, trait = {{trait}}) |> 
        mutate(trait = ifelse(P1P2 %in% test_set, -999, trait))
      
      glue_directive <- glue("
        $ANALYSE 11 9 0 0
        
        $DATA  ASCII (15, 1, -999) data_frame.dat
         
        $VARIABLE
        {glue_collapse(colnames(train_df),  sep = ' ')}
        
        $MODEL
        1
        0
        1 0 13 0 2 3 4 5 6 7 8 9 10 11 12 13
        10           1 2+2 3 4+4 5  6+6   7 
        0
        0
      
        $VAR_STR 2 GREL ASCII GRM1.txt
        $VAR_STR 3 GREL ASCII GRM2.txt
        
        $PRIOR prior.txt
      
        ")
      
      fm <- run_dmu(y = train_df, 
                    trait = trait, 
                    G1 = G_inv_parents, 
                    G2 = G_inv_crosses, 
                    dmu4_prior = out_full,
                    save_lst = F, 
                    dir = glue_directive)
      
      gca <- fm$sol |> 
        filter(V1 == 3 & V3 == 1 & V4 == 2) |> 
        dplyr::select(parent = V5, gca = V8)
      sca <- fm$sol |> 
        filter(V1 == 3 & V3 == 1 & V4 == 3) |> 
        dplyr::select(P1P2 = V5, n_obs = V6,  sca = V8) |> 
        mutate(P1P2 = as.character(P1P2))
      
      solutions <- adj.means.final |> 
        dplyr::select(P1, P2, P1P2, {{trait}}) |> 
        left_join(gca, by = c("P1" = "parent")) |> 
        rename(gca_P1 = gca) |> 
        left_join(gca, by = c("P2" = "parent")) |> 
        rename(gca_P2 = gca) |> 
        left_join(sca, by = "P1P2") |> 
        mutate(cross_pred = gca_P1 + gca_P2 + sca) |> 
        filter(n_obs == 0)
      
      c <- cor(solutions[ , trait], solutions$cross_pred)
      return(data.frame(trait = trait, rep = i, fold = j, corr = c))
      
    }
    
    out_all <- rbind(out_all, out)
    
  }
  
}
# Stop cluster
stopCluster(cl)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Model Fitting (Reduced LMM: GCA only) ----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Reduced model ----
trait_list <- c("sumDMY", "ADF","ADL", "DMDig", "NDF", "NDFD", "Prot", "WSC" )
out_reduced <- list() # No SCA and SCA interactions
for (i in 1:length(trait_list)) {
  
  trait = trait_list[i]
  
  glue_directive <- glue("
  $ANALYSE 1 1 0 0
  
  $DATA  ASCII (15, 8, -999) data_frame.dat
   
  $VARIABLE
  {glue_collapse(colnames(pheno_dmu),  sep = ' ')}


  $MODEL
  1
  0
  {i} 0 10 0 2 3 4 5 6 8 9 11 12
  7              1 2+2 3+3 4+4
  0
  0

  $VAR_STR 2 GREL ASCII GRM1.txt
  $RESIDUALS ASCII

  $DMUAI
  10
  1.0d-7
  1.0d-6
  1
  0
  0
  ")
  
  pheno_in <- pheno_dmu |> 
    filter(!!sym(trait) != -999)
  
  fm <- run_dmu(y = pheno_in, 
                trait = trait, 
                G1 = G_inv_parents,
                save_lst = T, 
                dir = glue_directive)
  out_reduced[[trait]] <- fm
  
}
saveRDS(out_reduced, "../results/reduced_model_outputs.RDS")
if (!exists("out_reduced") || is.null(out_reduced) || identical(out_reduced, list())) {
  out_reduced <- readRDS("../results/reduced_model_outputs.RDS")
}

### Cross-validation ----
# Make repeated k-folds
n_folds = 5
n_reps = 10
set.seed(1234) # Same seed as before to make the CVs more comparable
fold_lists <- replicate(n_reps, {
  pheno_dmu %>%
    distinct(P1P2) %>%
    sample_n(size = n(), replace = FALSE) %>%
    mutate(fold = cut(row_number(), breaks = n_folds, labels = FALSE)) %>%
    group_split(fold) %>%
    map(~ .x$P1P2)
}, simplify = FALSE)

# Run CV
cl <- makeCluster(4)
registerDoParallel(cl)
out_all_reduced <- NULL
for (i in 1:n_reps) {
  for (j in 1:n_folds) {
    
    test_set <- fold_lists[[i]][[j]]
    
    out <- foreach(k = 1:8, .packages = c("dplyr", "data.table", "glue", "lme4", "tibble", "stringr"), .combine = "rbind") %dopar% {
      
      trait <- trait_list[k]
      
      train_df <- pheno_dmu |>
        dplyr::select(1:15, trait = {{trait}}) |> 
        mutate(trait = ifelse(P1P2 %in% test_set, -999, trait))
      
      glue_directive <- glue("
        $ANALYSE 11 9 0 0
        
        $DATA  ASCII (15, 1, -999) data_frame.dat
         
        $VARIABLE
        {glue_collapse(colnames(train_df),  sep = ' ')}
        
        $MODEL
        1
        0
        1 0 10 0 2 3 4 5 6 8 9 11 12
        7            1 2+2 3+3 4+4
        0
        0
      
        $VAR_STR 2 GREL ASCII GRM1.txt
        
        $PRIOR prior.txt
        ")
      
      fm <- run_dmu(y = train_df, 
                    trait = trait, 
                    G1 = G_inv_parents, 
                    dmu4_prior = out_reduced,
                    save_lst = F, 
                    dir = glue_dir_reduced)
      
      gca <- fm$sol |> 
        filter(V1 == 3 & V3 == 1 & V4 == 2) |> 
        dplyr::select(parent = V5, gca = V8)
      
      solutions <- adj.means.final |> 
        dplyr::select(P1, P2, P1P2, {{trait}}) |> 
        left_join(gca, by = c("P1" = "parent")) |> 
        rename(gca_P1 = gca) |> 
        left_join(gca, by = c("P2" = "parent")) |> 
        rename(gca_P2 = gca) |> 
        mutate(cross_pred = gca_P1 + gca_P2) |> 
        filter(P1P2 %in% test_set)
      
      c <- cor(solutions[ , trait], solutions$cross_pred)
      return(data.frame(trait = trait, rep = i, fold = j, corr = c))
      
    }
    
    out_all_reduced <- rbind(out_all_reduced, out)
    
  }
  
}
# Stop cluster
stopCluster(cl)
# Save CV results to the disk 
out_all$model <- "complete"
out_all_reduced$model <- "reduced"
cv_full_data_final <- rbind(out_all, out_all_reduced)
write.table(cv_full_data_final, "../results/cv_full_data_final.csv", sep = ",", row.names = F, quote = F)
cv_full_data_final <- read.table("../results/cv_full_data_final.csv", h = T, sep = ",")

### Plot CV ----
cv_full_data_final |> 
  ggplot(aes(x = trait, y = corr)) +
  geom_violin()
colours <- c("dodgerblue2", "darkorange")

cv_full_data_final |> 
  mutate(model = ifelse(model == "complete", "LMM2: GCA+SCA", "LMM1: GCA"),
         trait = ifelse(trait == "sumDMY", "DMY", trait)) |> 
  ggplot(aes(x = fct_inorder(trait), y = corr, fill = model)) +
  introdataviz::geom_split_violin(alpha = .6) +
  geom_boxplot(width = .2, alpha = .6, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175), size = .1) +
  scale_fill_manual(values = colours, name = "Parental model") +
  scale_y_continuous(name = "Predictive ability", n.breaks = 15) +
  xlab("Trait") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 11),
    panel.grid = element_blank(),
    axis.ticks.x = element_line(linewidth = .5),
    axis.ticks.y = element_line(linewidth = .5)
  ) -> p4
p4
ggsave(plot = p4, filename = "../figures/prediction_acc_full_dataset.png", width = 18, height = 7, units = "cm")

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Second Cross-Validation Scenario (Relatedness T0, T1, T2) ----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Get pedigree structure
pedigree <- pheno_dmu |> 
  dplyr::select(P1, P2, P1P2) |> 
  group_by(P1P2) |> 
  slice_head(n = 1) |> 
  as.data.frame() 

# Run CV
n_rep <- 25 # Number of repetitions of the random cross-validation scheme
cl <- makeCluster(4)
registerDoParallel(cl)
trait_list <- c("sumDMY", "ADF","ADL", "DMDig", "NDF", "NDFD", "Prot", "WSC" )
relationship_cv <- NULL
for (int in 1:n_rep) {
  
  # T2
  set.seed(int)
  spl1 <- sample(pedigree$P1P2, 30)
  test.t2 <- pedigree |> 
    filter(P1P2 %in% spl1)
  train.t2 <- pedigree |> 
    filter(!(P1P2 %in% test.t2$P1P2)) |> 
    filter(!(P1 %in% c(test.t2$P1, test.t2$P2) & !(P2 %in% c(test.t2$P1, test.t2$P2)))) |> 
    slice_sample(n = 202, replace = F) |> 
    pull(P1P2)
  test.t2 <- test.t2 |> pull(P1P2)
  
  # T1
  set.seed(int)
  spl2 <- sample(1:104, 27) # 27 to have a training set of around 200 families
  at_least_one <- pedigree |> 
    filter(P1 %in% spl2 | P2 %in% spl2)
  train.t1 <- pedigree |> 
    filter(!(P1P2 %in% at_least_one$P1P2)) |> 
    pull(P1P2)
  test.t1 <- at_least_one |> 
    filter(!(P1 %in% spl2) | !(P2 %in% spl2)) |> 
    slice_sample(n = 30, replace = F) |> 
    pull(P1P2)
  
  # T0
  set.seed(int)
  spl3 <- sample(1:104, 28) # 28 to have a training set of around 200 families
  test.t0 <- pedigree |> 
    filter(P1 %in% spl3 & P2 %in% spl3) |> 
    pull(P1P2)
  train.t0 <- pedigree |> 
    filter(!(P1 %in% spl3) & !(P2 %in% spl3)) |> 
    pull(P1P2)
  
  for (cv in c("t2", "t1", "t0")) {
    
    # Training set
    list_families_train <- get(paste("train", cv, sep = "."))
    list_families_test <- get(paste("test", cv, sep = "."))
    
    out_par <- foreach(k = 1:8, .packages = c("dplyr", "data.table", "glue", "tibble"), .combine = "rbind") %dopar% {
      
      trait <- trait_list[k]
      
      train_df <- pheno_dmu |>
        dplyr::select(1:15, trait = {{trait}}) |> 
        filter(P1P2 %in% list_families_train)
      
      # %%%%%%%%%%%%%%%%%% Complete %%%%%%%%%%%%%%%%%%%
      
      glue_dir_complete <- glue("
        $ANALYSE 11 9 0 0
        
        $DATA  ASCII (15, 1, -999) data_frame.dat
         
        $VARIABLE
        {glue_collapse(colnames(train_df),  sep = ' ')}
        
        $MODEL
        1
        0
        1 0 13 0 2 3 4 5 6 7 8 9 10 11 12 13
        10           1 2+2 3 4+4 5  6+6   7 
        0
        0
      
        $VAR_STR 2 GREL ASCII GRM1.txt
        $VAR_STR 3 GREL ASCII GRM2.txt
        
        $PRIOR prior.txt
      
        ")
      
      fm_complete <- run_dmu(y = train_df, 
                             trait = trait, 
                             G1 = G_inv_parents, 
                             G2 = G_inv_crosses, 
                             dmu4_prior = out_full,
                             save_lst = F, 
                             dir = glue_dir_complete)
      
      gca <- fm_complete$sol |> 
        filter(V1 == 3 & V3 == 1 & V4 == 2) |> 
        dplyr::select(parent = V5, gca = V8)
      sca <- fm_complete$sol |> 
        filter(V1 == 3 & V3 == 1 & V4 == 3) |> 
        dplyr::select(P1P2 = V5, n_obs = V6,  sca = V8) |> 
        mutate(P1P2 = as.character(P1P2))
      
      solutions <- adj.means.final |> 
        dplyr::select(P1, P2, P1P2, {{trait}}) |> 
        left_join(gca, by = c("P1" = "parent")) |> 
        rename(gca_P1 = gca) |> 
        left_join(gca, by = c("P2" = "parent")) |> 
        rename(gca_P2 = gca) |> 
        left_join(sca, by = "P1P2") |> 
        mutate(cross_pred = gca_P1 + gca_P2 + sca) |> 
        filter(P1P2 %in% list_families_test)
      
      c <- cor(solutions[, trait], solutions$cross_pred)
      df_complete <- data.frame(model = "complete", CV_scheme = cv, trait = trait, replicate = int, corr = c)
      
      # %%%%%%%%%%%%%%%%%% Reduced %%%%%%%%%%%%%%%%%%%
      glue_dir_reduced <- glue("
        $ANALYSE 11 9 0 0
        
        $DATA  ASCII (15, 1, -999) data_frame.dat
         
        $VARIABLE
        {glue_collapse(colnames(train_df),  sep = ' ')}
        
        $MODEL
        1
        0
        1 0 10 0 2 3 4 5 6 8 9 11 12
        7            1 2+2 3+3 4+4
        0
        0
      
        $VAR_STR 2 GREL ASCII GRM1.txt
        
        $PRIOR prior.txt
      
        ")
      
      fm_reduced <- run_dmu(y = train_df, 
                            trait = trait, 
                            G1 = G_inv_parents, 
                            dmu4_prior = out_reduced,
                            save_lst = F, 
                            dir = glue_dir_reduced)
      
      gca_red <- fm_reduced$sol |> 
        filter(V1 == 3 & V3 == 1 & V4 == 2) |> 
        dplyr::select(parent = V5, gca = V8)
      
      solutions_red <- adj.means.final |> 
        dplyr::select(P1, P2, P1P2, {{trait}}) |>
        left_join(gca_red, by = c("P1" = "parent")) |> 
        rename(gca_P1 = gca) |> 
        left_join(gca_red, by = c("P2" = "parent")) |> 
        rename(gca_P2 = gca) |> 
        mutate(cross_pred = gca_P1 + gca_P2) |> 
        filter(P1P2 %in% list_families_test)
      
      c_red <- cor(solutions_red[,trait], solutions_red$cross_pred)
      df_reduced <- data.frame(model = "reduced", CV_scheme = cv, trait = trait, replicate = int, corr = c_red)
      
      return(rbind(df_complete, df_reduced))
      
      Sys.sleep(.5)
      
    }
    
    relationship_cv <- rbind(relationship_cv, out_par)
    
  }
  
}
# Stop cluster
stopCluster(cl)
# Write to the disk
write.table(relationship_cv, "../results/cv_t2t1t0.csv", sep = ",", row.names = F, quote = F)
#relationship_cv <- read.table("../results/cv_t2t1t0.csv", h = T, sep = ",")

## Plot CV ----
stat_df <- relationship_cv |> 
  group_by(trait, CV_scheme, model) |> 
  summarise(corr_mean = mean(corr), sd = sd(corr)) |> 
  mutate(trait = ifelse(trait == "sumDMY", "DMY", trait))

relationship_cv |> 
  mutate(trait = ifelse(trait == "sumDMY", "DMY", trait)) |>
  ggplot(aes(x = CV_scheme, y = corr, fill = CV_scheme)) +
  geom_point(aes(shape = model), position = position_jitterdodge(jitter.width = 0.1, seed = 123), alpha = .2, size = 1) +
  geom_boxplot(aes(group = interaction(model, CV_scheme)), color = "black", outliers = F, alpha = .6, linewidth = .3) +
  scale_fill_manual(values = c('darkgray',"dodgerblue2", "darkorange"), name = "CV scheme", labels = c("T0", "T1", "T2")) +
  scale_shape_manual(values = c(1,2), labels = c("LMM1: GCA", "LMM2: GCA+SCA"), name = "Parental model") +
  xlab("Trait") +
  scale_y_continuous(name = "Predictive ability", n.breaks = 15) +
  guides(shape = guide_legend(override.aes = list(size = 1))) +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  facet_wrap(~ fct_inorder(trait), ncol = 8, strip.position = "bottom") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 10),
    panel.grid = element_blank(),
    axis.ticks.y = element_line(linewidth = .5),
    axis.text.x = element_blank(),
    legend.title=element_text(size=10), 
    legend.text=element_text(size=9),
    strip.text = element_text(size = 12)
  )-> p5
p5
ggsave(plot = p5, filename = "../figures/prediction_acc_relationship.png", width = 18, height = 7, units = "cm")

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Correlations & Nitrogen Effect Figures ----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X <- adj.means.final[ , -c(1:3)]
Xc <- cor(X)
Xp <- pcor(X)$estimate
corr_m <- matrix(NA, ncol = ncol(X), nrow = ncol(X))
corr_m[upper.tri(corr_m)] <- Xc[upper.tri(corr_m)]
corr_m[lower.tri(corr_m)] <- Xp[lower.tri(corr_m)]
diag(corr_m) <- NA
colnames(corr_m) = rownames(corr_m) = colnames(X)

melted_cormat <- reshape2::melt(corr_m) |> 
  mutate(
    Var1 = recode(Var1, sumDMY = "DMY"),
    Var2 = recode(Var2, sumDMY = "DMY")
  ) |> 
  mutate(value = round(value, 2))
# Heatmap
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "black")+
  scale_fill_gradient2(low = "#d73027", high = "#1a9850", mid = "#ffffbf", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", na.value = 'white',
                       name="Correlation\ncoefficient") +
  geom_text(aes(label = value), color = "black", size = 2) +
  xlab("Simple correlation") +
  ylab("Partial correlation") +
  theme_minimal() +
  coord_fixed() +
  theme(
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", colour = "white"), 
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.title.y = element_text(hjust = 1),
    axis.title.x = element_text(hjust = 1)
  ) -> p6
p6
ggsave(plot = p6, filename = "../figures/Figure_corr_plus_BLUE/correlation.svg", width = 10, height = 12, units = "cm", dpi = 500)

# Nitrogen effect Fig
env_BLUE <- read.table("../results/E_BLUEs_of_environments.txt", h = T, sep = "\t") |> 
  mutate(trait = fct_inorder(trait))

range_data <- env_BLUE %>%
  group_by(trait) |> 
  summarize(
    min_val = min(blue - se, na.rm = TRUE) * .9,
    max_val = max(blue + se, na.rm = TRUE) * 1.1
  ) |> 
  as.data.frame()

env_BLUE |> 
  mutate(manag = as.factor(manag),
         code = as.factor(code)) |> 
  ggplot(aes(x = manag, y = blue, group = yl)) +
  geom_pointrange(aes(ymin = blue-se, ymax = blue + se, colour = yl), position = position_dodge(.5), size = .2) +
  scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")) +
  labs(color = "Year - Location", x = "Nitrogen managment", y = "BLUE") +
  theme_bw() +
  theme(legend.position = c(.85, 0.15),
        legend.key.size = unit(.7, 'cm'),
        legend.title = element_text(colour = "black", size = 18),
        legend.text = element_text(colour = "black", size = 18),
        legend.spacing.x = unit(0.1, 'cm'),
        axis.text.x = element_text(size=16, colour = 'black'), 
        axis.title.x = element_text(size= 18), 
        panel.grid = element_blank(),
        axis.text.y = element_text(size=13, colour = 'black'), 
        axis.title.y = element_text(size= 18),  
        panel.spacing = unit(1, "lines"),
        strip.background = element_rect(color = NA, fill = NA),
        strip.text.x = element_text(size = 16, colour = "black", face = 'bold', hjust = 0)) +
  facet_wrap( ~ trait, scales = "free") -> p7
p7 + ggh4x::facetted_pos_scales(
  y = list(trait == "sumDMY" ~ scale_y_continuous(limits = c(range_data[1,2], range_data[1,3]), n.breaks = 6),
           trait == "ADF" ~ scale_y_continuous(limits = c(range_data[2,2], range_data[2,3]), n.breaks = 6),
           trait == "ADL" ~ scale_y_continuous(limits = c(range_data[3,2], range_data[3,3]), n.breaks = 6),
           trait == "DMDig" ~ scale_y_continuous(limits = c(range_data[4,2], range_data[4,3]), n.breaks = 6),
           trait == "NDF" ~ scale_y_continuous(limits = c(range_data[5,2], range_data[5,3]), n.breaks = 6),
           trait == "NDFD" ~ scale_y_continuous(limits = c(range_data[6,2], range_data[6,3]), n.breaks = 6),
           trait == "Prot" ~ scale_y_continuous(limits = c(range_data[7,2], range_data[7,3]), n.breaks = 6),
           trait == "WSC" ~ scale_y_continuous(limits = c(range_data[8,2], range_data[8,3]), n.breaks = 6))
) -> p8
p8

ggsave('../figures/BLUEs.png', plot = p8, width = 20, height = 15, units = 'cm')


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GWAS analysis ----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freqAvg <- colMeans(GD_tmp, na.rm=T)
p <- freqAvg / 2
M <- t(GD_tmp) - 1         # Re-code matrix as -1, 0, 1
P <- 2 * (p - 0.5)         # Deviation from 0.5
W <- M - P
WtW = t(W) %*% W
d <- 2 * sum(p * (1 - p))  # Scaling
G <- WtW / d               # GRM

# Load solutions from the parental model
if (!exists("out_full") || is.null(out_full) || identical(out_full, list())) {
  out_full <- readRDS("../results/full_model_outputs.RDS")
}

# Load physical map
map <- data.table::fread("../data/greenSelectParents_QUAL20_MQ10_maf005_miss50_imputed_GM.csv")

# Run and plot GWAS
source("https://raw.githubusercontent.com/YinLiLin/CMplot/master/R/CMplot.r") # Function for GWAS plotting
dir_path <- "../figures/gwas"
if(!dir.exists(dir_path)) dir.create(dir_path)
here <- getwd()
setwd(dir_path)
trait_list <- c("sumDMY", "ADF","ADL", "DMDig", "NDF", "NDFD", "Prot", "WSC" )
gwass_all <- data.frame()
trait_GCA <- data.frame()
for (i in trait_list) {
  gca_hat <- out_full |> 
    pluck(i, "sol") |> 
    filter(V1 == 3 & V3 == 1 & V4 == 2) |> 
    dplyr::select(parent = V5, gca = V8) |> 
    column_to_rownames("parent") |> 
    as.matrix()
  gca_hat <- gca_hat[match(colnames(M), rownames(gca_hat)), , drop = FALSE]
  
  if (i == "sumDMY") {
    trait_GCA <- gca_hat |> as.data.frame() |> rename("{i}" := gca) |>  rownames_to_column("parent")
  } else {
    trait_GCA <- trait_GCA |> left_join(gca_hat |> as.data.frame() |> rename("{i}" := gca) |> rownames_to_column("parent"), by = "parent")
  }
  
  # Back-solving for SNP effects
  backSolve <- 1 / d * W %*% solve(G) %*% gca_hat
  
  # Approximate p-value based on a t-distribution
  pvalBlup <- 2 * pt(-abs(backSolve / sd(backSolve)), df = length(backSolve) - 1)
  gwas <- data.frame(map, pval = pvalBlup[,1])
  
  # Threshold
  BF <- 0.05 / nrow(M) # Bonferroni
  p_values_fdr <- max(gwas$pval[p.adjust(gwas$pval, method = "fdr") <= 0.05]) # FDR
  
  
  # Significant SNPs
  SNPs <- list(
    gwas$Name[gwas$pval <= p_values_fdr]
  )
  
  # Plot Manhattan
  if (FALSE) {
    CMplot(gwas,
           type = "p", 
           plot.type = "m",
           col=c("dodgerblue2","grey60"),
           LOG10 = TRUE,
           threshold=c(BF,p_values_fdr),
           threshold.lty=c(1,2),
           threshold.lwd=c(1,2),
           threshold.col=c("red","blue"),
           file = "jpg",
           file.name = i,
           dpi = 300,
           file.output = TRUE,
           width = 10,
           height = 5,
           highlight = SNPs,
           highlight.text=SNPs
    )
  }
  
  if (i == "sumDMY") {
    SNPs_all <- gwas$Name[gwas$pval <= BF]
    gwass_all <- gwas |> 
      rename(!!i := pval)
    rownames(gwass_all) <- NULL
  } else {
    SNPs_all <- append(SNPs_all, gwas$Name[gwas$pval <= BF])
    gwass_all <- gwass_all |> 
      left_join(gwas |> dplyr::select(Name, !!i := pval), by = "Name")
  }
  
  
}

CMplot(gwass_all, 
       plot.type="m",
       multraits=TRUE,
       threshold=BF,
       bin.size=1e6,
       chr.den.col=c("#01665e", "#f5f5f5", "#8c510a"),
       col="#d9d9d9",
       signal.col=c("#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0", "#f0027f", "#bf5b17", "#666666"),
       amplify=TRUE,
       signal.cex=1, 
       threshold.col="black",
       threshold.lwd=.5,
       threshold.lty=1,
       points.alpha=200,
       file="jpg",
       file.name="multi_trait",
       dpi=300,
       file.output=TRUE,
       verbose=TRUE,
       legend.ncol=1,
       legend.pos="left")
setwd(here)
