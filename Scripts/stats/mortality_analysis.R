# mortality_statistics.R
# UVBR + temperature effects on embryo survival
# ONLY STATISTICAL ANALYSES

rm(list = ls())

# =========================
# 1. PACKAGES
# =========================
required_packages <- c("dplyr", "lme4", "emmeans")

to_install <- required_packages[!required_packages %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) {
  install.packages(to_install)
}

invisible(lapply(required_packages, library, character.only = TRUE))

# =========================
# 2. LOAD DATA
# =========================
ruta <- "C:/Users/KGH/OneDrive - ut.edu.co/DCB/Paper UVB+temp/Herpetological Conservation and Biology/New submission/mortality_counts.csv"

file.exists(ruta)

if (!file.exists(ruta)) {
  stop("The file was not found. Please check the path in 'ruta'.")
}

mortality <- read.csv(ruta, stringsAsFactors = FALSE)

# Optional quick check
str(mortality)
head(mortality)

# =========================
# 3. PREPARE VARIABLES
# =========================
mortality <- mortality %>%
  mutate(
    species = factor(
      species,
      levels = c(
        "Boana_platanera",
        "Engystomops_pustulosus",
        "Rhinella_horribilis"
      )
    ),
    uvbr = factor(uvbr, levels = c(0, 25, 50)),
    temperature = factor(temperature, levels = c(23, 28, 33)),
    replicate = factor(replicate),
    clutch = factor(clutch),
    survival = 1 - mortality
  )

str(mortality)
levels(mortality$species)
levels(mortality$uvbr)
levels(mortality$temperature)
levels(mortality$replicate)

# =========================
# 4. CHECK DESIGN AND OBSERVED MORTALITY
# =========================
design_check <- mortality %>%
  count(species, uvbr, temperature, replicate)

design_check

summary_mort <- mortality %>%
  group_by(species, uvbr, temperature) %>%
  summarise(
    n = n(),
    dead = sum(mortality),
    alive = sum(survival),
    mort_rate = mean(mortality),
    surv_rate = mean(survival),
    .groups = "drop"
  )

summary_mort

zero_uv_check <- mortality %>%
  filter(uvbr == "0") %>%
  group_by(species, temperature) %>%
  summarise(
    n = n(),
    dead = sum(mortality),
    alive = sum(survival),
    .groups = "drop"
  )

zero_uv_check

# =========================
# 5. PRIMARY INFERENTIAL MODEL
#    Excluding UVBR = 0 because all embryos survived
# =========================
mort_model <- mortality %>%
  filter(uvbr != "0") %>%
  droplevels()

summary(mort_model)
levels(mort_model$uvbr)

model_full <- glmer(
  survival ~ species * uvbr * temperature + (1 | replicate),
  family = binomial(link = "logit"),
  data = mort_model,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

summary(model_full)

# =========================
# 6. MODEL DIAGNOSTICS
# =========================
singular_result <- isSingular(model_full, tol = 1e-4)
singular_result

VarCorr(model_full)

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  ratio <- Pearson.chisq / rdf
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  
  data.frame(
    chisq = Pearson.chisq,
    ratio = ratio,
    rdf = rdf,
    p_value = pval
  )
}

overdisp_results <- overdisp_fun(model_full)
overdisp_results

model_no3 <- update(model_full, . ~ . - species:uvbr:temperature)
lrt_three_way <- anova(model_full, model_no3, test = "Chisq")
lrt_three_way

# =========================
# 7. POST HOC COMPARISONS
# =========================
emm_temp <- emmeans(
  model_full,
  ~ temperature | species * uvbr,
  type = "response"
)

emm_temp

temp_pairs <- as.data.frame(
  pairs(emm_temp, adjust = "tukey")
)

temp_pairs

emm_uv <- emmeans(
  model_full,
  ~ uvbr | species * temperature,
  type = "response"
)

emm_uv

uv_pairs <- as.data.frame(
  pairs(emm_uv, adjust = "holm")
)

uv_pairs

# =========================
# 8. PREDICTED VALUES FROM MODEL
# =========================
emm_cells <- emmeans(
  model_full,
  ~ species * uvbr * temperature,
  type = "response"
)

emm_cells_df <- as.data.frame(emm_cells)

emm_cells_df <- emm_cells_df %>%
  mutate(
    pred_mortality = 1 - prob,
    mort_LCL = 1 - asymp.UCL,
    mort_UCL = 1 - asymp.LCL
  )

emm_cells_df

# =========================
# 9. OPTIONAL MODEL INCLUDING UVBR = 0
# =========================
model_with0 <- try(
  glmer(
    survival ~ species * uvbr * temperature + (1 | replicate),
    family = binomial(link = "logit"),
    data = mortality,
    control = glmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 2e5)
    )
  ),
  silent = TRUE
)

if (!inherits(model_with0, "try-error")) {
  summary(model_with0)
}
