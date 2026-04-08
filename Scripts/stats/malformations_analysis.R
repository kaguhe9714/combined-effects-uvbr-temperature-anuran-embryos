# malformations_analysis.R
# UVBR + temperature effects on embryo malformations


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
ruta <- "C:/Users/KGH/OneDrive - ut.edu.co/DCB/Paper UVB+temp/Herpetological Conservation and Biology/New submission/malformations_counts.csv"

file.exists(ruta)

if (!file.exists(ruta)) {
  stop("The file was not found. Please check the path in 'ruta'.")
}

malf <- read.csv(ruta, stringsAsFactors = FALSE)


# =========================
# 3. CHECK DATA
# =========================
str(malf)
head(malf)
names(malf)
dim(malf)

sort(unique(malf$uvbr))
sort(unique(malf$temperature))
table(malf$uvbr, useNA = "ifany")
table(malf$temperature, useNA = "ifany")


# =========================
# 4. PREPARE VARIABLES
# =========================
malf <- malf %>%
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
    malformation = as.integer(malformation)
  )

str(malf)


# =========================
# 5. INFERENTIAL DATASET
# =========================
# UVBR = 0 is excluded from the inferential model
# to avoid complete separation, but can still be used in figures

malf_model <- malf %>%
  filter(uvbr != "0") %>%
  droplevels()

summary(malf_model)
levels(malf_model$uvbr)


# =========================
# 6. MAIN GLMM
# =========================
model_malf_full <- glmer(
  malformation ~ species * uvbr * temperature + (1 | replicate),
  family = binomial(link = "logit"),
  data = malf_model,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

summary(model_malf_full)


# =========================
# 7. MODEL CHECKS
# =========================
isSingular(model_malf_full, tol = 1e-4)
VarCorr(model_malf_full)

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  ratio <- Pearson.chisq / rdf
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  
  c(
    chisq = Pearson.chisq,
    ratio = ratio,
    rdf = rdf,
    p = pval
  )
}

overdisp_fun(model_malf_full)


# =========================
# 8. TEST THREE-WAY INTERACTION
# =========================
model_malf_no3 <- update(model_malf_full, . ~ . - species:uvbr:temperature)

anova(model_malf_full, model_malf_no3, test = "Chisq")


# =========================
# 9. IF RANDOM EFFECT IS NEGLIGIBLE, FIT GLM
# =========================
model_malf_glm <- glm(
  malformation ~ species * uvbr * temperature,
  family = binomial(link = "logit"),
  data = malf_model
)

summary(model_malf_glm)

model_malf_glm_no3 <- update(model_malf_glm, . ~ . - species:uvbr:temperature)

anova(model_malf_glm_no3, model_malf_glm, test = "Chisq")


# =========================
# 10. POST HOC: TEMPERATURE WITHIN SPECIES × UVBR
# =========================
emm_malf_temp <- emmeans(
  model_malf_glm,
  ~ temperature | species * uvbr,
  type = "response"
)

emm_malf_temp

malf_temp_pairs <- pairs(emm_malf_temp, adjust = "tukey")
malf_temp_pairs


# =========================
# 11. POST HOC: UVBR WITHIN SPECIES × TEMPERATURE
# =========================
emm_malf_uv <- emmeans(
  model_malf_glm,
  ~ uvbr | species * temperature,
  type = "response"
)

emm_malf_uv

malf_uv_pairs <- pairs(emm_malf_uv, adjust = "holm")
malf_uv_pairs


# =========================
# 12. PREDICTED PROBABILITIES
# =========================
emm_malf_cells <- emmeans(
  model_malf_glm,
  ~ species * uvbr * temperature,
  type = "response"
)

emm_malf_cells_df <- as.data.frame(emm_malf_cells)

emm_malf_cells
emm_malf_cells_df
