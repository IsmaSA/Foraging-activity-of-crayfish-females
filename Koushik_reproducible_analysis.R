## ============================================================================
##  Koushik et al. — "An overlooked reproductive stage: Foraging activity of
##  crayfish females with glair glands" (Inland Waters)
##
##  FULL REPRODUCIBLE ANALYSIS SCRIPT
##  Reproduces every result in the manuscript on the complete dataset and adds
##  the new body-size analysis (Table S6) requested by Reviewer 2 (comment 11),
##  the beta-binomial refit + dispersion (comments 5 & 15) and the Figure 3 fix
##  (comment 12). All tables and figures are written to ./outputs_repro/.
##
##  Data files expected in the working directory:
##    - "Tonda_student.xlsx"                              (sheet "Ismael", "Experiment 2")
##    - "FR REPRODUCTION STATUS OF FEMALES complete.xlsx" (sheets "Experiment 1", "Experiment 2")
##
##  R 4.1.3. Package versions this script was written against:
##    readxl 1.4.3 · dplyr 1.1.4 · tidyr 1.3.1 · stringr 1.5.1 · ggplot2 3.5.1
##    patchwork 1.2.0 · frair 0.5.100 · glmmTMB 1.1.9 · DHARMa 0.4.6
##    emmeans 1.11.0 · car 3.1-2 · broom 1.0.6 · broom.mixed 0.2.9.5 · writexl 1.5.0
## ============================================================================

## ---- 0. Reproducibility preamble -------------------------------------------
rm(list = ls())
set.seed(1234)                       # governs frair bootstraps & DHARMa simulation

pkgs <- c("readxl","dplyr","tidyr","stringr","ggplot2","patchwork","frair",
          "glmmTMB","DHARMa","emmeans","car","broom","broom.mixed","writexl")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

## Edit this if the workbooks live elsewhere:
data_dir <- "."
setwd(data_dir)
out_dir  <- file.path(data_dir, "outputs_repro")
dir.create(out_dir, showWarnings = FALSE)

f_fr  <- "Tonda_student.xlsx"                                  # functional-response & partial data
f_raw <- "FR REPRODUCTION STATUS OF FEMALES complete.xlsx"     # per-individual CL / weight / killed

cat("R version:", R.version.string, "\n")
print(sapply(pkgs, function(p) as.character(packageVersion(p))))

## Colour palette used throughout (colour-blind friendly, consistent with the MS)
col_grp <- c("Non-reproductive" = "deepskyblue3",
             "Glair glands"     = "chocolate3",
             "Egg-carrying"     = "forestgreen")

## ============================================================================
##  EXPERIMENT 1 — FUNCTIONAL RESPONSE  (Table 1, Figure 1, Figure 2)
## ============================================================================
## Wide FR sheet: N in cols 1-11 (Total), Glair in 12-20 (Total1), Eggs in 21-29 (Total2).
fr <- read_excel(f_fr, sheet = "Ismael")
fr[is.na(fr)] <- 0
df_N <- fr[, c(1:11)]          # Non-reproductive  -> response column "Total"
df_G <- fr[, c(2, 12:20)]      # Glair glands      -> response column "Total1"
df_E <- fr[, c(2, 21:29)]      # Egg-carrying      -> response column "Total2"

## 1a. Functional-response TYPE test (phi first-order term; Table 1 "Mean"/"SE")
type_N <- frair_test(Total  ~ `Prey density`, data = df_N)
type_G <- frair_test(Total1 ~ `Prey density`, data = df_G)
type_E <- frair_test(Total2 ~ `Prey density`, data = df_E)

## 1b. Fit Rogers' random-predator equation (T fixed = 1; rate expressed per exposure)
fit_N <- frair_fit(Total  ~ `Prey density`, data = df_N, response = "rogersII",
                   start = list(a = 1, h = 0.1), fixed = list(T = 1))
fit_G <- frair_fit(Total1 ~ `Prey density`, data = df_G, response = "rogersII",
                   start = list(a = 1, h = 0.1), fixed = list(T = 1))
fit_E <- frair_fit(Total2 ~ `Prey density`, data = df_E, response = "rogersII",
                   start = list(a = 1, h = 0.1), fixed = list(T = 1))

## 1c. Non-parametric bootstrap CIs for a and h (999 resamples; seed set above)
boot_N <- frair_boot(fit_N, nboot = 999)
boot_G <- frair_boot(fit_G, nboot = 999)
boot_E <- frair_boot(fit_E, nboot = 999)

## 1d. Assemble Table 1 (a, h, 1/hT, FRR) with bootstrap CIs
## Robust extraction: read the bootstrap replicate matrix directly rather than
## indexing confint()'s printed object (whose structure varies across frair versions).
fr_params <- function(fit, boot, group) {
  co <- boot$coefficients                       # original-data fit: named a, h, T
  a  <- as.numeric(co[["a"]]); h <- as.numeric(co[["h"]])
  T  <- if ("T" %in% names(co)) as.numeric(co[["T"]]) else 1
  bc  <- as.data.frame(boot$bootcoefs)          # 999 x {a,h,T} bootstrap replicates
  aci <- quantile(bc$a, c(.025, .975), na.rm = TRUE)   # percentile 95% CI
  hci <- quantile(bc$h, c(.025, .975), na.rm = TRUE)
  data.frame(
    Group   = group,
    a       = round(a, 3),
    a_CI    = paste0(round(aci[1], 3), "-", round(aci[2], 3)),
    h       = round(h, 3),
    h_CI    = paste0(round(hci[1], 3), "-", round(hci[2], 3)),
    max_feed_1hT = round(1/(h*T), 3),     # maximum feeding rate
    FRR     = round(a/h, 3),
    row.names = NULL)
}
## frair's own BCa intervals, printed for reference (wrapped so it cannot halt the run):
print_bca <- function(boot, group) {
  cat("\n[BCa 95% CI] ", group, "\n"); try(print(confint(boot, citypes = "bca")), silent = TRUE)
}
tableS1_FR <- bind_rows(
  fr_params(fit_N, boot_N, "Non-reproductive"),
  fr_params(fit_G, boot_G, "Glair glands"),
  fr_params(fit_E, boot_E, "Egg-carrying"))
cat("\n== Table 1: functional-response parameters (percentile bootstrap CI) ==\n"); print(tableS1_FR)
## frair's native BCa intervals for reference (match those quoted in the MS):
for (b in list(list(boot_N,"Non-reproductive"), list(boot_G,"Glair glands"), list(boot_E,"Egg-carrying")))
  print_bca(b[[1]], b[[2]])

## 1e. Figure 1 — fitted FR curves with data
svg(file.path(out_dir, "Figure1_FR_curves.svg"), width = 7, height = 5)
plot(1, type = "n", xlim = c(0, 26), ylim = c(0, 15),
     xlab = "Prey density", ylab = "No. of prey killed", xaxt = "n", font.lab = 2)
axis(1, at = c(1,3,6,12,18,24))
drawpoly(boot_N, col = "deepskyblue1", border = NA, tozero = TRUE)
drawpoly(boot_G, col = "chocolate1",   border = NA, tozero = TRUE)
drawpoly(boot_E, col = "lightgreen",   border = NA, tozero = TRUE)
lines(fit_N, col = "deepskyblue3", lwd = 3, lty = 1)
lines(fit_G, col = "chocolate3",   lwd = 3, lty = 2)
lines(fit_E, col = "forestgreen",  lwd = 3, lty = 3)
legend(1, 14, legend = names(col_grp), col = col_grp, lty = 1:3, cex = 0.8)
dev.off()

## 1f. Figure 2 — attack rate (a) and handling time (h) with bootstrap CIs
ci_tab <- function(boot, par, group) {
  bc <- as.data.frame(boot$bootcoefs)
  q  <- quantile(bc[[par]], c(.025, .975), na.rm = TRUE)
  data.frame(Group = group, point = as.numeric(boot$coefficients[[par]]),
             lower = as.numeric(q[1]), upper = as.numeric(q[2]))
}
a_df <- bind_rows(ci_tab(boot_N,"a","Non-reproductive"),
                  ci_tab(boot_G,"a","Glair glands"),
                  ci_tab(boot_E,"a","Egg-carrying"))
h_df <- bind_rows(ci_tab(boot_N,"h","Non-reproductive"),
                  ci_tab(boot_G,"h","Glair glands"),
                  ci_tab(boot_E,"h","Egg-carrying"))
lev <- c("Non-reproductive","Glair glands","Egg-carrying")
a_df$Group <- factor(a_df$Group, lev); h_df$Group <- factor(h_df$Group, lev)
p_a <- ggplot(a_df, aes(Group, point, colour = Group)) +
  geom_point(size = 2) + geom_errorbar(aes(ymin = lower, ymax = upper), width = .2) +
  scale_colour_manual(values = col_grp) + ylab("Attack rate (a) ± 95% CI") + xlab("") +
  theme_bw(base_size = 13) + theme(legend.position = "none")
p_h <- ggplot(h_df, aes(Group, point, colour = Group)) +
  geom_point(size = 2) + geom_errorbar(aes(ymin = lower, ymax = upper), width = .2) +
  scale_colour_manual(values = col_grp) + ylab("Handling time (h) ± 95% CI") + xlab("") +
  theme_bw(base_size = 13) + theme(legend.position = "none")
ggsave(file.path(out_dir,"Figure2_a_h.svg"), p_a + p_h, width = 9, height = 4)

## ============================================================================
##  EXPERIMENT 1 — PARTIAL CONSUMPTION  (Table S2; zero-inflated binomial GLM)
## ============================================================================
pc <- read_excel(f_fr, sheet = "Ismael"); pc[is.na(pc)] <- 0
d1 <- pc[, 1:10];  d1 <- d1[,-1]
colnames(d1) <- c("Prey.density","Code","Alive","Killed","p75","p50","p25","p0","Total")
d2 <- pc[, c(2,11:18)]; colnames(d2) <- c("Prey.density","NA","Code","Alive","Killed","p75","p50","p25","p0")
d3 <- pc[, c(2,19:27)]; colnames(d3) <- c("Prey.density","NA1","NA2","Code","Alive","Killed","p75","p50","p25","p0")
d1 <- d1[,-9]; d2 <- d2[,-2]; d3 <- d3[,-c(2,3)]
pcd <- rbind(d1, d2, d3)
## NOTE: original koushik.R summed p25 twice; corrected here to sum each partial class once.
pcd$Partial <- pcd$p75 + pcd$p50 + pcd$p25
pcd$Full    <- pcd$Killed + pcd$Partial
pcd$Code    <- factor(pcd$Code)

m_pc1 <- glmmTMB(cbind(Partial, Full) ~ Code * Prey.density, ziformula = ~1,
                 data = pcd, family = binomial(link = "logit"))
tableS2 <- broom.mixed::tidy(m_pc1)
anova_S2 <- car::Anova(m_pc1, type = 3)
emm_S2  <- emmeans(m_pc1, pairwise ~ Code)
cat("\n== Table S2: Exp 1 partial consumption ==\n"); print(tableS2); print(anova_S2)

## ============================================================================
##  EXPERIMENT 2 — TOTAL PREY KILLED  (Table S4, Figure 4)
##  Comments 5 & 15: bounded count -> binomial(12,p); check dispersion; beta-binomial.
## ============================================================================
## Per-individual killed / CL from the raw workbook (blocks G,E,N in Experiment 2 sheet).
pull_e2 <- function(start_col, group) {
  raw <- read_excel(f_raw, sheet = "Experiment 2", col_names = FALSE, skip = 1)
  tibble(group = group,
         Code  = as.character(raw[[start_col]]),
         CL    = suppressWarnings(as.numeric(raw[[start_col + 3]])),
         W     = suppressWarnings(as.numeric(raw[[start_col + 4]])),
         Alive = suppressWarnings(as.numeric(raw[[start_col + 6]]))) %>%
    filter(!is.na(Code), !is.na(Alive))
}
e2 <- bind_rows(pull_e2(1,"Glair glands"), pull_e2(14,"Egg-carrying"), pull_e2(27,"Non-reproductive")) %>%
  mutate(shelter = factor(ifelse(str_detect(Code, "\\+"), "Shelter", "No shelter"),
                          levels = c("No shelter","Shelter")),
         type    = factor(group, levels = c("Egg-carrying","Glair glands","Non-reproductive")),
         killed  = 12 - Alive) %>%
  filter(killed >= 0, killed <= 12)
cat("\n== Experiment 2 sample sizes ==\n"); print(table(e2$type, e2$shelter))

## Binomial (what the revised MS currently describes)
m_bin <- glmmTMB(cbind(killed, 12 - killed) ~ shelter * type,
                 data = e2, family = binomial(link = "logit"))
## Pearson dispersion = sum of squared Pearson residuals / residual df  (CORRECT wording; comment 15)
pearson_disp <- sum(residuals(m_bin, type = "pearson")^2) / df.residual(m_bin)
cat(sprintf("\nBinomial Pearson dispersion = %.2f\n", pearson_disp))   # ~2.25 -> over-dispersed
print(testDispersion(simulateResiduals(m_bin, n = 1000)))

## Beta-binomial (respects the ceiling of 12 AND extra-binomial variation)
m_bb <- glmmTMB(cbind(killed, 12 - killed) ~ shelter * type,
                data = e2, family = betabinomial(link = "logit"))
cat("\nAIC comparison (binomial vs beta-binomial):\n"); print(AIC(m_bin, m_bb))
tableS4     <- broom.mixed::tidy(m_bb)
emm_shelter <- emmeans(m_bb, pairwise ~ shelter | type, type = "response")
emm_type    <- emmeans(m_bb, pairwise ~ type | shelter, type = "response")
cat("\n== Table S4: Exp 2 total killed (beta-binomial) ==\n"); print(tableS4)
print(emm_shelter); print(emm_type)

## Figure 4 — killed by group x shelter
fig4_df <- e2 %>% group_by(type, shelter) %>%
  summarise(mean = mean(killed), se = sd(killed)/sqrt(n()), .groups = "drop")
p4 <- ggplot(fig4_df, aes(type, mean, fill = type)) +
  geom_col(colour = "black", linewidth = .2) +
  geom_errorbar(aes(ymin = pmax(0, mean - se), ymax = mean + se), width = .2) +
  facet_wrap(~shelter, labeller = as_labeller(c("Shelter"="Shelter present","No shelter"="Shelter absent"))) +
  scale_fill_manual(values = col_grp) + labs(y = "No. of prey killed", x = NULL) +
  theme_bw(base_size = 13) + theme(legend.position = "none",
                                   strip.text = element_text(face = "bold"))
ggsave(file.path(out_dir,"Figure4_killed_shelter.svg"), p4, width = 8, height = 4)

## ============================================================================
##  EXPERIMENT 2 — PARTIAL CONSUMPTION  (Table S5; zero-inflated binomial GLM)
## ============================================================================
pc2 <- read_excel(f_fr, sheet = "Experiment 2")
pc2 <- pc2[, -c(2:6,13,15:19,26,28:32)]
q1 <- pc2[,1:7];  q1 <- q1[-1,]; colnames(q1) <- c("Code","Alive","Killed","p75","p50","p25","p0")
q2 <- pc2[,8:14]; q2 <- q2[-1,]; colnames(q2) <- c("Code","Alive","Killed","p75","p50","p25","p0")
q3 <- pc2[,15:21];q3 <- q3[-1,]; colnames(q3) <- c("Code","Alive","Killed","p75","p50","p25","p0")
pc2d <- rbind(q1,q2,q3) %>% filter(!is.na(Code))
pc2d <- pc2d %>% mutate(across(c(p75,p50,p25,p0,Killed), ~ as.numeric(replace_na(as.character(.),"0"))),
                        Partial = p75 + p50 + p25,
                        Full    = Killed + Partial,
                        Group   = substr(Code, 1, 1))
m_pc2 <- glmmTMB(cbind(Partial, Full) ~ Group, ziformula = ~1,
                 data = pc2d, family = binomial(link = "logit"))
tableS5  <- broom.mixed::tidy(m_pc2)
anova_S5 <- car::Anova(m_pc2, type = 3)
cat("\n== Table S5: Exp 2 partial consumption ==\n"); print(tableS5); print(anova_S5)

## ============================================================================
##  NEW — Reviewer 2, comment 11:  BODY SIZE by reproductive group (Table S6)
##  Uses the FULL recorded set of individuals (every row with a Code and a CL).
##  Set analysed_only <- TRUE to restrict to the 72/group feeding-analysis set.
## ============================================================================
analysed_only <- FALSE
pull_block <- function(sheet, start_col, group) {
  raw <- read_excel(f_raw, sheet = sheet, col_names = FALSE, skip = 1)
  tibble(group = group,
         Code  = as.character(raw[[start_col]]),
         prey  = suppressWarnings(as.numeric(raw[[start_col + 1]])),
         CL    = suppressWarnings(as.numeric(raw[[start_col + 3]])),
         W     = suppressWarnings(as.numeric(raw[[start_col + 4]])),
         Alive = suppressWarnings(as.numeric(raw[[start_col + 6]]))) %>%
    filter(!is.na(Code), !is.na(CL)) %>%
    mutate(W = ifelse(W > 0.3 & W < 4, W, NA_real_))       # drop impossible weights (stray cells)
}
## Experiment 1: blocks N(1), G(14), E(27)
e1_size <- bind_rows(pull_block("Experiment 1", 1,"Non-reproductive"),
                     pull_block("Experiment 1",14,"Glair glands"),
                     pull_block("Experiment 1",27,"Egg-carrying")) %>%
  mutate(killed = prey - Alive)
## Experiment 2: blocks G(1), E(14), N(27)
e2_size <- bind_rows(pull_block("Experiment 2", 1,"Glair glands"),
                     pull_block("Experiment 2",14,"Egg-carrying"),
                     pull_block("Experiment 2",27,"Non-reproductive")) %>%
  mutate(shelter = ifelse(str_detect(Code,"\\+"),"Shelter","No shelter"),
         killed  = 12 - Alive)
if (analysed_only) {
  e1_size <- e1_size %>% group_by(group) %>% slice(1:72) %>% ungroup()
  e2_size <- e2_size %>% group_by(group, shelter) %>% slice(1:26) %>% ungroup()
}
glev <- c("Egg-carrying","Glair glands","Non-reproductive")
e1_size$group <- factor(e1_size$group, glev); e2_size$group <- factor(e2_size$group, glev)

## Table S6 — descriptive statistics
tableS6 <- bind_rows(
  e1_size %>% group_by(Experiment = "Experiment 1", Group = group) %>%
    summarise(n = n(),
              CL_mean = mean(CL), CL_sd = sd(CL), CL_min = min(CL), CL_max = max(CL),
              W_mean = mean(W, na.rm = TRUE), W_sd = sd(W, na.rm = TRUE), .groups = "drop"),
  e2_size %>% group_by(Experiment = "Experiment 2", Group = group) %>%
    summarise(n = n(),
              CL_mean = mean(CL), CL_sd = sd(CL), CL_min = min(CL), CL_max = max(CL),
              W_mean = mean(W, na.rm = TRUE), W_sd = sd(W, na.rm = TRUE), .groups = "drop")
) %>% mutate(across(where(is.numeric), ~round(., 3)))
cat("\n== Table S6: body size by reproductive group ==\n"); print(tableS6, n = Inf)

## Formal tests — do groups differ in size?
cat("\n-- Experiment 1 CL --\n")
print(summary(aov(CL ~ group, data = e1_size)))
print(kruskal.test(CL ~ group, data = e1_size))
print(TukeyHSD(aov(CL ~ group, data = e1_size)))
cat("\n-- Experiment 1 weight --\n"); print(summary(aov(W ~ group, data = e1_size)))
cat("\n-- Experiment 2 CL (group x shelter) --\n")
print(car::Anova(lm(CL ~ group * shelter, data = e2_size), type = 2))

## Confound test — does the reproductive-stage effect survive CL adjustment?
e1m <- e1_size %>% filter(!is.na(killed), killed >= 0, killed <= prey)
b1_base <- glmmTMB(cbind(killed, prey - killed) ~ group + scale(prey),
                   data = e1m, family = binomial)
b1_cov  <- glmmTMB(cbind(killed, prey - killed) ~ group + scale(prey) + scale(CL),
                   data = e1m, family = binomial)
cat("\n-- Exp 1 feeding model without / with CL covariate --\n")
print(broom.mixed::tidy(b1_base)); print(broom.mixed::tidy(b1_cov))
cat("Pearson r (CL, killed), Exp 1: ",
    round(cor(e1m$CL, e1m$killed, use = "complete.obs"), 3), "\n")

e2m <- e2_size %>% mutate(type = factor(group, glev),
                          shelter = factor(shelter, c("No shelter","Shelter")))
b2_base <- glmmTMB(cbind(killed, 12 - killed) ~ type * shelter,
                   data = e2m, family = binomial)
b2_cov  <- glmmTMB(cbind(killed, 12 - killed) ~ type * shelter + scale(CL),
                   data = e2m, family = binomial)
cat("\n-- Exp 2 feeding model without / with CL covariate --\n")
print(broom.mixed::tidy(b2_cov))

## ============================================================================
##  FIGURE 3 — partial consumption across densities (comment 12: whiskers >= 0)
## ============================================================================
grab <- function(sheet, rows) read_xlsx(f_raw, sheet = sheet)[rows, c(1:3,13,14:16,26:29,39)]
ex1 <- grab("Experiment 1", 1:72)[,-c(3,7,11)]
ex2 <- grab("Experiment 2", 1:54)[,-c(2,3,6,7,10,11)]; ex2$Prey <- 12
colnames(ex1) <- c("Code_N","Prey_N","Partial_N","Code_G","Prey_G","Partial_G","Code_E","Prey_E","Partial_E")
long1 <- ex1 %>% pivot_longer(everything(), names_to = c(".value","Group"), names_sep = "_") %>%
  transmute(Group, Prey, Partial) %>% mutate(Experiment = "Experiment 1")
summ <- long1 %>% group_by(Group, Prey) %>%
  summarise(mean_partial = mean(Partial, na.rm = TRUE),
            sd_partial   = sd(Partial, na.rm = TRUE), .groups = "drop")
p3 <- ggplot(summ, aes(Prey, mean_partial)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(ymin = pmax(0, mean_partial - sd_partial),   # <-- floored at zero
                      ymax =        mean_partial + sd_partial,
                      colour = Group),
                  position = position_dodge(width = 2)) +
  scale_colour_manual(values = c("E"="forestgreen","N"="deepskyblue3","G"="chocolate3")) +
  scale_x_continuous(breaks = c(1,3,6,12,18,24)) +
  ylab("Prey partially consumed") + xlab("Prey density") +
  theme_bw(base_size = 14) + theme(axis.title = element_text(face = "bold"))
ggsave(file.path(out_dir,"Figure3_partial_consumption.svg"), p3, width = 7, height = 5)

## ============================================================================
##  WRITE ALL TABLES + SESSION INFO
## ============================================================================
writexl::write_xlsx(
  list("Table1_FR"  = tableS1_FR,
       "TableS2_PC1"= as.data.frame(tableS2),
       "TableS4_E2" = as.data.frame(tableS4),
       "TableS5_PC2"= as.data.frame(tableS5),
       "TableS6_size" = as.data.frame(tableS6)),
  path = file.path(out_dir, "Reproduced_tables.xlsx"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("\nDONE. Tables + figures written to: ", normalizePath(out_dir), "\n")
