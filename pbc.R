# ────────────────────────────────────────────────────────────────
# Landmark-based dynamic prediction of mean residual lifetime (MRL)
# Using pseudo-observations and IPCW in the pbc2 dataset
# ────────────────────────────────────────────────────────────────

# Load required packages
library(JM)         # pbc2 dataset
library(pseudo)     # pseudo-observation calculation
library(dynpred)    # landmarking utilities
library(gee)        # GEE models
library(survival)   # needed for survfit in IPCW

# ────────────────────────────────────────────────────────────────
# Data preparation
# ────────────────────────────────────────────────────────────────
data(pbc2)

data <- data.frame(
  id          = pbc2$id,
  time        = pbc2$years,
  year        = pbc2$year,
  age         = pbc2$age,
  serBilir    = pbc2$serBilir,
  albumin     = pbc2$albumin,
  prothrombin = pbc2$prothrombin,
  histologic  = pbc2$histologic
)

data$delta      <- as.numeric(pbc2$status == "dead")
data$drug       <- as.numeric(pbc2$drug  == "D-penicil")
data$sex        <- as.numeric(pbc2$sex   == "female")
data$edema      <- as.numeric(pbc2$edema != "No edema")

# Truncation time
tau <- 12

data$Y         <- pmin(data$time, tau)
data$delta_L   <- ifelse(tau <= data$time, 1, 0)
data$delta_star <- data$delta + data$delta_L * (1 - data$delta)

# ────────────────────────────────────────────────────────────────
# Landmark settings
# ────────────────────────────────────────────────────────────────
sL  <- 6                # maximum follow-up time considered
LMs <- seq(0, sL, by = 0.2)
nsl <- length(LMs)

fixed   <- c("drug", "sex", "age")
varying <- c("serBilir", "edema", "albumin", "prothrombin", "histologic")

# ────────────────────────────────────────────────────────────────
# IPCW weight function
# ────────────────────────────────────────────────────────────────
ipcw_computation <- function(df, t) {
  present_ind  <- as.integer(df$Y > t)
  censored_ind <- 1 - df$delta_star
  
  kmf <- survfit(Surv(df$Y, censored_ind) ~ 1)
  kmf.summ <- summary(kmf, times = df$Y)
  kmf_df <- data.frame(days = kmf.summ[[2]], prob = kmf.summ[[6]])
  kmf_df <- kmf_df[match(df$Y, kmf_df$days), ]
  
  ipcw <- (df$delta_star * present_ind) / kmf_df$prob
  ipcw <- ifelse(df$delta_star == 0, 0, ipcw)
  return(ipcw)
}

# ────────────────────────────────────────────────────────────────
# Create landmark super dataset
# ────────────────────────────────────────────────────────────────
LMdata <- NULL
for (j in seq_along(LMs)) {
  LM <- cutLM(data,
              outcome = list(time = "Y", status = "delta"),
              LM = LMs[j], horizon = tau,
              covs = list(fixed = c(fixed, "delta_star"), varying = varying),
              format = "long", id = "id", rtime = "year", right = FALSE
  )
  
  LM$residual <- LM$Y - LMs[j]
  
  if (sum(LM$delta == 1) == 0) next
  
  LM$ipcw   <- ipcw_computation(df = LM, t = LMs[j])
  LM$pseudo <- pseudomean(time = LM$residual, event = LM$delta, tmax = tau - LMs[j])
  
  LMdata <- rbind(LMdata, LM)
}

LMdata <- LMdata[order(LMdata$id), ]
LMdata$LM_st <- (LMdata$LM - mean(LMs)) / sL

# ────────────────────────────────────────────────────────────────
# Create time interactions and LM polynomial terms
# ────────────────────────────────────────────────────────────────
define <- function(data, cov, time) {
  for (i in 1:length(cov)) {
    s <- paste("data$", cov[i], ".t", 0:2, "<-data$", cov[i],
               c("", "*time", "*time^2"), sep = "")
    eval(parse(text = s))
  }
  data$LM1 <- time
  data$LM2 <- time^2
  return(data)
}

LMdata <- define(data = LMdata, cov = c(fixed, varying), time = LMdata$LM_st)

# ────────────────────────────────────────────────────────────────
# Dynamic models (super models)
# ────────────────────────────────────────────────────────────────
dyn_po <- geeglm(
  pseudo ~ drug.t0 + sex.t0 + age.t0 +
    serBilir.t0 + serBilir.t1 +
    edema.t0 + albumin.t0 +
    prothrombin.t0 + prothrombin.t1 + prothrombin.t2 +
    histologic.t0 + histologic.t1 + LM1 + LM2,
  data = LMdata, id = id,
  scale.fix = FALSE, family = gaussian, corstr = "independence"
)

dyn_ipcw <- geeglm(
  residual ~ drug.t0 + sex.t0 + age.t0 +
    serBilir.t0 + serBilir.t1 +
    edema.t0 + albumin.t0 + albumin.t1 + albumin.t2 +
    prothrombin.t0 + prothrombin.t1 + prothrombin.t2 +
    histologic.t0 + LM1 + LM2,
  data = LMdata, id = id, weights = ipcw,
  scale.fix = FALSE, family = gaussian, corstr = "independence"
)

# Coefficient summary
round(cbind(
  mean = dyn_po$geese$beta,
  SD   = sqrt(diag(dyn_po$geese$vbeta)),
  Z    = dyn_po$geese$beta / sqrt(diag(dyn_po$geese$vbeta)),
  PVal = 2 - 2 * pnorm(abs(dyn_po$geese$beta / sqrt(diag(dyn_po$geese$vbeta))))
), 3)

round(cbind(
  mean = dyn_ipcw$geese$beta,
  SD   = sqrt(diag(dyn_ipcw$geese$vbeta)),
  Z    = dyn_ipcw$geese$beta / sqrt(diag(dyn_ipcw$geese$vbeta)),
  PVal = 2 - 2 * pnorm(abs(dyn_ipcw$geese$beta / sqrt(diag(dyn_ipcw$geese$vbeta))))
), 3)

# ────────────────────────────────────────────────────────────────
# Final formulas after variable selection
# ────────────────────────────────────────────────────────────────
fpseudo <- pseudo ~ drug.t0 + sex.t0 + age.t0 + serBilir.t0 + serBilir.t1 +
  edema.t0 + albumin.t0 +
  prothrombin.t0 + prothrombin.t2 +
  histologic.t0 + histologic.t1 + LM1 + LM2

fipcw <- residual ~ drug.t0 + sex.t0 + age.t0 + serBilir.t0 + edema.t0 +
  albumin.t0 + albumin.t1 + albumin.t2 +
  prothrombin.t0 + prothrombin.t1 + prothrombin.t2 +
  histologic.t0 + LM1 + LM2

# ────────────────────────────────────────────────────────────────
# Static dataset preparation (baseline)
# ────────────────────────────────────────────────────────────────
dex <- duplicated(data$id)
data_MRL <- data[!dex, ]

static <- NULL
for (j in seq(along = LMs)) {
  dt <- data_MRL[data_MRL$Y > LMs[j], ]
  if (sum(dt$delta == 1) == 0) next
  
  dt$LM       <- LMs[j]
  dt$residual <- dt$Y - LMs[j]
  dt$ipcw     <- ipcw_computation(df = dt, t = LMs[j])
  dt$pseudo   <- pseudomean(time = dt$residual, event = dt$delta, tmax = tau - LMs[j])
  
  static <- rbind(static, dt)
}

# ────────────────────────────────────────────────────────────────
# Static models at landmark = 0
# ────────────────────────────────────────────────────────────────
data_MRL_0 <- static[static$LM == 0, ]

ffpseudo <- pseudo ~ drug + sex + age + serBilir + edema + albumin + prothrombin + histologic
ffipcw   <- residual ~ drug + sex + age + serBilir + edema + albumin + prothrombin + histologic

sta_po <- geeglm(ffpseudo, data = data_MRL_0, id = id,
                 scale.fix = FALSE, family = gaussian, corstr = "independence")

sta_ipcw <- geeglm(ffipcw, data = data_MRL_0, id = id, weights = ipcw,
                   scale.fix = FALSE, family = gaussian, corstr = "independence")

# Coefficient summary for static models
round(cbind(
  mean = sta_po$geese$beta,
  SD   = sqrt(diag(sta_po$geese$vbeta)),
  Z    = sta_po$geese$beta / sqrt(diag(sta_po$geese$vbeta)),
  PVal = 2 - 2 * pnorm(abs(sta_po$geese$beta / sqrt(diag(sta_po$geese$vbeta))))
), 3)

round(cbind(
  mean = sta_ipcw$geese$beta,
  SD   = sqrt(diag(sta_ipcw$geese$vbeta)),
  Z    = sta_ipcw$geese$beta / sqrt(diag(sta_ipcw$geese$vbeta)),
  PVal = 2 - 2 * pnorm(abs(sta_ipcw$geese$beta / sqrt(diag(sta_ipcw$geese$vbeta))))
), 3)

# ────────────────────────────────────────────────────────────────
# Plot coefficient evolution over landmark times
# ────────────────────────────────────────────────────────────────
plot_coef <- function(x, start, ylimit, mintitle, pl = TRUE,
                      legendx, legendy, legendtxt) {
  bet <- matrix(x, length(LMs))
  stop <- start + ncol(bet) - 1
  
  LMsmooth <- data.frame(
    LMs   = LMs,
    dRMST = as.numeric(bet %*% dyn_po$geese$beta[start:stop])
  )
  
  se <- sqrt(diag(bet %*% dyn_po$geese$vbeta[start:stop, start:stop] %*% t(bet)))
  LMsmooth$lower <- LMsmooth$dRMST - qnorm(0.975) * se
  LMsmooth$upper <- LMsmooth$dRMST + qnorm(0.975) * se
  
  par(xaxs = "i", yaxs = "i", mar = c(2, 2.5, 2, 0.5))
  plot(LMsmooth$LMs, LMsmooth$dRMST, type = "l", lwd = 3,
       xlim = c(min(LMs), max(LMs)), ylim = ylimit,
       bty = "l", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  
  axis(2, las = 1, pos = min(LMs), cex.axis = 0.8, tcl = -0.2, hadj = 0.3, lwd = 1)
  axis(1, las = 1, pos = ylimit[1], cex.axis = 0.8, tcl = -0.2, padj = -2, lwd = 1)
  
  title(main = list(mintitle, cex = 1), line = 0.7)
  title(xlab = "Prediction time in year (s)", cex.lab = 0.8, line = 0.7)
  title(ylab = "Difference in MRL", cex.lab = 0.8, line = 1.3)
  
  lines(LMsmooth$LMs, LMsmooth$lower, lty = 2, col = "darkgray", lwd = 1)
  lines(LMsmooth$LMs, LMsmooth$upper, lty = 2, col = "darkgray", lwd = 1)
  
  if (pl) abline(h = 0, lwd = 1, lty = 3, col = "lightgray")
  legend(legendx, legendy, legendtxt, cex = 1, bty = "n", lwd = 1)
}

# Extract model coefficients and variance-covariance matrix
bet <- dyn_po$geese$beta
sig <- dyn_po$geese$vbeta

# Create basis functions for time-varying effects
a <- rep(1, length(LMs)) # Intercept term
b <- (LMs - mean(LMs)) / sL # Linear term (scaled)
c <- ((LMs - mean(LMs)) / sL)^2 # Quadratic term
ab <- c(a, b)
ac <- c(a, c)
abc <- c(a, b, c)

# Example usage
cbind(round(dyn_po$coefficients, 3), 1:length(dyn_po$coefficients))
plot_coef(ab, start = 12, ylimit = c(-2, 0.5), mintitle = "Histologic",
          legendx = 0, legendy = 0.5, legendtxt = "histologic stage")

# ────────────────────────────────────────────────────────────────
# Static prediction for a single patient
# ────────────────────────────────────────────────────────────────
sta_pred <- function(patientid) {
  tt <- seq(0, sL, by = 0.2)
  data_pre <- data[!duplicated(data$id), ]
  data_pre <- data_pre[data_pre$id == patientid, ]
  
  sta_MRL <- data.frame(array(NA, dim = c(length(tt), 3)))
  
  for (i in 1:length(tt)) {
    LMdata_i <- static[static$LM == LMs[i], ]
    sta_po <- gee(ffpseudo, data = LMdata_i, id = id,
                  scale.fix = FALSE, family = gaussian, corstr = "independence")
    
    a <- paste("data_pre$", sta_po$xnames[-1], sep = "", collapse = ",")
    b <- paste("matrix(c(rep(1,1),", a, "),1,length(sta_po$xnames))")
    mm <- eval(parse(text = b))
    
    sta_MRL[i, 1] <- as.numeric(sta_po$coef %*% t(mm))
    
    df <- sta_po$nobs - length(sta_po$coef)
    tfrac <- qt(0.025, df)
    ip <- mm %*% sta_po$robust.variance %*% t(mm)
    hwid <- tfrac * sqrt(ip)
    
    sta_MRL[i, 2:3] <- sta_MRL[i, 1] + c(hwid, -hwid)
  }
  return(cbind(tt, sta_MRL))
}

# ────────────────────────────────────────────────────────────────
# Dynamic prediction for a single patient
# ────────────────────────────────────────────────────────────────
dyn_pred <- function(patientid) {
  tt <- seq(0, sL, by = 0.2)
  data_pre <- LMdata[LMdata$id == patientid, ]
  inter <- findInterval(tt, data_pre$year)
  
  dyn_MRL <- data.frame(array(NA, dim = c(length(tt), 3)))
  data_p <- NULL
  
  for (i in 1:length(tt)) {
    data_pred <- data_pre[inter[i], ]
    data_pred$LM_st <- (tt[i] - mean(tt)) / sL
    data_p <- rbind(data_p, data_pred)
  }
  
  dyn_MRL[, 1] <- predict(dyn_po, data_p)
  
  df <- nrow(LMdata) - length(dyn_po$geese$beta)
  tfrac <- qt(0.025, df)
  
  data_pp <- data.frame(
    rep(1, nrow(data_p)),
    data_p$drug.t0, data_p$sex.t0, data_p$age.t0,
    data_p$serBilir.t0, data_p$serBilir.t1,
    data_p$edema.t0, data_p$albumin.t0,
    data_p$prothrombin.t0, data_p$prothrombin.t1, data_p$prothrombin.t2,
    data_p$histologic.t0, data_p$histologic.t1,
    data_p$LM1, data_p$LM2
  )
  
  data_pp <- as.matrix(data_pp)
  ip <- diag(data_pp %*% dyn_po$geese$vbeta %*% t(data_pp))
  hwid <- tfrac * sqrt(ip)
  dyn_MRL[, 2:3] <- dyn_MRL[, 1] + c(hwid, -hwid)
  
  return(cbind(tt, dyn_MRL))
}

# Example: Patient 113
num <- 113
sta_MRL <- sta_pred(patientid = num)
dyn_MRL <- dyn_pred(patientid = num)

# Plot comparison for one patient
plot(NA, type = "n", xlim = c(0, 3), ylim = c(0, 10),
     bty = "l", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
points(sta_MRL[,1], sta_MRL[,2], pch = 19, col = "#218365", cex = 1)
segments(sta_MRL[,1], sta_MRL[,4], sta_MRL[,1], sta_MRL[,3], col = "darkgray", lwd = 1.5)
lines(dyn_MRL[,1], dyn_MRL[,2], lwd = 2, col = "#f3990b")
lines(dyn_MRL[,1], dyn_MRL[,3], lwd = 1, lty = 2, col = "#f3990b")
lines(dyn_MRL[,1], dyn_MRL[,4], lwd = 1, lty = 2, col = "#f3990b")

axis(1, las = 1, pos = 0, cex.axis = 0.8, tcl = -0.2, padj = -2, lwd = 1)
axis(2, las = 1, pos = 0, cex.axis = 0.8, tcl = -0.2, hadj = 0.3, lwd = 1)

title(main = list("Patient A", cex = 1), line = 0.7)
title(xlab = "Prediction time (s)", cex.lab = 0.8, line = 0.7)
title(ylab = "MRL in years", cex.lab = 0.8, line = 1.3)

legend("topright", c("Static fixed model", "Dynamic supermodel"),
       col = c("#218365", "#f3990b"), pch = c(19, NA), lty = c(NA, 1),
       lwd = c(NA, 3), bty = "n", xpd = TRUE, cex = 0.75, pt.cex = 1)

# ────────────────────────────────────────────────────────────────
# Performance measures
# ────────────────────────────────────────────────────────────────
MRL_cindex <- function(model, data) {
  nt <- length(data$Y)
  ord <- order(data$Y, -data$delta)
  Y    <- data$Y[ord]
  delta <- data$delta[ord]
  dynpred.MRL <- predict(model, data)[ord]
  
  wh <- which(delta == 1)
  wh <- wh[wh <= nt - 1]
  
  total <- con <- 0
  for (m in wh) {
    for (n in (m+1):nt) {
      if (Y[n] > Y[m]) {
        total <- total + 1
        if (dynpred.MRL[n] > dynpred.MRL[m]) con <- con + 1
        if (dynpred.MRL[n] == dynpred.MRL[m]) con <- con + 0.5
      }
    }
  }
  return(con / total)
}

MRL_PE <- function(model, data, Y) {
  pre_MRL <- predict(model, data)
  numerator   <- sum(data$ipcw * abs(Y - pre_MRL))
  denominator <- sum(data$ipcw)
  return(numerator / denominator)
}

# ────────────────────────────────────────────────────────────────
# Monte-Carlo validation (dynamic vs static performance)
# ────────────────────────────────────────────────────────────────
u <- 200

dcii <- dcip <- dpei <- dpep <- matrix(NA, u, nsl)   # dynamic
scii <- scip <- spei <- spep <- matrix(NA, u, nsl)   # static

success_count <- 0
m <- 1
max_attempts <- 1000
attempt <- 0

while (success_count < u && attempt < max_attempts) {
  result <- tryCatch({
    set.seed(20250709 + 10 * m)
    index <- sample(2, nrow(data_MRL), replace = TRUE, prob = c(0.7, 0.3))
    trainid <- data_MRL$id[index == 1]
    testid  <- data_MRL$id[index == 2]
    
    dynamic_train <- LMdata[LMdata$id %in% trainid, ]
    dynamic_test  <- LMdata[LMdata$id %in% testid, ]
    static_train  <- static[static$id %in% trainid, ]
    static_test   <- static[static$id %in% testid, ]
    
    # Train dynamic models on bootstrap training set
    dynamic_train_data <- NULL
    for (j in seq_along(LMs)) {
      if (!any(dynamic_train$LM == LMs[j])) next
      dynamic_train_j <- dynamic_train[dynamic_train$LM == LMs[j], ]
      if (sum(dynamic_train_j$delta == 1) == 0) next
      
      dynamic_train_j$ipcw <- ipcw_computation(dynamic_train_j, LMs[j])
      dynamic_train_j$pseudo <- pseudomean(
        dynamic_train_j$residual, dynamic_train_j$delta, tau - LMs[j]
      )
      dynamic_train_data <- rbind(dynamic_train_data, dynamic_train_j)
    }
    
    dyn_ipcw_boot <- geeglm(fipcw, data = dynamic_train_data, id = id, weights = ipcw,
                            scale.fix = FALSE, family = gaussian, corstr = "independence")
    
    dyn_po_boot <- geeglm(fpseudo, data = dynamic_train_data, id = id,
                          scale.fix = FALSE, family = gaussian, corstr = "independence")
    
    if (is.null(dyn_ipcw_boot) || is.null(dyn_po_boot)) stop("model failed")
    
    # Evaluate on test set at each landmark
    for (k in seq_along(LMs)) {
      if (!any(dynamic_test$LM == LMs[k])) next
      dynamic_test_k <- dynamic_test[dynamic_test$LM == LMs[k], ]
      dynamic_test_k$ipcw <- ipcw_computation(dynamic_test_k, LMs[k])
      
      dcii[m, k] <- MRL_cindex(dyn_ipcw_boot, dynamic_test_k)
      dcip[m, k] <- MRL_cindex(dyn_po_boot,   dynamic_test_k)
      dpei[m, k] <- MRL_PE(dyn_ipcw_boot, dynamic_test_k, dynamic_test_k$Y - LMs[k])
      dpep[m, k] <- MRL_PE(dyn_po_boot,   dynamic_test_k, dynamic_test_k$Y - LMs[k])
    }
    
    # Static models per landmark
    for (j in seq_along(LMs)) {
      if (!any(static_train$LM == LMs[j])) next
      static_train_j <- static_train[static_train$LM == LMs[j], ]
      if (sum(static_train_j$delta == 1) == 0) next
      
      static_train_j$pseudo <- pseudomean(static_train_j$residual, static_train_j$delta, tau - LMs[j])
      static_train_j$ipcw   <- ipcw_computation(static_train_j, LMs[j])
      
      sta_ipcw_boot <- geeglm(ffipcw, data = static_train_j, id = id, weights = ipcw,
                              scale.fix = FALSE, family = gaussian, corstr = "independence")
      
      sta_po_boot <- geeglm(ffpseudo, data = static_train_j, id = id,
                            scale.fix = FALSE, family = gaussian, corstr = "independence")
      
      if (is.null(sta_ipcw_boot) || is.null(sta_po_boot)) next
      
      static_test_j <- static_test[static_test$LM == LMs[j], ]
      static_test_j$ipcw <- ipcw_computation(static_test_j, LMs[j])
      
      scii[m, j] <- MRL_cindex(sta_ipcw_boot, static_test_j)
      scip[m, j] <- MRL_cindex(sta_po_boot,   static_test_j)
      spei[m, j] <- MRL_PE(sta_ipcw_boot, static_test_j, static_test_j$Y - LMs[j])
      spep[m, j] <- MRL_PE(sta_po_boot,   static_test_j, static_test_j$Y - LMs[j])
    }
    
    success_count <- success_count + 1
    m <- m + 1
    attempt <- 0
    TRUE
  }, error = function(e) FALSE)
  
  if (!result) attempt <- attempt + 1
}

# Average performance over bootstrap replicates
dcii <- apply(dcii, 2, mean, na.rm = TRUE)
dcip <- apply(dcip, 2, mean, na.rm = TRUE)
dpei <- apply(dpei, 2, mean, na.rm = TRUE)
dpep <- apply(dpep, 2, mean, na.rm = TRUE)

scii <- apply(scii, 2, mean, na.rm = TRUE)
scip <- apply(scip, 2, mean, na.rm = TRUE)
spei <- apply(spei, 2, mean, na.rm = TRUE)
spep <- apply(spep, 2, mean, na.rm = TRUE)

CI <- cbind(LMs, dcii, dcip, dpei, dpep, scii, scip, spei, spep)

# ────────────────────────────────────────────────────────────────
# Final performance plots
# ────────────────────────────────────────────────────────────────
# C-index plot
plot(CI[,1], CI[,2], lwd = 2, lty = 3, col = "#54acce", type = "l", bty = "l",
     xlim = c(0, sL), ylim = c(0.45, 0.9), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
axis(1, las = 1, pos = 0.45, cex.axis = 0.8, tcl = -0.2, padj = -2, lwd = 1)
axis(2, las = 1, pos = 0, cex.axis = 0.8, tcl = -0.2, hadj = 0.3, lwd = 1)
title(xlab = "Prediction time (s)", cex.lab = 0.8, line = 0.7)
title(ylab = "C-index", cex.lab = 0.8, line = 1.3)

lines(CI[,1], CI[,3], lwd = 2, lty = 4, col = "#f3990b")
lines(CI[,1], CI[,6], lwd = 2, lty = 1, col = "#ffe93f")
lines(CI[,1], CI[,7], lwd = 2, lty = 2, col = "#218365")

legend("bottomleft", c("Dynamic IPCW", "Dynamic PO", "Static IPCW", "Static PO"),
       col = c("#54acce", "#f3990b", "#ffe93f", "#218365"),
       lty = c(3,4,1,2), bty = "n", lwd = 2, xpd = TRUE, cex = 0.75)

# Prediction error plot
plot(CI[,1], CI[,4], lwd = 2, lty = 3, col = "#54acce", type = "l", bty = "l",
     xlim = c(0, sL), ylim = c(0, 3), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
axis(1, las = 1, pos = 0, cex.axis = 0.8, tcl = -0.2, padj = -2, lwd = 1)
axis(2, las = 1, pos = 0, cex.axis = 0.8, tcl = -0.2, hadj = 0.3, lwd = 1)
title(xlab = "Prediction time (s)", cex.lab = 0.8, line = 0.7)
title(ylab = "Prediction error", cex.lab = 0.8, line = 1.3)

lines(CI[,1], CI[,5], lwd = 2, lty = 4, col = "#f3990b")
lines(CI[,1], CI[,8], lwd = 2, lty = 1, col = "#ffe93f")
lines(CI[,1], CI[,9], lwd = 2, lty = 2, col = "#218365")

legend("bottomleft", c("Dynamic IPCW", "Dynamic PO", "Static IPCW", "Static PO"),
       col = c("#54acce", "#f3990b", "#ffe93f", "#218365"),
       lty = c(3,4,1,2), bty = "n", lwd = 2, xpd = TRUE, cex = 0.75)