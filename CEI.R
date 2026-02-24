# =====================================================
# Bản essay1 chốt 
# EMH pipeline MASTER (FX + STOCK) -- v1 + SCALE TESTS
# - Input: daily OHLC (FX hoặc stock index)
# - Output:
#   + df_month: 10/12 sub-indices EMH theo (id, ym)
#   + df_emh  : composite EMH (pc, ave) theo (id, ym)
#   + In-console: Kiểm định thang đo (alpha, omega, convergent)
#
# Tham số:
#   DATA_MODE       = "FX" | "STOCK"
#   USE_LOCAL_HURST = TRUE  -> Hurst local 252d
#   USE_CLV_GAP     = FALSE -> 10 chỉ số; TRUE -> 12 chỉ số (thêm CLV & GAP)
# =====================================================

# ---- 0) SETUP & PARAMS ----
DATA_MODE       <- "FX"   # "FX" hoặc "STOCK"
USE_LOCAL_HURST <- TRUE
USE_CLV_GAP     <- FALSE    # FALSE => 10 chỉ số; TRUE => 12 chỉ số

need <- c("readxl","dplyr","tidyr","lubridate","zoo","e1071","psych")
to_install <- setdiff(need, rownames(installed.packages()))
if(length(to_install)) install.packages(to_install, dependencies = TRUE)

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr)
  library(lubridate); library(zoo); library(e1071)
  library(psych)
})

# ---- 1) LOAD DATA (FX vs STOCK) ----
if (DATA_MODE == "FX") {
  cat(">>> Đọc dữ liệu FX...\n")
  excel_path_fx <- "C:/Users/trong/OneDrive/HUYNHTT/Thesis 2025/Paper/Measure EMH/FX ASEAN.xlsx"
  df <- read_excel(excel_path_fx, sheet = "full_exchange") %>%
    rename_with(~ trimws(.x)) %>% rename_with(tolower)
  
} else if (DATA_MODE == "STOCK") {
  cat(">>> Đọc dữ liệu stock index...\n")
  excel_path_stock <- "C:/Users/trong/OneDrive/HUYNHTT/Thesis 2025/Data/stock asean/index_asean.xlsx"
  df <- read_excel(excel_path_stock) %>%
    rename_with(~ trimws(.x)) %>% rename_with(tolower)
} else {
  stop("DATA_MODE phải là 'FX' hoặc 'STOCK'")
}

# Chuẩn hoá tên cột
if ("date" %in% names(df) && !("dates" %in% names(df))) df <- df %>% rename(dates = date)
if (!("id" %in% names(df))) {
  cand <- intersect(c("country","market","ticker","symbol","code"), names(df))
  if (length(cand) >= 1) df <- df %>% rename(id = all_of(cand[1])) else stop("Không tìm được cột id/market.")
}
stopifnot(all(c("close","open","high","low","dates","id") %in% names(df)))

# Ép kiểu & tạo ym
if (!inherits(df$dates, "Date")) {
  df$dates <- if (is.numeric(df$dates)) as.Date(df$dates, origin="1899-12-30") else as.Date(df$dates)
}
if (!("ym" %in% names(df))) {
  df <- df %>% mutate(ym = year(dates)*100 + month(dates))
}

df <- df %>%
  mutate(
    id = as.character(id),
    ym = as.integer(ym)
  ) %>%
  arrange(id, dates)

cat(">>> Số thị trường (id):", dplyr::n_distinct(df$id), "\n")

# ---- 2) HELPERS ----
safe_div <- function(a,b) ifelse(is.finite(a)&is.finite(b)&b!=0, a/b, NA_real_)

ac_abs_sum <- function(x, p=5){
  x <- x[is.finite(x)]
  if(length(x) < p+2) return(NA_real_)
  as <- stats::acf(x, plot=FALSE, lag.max=p, na.action=na.pass)$acf[-1]
  sum(abs(as), na.rm=TRUE)
}

lb_inEff <- function(x, lag=10){
  x <- x[is.finite(x)]
  if(length(x) < lag+2) return(NA_real_)
  p <- try(stats::Box.test(x, type="Ljung-Box", lag=lag)$p.value, silent=TRUE)
  if(inherits(p,"try-error")) return(NA_real_)
  1 - as.numeric(p)
}

variance_ratio <- function(r, k){
  r <- as.numeric(r)
  n <- length(r)
  if(n < k+2 || sum(is.finite(r)) < k+2) return(NA_real_)
  mu   <- mean(r, na.rm=TRUE)
  r_dm <- r - mu
  var1 <- var(r_dm, na.rm=TRUE)
  if(!is.finite(var1) || var1 <= 0) return(NA_real_)
  rk <- zoo::rollapply(r, width=k, FUN=sum, align="right", partial=FALSE)
  rk <- as.numeric(rk) - k*mu
  vark <- var(rk, na.rm=TRUE)
  vark / (k * var1)
}

tail_freq <- function(r, k=2){
  r <- r[is.finite(r)]
  if(length(r) < 10) return(NA_real_)
  s <- sd(r, na.rm=TRUE)
  if(!is.finite(s) || s==0) return(NA_real_)
  mean(abs(r) > k*s, na.rm=TRUE)
}

gk_var <- function(O,H,L,C){
  u <- log(H/O); d <- log(L/O); cc <- log(C/O)
  0.511*(u-d)^2 - 0.019*(cc*(u+d)) - 0.383*(cc^2)
}

cs_spread_daily <- function(H, L){
  n <- length(H); if(n < 2) return(rep(NA_real_, n))
  H <- as.numeric(H); L <- as.numeric(L)
  H[H <= 0] <- NA_real_; L[L <= 0] <- NA_real_
  hl <- log(H/L); hl[!is.finite(hl)] <- NA_real_
  hl1 <- dplyr::lag(hl)
  beta <- (hl^2 + hl1^2)/2
  Hlag <- dplyr::lag(H); Llag <- dplyr::lag(L)
  ratio <- ifelse(is.finite(H)&is.finite(L)&is.finite(Hlag)&is.finite(Llag)&
                    pmax(H,Hlag) > 0 & pmin(L,Llag) > 0,
                  pmax(H,Hlag)/pmin(L,Llag), NA_real_)
  gamma <- log(ratio)^2
  k <- 3 - 2*sqrt(2)
  alpha <- suppressWarnings((sqrt(2*beta) - sqrt(beta))/k - sqrt(gamma/k))
  alpha[!is.finite(alpha)] <- NA_real_
  S <- 2*(exp(alpha)-1)/(1+exp(alpha))
  S[!is.finite(S)] <- NA_real_
  pmax(pmin(S,1),0)
}

safe_mean <- function(x){
  x <- x[is.finite(x)]
  if(length(x)==0) NA_real_ else mean(x)
}

clip3 <- function(z) pmin(pmax(z, -3), 3)

zscore_full <- function(x){
  m <- mean(x, na.rm=TRUE); s <- sd(x, na.rm=TRUE)
  if(!is.finite(s) || s==0) return(rep(NA_real_, length(x)))
  (x - m)/s
}

# ---- 3) DAILY TRANSFORMS ----
df_daily <- df %>%
  group_by(id) %>%
  mutate(
    r_cc   = log(close) - dplyr::lag(log(close)),
    r_on   = log(open)  - dplyr::lag(log(close)),
    r_id   = log(close) - log(open),
    rel_rg_pc = safe_div(high - low, dplyr::lag(close)),
    CLV    = safe_div(close - open, pmax(high - low, 1e-12)),
    gap    = safe_div(abs(open - dplyr::lag(close)), pmax(high - low, 1e-12)),
    gkv    = gk_var(open, high, low, close),
    cs     = cs_spread_daily(high, low)
  ) %>%
  ungroup()

# ---- 4) LOCAL HURST (tùy chọn) ----
get_H_one <- function(x){
  x <- x[is.finite(x)]
  if(length(x) < 60) return(NA_real_)
  if(requireNamespace("pracma", quietly = TRUE)){
    out <- try(pracma::hurstexp(x, display=FALSE), silent=TRUE)
    if(!inherits(out,"try-error")){
      if(!is.null(out$Hs) && is.finite(out$Hs)) return(out$Hs)
      if(!is.null(out$Ha) && is.finite(out$Ha)) return(out$Ha)
    }
  }
  # fallback: ACF(1) heuristic
  r1 <- suppressWarnings(acf(x, plot=FALSE, lag.max=1, na.action=na.pass)$acf[2])
  if(!is.finite(r1)) return(NA_real_)
  H <- 0.5 + 0.5*max(min(r1,1),-1)
  max(min(H,1),0)
}

if (USE_LOCAL_HURST) {
  cat(">>> Tính Hurst LOCAL 252d...\n")
  win_days <- 252
  min_days <- 120
  
  calc_local_H_daily <- function(dd){
    dd <- dd %>% arrange(dates)
    r  <- dd$r_cc; n <- length(r)
    Ht <- rep(NA_real_, n)
    if(n < min_days) return(dd %>% mutate(H = NA_real_, inef_hurst_d = NA_real_))
    for(i in seq_len(n)){
      i2 <- i; i1 <- max(1, i2 - win_days + 1)
      w  <- r[i1:i2]
      if(sum(is.finite(w)) >= min_days) Ht[i] <- get_H_one(w) else Ht[i] <- NA_real_
    }
    dd %>% mutate(H = Ht, inef_hurst_d = ifelse(is.finite(H), abs(H - 0.5), NA_real_))
  }
  
  df_H_daily <- df_daily %>%
    group_by(id) %>%
    group_modify(~ calc_local_H_daily(.x)) %>%
    ungroup()
  
  df_hm <- df_H_daily %>%
    mutate(ym = year(dates)*100 + month(dates)) %>%
    group_by(id, ym) %>%
    summarise(
      inef_hurst_m = mean(inef_hurst_d, na.rm=TRUE),
      nH_days      = sum(is.finite(inef_hurst_d)),
      .groups="drop"
    )
} else {
  df_hm <- NULL
}

# ---- 5) MONTHLY SUB-INDICES (12m window) ----

df_month <- df %>% distinct(id, ym) %>% arrange(id, ym)

# bộ 10 chỉ số core
sub_names_core <- c(
  "inef_autocorr","inef_vr","inef_vc","inef_rvp",
  "inef_hurst","inef_skew","inef_kurt","inef_rdi",
  "inef_tail","inef_rcvm"
)

# nếu muốn thêm CLV & GAP (12 chỉ số), bật USE_CLV_GAP = TRUE
if (USE_CLV_GAP) {
  sub_names <- c(sub_names_core, "inef_clv","inef_gap")
} else {
  sub_names <- sub_names_core
}

for(nm in sub_names) df_month[[nm]] <- NA_real_

calc_subindices_window <- function(d){
  r <- d$r_cc
  res <- list()
  
  # 1) Autocorr inefficiency
  res$inef_autocorr <- mean(c(ac_abs_sum(r,5), lb_inEff(r,10)), na.rm=TRUE)
  
  # 2) Variance Ratio
  vrs <- sapply(c(2,5,10), function(k) variance_ratio(r,k))
  res$inef_vr <- if(all(!is.finite(vrs))) NA_real_ else sum(abs(vrs-1), na.rm=TRUE)
  
  # 3) Vol clustering
  res$inef_vc <- mean(c(ac_abs_sum(abs(r),5), ac_abs_sum(r^2,5)), na.rm=TRUE)
  
  # 4) Range-vol predictability
  sig_rg <- sqrt(pmax(0, d$gkv))
  res$inef_rvp <- ac_abs_sum(sig_rg,5)
  
  # 5) Hurst (placeholder, sẽ ghi đè bởi local H nếu có)
  res$inef_hurst <- NA_real_
  
  # 6) Skew
  res$inef_skew <- if(sum(is.finite(r))>=10) abs(e1071::skewness(r, na.rm=TRUE, type=2)) else NA_real_
  
  # 7) Kurtosis excess
  if(sum(is.finite(r))>=10){
    kx <- e1071::kurtosis(r, na.rm=TRUE, type=2)
    res$inef_kurt <- max(kx - 3, 0)
  } else res$inef_kurt <- NA_real_
  
  # 8) Return decomposition (ON vs ID)
  mdiff <- abs(safe_mean(d$r_on) - safe_mean(d$r_id))
  vrat  <- safe_div(var(d$r_on, na.rm=TRUE), var(d$r_id, na.rm=TRUE))
  res$inef_rdi <- mdiff + ifelse(is.finite(vrat)&&vrat>0, abs(log(vrat)), NA_real_)
  
  # 9) Tail: range + exceedance
  res$inef_tail <- safe_mean(d$rel_rg_pc) + tail_freq(r,2)
  
  # 10) RCVM
  s_rg <- sqrt(mean(pmax(0, d$gkv), na.rm=TRUE))
  s_cc <- sd(r, na.rm=TRUE)
  res$inef_rcvm <- if(is.finite(s_rg)&&is.finite(s_cc)&&s_rg>0&&s_cc>0) abs(log(s_rg)-log(s_cc)) else NA_real_
  
  # 11,12) CLV & GAP (optional)
  if (USE_CLV_GAP) {
    res$inef_clv <- mean(c(abs(safe_mean(d$CLV)), ac_abs_sum(d$CLV,5)), na.rm=TRUE)
    res$inef_gap <- mean(c(safe_mean(d$gap),      ac_abs_sum(d$gap,5)), na.rm=TRUE)
  }
  
  res
}

cat(">>> Tính sub-indices 12m rolling window...\n")
out_list <- vector("list", nrow(df_month))
for(i in seq_len(nrow(df_month))){
  this_id <- df_month$id[i]
  this_ym <- df_month$ym[i]
  lb      <- this_ym - 100L   # ~12m (phía trước)
  
  d_sub <- df_daily %>%
    filter(id == this_id, ym > lb, ym <= this_ym)
  
  if (nrow(d_sub) < 30) {
    out_list[[i]] <- c(list(id=this_id, ym=this_ym),
                       setNames(as.list(rep(NA_real_, length(sub_names))), sub_names))
  } else {
    out_list[[i]] <- c(list(id=this_id, ym=this_ym), calc_subindices_window(d_sub))
  }
}

df_month <- bind_rows(out_list)

# GHI ĐÈ inef_hurst = local H nếu dùng
if (USE_LOCAL_HURST && !is.null(df_hm)) {
  df_month <- df_month %>%
    left_join(df_hm, by = c("id","ym")) %>%
    mutate(inef_hurst = ifelse(is.finite(inef_hurst_m), inef_hurst_m, inef_hurst)) %>%
    select(-inef_hurst_m, -nH_days)
}

# ---- 6) COMPOSITE EMH = PCA GLOBAL (panel pooled) ----
cat(">>> PCA pooled (panel) để tạo EMH composite...\n")

df_month_clean <- df_month %>% arrange(id, ym)

sub_cols <- sub_names   # 10 hoặc 12 tuỳ USE_CLV_GAP

# z-score toàn panel + clip 3-sigma
dfZ <- df_month_clean %>%
  mutate(across(all_of(sub_cols), zscore_full, .names="{.col}_z")) %>%
  mutate(across(ends_with("_z"), clip3))

sub_cols_z <- paste0(sub_cols, "_z")
Z <- as.matrix(dfZ[, sub_cols_z])

# bỏ cột gần hằng số
keep_cols <- which(apply(Z, 2, function(col) sd(col, na.rm=TRUE) > 1e-8))
Zk <- Z[, keep_cols, drop=FALSE]

# fit PCA trên các hàng complete
keep_rows <- stats::complete.cases(Zk)
pc <- prcomp(Zk[keep_rows, ], center=FALSE, scale.=FALSE)

# định hướng PC1 cùng chiều với mean(Z)
sc_tr <- pc$x[,1]
mbar  <- rowMeans(Zk[keep_rows, ], na.rm=TRUE)
r     <- suppressWarnings(cor(sc_tr, mbar, use="pairwise.complete.obs"))
sgn   <- ifelse(is.finite(r) && r < 0, -1, 1)
v1    <- pc$rotation[,1] * sgn

# score toàn bộ
pc_all <- rep(NA_real_, nrow(dfZ))
rows_ok <- which(complete.cases(Zk))
pc_all[rows_ok] <- as.numeric(Zk[rows_ok, , drop=FALSE] %*% v1)

# ave_z: trung bình 10–12 z-score pooled
df_emh <- dfZ %>%
  mutate(
    pc  = pc_all,
    ave = ifelse(complete.cases(dplyr::select(., all_of(sub_cols_z))),
                 rowMeans(dplyr::select(., all_of(sub_cols_z)), na.rm=FALSE),
                 NA_real_)
  ) %>%
  select(id, ym, pc, ave)

cat(">>> Done EMH composite. Số quan sát có pc:", sum(is.finite(df_emh$pc)), "\n")

# =====================================================
# 8) SCALE TESTS: RELIABILITY + CONVERGENT VALIDITY
# =====================================================

cat("\n===== [SCALE TESTS] =====\n")

# ---- 8.1 RELIABILITY (POOLED PANEL) ----
cat(">>> [1] Reliability pooled (Cronbach alpha, Omega total)\n")

X_all <- dfZ %>%
  select(all_of(sub_cols_z)) %>%
  as.matrix()

X_cc <- X_all[stats::complete.cases(X_all), , drop = FALSE]

if (nrow(X_cc) < 30 || ncol(X_cc) < 2) {
  cat("Không đủ complete cases để ước lượng reliability (pooled).\n")
} else {
  alpha_all <- psych::alpha(X_cc, check.keys = TRUE)
  alpha_val <- as.numeric(alpha_all$total$raw_alpha)
  
  omega_val <- NA_real_
  omega_obj <- try(psych::omega(as.data.frame(X_cc),
                                nfactors = 1, plot = FALSE),
                   silent = TRUE)
  if (!inherits(omega_obj,"try-error") && !is.null(omega_obj$omega.tot)) {
    omega_val <- as.numeric(omega_obj$omega.tot)
  }
  
  cat(sprintf("Cronbach alpha (pooled, check.keys=TRUE): %.3f\n", alpha_val))
  if (is.finite(omega_val)) {
    cat(sprintf("Omega total (pooled):                  %.3f\n", omega_val))
  } else {
    cat("Omega total (pooled):                      NA (omega failed)\n")
  }
  
  cat("\nItem keying theo psych::alpha (1 = cùng chiều, -1 = đảo chiều):\n")
  print(alpha_all$keys)
}

# ---- 8.2 RELIABILITY BY ID (Cronbach alpha) ----
cat("\n>>> [1b] Reliability theo từng id (Cronbach alpha, check.keys=TRUE)\n")

rel_by_id <- dfZ %>%
  arrange(id, ym) %>%
  group_by(id) %>%
  group_modify(~{
    Xi <- .x %>% select(all_of(sub_cols_z)) %>% as.matrix()
    Xi <- Xi[stats::complete.cases(Xi), , drop = FALSE]
    if (nrow(Xi) < 20 || ncol(Xi) < 2) {
      tibble(alpha_id = NA_real_, n_cc = nrow(Xi))
    } else {
      ai <- psych::alpha(Xi, check.keys = TRUE)
      tibble(alpha_id = as.numeric(ai$total$raw_alpha),
             n_cc     = nrow(Xi))
    }
  }) %>%
  ungroup()

print(rel_by_id)

# ---- 8.3 CONVERGENT VALIDITY ----
cat("\n>>> [2] Convergent validity: pc vs ave, vs từng subindex\n")

# (a) Corr(pc, ave)
conv_pc_ave <- df_emh %>%
  filter(is.finite(pc), is.finite(ave)) %>%
  summarise(cor_pc_ave = cor(pc, ave, use = "pairwise.complete.obs")) %>%
  pull(cor_pc_ave)

cat(sprintf("Corr(pc, ave) toàn panel: %.3f\n", conv_pc_ave))

# (b) Corr(pc) / Corr(ave) với từng sub-index z-pooled
df_for_corr <- df_emh %>%
  inner_join(dfZ %>% select(id, ym, all_of(sub_cols_z)),
             by = c("id","ym"))

corr_pc_sub <- tibble(
  sub_index = sub_cols,
  corr      = sapply(sub_cols_z, function(v) {
    cor(df_for_corr$pc, df_for_corr[[v]], use = "pairwise.complete.obs")
  })
) %>%
  arrange(desc(abs(corr)))

corr_pc_sub_print <- corr_pc_sub %>%
  mutate(corr = round(corr, 3))

cat("\nCorr(pc) với từng sub-index (z-score pooled):\n")
print(corr_pc_sub_print)

corr_ave_sub <- tibble(
  sub_index = sub_cols,
  corr      = sapply(sub_cols_z, function(v) {
    cor(df_for_corr$ave, df_for_corr[[v]], use = "pairwise.complete.obs")
  })
) %>%
  arrange(desc(abs(corr)))

corr_ave_sub_print <- corr_ave_sub %>%
  mutate(corr = round(corr, 3))

cat("\nCorr(ave) với từng sub-index (z-score pooled):\n")
print(corr_ave_sub_print)

cat("\n===== [SCALE TESTS DONE] =====\n")

# ---- 9) LƯU OUTPUT ----
out_dir <- "C:/Users/trong/OneDrive/HUYNHTT/Year 2025/EMH"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

tag_asset <- ifelse(DATA_MODE=="FX","fx","stock")
tag_idx   <- ifelse(USE_CLV_GAP,"12idx","10idx")
tag_hurst <- ifelse(USE_LOCAL_HURST,"Hlocal","Hnone")

file_month <- file.path(out_dir, paste0("emh_subindices_", tag_asset, "_", tag_idx, "_", tag_hurst, ".csv"))
file_emh   <- file.path(out_dir, paste0("emh_composite_", tag_asset, "_", tag_idx, "_", tag_hurst, ".csv"))

write.csv(df_month, file_month, row.names=FALSE)
write.csv(df_emh,   file_emh,   row.names=FALSE)

cat("Saved sub-indices to: ", file_month, "\n")
cat("Saved composite EMH to:", file_emh, "\n")










# ================================
# A) PREPARE DATA FOR EFA
# ================================
X_efa <- dfZ %>%
  select(all_of(sub_cols_z)) %>%
  as.matrix()

# complete cases (giống logic alpha pooled)
X_efa_cc <- X_efa[complete.cases(X_efa), ]

cat("N obs (EFA):", nrow(X_efa_cc),
    " | N items:", ncol(X_efa_cc), "\n")

# ================================
# B) KMO & BARTLETT
# ================================
KMO_res <- psych::KMO(X_efa_cc)
Bart_res <- cortest.bartlett(X_efa_cc)

cat("\n--- KMO ---\n")
print(KMO_res$MSA)      # overall + từng biến

cat("\n--- Bartlett test ---\n")
print(Bart_res)

# ================================
# C) PARALLEL ANALYSIS
# ================================

#fa.parallel(  X_efa_cc,   fa       = "fa",   fm       = "ml",  n.iter   = 100,   error.bars = FALSE,  show.legend = TRUE,   main = "Parallel Analysis (ML)")


# ================================
# D) EFA: 5 FACTORS (ML + OBLIMIN)
# ================================
efa_5f <- psych::fa(
  X_efa_cc,
  nfactors = 5,
  rotate   = "oblimin",
  fm       = "ml"
)

print(efa_5f$loadings, cutoff = 0.30)

# ================================
# E) RELIABILITY BY FACTOR BLOCK
# ================================

blocks <- list(
  Moment   = c("inef_skew_z","inef_kurt_z","inef_tail_z"),
  Depend1  = c("inef_vr_z","inef_hurst_z","inef_autocorr_z"),
  Depend2  = c("inef_vc_z","inef_rvp_z")
)

for (b in names(blocks)) {
  Xb <- X_efa_cc[, blocks[[b]]]
  cat("\n---", b, "---\n")
  
  a <- psych::alpha(Xb, check.keys = TRUE)
  o <- psych::omega(as.data.frame(Xb), nfactors = 1, plot = FALSE)
  
  cat("Alpha :", round(a$total$raw_alpha, 3), "\n")
  cat("Omega :", round(o$omega.tot, 3), "\n")
}


# ================================
# F) FACTOR SCORES (THOMSON)
# ================================
fs <- psych::factor.scores(
  X_efa_cc,
  efa_5f,
  method = "Thurstone"
)$scores

fs <- as.data.frame(scale(fs))  # z-score
# ================================
# G) SIGN ORIENTATION
# ================================
for (j in seq_len(ncol(fs))) {
  ref <- rowMeans(X_efa_cc[, abs(efa_5f$loadings[,j]) > 0.3, drop=FALSE])
  sgn <- sign(cor(fs[,j], ref, use="pairwise.complete.obs"))
  fs[,j] <- fs[,j] * sgn
}

# ================================
# H) SECOND-STAGE PCA (CEI)
# ================================
pca2 <- prcomp(fs, center = FALSE, scale. = FALSE)

summary(pca2)   # xem EVR
cei_efa <- pca2$x[,1]

cat("Explained variance PC1:",
    round(summary(pca2)$importance[2,1]*100, 1), "%\n")




########## thêm biến thời gian
df_emh <- df_emh %>%
  mutate(
    y = ym %/% 100,
    m = ym %% 100
  )

y0 <- min(df_emh$y)
m0 <- df_emh$m[df_emh$y == y0][which.min(df_emh$ym)]

df_emh <- df_emh %>%
  mutate(
    t_seq = 12 * (y - y0) + (m - m0)
  ) %>%
  select(-y, -m)


##### Kiểm định AMH
df <- subset(df_emh,ym>200900 & ym < 202508)
library(dplyr)

desc_pc_id <- df %>%
  group_by(id) %>%
  summarise(
    obs    = sum(!is.na(pc)),
    min    = min(pc, na.rm = TRUE),
    P25    = quantile(pc, 0.25, na.rm = TRUE),
    median = median(pc, na.rm = TRUE),
    P75    = quantile(pc, 0.75, na.rm = TRUE),
    max    = max(pc, na.rm = TRUE),
    mean   = mean(pc, na.rm = TRUE),
    std    = sd(pc, na.rm = TRUE),
    .groups = "drop"
  )

# Save to CSV
write.csv(desc_pc_id, "desc_pc_by_id.csv", row.names = FALSE)

library(dplyr)

desc_pc_id_period <- df %>%
  mutate(
    period = ifelse(ym < 201400, "Pre-2014", "Post-2014")
  ) %>%
  group_by(id, period) %>%
  summarise(
    obs    = sum(!is.na(pc)),
    min    = min(pc, na.rm = TRUE),
    P25    = quantile(pc, 0.25, na.rm = TRUE),
    median = median(pc, na.rm = TRUE),
    P75    = quantile(pc, 0.75, na.rm = TRUE),
    max    = max(pc, na.rm = TRUE),
    mean   = mean(pc, na.rm = TRUE),
    std    = sd(pc, na.rm = TRUE),
    .groups = "drop"
  )

# Save CSV
write.csv(desc_pc_id_period,
          "desc_pc_by_id_pre_post_2014.csv",
          row.names = FALSE)
library(car)

leveneTest(pc ~ factor(id), data = df)

welch_res <- oneway.test(pc ~ factor(id),
                         data = df,
                         var.equal = FALSE)




welch_res

kruskal_res <- kruskal.test(pc ~ factor(id), data = df)
kruskal_res

library(FSA)

dunnTest(pc ~ factor(id),
         data = df,
         method = "bonferroni")




# =====================================================
# DESCRIPTIVE STATS + VISUALIZATION (Data in Brief style)
# - Input: df_month (sub-indices), df_emh (composite pc)
# - Output:
#   + Tables (csv): CEI summary by id, sub-index summary, corr matrix
#   + Figures (png): CEI time series, CEI distribution by id, average CEI by id
# =====================================================

need2 <- c("dplyr","tidyr","lubridate","ggplot2","scales")
to_install2 <- setdiff(need2, rownames(installed.packages()))
if(length(to_install2)) install.packages(to_install2, dependencies = TRUE)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(lubridate)
  library(ggplot2); library(scales)
})

# ---------- 0) Helpers ----------
ym_to_date <- function(ym){
  # ym = YYYYMM -> Date (first day of month)
  y <- ym %/% 100
  m <- ym %% 100
  as.Date(sprintf("%04d-%02d-01", y, m))
}

# If you want to restrict sample (optional)
# df_emh <- df_emh %>% filter(ym >= 200001, ym <= 202508)

# Ensure time variable
df_emh2 <- df_emh %>%
  mutate(date = ym_to_date(ym)) %>%
  arrange(id, date)

# Sub-index columns (10 core)
sub_cols_10 <- c(
  "inef_autocorr","inef_vr","inef_vc","inef_rvp",
  "inef_hurst","inef_skew","inef_kurt","inef_rdi",
  "inef_tail","inef_rcvm"
)

# If your df_month may contain extra cols (CLV/GAP), you can keep only the 10 core:
df_month2 <- df_month %>%
  select(any_of(c("id","ym", sub_cols_10))) %>%
  mutate(date = ym_to_date(ym)) %>%
  arrange(id, date)

# ---------- 1) Table: Descriptive stats of CEI by market ----------
tab_cei_by_id <- df_emh2 %>%
  group_by(id) %>%
  summarise(
    n   = sum(is.finite(pc)),
    min = min(pc, na.rm = TRUE),
    p25 = quantile(pc, 0.25, na.rm = TRUE),
    med = median(pc, na.rm = TRUE),
    p75 = quantile(pc, 0.75, na.rm = TRUE),
    max = max(pc, na.rm = TRUE),
    mean = mean(pc, na.rm = TRUE),
    sd   = sd(pc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean)

print(tab_cei_by_id)

# ---------- 2) Table: Descriptive stats of sub-indices (pooled) ----------
tab_sub_pooled <- df_month2 %>%
  select(all_of(sub_cols_10)) %>%
  pivot_longer(everything(), names_to = "sub_index", values_to = "value") %>%
  summarise(
    n    = sum(is.finite(value)),
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value, na.rm = TRUE),
    min  = min(value, na.rm = TRUE),
    p25  = quantile(value, 0.25, na.rm = TRUE),
    med  = median(value, na.rm = TRUE),
    p75  = quantile(value, 0.75, na.rm = TRUE),
    max  = max(value, na.rm = TRUE),
    .by  = sub_index
  ) %>%
  arrange(desc(sd))

print(tab_sub_pooled)

# ---------- 3) Table: Correlation matrix of sub-indices (pooled) ----------
# (use complete cases for a clean correlation matrix)
X_sub <- df_month2 %>% select(all_of(sub_cols_10))
X_cc  <- X_sub[complete.cases(X_sub), , drop=FALSE]

corr_sub <- cor(X_cc, use = "pairwise.complete.obs")
print(round(corr_sub, 3))

# Optional: correlation between CEI and sub-indices (needs join on id, ym)
df_join <- df_emh %>%
  select(id, ym, pc) %>%
  inner_join(df_month %>% select(id, ym, any_of(sub_cols_10)), by = c("id","ym"))

corr_pc_sub <- tibble(
  sub_index = sub_cols_10,
  corr_pc   = sapply(sub_cols_10, function(v) cor(df_join$pc, df_join[[v]], use="pairwise.complete.obs"))
) %>% arrange(desc(abs(corr_pc)))

print(corr_pc_sub)

# ---------- 4) Figures ----------
# Output directory
fig_dir <- file.path(out_dir, "dib_figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# 4.1 Time-series plot of CEI (6 markets)
p_ts <- ggplot(df_emh2, aes(x = date, y = pc, group = id)) +
  geom_line(linewidth = 0.6, alpha = 0.9) +
  facet_wrap(~ id, ncol = 2, scales = "free_y") +
  labs(
    x = NULL, y = "Composite Efficiency Index (CEI, PCA score)",
    title = "Time-varying Composite Efficiency Index (CEI) in ASEAN FX markets"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

print(p_ts)
ggsave(filename = file.path(fig_dir, "Fig1_CEI_timeseries_facets.png"),
       plot = p_ts, width = 10, height = 7, dpi = 300)

# 4.2 Distribution by market (boxplot)
p_box <- ggplot(df_emh2, aes(x = reorder(id, pc, median, na.rm = TRUE), y = pc)) +
  geom_boxplot(outlier.alpha = 0.25) +
  coord_flip() +
  labs(
    x = NULL, y = "CEI (PCA score)",
    title = "Distribution of CEI by market"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

print(p_box)
ggsave(filename = file.path(fig_dir, "Fig2_CEI_boxplot_by_market.png"),
       plot = p_box, width = 9, height = 6, dpi = 300)

# 4.3 Average CEI by market (bar + CI optional)
tab_cei_mean <- df_emh2 %>%
  group_by(id) %>%
  summarise(
    mean = mean(pc, na.rm = TRUE),
    se   = sd(pc, na.rm = TRUE)/sqrt(sum(is.finite(pc))),
    .groups = "drop"
  ) %>% arrange(mean)

p_bar <- ggplot(tab_cei_mean, aes(x = reorder(id, mean), y = mean)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se), linewidth = 0.5) +
  coord_flip() +
  labs(
    x = NULL, y = "Mean CEI (± 95% CI)",
    title = "Average CEI by market"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

print(p_bar)
ggsave(filename = file.path(fig_dir, "Fig3_CEI_mean_by_market.png"),
       plot = p_bar, width = 9, height = 6, dpi = 300)

# ---------- 5) Save tables (csv) ----------
write.csv(tab_cei_by_id,  file.path(fig_dir, "Table_CEI_desc_by_market.csv"), row.names = FALSE)
write.csv(tab_sub_pooled, file.path(fig_dir, "Table_SubIndices_desc_pooled.csv"), row.names = FALSE)
write.csv(round(corr_sub, 4), file.path(fig_dir, "Table_SubIndices_corr_matrix.csv"))
write.csv(corr_pc_sub, file.path(fig_dir, "Table_Corr_CEI_vs_SubIndices.csv"), row.names = FALSE)

cat("\n[Data in Brief outputs saved]\n")
cat("Tables & figures folder: ", fig_dir, "\n")

pca_sum <- summary(pc)
evr1 <- pca_sum$importance[2,1]   # Proportion of Variance of PC1
cat("Explained variance PC1 (pooled PCA):", round(evr1*100, 1), "%\n")

