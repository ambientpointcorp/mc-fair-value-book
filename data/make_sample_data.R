# ============================================================================
# Sample dataset generator — Monte Carlo Fair Value Analysis for Loan Pools
# ============================================================================
# Expands the anonymized pool extract (sample_loan_tape.csv) into the two
# estimation datasets used throughout Parts II-III:
#
#   sample_loan_panel.csv    one row per loan, survival format (2008-2025)
#   sample_monthly_panel.csv one row per loan-month, 2018-01 .. 2025-11
#
# Loan characteristics are bootstrapped from the tape rows; performance
# histories are simulated from seasoning/incentive/vintage hazards whose
# levels are calibrated to the aggregate behavior documented in the SBA
# 7(a) FOIA dataset (seasoning ramp over 12-18 months, long-run CPR near
# 10%, lifetime defaults in the low single digits with elevated
# 2008-2009 and 2020 vintages). Everything is seeded and deterministic.
#
# Usage: Rscript make_sample_data.R   (from the code/data directory)
# ============================================================================

set.seed(7031)
n_loans   <- 1500
report_ym <- c(2025, 11)

tape <- read.csv("sample_loan_tape.csv", stringsAsFactors = FALSE)

# ── Approximate public rate/macro history (annual averages) ────────────────
yrs   <- 2008:2025
prime <- c(5.1, 3.25, 3.25, 3.25, 3.25, 3.25, 3.25, 3.3, 3.5, 4.1,
           4.9, 5.3, 3.5, 3.25, 4.9, 8.2, 8.5, 7.6) / 100
unemp <- c(5.8, 9.3, 9.6, 8.9, 8.1, 7.4, 6.2, 5.3, 4.9, 4.4,
           3.9, 3.7, 8.1, 5.4, 3.6, 3.6, 4.0, 4.2) / 100
names(prime) <- names(unemp) <- yrs

sector_grp_map <- c(
  "72" = "Accommodation and Food", "44" = "Retail", "45" = "Retail",
  "62" = "Health Care", "81" = "Services", "56" = "Services",
  "54" = "Services", "23" = "Construction", "11" = "Agriculture"
)
grp_of <- function(s) ifelse(is.na(sector_grp_map[s]), "Other",
                             sector_grp_map[s])

# ── Bootstrap loan characteristics from the tape ───────────────────────────
idx  <- sample(nrow(tape), n_loans, replace = TRUE)
base <- tape[idx, c("state", "sector", "term_months", "loan_spread",
                    "orig_balance")]
base$loan_spread <- pmin(0.0375, pmax(0.005,
  base$loan_spread + sample(c(-0.0025, 0, 0.0025), n_loans, TRUE)))
base$orig_balance <- round(base$orig_balance *
  exp(rnorm(n_loans, 0, 0.35)) / 100) * 100
base$orig_balance <- pmin(5e6, pmax(5e4, base$orig_balance))

vint_w <- c(rep(1, 2), rep(1.2, 4), rep(1.6, 6), rep(2.4, 5))  # 2008..2024
orig_year  <- sample(2008:2024, n_loans, TRUE, prob = vint_w)
orig_month <- sample(1:12, n_loans, TRUE)

loans <- data.frame(
  loan_id      = sprintf("S%04d", seq_len(n_loans)),
  year         = orig_year,
  orig_month   = orig_month,
  state        = base$state,
  sector       = as.character(base$sector),
  sector_grp   = grp_of(as.character(base$sector)),
  term_months  = base$term_months,
  spread       = base$loan_spread,
  orig_balance = base$orig_balance,
  stringsAsFactors = FALSE
)

# ── Hazards ─────────────────────────────────────────────────────────────────
smm_at <- function(age, inc, grp) {
  ramp <- pmin(age / 18, 1)
  mult <- ifelse(grp == "Accommodation and Food", 1.10,
          ifelse(grp == "Retail", 1.00, 0.95))
  pmin(0.06, 0.009 * ramp * mult * exp(pmin(2, 25 * inc)))
}
mdr_at <- function(age, year, grp, st_unemp) {
  shape <- dlnorm(pmax(age, 1), log(36), 0.75) / dlnorm(36, log(36), 0.75)
  vint  <- ifelse(year %in% 2008:2009, 2.5, ifelse(year == 2020, 1.4, 1))
  mult  <- ifelse(grp == "Accommodation and Food", 1.25, 1)
  pmin(0.03, 0.0025 * shape * vint * mult * exp(15 * (st_unemp - 0.05)))
}

# ── Simulate monthly performance to the reporting date ─────────────────────
months_to_report <- function(y, m) (report_ym[1] - y) * 12 + (report_ym[2] - m)

# Coupons float on Prime (coupon_t = prime_t + spread) and payments
# re-amortize at each step, so the refinance incentive -- the floating
# coupon against a fixed-rate alternative -- RISES as Prime rises,
# matching the variable-rate prepayment response documented in ch.12.
sim_one <- function(i) {
  L <- loans[i, ]
  n_obs <- min(L$term_months, months_to_report(L$year, L$orig_month))
  if (n_obs < 1) return(NULL)
  bal <- L$orig_balance
  st_jit <- rnorm(1, 0, 0.004)
  out <- vector("list", n_obs)
  for (a in seq_len(n_obs)) {
    cal_y <- L$year + (L$orig_month + a - 1) %/% 12
    cal_m <- (L$orig_month + a - 1) %% 12 + 1
    yk    <- as.character(min(cal_y, 2025))
    stu   <- max(0.025, unemp[yk] + st_jit)
    cpn   <- prime[yk] + L$spread                 # floating coupon
    fixed_alt <- 0.04 + 0.5 * prime[yk]           # fixed-rate refi proxy
    inc   <- cpn - fixed_alt
    smm   <- smm_at(a, inc, L$sector_grp)
    mdr   <- mdr_at(a, L$year, L$sector_grp, stu)
    u     <- runif(1)
    ev    <- if (u < mdr) 2L else if (u < mdr + smm) 1L else 0L
    r     <- cpn / 12
    rem   <- L$term_months - a + 1
    pay   <- bal * r / (1 - (1 + r)^(-rem))
    int   <- bal * r
    bal   <- max(0, bal - (pay - int))
    out[[a]] <- c(a, cal_y, cal_m, bal, smm, mdr, inc, stu, ev, cpn)
    if (ev > 0L || bal <= 0) break
  }
  m <- do.call(rbind, out)
  colnames(m) <- c("age", "cal_y", "cal_m", "balance", "smm_true",
                   "mdr_true", "inc", "st_unemp", "ev", "cpn")
  cbind(loan = i, m)
}

hist_all <- do.call(rbind, lapply(seq_len(n_loans), sim_one))
hist_all <- as.data.frame(hist_all)

# ── Loan-level survival dataset ─────────────────────────────────────────────
last <- hist_all[!duplicated(hist_all$loan, fromLast = TRUE), ]
L    <- loans[last$loan, ]
window_start <- 2015
entry <- pmax(0, (window_start - L$year) * 12 - L$orig_month + 1)

loan_panel <- data.frame(
  loan_id       = L$loan_id,
  year          = L$year,
  sector        = L$sector,
  sector_grp    = L$sector_grp,
  state         = L$state,
  coupon        = round(last$cpn, 5),
  coupon_rate   = round(last$cpn, 5),
  term_months   = L$term_months,
  orig_balance  = L$orig_balance,
  balance       = pmax(round(last$balance), 1),
  balance_log   = log(pmax(round(last$balance), 1)),
  loan_age      = last$age,
  time          = last$age,
  age_months    = last$age,
  end_month     = last$age,
  end_age       = last$age,
  entry_age     = pmin(entry, last$age - 1),
  exit_age      = last$age,
  age_at_start  = pmin(entry, last$age - 1),
  prepay_event  = as.integer(last$ev == 1),
  prepay_flag   = as.integer(last$ev == 1),
  default_event = as.integer(last$ev == 2),
  event         = as.integer(last$ev > 0),
  event_type    = last$ev,
  status        = c("CURRENT", "PAID_IN_FULL", "CHARGE_OFF")[last$ev + 1],
  state_unemp   = round(last$st_unemp, 4),
  stringsAsFactors = FALSE
)
loan_panel$entry_age <- pmax(0, loan_panel$entry_age)
write.csv(loan_panel, "sample_loan_panel.csv", row.names = FALSE)

# ── Loan-month panel, 2018-01 .. 2025-11 ───────────────────────────────────
mp <- hist_all[hist_all$cal_y >= 2018, ]
Lm <- loans[mp$loan, ]
monthly_panel <- data.frame(
  loan_id      = Lm$loan_id,
  month        = as.Date(sprintf("%d-%02d-01", mp$cal_y, mp$cal_m)),
  month_id     = (mp$cal_y - 2018) * 12 + mp$cal_m,
  year         = mp$cal_y,
  age          = mp$age,
  loan_age     = mp$age,
  coupon       = round(mp$cpn, 5),
  prime        = as.numeric(prime[as.character(pmin(mp$cal_y, 2025))]),
  refi_incentive = round(mp$inc, 5),
  macro_index  = round(as.numeric(scale(mp$st_unemp)), 4),
  sector       = Lm$sector,
  state        = Lm$state,
  balance      = pmax(round(mp$balance), 1),
  balance_log  = log(pmax(round(mp$balance), 1)),
  prepay_event = as.integer(mp$ev == 1),
  smm_event    = as.integer(mp$ev == 1),
  event_type   = c("alive", "prepay", "default")[mp$ev + 1],
  smm          = round(pmax(1e-5, mp$smm_true *
                              exp(rnorm(nrow(mp), 0, 0.35))), 6),
  stringsAsFactors = FALSE
)
monthly_panel$log_smm <- round(log(monthly_panel$smm), 5)
write.csv(monthly_panel, "sample_monthly_panel.csv", row.names = FALSE)

cat("loan_panel:", nrow(loan_panel), "loans |",
    "prepaid", sum(loan_panel$prepay_event),
    "| defaulted", sum(loan_panel$default_event), "\n")
cat("monthly_panel:", nrow(monthly_panel), "loan-months |",
    "monthly prepay rate", round(mean(monthly_panel$prepay_event), 4), "\n")
