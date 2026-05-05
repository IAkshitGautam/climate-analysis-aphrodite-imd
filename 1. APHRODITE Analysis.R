# ======================================================================
# APHRODITE TEMPERATURE & PRECIPITATION & RSN ANALYSIS
# ======================================================================

# ----------------------------------------------------------------------
# 0. Packages
# ----------------------------------------------------------------------
pkgs <- c("terra", "trend", "ggplot2", "dplyr", "patchwork", "scales", "sf", "viridis", "R.utils")
new <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))
if (!requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils")
library(R.utils)

# ----------------------------------------------------------------------
# 1. Configuration
# ----------------------------------------------------------------------
base_dir <- "D:/RD539---Assignment"

cfg <- list(
  temp_raw_dir      = file.path(base_dir, "1. Data/APHRO_MA_TAVE_025deg_V1808/"),
  prcp_v1101_dir    = file.path(base_dir, "1. Data/APHRO_MA_025deg_V1101/"),
  prcp_v1101ex_dir  = file.path(base_dir, "1. Data/APHRO_MA_025deg_V1101_EXR1/"),
  shp_path          = file.path(base_dir, "madhya_pradesh.shp"),
  temp_out_dir      = file.path(base_dir, "Output_1-Temperature/"),
  prcp_out_dir      = file.path(base_dir, "Output_2-Precipitation/"),
  years             = 1981:2015,
  base_years        = 1981:2010,
  years_v1101       = 1981:2007,
  years_v1101ex     = 2008:2015,
  nx                = 360L,
  ny                = 280L,
  lon               = seq(60.125, 149.875, by = 0.25),
  lat               = seq(-14.875, 54.875, by = 0.25),
  miss              = -99.9,
  temp_tmpl         = "APHRO_MA_TAVE_025deg_V1808.%d",
  v1101_tmpl        = "APHRO_MA_025deg_V1101.%d",
  v1101ex_tmpl      = "APHRO_MA_025deg_V1101_EXR1.%d",
  cores             = max(1L, tryCatch(parallel::detectCores(logical = FALSE) - 1L, error = function(e) 1L))
)

force_rebuild <- FALSE
MON <- month.abb

for (d in c("nc", "maps", "timeseries", "trends", "variability", "qa")) {
  dir.create(file.path(cfg$temp_out_dir, d), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(cfg$prcp_out_dir, d), recursive = TRUE, showWarnings = FALSE)
}

terra::terraOptions(progress = 1, threads = cfg$cores)

# ----------------------------------------------------------------------
# 2. Generic helpers
# ----------------------------------------------------------------------
monthly_names <- function(var, years) {
  unlist(lapply(years, function(yr) sprintf("%s_%04d_%02d", var, yr, 1:12)), use.names = FALSE)
}

fix_layer_names <- function(stk, var, years) {
  exp_names <- monthly_names(var, years)
  if (terra::nlyr(stk) == length(exp_names)) names(stk) <- exp_names
  stk
}

month_days_year <- function(yr) {
  leap <- ((yr %% 4 == 0 & yr %% 100 != 0) | yr %% 400 == 0)
  c(31, if (leap) 29 else 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
}

mat_to_rast <- function(mat, cfg) {
  r_mat <- t(mat)[cfg$ny:1, ]
  ext <- terra::ext(
    min(cfg$lon) - 0.125, max(cfg$lon) + 0.125,
    min(cfg$lat) - 0.125, max(cfg$lat) + 0.125
  )
  r <- terra::rast(nrows = cfg$ny, ncols = cfg$nx, ext = ext, crs = "EPSG:4326")
  terra::values(r) <- c(t(r_mat))
  r
}

rast_df <- function(r) {
  df <- as.data.frame(r, xy = TRUE)
  names(df)[3] <- "value"
  df
}

monthly_layers <- function(stk, var, years, m) {
  idx <- match(sprintf("%s_%04d_%02d", var, years, m), names(stk))
  if (anyNA(idx)) stop(sprintf("Missing monthly layers for %s month %02d.", var, m))
  stk[[idx]]
}

annual_sum_from_monthly <- function(monthly_stack, var, years) {
  out <- lapply(years, function(yr) {
    idx <- match(sprintf("%s_%04d_%02d", var, yr, 1:12), names(monthly_stack))
    if (anyNA(idx)) stop(sprintf("Annual aggregation failed for %s %d.", var, yr))
    r <- sum(monthly_stack[[idx]], na.rm = TRUE)
    names(r) <- as.character(yr)
    r
  })
  combine_spatrasters(out)
}

annual_weighted_from_monthly <- function(monthly_stack, var, years) {
  out <- lapply(years, function(yr) {
    idx <- match(sprintf("%s_%04d_%02d", var, yr, 1:12), names(monthly_stack))
    if (anyNA(idx)) stop(sprintf("Annual aggregation failed for %s %d.", var, yr))
    w <- month_days_year(yr)
    layers <- lapply(seq_along(idx), function(i) monthly_stack[[idx[i]]] * w[i])
    r <- Reduce(`+`, layers) / sum(w)
    names(r) <- as.character(yr)
    r
  })
  combine_spatrasters(out)
}

monthly_climatology <- function(monthly_stack, var, years) {
  out <- lapply(1:12, function(m) {
    idx <- match(sprintf("%s_%04d_%02d", var, years, m), names(monthly_stack))
    if (anyNA(idx)) stop(sprintf("Missing baseline layers for %s month %02d.", var, m))
    r <- mean(monthly_stack[[idx]], na.rm = TRUE)
    names(r) <- MON[m]
    r
  })
  combine_spatrasters(out)
}

annual_climatology <- function(month_clim, kind = c("temp", "prcp", "rsn")) {
  kind <- match.arg(kind)
  if (kind == "prcp") {
    terra::app(month_clim, function(x) sum(x, na.rm = TRUE))
  } else {
    terra::mean(month_clim, na.rm = TRUE)
  }
}

safe_cv <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (is.na(m) || abs(m) < 1e-8 || is.na(s)) return(NA_real_)
  100 * s / abs(m)
}

mk_sen_pixel <- function(x) {
  if (sum(!is.na(x)) < 5) return(c(NA_real_, NA_real_, NA_real_))
  mk <- tryCatch(trend::mk.test(x), error = function(e) NULL)
  sl <- tryCatch(trend::sens.slope(x), error = function(e) NULL)
  if (is.null(mk) || is.null(sl)) return(c(NA_real_, NA_real_, NA_real_))
  c(unname(mk$statistic), mk$p.value, unname(sl$estimates))
}

compute_trend <- function(stk) {
  r3 <- terra::app(stk, mk_sen_pixel, cores = cfg$cores)
  names(r3) <- c("mk_z", "pvalue", "slope")
  r3
}

trend_summary <- function(tr, label, alpha = 0.05) {
  pv <- terra::values(tr[["pvalue"]], na.rm = FALSE)
  sl <- terra::values(tr[["slope"]],  na.rm = FALSE)
  ok <- !is.na(pv)
  n_ok <- sum(ok)
  n_sig <- sum(pv[ok] < alpha, na.rm = TRUE)
  sig_sl <- sl[ok & !is.na(sl) & pv < alpha]
  data.frame(
    Variable     = label,
    N_pixels     = n_ok,
    Sig_pct      = round(if (n_ok == 0) NA_real_ else 100 * n_sig / n_ok, 1),
    Pos_sig_pct  = round(if (n_sig == 0) NA_real_ else 100 * sum(sig_sl > 0, na.rm = TRUE) / n_sig, 1),
    Median_slope = round(median(sl[ok], na.rm = TRUE), 4),
    stringsAsFactors = FALSE
  )
}

combine_spatrasters <- function(x) {
  x <- Filter(function(z) inherits(z, "SpatRaster") && terra::nlyr(z) == 1, x)
  if (!length(x)) stop("No valid SpatRaster layers to combine.")
  Reduce(c, x)
}

symmetric_limits <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(-1, 1))
  m <- max(abs(range(x, na.rm = TRUE)))
  if (!is.finite(m) || m == 0) m <- 1
  c(-m, m)
}

# ----------------------------------------------------------------------
# 3. Maps
# ----------------------------------------------------------------------
base_map <- function(df, mp_sf) {
  ggplot() +
    geom_raster(data = df, aes(x, y, fill = value), na.rm = TRUE) +
    geom_sf(data = mp_sf, fill = NA, color = "black", linewidth = 0.35) +
    coord_sf(expand = FALSE) +
    labs(x = "Longitude", y = "Latitude") +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.title = element_text(color = "black"),
      axis.text  = element_text(color = "black", size = 8),
      axis.ticks = element_line(color = "black"),
      axis.line  = element_line(color = "black", linewidth = 0.35),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text  = element_text(size = 10),
      legend.key.height = grid::unit(0.45, "cm"),
      legend.key.width  = grid::unit(0.35, "cm")
    )
}

fill_temp <- function(limits, name = "°C") {
  scale_fill_distiller(
    palette = "YlOrRd", direction = 1,
    limits = limits, oob = scales::squish,
    na.value = "grey90", name = name,
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5,
                           barheight = grid::unit(4.8, "cm"),
                           barwidth  = grid::unit(0.4, "cm"))
  )
}

fill_prcp <- function(limits, name = "mm") {
  scale_fill_distiller(
    palette = "YlGnBu", direction = 1,
    limits = limits, oob = scales::squish,
    na.value = "grey90", name = name,
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5,
                           barheight = grid::unit(4.8, "cm"),
                           barwidth  = grid::unit(0.4, "cm"))
  )
}

fill_rsn <- function(limits, name = "RSN (%)") {
  scale_fill_viridis_c(
    option = "C", limits = limits, oob = scales::squish,
    na.value = "grey90", name = name,
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5,
                           barheight = grid::unit(4.8, "cm"),
                           barwidth  = grid::unit(0.4, "cm"))
  )
}

fill_div <- function(limits, name = "") {
  scale_fill_gradient2(
    low = "firebrick", mid = "white", high = "royalblue", midpoint = 0,
    limits = limits, oob = scales::squish,
    na.value = "grey90", name = name,
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5,
                           barheight = grid::unit(4.8, "cm"),
                           barwidth  = grid::unit(0.45, "cm"))
  )
}

fill_pval <- function(limits) {
  scale_fill_viridis_c(
    option = "D", direction = -1,
    limits = limits, oob = scales::squish,
    na.value = "grey90", name = "p-value",
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5,
                           barheight = grid::unit(4.8, "cm"),
                           barwidth  = grid::unit(0.4, "cm"))
  )
}

fill_cv <- function(limits) {
  scale_fill_viridis_c(
    option = "C", limits = limits, oob = scales::squish,
    na.value = "grey90", name = "CV (%)",
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5,
                           barheight = grid::unit(4.8, "cm"),
                           barwidth  = grid::unit(0.4, "cm"))
  )
}

assemble_panel <- function(plots, title, ncol = 4) {
  p <- patchwork::wrap_plots(plots, ncol = ncol, guides = "collect")
  p <- p + plot_annotation(title = title, theme = theme(plot.title = element_text(size = 14, face = "bold")))
  p <- p & theme(legend.position = "right")
  print(p)
  invisible(p)
}

# ----------------------------------------------------------------------
# 4. Reading APHRODITE daily binary files (field1 + RSN)
# ----------------------------------------------------------------------
read_aphro_pair_year <- function(filepath, cfg, mp_vect, kind = c("temp", "prcp")) {
  kind <- match.arg(kind)
  stopifnot(file.exists(filepath))
  
  yr <- as.integer(sub(".*\\.(\\d{4})$", "\\1", filepath))
  nday <- if ((yr %% 4 == 0 & yr %% 100 != 0) | yr %% 400 == 0) 366L else 365L
  dates <- seq(as.Date(sprintf("%d-01-01", yr)), by = "day", length.out = nday)
  mon <- as.integer(format(dates, "%m"))
  
  n_cell <- cfg$nx * cfg$ny
  n_per_rec <- 2L * n_cell  # field1 + rsn
  
  f_sum <- array(0, dim = c(cfg$nx, cfg$ny, 12))
  f_cnt <- array(0L, dim = c(cfg$nx, cfg$ny, 12))
  
  rsn_sum <- array(0, dim = c(cfg$nx, cfg$ny, 12))
  rsn_cnt <- array(0L, dim = c(cfg$nx, cfg$ny, 12))
  
  daily_field1 <- numeric(nday)
  daily_rsn <- numeric(nday)
  
  con <- if (grepl("\\.gz$", filepath)) gzfile(filepath, "rb") else file(filepath, "rb")
  on.exit(close(con), add = TRUE)
  
  for (d in seq_len(nday)) {
    raw <- readBin(con, what = "numeric", n = n_per_rec, size = 4L, endian = "little")
    if (length(raw) < n_per_rec) stop("Unexpected end of file while reading: ", filepath)
    
    f1  <- raw[seq_len(n_cell)]
    rsn <- raw[n_cell + seq_len(n_cell)]
    
    f1[f1 <= cfg$miss + 0.01]   <- NA_real_
    rsn[rsn <= cfg$miss + 0.01] <- NA_real_
    
    mat1 <- matrix(f1, nrow = cfg$nx, ncol = cfg$ny)
    matr <- matrix(rsn, nrow = cfg$nx, ncol = cfg$ny)
    
    r1_full <- mat_to_rast(mat1, cfg)
    rr_full <- mat_to_rast(matr, cfg)
    
    r1_mp <- terra::crop(r1_full, mp_vect) |> terra::mask(mp_vect)
    rr_mp <- terra::crop(rr_full, mp_vect) |> terra::mask(mp_vect)
    
    vals1 <- terra::values(r1_mp, mat = FALSE)
    valsr <- terra::values(rr_mp, mat = FALSE)
    
    daily_field1[d] <- if (all(is.na(vals1))) NA_real_ else mean(vals1, na.rm = TRUE)
    daily_rsn[d]    <- if (all(is.na(valsr))) NA_real_ else mean(valsr, na.rm = TRUE)
    
    m <- mon[d]
    valid1 <- !is.na(mat1)
    validr <- !is.na(matr)
    
    if (kind == "prcp") {
      f_sum[, , m] <- f_sum[, , m] + ifelse(valid1, mat1, 0)
    } else {
      f_sum[, , m] <- f_sum[, , m] + ifelse(valid1, mat1, 0)
      f_cnt[, , m] <- f_cnt[, , m] + as.integer(valid1)
    }
    
    rsn_sum[, , m] <- rsn_sum[, , m] + ifelse(validr, matr, 0)
    rsn_cnt[, , m] <- rsn_cnt[, , m] + as.integer(validr)
  }
  
  monthly_field1 <- array(NA_real_, c(cfg$nx, cfg$ny, 12))
  monthly_rsn    <- array(NA_real_, c(cfg$nx, cfg$ny, 12))
  
  for (m in 1:12) {
    if (kind == "prcp") {
      monthly_field1[, , m] <- f_sum[, , m]
    } else {
      monthly_field1[, , m] <- ifelse(f_cnt[, , m] > 0, f_sum[, , m] / f_cnt[, , m], NA_real_)
    }
    monthly_rsn[, , m] <- ifelse(rsn_cnt[, , m] > 0, rsn_sum[, , m] / rsn_cnt[, , m], NA_real_)
  }
  
  daily_df <- data.frame(
    date = dates,
    year = as.integer(format(dates, "%Y")),
    month = as.integer(format(dates, "%m")),
    day = as.integer(format(dates, "%d")),
    doy = as.integer(format(dates, "%j")),
    field1 = daily_field1,
    rsn_mp = daily_rsn
  )
  if (kind == "temp") names(daily_df)[names(daily_df) == "field1"] <- "tmean_mp"
  if (kind == "prcp") names(daily_df)[names(daily_df) == "field1"] <- "prcp_mp"
  
  list(monthly_field1 = monthly_field1, monthly_rsn = monthly_rsn, daily = daily_df, year = yr)
}

# ----------------------------------------------------------------------
# 5. Caching helpers
# ----------------------------------------------------------------------
check_mp_values <- function(r_full, mp_vect, var_label) {
  r_mp <- terra::crop(r_full, mp_vect) |> terra::mask(mp_vect)
  mn <- terra::global(r_mp, "mean", na.rm = TRUE)[1, 1]
  mn2 <- terra::global(r_mp, "min",  na.rm = TRUE)[1, 1]
  mx  <- terra::global(r_mp, "max",  na.rm = TRUE)[1, 1]
  message(sprintf("  [DIAG] %s -> mean=%.2f  min=%.2f  max=%.2f", var_label, mn, mn2, mx))
}

load_stack_if_valid <- function(path, expected_names = NULL) {
  if (!file.exists(path)) return(NULL)
  r <- tryCatch(terra::rast(path), error = function(e) NULL)
  if (is.null(r)) return(NULL)
  if (!is.null(expected_names)) {
    if (terra::nlyr(r) != length(expected_names)) return(NULL)
    names(r) <- expected_names
  }
  r
}

load_or_build_product <- function(spec, cfg, kind = c("temp", "prcp"), force_rebuild = FALSE) {
  kind <- match.arg(kind)
  
  out_dir   <- spec$out_dir
  prefix    <- spec$prefix
  years     <- spec$years
  file_fun  <- spec$file_fun
  value_col <- spec$value_col
  
  main_cache <- file.path(out_dir, "nc", sprintf("%s_monthly_MP.tif", prefix))
  rsn_cache  <- file.path(out_dir, "nc", "rsn_monthly_MP.tif")
  daily_csv  <- file.path(out_dir, "timeseries", sprintf("%s_daily_mean_mp.csv", prefix))
  
  if (force_rebuild) {
    if (file.exists(main_cache)) unlink(main_cache)
    if (file.exists(rsn_cache)) unlink(rsn_cache)
    if (file.exists(daily_csv)) unlink(daily_csv)
  }
  
  main_names <- monthly_names(prefix, years)
  rsn_names  <- monthly_names("rsn", years)
  
  if (file.exists(main_cache) && file.exists(rsn_cache) && file.exists(daily_csv)) {
    main_stk <- load_stack_if_valid(main_cache, main_names)
    rsn_stk  <- load_stack_if_valid(rsn_cache, rsn_names)
    daily_df <- tryCatch(read.csv(daily_csv), error = function(e) NULL)
    
    if (!is.null(main_stk) && !is.null(rsn_stk) && !is.null(daily_df)) {
      daily_df$date <- as.Date(daily_df$date)
      if (value_col %in% names(daily_df) && "rsn_mp" %in% names(daily_df)) {
        daily_df$month_name <- factor(MON[daily_df$month], levels = MON)
        return(list(main_stack = main_stk, rsn_stack = rsn_stk, daily = daily_df))
      }
    }
  }
  
  mp <- terra::vect(cfg$shp_path) |> terra::project("EPSG:4326")
  
  main_layers <- list()
  rsn_layers   <- list()
  daily_list   <- list()
  diag_done    <- FALSE
  
  for (yr in years) {
    f <- file_fun(yr)
    fgz <- paste0(f, ".gz")
    rm_tmp <- FALSE
    
    if (!file.exists(f) && file.exists(fgz)) {
      message(sprintf("  Decompressing %s %d ...", prefix, yr))
      tmp <- tempfile(fileext = sprintf(".%d", yr))
      R.utils::gunzip(fgz, destname = tmp, remove = FALSE, overwrite = TRUE)
      f <- tmp
      rm_tmp <- TRUE
    }
    
    if (!file.exists(f)) {
      warning("Missing file: ", f)
      next
    }
    
    message(sprintf("  Reading %s %d ...", prefix, yr))
    res <- read_aphro_pair_year(f, cfg, mp, kind = kind)
    daily_list[[length(daily_list) + 1]] <- res$daily
    
    if (rm_tmp) unlink(f)
    
    for (m in 1:12) {
      r_main_full <- mat_to_rast(res$monthly_field1[, , m], cfg)
      r_rsn_full  <- mat_to_rast(res$monthly_rsn[, , m], cfg)
      
      if (!diag_done && m == 1) {
        check_mp_values(r_main_full, mp, sprintf("%s %d Jan", prefix, yr))
        check_mp_values(r_rsn_full,  mp, sprintf("rsn %d Jan", yr))
        diag_done <- TRUE
      }
      
      r_main <- terra::crop(r_main_full, mp) |> terra::mask(mp)
      r_rsn  <- terra::crop(r_rsn_full,  mp) |> terra::mask(mp)
      
      nm_main <- sprintf("%s_%04d_%02d", prefix, yr, m)
      nm_rsn  <- sprintf("rsn_%04d_%02d", yr, m)
      
      names(r_main) <- nm_main
      names(r_rsn)  <- nm_rsn
      
      main_layers[[nm_main]] <- r_main
      rsn_layers[[nm_rsn]]   <- r_rsn
    }
  }
  
  if (length(main_layers) == 0) stop("No layers created for ", prefix)
  
  main_stk <- combine_spatrasters(main_layers)
  rsn_stk  <- combine_spatrasters(rsn_layers)
  
  names(main_stk) <- names(main_layers)
  names(rsn_stk)  <- names(rsn_layers)
  
  main_stk <- main_stk[[intersect(main_names, names(main_stk))]]
  rsn_stk  <- rsn_stk[[intersect(rsn_names, names(rsn_stk))]]
  
  terra::writeRaster(main_stk, main_cache, overwrite = TRUE)
  terra::writeRaster(rsn_stk, rsn_cache, overwrite = TRUE)
  
  daily_df <- dplyr::bind_rows(daily_list) |>
    arrange(date) |>
    mutate(month_name = factor(MON[month], levels = MON))
  
  if (!(value_col %in% names(daily_df))) names(daily_df)[names(daily_df) == "field1"] <- value_col
  write.csv(daily_df, daily_csv, row.names = FALSE)
  
  list(main_stack = main_stk, rsn_stack = rsn_stk, daily = daily_df)
}

# ----------------------------------------------------------------------
# 6. Base plotting helpers
# ----------------------------------------------------------------------
sen_trend_band <- function(years, values) {
  ok <- is.finite(years) & is.finite(values)
  years <- years[ok]
  values <- values[ok]
  if (length(values) < 5) return(NULL)
  
  sl <- tryCatch(as.numeric(trend::sens.slope(values)$estimates), error = function(e) NA_real_)
  if (!is.finite(sl)) return(NULL)
  
  intercept <- median(values - sl * years, na.rm = TRUE)
  fit <- intercept + sl * years
  resid <- values - fit
  
  r1 <- tryCatch(stats::acf(resid, plot = FALSE, lag.max = 1, na.action = na.pass)$acf[2], error = function(e) NA_real_)
  if (!is.finite(r1)) r1 <- 0
  
  n <- length(resid)
  neff <- round(n * (1 - r1) / (1 + r1))
  neff <- max(3, min(n, neff))
  
  se <- stats::sd(resid, na.rm = TRUE) / sqrt(neff)
  crit <- stats::qt(0.975, df = max(1, neff - 2))
  
  data.frame(year = years, fit = fit, lwr = fit - crit * se, upr = fit + crit * se, slope = sl, neff = neff)
}

plot_daily_series <- function(df, value_col, roll_col = NULL, title, ylab) {
  plot(df$date, df[[value_col]], type = "l", col = "grey50", lwd = 1,
       main = title, xlab = "", ylab = ylab)
  
  if (!is.null(roll_col) && roll_col %in% names(df)) {
    lines(df$date, df[[roll_col]], lwd = 2)
  }
  invisible(NULL)
}

plot_annual_trend_series <- function(df, year_col, value_col, title, ylab) {
  years <- df[[year_col]]
  vals <- df[[value_col]]
  band <- sen_trend_band(years, vals)
  
  ylim <- range(c(vals, if (!is.null(band)) c(band$lwr, band$upr) else NULL), na.rm = TRUE)
  plot(years, vals, type = "n", main = title, xlab = "Year", ylab = ylab, ylim = ylim)
  
  if (!is.null(band)) {
    polygon(c(band$year, rev(band$year)), c(band$lwr, rev(band$upr)),
            col = adjustcolor("grey60", 0.25), border = NA)
    lines(band$year, band$fit, lwd = 2)
  }
  
  lines(years, vals, lwd = 1.2)
  points(years, vals, pch = 16, cex = 0.7)
  invisible(NULL)
}

plot_annual_anomaly_bars <- function(df, year_col, anom_col, title, ylab) {
  years <- df[[year_col]]
  anom <- df[[anom_col]]
  cols <- ifelse(anom >= 0, "royalblue", "firebrick")
  
  ylim <- range(c(anom, 0), na.rm = TRUE)
  plot(years, anom, type = "n", main = title, xlab = "Year", ylab = ylab, ylim = ylim)
  abline(h = 0, lwd = 1)
  
  for (i in seq_along(years)) {
    rect(years[i] - 0.35, 0, years[i] + 0.35, anom[i], col = cols[i], border = cols[i])
  }
  
  lines(years, anom, lwd = 0.8)
  points(years, anom, pch = 16, cex = 0.6)
  invisible(NULL)
}

plot_monthly_anomalies_base <- function(df, year_col, anom_col, title) {
  op <- par(mfrow = c(3, 4), mar = c(3, 3, 2.2, 0.8), oma = c(0, 0, 3, 0))
  on.exit(par(op), add = TRUE)
  
  lim <- max(abs(df[[anom_col]]), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 1
  
  for (m in 1:12) {
    d <- df[df$month == m, ]
    plot(d[[year_col]], d[[anom_col]], type = "l", main = MON[m], xlab = "", ylab = "",
         ylim = c(-lim, lim))
    abline(h = 0, lwd = 1)
  }
  mtext(title, outer = TRUE, font = 2, cex = 1.1)
  invisible(NULL)
}

plot_boxplot_base <- function(df, x_col, y_col, title, ylab) {
  boxplot(df[[y_col]] ~ df[[x_col]], notch = FALSE, las = 2, col = "grey80",
          main = title, xlab = "", ylab = ylab, outline = FALSE)
  invisible(NULL)
}

# ----------------------------------------------------------------------
# 7. Generic analysis runner
# ----------------------------------------------------------------------
run_product_analysis <- function(spec, cfg, kind = c("temp", "prcp"), force_rebuild = FALSE) {
  kind <- match.arg(kind)
  
  out_dir   <- spec$out_dir
  years     <- spec$years
  prefix    <- spec$prefix
  value_col <- spec$value_col
  value_unit <- spec$value_unit
  annual_unit <- spec$annual_unit
  
  message(sprintf("\n>>> Building/loading %s products <<<", prefix))
  obj <- load_or_build_product(spec, cfg, kind = kind, force_rebuild = force_rebuild)
  
  main_stk <- fix_layer_names(obj$main_stack, prefix, years)
  rsn_stk  <- fix_layer_names(obj$rsn_stack, "rsn", years)
  daily_df  <- obj$daily
  mp_sf     <- sf::st_read(cfg$shp_path, quiet = TRUE) |> sf::st_transform(4326)
  
  # ---------- QA ----------
  qa_dir <- file.path(out_dir, "qa")
  dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)
  
  expected_dates <- seq(as.Date(sprintf("%d-01-01", min(years))),
                        as.Date(sprintf("%d-12-31", max(years))), by = "day")
  
  qa_report <- data.frame(
    item = c("monthly layers in stack", "daily rows", "duplicate dates", "missing dates in span", "NA values in daily series", "baseline years present"),
    value = c(
      terra::nlyr(main_stk),
      nrow(daily_df),
      sum(duplicated(daily_df$date)),
      length(setdiff(expected_dates, daily_df$date)),
      sum(is.na(daily_df[[value_col]])),
      all(cfg$base_years %in% unique(daily_df$year))
    ),
    stringsAsFactors = FALSE
  )
  write.csv(qa_report, file.path(qa_dir, sprintf("%s_data_health_check.csv", prefix)), row.names = FALSE)
  
  if (sum(duplicated(daily_df$date)) != 0) stop("Duplicate dates in ", prefix, " daily CSV.")
  if (sum(is.na(daily_df[[value_col]])) != 0) stop("NA values in ", prefix, " daily series.")
  
  # ---------- Daily preprocessing ----------
  daily_df <- daily_df |>
    arrange(date) |>
    mutate(
      month_name = factor(MON[month], levels = MON),
      roll30     = as.numeric(stats::filter(.data[[value_col]], rep(1 / 30, 30), sides = 2)),
      roll30_rsn = as.numeric(stats::filter(rsn_mp, rep(1 / 30, 30), sides = 2))
    )
  
  # ---------- State summaries ----------
  monthly_state <- daily_df |>
    group_by(year, month, month_name) |>
    summarise(
      mean_val = if (kind == "prcp") sum(.data[[value_col]], na.rm = TRUE) else mean(.data[[value_col]], na.rm = TRUE),
      sd_val   = sd(.data[[value_col]], na.rm = TRUE),
      min_val  = min(.data[[value_col]], na.rm = TRUE),
      max_val  = max(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    )
  
  annual_state <- daily_df |>
    group_by(year) |>
    summarise(
      mean_val = if (kind == "prcp") sum(.data[[value_col]], na.rm = TRUE) else mean(.data[[value_col]], na.rm = TRUE),
      sd_val   = sd(.data[[value_col]], na.rm = TRUE),
      min_val  = min(.data[[value_col]], na.rm = TRUE),
      max_val  = max(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    )
  
  rsn_monthly_state <- daily_df |>
    group_by(year, month, month_name) |>
    summarise(mean_rsn = mean(rsn_mp, na.rm = TRUE), sd_rsn = sd(rsn_mp, na.rm = TRUE), .groups = "drop")
  
  rsn_annual_state <- daily_df |>
    group_by(year) |>
    summarise(mean_rsn = mean(rsn_mp, na.rm = TRUE), sd_rsn = sd(rsn_mp, na.rm = TRUE), .groups = "drop")
  
  baseline_month_clim_state <- monthly_state |>
    filter(year %in% cfg$base_years) |>
    group_by(month, month_name) |>
    summarise(clim_val = mean(mean_val, na.rm = TRUE), .groups = "drop") |>
    arrange(month)
  
  baseline_annual_mean <- annual_state |>
    filter(year %in% cfg$base_years) |>
    summarise(clim_val = mean(mean_val, na.rm = TRUE)) |>
    pull(clim_val)
  
  baseline_rsn_month <- rsn_monthly_state |>
    filter(year %in% cfg$base_years) |>
    group_by(month, month_name) |>
    summarise(clim_rsn = mean(mean_rsn, na.rm = TRUE), .groups = "drop") |>
    arrange(month)
  
  baseline_rsn_annual <- rsn_annual_state |>
    filter(year %in% cfg$base_years) |>
    summarise(clim_rsn = mean(mean_rsn, na.rm = TRUE)) |>
    pull(clim_rsn)
  
  monthly_state <- monthly_state |>
    left_join(baseline_month_clim_state, by = c("month", "month_name")) |>
    mutate(anom = mean_val - clim_val)
  
  annual_state <- annual_state |>
    mutate(clim_val = baseline_annual_mean, anom = mean_val - clim_val)
  
  rsn_monthly_state <- rsn_monthly_state |>
    left_join(baseline_rsn_month, by = c("month", "month_name")) |>
    mutate(anom = mean_rsn - clim_rsn)
  
  rsn_annual_state <- rsn_annual_state |>
    mutate(clim_rsn = baseline_rsn_annual, anom = mean_rsn - clim_rsn)
  
  # Save state tables
  write.csv(daily_df, file.path(out_dir, "timeseries", sprintf("%s_daily_mean_mp.csv", prefix)), row.names = FALSE)
  write.csv(monthly_state, file.path(out_dir, "timeseries", sprintf("%s_state_monthly.csv", prefix)), row.names = FALSE)
  write.csv(annual_state, file.path(out_dir, "timeseries", sprintf("%s_state_annual.csv", prefix)), row.names = FALSE)
  write.csv(rsn_annual_state, file.path(out_dir, "timeseries", "rsn_state_annual.csv"), row.names = FALSE)
  
  # ---------- Daily graph ----------
  plot_daily_series(
    daily_df, value_col, "roll30",
    sprintf("Daily Mean %s — Madhya Pradesh", prefix),
    value_unit
  )
  
  # ---------- Climatology maps for main variable ----------
  base_idx <- match(monthly_names(prefix, cfg$base_years), names(main_stk))
  if (anyNA(base_idx)) stop("Baseline monthly layers missing in ", prefix, " stack.")
  base_stk <- main_stk[[base_idx]]
  names(base_stk) <- monthly_names(prefix, cfg$base_years)
  
  month_clim <- monthly_climatology(base_stk, prefix, cfg$base_years)
  annual_clim <- annual_climatology(month_clim, kind = kind)
  
  month_lim <- range(terra::values(month_clim, na.rm = TRUE), na.rm = TRUE)
  annual_lim <- range(terra::values(annual_clim, na.rm = TRUE), na.rm = TRUE)
  
  p_month_clim <- assemble_panel(
    lapply(1:12, function(m) {
      p <- base_map(rast_df(month_clim[[m]]), mp_sf)
      p <- p + (if (kind == "prcp") fill_prcp(month_lim, "mm/month") else fill_temp(month_lim, value_unit))
      p + labs(title = sprintf("%s\nState mean = %.2f %s",
                               MON[m],
                               baseline_month_clim_state$clim_val[match(m, baseline_month_clim_state$month)],
                               if (kind == "prcp") "mm/month" else value_unit))
    }),
    sprintf("Monthly %s Climatology (1981–2010) — Madhya Pradesh", prefix),
    ncol = 4
  )
  
  p_annual_clim <- base_map(rast_df(annual_clim), mp_sf) +
    (if (kind == "prcp") fill_prcp(annual_lim, annual_unit) else fill_temp(annual_lim, value_unit)) +
    labs(title = sprintf("Mean Annual %s (1981–2010)\nState mean = %.2f %s", prefix, baseline_annual_mean, annual_unit))
  print(p_annual_clim)
  
  terra::writeRaster(month_clim, file.path(out_dir, "nc", sprintf("%s_monthly_climatology_1981_2010.tif", prefix)), overwrite = TRUE)
  terra::writeRaster(annual_clim, file.path(out_dir, "nc", sprintf("%s_annual_climatology_1981_2010.tif", prefix)), overwrite = TRUE)
  
  # ---------- RSN annual climatology map + annual RSN series ----------
  rsn_base_idx <- match(monthly_names("rsn", cfg$base_years), names(rsn_stk))
  if (anyNA(rsn_base_idx)) stop("Baseline RSN layers missing in stack.")
  rsn_base_stk <- rsn_stk[[rsn_base_idx]]
  names(rsn_base_stk) <- monthly_names("rsn", cfg$base_years)
  
  rsn_month_clim <- monthly_climatology(rsn_base_stk, "rsn", cfg$base_years)
  rsn_annual_clim <- annual_climatology(rsn_month_clim, kind = "rsn")
  rsn_lim <- range(terra::values(rsn_annual_clim, na.rm = TRUE), na.rm = TRUE)
  rsn_lim <- c(max(0, floor(rsn_lim[1])), min(100, ceiling(rsn_lim[2])))
  
  p_rsn_annual <- base_map(rast_df(rsn_annual_clim), mp_sf) + fill_rsn(rsn_lim) +
    labs(title = sprintf("Mean Annual RSN (1981–2010) — %s", prefix))
  print(p_rsn_annual)
  
  terra::writeRaster(rsn_annual_clim, file.path(out_dir, "nc", "rsn_annual_climatology_1981_2010.tif"), overwrite = TRUE)
  
  plot_annual_trend_series(
    rsn_annual_state, "year", "mean_rsn",
    sprintf("Annual Mean RSN Time Series — %s", prefix),
    "RSN (%)"
  )
  
  # ---------- Trend analysis for main variable ----------
  annual_stack <- if (kind == "prcp") annual_sum_from_monthly(main_stk, prefix, years) else annual_weighted_from_monthly(main_stk, prefix, years)
  main_ann_trend <- compute_trend(annual_stack)
  main_month_trend <- lapply(1:12, function(m) compute_trend(monthly_layers(main_stk, prefix, years, m)))
  
  month_slope_vals <- unlist(lapply(main_month_trend, function(tr) terra::values(tr[["slope"]], na.rm = TRUE)))
  annual_slope_vals <- terra::values(main_ann_trend[["slope"]], na.rm = TRUE)
  month_mkz_vals <- unlist(lapply(main_month_trend, function(tr) terra::values(tr[["mk_z"]], na.rm = TRUE)))
  annual_mkz_vals <- terra::values(main_ann_trend[["mk_z"]], na.rm = TRUE)
  month_p_vals <- unlist(lapply(main_month_trend, function(tr) terra::values(tr[["pvalue"]], na.rm = TRUE)))
  annual_p_vals <- terra::values(main_ann_trend[["pvalue"]], na.rm = TRUE)
  
  month_slope_lim <- symmetric_limits(month_slope_vals)
  annual_slope_lim <- symmetric_limits(annual_slope_vals)
  mkz_lim <- symmetric_limits(c(month_mkz_vals, annual_mkz_vals))
  
  month_p_lim <- range(month_p_vals, na.rm = TRUE)
  annual_p_lim <- range(annual_p_vals, na.rm = TRUE)
  if (!all(is.finite(month_p_lim))) month_p_lim <- c(0, 1)
  if (!all(is.finite(annual_p_lim))) annual_p_lim <- c(0, 1)
  
  p_ann_trend <- {
    p1 <- base_map(rast_df(main_ann_trend[["slope"]]), mp_sf) +
      fill_div(annual_slope_lim, if (kind == "temp") "°C/yr" else "mm/yr") +
      labs(title = "Sen Slope")
    p2 <- base_map(rast_df(main_ann_trend[["mk_z"]]), mp_sf) +
      fill_div(mkz_lim, "MK Z") +
      labs(title = "Mann-Kendall Z")
    p3 <- base_map(rast_df(main_ann_trend[["pvalue"]]), mp_sf) +
      fill_pval(annual_p_lim) +
      labs(title = "p-value")
    (p1 | p2 | p3) + plot_annotation(title = sprintf("Annual %s Trend (1981–2015) — Madhya Pradesh", prefix))
  }
  print(p_ann_trend)
  
  p_month_slope <- assemble_panel(
    lapply(1:12, function(m) base_map(rast_df(main_month_trend[[m]][["slope"]]), mp_sf) + fill_div(month_slope_lim, if (kind == "temp") "°C/yr" else "mm/yr") + labs(title = MON[m])),
    sprintf("Monthly %s Sen Slope (1981–2015)", prefix),
    ncol = 4
  )
  
  p_month_mkz <- assemble_panel(
    lapply(1:12, function(m) base_map(rast_df(main_month_trend[[m]][["mk_z"]]), mp_sf) + fill_div(mkz_lim, "MK Z") + labs(title = MON[m])),
    sprintf("Monthly %s Mann-Kendall Z (1981–2015)", prefix),
    ncol = 4
  )
  
  p_month_p <- assemble_panel(
    lapply(1:12, function(m) base_map(rast_df(main_month_trend[[m]][["pvalue"]]), mp_sf) + fill_pval(month_p_lim) + labs(title = MON[m])),
    sprintf("Monthly %s p-value (1981–2015)", prefix),
    ncol = 4
  )
  
  plot_annual_trend_series(
    annual_state, "year", "mean_val",
    sprintf("State Mean Annual %s Time Series — Madhya Pradesh", prefix),
    annual_unit
  )
  
  # trend summary
  trend_summary_df <- rbind(
    trend_summary(main_ann_trend, sprintf("Annual %s", prefix)),
    do.call(rbind, lapply(1:12, function(m) trend_summary(main_month_trend[[m]], paste(prefix, MON[m]))))
  )
  write.csv(trend_summary_df, file.path(out_dir, "trends", "trend_significance_summary.csv"), row.names = FALSE)
  
  # ---------- Variability (main variable only) ----------
  annual_cv <- terra::app(annual_stack, safe_cv)
  names(annual_cv) <- "annual_cv"
  month_cv <- lapply(1:12, function(m) terra::app(monthly_layers(main_stk, prefix, years, m), safe_cv))
  
  month_cv_vals <- unlist(lapply(month_cv, function(r) terra::values(r, na.rm = TRUE)))
  annual_cv_vals <- terra::values(annual_cv, na.rm = TRUE)
  month_cv_lim <- c(0, max(month_cv_vals, na.rm = TRUE))
  annual_cv_lim <- c(0, max(annual_cv_vals, na.rm = TRUE))
  if (!all(is.finite(month_cv_lim))) month_cv_lim <- c(0, 1)
  if (!all(is.finite(annual_cv_lim))) annual_cv_lim <- c(0, 1)
  
  p_cv_annual <- base_map(rast_df(annual_cv), mp_sf) +
    fill_cv(annual_cv_lim) +
    labs(title = sprintf("Annual %s Variability Hotspots (CV %%) — Madhya Pradesh", prefix))
  print(p_cv_annual)
  
  p_cv_month <- assemble_panel(
    lapply(1:12, function(m) base_map(rast_df(month_cv[[m]]), mp_sf) + fill_cv(month_cv_lim) + labs(title = MON[m])),
    sprintf("Monthly %s Variability Hotspots (CV %%) — Madhya Pradesh", prefix),
    ncol = 4
  )
  
  terra::writeRaster(annual_cv, file.path(out_dir, "nc", sprintf("%s_annual_cv_1981_2015.tif", prefix)), overwrite = TRUE)
  
  # Monthly boxplot only
  plot_boxplot_base(
    monthly_state, "month_name", "mean_val",
    sprintf("Boxplot of Monthly %s — Madhya Pradesh", prefix),
    if (kind == "temp") "°C" else "mm/month"
  )
  
  # Monthly and annual anomaly time series
  plot_monthly_anomalies_base(
    monthly_state, "year", "anom",
    sprintf("Monthly %s Anomaly Time Series (Baseline: 1981–2010)", prefix)
  )
  
  plot_annual_anomaly_bars(
    annual_state, "year", "anom",
    sprintf("Annual %s Anomaly Time Series (Baseline: 1981–2010)", prefix),
    if (kind == "temp") "°C anomaly" else "mm/year anomaly"
  )
  
  write.csv(monthly_state, file.path(out_dir, "variability", sprintf("monthly_%s_anomaly.csv", prefix)), row.names = FALSE)
  write.csv(annual_state, file.path(out_dir, "variability", sprintf("annual_%s_anomaly.csv", prefix)), row.names = FALSE)
  
  # Export core rasters
  terra::writeRaster(annual_stack, file.path(out_dir, "nc", sprintf("%s_annual_%s_1981_2015.tif", prefix, if (kind == "prcp") "total" else "mean")), overwrite = TRUE)
  terra::writeRaster(main_ann_trend, file.path(out_dir, "nc", sprintf("%s_annual_trend_1981_2015.tif", prefix)), overwrite = TRUE)
  
  for (m in 1:12) {
    terra::writeRaster(main_month_trend[[m]], file.path(out_dir, "nc", sprintf("%s_monthly_trend_%02d.tif", prefix, m)), overwrite = TRUE)
    terra::writeRaster(month_cv[[m]], file.path(out_dir, "nc", sprintf("%s_monthly_cv_%02d.tif", prefix, m)), overwrite = TRUE)
  }
  
  invisible(list(main_stack = main_stk, rsn_stack = rsn_stk, daily = daily_df))
}

# ----------------------------------------------------------------------
# 8. Product specs
# ----------------------------------------------------------------------
temp_spec <- list(
  prefix    = "temp",
  out_dir   = cfg$temp_out_dir,
  years     = cfg$years,
  value_col = "tmean_mp",
  value_unit = "°C",
  annual_unit = "°C",
  file_fun = function(yr) file.path(cfg$temp_raw_dir, sprintf(cfg$temp_tmpl, yr))
)

prcp_spec <- list(
  prefix    = "prcp",
  out_dir   = cfg$prcp_out_dir,
  years     = cfg$years,
  value_col = "prcp_mp",
  value_unit = "mm/day",
  annual_unit = "mm/year",
  file_fun = function(yr) {
    if (yr %in% cfg$years_v1101) {
      file.path(cfg$prcp_v1101_dir, sprintf(cfg$v1101_tmpl, yr))
    } else {
      file.path(cfg$prcp_v1101ex_dir, sprintf(cfg$v1101ex_tmpl, yr))
    }
  }
)

# ----------------------------------------------------------------------
# 9. Run analyses (loads caches if available)
# ----------------------------------------------------------------------
run_product_analysis(temp_spec, cfg, kind = "temp", force_rebuild = force_rebuild)
run_product_analysis(prcp_spec, cfg, kind = "prcp", force_rebuild = force_rebuild)

message("\n=== Complete. Temperature output: ", cfg$temp_out_dir, " ===")
message("=== Complete. Precipitation output: ", cfg$prcp_out_dir, " ===")
