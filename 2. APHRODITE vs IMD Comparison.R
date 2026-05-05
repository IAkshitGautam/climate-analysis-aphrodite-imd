# ==============================================================================
# IMD + APHRODITE Validation and Comparison Script
# Madhya Pradesh | 1981–2015
#
# Uses the APHRODITE outputs:
#   Output_1-Temperature/nc/temp_monthly_MP.tif
#   Output_2-Precipitation/nc/prcp_monthly_MP.tif
#   Output_1-Temperature/timeseries/temp_daily_mean_mp.csv
#   Output_2-Precipitation/timeseries/prcp_daily_mean_mp.csv
# ==============================================================================

# ==============================================================================
# 0. PACKAGES
# ==============================================================================
pkgs <- c("terra", "trend", "ggplot2", "dplyr", "patchwork", "scales", "sf", "viridis")
new <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))
terra::terraOptions(progress = 1, threads = max(1L, parallel::detectCores(logical = FALSE) - 1L))

# ==============================================================================
# 1. USER CONFIGURATION
# ==============================================================================
base_dir <- "D:/RD539---Assignment"

cfg <- list(
  # IMD files
  imd_rain_nc = file.path(base_dir, "1. Data/IMD/imdlib_rain_1981-01-01_to_2015-12-31_polygon.nc"),
  imd_tmax_nc = file.path(base_dir, "1. Data/IMD/imdlib_tmax_1981-01-01_to_2015-12-31_polygon.nc"),
  imd_tmin_nc = file.path(base_dir, "1. Data/IMD/imdlib_tmin_1981-01-01_to_2015-12-31_polygon.nc"),
  
  # APHRODITE outputs from the merged script
  aphro_rain_stack = file.path(base_dir, "Output_2-Precipitation/nc/prcp_monthly_MP.tif"),
  aphro_temp_stack = file.path(base_dir, "Output_1-Temperature/nc/temp_monthly_MP.tif"),
  aphro_rain_daily  = file.path(base_dir, "Output_2-Precipitation/timeseries/prcp_daily_mean_mp.csv"),
  aphro_temp_daily  = file.path(base_dir, "Output_1-Temperature/timeseries/temp_daily_mean_mp.csv"),
  
  # Boundary
  shp_path = file.path(base_dir, "madhya_pradesh.shp"),
  
  # Output folder for tables only
  out_dir = file.path(base_dir, "FINAL_IMD_APHRO_Validation_Output"),
  
  years = 1981:2015,
  base_years = 1981:2010
)

for (d in c("timeseries", "stats")) {
  dir.create(file.path(cfg$out_dir, d), recursive = TRUE, showWarnings = FALSE)
}

MON <- month.abb

# ==============================================================================
# 2. SMALL HELPERS
# ==============================================================================
need_file <- function(path, label = basename(path)) {
  if (!file.exists(path)) stop("Missing required file: ", label, "\n  ", path)
}

monthly_names <- function(years) {
  sprintf("%04d_%02d", rep(years, each = 12), rep(1:12, times = length(years)))
}

month_days_year <- function(yr) {
  leap <- ((yr %% 4 == 0 & yr %% 100 != 0) | yr %% 400 == 0)
  c(31, if (leap) 29 else 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
}

theme_ts <- function() {
  theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, color = "grey30", hjust = 0.5, margin = margin(b = 10)),
      axis.title = element_text(size = 13, face = "bold"),
      axis.text  = element_text(size = 11, color = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey85", linetype = "dotted"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text  = element_text(size = 11),
      legend.position = "right"
    )
}

theme_map <- function() {
  theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 12.5, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10.5, color = "grey30", hjust = 0.5, margin = margin(b = 8)),
      axis.title = element_blank(),
      axis.text  = element_text(color = "grey40", size = 9),
      panel.grid.major = element_line(color = "grey90", linetype = "dashed"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text  = element_text(size = 10),
      legend.position = "right"
    )
}

assemble_panel <- function(plots, title, subtitle = NULL, ncol = 4) {
  patchwork::wrap_plots(plots, ncol = ncol, guides = "collect") +
    patchwork::plot_annotation(
      title = title,
      subtitle = subtitle,
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11.5, color = "grey30", hjust = 0.5),
        legend.position = "right"
      )
    )
}

show_plot <- function(p) {
  print(p)
  invisible(p)
}

rast_df <- function(r) {
  df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  names(df)[3] <- "value"
  df
}

clip_mask_rast <- function(r, mp_vect) {
  terra::mask(terra::crop(r, mp_vect, snap = "out"), mp_vect)
}

plot_map <- function(r, mp_sf, fill_scale, title = NULL, subtitle = NULL) {
  v <- terra::vect(mp_sf)
  rr <- clip_mask_rast(r, v)
  ggplot() +
    geom_raster(data = rast_df(rr), aes(x = x, y = y, fill = value)) +
    fill_scale +
    geom_sf(data = mp_sf, fill = NA, color = "black", linewidth = 0.5) +
    coord_sf(expand = FALSE) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_map()
}

fill_temp <- function(limits = NULL, name = "°C") {
  scale_fill_distiller(
    palette = "YlOrRd", direction = 1,
    limits = limits, oob = scales::squish,
    na.value = "grey90", name = name,
    guide = guide_colorbar(
      title.position = "top", title.hjust = 0.5,
      barheight = grid::unit(10.0, "cm"),
      barwidth  = grid::unit(1.2, "cm")
    )
  )
}

fill_prcp <- function(limits = NULL, name = "mm") {
  scale_fill_distiller(
    palette = "YlGnBu", direction = 1,
    limits = limits, oob = scales::squish,
    na.value = "grey90", name = name,
    guide = guide_colorbar(
      title.position = "top", title.hjust = 0.5,
      barheight = grid::unit(10.0, "cm"),
      barwidth  = grid::unit(1.2, "cm")
    )
  )
}

fill_div <- function(limits = NULL, name = "", reverse = FALSE) {
  if (reverse) {
    low_col <- "firebrick"
    high_col <- "royalblue"
  } else {
    low_col <- "royalblue"
    high_col <- "firebrick"
  }
  scale_fill_gradient2(
    low = low_col, mid = "white", high = high_col, midpoint = 0,
    limits = limits, oob = scales::squish, na.value = "grey90", name = name,
    guide = guide_colorbar(
      title.position = "top", title.hjust = 0.5,
      barheight = grid::unit(10.0, "cm"),
      barwidth  = grid::unit(1.2, "cm")
    )
  )
}

fill_pval <- function() {
  scale_fill_stepsn(
    colours = c("#1a9641", "#a6d96a", "#ffffbf", "#fdae61", "#d7191c"),
    breaks  = c(0, 0.01, 0.05, 0.10, 0.20, 1),
    limits  = c(0, 1),
    name    = "p-value",
    na.value = "grey90",
    guide = guide_colorbar(
      title.position = "top", title.hjust = 0.5,
      barheight = grid::unit(10.0, "cm"),
      barwidth  = grid::unit(1.2, "cm")
    )
  )
}

extract_dates_from_filename <- function(file) {
  b <- basename(file)
  m <- regexec("(\\d{4}-\\d{2}-\\d{2})_to_(\\d{4}-\\d{2}-\\d{2})", b)
  x <- regmatches(b, m)[[1]]
  if (length(x) == 3) seq(as.Date(x[2]), as.Date(x[3]), by = "day") else NULL
}

assign_or_check_dates <- function(r, file, expected_n = NULL) {
  d <- terra::time(r)
  if (!is.null(d) && length(d) == terra::nlyr(r) && !all(is.na(d))) return(as.Date(d))
  d <- extract_dates_from_filename(file)
  if (!is.null(d) && length(d) == terra::nlyr(r)) {
    terra::time(r) <- d
    return(d)
  }
  if (!is.null(expected_n) && terra::nlyr(r) == expected_n) {
    d <- seq(as.Date("1981-01-01"), by = "day", length.out = expected_n)
    terra::time(r) <- d
    return(d)
  }
  stop("Could not assign dates for: ", file)
}

read_nc_raster <- function(ncfile, var_pattern = NULL) {
  s <- try(terra::sds(ncfile), silent = TRUE)
  if (!inherits(s, "try-error")) {
    if (length(s) > 1) {
      nm <- names(s)
      if (!is.null(var_pattern)) {
        idx <- grep(var_pattern, nm, ignore.case = TRUE)
        if (length(idx) > 0) return(terra::rast(s[[idx[1]]]))
      }
      return(terra::rast(s[[1]]))
    }
  }
  terra::rast(ncfile)
}

select_month_layers <- function(stk, m) {
  idx <- which(substr(names(stk), 6, 7) == sprintf("%02d", m))
  if (length(idx) == 0) stop("No monthly layers found for month ", m)
  stk[[idx]]
}

daily_to_monthly_annual <- function(daily, dates, kind = c("temp", "prcp")) {
  kind <- match.arg(kind)
  year  <- as.integer(format(dates, "%Y"))
  month <- as.integer(format(dates, "%m"))
  ym <- sprintf("%04d_%02d", year, month)
  
  monthly <- terra::tapp(daily, factor(ym, levels = unique(ym)),
                         fun = if (kind == "prcp") sum else mean, na.rm = TRUE)
  names(monthly) <- unique(ym)
  
  annual <- terra::tapp(daily, factor(year, levels = sort(unique(year))),
                        fun = if (kind == "prcp") sum else mean, na.rm = TRUE)
  names(annual) <- as.character(sort(unique(year)))
  
  list(monthly = monthly, annual = annual)
}

monthly_climatology <- function(monthly_stack, years) {
  out <- lapply(1:12, function(m) {
    idx <- match(sprintf("%04d_%02d", years, m), names(monthly_stack))
    if (anyNA(idx)) stop("Missing baseline layers for month ", m)
    r <- monthly_stack[[idx[1]]]
    if (length(idx) > 1) {
      for (i in 2:length(idx)) r <- r + monthly_stack[[idx[i]]]
      r <- r / length(idx)
    }
    names(r) <- MON[m]
    r
  })
  terra::rast(out)
}

annual_climatology <- function(monthly_clim, kind = c("temp", "prcp")) {
  kind <- match.arg(kind)
  if (kind == "prcp") {
    terra::app(monthly_clim, sum, na.rm = TRUE)
  } else {
    w <- month_days_year(2001)
    r <- monthly_clim[[1]] * w[1]
    for (i in 2:12) r <- r + monthly_clim[[i]] * w[i]
    r / sum(w)
  }
}

annual_weighted_from_monthly <- function(monthly_stack, years) {
  out <- lapply(years, function(yr) {
    idx <- match(sprintf("%04d_%02d", yr, 1:12), names(monthly_stack))
    if (anyNA(idx)) stop("Missing monthly layers for year ", yr)
    w <- month_days_year(yr)
    r <- monthly_stack[[idx[1]]] * w[1]
    for (i in 2:12) r <- r + monthly_stack[[idx[i]]] * w[i]
    r <- r / sum(w)
    names(r) <- as.character(yr)
    r
  })
  terra::rast(out)
}

annual_sum_from_monthly <- function(monthly_stack, years) {
  out <- lapply(years, function(yr) {
    idx <- match(sprintf("%04d_%02d", yr, 1:12), names(monthly_stack))
    if (anyNA(idx)) stop("Missing monthly layers for year ", yr)
    r <- monthly_stack[[idx[1]]]
    for (i in 2:12) r <- r + monthly_stack[[idx[i]]]
    names(r) <- as.character(yr)
    r
  })
  terra::rast(out)
}

mk_sen_pixel <- function(x) {
  if (sum(!is.na(x)) < 5) return(c(NA_real_, NA_real_, NA_real_))
  mk <- tryCatch(trend::mk.test(x), error = function(e) NULL)
  sl <- tryCatch(trend::sens.slope(x), error = function(e) NULL)
  if (is.null(mk) || is.null(sl)) return(c(NA_real_, NA_real_, NA_real_))
  c(unname(mk$statistic), mk$p.value, unname(sl$estimates))
}

compute_trend <- function(stk) {
  r3 <- terra::app(stk, mk_sen_pixel, cores = 1)
  names(r3) <- c("mk_z", "pvalue", "slope")
  r3
}

trend_summary <- function(tr, label, alpha = 0.05) {
  pv <- terra::values(tr[["pvalue"]], na.rm = FALSE)
  sl <- terra::values(tr[["slope"]], na.rm = FALSE)
  ok <- !is.na(pv)
  n_ok <- sum(ok)
  n_sig <- sum(pv[ok] < alpha, na.rm = TRUE)
  sig_sl <- sl[ok & !is.na(sl) & pv < alpha]
  data.frame(
    Variable = label,
    N_pixels = n_ok,
    Sig_pct = round(if (n_ok == 0) NA_real_ else 100 * n_sig / n_ok, 1),
    Pos_sig_pct = round(if (n_sig == 0) NA_real_ else 100 * sum(sig_sl > 0, na.rm = TRUE) / n_sig, 1),
    Median_slope = round(median(sl[ok], na.rm = TRUE), 4),
    stringsAsFactors = FALSE
  )
}

pixel_correlation <- function(x, y) {
  n <- terra::nlyr(x)
  cxy <- c(x, y)
  terra::app(cxy, fun = function(v) {
    a <- v[1:n]
    b <- v[(n + 1):(2 * n)]
    ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < 10) return(NA_real_)
    cor(a[ok], b[ok], use = "complete.obs")
  })
}

pixel_bias <- function(x, y) {
  n <- terra::nlyr(x)
  cxy <- c(x, y)
  terra::app(cxy, fun = function(v) {
    a <- v[1:n]
    b <- v[(n + 1):(2 * n)]
    ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < 1) return(NA_real_)
    mean(a[ok] - b[ok], na.rm = TRUE)
  })
}

series_metrics <- function(obs, mod) {
  ok <- is.finite(obs) & is.finite(mod)
  if (sum(ok) < 3) {
    return(data.frame(n = sum(ok), bias = NA_real_, mae = NA_real_, rmse = NA_real_, cor = NA_real_))
  }
  data.frame(
    n = sum(ok),
    bias = mean(mod[ok] - obs[ok]),
    mae  = mean(abs(mod[ok] - obs[ok])),
    rmse = sqrt(mean((mod[ok] - obs[ok])^2)),
    cor  = cor(obs[ok], mod[ok]),
    stringsAsFactors = FALSE
  )
}

stack_state_series <- function(stk) {
  data.frame(
    name = names(stk),
    year = as.integer(substr(names(stk), 1, 4)),
    month = as.integer(substr(names(stk), 6, 7)),
    month_name = factor(MON[as.integer(substr(names(stk), 6, 7))], levels = MON),
    value = terra::global(stk, "mean", na.rm = TRUE)[, 1],
    stringsAsFactors = FALSE
  )
}

annual_state_series <- function(stk) {
  data.frame(
    year = as.integer(names(stk)),
    value = terra::global(stk, "mean", na.rm = TRUE)[, 1],
    stringsAsFactors = FALSE
  )
}

daily_state_from_stack <- function(daily_stack, mp_sf) {
  v <- terra::vect(mp_sf)
  r <- terra::mask(terra::crop(daily_stack, v), v)
  as.numeric(terra::global(r, "mean", na.rm = TRUE)[, 1])
}

load_daily_series <- function(csv_path, value_col) {
  need_file(csv_path)
  df <- read.csv(csv_path)
  if (!("date" %in% names(df))) stop("CSV missing date column: ", csv_path)
  if (!(value_col %in% names(df))) stop("CSV missing value column ", value_col, ": ", csv_path)
  df$date <- as.Date(df$date)
  df <- dplyr::arrange(df, date)
  if (!("year" %in% names(df))) df$year <- as.integer(format(df$date, "%Y"))
  if (!("month" %in% names(df))) df$month <- as.integer(format(df$date, "%m"))
  if (!("month_name" %in% names(df))) df$month_name <- factor(MON[df$month], levels = MON)
  df
}

monthly_and_annual_state <- function(df, value_col, kind = c("temp", "prcp")) {
  kind <- match.arg(kind)
  
  monthly <- df |>
    dplyr::group_by(year, month, month_name) |>
    dplyr::summarise(
      value = if (kind == "prcp") sum(.data[[value_col]], na.rm = TRUE) else mean(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    )
  
  annual <- df |>
    dplyr::group_by(year) |>
    dplyr::summarise(
      value = if (kind == "prcp") sum(.data[[value_col]], na.rm = TRUE) else mean(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    )
  
  list(monthly = monthly, annual = annual)
}

add_baseline_anomaly <- function(monthly_df, annual_df, base_years) {
  base_month <- monthly_df |>
    dplyr::filter(year %in% base_years) |>
    dplyr::group_by(month, month_name) |>
    dplyr::summarise(clim = mean(value, na.rm = TRUE), .groups = "drop")
  
  base_ann <- annual_df |>
    dplyr::filter(year %in% base_years) |>
    dplyr::summarise(clim = mean(value, na.rm = TRUE)) |>
    dplyr::pull(clim)
  
  monthly_df <- dplyr::left_join(monthly_df, base_month, by = c("month", "month_name")) |>
    dplyr::mutate(anom = value - clim)
  
  annual_df <- annual_df |>
    dplyr::mutate(clim = base_ann, anom = value - clim)
  
  list(monthly = monthly_df, annual = annual_df)
}

# ==============================================================================
# 3. LOAD BOUNDARY
# ==============================================================================
need_file(cfg$shp_path)
mp_sf <- sf::st_read(cfg$shp_path, quiet = TRUE) |> sf::st_transform(4326)
mp_vect <- terra::vect(mp_sf) |> terra::project("EPSG:4326")

# ==============================================================================
# 4. LOAD IMD DAILY DATA
# ==============================================================================
message("\n>>> Loading IMD data <<<")
need_file(cfg$imd_rain_nc)
need_file(cfg$imd_tmax_nc)
need_file(cfg$imd_tmin_nc)

imd_rain <- read_nc_raster(cfg$imd_rain_nc, "rain|precip|prcp")
imd_tmax <- read_nc_raster(cfg$imd_tmax_nc, "tmax|max")
imd_tmin <- read_nc_raster(cfg$imd_tmin_nc, "tmin|min")

imd_rain_dates <- assign_or_check_dates(imd_rain, cfg$imd_rain_nc)
imd_tmax_dates <- assign_or_check_dates(imd_tmax, cfg$imd_tmax_nc, expected_n = length(seq(as.Date("1981-01-01"), as.Date("2015-12-31"), by = "day")))
imd_tmin_dates <- assign_or_check_dates(imd_tmin, cfg$imd_tmin_nc, expected_n = length(seq(as.Date("1981-01-01"), as.Date("2015-12-31"), by = "day")))

terra::time(imd_rain) <- imd_rain_dates
terra::time(imd_tmax) <- imd_tmax_dates
terra::time(imd_tmin) <- imd_tmin_dates
imd_tmean <- (imd_tmax + imd_tmin) / 2
terra::time(imd_tmean) <- imd_tmax_dates

message("\n>>> Aggregating IMD daily data <<<")
imd_rain_ag  <- daily_to_monthly_annual(imd_rain,  imd_rain_dates, kind = "prcp")
imd_tmax_ag  <- daily_to_monthly_annual(imd_tmax,  imd_tmax_dates, kind = "temp")
imd_tmin_ag  <- daily_to_monthly_annual(imd_tmin,  imd_tmin_dates, kind = "temp")
imd_tmean_ag <- daily_to_monthly_annual(imd_tmean, imd_tmax_dates, kind = "temp")

imd_rain_monthly  <- clip_mask_rast(imd_rain_ag$monthly,  mp_vect)
imd_rain_annual   <- clip_mask_rast(imd_rain_ag$annual,   mp_vect)
imd_tmean_monthly <- clip_mask_rast(imd_tmean_ag$monthly, mp_vect)
imd_tmean_annual  <- clip_mask_rast(imd_tmean_ag$annual,  mp_vect)
imd_tmax_monthly  <- clip_mask_rast(imd_tmax_ag$monthly,  mp_vect)
imd_tmax_annual   <- clip_mask_rast(imd_tmax_ag$annual,   mp_vect)
imd_tmin_monthly  <- clip_mask_rast(imd_tmin_ag$monthly,  mp_vect)
imd_tmin_annual   <- clip_mask_rast(imd_tmin_ag$annual,   mp_vect)

# ==============================================================================
# 5. LOAD MERGED APHRODITE PRODUCTS
# ==============================================================================
message("\n>>> Loading merged APHRODITE outputs <<<")
need_file(cfg$aphro_rain_stack)
need_file(cfg$aphro_temp_stack)
need_file(cfg$aphro_rain_daily)
need_file(cfg$aphro_temp_daily)

aphro_rain <- terra::rast(cfg$aphro_rain_stack)
aphro_temp <- terra::rast(cfg$aphro_temp_stack)
aphro_rain_daily <- load_daily_series(cfg$aphro_rain_daily, "prcp_mp")
aphro_temp_daily <- load_daily_series(cfg$aphro_temp_daily, "tmean_mp")

if (terra::nlyr(aphro_rain) == length(monthly_names(cfg$years))) names(aphro_rain) <- monthly_names(cfg$years)
if (terra::nlyr(aphro_temp) == length(monthly_names(cfg$years))) names(aphro_temp) <- monthly_names(cfg$years)

imd_rain_template  <- imd_rain_annual[[1]]
imd_tmean_template <- imd_tmean_annual[[1]]
aphro_rain_al <- clip_mask_rast(terra::resample(aphro_rain, imd_rain_template, method = "bilinear"), mp_vect)
aphro_temp_al <- clip_mask_rast(terra::resample(aphro_temp, imd_tmean_template, method = "bilinear"), mp_vect)

aphro_rain_annual <- annual_sum_from_monthly(aphro_rain_al, cfg$years)
aphro_temp_annual <- annual_weighted_from_monthly(aphro_temp_al, cfg$years)

# ==============================================================================
# 6. IMD AND APHRO STATE SERIES
# ==============================================================================
message("\n>>> Building state time series <<<")
imd_rain_daily_mp  <- daily_state_from_stack(imd_rain, mp_sf)
imd_tmean_daily_mp <- daily_state_from_stack(imd_tmean, mp_sf)

imd_rain_daily_df <- data.frame(
  date = imd_rain_dates,
  year = as.integer(format(imd_rain_dates, "%Y")),
  month = as.integer(format(imd_rain_dates, "%m")),
  month_name = factor(MON[as.integer(format(imd_rain_dates, "%m"))], levels = MON),
  rain_mm = imd_rain_daily_mp
)

imd_tmean_daily_df <- data.frame(
  date = imd_tmax_dates,
  year = as.integer(format(imd_tmax_dates, "%Y")),
  month = as.integer(format(imd_tmax_dates, "%m")),
  month_name = factor(MON[as.integer(format(imd_tmax_dates, "%m"))], levels = MON),
  tmean_c = imd_tmean_daily_mp
)

imd_rain_state  <- monthly_and_annual_state(imd_rain_daily_df,  "rain_mm",  kind = "prcp")
imd_tmean_state <- monthly_and_annual_state(imd_tmean_daily_df, "tmean_c",  kind = "temp")

aphro_rain_state  <- monthly_and_annual_state(aphro_rain_daily,  "prcp_mp",  kind = "prcp")
aphro_temp_state  <- monthly_and_annual_state(aphro_temp_daily,  "tmean_mp", kind = "temp")

imd_rain_state  <- add_baseline_anomaly(imd_rain_state$monthly,  imd_rain_state$annual,  cfg$base_years)
imd_tmean_state <- add_baseline_anomaly(imd_tmean_state$monthly, imd_tmean_state$annual, cfg$base_years)
aphro_rain_state <- add_baseline_anomaly(aphro_rain_state$monthly, aphro_rain_state$annual, cfg$base_years)
aphro_temp_state <- add_baseline_anomaly(aphro_temp_state$monthly, aphro_temp_state$annual, cfg$base_years)

write.csv(imd_rain_state$annual,  file.path(cfg$out_dir, "timeseries", "imd_rain_annual_state.csv"), row.names = FALSE)
write.csv(imd_rain_state$monthly, file.path(cfg$out_dir, "timeseries", "imd_rain_monthly_state.csv"), row.names = FALSE)
write.csv(imd_tmean_state$annual, file.path(cfg$out_dir, "timeseries", "imd_tmean_annual_state.csv"), row.names = FALSE)
write.csv(imd_tmean_state$monthly,file.path(cfg$out_dir, "timeseries", "imd_tmean_monthly_state.csv"), row.names = FALSE)

write.csv(aphro_rain_state$annual,  file.path(cfg$out_dir, "timeseries", "aphro_rain_annual_state.csv"), row.names = FALSE)
write.csv(aphro_rain_state$monthly, file.path(cfg$out_dir, "timeseries", "aphro_rain_monthly_state.csv"), row.names = FALSE)
write.csv(aphro_temp_state$annual,  file.path(cfg$out_dir, "timeseries", "aphro_tmean_annual_state.csv"), row.names = FALSE)
write.csv(aphro_temp_state$monthly, file.path(cfg$out_dir, "timeseries", "aphro_tmean_monthly_state.csv"), row.names = FALSE)

# ==============================================================================
# 7. IMD CLIMATOLOGY MAPS
# ==============================================================================
message("\n>>> IMD climatology maps <<<")
imd_rain_month_clim  <- monthly_climatology(imd_rain_monthly,  cfg$base_years)
imd_rain_annual_clim <- annual_climatology(imd_rain_month_clim, "prcp")
imd_tmean_month_clim <- monthly_climatology(imd_tmean_monthly, cfg$base_years)
imd_tmean_annual_clim<- annual_climatology(imd_tmean_month_clim, "temp")

imd_tmax_month_clim  <- monthly_climatology(imd_tmax_monthly,  cfg$base_years)
imd_tmin_month_clim  <- monthly_climatology(imd_tmin_monthly,  cfg$base_years)

prcp_mon_lim  <- range(terra::values(imd_rain_month_clim,  na.rm = TRUE), na.rm = TRUE)
tmean_mon_lim <- range(terra::values(imd_tmean_month_clim, na.rm = TRUE), na.rm = TRUE)
tmax_mon_lim  <- range(terra::values(imd_tmax_month_clim,  na.rm = TRUE), na.rm = TRUE)
tmin_mon_lim  <- range(terra::values(imd_tmin_month_clim,  na.rm = TRUE), na.rm = TRUE)

p_rain_month <- assemble_panel(
  lapply(1:12, function(m) plot_map(imd_rain_month_clim[[m]], mp_sf, fill_prcp(prcp_mon_lim, "mm/mo"), MON[m])),
  title = "IMD Monthly Precipitation Normal",
  subtitle = "Baseline: 1981–2010 | Unified scale across months",
  ncol = 4
)
show_plot(p_rain_month)

p_rain_annual <- plot_map(imd_rain_annual_clim, mp_sf, fill_prcp(NULL, "mm/yr"),
                          title = "IMD Mean Annual Precipitation Normal", subtitle = "Baseline: 1981–2010")
show_plot(p_rain_annual)

p_tmean_month <- assemble_panel(
  lapply(1:12, function(m) plot_map(imd_tmean_month_clim[[m]], mp_sf, fill_temp(tmean_mon_lim, "°C"), MON[m])),
  title = "IMD Monthly Mean Temperature Normal",
  subtitle = "Baseline: 1981–2010 | Unified scale across months",
  ncol = 4
)
show_plot(p_tmean_month)

p_tmean_annual <- plot_map(imd_tmean_annual_clim, mp_sf, fill_temp(NULL, "°C"),
                           title = "IMD Mean Annual Temperature Normal", subtitle = "Baseline: 1981–2010")
show_plot(p_tmean_annual)

p_tmax_month <- assemble_panel(
  lapply(1:12, function(m) plot_map(imd_tmax_month_clim[[m]], mp_sf, fill_temp(tmax_mon_lim, "°C"), MON[m])),
  title = "IMD Monthly Tmax Normal (1981–2010)", subtitle = "Unified scale across months", ncol = 4
)
show_plot(p_tmax_month)

p_tmin_month <- assemble_panel(
  lapply(1:12, function(m) plot_map(imd_tmin_month_clim[[m]], mp_sf, fill_temp(tmin_mon_lim, "°C"), MON[m])),
  title = "IMD Monthly Tmin Normal (1981–2010)", subtitle = "Unified scale across months", ncol = 4
)
show_plot(p_tmin_month)

# ==============================================================================
# 8. IMD TREND ANALYSIS
# ==============================================================================
message("\n>>> IMD trend analysis <<<")
imd_rain_ann_trend  <- compute_trend(imd_rain_annual)
imd_tmean_ann_trend <- compute_trend(imd_tmean_annual)

imd_rain_mon_trend  <- lapply(1:12, function(m) compute_trend(select_month_layers(imd_rain_monthly,  m)))
imd_tmean_mon_trend <- lapply(1:12, function(m) compute_trend(select_month_layers(imd_tmean_monthly, m)))

all_rain_slopes  <- unlist(lapply(imd_rain_mon_trend,  function(tr) terra::values(tr[["slope"]], na.rm = TRUE)))
all_tmean_slopes <- unlist(lapply(imd_tmean_mon_trend, function(tr) terra::values(tr[["slope"]], na.rm = TRUE)))
rain_slope_lim   <- c(-quantile(abs(all_rain_slopes),  0.99, na.rm = TRUE), quantile(abs(all_rain_slopes),  0.99, na.rm = TRUE))
tmean_slope_lim  <- c(-quantile(abs(all_tmean_slopes), 0.99, na.rm = TRUE), quantile(abs(all_tmean_slopes), 0.99, na.rm = TRUE))

p_rain_trend <- {
  p1 <- plot_map(imd_rain_ann_trend[["slope"]], mp_sf, fill_div(NULL, "mm/yr", reverse = TRUE),
                 "Sen's Slope", "Precipitation: positive = wetter, negative = drier")
  p2 <- plot_map(imd_rain_ann_trend[["mk_z"]], mp_sf, fill_div(NULL, "MK Z", reverse = TRUE),
                 "Mann-Kendall Z", "Trend direction & strength")
  p3 <- plot_map(imd_rain_ann_trend[["pvalue"]], mp_sf, fill_pval(),
                 "p-value", "Statistical significance")
  (p1 | p2 | p3) + patchwork::plot_annotation(
    title = "IMD Annual Precipitation Trend",
    subtitle = "Period: 1981–2015",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 12, hjust = 0.5))
  )
}
show_plot(p_rain_trend)

p_tmean_trend <- {
  p1 <- plot_map(imd_tmean_ann_trend[["slope"]], mp_sf, fill_div(NULL, "°C/yr", reverse = FALSE),
                 "Sen's Slope", "Temperature: positive = warmer, negative = cooler")
  p2 <- plot_map(imd_tmean_ann_trend[["mk_z"]], mp_sf, fill_div(NULL, "MK Z", reverse = FALSE),
                 "Mann-Kendall Z", "Trend direction & strength")
  p3 <- plot_map(imd_tmean_ann_trend[["pvalue"]], mp_sf, fill_pval(),
                 "p-value", "Statistical significance")
  (p1 | p2 | p3) + patchwork::plot_annotation(
    title = "IMD Annual Mean Temperature Trend",
    subtitle = "Period: 1981–2015",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 12, hjust = 0.5))
  )
}
show_plot(p_tmean_trend)

p_rain_monthly_slope <- assemble_panel(
  lapply(1:12, function(m) plot_map(imd_rain_mon_trend[[m]][["slope"]], mp_sf, fill_div(rain_slope_lim, "mm/yr", reverse = TRUE), MON[m])),
  title = "IMD Monthly Precipitation Sen's Slope",
  subtitle = "Period: 1981–2015 | Unified scale across months",
  ncol = 4
)
show_plot(p_rain_monthly_slope)

p_tmean_monthly_slope <- assemble_panel(
  lapply(1:12, function(m) plot_map(imd_tmean_mon_trend[[m]][["slope"]], mp_sf, fill_div(tmean_slope_lim, "°C/yr", reverse = FALSE), MON[m])),
  title = "IMD Monthly Mean Temperature Sen's Slope",
  subtitle = "Period: 1981–2015 | Unified scale across months",
  ncol = 4
)
show_plot(p_tmean_monthly_slope)

trend_sum <- dplyr::bind_rows(
  trend_summary(imd_rain_ann_trend, "IMD Rain Annual"),
  dplyr::bind_rows(lapply(1:12, function(m) trend_summary(imd_rain_mon_trend[[m]], paste("IMD Rain", MON[m])))),
  trend_summary(imd_tmean_ann_trend, "IMD Tmean Annual"),
  dplyr::bind_rows(lapply(1:12, function(m) trend_summary(imd_tmean_mon_trend[[m]], paste("IMD Tmean", MON[m]))))
)
write.csv(trend_sum, file.path(cfg$out_dir, "stats", "imd_trend_significance_summary.csv"), row.names = FALSE)

# ==============================================================================
# 9. COMPARISON METRICS AND MAPS
# ==============================================================================
message("\n>>> APHRODITE vs IMD comparisons <<<")
rain_ann_cor   <- pixel_correlation(aphro_rain_annual, imd_rain_annual)
rain_ann_bias  <- pixel_bias(aphro_rain_annual, imd_rain_annual)
tmean_ann_cor  <- pixel_correlation(aphro_temp_annual, imd_tmean_annual)
tmean_ann_bias <- pixel_bias(aphro_temp_annual, imd_tmean_annual)

p_rain_ann_cor <- plot_map(
  rain_ann_cor, mp_sf,
  scale_fill_viridis_c(limits = c(-1, 1), oob = scales::squish, name = "Pearson r",
                       guide = guide_colorbar(barheight = grid::unit(10, "cm"), barwidth = grid::unit(1.2, "cm"))),
  title = "Annual Precipitation Correlation",
  subtitle = "APHRODITE vs IMD (1981–2015)"
)
show_plot(p_rain_ann_cor)

p_tmean_ann_cor <- plot_map(
  tmean_ann_cor, mp_sf,
  scale_fill_viridis_c(limits = c(-1, 1), oob = scales::squish, name = "Pearson r",
                       guide = guide_colorbar(barheight = grid::unit(10, "cm"), barwidth = grid::unit(1.2, "cm"))),
  title = "Annual Temperature Correlation",
  subtitle = "APHRODITE vs IMD (1981–2015)"
)
show_plot(p_tmean_ann_cor)

p_rain_ann_bias <- plot_map(
  rain_ann_bias, mp_sf, fill_div(NULL, "mm Bias", reverse = TRUE),
  title = "Annual Precipitation Bias",
  subtitle = "APHRODITE - IMD (1981–2015)"
)
show_plot(p_rain_ann_bias)

p_tmean_ann_bias <- plot_map(
  tmean_ann_bias, mp_sf, fill_div(NULL, "°C Bias", reverse = FALSE),
  title = "Annual Temperature Bias",
  subtitle = "APHRODITE - IMD (1981–2015)"
)
show_plot(p_tmean_ann_bias)

rain_month_cor <- lapply(1:12, function(m) {
  pixel_correlation(select_month_layers(aphro_rain_al, m), select_month_layers(imd_rain_monthly, m))
})
tmean_month_cor <- lapply(1:12, function(m) {
  pixel_correlation(select_month_layers(aphro_temp_al, m), select_month_layers(imd_tmean_monthly, m))
})

p_rain_month_cor <- assemble_panel(
  lapply(1:12, function(m) {
    plot_map(rain_month_cor[[m]], mp_sf,
             scale_fill_viridis_c(limits = c(-1, 1), oob = scales::squish, name = "r",
                                  guide = guide_colorbar(barheight = grid::unit(10, "cm"), barwidth = grid::unit(1.2, "cm"))),
             MON[m])
  }),
  title = "Monthly Precipitation Correlation Maps",
  subtitle = "APHRODITE vs IMD (1981–2015)",
  ncol = 4
)
show_plot(p_rain_month_cor)

p_tmean_month_cor <- assemble_panel(
  lapply(1:12, function(m) {
    plot_map(tmean_month_cor[[m]], mp_sf,
             scale_fill_viridis_c(limits = c(-1, 1), oob = scales::squish, name = "r",
                                  guide = guide_colorbar(barheight = grid::unit(10, "cm"), barwidth = grid::unit(1.2, "cm"))),
             MON[m])
  }),
  title = "Monthly Temperature Correlation Maps",
  subtitle = "APHRODITE vs IMD (1981–2015)",
  ncol = 4
)
show_plot(p_tmean_month_cor)

# ==============================================================================
# 10. STATE SCATTER PLOTS
# ==============================================================================
message("\n>>> Scatter / boxplot comparisons <<<")
imd_rain_ann_state2   <- annual_state_series(imd_rain_annual)
aphro_rain_ann_state  <- annual_state_series(aphro_rain_annual)
imd_tmean_ann_state2  <- annual_state_series(imd_tmean_annual)
aphro_tmean_ann_state <- annual_state_series(aphro_temp_annual)

scatter_annual_df <- dplyr::inner_join(
  dplyr::rename(imd_rain_ann_state2, imd_rain = value),
  dplyr::rename(aphro_rain_ann_state, aphro_rain = value),
  by = "year"
) |>
  dplyr::inner_join(
    dplyr::rename(imd_tmean_ann_state2, imd_tmean = value),
    by = "year"
  ) |>
  dplyr::inner_join(
    dplyr::rename(aphro_tmean_ann_state, aphro_tmean = value),
    by = "year"
  )

p_scatter_rain_ann <- ggplot(scatter_annual_df, aes(imd_rain, aphro_rain)) +
  geom_point(size = 3, alpha = 0.7, color = "dodgerblue4") +
  geom_smooth(method = "lm", se = TRUE, color = "darkred", fill = "pink") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = "Annual Precipitation Scatter",
       subtitle = "APHRODITE vs IMD (1981–2015) | Linear trend and 1:1 line",
       x = "IMD (mm)", y = "APHRODITE (mm)") +
  theme_ts()
show_plot(p_scatter_rain_ann)

p_scatter_tmean_ann <- ggplot(scatter_annual_df, aes(imd_tmean, aphro_tmean)) +
  geom_point(size = 3, alpha = 0.7, color = "dodgerblue4") +
  geom_smooth(method = "lm", se = TRUE, color = "darkred", fill = "pink") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = "Annual Mean Temperature Scatter",
       subtitle = "APHRODITE vs IMD (1981–2015) | Linear trend and 1:1 line",
       x = "IMD (°C)", y = "APHRODITE (°C)") +
  theme_ts()
show_plot(p_scatter_tmean_ann)

imd_rain_monthly_state2   <- imd_rain_state$monthly
aphro_rain_monthly_state2 <- aphro_rain_state$monthly
imd_tmean_monthly_state2  <- imd_tmean_state$monthly
aphro_tmean_monthly_state2<- aphro_temp_state$monthly

rain_month_scatter <- dplyr::inner_join(
  dplyr::rename(imd_rain_monthly_state2, imd = value),
  dplyr::rename(aphro_rain_monthly_state2, aphro = value),
  by = c("year", "month", "month_name")
)
tmean_month_scatter <- dplyr::inner_join(
  dplyr::rename(imd_tmean_monthly_state2, imd = value),
  dplyr::rename(aphro_tmean_monthly_state2, aphro = value),
  by = c("year", "month", "month_name")
)

p_scatter_month_rain <- ggplot(rain_month_scatter, aes(imd, aphro)) +
  geom_point(alpha = 0.6, color = "dodgerblue4", size = 1.5) +
  geom_smooth(method = "lm", se = FALSE, color = "darkred") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~month_name, scales = "free") +
  labs(title = "Monthly Precipitation Scatter", subtitle = "APHRODITE vs IMD", x = "IMD (mm)", y = "APHRODITE (mm)") +
  theme_ts()
show_plot(p_scatter_month_rain)

p_scatter_month_tmean <- ggplot(tmean_month_scatter, aes(imd, aphro)) +
  geom_point(alpha = 0.6, color = "dodgerblue4", size = 1.5) +
  geom_smooth(method = "lm", se = FALSE, color = "darkred") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~month_name, scales = "free") +
  labs(title = "Monthly Temperature Scatter", subtitle = "APHRODITE vs IMD", x = "IMD (°C)", y = "APHRODITE (°C)") +
  theme_ts()
show_plot(p_scatter_month_tmean)

# ==============================================================================
# 11. BOXPLOTS AND SUMMARY METRICS
# ==============================================================================
annual_box_rain <- dplyr::bind_rows(
  data.frame(Product = "APHRODITE", Value = aphro_rain_ann_state$value),
  data.frame(Product = "IMD",       Value = imd_rain_ann_state2$value)
)
annual_box_temp <- dplyr::bind_rows(
  data.frame(Product = "APHRODITE", Value = aphro_tmean_ann_state$value),
  data.frame(Product = "IMD",       Value = imd_tmean_ann_state2$value)
)

p_box_rain_ann <- ggplot(annual_box_rain, aes(Product, Value, fill = Product)) +
  geom_boxplot(notch = TRUE, alpha = 0.7, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.5, color = "black") +
  scale_fill_viridis_d(option = "D", guide = "none") +
  labs(title = "Annual Precipitation Comparison", subtitle = "Distribution across 1981–2015", x = NULL, y = "Total Precipitation (mm)") +
  theme_ts()
show_plot(p_box_rain_ann)

p_box_temp_ann <- ggplot(annual_box_temp, aes(Product, Value, fill = Product)) +
  geom_boxplot(notch = TRUE, alpha = 0.7, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.5, color = "black") +
  scale_fill_viridis_d(option = "D", guide = "none") +
  labs(title = "Annual Temperature Comparison", subtitle = "Distribution across 1981–2015", x = NULL, y = "Mean Temperature (°C)") +
  theme_ts()
show_plot(p_box_temp_ann)

monthly_box_rain <- dplyr::bind_rows(
  dplyr::mutate(imd_rain_monthly_state2, Product = "IMD"),
  dplyr::mutate(aphro_rain_monthly_state2, Product = "APHRODITE")
)
monthly_box_tmean <- dplyr::bind_rows(
  dplyr::mutate(imd_tmean_monthly_state2, Product = "IMD"),
  dplyr::mutate(aphro_tmean_monthly_state2, Product = "APHRODITE")
)

p_box_rain_month <- ggplot(monthly_box_rain, aes(Product, value, fill = Product)) +
  geom_boxplot(notch = TRUE, alpha = 0.7, outlier.size = 0.8) +
  scale_fill_viridis_d(option = "D", guide = "none") +
  facet_wrap(~month_name, scales = "free_y", ncol = 4) +
  labs(title = "Monthly Precipitation Comparison", subtitle = "Distribution grouped by month (1981–2015)", x = NULL, y = "Precipitation (mm)") +
  theme_ts()
show_plot(p_box_rain_month)

p_box_temp_month <- ggplot(monthly_box_tmean, aes(Product, value, fill = Product)) +
  geom_boxplot(notch = TRUE, alpha = 0.7, outlier.size = 0.8) +
  scale_fill_viridis_d(option = "D", guide = "none") +
  facet_wrap(~month_name, scales = "free_y", ncol = 4) +
  labs(title = "Monthly Temperature Comparison", subtitle = "Distribution grouped by month (1981–2015)", x = NULL, y = "Temperature (°C)") +
  theme_ts()
show_plot(p_box_temp_month)

rain_stats_ann  <- series_metrics(imd_rain_ann_state2$value,  aphro_rain_ann_state$value)
tmean_stats_ann <- series_metrics(imd_tmean_ann_state2$value, aphro_tmean_ann_state$value)

rain_stats_month <- dplyr::bind_rows(lapply(1:12, function(m) {
  obs_mod <- dplyr::inner_join(
    dplyr::filter(imd_rain_monthly_state2, month == m) |> dplyr::rename(obs = value),
    dplyr::filter(aphro_rain_monthly_state2, month == m) |> dplyr::rename(mod = value),
    by = c("year", "month", "month_name")
  )
  out <- series_metrics(obs_mod$obs, obs_mod$mod)
  out$month <- MON[m]
  out
}))

tmean_stats_month <- dplyr::bind_rows(lapply(1:12, function(m) {
  obs_mod <- dplyr::inner_join(
    dplyr::filter(imd_tmean_monthly_state2, month == m) |> dplyr::rename(obs = value),
    dplyr::filter(aphro_tmean_monthly_state2, month == m) |> dplyr::rename(mod = value),
    by = c("year", "month", "month_name")
  )
  out <- series_metrics(obs_mod$obs, obs_mod$mod)
  out$month <- MON[m]
  out
}))

stats_all <- dplyr::bind_rows(
  dplyr::mutate(rain_stats_ann,  Variable = "Rain Annual"),
  dplyr::mutate(tmean_stats_ann, Variable = "Temperature Annual")
)
write.csv(stats_all,       file.path(cfg$out_dir, "stats", "annual_comparison_stats.csv"), row.names = FALSE)
write.csv(rain_stats_month, file.path(cfg$out_dir, "stats", "monthly_rain_comparison_stats.csv"), row.names = FALSE)
write.csv(tmean_stats_month,file.path(cfg$out_dir, "stats", "monthly_tmean_comparison_stats.csv"), row.names = FALSE)

# ==============================================================================
# 12. OPTIONAL: DISPLAY COMPARISON RASTERS WITHOUT SAVING
# ==============================================================================
show_plot(
  plot_map(clip_mask_rast(rain_ann_cor, mp_vect), mp_sf,
           scale_fill_viridis_c(limits = c(-1, 1), oob = scales::squish, name = "Pearson r"),
           title = "Annual Precipitation Correlation", subtitle = "APHRODITE vs IMD (1981–2015)")
)

show_plot(
  plot_map(clip_mask_rast(tmean_ann_cor, mp_vect), mp_sf,
           scale_fill_viridis_c(limits = c(-1, 1), oob = scales::squish, name = "Pearson r"),
           title = "Annual Temperature Correlation", subtitle = "APHRODITE vs IMD (1981–2015)")
)

show_plot(
  plot_map(clip_mask_rast(rain_ann_bias, mp_vect), mp_sf, fill_div(NULL, "mm Bias", reverse = TRUE),
           title = "Annual Precipitation Bias", subtitle = "APHRODITE - IMD (1981–2015)")
)

show_plot(
  plot_map(clip_mask_rast(tmean_ann_bias, mp_vect), mp_sf, fill_div(NULL, "°C Bias", reverse = FALSE),
           title = "Annual Temperature Bias", subtitle = "APHRODITE - IMD (1981–2015)")
)

message("\n=== Complete. Tables saved to: ", cfg$out_dir, " ===")


# ==============================================================================
# 13. PRINT FINAL STATISTICAL RESULTS
# ==============================================================================
message("\n================ FINAL COMPARISON STATISTICS ================\n")

cat("\n--- Annual comparison stats ---\n")
print(stats_all, row.names = FALSE)

cat("\n--- Rainfall monthly comparison stats ---\n")
print(rain_stats_month, row.names = FALSE)

cat("\n--- Temperature monthly comparison stats ---\n")
print(tmean_stats_month, row.names = FALSE)

cat("\n--- Key annual metrics ---\n")
cat(sprintf("Rainfall  : n=%d, bias=%.3f, MAE=%.3f, RMSE=%.3f, r=%.3f\n",
            rain_stats_ann$n, rain_stats_ann$bias, rain_stats_ann$mae,
            rain_stats_ann$rmse, rain_stats_ann$cor))

cat(sprintf("Temperature: n=%d, bias=%.3f, MAE=%.3f, RMSE=%.3f, r=%.3f\n",
            tmean_stats_ann$n, tmean_stats_ann$bias, tmean_stats_ann$mae,
            tmean_stats_ann$rmse, tmean_stats_ann$cor))

message("\n================ END OF STATISTICAL RESULTS ================\n")
