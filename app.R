# =============================================================================
#  TrypTrack (T_curv v2) - growth curves and doubling times for T. brucei cultures
#
#  Workflow
#    1. Experiment design -> generates a form to fill in inside the app
#    2. Data              -> fill in the form, or upload a CSV saved earlier
#    3. Plots             -> cumulative growth, daily segments, doubling time,
#                            growth relative to a control condition
#    4. Doubling time     -> log-linear fit over a chosen time window
#
#  Data format (one row per count):
#    cell_line             RNAi line / genotype (e.g. "MM")
#    condition             e.g. tet-, tet+
#    replicate             biological replicate or clone (e.g. MM1, MM2, MM3)
#    day                   day of the experiment (0, 1, 2, ...)
#    clock                 clock time of counting, e.g. 11:00 or 16:30
#                          (blank = same clock time as the first count)
#    count                 average count from the counting device
#                          (density = count x multiplier, default 10^4 cells/ml)
#    dilution_factor       dilution made right after this count
#                          (1 or blank = not diluted)
#    count_after_dilution  optional; only if you counted again after diluting.
#                          If given, it is used instead of dilution_factor.
#    time_h                optional; hours since the first count. Only used
#                          for rows without a day (e.g. older files).
#
#  Elapsed time = 24 x (day - first day) + (clock - clock of first count),
#  so counting at 11:00 one day and 16:00 the next gives a 29 h interval.
#
#  Old T_curv files (cell_line, condition, hours, cells_diluted, cells_24h)
#  are recognised and converted automatically.
# =============================================================================

library(shiny)
library(bslib)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rhandsontable)
library(zip)

COLS     <- c("cell_line", "condition", "replicate", "day", "clock", "count",
              "dilution_factor", "count_after_dilution", "time_h")
CHR_COLS <- c("cell_line", "condition", "replicate", "clock")
NUM_COLS <- c("day", "count", "dilution_factor", "count_after_dilution", "time_h")
APP_NAME <- "TrypTrack"       # change the app name here

CLOCK_CHOICES <- sprintf("%02d:%02d", rep(0:23, each = 2), c(0, 30))

# ----------------------------------------------------------------------------
#  Helpers
# ----------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

parse_list <- function(x) {
  if (is.null(x) || !nzchar(trimws(x))) return(character(0))
  v <- trimws(strsplit(x, ",")[[1]])
  v[nzchar(v)]
}

parse_num_list <- function(x) {
  v <- suppressWarnings(as.numeric(parse_list(x)))
  sort(unique(v[!is.na(v)]))
}

# Clock time -> hours after midnight. Accepts 11:00, 11.30, 11h30, 11,
# 11:00:00, 4pm, 4:30 pm. Returns NA for blank or unreadable values.
parse_clock <- function(x) {
  s <- tolower(gsub("\\s", "", as.character(x)))
  s[is.na(s)] <- ""
  pm <- grepl("pm$", s)
  am <- grepl("am$", s)
  s  <- sub("(am|pm)$", "", s)
  ok <- grepl("^[0-9]{1,2}([.:h][0-9]{1,2}(:[0-9]{1,2})?)?$", s)
  h  <- suppressWarnings(as.numeric(sub("^([0-9]{1,2}).*$", "\\1", s)))
  m  <- ifelse(grepl("[.:h]", s),
               suppressWarnings(as.numeric(sub("^[0-9]{1,2}[.:h]([0-9]{1,2}).*$", "\\1", s))),
               0)
  h  <- ifelse(pm & h < 12, h + 12, ifelse(am & h == 12, 0, h))
  ifelse(ok & !is.na(h) & h < 24 & !is.na(m) & m < 60, h + m / 60, NA_real_)
}

clock_invalid <- function(x) {
  x <- as.character(x)
  !is.na(x) & nzchar(trimws(x)) & is.na(parse_clock(x))
}

# Make sure a data frame has the expected columns and types
normalize_df <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  for (cc in setdiff(COLS, names(df))) df[[cc]] <- NA
  df <- df[, COLS]
  for (cc in CHR_COLS) {
    v <- as.character(df[[cc]])
    v[is.na(v)] <- ""
    df[[cc]] <- trimws(v)
  }
  for (cc in NUM_COLS) df[[cc]] <- suppressWarnings(as.numeric(df[[cc]]))
  df$replicate[df$replicate == ""] <- "1"
  rownames(df) <- NULL
  df
}

# Elapsed hours since the first count of each culture
add_elapsed_time <- function(df) {
  df %>%
    mutate(.clock_h = parse_clock(clock),
           .order = if_else(is.na(day), time_h / 24, day)) %>%
    group_by(cell_line, condition, replicate) %>%
    arrange(.order, .by_group = TRUE) %>%
    mutate(
      .ref   = if (any(!is.na(.clock_h))) .clock_h[!is.na(.clock_h)][1] else 0,
      .ch    = if_else(is.na(.clock_h), .ref, .clock_h),
      .day0  = suppressWarnings(min(day, na.rm = TRUE)),
      .ch0   = .ch[!is.na(day)][1],
      elapsed_h = if_else(!is.na(day),
                          24 * (day - .day0) + (.ch - .ch0),
                          time_h)
    ) %>%
    ungroup() %>%
    select(-starts_with("."))
}

# Old T_curv format: one row per 24 h interval, with cells_diluted (start of
# interval) and cells_24h (count at 'hours')
convert_old <- function(df) {
  df %>%
    mutate(across(c(hours, cells_diluted, cells_24h), as.numeric)) %>%
    group_by(cell_line, condition) %>%
    arrange(hours, .by_group = TRUE) %>%
    mutate(count_after_dilution = lead(cells_diluted)) %>%
    ungroup() %>%
    transmute(cell_line, condition, replicate = "1",
              day = hours / 24, clock = "", count = cells_24h,
              dilution_factor = round(cells_24h / count_after_dilution, 3),
              count_after_dilution) %>%
    normalize_df()
}

read_any_csv <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                 strip.white = TRUE, fileEncoding = "UTF-8-BOM")
  names(df) <- trimws(tolower(names(df)))
  if (all(c("cells_diluted", "cells_24h", "hours") %in% names(df))) {
    return(list(data = convert_old(df), old = TRUE))
  }
  if (!"time_h" %in% names(df) && "hours" %in% names(df))
    names(df)[names(df) == "hours"] <- "time_h"
  missing <- setdiff(c("cell_line", "condition", "count"), names(df))
  if (!any(c("day", "time_h") %in% names(df))) missing <- c(missing, "day")
  if (length(missing))
    stop("Missing column(s): ", paste(missing, collapse = ", "))
  list(data = normalize_df(df), old = FALSE)
}

# Empty form from an experiment design
make_template <- function(lines, conds, reps, days, dfs, clock = "", start = NA) {
  g <- expand.grid(replicate = reps, condition = conds, cell_line = lines,
                   day = days, stringsAsFactors = FALSE)
  g <- g[order(g$day, match(g$cell_line, lines),
               match(g$condition, conds), match(g$replicate, reps)), ]
  g$clock <- clock
  g$count <- NA_real_
  if (!is.na(start)) g$count[g$day == min(days)] <- start
  g$dilution_factor <- unname(dfs[g$condition])
  g$dilution_factor[g$day == min(days)] <- 1      # first count: set-up
  g$dilution_factor[g$day == max(days)] <- NA     # last count: no dilution
  g$count_after_dilution <- NA_real_
  normalize_df(g)
}

# Core calculation ------------------------------------------------------------
#   start_after = density right after dilution
#               = count_after_dilution, or count / dilution_factor, or count
#   growth_i    = count_i / start_after_(i-1)
#   cumulative  = first count x product of growth factors
#   doubling time of an interval = dt * ln2 / ln(growth)
process_data <- function(df, mult = 1e4) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df <- normalize_df(df) %>%
    filter(nzchar(cell_line), nzchar(condition), !is.na(day) | !is.na(time_h))
  if (nrow(df) == 0) return(NULL)
  df %>%
    add_elapsed_time() %>%
    mutate(time_h = elapsed_h) %>%
    select(-elapsed_h) %>%
    filter(!is.na(time_h)) %>%
    group_by(cell_line, condition, replicate) %>%
    arrange(time_h, .by_group = TRUE) %>%
    mutate(
      start_after = case_when(
        !is.na(count_after_dilution) ~ count_after_dilution,
        !is.na(dilution_factor) & dilution_factor > 0 ~ count / dilution_factor,
        TRUE ~ count),
      interval_h        = time_h - lag(time_h),
      growth_factor     = count / lag(start_after),
      cum_count         = first(count) *
                            cumprod(if_else(row_number() == 1, 1, growth_factor)),
      doubling_time_h   = if_else(growth_factor > 1,
                                  interval_h * log(2) / log(growth_factor),
                                  NA_real_),
      doublings_per_24h = log2(growth_factor) * 24 / interval_h,
      point             = row_number()
    ) %>%
    ungroup() %>%
    mutate(density_per_ml      = count * mult,
           cum_density_per_ml  = cum_count * mult,
           after_dilution_per_ml = start_after * mult)
}

check_data <- function(df) {
  msgs <- character(0)
  if (is.null(df) || nrow(df) == 0) return(msgs)
  d <- normalize_df(df)
  where <- function(x) sprintf("%s / %s / %s, day %s", x$cell_line, x$condition,
                               x$replicate, ifelse(is.na(x$day), "?", x$day))
  badclk <- d[clock_invalid(d$clock), ]
  if (nrow(badclk))
    msgs <- c(msgs, sprintf("Clock time not understood ('%s'): %s - use e.g. 11:00",
                            badclk$clock, where(badclk)))
  nodate <- d %>% filter(nzchar(cell_line), is.na(day), is.na(time_h))
  if (nrow(nodate)) msgs <- c(msgs, sprintf("%d row(s) without a day are ignored.", nrow(nodate)))
  d <- d %>% filter(nzchar(cell_line), !is.na(day) | !is.na(time_h)) %>% add_elapsed_time()
  dup <- d %>% count(cell_line, condition, replicate, elapsed_h) %>% filter(n > 1)
  if (nrow(dup))
    msgs <- c(msgs, sprintf("Two counts at the same time (%s h): %s / %s / %s",
                            round(dup$elapsed_h, 2), dup$cell_line, dup$condition, dup$replicate))
  gaps <- d %>%
    group_by(cell_line, condition, replicate) %>% arrange(elapsed_h, .by_group = TRUE) %>%
    mutate(later = rev(cumsum(rev(!is.na(count)))) - !is.na(count)) %>%
    filter(is.na(count), later > 0) %>% ungroup()
  if (nrow(gaps))
    msgs <- c(msgs, sprintf("Missing count before later counts (curve stops here): %s",
                            where(gaps)))
  bad <- d %>% filter(!is.na(count_after_dilution), !is.na(count),
                      count_after_dilution > count)
  if (nrow(bad))
    msgs <- c(msgs, sprintf("count_after_dilution > count: %s", where(bad)))
  if (any(d$count <= 0, na.rm = TRUE)) msgs <- c(msgs, "Some counts are zero or negative.")
  if (any(d$dilution_factor < 1, na.rm = TRUE))
    msgs <- c(msgs, "Some dilution factors are < 1 (a 10x dilution should be entered as 10).")
  msgs
}

# Mean / SD over replicates, matched by time point index
summarise_reps <- function(d, var, by = c("cell_line", "condition", "point")) {
  d %>% filter(is.finite(.data[[var]])) %>%
    group_by(across(all_of(by))) %>%
    summarise(time_h = mean(time_h), mean = mean(.data[[var]]),
              sd = sd(.data[[var]]), n = n(), .groups = "drop")
}

fit_doubling <- function(t, y) {
  ok <- is.finite(t) & is.finite(y) & y > 0
  if (sum(ok) < 2) return(c(NA, NA, sum(ok)))
  f <- lm(log2(y[ok]) ~ t[ok])
  slope <- unname(coef(f)[2])                      # doublings per hour
  r2 <- if (sum(ok) > 2) summary(f)$r.squared else NA
  c(slope, r2, sum(ok))
}

okabe <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
           "#0072B2", "#D55E00", "#CC79A7")

cond_colors <- function(levels, scheme) {
  n <- length(levels)
  base <- switch(scheme,
    "Black/Red"   = c("black", "red", "#0072B2", "#009E73", "#E69F00",
                      "#CC79A7", "#56B4E9", "#D55E00"),
    "Blue/Orange" = c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
                      "#56B4E9", "black", "#F0E442"),
    "Okabe-Ito"   = okabe,
    "Viridis"     = scales::viridis_pal(end = 0.85)(max(n, 2)))
  setNames(rep(base, length.out = n), levels)
}

log_labels <- scales::trans_format("log10", scales::math_format(10^.x))

# Example data: one RNAi line (MM), three independent clones, tet-/tet+.
# Only dilution factors are given; densities after dilution are calculated.
example_data <- local({
  old <- data.frame(
    clone = rep(c("MM1", "MM1", "MM2", "MM2", "MM3", "MM3"), 4),
    condition = rep(c("tet-", "tet+"), 12),
    hours = rep(c(0, 24, 48, 72), each = 6),
    cells_diluted = c(1, 1, 1, 1, 1, 1,
                      13.3, 13.3, 11.44285714, 11.44285714, 12.4, 12.4,
                      11.82857143, 12.42857143, 10.25, 10.2, 9.25, 9.85,
                      16.3, 21.5, 11.35714286, 26, 10.37142857, 10.85714286),
    cells_24h = c(13.3, 13.3, 11.44285714, 11.44285714, 12.4, 12.4,
                  138, 145, 102.5, 102, 92.5, 98.5,
                  163, 107.5, 53, 26.5, 60.5, 38,
                  164.5, 47.5, 90, 25.5, 112, 20.5),
    stringsAsFactors = FALSE)
  old$cell_line <- old$clone
  d <- convert_old(old)
  d$replicate <- d$cell_line
  d$cell_line <- "MM"
  d$dilution_factor <- round(d$dilution_factor, 2)
  d$count_after_dilution <- NA_real_
  d$clock <- "10:00"
  d[order(d$day, d$condition, d$replicate), ] %>% normalize_df()
})

# ----------------------------------------------------------------------------
#  UI
# ----------------------------------------------------------------------------

ui <- navbarPage(
  title = tags$img(src = "logo.svg", height = "46px", alt = APP_NAME,
                   style = paste("background: #fff; border-radius: 8px;",
                                 "padding: 2px 8px; margin: -10px 0;")),
  windowTitle = APP_NAME,
  header = tags$head(tags$link(rel = "icon", type = "image/svg+xml", href = "tryp.svg")),
  id = "nav",
  theme = bs_theme(version = 5, bootswatch = "flatly"),

  # --- 1. Experiment design ---------------------------------------------------
  tabPanel("1. Experiment design",
    sidebarLayout(
      sidebarPanel(width = 4,
        numericInput("n_lines", "Number of cell lines", 1, min = 1, max = 20),
        textInput("line_names", "Cell line names (comma separated, optional)",
                  placeholder = "e.g. MM, POLx, WT"),
        textInput("conditions", "Conditions (comma separated)", "tet-, tet+"),
        numericInput("n_reps", "Biological replicates / clones per line",
                     3, min = 1, max = 20),
        textInput("rep_names", "Replicate names (comma separated, optional)",
                  placeholder = "e.g. MM1, MM2, MM3"),
        hr(),
        radioButtons("schedule", "Counting and dilution",
                     c("Every day" = "1", "Every 2 days" = "2",
                       "Custom days" = "custom")),
        conditionalPanel("input.schedule != 'custom'",
          numericInput("n_days", "Length of experiment (days)", 4, min = 1, max = 60)),
        conditionalPanel("input.schedule == 'custom'",
          textInput("custom_days", "Count days (comma separated)",
                    "0, 1, 2, 3, 6, 7")),
        textInput("default_clock", "Usual counting time", "10:00"),
        helpText("Pre-filled in the form. On days you counted at another time, ",
                 "just change the clock time of those rows (e.g. 16:00) - the ",
                 "longer or shorter interval is taken into account."),
        uiOutput("df_inputs"),
        numericInput("start_count", "Starting count on day 0 (optional, same for all)",
                     NA, min = 0),
        hr(),
        actionButton("to_form", "Create the form",
                     class = "btn-primary", icon = icon("table")),
        helpText("Replaces the table in '2. Data'. Save the table there first ",
                 "if it holds data you want to keep.")
      ),
      mainPanel(width = 8,
        h4("Form preview"),
        textOutput("template_info"),
        br(),
        tableOutput("template_preview")
      )
    )
  ),

  # --- 2. Data ----------------------------------------------------------------
  tabPanel("2. Data",
    sidebarLayout(
      sidebarPanel(width = 3,
        h5("Counting time per day"),
        helpText("Change the time here to set it for all cultures counted ",
                 "that day. Single rows can still be changed in the table."),
        rHandsontableOutput("day_clock"),
        hr(),
        fileInput("file", "Upload a saved CSV", accept = ".csv"),
        helpText("A table saved from this app, or an old T_curv file."),
        downloadButton("dl_data", "Save table as CSV"),
        helpText("Save after each count and upload the file the next day ",
                 "to continue - everything stays editable."),
        hr(),
        h5("Add a count day"),
        fluidRow(
          column(6, numericInput("new_day", "Day", NA, min = 0)),
          column(6, textInput("new_clock", "Clock", placeholder = "10:00"))
        ),
        helpText("Empty day = the day after the last one. Empty clock = same ",
                 "time as the last day."),
        actionButton("add_day", "Add day for all cultures", icon = icon("plus")),
        br(),
        actionButton("load_example", "Load example data", icon = icon("flask"),
                     class = "btn-sm btn-outline-secondary"),
        hr(),
        numericInput("mult", "Cells/ml = count x", 1e4, min = 1),
        hr(),
        helpText(HTML(
          "<b>day</b>: day of the experiment<br>",
          "<b>clock</b>: time of counting, e.g. 11:00 or 16:30 ",
          "(empty = same as the first count)<br>",
          "<b>count</b>: average count from the counting device<br>",
          "<b>dilution_factor</b>: dilution made right after this count ",
          "(1 or empty = not diluted; last day empty)<br>",
          "<b>count_after_dilution</b>: optional, only if you counted again ",
          "after diluting<br>",
          "Right-click the table to add or remove rows."))
      ),
      mainPanel(width = 9,
        uiOutput("data_messages"),
        rHandsontableOutput("hot", height = "600px"),
        br(),
        h5("Hours since first count (calculated from day + clock)"),
        tableOutput("time_preview")
      )
    )
  ),

  # --- 3. Plots ---------------------------------------------------------------
  tabPanel("3. Plots",
    sidebarLayout(
      sidebarPanel(width = 3,
        radioButtons("mode", "Replicates",
                     c("Mean +/- SD" = "mean", "Separate (one panel each)" = "sep")),
        conditionalPanel("input.mode == 'mean'",
          checkboxInput("show_reps", "Show individual replicates", TRUE)),
        uiOutput("line_select"),
        uiOutput("cond_select"),
        uiOutput("rep_select"),
        radioButtons("time_unit", "Time axis", c("Hours" = "h", "Days" = "d"),
                     inline = TRUE),
        selectInput("color_scheme", "Colours",
                    c("Black/Red", "Blue/Orange", "Okabe-Ito", "Viridis")),
        hr(),
        conditionalPanel("input.plot_tab == 'Daily segments'",
          checkboxInput("seg_log", "Log scale", FALSE)),
        conditionalPanel("input.plot_tab == 'Growth rate per interval'",
          radioButtons("rate_metric", "Show",
                       c("Doubling time (h)" = "dt",
                         "Doublings per 24 h" = "dpd"))),
        conditionalPanel("input.plot_tab == 'Relative growth'",
          uiOutput("ref_select")),
        hr(),
        fluidRow(
          column(4, numericInput("plot_width", "Width (in)", 8, min = 3, max = 30)),
          column(4, numericInput("plot_height", "Height (in)", 5, min = 3, max = 30)),
          column(4, selectInput("plot_format", "Format", c("pdf", "png")))
        ),
        downloadButton("dl_plot", "Download this plot"),
        br(), br(),
        downloadButton("dl_all", "Download all plots + tables (zip)")
      ),
      mainPanel(width = 9,
        tabsetPanel(id = "plot_tab",
          tabPanel("Cumulative growth", plotOutput("p_cum", height = "550px")),
          tabPanel("Daily segments", plotOutput("p_seg", height = "550px")),
          tabPanel("Growth rate per interval", plotOutput("p_rate", height = "550px")),
          tabPanel("Relative growth", plotOutput("p_rel", height = "550px"))
        )
      )
    )
  ),

  # --- 4. Doubling time -------------------------------------------------------
  tabPanel("4. Doubling time",
    sidebarLayout(
      sidebarPanel(width = 3,
        uiOutput("fit_window_ui"),
        helpText("Doubling time from a straight-line fit of log2(cumulative ",
                 "cells) vs time inside this window. Use e.g. the window ",
                 "before and after the RNAi phenotype starts."),
        hr(),
        downloadButton("dl_fit", "Download fit table (CSV)"),
        br(), br(),
        downloadButton("dl_processed", "Download processed data (CSV)")
      ),
      mainPanel(width = 9,
        h4("Mean over replicates"),
        tableOutput("fit_summary"),
        h4("Per replicate"),
        tableOutput("fit_table"),
        h4("Per interval"),
        tableOutput("interval_table")
      )
    )
  ),

  tabPanel("Help",
    fluidRow(column(8, offset = 1,
      h3("How it works"),
      p("Cultures are counted and diluted repeatedly. For each count the app ",
        "works out the density right after dilution (count / dilution factor, ",
        "or the count after dilution if you measured it). Growth in an interval ",
        "is the next count divided by that density."),
      tags$ul(
        tags$li("Cumulative cells = first count x product of all growth factors ",
                "(the curve you would see without dilution)."),
        tags$li("Doubling time of an interval = interval length x ln 2 / ln(growth)."),
        tags$li("Doublings per 24 h = log2(growth) x 24 / interval length."),
        tags$li("Relative growth = cumulative cells of a condition / cumulative ",
                "cells of the reference condition, same line and replicate."),
        tags$li("Mean +/- SD: replicates are matched by time point number, so ",
                "slightly different count times are averaged.")
      ),
      p("Times: each count has a day and a clock time. Elapsed time is ",
        "24 x (day - first day) + (clock time - clock time of the first count), ",
        "so counting at 11:00 on day 1 and 16:00 on day 2 gives a 29 h interval, ",
        "and growth in that interval is scaled to its real length. The ",
        "calculated hours are shown under the table in '2. Data'."),
      p("Intervals do not have to be equal: a 2-day interval with a 100x ",
        "dilution is handled the same way as a 1-day interval with 10x."),
      p("Day-by-day use: create the form, fill in today's counts, click ",
        "'Save table as CSV'. Next time, upload that file in '2. Data' and ",
        "continue. 'Add a count day' extends the experiment; any cell can ",
        "be corrected later, and rows can be added or removed by right-click."),
      p("Old T_curv CSV files (cells_diluted / cells_24h) are converted ",
        "automatically when uploaded.")
    ))
  )
)

# ----------------------------------------------------------------------------
#  Server
# ----------------------------------------------------------------------------

server <- function(input, output, session) {

  # ===== 1. Design ============================================================

  conds_design <- reactive({
    v <- parse_list(input$conditions)
    if (!length(v)) v <- "control"
    unique(v)
  })

  output$df_inputs <- renderUI({
    def <- if (input$schedule == "2") 100 else 10
    tagList(
      tags$label("Default dilution factor per condition"),
      lapply(seq_along(conds_design()), function(i) {
        numericInput(paste0("dfac_", i), conds_design()[i], def, min = 1)
      }),
      helpText("Editable for every single culture and day in the form.")
    )
  })

  design_template <- reactive({
    n_l <- max(1, as.integer(input$n_lines %||% 1))
    n_r <- max(1, as.integer(input$n_reps %||% 1))
    lines <- parse_list(input$line_names)
    if (length(lines) < n_l)
      lines <- c(lines, paste0("Line", (length(lines) + 1):n_l))
    lines <- lines[seq_len(n_l)]
    reps <- parse_list(input$rep_names)
    if (length(reps) < n_r)
      reps <- c(reps, paste0("R", (length(reps) + 1):n_r))
    reps <- reps[seq_len(n_r)]
    conds <- conds_design()
    days <- if (input$schedule == "custom") {
      parse_num_list(input$custom_days)
    } else {
      seq(0, max(1, input$n_days %||% 1), by = as.numeric(input$schedule))
    }
    validate(need(length(days) >= 2, "Need at least two count days."))
    clk <- trimws(input$default_clock %||% "")
    validate(need(!clock_invalid(clk), "Usual counting time not understood - use e.g. 10:00."))
    dfs <- sapply(seq_along(conds), function(i) {
      v <- input[[paste0("dfac_", i)]]
      if (is.null(v) || is.na(v)) 10 else v
    })
    names(dfs) <- conds
    make_template(unique(lines), conds, unique(reps), days, dfs, clock = clk,
                  start = input$start_count %||% NA)
  })

  output$template_info <- renderText({
    t <- design_template()
    sprintf("%d rows: %d cell line(s) x %d condition(s) x %d replicate(s) x %d count days (day %s).",
            nrow(t), n_distinct(t$cell_line), n_distinct(t$condition),
            n_distinct(t$replicate), n_distinct(t$day),
            paste(sort(unique(t$day)), collapse = ", "))
  })

  output$template_preview <- renderTable({
    head(design_template() %>% select(-time_h), 40)
  }, na = "", digits = 2)

  # ===== 2. Data ==============================================================

  rv <- reactiveVal(example_data)        # current data table
  table_version <- reactiveVal(0)        # bump to redraw the editable table

  set_data <- function(df) {
    rv(normalize_df(df))
    table_version(table_version() + 1)
  }

  observeEvent(input$to_form, {
    set_data(design_template())
    updateNavbarPage(session, "nav", selected = "2. Data")
    showNotification("Form created. Enter counts in the 'count' column and change clock times where needed.",
                     type = "message", duration = 8)
  })

  observeEvent(input$load_example, {
    set_data(example_data)
    showNotification("Example data loaded (line MM, clones MM1-3).", type = "message")
  })

  observeEvent(input$file, {
    res <- tryCatch(read_any_csv(input$file$datapath),
                    error = function(e) {
                      showNotification(paste("Could not read file:", conditionMessage(e)),
                                       type = "error", duration = 10)
                      NULL
                    })
    req(res)
    set_data(res$data)
    if (res$old)
      showNotification(paste("Old T_curv format converted. All rows got replicate '1';",
                             "edit cell_line / replicate if the lines are clones of one line."),
                       type = "warning", duration = 12)
  })

  output$hot <- renderRHandsontable({
    table_version()
    d <- isolate(rv())
    if (all(is.na(d$time_h))) d$time_h <- NULL     # only shown for hour-based files
    rhandsontable(d, rowHeaders = NULL, stretchH = "all", height = 600,
                  useTypes = TRUE) %>%
      hot_cols(columnSorting = FALSE) %>%
      hot_col(intersect(NUM_COLS, names(d)), format = "0.[000]") %>%
      hot_col("clock", type = "autocomplete", source = CLOCK_CHOICES,
              strict = FALSE) %>%
      hot_context_menu(allowRowEdit = TRUE, allowColEdit = FALSE)
  })

  # --- counting time per day ---
  day_clock_df <- reactive({
    d <- rv()
    d <- d[!is.na(d$day), ]
    if (!nrow(d)) return(NULL)
    d %>%
      group_by(day) %>%
      summarise(clock = names(sort(table(clock), decreasing = TRUE))[1],
                note  = if (n_distinct(clock) > 1) "rows differ" else "",
                .groups = "drop") %>%
      arrange(day) %>%
      as.data.frame()
  })

  output$day_clock <- renderRHandsontable({
    dc <- day_clock_df()
    req(dc)
    rhandsontable(dc, rowHeaders = NULL, stretchH = "all", useTypes = TRUE) %>%
      hot_col("day", readOnly = TRUE, format = "0.[0]") %>%
      hot_col("note", readOnly = TRUE) %>%
      hot_col("clock", type = "autocomplete", source = CLOCK_CHOICES,
              strict = FALSE) %>%
      hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE)
  })

  observeEvent(input$day_clock, {
    new <- tryCatch(hot_to_r(input$day_clock), error = function(e) NULL)
    old <- isolate(day_clock_df())
    req(new, old, nrow(new) == nrow(old))
    new$clock <- trimws(as.character(new$clock))
    new$clock[is.na(new$clock)] <- ""
    changed <- which(new$clock != old$clock)
    if (!length(changed)) return()
    bad <- clock_invalid(new$clock[changed])
    if (any(bad))
      showNotification("Clock time not understood - use e.g. 11:00 or 16:30.",
                       type = "warning")
    d <- rv()
    for (i in changed[!bad]) d$clock[!is.na(d$day) & d$day == old$day[i]] <- new$clock[i]
    set_data(d)
  })

  observeEvent(input$hot, {
    d <- tryCatch(hot_to_r(input$hot), error = function(e) NULL)
    if (!is.null(d)) rv(normalize_df(d))
  })

  output$data_messages <- renderUI({
    m <- check_data(rv())
    if (!length(m)) return(NULL)
    div(class = "alert alert-warning",
        tags$b("Check your data:"), tags$ul(lapply(head(m, 15), tags$li)))
  })

  output$time_preview <- renderTable({
    p <- processed()
    req(p)
    p %>%
      arrange(day, time_h) %>%
      mutate(col = if_else(is.na(day), paste0("count ", point), paste0("day ", day))) %>%
      select(cell_line, condition, replicate, col, time_h) %>%
      pivot_wider(names_from = col, values_from = time_h,
                  values_fn = function(x) x[1]) %>%
      arrange(cell_line, condition, replicate)
  }, na = "", digits = 1)

  observeEvent(input$add_day, {
    d <- rv()
    keys <- c("cell_line", "condition", "replicate")
    cultures <- d %>% filter(nzchar(cell_line), nzchar(condition)) %>% distinct(across(all_of(keys)))
    if (!nrow(cultures)) {
      showNotification("No cultures in the table yet - create the form first.", type = "warning")
      return()
    }
    days <- d$day[!is.na(d$day)]
    nd <- if (is.na(input$new_day)) (if (length(days)) max(days) + 1 else 0) else input$new_day
    clk <- trimws(input$new_clock %||% "")
    if (!nzchar(clk) && length(days)) {
      last <- d$clock[!is.na(d$day) & d$day == max(days)]
      clk <- names(sort(table(last), decreasing = TRUE))[1]
    }
    if (clock_invalid(clk)) {
      showNotification("Clock time not understood - use e.g. 10:00.", type = "error")
      return()
    }
    have <- d %>% filter(!is.na(day), day == nd) %>% distinct(across(all_of(keys)))
    new <- cultures %>% anti_join(have, by = keys)
    if (!nrow(new)) {
      showNotification(sprintf("Day %s already exists for all cultures.", nd), type = "warning")
      return()
    }
    # The previous count of each culture was the last one so far and may have
    # no dilution factor - give it the culture's last used dilution factor.
    prev <- d %>% filter(!is.na(day), day < nd) %>%
      group_by(across(all_of(keys))) %>%
      summarise(prev_day = max(day),
                df_used = {
                  v <- dilution_factor[order(day)]
                  v <- v[!is.na(v) & v > 1]
                  if (length(v)) v[length(v)] else NA_real_
                }, .groups = "drop")
    d <- d %>% left_join(prev, by = keys)
    fill <- !is.na(d$day) & !is.na(d$prev_day) & d$day == d$prev_day &
      is.na(d$dilution_factor) & !is.na(d$df_used)
    d$dilution_factor[fill] <- d$df_used[fill]
    d$prev_day <- NULL
    d$df_used <- NULL
    new <- new %>% mutate(day = nd, clock = clk)
    d <- bind_rows(d, new) %>% normalize_df()
    d <- d[order(is.na(d$day), d$day), ]          # stable: keeps culture order
    set_data(d)
    updateNumericInput(session, "new_day", value = NA)
    msg <- sprintf("Day %s added for %d culture(s).", nd, nrow(new))
    if (any(fill)) msg <- paste(msg, sprintf(
      "Dilution factor of the previous day filled in for %d row(s) - please check.", sum(fill)))
    showNotification(msg, type = "message", duration = 8)
  })

  output$dl_data <- downloadHandler(
    filename = function() paste0("T_curv_data_", Sys.Date(), ".csv"),
    content = function(file) {
      d <- rv()
      if (all(is.na(d$time_h))) d$time_h <- NULL
      write.csv(d, file, row.names = FALSE, na = "")
    }
  )

  # ===== Processing ===========================================================

  processed <- reactive({
    mult <- if (is.na(input$mult) || input$mult <= 0) 1 else input$mult
    process_data(rv(), mult)
  })

  cond_levels <- reactive({
    d <- rv()
    unique(d$condition[nzchar(d$condition)])
  })

  output$line_select <- renderUI({
    req(processed())
    ch <- unique(processed()$cell_line)
    checkboxGroupInput("lines", "Cell lines", ch, selected = ch)
  })
  output$cond_select <- renderUI({
    req(processed())
    ch <- cond_levels()
    checkboxGroupInput("conds", "Conditions", ch, selected = ch, inline = TRUE)
  })
  output$rep_select <- renderUI({
    req(processed())
    ch <- unique(processed()$replicate)
    checkboxGroupInput("reps", "Replicates", ch, selected = ch, inline = TRUE)
  })
  output$ref_select <- renderUI({
    selectInput("ref_cond", "Reference condition (= 100 %)", cond_levels())
  })

  for (o in c("line_select", "cond_select", "rep_select", "ref_select"))
    outputOptions(output, o, suspendWhenHidden = FALSE)

  dsel <- reactive({
    req(processed(), input$lines, input$conds, input$reps)
    processed() %>%
      filter(cell_line %in% input$lines, condition %in% input$conds,
             replicate %in% input$reps) %>%
      mutate(condition = factor(condition, levels = cond_levels()))
  })

  tdiv   <- reactive(if (input$time_unit == "d") 24 else 1)
  xlabel <- reactive(if (input$time_unit == "d") "Time (days)" else "Time (hours)")
  cols   <- reactive(cond_colors(cond_levels(), input$color_scheme))

  facet_layer <- function() {
    if (input$mode == "mean") facet_wrap(~cell_line)
    else facet_wrap(~cell_line + replicate)
  }

  base_theme <- function() {
    list(scale_colour_manual(values = cols()),
         scale_fill_manual(values = cols()),
         theme_bw(base_size = 13),
         theme(strip.background = element_rect(fill = "grey92"),
               legend.position = "right"))
  }

  # ===== Plot: cumulative growth ==============================================

  plot_cum <- reactive({
    d <- dsel() %>% filter(is.finite(cum_density_per_ml))
    validate(need(nrow(d) > 0, "No complete data to plot yet."))
    td <- tdiv()
    ylab <- "Cumulative cells/ml (corrected for dilution)"
    if (input$mode == "mean") {
      s <- summarise_reps(d, "cum_density_per_ml")
      w <- 0.02 * diff(range(s$time_h / td)) + 1e-9
      p <- ggplot(s, aes(time_h / td, mean, colour = condition))
      if (isTRUE(input$show_reps))
        p <- p + geom_line(data = d, aes(time_h / td, cum_density_per_ml,
                                         group = interaction(condition, replicate)),
                           alpha = 0.25, linewidth = 0.4) +
          geom_point(data = d, aes(time_h / td, cum_density_per_ml),
                     alpha = 0.3, size = 1.3)
      p <- p + geom_line(linewidth = 0.8) + geom_point(size = 2.3) +
        geom_errorbar(aes(ymin = pmax(mean - sd, mean * 0.05), ymax = mean + sd),
                      width = w, na.rm = TRUE)
    } else {
      p <- ggplot(d, aes(time_h / td, cum_density_per_ml, colour = condition)) +
        geom_line(linewidth = 0.7) + geom_point(size = 2)
    }
    p + facet_layer() + scale_y_log10(labels = log_labels) +
      labs(x = xlabel(), y = ylab, colour = NULL,
           title = "Growth corrected for dilution",
           subtitle = if (input$mode == "mean") "mean +/- SD of replicates" else NULL) +
      base_theme()
  })
  output$p_cum <- renderPlot(plot_cum())

  # ===== Plot: daily segments (real densities) ================================

  plot_seg <- reactive({
    d <- dsel()
    td <- tdiv()
    if (input$mode == "mean") {
      grp <- c("cell_line", "condition")
      pts <- d %>% filter(is.finite(density_per_ml)) %>%
        group_by(cell_line, condition, point) %>%
        summarise(time_h = mean(time_h),
                  count = mean(density_per_ml),
                  start = mean(after_dilution_per_ml), .groups = "drop")
    } else {
      grp <- c("cell_line", "condition", "replicate")
      pts <- d %>% filter(is.finite(density_per_ml)) %>%
        transmute(cell_line, condition, replicate, point, time_h,
                  count = density_per_ml, start = after_dilution_per_ml)
    }
    validate(need(nrow(pts) > 0, "No counts to plot yet."))
    pts <- pts %>% group_by(across(all_of(grp))) %>% arrange(point, .by_group = TRUE)
    grow <- pts %>%
      mutate(x = lag(time_h), y = lag(start), xend = time_h, yend = count,
             step = "growth", gap = point - lag(point)) %>%
      filter(!is.na(x), gap == 1) %>% ungroup()
    dil <- pts %>% filter(is.finite(start), abs(start - count) > 1e-9 * count) %>%
      mutate(x = time_h, y = count, xend = time_h, yend = start,
             step = "dilution") %>% ungroup()
    segs <- bind_rows(grow, dil)
    after <- pts %>% filter(is.finite(start), abs(start - count) > 1e-9 * count) %>% ungroup()

    p <- ggplot() +
      geom_segment(data = segs,
                   aes(x = x / td, y = y, xend = xend / td, yend = yend,
                       colour = condition, linetype = step), linewidth = 0.6) +
      geom_point(data = ungroup(pts), aes(time_h / td, count, colour = condition),
                 size = 2) +
      geom_point(data = after, aes(time_h / td, start, colour = condition),
                 size = 1.8, shape = 1) +
      scale_linetype_manual(values = c(growth = "solid", dilution = "dashed"),
                            name = NULL) +
      facet_layer() +
      labs(x = xlabel(), y = "Cells/ml", colour = NULL,
           title = "Measured densities and dilutions",
           subtitle = "filled = counted, open = after dilution (calculated or counted)") +
      base_theme()
    if (isTRUE(input$seg_log)) p + scale_y_log10(labels = log_labels)
    else p + scale_y_continuous(labels = scales::label_scientific())
  })
  output$p_seg <- renderPlot(plot_seg())

  # ===== Plot: growth rate per interval =======================================

  plot_rate <- reactive({
    var <- if (input$rate_metric == "dt") "doubling_time_h" else "doublings_per_24h"
    ylab <- if (var == "doubling_time_h") "Doubling time (h)" else "Doublings per 24 h"
    d <- dsel() %>% filter(point > 1, is.finite(.data[[var]]))
    validate(need(nrow(d) > 0, "No complete intervals yet."))
    td <- tdiv(); unit <- if (td == 24) "d" else "h"
    # label intervals by mean start/end time of each interval number
    labs_df <- dsel() %>% group_by(point) %>% summarise(t = mean(time_h), .groups = "drop") %>%
      arrange(point) %>%
      mutate(label = paste0(signif(lag(t) / td, 3), "-", signif(t / td, 3), " ", unit)) %>%
      filter(point > 1)
    d <- d %>% left_join(labs_df %>% select(point, label), by = "point") %>%
      mutate(label = factor(label, levels = unique(labs_df$label)))
    dodge <- position_dodge(width = 0.8)
    if (input$mode == "mean") {
      s <- d %>% group_by(cell_line, condition, label) %>%
        summarise(mean = mean(.data[[var]]), sd = sd(.data[[var]]), n = n(),
                  .groups = "drop")
      p <- ggplot(s, aes(label, mean, fill = condition)) +
        geom_col(position = dodge, width = 0.75, alpha = 0.85) +
        geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                      position = dodge, width = 0.25, na.rm = TRUE)
      if (isTRUE(input$show_reps))
        p <- p + geom_point(data = d, aes(label, .data[[var]], fill = condition),
                            shape = 21, colour = "grey20",
                            position = position_jitterdodge(jitter.width = 0.1,
                                                            dodge.width = 0.8),
                            size = 1.6, show.legend = FALSE)
    } else {
      p <- ggplot(d, aes(label, .data[[var]], fill = condition)) +
        geom_col(position = dodge, width = 0.75, alpha = 0.85)
    }
    cap <- if (var == "doubling_time_h")
      "Intervals without growth have no doubling time and are not shown." else NULL
    p + facet_layer() +
      labs(x = "Interval", y = ylab, fill = NULL, title = "Growth rate per interval",
           caption = cap) +
      base_theme() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
  })
  output$p_rate <- renderPlot(plot_rate())

  # ===== Plot: relative growth ================================================

  plot_rel <- reactive({
    ref <- input$ref_cond
    validate(need(!is.null(ref) && length(cond_levels()) > 1,
                  "Needs at least two conditions (e.g. tet- and tet+)."))
    d <- dsel()
    validate(need(ref %in% d$condition, "Reference condition is not selected."))
    r <- d %>% filter(is.finite(cum_count)) %>%
      group_by(cell_line, replicate, point) %>%
      mutate(ref_val = cum_count[condition == ref][1]) %>%
      ungroup() %>%
      filter(condition != ref, is.finite(ref_val)) %>%
      mutate(rel = 100 * cum_count / ref_val)
    validate(need(nrow(r) > 0, "No matching reference data."))
    td <- tdiv()
    if (input$mode == "mean") {
      s <- summarise_reps(r, "rel")
      w <- 0.02 * diff(range(s$time_h / td)) + 1e-9
      p <- ggplot(s, aes(time_h / td, mean, colour = condition))
      if (isTRUE(input$show_reps))
        p <- p + geom_line(data = r, aes(time_h / td, rel,
                                         group = interaction(condition, replicate)),
                           alpha = 0.25, linewidth = 0.4)
      p <- p + geom_line(linewidth = 0.8) + geom_point(size = 2.3) +
        geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = w, na.rm = TRUE)
    } else {
      p <- ggplot(r, aes(time_h / td, rel, colour = condition)) +
        geom_line(linewidth = 0.7) + geom_point(size = 2)
    }
    p + geom_hline(yintercept = 100, linetype = "dotted") +
      facet_layer() +
      labs(x = xlabel(), y = paste0("Cumulative growth (% of ", ref, ")"),
           colour = NULL, title = paste("Growth relative to", ref)) +
      base_theme()
  })
  output$p_rel <- renderPlot(plot_rel())

  # ===== 4. Doubling time fit ================================================

  output$fit_window_ui <- renderUI({
    req(processed())
    rng <- range(processed()$time_h, na.rm = TRUE)
    sliderInput("fit_window", "Fit window (hours)", min = rng[1], max = rng[2],
                value = rng, step = 1)
  })

  fit_tab <- reactive({
    req(processed(), input$fit_window)
    w <- input$fit_window
    processed() %>%
      filter(time_h >= w[1] - 1e-9, time_h <= w[2] + 1e-9) %>%
      group_by(cell_line, condition, replicate) %>%
      summarise(f = list(fit_doubling(time_h, cum_count)), .groups = "drop") %>%
      mutate(slope = sapply(f, `[`, 1),
             r_squared = sapply(f, `[`, 2),
             n_points = sapply(f, `[`, 3),
             doubling_time_h = if_else(slope > 0, 1 / slope, NA_real_),
             doublings_per_24h = 24 * slope) %>%
      select(-f, -slope) %>%
      mutate(condition = factor(condition, levels = cond_levels())) %>%
      arrange(cell_line, condition, replicate) %>%
      mutate(condition = as.character(condition))
  })

  fit_summary <- reactive({
    fit_tab() %>%
      mutate(condition = factor(condition, levels = cond_levels())) %>%
      group_by(cell_line, condition) %>%
      summarise(n_replicates = sum(is.finite(doublings_per_24h)),
                doubling_time_h_mean = mean(doubling_time_h, na.rm = TRUE),
                doubling_time_h_sd = sd(doubling_time_h, na.rm = TRUE),
                doublings_per_24h_mean = mean(doublings_per_24h, na.rm = TRUE),
                doublings_per_24h_sd = sd(doublings_per_24h, na.rm = TRUE),
                .groups = "drop") %>%
      mutate(condition = as.character(condition))
  })

  interval_tab <- reactive({
    req(processed())
    processed() %>% filter(point > 1) %>%
      transmute(cell_line, condition, replicate,
                from_h = time_h - interval_h, to_h = time_h,
                growth_factor, doubling_time_h, doublings_per_24h)
  })

  output$fit_summary <- renderTable(fit_summary(), na = "-", digits = 2)
  output$fit_table   <- renderTable(fit_tab(), na = "-", digits = 2)
  output$interval_table <- renderTable(interval_tab(), na = "-", digits = 2)

  output$dl_fit <- downloadHandler(
    filename = function() paste0("T_curv_doubling_time_", Sys.Date(), ".csv"),
    content = function(file) {
      write.csv(bind_rows(
        fit_tab() %>% mutate(level = "replicate", .before = 1),
        fit_summary() %>% mutate(level = "mean", .before = 1)),
        file, row.names = FALSE, na = "")
    }
  )
  output$dl_processed <- downloadHandler(
    filename = function() paste0("T_curv_processed_", Sys.Date(), ".csv"),
    content = function(file) write.csv(processed(), file, row.names = FALSE, na = "")
  )

  # ===== Plot downloads =======================================================

  current_plot <- reactive({
    switch(input$plot_tab,
           "Cumulative growth" = plot_cum(),
           "Daily segments" = plot_seg(),
           "Growth rate per interval" = plot_rate(),
           "Relative growth" = plot_rel())
  })

  save_plot <- function(file, p, fmt) {
    ggsave(file, p, width = input$plot_width, height = input$plot_height,
           dpi = 300, device = fmt)
  }

  output$dl_plot <- downloadHandler(
    filename = function() paste0("T_curv_", gsub(" ", "_", tolower(input$plot_tab)),
                                 ".", input$plot_format),
    content = function(file) save_plot(file, current_plot(), input$plot_format)
  )

  output$dl_all <- downloadHandler(
    filename = function() paste0("T_curv_results_", Sys.Date(), ".zip"),
    content = function(file) {
      dir <- file.path(tempdir(), paste0("tcurv_", as.integer(Sys.time())))
      dir.create(dir, showWarnings = FALSE)
      fmt <- input$plot_format
      plots <- list(cumulative_growth = function() plot_cum(),
                    daily_segments = function() plot_seg(),
                    growth_rate = function() plot_rate(),
                    relative_growth = function() plot_rel())
      files <- character(0)
      for (nm in names(plots)) {
        p <- tryCatch(plots[[nm]](), error = function(e) NULL)
        if (!is.null(p)) {
          f <- file.path(dir, paste0(nm, ".", fmt))
          ok <- tryCatch({ save_plot(f, p, fmt); TRUE }, error = function(e) FALSE)
          if (ok) files <- c(files, f)
        }
      }
      tabs <- list(data = rv(), processed = processed(),
                   doubling_time_fit = fit_tab(), doubling_time_mean = fit_summary(),
                   intervals = interval_tab())
      for (nm in names(tabs)) {
        f <- file.path(dir, paste0(nm, ".csv"))
        write.csv(tabs[[nm]], f, row.names = FALSE, na = "")
        files <- c(files, f)
      }
      zip::zip(file, files = basename(files), root = dir)
    }
  )
}

shinyApp(ui, server)
