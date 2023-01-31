#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(vroom))

parse_args <- function(args) {
  # create parser object
  parser <- ArgumentParser(prog='Define sweep regions using smooth splines')
  parser$add_argument('--do_tests', action='store_true', help='Perform units tests')
  subparsers = parser$add_subparsers()
  
  parser_neutral = subparsers$add_parser("neutral_cutoff", 
  help="Define cutoff using null/neutral data and quantile value.")
  
  ####################
  ## neutral_cutoff ##
  ####################
  
  parser_neutral$add_argument(
    "-q",
    "--quantile",
    type = "double",
    help = "A quantile probabilitiy to define outliers."
  )
  parser_neutral$add_argument(
    "-f",
    "--data_frame",
    type = "character",
    help = "Name of the file to treat the input data from."
  )
  parser_neutral$add_argument(
    "-R",
    "--input_raisd2",
    type = "logical",
    default = FALSE,
    help = "Option to treat input as 2 column file from RAiSD. Generated when RAiSD input is from the ms simulator"
  )
  parser_neutral$add_argument(
    "-d",
    "--delimeter",
    type = "character",
    default = ",",
    help = "Characeter used to seperate fields on input."
  )
  parser_neutral$add_argument(
    "-s",
    "--skip_rows",
    type = "integer",
    default = 0,
    help = "Number of rows to skip when reading input."
  )
  parser_neutral$add_argument(
    "-b",
    "--skip_character",
    type = "character",
    default = "#",
    help = "Skips past lines that starts with this character."
  )
  parser_neutral$add_argument(
    "-p",
    "--positions",
    type = "character",
    help = "Name of column to use as positions"
  )
  parser_neutral$add_argument(
    "-v",
    "--values",
    type = "character",
    help = "Name of column to use as values of interest. Must match length of positions argument. Missing values are dropped."
  )
  
  ###################
  ## sweep_regions ##
  ###################
  parser_regions = subparsers$add_parser('sweep_regions', help='Define sweep regions based on various user defined cutoffs')
  
  parser_regions$add_argument(
    "-c",
    "--chromosome",
    type = "character",
    default = "chrom_NAN",
    help = "Chromosome data is taken from. Passed value is overridden if --input_raisd8 TRUE."
  )
  parser_regions$add_argument(
    "-f",
    "--data_frame",
    type = "character",
    help = "Name of the file to treat the input data from."
  )
  parser_regions$add_argument(
    "-L",
    "--truth_data",
    type = "character",
    help = "Tab delimeted file containing columns labeled: 'position' -- integer postions along input chromosome; and 'type' -- either 'sweep' or 'neutral' depending on state at the position"
  )
  parser_regions$add_argument(
    "-R",
    "--input_raisd8",
    type = "logical",
    default = FALSE,
    help = "Option to treat input as 8 column file from RAiSD."
  )
  parser_regions$add_argument(
    "-d",
    "--delimeter",
    type = "character",
    default = ",",
    help = "Characeter used to seperate fields on input."
  )
  parser_regions$add_argument(
    "-s",
    "--skip_rows",
    type = "integer",
    default = 0,
    help = "Number of rows to skip when reading input."
  )
  parser_regions$add_argument(
    "-b",
    "--skip_character",
    type = "character",
    default = "#",
    help = "Skips past lines that starts with this character."
  )
  parser_regions$add_argument(
    "-p",
    "--positions",
    type = "character",
    help = "Name of column to use as positions"
  )
  parser_regions$add_argument(
    "-v",
    "--values",
    type = "character",
    help = "Name of column to use as values of interest. Must match length of positions argument. Rows with missing values are dropped."
  )
  parser_regions$add_argument(
    "-S",
    "--merge_size",
    type = "double",
    help = "The distance between identifed regions to determine if they should be treated as a single region or distinct."
  )
  parser_regions$add_argument(
    "-C",
    "--cutoff",
    type = "double",
    help = "The cutoff value used to define regions. Ideally chosen using neutral simulations."
  )
  parser_regions$add_argument(
    "-m",
    "--min_size",
    type = "double",
    help = "The minium size an outlier region can be classified as a sweep"
  )
  parser_regions$add_argument(
    "-M",
    "--max_size",
    type = "double",
    help = "The maximum size an outlier region can be classified as a sweep"
  )
  parser_regions$add_argument(
    "-F",
    "--full_output",
    type = "logical",
    default = FALSE,
    help = "Option to write the entire output, rather than just one line per region. Good for plotting and troubleshooting."
  )
  parser_regions$add_argument(
    "-o",
    "--out_file",
    help = "Name of the outfile to write the merged regions to."
  )
  parser$parse_args()
}

### FUNCTIONS ###
make_sweep_regions <- function(chromosome,
                               positions,
                               mu_vector,
                               minimum_size,
                               maximum_size,
                               cutoff,
                               merge_size,
                               full_dataframe = FALSE) {
  if (class(positions) != "numeric") {
    stop("postions must be of class numeric")
  }
  if (class(mu_vector) != "numeric") {
    stop("mu_vector must be of class numeric")
  }
  if (length(positions) != length(mu_vector)) {
    stop("positions and mu_vector inputs must be the same length")
  }
  
  #fit spline, filter values above cutoff, merge into regions, keep regions of desired size
  spline_df <- smooth.spline(positions, mu_vector) %>%
    do(data.frame(positions = .$x, mu_vector = .$y), .) %>%
    filter(mu_vector > cutoff)
  
  if (nrow(spline_df) > 1) {
    spline_df <- spline_df %>%
      mutate(
        outlier_dists = c(0, positions[2:n()] - positions[1:(n() - 1)]),
        sweep_group = cumsum(outlier_dists > merge_size)
      ) %>%
      group_by(sweep_group) %>%
      mutate(
        chromosome = chromosome,
        start = min(positions),
        end = max(positions),
        sweep_size = end - start
      ) %>%
      filter(between(sweep_size, minimum_size, maximum_size)) %>%
      ungroup()
    
    if (!full_dataframe & nrow(spline_df) != 0) {
      spline_df <- spline_df %>%
        select(chromosome, sweep_group, start, end, sweep_size) %>%
        distinct() %>%
        ungroup()
    }
    return(spline_df)
  } else {
    tibble()
  }
}

position_checks <- function(true_positions, true_type) {
  if (mean(unique(true_type) %in% c("sweep", "neutral")) != 1) {
    stop("elements in type must be 'sweep' or 'neutral' only.")
  }
  
  if (sum(table(true_positions) > 1) != 0) {
    stop("true_positions must be unique")
  }
  
  return(TRUE)
}

#if true sweep positions are available, label inferred sweeps and TP/FP
label_TF_positives <-
  function(sweep_df, true_positions, true_type) {
    position_checks(true_positions, true_type)
    
    true_positions <- true_positions[true_type == "sweep"]
    true_type <- true_type[true_type == "sweep"]
    
    sweep_df %>%
      rowwise() %>%
      mutate(true_covered = sum(between(true_positions, start, end))) %>%
      mutate(
        truth_class = case_when(
          true_covered < 1 ~ "false_positive",
          true_covered >= 1 ~ "true_positive"
        )
      ) %>%
      ungroup()
  }

label_TF_negatives <-
  function(sweep_df, true_positions, true_type) {
    if (nrow(sweep_df) < 1) {
      sweep_df <- tibble(chromosome = NA,
                         start = NA,
                         end = NA)
    }
    
    chromosome <- sweep_df$chromosome[1]
    position_checks(true_positions, true_type)
    negative_df <- purrr::map2_dfr(true_positions, true_type,
                                   function(positions, type) {
                                     found_at <-
                                       which(sweep_df$start <= positions &
                                               sweep_df$end >= positions)
                                     out_tibble <-
                                       tibble(
                                         chromosome = chromosome,
                                         start = positions,
                                         end = positions,
                                         sweep_size = 0,
                                         type = type,
                                         truth_class = NA
                                       )
                                     if (length(found_at) == 0 &
                                         type == 'sweep') {
                                       return(mutate(out_tibble, truth_class = 'false_negative'))
                                     }
                                     if (length(found_at) == 0 &
                                         type == 'neutral') {
                                       return(mutate(out_tibble, truth_class = 'true_negative'))
                                     }
                                   })
    bind_rows(sweep_df, negative_df)
  }



testing_suite <- function() {
  #test setup
  #make up some example data using a GP
  set.seed(1334)
  n_points <- 500
  positions <- seq(1, n_points, length.out = n_points)
  gamma <- 10
  sq_exp <- as.matrix(exp(-dist(positions) ^ 2 / (2 * gamma ^ 2)))
  pred_mu <-
    suppressWarnings(mvtnorm::rmvnorm(n = 1, mean = rep(0, n_points), sq_exp))
  pred_mu <- pred_mu - min(pred_mu)
  mu_obs <- rnorm(n_points, pred_mu, 1)
  chr <- "chr1"
  min_s <- 10
  max_s <- 1e10
  merge_s <- 10
  cut <- 3
  
  out_full <- make_sweep_regions(
    chromosome = chr,
    positions = positions,
    mu_vector = mu_obs,
    minimum_size = min_s,
    maximum_size = max_s,
    merge_size = merge_s,
    cutoff = cut,
    full_dataframe = TRUE
  )
  
  outt <- make_sweep_regions(
    chromosome = chr,
    positions = positions,
    mu_vector = mu_obs,
    minimum_size = min_s,
    maximum_size = max_s,
    merge_size = merge_s,
    cutoff = cut,
    full_dataframe = FALSE
  )
  
  test_that("true_positive", {
    true_positions <- c(112, 485)
    true_labels <- c("sweep", "sweep")
    result_df <- outt %>%
      label_TF_positives(true_positions, true_labels) %>%
      label_TF_negatives(true_positions, true_labels)
    prop_tp <- mean(result_df$truth_class == "true_positive")
    expect_equal(prop_tp, 1)
  })
  
  test_that("false_negative", {
    true_positions <- c(10)
    true_labels <- c("sweep")
    result_df <- outt %>%
      label_TF_positives(true_positions, true_labels) %>%
      label_TF_negatives(true_positions, true_labels)
    fn <- filter(result_df, start == 10) %>% pull(truth_class)
    expect_equal(fn, "false_negative")
  })
  
  test_that("bad_positions", {
    bad_f <-
      function()
        make_sweep_regions(chromosome = "chr1", positions = c("a"))
    expect_error(bad_f(), "postions must be of class numeric")
  })
  
  test_that("bad_lengths", {
    bad_f <-
      function()
        make_sweep_regions(
          chromosome = "chr1",
          positions = (1),
          mu_vector = c(1, 2)
        )
    expect_error(bad_f(),
                 "positions and mu_vector inputs must be the same length")
  })
  
}

### RUN CODE! ###
main <- function() {
  cli_args = commandArgs(trailingOnly = TRUE)
  args <- parse_args(cli_args)
  
  #do tests or not
  if (args$do_tests) {
    library(testthat)
    testing_suite()
    do_tests <- TRUE
    return("Unit tests done.")
  }
  
  if(cli_args[1] == "sweep_regions"){
    #get the inout columns
    if (args$input_raisd8) {
      colnames <- c("chromosome", "mid", "start", "end", "mu1", "mu2","mu3", "mu") 
      sweep_data <- vroom::vroom(
        file = args$data_frame,
        col_names = colnames,
        delim = args$delimeter,
        skip = args$skip_rows,
        comment = args$skip_character
      ) %>% drop_na()
      positions_vector <- sweep_data$mid
      values_vector <- sweep_data$mu
      chrom <-  sweep_data$chromosome[1]
    } else {
      sweep_data <- vroom::vroom(
        file = args$data_frame,
        delim = args$delimeter,
        skip = args$skip_rows,
        comment = args$skip_character
      ) %>% drop_na()
      chrom <- args$chromosome
      positions_vector <- sweep_data[[args$positions]]
      values_vector <- sweep_data[[args$values]]
    }
    
    out_df <- make_sweep_regions(
      chromosome = chrom,
      positions =  positions_vector,
      mu_vector =  values_vector,
      minimum_size = args$min_size,
      maximum_size = args$max_size,
      merge_size = args$merge_size,
      cutoff = args$cutoff,
      full_dataframe = args$full_output
    )
    
    
    if(!args$full_out && !is.null(args$truth_data)){
      true_df <- read.table(args$truth_data, header = TRUE)
      out_df <- out_df %>%
        label_TF_positives(true_df$positions, true_df$type) %>%
        label_TF_negatives(true_df$positions, true_df$type)
    }
    
    
      #bedfile output
      bedname <- paste0(strsplit(args$out_file, "\\.")[[1]][1], ".bed")
      if(nrow(out_df > 0)){
        out_df %>% 
          mutate(start = start-1) %>% 
          select(chromosome, start, end) %>% 
          write.table(
            file = bedname, 
            sep = "\t", 
            quote = FALSE, 
            col.names = FALSE, 
            row.names = FALSE
          )  
      } else {
        warning("No sweep regions found.")
        out_df %>% 
        write.table(
          file = bedname, 
          sep = "\t", 
          quote = FALSE, 
          col.names = FALSE, 
          row.names = FALSE
        )
      }
      
      #redundant full output
      write.csv(out_df, file = args$out_file)  
    
  }
  
  if(cli_args[1] == "neutral_cutoff"){
    if (args$input_raisd2) {
      colnames <- c("position", "mu") 
      neutral_data <- vroom::vroom(
        file = args$data_frame,
        col_names = colnames,
        delim = args$delimeter,
        skip = 1,
        comment = args$skip_character
      ) %>% drop_na()
      positions_vector <- neutral_data$position
      values_vector <- neutral_data$mu
    } else {
      neutral_data <- vroom::vroom(
        file = args$data_frame,
        delim = args$delimeter,
        skip = args$skip_rows,
        comment = args$skip_character
      ) %>% drop_na()
      positions_vector <- neutral_data[[args$positions]]
      values_vector <- neutral_data[[args$values]]
    }
    
    smooth_model <- smooth.spline(positions_vector, values_vector)
    q_out <- quantile(smooth_model$y, args$quantile)
    cat(paste0("neutral cutoff: ", round(q_out, 5), "\n"))
    
  }
  
}

if (!interactive()) {
  main()
}
