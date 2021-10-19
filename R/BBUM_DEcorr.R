#' FDR correction of secondary effects and significance calling after DESeq2
#'
#' Process DESeq2 data tables. Out of two subsets of data rows, one with
#'   primary effects and one without, \code{BBUM_DEcorr} models the p values to
#'   distinguish secondary effects signal from primary effects signal using the
#'   BBUM model. It then calls significant genes within the signal set after
#'   correcting for FDR of **both** null **and** secondary effects in multiple
#'   testing.
#'
#' @param df.deseq The output DESeqResults object from DESeq2, after running
#'   results(). This object should normally be of \code{data.frame} subclass.
#'   The \code{data.frame} could also have been further processed before calling
#'   \code{BBUM_DEcorr}, as long as the \code{pvalue} column is unchanged.
#' @param classCol A \code{str} of the name of the column that indicates whether
#'   a data point belongs in the signal set or the background set.
#' @param geneName A \code{str} of the name of the column that contains the
#'   names of genes, e.g. \code{"FeatureID"}. Leave as default (\code{NULL}) if
#'   the column has not been changed (as row names, or \code{"row"} if
#'   \code{tidy=T} in \code{results()}).
#' @param pBBUM.alpha Cutoff level of \code{pBBUM} for significance calling.
#' @param excluded A \code{str} vector of all gene names that did **not** meet
#'   some user-determined prior filtering criteria (e.g. read cutoff). If
#'   \code{excluded} is not supplied, all genes with valid p values will be used
#'   for analysis.
#' @param outliers A \code{str} vector of all (additional) genes names that
#'   should be regarded as outliers among the background set.
#' @inheritParams BBUM_corr
#'
#' @details The output \code{data.frame} preserves all original rows in order,
#'   all original columns in order, and row names if present.
#' @details \code{classCol}: For example, in the case where positive log fold
#'   changes represent the signal class, and negative log fold changes represent
#'   the background class, all positive fold change points have
#'   \code{classCol == T}, and all negative fold change points have
#'   \code{classCol == F}.
#' @details \code{pBBUM} represents the expected overall FDR level if the cutoff
#'   were set at that particular p value. This is similar to the interpretation
#'   of p values corrected through the typical \code{p.adjust(method = "fdr")}.
#' @details \code{BBUM_DEcorr} functions properly only with p values filtered
#'   for poor quality data points, e.g. low read counts. Filtering has either
#'   been done manually in prior, through DESeq2's \code{independentFiltering}
#'   in \code{results()}, or through the \code{excluded} argument. Excluding a
#'   point via \code{excluded} is equivalent to filtering it out of the
#'   data.frame to begin with.
#' @details Outliers (\code{outliers}) are included in the graphs from
#'   \code{BBUM_plot()} and have p values associated. Excluded points
#'   (\code{excluded}) are removed from analysis altogether.
#' @details Default starts for BBUM fitting are implemented. If additional
#'   starts should included, or only custom starts should be considered, make
#'   use of \code{add_starts} and/or \code{only_start} arguments.
#' @details If more than one start achieved the identical maximum likelihood,
#'   A random start is chosen among them.
#' @details Automatic outlier detection relies on the model fitting a value of
#'   \code{r > 1}. For benchmarking of the trimming strategy, see
#'   Wang & Bartel, 2021.
#' @details Adding too many starts or allowing too much outlier trimming can
#'   increase computation time.
#'
#' @return A \code{data.frame} based on the input results table, with added
#'   columns from BBUM correction and significance calling:
#' * \code{geneName}: Copy of extracted gene names (\code{str})
#' * \code{excluded}: Whether the point is considered failing the cutoff and
#'   excluded in \code{excluded} (\code{bool})
#' * \code{outlier}: Whether the point is considered an outlier (\code{bool})
#' * \code{BBUM.class}: Whether the data point is considered in the singal set
#'   instead of the background set (\code{bool})
#' * \code{BBUM.l}, \code{BBUM.a}, \code{BBUM.th}, \code{BBUM.r}: Coefficients
#'   from BBUM model fit (\code{float})
#' * \code{pBBUM}: BBUM-transformed, FDR-adusted p value (\code{float})
#' * \code{BBUM.hit}: Whether a gene is a hit after signficance calling
#'   (\code{bool})
#' * \code{BBUM.fct}: Factor summarizing status as hits, non-hits, or outliers
#'   (\code{factor})
#'
#' @examples
#' \dontrun{
#' # Typical default run
#' # DESeq2
#' dds = DESeqDataSetFromMatrix(countData = cts, ...) # DESeq2
#' dds = DESeq(dds)
#' res = results(dds) %>%
#'   as.data.frame() %>%
#'   mutate(FCup = log2FoldChange > 0)
#'
#' # Call bbum function
#' res.BBUMcorr = BBUM_DEcorr(
#'   df.deseq = res,
#'   classCol = "FCup"
#'   )
#'
#' # Or with some customized behavior
#' res.BBUMcorr = BBUM_DEcorr(
#'   df.deseq = res,
#'   classCol = "FCup",
#'   geneName = "miRNAs",
#'   pBBUM.alpha = 0.01,
#'   excluded = c("hsa-miR-124-5p", "hsa-miR-1-3p"),
#'   outliers = c("hsa-miR-21-5p"),
#'   add_starts = list(c(lambda = 0.9, a = 0.6, theta = 0.1, r = 0.1))
#'   )
#' }
#'
#' @export
BBUM_DEcorr = function(
  df.deseq,
  classCol,
  geneName = NULL,
  pBBUM.alpha = 0.05,
  excluded = c(),
  outliers = c(),
  add_starts = list(), only_start = FALSE,
  limits = list(),
  auto_outliers = TRUE, rthres = 1,
  rtrimmax = 0.05, atrimmax = 10,
  quiet = FALSE
) {

  # Check, collect, and process data.frame input ----
  if(!quiet){print(">> Starting BBUM correction.")}
  df.res = tibble::as_tibble(df.deseq, rownames = NA, .name_repair = "minimal")

  # First handle whatever the name column is in the provided input, for usage convenience
  if(!is.null(geneName) && !(geneName %in% colnames(df.deseq))){
    stop("Provided geneName col not in input data.frame!")
    }
  restore.rownames = FALSE
  if(is.null(geneName)) {
    if("row" %in% colnames(df.res)){
      # The tidy=T case
      df.res = dplyr::mutate(df.res, geneName = as.character(row))
    } else {
      # Default
      restore.rownames = TRUE
      df.res = tibble::rownames_to_column(df.res, "geneName")
    }
  } else {
    # User-provided
    df.res = dplyr::mutate(df.res, geneName = as.character(!!dplyr::sym(geneName)))
  }

  # Now handle whatever the BBUM class column is in the provided input, for usage convenience
  if(!(classCol %in% colnames(df.res))){ stop("Provided class col not in input data.frame!") }
  df.res = dplyr::mutate(df.res, BBUM.class = !!dplyr::sym(classCol))

  # Process read cutoff and other filtering
  allGeneNames = as.character(dplyr::pull(df.res, geneName))
  meets.cutoff = allGeneNames[!allGeneNames %in% excluded]
  if(length(meets.cutoff) < 10){ stop("Number of genes meeting cutoff is too small!") }

  # Process outlier filter
  outliers = outliers[outliers %in% meets.cutoff]
  outliersN = length(outliers)
  if(!quiet && outliersN > 0){
    print(paste0(">> Outliers manually masked (# = ", outliersN, "):"))
    print(as.character(outliers))
  }

  # Final polishing
  df.dat = df.res %>%
    dplyr::mutate(
      excluded = !(as.character(geneName) %in% meets.cutoff),
      outlier  =   as.character(geneName) %in% outliers
    )

  # Model BBUM on p values ----
  BBUM.out = df.dat %>%
    BBUM_apply(
      add_starts = add_starts,
      only_start = only_start,
      limits = limits,
      auto_outliers = auto_outliers,
      rthres = rthres,
      rtrimmax = rtrimmax,
      atrimmax = atrimmax,
      quiet = quiet
    )

  # Call hits ----
  df.bbum = BBUM.out$data %>%
    dplyr::mutate(
      BBUM.hits = !outlier &
        dplyr::if_else(is.na(pBBUM), FALSE, pBBUM < pBBUM.alpha) &  # pBBUM criterion
        BBUM.class, # Call primary effect hits
      BBUM.fct = factor(
        dplyr::if_else(
          BBUM.hits,
          "hit",
          dplyr::if_else(outlier,
                         "outlier",
                         "none"
          )
        ), levels = c("none","hit","outlier"))  # summaritive factor for categories
    ) %>%
    dplyr::ungroup()

  # Benchmarking ----
  if(!quiet){
    # Print hits list
    hits = df.bbum %>%
      dplyr::filter(BBUM.hits) %>%
      dplyr::arrange(pBBUM) %>%
      dplyr::pull(geneName) %>%
      unname()
    c_hits = length(hits)
    print(paste0(">> Significant genes (# = ", c_hits, "):"))
    if(c_hits > 0){
      print(as.character(hits))
    }

    # Check hit rate and in comparison with theoretical hit rate from the fitted distribution
    c_N_pos = nrow(df.bbum %>% dplyr::filter(!excluded, !outlier, BBUM.class))
      # number of points in sample class
    # Hit counts calculations
    c_hits.empi = c_hits * (1-pBBUM.alpha)
    c_hits.theo = BBUM_expectedhits(
      BBUM.out,
      pBBUM.alpha,
      n = c_N_pos
    )
    print(paste0(">> Empirical true-hit count estimate: ",
                 round(c_hits.empi , 1),
                 ", theoretical true-hit count estimate: ",
                 round(c_hits.theo, 1)
    ))
  }

  # Restore how the input df was ----
  if(restore.rownames){
    df.bbum = df.bbum %>%
      dplyr::mutate(geneName.rowname = geneName) %>%
      tibble::column_to_rownames("geneName.rowname")
  }

  # Done ----
  return(df.bbum)
}
