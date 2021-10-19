#' Apply BBUM fitting and correction onto a data.frame
#'
#' Wrapper/intermediate function between \code{BBUM_corr} and
#'   \code{BBUM_DEcorr}. For internal use of the package.
#'
#' @param df \code{data.frame} input
#' @inheritParams BBUM_corr
#'
BBUM_apply = function(
  df,
  add_starts, only_start,
  limits,
  auto_outliers, rthres,
  rtrimmax, atrimmax,
  quiet
) {

  # Check input columns ----
  if (!all(c(
    "geneName", "BBUM.class", "excluded", "outlier", "pvalue"
  ) %in% colnames(df))) {
    stop("Missing cols in df!")
  }

  # Model BBUM on negative fold change miRs ----
  # Get appropriate inputs
  df.filtered = df %>%
    dplyr::mutate(pvalue.filtered = dplyr::if_else(
      !excluded & !outlier,
      pvalue, NA_real_
    ))  # exclude ones that didn't make the cut
  BBUM_input_p = df.filtered$pvalue.filtered
  BBUM_input_s = df.filtered$BBUM.class
  BBUM.output = BBUM_corr(
    pvals = BBUM_input_p,
    signal_set = BBUM_input_s,
    add_starts = add_starts,
    only_start = only_start,
    limits = limits,
    auto_outliers = auto_outliers,
    rthres = rthres,
    rtrimmax = rtrimmax,
    atrimmax = atrimmax,
    quiet = quiet
  )

  # Identify additional outliers flagged ----
  trimmedN = BBUM.output$outlier_trim
  trimmed_genes = df.filtered %>%
    dplyr::filter(!BBUM.class) %>%  # bg class
    dplyr::arrange(pvalue.filtered) %>%
    dplyr::slice(0:trimmedN) %>%
    dplyr::pull(geneName)
  if(auto_outliers && trimmedN > 0 && !quiet){
    print(paste0(">> Additional trimmed outliers:"))
    print(as.character(trimmed_genes))
  }

  # Attach new cols with coefs ----
  coefs = BBUM.output$estim
  df.out = df %>%
    dplyr::ungroup() %>%
    dplyr::mutate(outlier = (outlier | geneName %in% trimmed_genes)) %>%  # include new outliers
    dplyr::mutate(
      BBUM.l  = coefs$lambda,
      BBUM.a  = coefs$a,
      BBUM.th = coefs$theta,
      BBUM.r  = coefs$r,            # put coefs into the data table
      pBBUM   = BBUM.output$pBBUMs  # BBUM correction with FDR correction
    )

  # Output ----
  list(data = df.out,
       BBUM.l  = coefs$lambda,
       BBUM.a  = coefs$a,
       BBUM.th = coefs$theta,
       BBUM.r  = coefs$r
  )
}
