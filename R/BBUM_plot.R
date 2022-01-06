#' QC and results graphs for the BBUM-processed results data frame
#'
#' Useful graphs for checking and viewing the results of BBUM correction and
#'   significance calling.
#'
#' @param df.bbum The data.frame output of \code{BBUM_DEcorr()}, with any
#'   additional columns.
#' @param option Option of graph to plot. If a vector of length > 1 is provided,
#'   only the first element is used. Ignores case.
#' @param expressionCol A \code{str} of the name of the column used for the x
#'   axis of the MA graph. (Only for \code{option = "MA"}.)
#' @inheritParams BBUM_DEcorr
#'
#' @details The argument \code{expressionCol} allows plotting the MA graph
#'   against a specified column as the x axis (expression level). For instance,
#'   some may prefer to plot the mean normalized expression in control
#'   experiments only, rather than the default \code{"baseMean"} of DESeq2.
#' @details Graph \code{option}s are:
#' * \code{MA}: MA plot (\code{log2FoldChange} against \code{expressionCol})
#' * \code{volcano}: Volcano plot (\code{-log10(pvalue)} against
#'   \code{log2FoldChange})
#' * \code{hist}: p value histogram separated into signal and background set
#'   points, with BBUM model overlaid. Background histogram is normalized by a
#'   factor of \code{1 - theta} to account for the lack of primary effects for
#'   comparison.
#' * \code{ecdf}: p value ECDF separated into signal and background set
#'   points, with BBUM model overlaid.
#' * \code{ecdf_log}: ECDF in log scale for the p values, which helps to focus
#'   on the left-tail.
#' * \code{ecdf.corr}, \code{ecdf_log.corr}: Plots the \code{pBBUM} values
#'   instead to evaluate the FDR-corrected p values.
#' * \code{pp}: P-P plot to evaluate the goodness of fit.
#' * \code{pcorr}: Plot of p values from raw values to BBUM-FDR_corrected
#'   values, by data set. This plot is helpful for evaluating the correction of
#'   individual p values through the BBUM algorithm.
#' * \code{symm}: Modified symmetry plot of \code{-log10(p)} values, excluding
#'   hits, with \code{-log10(0)} as the mid-point instead of the median. Uses
#'   subsampling to account for the different number of points in the signal and
#'   background sets.
#' @details The most critical region of BBUM distribution for an appropriate
#'   correction for secondary effects is the "left-tail" around 0, where both
#'   primary and secondary beta components peak. An ECDF graph in log scale
#'   allows emphasis and better visualization of this region.
#' @details ECDF graphs are overlaid on the \code{x = y} diagonal line, which
#'   represents the uniform/null-only i.e. no secondary effects case.
#' @details Because the peak near \code{p = 0} is the most informative region
#'   for p values correction, a P-P plot is more appropriate to assess
#'   goodness-of-fit of BBUM models than a Q-Q plot.
#' @details Plot \code{symm} is for the validation of the assumption that the
#'   signal and background sets have roughly similar background (null and
#'   secondary effects) distributions of p values. As excluding hits does not
#'   exclude the false-negative region, there is still an expected discrepancy
#'   at low p values. The implemented color gradient attempts to reflect
#'   this expected up-deviation from the diagonal line when the fraction of
#'   remaining primary effects is large. Empirically, distributions that do not
#'   deviate from the \code{+/- log(10)} dashed lines when the expected primary
#'   effects fractions is low are symmetrical enough for accurate BBUM
#'   correction.
#'
#' @return \code{ggplot2} plot object.
#'
#' @examples
#' \dontrun{
#' BBUM_plot(df.bbum = res.BBUMcorr,
#'           option = "ecdf_log",
#'           expressionCol = "WTmean",
#'           pBBUM.alpha = 0.01)
#' }
#'
#' @export
BBUM_plot = function(
  df.bbum,
  option = c(
    "MA", "volcano",
    "hist",
    "ecdf", "ecdf_log",
    "ecdf.corr", "ecdf_log.corr",
    "pp",
    "pcorr",
    "symm"
  ),
  expressionCol = "baseMean",
  pBBUM.alpha = 0.05
) {

  # Set up theme ----
  customtheme = ggplot2::theme(
    axis.line = ggplot2::element_line(colour = "black"),
    axis.ticks = ggplot2::element_line(colour = "black"),
    axis.text = ggplot2::element_text(color = "black", size = 12),
    axis.title = ggplot2::element_text(color = "black", size = 12)
  )

  # Set up inputs ----
  plot.option = tolower(option[1])
  df.bbum = df.bbum %>%
    ungroup() %>%
    dplyr::filter(!excluded)
  # BBUM.l  = df.bbum$BBUM.l  %>% unique()
  # BBUM.a  = df.bbum$BBUM.a  %>% unique()
  BBUM.th = df.bbum$BBUM.th %>% unique()
  # BBUM.r  = df.bbum$BBUM.r  %>% unique()
  bbum.model.graph = tibble::tibble(
    p = sort(c(10^seq(-300,-3,0.05), seq(1E-3,1,1E-3)))
    ) %>%
    tidyr::crossing(df.bbum %>%
                      dplyr::select(BBUM.l, BBUM.a, BBUM.th, BBUM.r) %>%
                      dplyr::distinct()) %>%
    dplyr::mutate(dbum.model  = dbum(p, BBUM.l, BBUM.a),
                  pbum.model  = pbum(p, BBUM.l, BBUM.a),
                  dbbum.model = dbbum(p, BBUM.l, BBUM.a, BBUM.th, BBUM.r),
                  pbbum.model = pbbum(p, BBUM.l, BBUM.a, BBUM.th, BBUM.r)
    )

  FoldChange.lim = df.bbum$log2FoldChange %>%
                      unname() %>% abs() %>% max() %>% ceiling()
  topp = df.bbum %>%
    dplyr::filter(is.finite(pvalue), pvalue > 1E-300) %>%
    dplyr::group_by(BBUM.class) %>%
    dplyr::arrange(pvalue, .by_group = TRUE) %>%
    dplyr::slice(2) %>%
    dplyr::ungroup()
  down_lim = 10^(topp %>% dplyr::pull(pvalue) %>% log10() %>% mean())
  down_lim.trans = 10^(topp %>% dplyr::pull(pBBUM) %>% log10() %>% mean())
  pdir.lim = df.bbum %>%
    dplyr::filter(!BBUM.class) %>% dplyr::pull(pvalue) %>%
    stats::quantile(., probs = 0.01, na.rm = T, names = F) %>%
    -log10(.) %>% ceiling(.)+1
  pdir.breaks = seq(-350, 350,
                    dplyr::if_else(pdir.lim >= 10,
                                   ceiling((pdir.lim/3)/10)*10,
                                   2))
  # Plots ----
  if(plot.option == "ma"){

    ## Fold-change i.e. MA graph ----
    return(df.bbum %>%
      dplyr::arrange(BBUM.fct, abs(log2FoldChange)) %>%
      ggplot2::ggplot(ggplot2::aes(x = !!dplyr::sym(expressionCol),
                                   y = log2FoldChange,
                                   color = BBUM.fct)) +
      ggplot2::geom_hline(yintercept = 0, color = "black", alpha = 0.5,
                          size = 0.5, linetype = "dashed") +
      ggplot2::geom_point(alpha = 0.75, size = 1, shape = 16) +
      ggplot2::scale_color_manual(
        breaks = c("none","hit","outlier"),
        values = c(
          "gray80",
          "red3",
          "gray25"
        )) +
      ggplot2::scale_x_continuous(trans = "log10", breaks = 10^seq(0,100,2)) +
      ggplot2::scale_y_continuous(breaks = seq(-100,100,2)) +
      ggplot2::coord_cartesian(ylim = c(-FoldChange.lim,FoldChange.lim),
                               xlim = c(1, NA)) +
      ggplot2::labs(y = "Fold change (log2)", x = "Expression",
                    title = "MA plot", color = "Gene category") +
      ggplot2::theme_classic(base_size = 12) + customtheme
    )

  } else if(plot.option == "volcano"){

    ## Volcano plot ----
    return(df.bbum %>%
      dplyr::arrange(BBUM.fct, abs(log2FoldChange)) %>%
      ggplot2::ggplot(ggplot2::aes(x = log2FoldChange,
                                   y = -log10(pvalue),
                                   color = BBUM.fct)) +
      ggplot2::geom_vline(xintercept = 0, color = "black", alpha = 0.5,
                         size = 0.5) +
      ggplot2::geom_hline(yintercept = 0, color = "black", alpha = 0.5,
                          size = 0.5) +
      ggplot2::geom_point(alpha = 0.75, size = 1, shape = 16) +
      ggplot2::scale_color_manual(
       breaks = c("none","hit","outlier"),
       values = c(
         "gray80",
         "red3",
         "gray25"
       )) +
      ggplot2::scale_x_continuous(breaks = seq(-100,100,2)) +
      ggplot2::coord_cartesian(xlim = c(-FoldChange.lim,FoldChange.lim)) +
      ggplot2::labs(x = "Fold change (log2)", y = "-log10(p value)",
                   title = "Volcano plot", color = "Gene category") +
      ggplot2::theme_classic(base_size = 12) + customtheme
    )

  } else if(plot.option == "hist") {

    ## Histogram ----
    binN = (length(df.bbum$pBBUM[!is.na(df.bbum$pBBUM)])/2)^(1/3)*10
      # 10 times more precise than traditional rule of "N^(1/3)"
    binwidth = 1/binN
    df.bbum_plot_bin = df.bbum %>%
      dplyr::mutate(pvalue.binned.raw = ggplot2::cut_interval(pvalue, n = binN),
                    pvalue.binned = as.numeric(
                      gsub("(^[\\(\\[])|(,.*)", "", pvalue.binned.raw)
                        # extract number from range text
                    )
      ) %>%
      dplyr::group_by(BBUM.class, pvalue.binned.raw, pvalue.binned) %>%
      dplyr::tally(.) %>%
      dplyr::group_by(BBUM.class) %>%
      dplyr::mutate(freq = n/sum(n)/binwidth) %>%
      dplyr::mutate(freq = dplyr::if_else(BBUM.class, freq, freq * (1-BBUM.th)))
    return(df.bbum_plot_bin %>%
      ggplot2::ggplot(ggplot2::aes(
        x = pvalue.binned,
        y = freq,
        color = factor(BBUM.class, levels = c(TRUE,FALSE)))) +
      ggplot2::geom_hline(yintercept = 0, color = "gray60", alpha = 0.75,
                          size = 0.5) +
      ggplot2::geom_vline(xintercept = c(0,1), color = "gray60", alpha = 0.75,
                          size = 0.5) +
      ggplot2::geom_step(stat = "identity", alpha = 0.75, size = 1) +
      ggplot2::geom_line(data = bbum.model.graph, inherit.aes = F,
                         ggplot2::aes(x = p, y = dbum.model*(1-BBUM.th)),
                         color = "goldenrod3",
                         alpha = 0.5, size = 0.5) +
      ggplot2::geom_line(data = bbum.model.graph, inherit.aes = FALSE,
                         ggplot2::aes(x = p, y = dbbum.model),
                         color = "turquoise4",
                         alpha = 0.5, size = 0.5) +
      ggplot2::scale_color_manual(breaks = c(FALSE, TRUE),
                                  values = c("goldenrod3","turquoise4"),
                                  labels = c("Background", "Signal")
      ) +
      ggplot2::labs(x = "p value", y = "Density of probability",
                    title = "Histogram of p values", color = "Data set") +
      ggplot2::coord_cartesian(ylim = c(0, max(df.bbum_plot_bin$freq)*1.2)) +
      ggplot2::theme_classic(base_size = 12) + customtheme
    )

  } else if(plot.option == "ecdf") {

    ## ECDF ----
    return(df.bbum %>%
      ggplot2::ggplot(ggplot2::aes(
        x = pvalue,
        color = factor(BBUM.class, levels = c(TRUE,FALSE)))) +
      ggplot2::geom_hline(yintercept = 0, color = "gray60",
                          alpha = 0.75, size = 0.5) +
      ggplot2::geom_vline(xintercept = c(0,1), color = "gray60",
                          alpha = 0.75, size = 0.5) +
      ggplot2::geom_abline(slope = 1, intercept = 0, color = "gray60",
                           alpha = 0.75, size = 0.7, linetype = "dashed") +
      ggplot2::geom_step(stat = "ecdf", alpha = 0.75, size = 0.5) +
      ggplot2::geom_line(data = bbum.model.graph, inherit.aes = F,
                         ggplot2::aes(x = p, y = pbum.model),
                         color = "goldenrod3",
                         alpha = 0.5, size = 0.7) +
      ggplot2::geom_line(data = bbum.model.graph, inherit.aes = FALSE,
                         ggplot2::aes(x = p, y = pbbum.model),
                         color = "turquoise4",
                         alpha = 0.5, size = 0.7) +
      ggplot2::scale_color_manual(breaks = c(FALSE, TRUE),
                                  values = c("goldenrod3","turquoise4"),
                                  labels = c("Background", "Signal")
      ) +
      ggplot2::labs(x = "p value", y = "Cumulative frequency",
                    title = "ECDF of p values", color = "Data set") +
      ggplot2::theme_classic(base_size = 12) + customtheme
    )

  } else if(plot.option == "ecdf_log") {

    ## ECDF, in log ----
    return(df.bbum %>%
      ggplot2::ggplot(ggplot2::aes(
        x = pvalue,
        color = factor(BBUM.class, levels = c(TRUE,FALSE)))) +
      ggplot2::geom_hline(yintercept = 0, color = "gray60",
                          alpha = 0.75, size = 0.5) +
      ggplot2::geom_vline(xintercept = c(1), color = "gray60",
                          alpha = 0.75, size = 0.5) +
      ggplot2::geom_line(data = tibble::tibble(
        pvalue     = 10^(seq(-10,0,0.01)),
        pvalue.lin = 10^(seq(-10,0,0.01))
      ), ggplot2::aes(y = pvalue.lin), color = "gray60", alpha = 0.75,
      size = 0.7, linetype = "dashed") +
      ggplot2::geom_step(stat = "ecdf", alpha = 0.75, size = 0.5) +
      ggplot2::geom_line(data = bbum.model.graph, inherit.aes = F,
                         ggplot2::aes(x = p, y = pbum.model),
                         color = "goldenrod3",
                         alpha = 0.5, size = 0.7) +
      ggplot2::geom_line(data = bbum.model.graph, inherit.aes = FALSE,
                         ggplot2::aes(x = p, y = pbbum.model),
                         color = "turquoise4",
                         alpha = 0.5, size = 0.7) +
      ggplot2::scale_color_manual(breaks = c(FALSE, TRUE),
                                  values = c("goldenrod3","turquoise4"),
                                  labels = c("Background", "Signal")
      ) +
      ggplot2::scale_x_continuous(trans = "log10") +
      ggplot2::coord_cartesian(xlim = c(down_lim,1)) +
      ggplot2::labs(x = "p value", y = "Cumulative frequency",
                    title = "ECDF of p values, in log", color = "Data set") +
      ggplot2::theme_classic(base_size = 12) + customtheme
    )

  } else if(plot.option == "ecdf.corr") {

    ## ECDF, transformed ----
    return(df.bbum %>%
      ggplot2::ggplot(ggplot2::aes(
        x = pBBUM,
        color = factor(BBUM.class, levels = c(TRUE,FALSE)))) +
      ggplot2::geom_hline(yintercept = 0, color = "gray60",
                          alpha = 0.75, size = 0.5) +
      ggplot2::geom_vline(xintercept = c(0,1), color = "gray60",
                          alpha = 0.75, size = 0.5) +
      ggplot2::geom_vline(xintercept = pBBUM.alpha, color = "salmon4",
                          alpha = 0.75, size = 0.5, linetype = "dashed") +
      ggplot2::geom_abline(slope = 1, intercept = 0, color = "gray60",
                           alpha = 0.75, size = 0.7, linetype = "dashed") +
      ggplot2::geom_step(stat = "ecdf", alpha = 0.75, size = 0.5) +
      ggplot2::scale_color_manual(breaks = c(FALSE, TRUE),
                                  values = c("goldenrod3","turquoise4"),
                                  labels = c("Background", "Signal")
      ) +
      ggplot2::labs(x = "pBBUM value", y = "Cumulative frequency",
                    title = "ECDF of BBUM-FDR-adjusted p values", color = "Data set") +
      ggplot2::theme_classic(base_size = 12) + customtheme
    )

  } else if(plot.option == "ecdf_log.corr") {

    ## ECDF, transformed, in log ----
    return(df.bbum %>%
      ggplot2::ggplot(ggplot2::aes(
        x = pBBUM,
        color = factor(BBUM.class, levels = c(TRUE,FALSE)))) +
      ggplot2::geom_hline(yintercept = 0, color = "gray60",
                          alpha = 0.75, size = 0.5) +
      ggplot2::geom_vline(xintercept = c(1), color = "gray60",
                          alpha = 0.75, size = 0.5) +
      ggplot2::geom_vline(xintercept = pBBUM.alpha, color = "salmon4",
                          alpha = 0.75, size = 0.5, linetype = "dashed") +
      ggplot2::geom_step(stat = "ecdf", alpha = 0.75, size = 0.5) +
      ggplot2::scale_color_manual(breaks = c(FALSE, TRUE),
                                  values = c("goldenrod3","turquoise4"),
                                  labels = c("Background", "Signal")
      ) +
      ggplot2::scale_x_continuous(trans = "log10") +
      ggplot2::coord_cartesian(xlim = c(down_lim.trans,1)) +
      ggplot2::labs(x = "pBBUM value", y = "Cumulative frequency",
                    title = "ECDF of BBUM-FDR-adjusted p values, in log", color = "Data set") +
      ggplot2::theme_classic(base_size = 12) + customtheme
    )

  } else if(plot.option == "pp") {

    ## P-P plot for goodness of fit ----

    return(df.bbum %>%
             dplyr::mutate(p = dplyr::if_else(
               BBUM.class,
               pbbum(q = pvalue,
                     lambda = BBUM.l,  a = BBUM.a,
                     theta  = BBUM.th, r = BBUM.r
               ),
               pbbum(q = pvalue,
                     lambda = BBUM.l, a = BBUM.a,
                     theta  = 0,      r = 1
               )
             )) %>%
             ggplot2::ggplot(ggplot2::aes(
               x = p,
               color = factor(BBUM.class, levels = c(TRUE,FALSE)))) +
             ggplot2::geom_hline(yintercept = 0, color = "gray60",
                                 alpha = 0.75, size = 0.5) +
             ggplot2::geom_vline(xintercept = c(0,1), color = "gray60",
                                 alpha = 0.75, size = 0.5) +
             ggplot2::geom_abline(slope = 1, intercept = 0, color = "gray60",
                                  alpha = 0.75, size = 0.7, linetype = "dashed") +
             ggplot2::geom_step(stat = "ecdf", alpha = 0.75, size = 0.5) +
             ggplot2::scale_color_manual(breaks = c(FALSE, TRUE),
                                         values = c("goldenrod3","turquoise4"),
                                         labels = c("Background", "Signal")
             ) +
             ggplot2::labs(x = "Theoretical cumulative distribution",
                           y = "Theoretical cumulative distribution",
                           title = "P-P plot for goodness-of-fit", color = "Data set") +
             ggplot2::theme_classic(base_size = 12) + customtheme
    )

  } else if(plot.option == "pcorr") {

    ## p values plot in two directions of two data sets ----
    return(df.bbum %>%
      dplyr::select(geneName, BBUM.fct, pvalue, pBBUM, BBUM.class) %>%
      tidyr::pivot_longer(cols = tidyr::starts_with("p"),
                          names_to = "p", values_to = "val") %>%
      dplyr::mutate(p = factor(p, levels = c("pvalue","pBBUM"))) %>%
      dplyr::mutate(p.dir = -log10(val) * dplyr::if_else(BBUM.class,+1,-1)) %>%
      dplyr::arrange(BBUM.fct) %>%
      ggplot2::ggplot(ggplot2::aes(y = p, x = p.dir,
                                   color = BBUM.fct,
                                   group = geneName
      )) +
      ggplot2::geom_path(alpha = 0.3, size = 0.25) +
      ggplot2::geom_point(alpha = 0.6, shape = 16) +
      ggplot2::scale_color_manual(
        breaks = c("none","hit","outlier"),
        values = c(
          "gray80",
          "red3",
          "gray25"
        )) +
      ggplot2::geom_vline(xintercept = 0, color = "black", size = 0.5,
                          linetype = "dashed") +
      ggplot2::geom_vline(xintercept = c(log10(pBBUM.alpha),
                                         -log10(pBBUM.alpha)),
                          color = "salmon1", alpha = 0.5, size = 0.5,
                          linetype = "dashed") +
      ggplot2::scale_x_continuous(breaks = pdir.breaks,
                                  labels = 10^-abs(pdir.breaks)) +
      ggplot2::scale_y_discrete(limits = rev) +
      ggplot2::coord_cartesian(xlim = c(-pdir.lim, pdir.lim*2)) +
      ggplot2::labs(x = "Value", y = "Statistic",
                    title = "Correction of p values",
                    color = "Gene category") +
      ggplot2::theme_classic(base_size = 12) + customtheme
    )

  } else if(plot.option == "symm") {

    ## Symmetrical assumption check graph, with -log-transformed p values ----
    pts_p = df.bbum %>%
      dplyr::filter(!is.na(pBBUM), !outlier) %>%
      dplyr::filter(!BBUM.hits) %>%
      dplyr::mutate(logp = -log10(pvalue))
    v_p = pts_p %>% dplyr::pull(logp) %>% sort()
    v_p.up   = pts_p %>% dplyr::filter( BBUM.class) %>%
      dplyr::pull(logp) %>% sort()
    v_p.down = pts_p %>% dplyr::filter(!BBUM.class) %>%
      dplyr::pull(logp) %>% sort()
    symmN = min( length(v_p.up), length(v_p.down) )
    symm.samples = tibble::tibble(sampleid = paste0("n",seq(1,20,1)))
    if(symmN == length(v_p.up)) {
      symm.samples = symm.samples %>%
        dplyr::rowwise() %>%
        dplyr::do(., tibble::tibble(
          sampleid = .$sampleid,
          v_p.up.sample = v_p.up,
          v_p.down.sample = sample(v_p.down, symmN, replace = FALSE) %>% sort()
        )
        )
    } else {
      symm.samples = symm.samples %>%
        dplyr::rowwise() %>%
        dplyr::do(., tibble::tibble(
          sampleid = .$sampleid,
          v_p.up.sample = sample(v_p.up, symmN, replace = FALSE) %>% sort(),
          v_p.down.sample = v_p.down
        )
        )
    }
    dt.coefs = df.bbum %>%
      dplyr::select(tidyr::starts_with("BBUM")) %>%
      dplyr::distinct()
    symm.samples = symm.samples %>%
      tidyr::crossing(., dt.coefs) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        p.up = 10^(-v_p.up.sample),
        primfrac.here = BBUM_primfrac(
          x = p.up,
          lambda = BBUM.l,
          a = BBUM.a,
          theta = BBUM.th,
          r = BBUM.r
        )
      ) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(primfrac.here)
    return(symm.samples %>%
      ggplot2::ggplot(ggplot2::aes(x = v_p.down.sample, y = v_p.up.sample,
                                   group = sampleid, color = primfrac.here)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, color = "black",
                           alpha = 0.5, size = 0.75, linetype = "dotted") +
      ggplot2::geom_abline(slope = 1, intercept = c(-1,1), color = "black",
                           alpha = 0.1, size = 0.75, linetype = "dotted") +
      ggplot2::geom_point(alpha = 0.1, size = 1, shape = 16) +
      ggplot2::scale_color_gradient(low = "gray75", high = "red4",
                                    limits = c(0,1)) +
      ggplot2::coord_fixed(
        xlim = c(0,3), ylim = c(0,3)
      ) +
      ggplot2::scale_x_continuous(breaks = 0:3,labels=10^-(0:3)) +
      ggplot2::scale_y_continuous(breaks = 0:3,labels=10^-(0:3)) +
      ggplot2::labs(x = "p value, background set", y = "p value, signal set",
                    title = "Symmetry plot of non-hits p values",
                    color = "Estimated fraction of primary effects") +
      ggplot2::theme_classic(base_size = 12) + customtheme
    )

  } else {
    stop("Plot option not available.")
  }
}
