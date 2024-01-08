
icompareCI <- function(y,
                       B = 1000L,
                       alpha = 0.05,
                       plotCI = TRUE,
                       digitsgini = 2L,
                       digitsvar = 4L,
                       cum.sums = NULL,
                       na.rm = TRUE,
                       precisionEL = 1e-4,
                       maxiterEL = 100L,
                       line.types = c(1L, 2L),
                       colors = c("red", "green"),
                       save.plot = FALSE
                       )
{
  gini <- upperlimit <- bc <- lowerlimit <- NULL
if (min(y) >= 0)
{
    SOL <- matrix(nrow = 20L, ncol = 4L)
    SOL[1L,]  <- unlist(igini(y = y, bias.correction = FALSE, interval = "zjackknife",     B = B, alpha = alpha, cum.sums = cum.sums, na.rm = na.rm))
    SOL[2L,]  <- unlist(igini(y = y, bias.correction = TRUE,  interval = "zjackknife",     B = B, alpha = alpha, cum.sums = cum.sums, na.rm = na.rm))
    SOL[3L,]  <- unlist(igini(y = y, bias.correction = FALSE, interval = "tjackknife",     B = B, alpha = alpha, cum.sums = cum.sums, na.rm = na.rm))
    SOL[4L,]  <- unlist(igini(y = y, bias.correction = TRUE,  interval = "tjackknife",     B = B, alpha = alpha, cum.sums = cum.sums, na.rm = na.rm))
    SOL[5L,]  <- unlist(igini(y = y, bias.correction = FALSE, interval = "zalinearization", B = B, alpha = alpha, cum.sums = cum.sums, na.rm = na.rm))
    SOL[6L,]  <- unlist(igini(y = y, bias.correction = TRUE,  interval = "zalinearization", B = B, alpha = alpha, cum.sums = cum.sums, na.rm = na.rm))
    SOL[7L,]  <- unlist(igini(y = y, bias.correction = FALSE, interval = "talinearization", B = B, alpha = alpha, cum.sums = cum.sums, na.rm = na.rm))
    SOL[8L,]  <- unlist(igini(y = y, bias.correction = TRUE,  interval = "talinearization", B = B, alpha = alpha, cum.sums = cum.sums, na.rm = na.rm))
    SOL[9L,]  <- unlist(igini(y = y, bias.correction = FALSE, interval = "zblinearization", B = B, alpha = alpha, cum.sums = cum.sums, na.rm = na.rm))
    SOL[10L,] <- unlist(igini(y = y, bias.correction = TRUE,  interval = "zblinearization", B = B, alpha = alpha, cum.sums = cum.sums, na.rm = na.rm))
    SOL[11L,] <- unlist(igini(y = y, bias.correction = FALSE, interval = "tblinearization", B = B, alpha = alpha, cum.sums = cum.sums, na.rm = na.rm))
    SOL[12L,] <- unlist(igini(y = y, bias.correction = TRUE,  interval = "tblinearization", B = B, alpha = alpha, cum.sums = cum.sums, na.rm = na.rm))
    SOL[13L,] <- unlist(igini(y = y, bias.correction = FALSE, interval = "pbootstrap",   B = B, alpha = alpha, cum.sums = cum.sums, na.rm = na.rm))
    SOL[14L,] <- unlist(igini(y = y, bias.correction = TRUE,  interval = "pbootstrap",  B = B, alpha = alpha, cum.sums = cum.sums, na.rm = na.rm))
    SOL[15L,] <- unlist(igini(y = y, bias.correction = FALSE, interval = "BCa",   B = B, alpha = alpha, na.rm = na.rm, cum.sums = cum.sums))
    SOL[16L,] <- unlist(igini(y = y, bias.correction = TRUE,  interval = "BCa",   B = B, alpha = alpha, na.rm = na.rm, cum.sums = cum.sums))
    SOL[17L,] <- unlist(igini(y = y, bias.correction = FALSE, interval = "ELchisq", B = B, alpha = alpha, na.rm = na.rm, cum.sums = cum.sums, precisionEL = precisionEL, maxiterEL = maxiterEL))
    SOL[18L,] <- unlist(igini(y = y, bias.correction = TRUE,  interval = "ELchisq", B = B, alpha = alpha, na.rm = na.rm, cum.sums = cum.sums, precisionEL = precisionEL, maxiterEL = maxiterEL))
    SOL[19L,] <- unlist(igini(y = y, bias.correction = FALSE, interval = "ELboot", B = B, alpha = alpha, na.rm = na.rm, cum.sums = cum.sums, precisionEL = precisionEL, maxiterEL = maxiterEL))
    SOL[20L,] <- unlist(igini(y = y, bias.correction = TRUE,  interval = "ELboot", B = B, alpha = alpha, na.rm = na.rm, cum.sums = cum.sums, precisionEL = precisionEL, maxiterEL = maxiterEL))

    SOLround <- SOL
    SOLround[,1L:3L] <- round(SOL[,1L:3L], digits=digitsgini)
    SOLround[,4L]    <- round(SOL[,4L],    digits=digitsvar)

    Lmethods <- c("zjackknife",      "zjackknife",       "tjackknife",      "tjackknife",
                  "zalinearization", "zalinearization",  "talinearization", "talinearization",
                  "zblinearization", "zblinearization",  "tblinearization", "tblinearization",
                  "pbootstrap", "pbootstrap", "BCa", "BCa", "ELchisq", "ELchisq", "ELboot", "ELboot")
    Lbc <- rep(c(FALSE, TRUE), 10L)
    output <- data.frame(interval = Lmethods, bc = Lbc, gini = rep(NA, 20L), lowerlimit = rep(NA,20L), upperlimit = rep(NA,20L), var.gini = rep(NA,20L))
    output[,3L:6L] <- SOLround

    if (plotCI | save.plot)
    {
      if (length(line.types) != 2L) stop("argument 'line.types' must have length = 2")
      if (length(colors) != 2L) stop("object 'colors' must have length = 2")

      output[,3L:6L] <- SOL
      output$interval <- factor(output$interval,
                                levels = c("zjackknife", "tjackknife",
                                           "zalinearization", "talinearization",
                                           "zblinearization", "tblinearization",
                                           "pbootstrap", "BCa", "ELchisq", "ELboot"))

      p <-
        ggplot2::ggplot(output, ggplot2::aes(x = rep(c(1,2), 10)*0.3, y = gini)) +
        ggplot2::geom_point() +
        ggplot2::facet_grid(. ~ interval, scales = "free", space = "free") +
        ggplot2::theme_bw() +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = lowerlimit, ymax = upperlimit,
                                            color = bc, linetype = bc),
                               width = .2, linewidth = 0.85) +
        ggplot2::ylab("Estimation of the Gini index") +
        ggplot2::xlab("intervals") +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank(),
                       panel.grid.minor.x = ggplot2::element_blank(),
                       panel.grid.major.x = ggplot2::element_blank(),
                       legend.position = "bottom") +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = 0.3)) +
        ggplot2::scale_color_manual(values = colors) +
        ggplot2::scale_linetype_manual(values = line.types)

      if (plotCI) plot(p)

      output[,3L:6L] <- SOLround
      if (save.plot) output <- list("base.CI" = output, "plot" = p)

    }
}
else
{
  warning("sample values must be non-negative")
  output <- NA
}
return(output)
}
