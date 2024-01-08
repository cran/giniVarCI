fcompareCI <- function(y,
                       w,
                       Pi = NULL,
                       Pij = NULL,
                       PiU,
                       alpha = 0.05,
                       B = 1000L,
                       digitsgini = 2L,
                       digitsvar = 4L,
                       na.rm = TRUE,
                       plotCI = TRUE,
                       line.types = c(1L, 2L, 4L),
                       colors = c("red", "green", "blue"),
                       shapes = c(8L, 4L, 3L),
                       save.plot = FALSE,
                       large.sample = FALSE)
  {
  x <- gini <- upperlimit <- method <- varformula <- lowerlimit <- NULL
if (min(y) >= 0)
{
  if (missing(PiU))
  {
    SOL <- matrix(nrow = 13L, ncol = 4L)
    SOL[1L,]  <- unlist(fgini(y = y, w = w, method = 2L, interval = "zjackknife",      Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "HT",  large.sample = large.sample))
    SOL[2L,]  <- unlist(fgini(y = y, w = w, method = 2L, interval = "zjackknife",      Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "SYG", large.sample = large.sample))
    SOL[3L,]  <- unlist(fgini(y = y, w = w, method = 4L, interval = "zjackknife",      Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "HT",  large.sample = large.sample))
    SOL[4L,]  <- unlist(fgini(y = y, w = w, method = 4L, interval = "zjackknife",      Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "SYG", large.sample = large.sample))
    SOL[5L,]  <- unlist(fgini(y = y, w = w, method = 5L, interval = "zjackknife",      Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "HT",  large.sample = large.sample))
    SOL[6L,]  <- unlist(fgini(y = y, w = w, method = 5L, interval = "zjackknife",      Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "SYG", large.sample = large.sample))
    SOL[7L,]  <- unlist(fgini(y = y, w = w, method = 2L, interval = "zalinearization", Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "HT",  large.sample = large.sample))
    SOL[8L,]  <- unlist(fgini(y = y, w = w, method = 2L, interval = "zalinearization", Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "SYG", large.sample = large.sample))
    SOL[9L,]  <- unlist(fgini(y = y, w = w, method = 4L, interval = "zblinearization", Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "HT",  large.sample = large.sample))
    SOL[10L,] <- unlist(fgini(y = y, w = w, method = 4L, interval = "zblinearization", Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "SYG", large.sample = large.sample))
    SOL[11L,] <- unlist(fgini(y = y, w = w, method = 2L, interval = "pbootstrap",      Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "SYG", large.sample = large.sample))
    SOL[12L,] <- unlist(fgini(y = y, w = w, method = 4L, interval = "pbootstrap",      Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "SYG", large.sample = large.sample))
    SOL[13L,] <- unlist(fgini(y = y, w = w, method = 5L, interval = "pbootstrap",      Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "SYG", large.sample = large.sample))

    SOLround <- SOL
    SOLround[,1L:3L] <- round(SOL[,1L:3L], digits=digitsgini)
    SOLround[,4L]    <- round(SOL[,4L],    digits=digitsvar)

    Lmethods <- c("zjackknife", "zjackknife", "zjackknife", "zjackknife", "zjackknife", "zjackknife",
                  "zalinearization", "zalinearization", "zblinearization", "zblinearization",
                  "pbootstrap", "pbootstrap", "pbootstrap")
    Emethods  <- c(2L,2L,4L,4L,5L,5L,2L,2L,4L,4L,2L,4L,5L)
    Formulas <- c( rep(c("HT", "SYG"), 5L), NA, NA, NA)
    output <- data.frame(interval=Lmethods, method = Emethods, varformula = Formulas, gini = rep(NA,13L), lowerlimit = rep(NA,13L), upperlimit = rep(NA,13L), var.gini = rep(NA,13L))
    output[,4L:7L] <- SOLround

    if (plotCI | save.plot)
    {
      shapes <- shapes[1:2]
      if (length(line.types) != 3L) stop("argument 'line.types' must have length 3")
      if (length(colors) != 3L) stop("argument 'colors' must have length 3")
      if (length(shapes) != 2L) stop("PiU is not missing, argument 'shapes' must have at least length 2")
      output[ ,4L:7L] <- SOL
      output$x <- c(1L:6L, 1L:2L, 1L:2L, 1L:3L)*.5
      output$method <- as.character(output$method)
      output$interval <- factor(output$interval,
                                levels = c("zjackknife", "zalinearization",
                                           "zblinearization", "pbootstrap"))
      out2 <- output[-c(11L:13L), ]

      p <-
        ggplot2::ggplot(output, ggplot2::aes(x = x, y = gini)) +
        ggplot2::geom_point() +
        ggplot2::facet_grid(. ~ interval, scales = "free", space = "free") +
        ggplot2::theme_bw() +
        ggplot2::geom_point(data = out2,
                            mapping = ggplot2::aes(x = x, y = upperlimit,
                                                   color = method, shape = varformula),
                            size = 4) +
        ggplot2::geom_point(data = out2,
                            mapping = ggplot2::aes(x = x, y = lowerlimit,
                                                   color = method, shape = varformula),
                            size = 4) +
        ggplot2::geom_segment(ggplot2::aes(x = x, xend = x,
                                           y = lowerlimit, yend = upperlimit,
                                           color = method, linetype = method),
                              linewidth = 0.85) +
        ggplot2::ylab("Estimation of the Gini index") +
        ggplot2::xlab("intervals") +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank(),
                       panel.grid.minor.x = ggplot2::element_blank(),
                       panel.grid.major.x = ggplot2::element_blank(),
                       legend.position = "bottom") +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = 0.3)) +
        ggplot2::scale_shape_manual(values = c(shapes, 3L), breaks = c("HT", "SYG")) +
        ggplot2::scale_color_manual(values = colors) +
        ggplot2::scale_linetype_manual(values = line.types) +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 0.1)))

      if (plotCI) suppressWarnings(plot(p))

      output[,4L:7L] <- SOLround
      if (save.plot) output <- list("base.CI" = output[, 1L:7L], "plot" = p)
    }

  }
  else
  {
    SOL <- matrix(nrow = 18L, ncol = 4L)
    SOL[1L,]  <- unlist(fgini(y = y, w = w, method = 2L, interval = "zjackknife",      Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "HT",  large.sample = large.sample))
    SOL[2L,]  <- unlist(fgini(y = y, w = w, method = 2L, interval = "zjackknife",      Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "SYG", large.sample = large.sample))
    SOL[3L,]  <- unlist(fgini(y = y, w = w, method = 2L, interval = "zjackknife",      Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "HR",  large.sample = large.sample))
    SOL[4L,]  <- unlist(fgini(y = y, w = w, method = 4L, interval = "zjackknife",      Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "HT",  large.sample = large.sample))
    SOL[5L,]  <- unlist(fgini(y = y, w = w, method = 4L, interval = "zjackknife",      Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "SYG", large.sample = large.sample))
    SOL[6L,]  <- unlist(fgini(y = y, w = w, method = 4L, interval = "zjackknife",      Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "HR",  large.sample = large.sample))
    SOL[7L,]  <- unlist(fgini(y = y, w = w, method = 5L, interval = "zjackknife",      Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "HT",  large.sample = large.sample))
    SOL[8L,]  <- unlist(fgini(y = y, w = w, method = 5L, interval = "zjackknife",      Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "SYG", large.sample = large.sample))
    SOL[9L,]  <- unlist(fgini(y = y, w = w, method = 5L, interval = "zjackknife",      Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "HR",  large.sample = large.sample))
    SOL[10L,] <- unlist(fgini(y = y, w = w, method = 2L, interval = "zalinearization", Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "HT",  large.sample = large.sample))
    SOL[11L,] <- unlist(fgini(y = y, w = w, method = 2L, interval = "zalinearization", Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "SYG", large.sample = large.sample))
    SOL[12L,] <- unlist(fgini(y = y, w = w, method = 2L, interval = "zalinearization", Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "HR",  large.sample = large.sample))
    SOL[13L,] <- unlist(fgini(y = y, w = w, method = 4L, interval = "zblinearization", Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "HT",  large.sample = large.sample))
    SOL[14L,] <- unlist(fgini(y = y, w = w, method = 4L, interval = "zblinearization", Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "SYG", large.sample = large.sample))
    SOL[15L,] <- unlist(fgini(y = y, w = w, method = 4L, interval = "zblinearization", Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "HR",  large.sample = large.sample))
    SOL[16L,] <- unlist(fgini(y = y, w = w, method = 2L, interval = "pbootstrap",      Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "SYG", large.sample = large.sample))
    SOL[17L,] <- unlist(fgini(y = y, w = w, method = 4L, interval = "pbootstrap",      Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "SYG", large.sample = large.sample))
    SOL[18L,] <- unlist(fgini(y = y, w = w, method = 5L, interval = "pbootstrap",      Pi = Pi, Pij = Pij, PiU = PiU, alpha = alpha, B = B, na.rm = na.rm, varformula = "SYG", large.sample = large.sample))

    SOLround <- SOL
    SOLround[,1L:3L] <- round(SOL[,1L:3L], digits = digitsgini)
    SOLround[,4L]    <- round(SOL[,4L],    digits = digitsvar)

    Lmethods <- c(rep("zjackknife", 9L),
                  rep("zalinearization", 3L),
                  rep("zblinearization", 3L),
                  rep("pbootstrap", 3L) )
    Emethods  <- c(2L,2L,2L, 4L,4L,4L, 5L,5L,5L, 2L,2L,2L, 4L,4L,4L, 2L,4L,5L)
    Formulas <- c( rep(c("HT", "SYG", "HR"), 5L), NA, NA, NA)
    output <- data.frame(interval=Lmethods, method = Emethods, varformula = Formulas, gini = rep(NA,18L), lowerlimit = rep(NA,18L), upperlimit = rep(NA,18L), var.gini = rep(NA,18L))
    output[,4L:7L] <- SOLround

    if (plotCI | save.plot)
    {
      if (length(line.types) != 3L) stop("argument 'line.types' must have length 3")
      if (length(colors) != 3L) stop("argument 'colors' must have length 3")
      if (length(shapes) != 3L) stop("argument 'shapes' must have length 3")
      output[ ,4L:7L] <- SOL
      output$x <- c(1L:9L, 1L:3L, 1L:3L, 1L:3L)*.5
      output$method <- as.character(output$method)
      output$interval <- factor(output$interval,
                                levels = c("zjackknife", "zalinearization",
                                           "zblinearization", "pbootstrap"))
      output$varformula <- factor(output$varformula,
                                levels = c("HT", "SYG", "HR"))
      out2 <- output[-c(16L:18L), ]

      p <-
        ggplot2::ggplot(output, ggplot2::aes(x = x, y = gini)) +
        ggplot2::geom_point() +
        ggplot2::facet_grid(. ~ interval, scales = "free", space = "free") +
        ggplot2::theme_bw() +
        ggplot2::geom_point(data = out2,
                            mapping = ggplot2::aes(x = x, y = upperlimit,
                                                   color = method, shape = varformula),
                            size = 4) +
        ggplot2::geom_point(data = out2,
                            mapping = ggplot2::aes(x = x, y = lowerlimit,
                                                   color = method, shape = varformula),
                            size = 4) +
        ggplot2::geom_segment(ggplot2::aes(x = x, xend = x,
                                           y = lowerlimit, yend = upperlimit,
                                           color = method, linetype = method),
                              linewidth = 0.85) +
        ggplot2::ylab("Estimation of the Gini index") +
        ggplot2::xlab("intervals") +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank(),
                       panel.grid.minor.x = ggplot2::element_blank(),
                       panel.grid.major.x = ggplot2::element_blank(),
                       legend.position = "bottom") +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = 0.3)) +
        ggplot2::scale_shape_manual(values = c(shapes, 3L), breaks = c("HT", "SYG", "HR")) +
        ggplot2::scale_color_manual(values = colors) +
        ggplot2::scale_linetype_manual(values = line.types) +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 0.1)))

      if (plotCI) suppressWarnings(plot(p))

      output[,4L:7L] <- SOLround
      if (save.plot) output <- list("base.CI" = output[, 1L:7L], "plot" = p)
    }

  }
}
else
{
  warning("sample values must be non-negative")
  output <- NA
}
return(output)
}
