# adapted from roahd::outliergram
prep_outliergram <- function (data, grid, p_check = 0.05,
          Fvalue = 1.5, adjust = FALSE, display = TRUE, xlab = NULL,
          ylab = NULL, main = NULL, ...) {
  MBD_data = NULL
  MEI_data = NULL
  fData <- roahd::fData(grid, data)
  N = fData$N
  grid = seq(fData$t0, fData$tP, length.out = fData$P)
  out = roahd:::.outliergram(fData, MBD_data = MBD_data, MEI_data = MEI_data,
                       p_check = p_check, Fvalue = Fvalue, shift = TRUE)

  a_0_2 = -2/(N * (N - 1))
  a_1 = 2 * (N + 1)/(N - 1)
  grid_1D <- seq(0, 1, l = 100)
  limits <-
    data.frame(MEI = grid_1D,
          upper = a_0_2 + a_1 * grid_1D + a_0_2 * N^2 * grid_1D^2,
          lower = a_0_2 + a_1 * grid_1D + a_0_2 * N^2 * grid_1D^2 - out$Q_d3 - Fvalue * out$IQR_d)
  return(c(out, limits = list(limits)))
}

plot_compare_real <- function(data, dataname, legend = FALSE) {
  grid <- seq(1:ncol(data))
  og <- prep_outliergram(data, grid)
  ms <- fdaoutlier::msplot(data, plot = FALSE)
  color <- factor(10*(1:nrow(data) %in% og$ID_SO) +  (1:nrow(data) %in% ms$outliers),
                  levels = c(0, 1, 10, 11),
                  labels = c("inlier", "MS outlier", "MBD-MEI outlier", "both outlier"))
  clr_scl <- scale_color_manual(
    values = c(scales::alpha("black", 0.1), scales::alpha(c("red", "orange", "purple"), 0.5)),
    limits = levels(color))
  p1 <- plot_funs(data, col = color) +  ggtitle(dataname) +
    clr_scl + xlab("t") + ylab("x(t)") + guides(color = guide_none()) +
    scale_x_continuous(breaks= NULL, minor_breaks = NULL)
  p2 <- ggplot(data.frame(MBD = og$MBD_data, MEI = og$MEI_data, color = color)) +
    geom_ribbon(data = tail(og, 1)[[1]], aes(x=MEI, ymin = lower, ymax = upper), fill = "darkgrey") +
    geom_point(aes(y = MBD, x = MEI, color = color)) +
    ggtitle("Outliergram") + clr_scl + guides(color = guide_legend(title = "")) +
    theme(legend.position = "bottom")
  p3 <- ggplot(data.frame(append(ms[3:4], list(color = color)))) +
    geom_point(aes(x = mean_outlyingness, y = var_outlyingness, color = color)) +
    xlab("MO") + ylab("VO") + ggtitle("Magnitude-Shape") + clr_scl +
    guides(color = guide_legend(title = "")) + theme(legend.position = "bottom")
  if (!legend) return(p1 + p3 + p2 + plot_layout(guides = "collect") & theme_bw(base_size = 10) & theme(legend.position = "none"))
  p1 + p3 + p2 + plot_layout(guides = "collect") & theme_bw(base_size = 10) & theme(legend.position = "bottom")
}

plt_outliergram_real <- list(
  ecg = plot_compare_real(out_ecg200,  "ECG"),
  octane = plot_compare_real(t(t(oct) - colMeans(oct)),  "Octane"),
  weather = plot_compare_real(weather, "Spanish weather"),
  tecator = plot_compare_real(t(t(tec) - colMeans(tec)), "Tecator"),
  wine = plot_compare_real(t(t(wine) - colMeans(wine)), "Wine", legend = TRUE))

