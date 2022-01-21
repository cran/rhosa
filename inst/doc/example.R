## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7
)

## ----setup--------------------------------------------------------------------
library(rhosa)

## ----init---------------------------------------------------------------------
set.seed(1)
f_1 <- 0.35
f_2 <- 0.2
D <- function(t) {
    omega_a <- runif(1, min = 0, max = 2 * pi)
    omega_b <- runif(1, min = 0, max = 2 * pi)
    omega_c <- runif(1, min = 0, max = 2 * pi)
    omega_d <- runif(1, min = 0, max = 2 * pi)
    wave_a <- function(t) cos(2 * pi * f_1 * t + omega_a)
    wave_b <- function(t) cos(2 * pi * f_2 * t + omega_b)
    wave_c <- function(t) cos(2 * pi * f_1 * t + omega_c)
    wave_d <- function(t) cos(2 * pi * f_2 * t + omega_d)
    curve_v <- function(t) wave_a(t) + wave_b(t) + wave_a(t) * wave_b(t)
    curve_w <- function(t) wave_c(t) + wave_d(t) + wave_c(t) * wave_b(t)
    data.frame(v = curve_v(t) + rnorm(length(t)),
               w = curve_w(t) + rnorm(length(t)))
}

## ----ts, fig.cap = "v and w."-------------------------------------------------
data <- D(seq_len(2048))
with(data, {
    plot(seq_len(100), head(v, 100), type = "l", col = "green", ylim = c(-3, 3), xlab = "t", ylab = "value")
    lines(seq_len(100), head(w, 100), col = "orange")
})

## ----psd, fig.cap = "Spectral density estimation via periodograms.", fig.show = "hold"----
with(data, {
    spectrum(v, main = "v", col = "green")
    spectrum(w, main = "w", col = "orange")
})

## ----bc-----------------------------------------------------------------------
x <- replicate(100, D(seq_len(128)), simplify = FALSE)
m_v <- do.call(cbind, Map(function(d) {d$v}, x))
m_w <- do.call(cbind, Map(function(d) {d$w}, x))

library(rhosa)

bc_v <- bicoherence(m_v, window_function = 'hamming')
bc_w <- bicoherence(m_w, window_function = 'hamming')

## ----plot_bicoherence---------------------------------------------------------
library(ggplot2)

plot_bicoherence <- function(bc) {
    ggplot(bc, aes(f1, f2)) +
        geom_raster(aes(fill = value)) +
        scale_fill_gradient(limits = c(0, 10)) +
        coord_fixed()
}

## ----viz_bc_v, fig.cap = "v's estimated magnitude-squared bicoherence."-------
plot_bicoherence(bc_v)

## ----viz_bc_w, fig.cap = "w's estimated magnitude-squared bicoherence."-------
plot_bicoherence(bc_w)

