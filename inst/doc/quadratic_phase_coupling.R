## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7
)

## ----setup--------------------------------------------------------------------
library(rhosa)

## ----sec1_definition----------------------------------------------------------
triple_lambda <- function(a, b) c(a, b, a + b)
lambda <- c(triple_lambda(0.55, 0.75), triple_lambda(0.6, 0.8))
x1 <- function(k) {
    set.seed(k)
    init_phi <- runif(5, min = 0, max = 2*pi)
    phi <- c(init_phi, init_phi[4] + init_phi[5])
    function(t) do.call(sum, Map(function(l, p) cos(l * t + p), lambda, phi))
}
observe <- function(f) {
    sapply(seq_len(256), f)
}
N1 <- 100
m1 <- do.call(cbind, Map(observe, Map(x1, seq_len(N1))))

## ----sec1_samples, fig.cap = "100 realizations of `x1`."----------------------
ith_sample <- function(i) {
    data.frame(i = i, t = seq_len(256), v = m1[,i])
}
r1 <- do.call(rbind, Map(ith_sample, seq_len(100)))

library(ggplot2)

ggplot(r1) +
    geom_line(aes(t, v)) +
    facet_wrap(vars(i)) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank())

## ----sec1_X1, fig.cap = "The 256 time points in the first realization of `x1`."----
plot(m1[,1], type = "l", ylim = c(-5, 5), xlab = "t", ylab = "")

## ----sec1_spec, fig.cap = "The spectrum estimation of `x1`.", fig.height = 4----
spectrum(m1)

## ----sec1_bc1-----------------------------------------------------------------
bc1 <- bicoherence(m1)

## ----sec1_heatmap_bc1, fig.cap = "`x1`'s estimated bicoherence.", fig.height = 4----
heatmap_bicoherence <- function(bc) {
    ggplot(bc) +
        geom_raster(aes(f1, f2, fill = value)) +
        coord_fixed() +
        scale_alpha(guide = "none")
}

heatmap_bicoherence(bc1)

## ----tc-----------------------------------------------------------------------
Fcoef1 <- 1.2
Fcoef2 <- 0.7
Fcoef3 <- 0.8

i1 <- function(x, p) {2 * asin(sin(Fcoef1 * x + p))}
i2 <- function(x, p) {ifelse(cos(Fcoef2 * x + p) >= 0, -1, 1)}
i3 <- function(x, p) {cos(Fcoef3 * x + p)}

Qcoef <- 0.3

tc <- function(k) {
    set.seed(k)
    ps <- runif(3, min = 0, max = 2*pi)
    function(x) {
        c1 <- i1(x, ps[1]) + rnorm(length(x), mean = 0, sd = 1)
        c2 <- i2(x, ps[2]) + rnorm(length(x), mean = 0, sd = 1)
        c3 <- Qcoef * c1 * c2 +
            i3(x, ps[3]) + rnorm(length(x), mean = 0, sd = 1)
        data.frame(c1, c2, c3)
    }
}

N2 <- 100

sample_tc <- function() {
    Map(function(f) {f(seq_len(256))}, Map(tc, seq_len(N2)))
}

c1_data_frame <- function(y) {
    do.call(cbind, Map(function(k) {y[[k]]$c1}, seq_len(N2)))
}

c2_data_frame <- function(y) {
    do.call(cbind, Map(function(k) {y[[k]]$c2}, seq_len(N2)))
}

c3_data_frame <- function(y) {
    do.call(cbind, Map(function(k) {y[[k]]$c3}, seq_len(N2)))
}

y1 <- sample_tc()
d1 <- c1_data_frame(y1)
d2 <- c2_data_frame(y1)
d3 <- c3_data_frame(y1)

## ----sec2_d1, fig.cap = "A sample path of `d1`.", fig.height = 4--------------
plot(d1[,1], type = "l", ylim = c(-5, 5), xlab = "t", ylab = "")

## ----sec2_d2, fig.cap = "A sample path of `d2`.", fig.height = 4--------------
plot(d2[,1], type = "l", ylim = c(-5, 5), xlab = "t", ylab = "")

## ----sec2_d3, fig.cap = "A sample path of `d3`.", fig.height = 4--------------
plot(d3[,1], type = "l", ylim = c(-5, 5), xlab = "t", ylab = "")

## ----sec2_spec1, fig.cap = "The spectrum estimation of C1.", fig.height = 4----
spectrum(d1)

## ----sec2_spec2, fig.cap = "The spectrum estimation of C2.", fig.height = 4----
spectrum(d2)

## ----sec2_spec3, fig.cap = "The spectrum estimation of C3.", fig.height = 4----
spectrum(d3)

## ----sec2_heatmap_bc3, fig.cap = "`d3`'s estimated bicoherence.", fig.height = 4----
bc3 <- bicoherence(d3)
heatmap_bicoherence(bc3)

## ----sec2_heatmap_cb123, fig.cap = "The estimated cross-bicoherence between C1, C2, and C3.", fig.height = 6----
cb123 <- cross_bicoherence(d1, d2, d3)

heatmap_cross_bicoherence <- function(cb) {
    ggplot(cb) +
        geom_raster(aes(f1, f2, fill = value)) +
        coord_fixed() +
        scale_alpha(guide = "none")
}

heatmap_cross_bicoherence(cb123)

## ----sec2_heatmap_cb312, fig.cap = "The estimated cross-bicoherencde between C3, C1, and C2.", fig.height = 6----
cb312 <- cross_bicoherence(d3, d1, d2)
heatmap_cross_bicoherence(cb312)

## ----sec3_weak_coupling, fig.cap = "The case of weak coupling.", fig.height = 6----
Qcoef <- 0.01

y2 <- sample_tc()
cb2 <- cross_bicoherence(c1_data_frame(y2), c2_data_frame(y2), c3_data_frame(y2))
heatmap_cross_bicoherence(cb2)

## ----sec3_too_high_frequency, fig.cap = "The case of nearly-Nyquist frequency of `f1`.", fig.height = 6----
Fcoef1 <- pi - 0.1
Fcoef2 <- 2.3
Fcoef3 <- 1.5
Qcoef <- 0.3

y3 <- sample_tc()
cb3 <- cross_bicoherence(c1_data_frame(y3), c2_data_frame(y3), c3_data_frame(y3))
heatmap_cross_bicoherence(cb3)

## ----sec3_undersampling, fig.cap = "The case of insufficient sample.", fig.height = 6----
Fcoef1 <- 1.2
Fcoef2 <- 0.7
Fcoef3 <- 0.8
Qcoef <- 0.3
N2 <- 10

y4 <- sample_tc()
cb4 <- cross_bicoherence(c1_data_frame(y4), c2_data_frame(y4), c3_data_frame(y4))
heatmap_cross_bicoherence(cb4)

