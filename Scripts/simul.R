library(gamlss)
library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(DiscreteDists)

theme_set(theme_minimal())

set.seed(65534L)

sample_sizes <- c(50, 100, 200, 400)

simBerG <- \(n){
  x1        <- runif(n)
  x2        <- runif(n)
  X         <- model.matrix(~ x1 + x2)
  beta      <- c(1, 1.2, .2)
  eta.mu    <- X%*%beta
  mu        <- exp(eta.mu)
  z1        <- runif(n)
  z2        <- runif(n)
  Z         <- model.matrix(~ z1 + z2)
  gamma     <- c(2, 1.5, 1.5)
  eta.sigma <- Z%*%gamma
  sigma     <- exp(eta.sigma)
  
  coresimBerG <- \(r){
    data.frame(x1, x2, z1, z2) |>
      mutate(y = rBerG(n, mu = mu, sigma = sigma)) |>
      {\(.data) gamlss(y ~ x1 + x2,
                       sigma.formula = ~ z1 + z2,
                       family        = BerG,
                       data          = .data,
                       trace         = FALSE) |>
          {\(.x) c(.x$mu.coefficients,
                   sigma = .x$sigma.coefficients)}()
      }()
  }
  
  rate <- rate_backoff(pause_base = 0, pause_min = 0, max_times = 15)
  result <- seq_len(1e4) |>
    map(insistently(coresimBerG, rate = rate), .progress = TRUE) |>
    reduce(bind_rows) |>
    mutate(n = n, .before = 1)
  
  return(result)
}



rate <- rate_backoff(pause_base = 0, pause_min = 0, max_times = 15)
result <- sample_sizes |>
  map(insistently(simBerG, rate = rate), .progress = FALSE) |>
  reduce(bind_rows)

result <- readRDS(file="result.rds")

result |>
  group_by(n) |>
  summarise(across(`(Intercept)`:sigma.z2, 
                   list(mean, sd), 
                   .names = "{.col}_{.fn}")
            ) |>
  pivot_longer(!n,
               names_to = c("parameters", ".value"),
               names_sep = "_") |>
  rename("mean" = `1`, "sd" = `2`) |>
  mutate(bias     = mean - c(1, 1.2, .2, 2, 1.5, 1.5),
         mse.sqrd = sqrt((1e4 - 1) * sd^2 / 1e4 + bias^2),
         .before = sd) |>
  print(n = Inf)

density.plots <- \(.n){
  p1 <- result |>
    filter(n == .n) |>
    ggplot(aes(x = `(Intercept)`, y = after_stat(density))) +
    geom_density() +
    geom_vline(xintercept = 1, lty = 2) +
    xlab(expression(beta[0]))
  
  p2 <- result |>
    filter(n == .n) |>
    ggplot(aes(x = x1, y = after_stat(density))) +
    geom_density() +
    geom_vline(xintercept = 1.2, lty = 2) +
    xlab(expression(beta[1]))
  
  p3 <- result |>
    filter(n == .n) |>
    ggplot(aes(x = x2, y = after_stat(density))) +
    geom_density() +
    geom_vline(xintercept = .2, lty = 2) +
    xlab(expression(beta[2]))
  
  p4 <- result |>
    filter(n == .n) |>
    ggplot(aes(x = `sigma.(Intercept)`, y = after_stat(density))) +
    geom_density() +
    geom_vline(xintercept = 2, lty = 2) +
    xlab(expression(gamma[0]))
  
  p5 <- result |>
    filter(n == .n) |>
    ggplot(aes(x = sigma.z1, y = after_stat(density))) +
    geom_density() +
    geom_vline(xintercept = 1.5, lty = 2) +
    xlab(expression(gamma[1]))
  
  p6 <- result |>
    filter(n == .n) |>
    ggplot(aes(x = sigma.z2, y = after_stat(density))) +
    geom_density() +
    geom_vline(xintercept = 1.5, lty = 2) +
    xlab(expression(gamma[2]))
  
  g <- (p1 + p2 + p3)/(p4 + p5 + p6) + plot_annotation(title = paste("n =", .n))
  #print(class(g))
  ggsave(file=paste0("Figs/res_simul_", .n), plot=g, device = "pdf")
  #plot(g)
  
}

sample_sizes |>
  map(density.plots)

# Graphs for beta^ versus n

res <- result %>% 
  drop_na() %>% 
  group_by(n) %>% 
  summarise(mean_b0=mean(`(Intercept)`),
            mean_b1=mean(x1),
            mean_b2=mean(x2),
            mean_g0=mean(`sigma.(Intercept)`),
            mean_g1=mean(sigma.z1),
            mean_g2=mean(sigma.z2),
            mse_b0=mean((1.0 - `(Intercept)`)^2), 
            mse_b1=mean((1.2 - x1)^2),
            mse_b2=mean((0.2 - x2)^2),
            mse_g0=mean((2.0 - `sigma.(Intercept)`)^2), 
            mse_g1=mean((1.5 - sigma.z1)^2),
            mse_g2=mean((1.5 - sigma.z2)^2)
            )

res

library(gridExtra)

# MSE -----------------------------------------------------
p1 <- ggplot(data=res, aes(x=n, y=mse_b0)) + 
  geom_line() + 
  #theme_bw() +
  labs(x="n", y=expression(MSE~hat(beta)[0]))

p2 <- ggplot(data=res, aes(x=n, y=mse_b1)) + 
  geom_line() + 
  #theme_bw() +
  labs(x="n", y=expression(MSE~hat(beta)[1]))

p3 <- ggplot(data=res, aes(x=n, y=mse_b2)) + 
  geom_line() + 
  #theme_bw() +
  labs(x="n", y=expression(MSE~hat(beta)[2]))

p4 <- ggplot(data=res, aes(x=n, y=mse_g0)) + 
  geom_line() + 
  #theme_bw() +
  labs(x="n", y=expression(MSE~hat(gamma)[0]))

p5 <- ggplot(data=res, aes(x=n, y=mse_g1)) + 
  geom_line() + 
  #theme_bw() +
  labs(x="n", y=expression(MSE~hat(gamma)[1]))

p6 <- ggplot(data=res, aes(x=n, y=mse_g2)) + 
  geom_line() + 
  #theme_bw() +
  labs(x="n", y=expression(MSE~hat(gamma)[2]))

mse <- grid.arrange(p1, p2, p3, p4, p5, p6, nrow=2, ncol=3)
mse
ggsave(filename="Figs/mse.pdf", 
       plot=mse, 
       width=12, height=7)
