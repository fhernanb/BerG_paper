
library(DiscreteDists)


# Figure pmf --------------------------------------------------------------

pdf("Figs/pmf.pdf", height=4, width=12)

# Parameters for the plots
mu_values_list <- list(c(0.8, 1, 1.3), c(0.3, 0.5, 1), c(2, 2.6, 3.3))
sigma_values <- c(0.5, 1, 3)
y_values <- 0:15

colores <- c("#fe4a49", "#fed766", "#009fb7")

# Set up the plotting area
par(mfrow=c(1, 3)) # Create a 3-row layout for the plots

# Loop through each combination of sigma and mu values
for (i in seq_along(sigma_values)) {
  sigma <- sigma_values[i]
  mu_values <- mu_values_list[[i]]
  
  # Initialize an empty plot
  plot(NULL, xlim=c(min(y_values), max(y_values)), ylim=c(0, 0.8),
       xlab="y", ylab="P(X=x)", 
       main=bquote(paste(sigma, "=", .(sigma)))
       )
  
  noise <- -0.2
  # Add lines for each mu value
  for (j in 1:3) {
    noise <- noise + 0.1
    mu <- mu_values[j]
    pmf_values <- sapply(y_values, dBerG, mu=mu, sigma=sigma)
    lines(y_values+noise, pmf_values, type="h", pch=16, lwd=4,
          col=colores[j])
  }
  
  # Add a legend
  legend("topright", 
         legend=paste("μ =", mu_values),
         col=colores, lty=1, lwd=4, bty="n")
}

dev.off()

