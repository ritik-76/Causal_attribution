library(ggplot2)

alpha <- 1
beta <- 4
c.phi <- 0.005

h1.si <- function(si){
  a1.phi + b1.phi * log( (si+0.5)^c.phi / (1- (si+0.5)^c.phi))
}

b1.phi <- 0.25 * (1/c.phi) * (1- (0.25)^ c.phi)
a1.phi <- -0.25 - (b1.phi * log( (0.25)^c.phi / (1- (0.25)^c.phi)))

g1.phi <- function(phi){
  ( 1 + exp(-(phi - a1.phi)/b1.phi)) ^(-1/c.phi) -0.5
}

prior1.phi <- function(phi) {
  (gamma(alpha + beta) / (gamma(alpha) * gamma(beta) * b1.phi * c.phi)) *
    (g1.phi(phi)+0.5) ^ (alpha + c.phi) * (0.5 - g1.phi(phi))^(beta - 1) *
    exp( -(phi-a1.phi) / b1.phi)
}


df <- data.frame(x = -0.5 + rbeta(8000, 1, 4))

p.si.plot <- ggplot(df, aes(x = h1.si(x))) +
  geom_density( color = "black",size=1, adjust=1.5) +
  labs(x = expression(psi),
       y = expression(pi(psi)))+
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )




library(ggplot2)


df <- data.frame(xi = seq(-0.5, 0.5, length.out = 500))

df$y1 <- df$xi
df$y2 <- h1.si(df$xi)     

p.fun.plot <- ggplot(df, aes(x = xi)) +
  geom_line(aes(y = y2), color = "black", linewidth = 1) +
  geom_line(aes(y = y1), color = "red", linewidth = 1) +
  
  # reference lines
  #geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
  geom_vline(xintercept = -0.25, color = "black", linewidth = 0.8) +
  
  labs(x = expression(xi),
       y = expression(f(xi))) +
  
  theme_bw(base_size = 14) +
  
  # equal grid spacing
  scale_x_continuous(
    breaks = seq(-0.5, 0.5, by = 0.25),
    minor_breaks = seq(-0.5, 0.5, by = 0.125)
  ) +
  scale_y_continuous(
    breaks = seq(-1.5, 2, by = 0.5),
    minor_breaks = seq(-1.5, 2, by = 0.25)
  ) +
  
  # box border + grid styling
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90"),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16)
  )