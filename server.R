library(shiny)
library(ggplot2)
library(Cairo)
library(tidyr)
library(dplyr)

shinyServer(function(input, output) {
  
  clamp <- function(x, a, b) max(a, min(x, b))
  prony.moments <- function(x, a) {
    n <- length(a)
    sapply(0:(2 * n - 1), function(k) sum(a * x ^ k))
  }
  prony.solve <- function(m) {
    n <- length(m) / 2
    u.A <- outer(1:n, 1:n, function(i, j) m[i - j + n] * (-1) ^ (j + 1))
    u.B <- m[(n + 1):(2 * n)]
    u <- solve(u.A, u.B)
    x <- sort(polyroot(c(rev(u), 1) * sapply(1:(length(u)+1), function(i) (-1)^i)))
    a.A <- t(sapply(0:(n - 1), function(i) x ^ i))
    a.B <- m[1:n]
    a <- solve(a.A, a.B)
    list(x = x, a = a)
  }
  prony.tr2d.solve <- function(m, delta) {
    m0 <- m[1]; m1 <- m[2]; m2 <- m[3]
    M0 <- Mod(m0); M1 <- Mod(m1); M2 <- Mod(m2)
    
    Delta.cos <- clamp((2 * M1 ^ 2 - M0 ^ 2 - M2 ^ 2) / (2 * M0 ^ 2 - 2 * M1 ^ 2), -1, 1)
    Delta <- acos(Delta.cos)
    D <- max(M0 ^ 2 - 4 * (M1 ^ 2 - M0 ^ 2) ^ 2 / (M2 ^ 2 + 3 * M0 ^ 2 - 4 * M1 ^ 2), 0)
    a.r <- (M0 - sqrt(D)) / 2
    b.r <- (M0 + sqrt(D)) / 2
    
    p <- a.r * Delta.cos + b.r
    q <- a.r * sin(Delta)
    theta.cos <- clamp((Re(m1) * p + Im(m1) * q) / (p ^ 2 + q ^ 2), -1, 1)
    theta <- -acos(theta.cos)
    phi <- theta + Delta  
    mu.r <- phi / (-2 * pi * delta)
    nu.r <- theta / (-2 * pi * delta)
    x.r <- exp(-2i * pi * delta * mu.r)
    y.r <- exp(-2i * pi * delta * nu.r)
    list(x = c(x.r, y.r), a = c(a.r, b.r))
  }
  generate.noise <- function(length, mod) {
    noise <- runif(length * 2) * mod
    noise <- noise[1:length] + noise[(length + 1):(length * 2)] * 1i
    noise
  }
  generate.noise.seed <- function(seed, length, mod) {
    set.seed(seed)
    generate.noise(length, mod)
  }
  
  output$mytable = renderDataTable({
    n <- 2
    mu <- c(input$mu1, input$mu2)
    a <- c(input$a1, input$a2)
    delta <- input$delta
    x <- sort(exp(-2i * pi * delta * mu))
    
    m <- prony.moments(x, a)
    noise <- generate.noise.seed(52777, 2 * n, input$noise)
    m.noise <- m + noise
    reconstruct <- prony.solve(m.noise)
    reconstruct.tr <- prony.tr2d.solve(m.noise, delta)
    pretty <- function(data) prettyNum(data, preserve.width = "common", drop0trailing = T, digits = 3, format = "f")
    
    df <- data.frame(
      name = c("Original system", "General", "Trigonometrical"),
      x1 = pretty(c(x[1], reconstruct$x[1], reconstruct.tr$x[1])),
      x2 = pretty(c(x[2], reconstruct$x[2], reconstruct.tr$x[2])),
      a1 = pretty(c(a[1], reconstruct$a[1], reconstruct.tr$a[1])),
      a2 = pretty(c(a[2], reconstruct$a[2], reconstruct.tr$a[2])),
      m0 = pretty(c(m[1], m.noise[1], m.noise[1])),
      m1 = pretty(c(m[2], m.noise[2], m.noise[2])),
      m2 = pretty(c(m[3], m.noise[3], m.noise[3])),
      m3 = pretty(c(m[4], m.noise[4], m.noise[4]))
    )
    
    df
  }, options = list(paging = FALSE, searching = FALSE))
  
  output$plot <- renderPlot({
    gauss <- function(a, x0, s2) { function(x) a * exp(-((x-x0)^2)/(2*s2))}
    substrRight <- function(x, n = 1){
      substr(x, nchar(x)-n+1, nchar(x))
    }
    
    mu <- c(input$mu1, input$mu2)
    a <- c(input$a1, input$a2)
    delta <- input$delta
    t <- seq(-5, 5, by = 0.01)
    f1 <- gauss(a[1], mu[1], 1)(t)
    f2 <- gauss(a[2], mu[2], 1)(t)
    f <- f1 + f2
    
    n <- 2
    x <- sort(exp(-2i * pi * delta * mu))
    m <- prony.moments(x, a)
    noise <- generate.noise.seed(52777, 2 * n, input$noise)
    m.noise <- m + noise
    reconstruct <- prony.solve(m.noise)
    reconstruct.tr <- prony.tr2d.solve(m.noise, delta)
    r.mu <- Arg(reconstruct$x) / (-2*pi*delta)
    r.mu.tr <- Arg(reconstruct.tr$x) / (-2*pi*delta)
    
    r.f1 <- gauss(Mod(reconstruct$a[1]), r.mu[1], 1)(t)
    r.f2 <- gauss(Mod(reconstruct$a[2]), r.mu[2], 1)(t)
    r.f <- r.f1 + r.f2
    r.f1.tr <- gauss(Mod(reconstruct.tr$a[1]), r.mu.tr[1], 1)(t)
    r.f2.tr <- gauss(Mod(reconstruct.tr$a[2]), r.mu.tr[2], 1)(t)
    r.f.tr <- r.f1.tr + r.f2.tr
    
    getType <- function(s) {
      if (substr(s, 1, 2) == "GO")
        "Original system"
      else if (substr(s, 1, 3) == "GRG")
        "General"
      else if (substr(s, 1, 3) == "GRT")
        "Trigonometrical"
      else
        "Unknown"
    }
    getType2 <- function(s) {
      s <- as.character(s)
      if (substrRight(s) == "1" || substrRight(s) == "2")
        "Signal"
      else
        "Composition"
    }
    getType <- Vectorize(getType)
    getType2 <- Vectorize(getType2)
    df <- data.frame(
      t, 
      GO1 = f1, GO2 = f2, GO = f, 
      GRG1 = r.f1, GRG2 = r.f2, GRG = r.f,
      GRT1 = r.f1.tr, GRT2 = r.f2.tr, GRT = r.f.tr) %>% 
      gather("name", "value", 2:10) %>%
      mutate(type = getType(name)) %>%
      mutate(type2 = getType2(name))
    plot <- ggplot(df, aes(t, value, group = name, colour = type, linetype = type2)) + geom_line()
    plot
  }, height=700)
}) 