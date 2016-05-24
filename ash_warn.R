options(warn = 2)
betahat <- c(-0.695, -0.497, -0.176, 0.115, -3.757, -0.752, 0.163, -0.308, 0.361)
sebetahat <- c(0.186, 0.221, 0.233, 0.436, 0.548, 0.307, 0.532, 0.997, 0.526)
ashr::ash.workhorse(betahat = betahat, sebetahat = sebetahat, df = 6) ## warnings
ashr::ash.workhorse(betahat = betahat, sebetahat = sebetahat, df = 7) ## works fine
