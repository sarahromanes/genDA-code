

                                                  ## Bernoulli Simulation ##

# Load Required Packages
library(genDA)
library(gllvm)
library(vegan)
library(ggplot2)
library(doSNOW)
library(tcltk)

# Set Up Parallel Process

TRIALS <-  500

cl <- makeSOCKcluster(40)
registerDoSNOW(cl)

pb <- txtProgressBar(max=TRIALS, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

results.bern <- foreach(i = 1:TRIALS,
                      .combine = rbind,
                      .packages = c("genDA", "gllvm", "mvtnorm", "vegan"),
                      .options.snow = opts) %dopar%{
                        
                        n_vals <- c(100, 50, 50, 50, 500)
                        m_vals <- c(50, 100, 500, 1000, 1000)
                        
                        
                        res <- data.frame(value = double(),
                                          parameter = character(),
                                          method = character(),
                                          metric = character(),
                                          n = integer(),
                                          m = integer(),
                                          stringsAsFactors = FALSE)
                        
                        for(s in 1:length(m_vals)){
                          
                          set.seed(i)
                          
                          n <- n_vals[s]
                          m <- m_vals[s]
                            
                            ## Simulate Data ## 
                            
                            m1 <- c(-2, 2)
                            m2 <- c(0, -1)
                            m3 <- c(0,  1)
                            
                            sigma <- diag(c(1,1))
                            
                            u1 <- mvtnorm::rmvnorm(0.5*n, m1, sigma)
                            u2 <- mvtnorm::rmvnorm(0.3*n, m2, sigma)
                            u3 <- mvtnorm::rmvnorm(0.2*n, m3, sigma)
                            
                            mU <- rbind(u1,u2,u3)
                            
                            mLambda <- matrix(0, m, 2)
                            mLambda[,1] <- seq(2, -2, length.out = m)
                            mLambda[,2] <- seq(1, -1, length.out = m)
                            mL <- t(mLambda)

                            vbeta0 <- rep(-1, m)
                            mB <- matrix(rep(1,m), 1, m)
                            
                            vc <- rbinom(n,1,0.5)
                            mX <- matrix(vc, n, 1)
                            
                            mEta.fixed  <- matrix(1,n,1)%*%vbeta0 + mX%*%mB
                            mEta.latent <- mU%*%mL
                            
                            mEta <- mEta.fixed + mEta.latent
                            
                            mY <- matrix(0,n,m)
                            for (j in 1:m) {
                              mY[,j] <- rbinom(n,1,1/(1+exp(-mEta[,j])))
                            }
                            
                            t1 <- system.time({
                              res_gllvm <- gllvm(y=mY, X= data.frame(mX), num.lv=2,family="binomial", method = "VA", sd.errors = TRUE)
                            })
                            
                            t2 <- system.time({
                              res_genDA <- genDA(y=mY, class = as.factor(mX), num.lv = 2, family="binomial", standard.errors = TRUE) 
                            })
      
                            
                            q <- qnorm(0.975)
                            
                   #   ------- ## mU ## ------- #
                            
                            # Procrustes error for latent variables
                            p1a <- procrustes(mU, res_gllvm$lvs)$ss
                            p1b <- procrustes(mU, res_genDA$lvs)$ss
                           
                            
                  #  ------- ## mL ## ------#
                            
                            # Procrustes error for mLambda
                            p2a <- procrustes(mLambda, res_gllvm$params$theta )$ss
                            p2b <- procrustes(mLambda, t(res_genDA$params$mL) )$ss
                            
                    # ------- ## mB ## ------  #  
                            
                            # Bias for mB
                            b3a <- mean(t(mB) - res_gllvm$params$Xcoef)
                            b3b <- mean(t(mB) - res_genDA$params$Xcoef)
                            
                            # RMSE for mB
                            r3a <- sqrt(mean((t(mB) - res_gllvm$params$Xcoef)^2))
                            r3b <- sqrt(mean((t(mB) - res_genDA$params$Xcoef)^2))
                            
                            # CI_width for mB
                            u_mB_gllvm = res_gllvm$params$Xcoef + res_gllvm$sd$Xcoef*q
                            l_mB_gllvm = res_gllvm$params$Xcoef - res_gllvm$sd$Xcoef*q
                            
                            u_mB_genDA = res_genDA$params$Xcoef +t(res_genDA$sd$Xcoef)*q
                            l_mB_genDA = res_genDA$params$Xcoef -t(res_genDA$sd$Xcoef)*q
                            
                            CI_3a <- mean(abs((u_mB_gllvm- l_mB_gllvm)))
                            CI_3b <- mean(abs((u_mB_genDA- l_mB_genDA)))
                            
                            # Coverage for mB
                            CO_3a <- mean( t(mB) >= l_mB_gllvm & t(mB) <= u_mB_gllvm)
                            CO_3b <- mean( t(mB) >= l_mB_genDA & t(mB) <= u_mB_genDA)
                            
                    # ------- ## vbeta0 ## ------ #                       
                            
                            # Bias for vbeta0
                            b4a <- mean(vbeta0 - res_gllvm$params$beta0)
                            b4b <- mean(vbeta0 - res_genDA$params$beta0)
                            
                            # RMSE for vbeta0
                            r4a <-  sqrt(mean((vbeta0 - res_gllvm$params$beta0)^2))
                            r4b <-  sqrt(mean((vbeta0 - res_genDA$params$beta0)^2))
                            
                            # CI_width for vbeta0
                            u_b0_gllvm = res_gllvm$params$beta0 + res_gllvm$sd$beta0*q
                            l_b0_gllvm = res_gllvm$params$beta0 - res_gllvm$sd$beta0*q
                            
                            u_b0_genDA = res_genDA$params$beta0 +res_genDA$sd$vbeta0*q
                            l_b0_genDA = res_genDA$params$beta0 -res_genDA$sd$vbeta0*q
                            
                            CI_4a <- mean(abs((u_b0_gllvm- l_b0_gllvm)))
                            CI_4b <- mean(abs((u_b0_genDA- l_b0_genDA)))
                            
                            # Coverage for vbeta0
                            CO_4a <- mean( vbeta0 >= l_b0_gllvm & vbeta0 <= u_b0_gllvm)
                            CO_4b <- mean( vbeta0 >= l_b0_genDA & vbeta0 <= u_b0_genDA)
                            
                            
                   # ----------------------------------- #   
                            
                   vals_gllvm  <- c(p1a, p2a, b3a, r3a, CI_3a, CO_3a, b4a, r4a, CI_4a, CO_4a, round(t1[3], 4))
                   vals_genDA  <- c(p1b, p2b, b3b, r3b, CI_3b, CO_3b, b4b, r4b, CI_4b, CO_4b, round(t2[3], 4))
                            
                   methods <- rep(c("gllvm", "genDA"), each=length(vals_gllvm))
                   params <- rep(c("mU", "mLambda", rep("mB", 4), rep("vbeta0", 4), "time"),2)
                   metric <- rep(c(rep("procrustes", 2), rep(c("Bias", "RMSE", "CI Width", "Coverage"), 2), "s"), 2)
                            
                   res_new <- data.frame(value = c(vals_gllvm, vals_genDA), parameter = params, method = methods, metric = metric, n = n, m = m)
                   
                   res <- rbind(res, res_new)
                            
                          
                    }
                        
                        return(res)
                      }

# Stop Distribution
stopCluster(cl)

# Save results to RData file
save(results.bern, file="Bernoulli_full.RData")