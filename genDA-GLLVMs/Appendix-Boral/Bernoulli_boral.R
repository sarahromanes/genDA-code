

                                                  ## Bernoulli Simulation ##

# Load Required Packages
library(vegan)
library(ggplot2)
library(doSNOW)
library(tcltk)
library(boral)

# Set Up Parallel Process

TRIALS <-  100


#This had to be run on my Windows Laptop -could not get it to parallelise on Debian system at USYD #

cl <- makeSOCKcluster(7)
registerDoSNOW(cl)

pb <- txtProgressBar(max=TRIALS, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

results.boral.bern <- foreach(i = 1:TRIALS,
                      .combine = rbind,
                      .packages = c("boral", "mvtnorm", "vegan"),
                      .options.snow = opts) %dopar%{
                        
                         n_vals <- c(100, 50)
                         m_vals <- c(50, 100)
                      
                        
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
                            
                            vtau <-rep(0, n)
                            vbeta0 <- rep(-1, m)
                            mB <- matrix(rep(1,m), 1, m)
                            
                            vc <- rbinom(n,1,0.5)
                            mX <- matrix(vc, n, 1)
                            
                            mEta.fixed  <- matrix(vtau)%*%matrix(1,1,m) + matrix(1,n,1)%*%vbeta0 + mX%*%mB
                            mEta.latent <- mU%*%mL
                            
                            mEta <- mEta.fixed + mEta.latent
                            
                            mY <- matrix(0,n,m)
                            for (j in 1:m) {
                              mY[,j] <- rbinom(n,1,1/(1+exp(-mEta[,j])))
                            }
                            
                            t1 <- system.time({
                              res_boral <- boral(y=mY, X= mX,family="binomial",lv.control = list(num.lv = 2))
                            })
                    
                            q <- qnorm(0.975)
                            
                   #   ------- ## mU ## ------- #
                            
                            # Procrustes error for latent variables
                            p1a <- procrustes(mU, res_boral$lv.mean)$ss
        
                  #  ------- ## mL ## ------#
                            
                            # Procrustes error for mLambda
                            p2a <- procrustes(mLambda, res_boral$lv.coefs.mean[, c(2,3)] )$ss
                 
                            
                    # ------- ## mB ## ------  #  
                            
                            # Bias for mB
                            b3a <- mean(t(mB) - res_boral$X.coefs.mean)
                         
                            # RMSE for mB
                            r3a <- sqrt(mean((t(mB) - res_boral$X.coefs.mean)^2))
          
                            
                            # CI_width for mB
                            u_mB_boral = res_boral$X.coefs.mean + res_boral$X.coefs.sd*q
                            l_mB_boral = res_boral$X.coefs.mean - res_boral$X.coefs.sd*q
                            
                            CI_3a <- mean(abs((u_mB_boral- l_mB_boral)))
                            
                            # Coverage for mB
                            CO_3a <- mean( t(mB) >= l_mB_boral & t(mB) <= u_mB_boral)
                            
                    # ------- ## vbeta0 ## ------ #                       
                            
                            # Bias for vbeta0
                            b4a <- mean(vbeta0 - res_boral$lv.coefs.mean[,1])
             
                            
                            # RMSE for vbeta0
                            r4a <-  sqrt(mean((vbeta0 - res_boral$lv.coefs.mean[,1])^2))
                        
                            
                            # CI_width for vbeta0
                            u_b0_boral = res_boral$lv.coefs.mean[,1] + res_boral$lv.coefs.sd[,1]*q
                            l_b0_boral = res_boral$lv.coefs.mean[,1]- res_boral$lv.coefs.sd[,1]*q
                            
                            
                            CI_4a <- mean(abs((u_b0_boral- l_b0_boral)))
                            
                            # Coverage for vbeta0
                            CO_4a <- mean( vbeta0 >= l_b0_boral & vbeta0 <= u_b0_boral)
                            
                            
                   # ----------------------------------- #   
                            
                   vals_boral  <- c(p1a, p2a, b3a, r3a, CI_3a, CO_3a, b4a, r4a, CI_4a, CO_4a, round(t1[3], 4))
                            
                   methods <- rep(c("boral"), each=length(vals_boral))
                   params <- rep(c("mU", "mLambda", rep("mB", 4), rep("vbeta0", 4), "time"),1)
                   metric <- rep(c(rep("procrustes", 2), rep(c("Bias", "RMSE", "CI Width", "Coverage"), 2), "s"), 1)
                            
                   res_new <- data.frame(value = c(vals_boral), parameter = params, method = methods, metric = metric, n = n, m = m)
                   
                   res <- rbind(res, res_new)
                            
                          
                    }
                        
                        return(res)
                      }

# Stop Distribution
stopCluster(cl)

# Save results to RData file
save(results.boral.bern, file="Bernoulli_boral.RData")