

                                                  ## NB Simulation ##

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

# Run the simulation 

results.NB <- foreach(i = 1:TRIALS,
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
                            
                            m1 <- c(-1, 1)
                            m2 <- c(2, 1.5)
                            m3 <- c(0.5, -1.5)
                            
                            sigma <- 0.5*diag(c(1,1))
                            
                            u1 <- mvtnorm::rmvnorm(0.4*n, m1, sigma)
                            u2 <- mvtnorm::rmvnorm(0.3*n, m2, sigma)
                            u3 <- mvtnorm::rmvnorm(0.3*n, m3, sigma)
                            
                            mU <- rbind(u1,u2,u3)
                            
                            mLambda <- matrix(0, m, 2)
                            mLambda[,1] <- runif(m, -2, 2)
                            mLambda[,2] <- runif(m, -2, 2)
                            mL <- t(mLambda)
                            
                            vbeta0 <- runif(m, -1, 1)
                            vphi <- rep(1, m)
                            
                            mEta.fixed  <- matrix(1,n,1)%*%vbeta0
                            mEta.latent <- mU%*%mL
                            
                            mEta <- mEta.fixed + mEta.latent
                            
                            mY <- matrix(0,n,m)
                            for (j in 1:m) {
                              mY[,j] <-  rnbinom(n, mu = exp(mEta[,j]), size = vphi[j])
                            }
                            
                            t1 <- system.time({
                              res_gllvm <- gllvm::gllvm(y=mY,num.lv=2,family="negative.binomial", method = "VA", row.eff = FALSE, sd.errors = TRUE)
                            })
                            
                            t2 <- system.time({
                              res_genDA <- genDA(y = mY, num.lv =2, family = "negative-binomial", standard.errors = TRUE)
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
                            
                            u_b0_genDA = res_genDA$params$beta0 +res_genDA$sd$beta0*q
                            l_b0_genDA = res_genDA$params$beta0 -res_genDA$sd$beta0*q
                            
                            CI_4a <- mean(abs((u_b0_gllvm- l_b0_gllvm)))
                            CI_4b <- mean(abs((u_b0_genDA- l_b0_genDA)))
                            
                            # Coverage for vbeta0
                            CO_4a <- mean( vbeta0 >= l_b0_gllvm & vbeta0 <= u_b0_gllvm)
                            CO_4b <- mean( vbeta0 >= l_b0_genDA & vbeta0 <= u_b0_genDA)
                            
                     #------- ## vphi ## ------  #                           
                            
                            # Bias for vphi
                            b5a <- mean(vphi - res_gllvm$params$inv.phi)
                            b5b <- mean(vphi - res_genDA$params$phi)
                            
                            # RMSE for vphi
                            r5a <- sqrt(mean((vphi - res_gllvm$params$inv.phi)^2))
                            r5b <- sqrt(mean((vphi - res_genDA$params$phi)^2))
                            
                            # CI_width for vphi
                            u_phi_gllvm = res_gllvm$params$inv.phi + res_gllvm$sd$inv.phi*q
                            l_phi_gllvm = res_gllvm$params$inv.phi - res_gllvm$sd$inv.phi*q
                            
                            u_phi_genDA = res_genDA$params$phi +res_genDA$sd$phi*q
                            l_phi_genDA = res_genDA$params$phi -res_genDA$sd$phi*q
                            
                            CI_5a <- mean(abs((u_phi_gllvm- l_phi_gllvm)))
                            CI_5b <- mean(abs((u_phi_genDA- l_phi_genDA)))
                            
                            # Coverage for vphi
                            CO_5a <- mean( vphi >= l_phi_gllvm & vphi <= u_phi_gllvm)
                            CO_5b <- mean( vphi >= l_phi_genDA & vphi <= u_phi_genDA)
                            
                   # ----------------------------------- #   
                            
                   vals_gllvm  <- c(p1a, p2a, b4a, r4a, CI_4a, CO_4a, b5a, r5a, CI_5a, CO_5a, round(t1[3], 4))
                   vals_genDA  <- c(p1b, p2b, b4b, r4b, CI_4b, CO_4b, b5b, r5b, CI_5b, CO_5b, round(t2[3], 4))
                            
                   methods <- rep(c("gllvm", "genDA"), each=length(vals_gllvm))
                   params <- rep(c("mU", "mLambda", rep("vbeta0", 4), rep("vphi", 4), "time"),2)
                   metric <- rep(c(rep("procrustes", 2), rep(c("Bias", "RMSE", "CI Width", "Coverage"), 2), "s"), 2)
                            
                   res_new <- data.frame(value = c(vals_gllvm, vals_genDA), parameter = params, method = methods, metric = metric, n = n, m = m)
                   
                   res <- rbind(res, res_new)

                    }
                        
                  return(res)
              }

# Stop Distribution
stopCluster(cl)
close(pb)

# Save results to RData file
save(results.NB, file="NB_full.RData")