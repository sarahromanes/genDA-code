

                                                  ## Mixed Simulation ##

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

results.mixed.genDA <- foreach(i = 1:TRIALS,
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
                            m3 <- c(0.5, -1.5)
                            
                            sigma <- 0.5*diag(c(1,1))
                            
                            u1 <- mvtnorm::rmvnorm(0.4*n, m1, sigma)
                            u3 <- mvtnorm::rmvnorm(0.6*n, m3, sigma)
                            
                            mU <- rbind(u1,u3)
                            
                            mLambda <- matrix(0, m, 2)
                            mLambda[,1] <- runif(m, -2, 2)
                            mLambda[,2] <- runif(m, -2, 2)
                            mL <- t(mLambda)
                            
                            vbeta0 <- rep(-1, m)
                            mB <- matrix(rep(1,m), 1, m)
                            vphi <- rep(1, m)
                            
                            vc <- rbinom(n,1,0.5)
                            mX <- matrix(vc, n, 1)
                            
                            
                            mEta.fixed  <- matrix(1,n,1)%*%vbeta0 +  mX%*%mB
                            mEta.latent <- mU%*%mL
                            
                            mEta <- mEta.fixed + mEta.latent
                            
                            splits = ceiling(m/3)
                            f1 <- rep("gaussian", splits)
                            f2 <- rep("poisson", splits)
                            f3 <- rep("binomial", (m-2*splits))
                            
                            f <- c(f1,f2,f3)
                            
                            mY <- matrix(0,n,m)
                            for (j in 1:m) {
                              
                              if(f[j]=="gaussian"){
                                mY[,j] <- rnorm(n,mEta[,j], vphi[j])
                              } else if (f[j]=="poisson"){
                                mY[,j] <- rpois(n,exp(mEta[,j]))
                              } else {
                                mY[,j] <-  rbinom(n,1,1/(1+exp(-mEta[,j])))
                              }
                            }
                            
                            
                            t2 <- system.time({
                              res_genDA <- genDA(y = mY, class = mX, num.lv =2, family = f, standard.errors = TRUE)
                            })
      
                            
                            q <- qnorm(0.975)
                            
                   #   ------- ## mU ## ------- #
                            
                            # Procrustes error for latent variables
                            p1a <- procrustes(mU, res_genDA$lvs)$ss
                     
                            
                  #  ------- ## mL ## ------#
                            
                            # Procrustes error for mLambda
                         
                            p2b <- procrustes(mLambda, t(res_genDA$params$mL) )$ss
                            
                    # ------- ## mB ## ------  #  
                            
                            # Bias for mB
               
                            b3b <- mean(t(mB) - res_genDA$params$Xcoef)
                            
                            # RMSE for mB
                        
                            r3b <- sqrt(mean((t(mB) - res_genDA$params$Xcoef)^2))
                            
                            # CI_width for mB
                    
                            
                            u_mB_genDA = res_genDA$params$Xcoef +t(res_genDA$sd$Xcoef)*q
                            l_mB_genDA = res_genDA$params$Xcoef -t(res_genDA$sd$Xcoef)*q
                            
                           
                            CI_3b <- mean(abs((u_mB_genDA- l_mB_genDA)))
                            
                            # Coverage for mB
                  
                            CO_3b <- mean( t(mB) >= l_mB_genDA & t(mB) <= u_mB_genDA)
                            
                    # ------- ## vbeta0 ## ------ #                       
                            
                            # Bias for vbeta0
                      
                            b4b <- mean(vbeta0 - res_genDA$params$beta0)
                            
                            # RMSE for vbeta0
                          
                            r4b <-  sqrt(mean((vbeta0 - res_genDA$params$beta0)^2))
                            
                            # CI_width for vbeta0
                  
                            u_b0_genDA = res_genDA$params$beta0 +res_genDA$sd$vbeta0*q
                            l_b0_genDA = res_genDA$params$beta0 -res_genDA$sd$vbeta0*q
                            
                         
                            CI_4b <- mean(abs((u_b0_genDA- l_b0_genDA)))
                            
                            # Coverage for vbeta0
                   
                            CO_4b <- mean( vbeta0 >= l_b0_genDA & vbeta0 <= u_b0_genDA)
                            
                     #------- ## vphi ## ------  #                           
                            
                            # Bias for vphi
                      
                            b5b <- mean(vphi[which(f=="gaussian")] - res_genDA$params$phi[which(f=="gaussian")] )
                            
                            # RMSE for vphi
                        
                            r5b <- sqrt(mean((vphi[which(f=="gaussian")] - res_genDA$params$phi[which(f=="gaussian")])^2))
                            
                            # CI_width for vphi
                          
                            u_phi_genDA = res_genDA$params$phi[which(f=="gaussian")] +res_genDA$sd$phi[which(f=="gaussian")]*q
                            l_phi_genDA = res_genDA$params$phi[which(f=="gaussian")] -res_genDA$sd$phi[which(f=="gaussian")]*q
                            
                          
                            CI_5b <- mean(abs((u_phi_genDA- l_phi_genDA)))
                            
                            # Coverage for vphi
                         
                            CO_5b <- mean( vphi[which(f=="gaussian")] >= l_phi_genDA & vphi[which(f=="gaussian")] <= u_phi_genDA)
                            
                   # ----------------------------------- #   
                            
              
                   vals_genDA  <- c(p1a, p2b, b3b, r3b, CI_3b, CO_3b, b4b, r4b, CI_4b, CO_4b, b5b, r5b, CI_5b, CO_5b, round(t2[3], 4))
                            
                   methods <- rep(c("genDA"), each=length(vals_genDA))
                   params <- rep(c("mU", "mLambda", rep("mB", 4), rep("vbeta0", 4), rep("vphi", 4), "time"),1)
                   metric <- rep(c(rep("procrustes", 2), rep(c("Bias", "RMSE", "CI Width", "Coverage"), 3), "s"), 1)
                            
                   res_new <- data.frame(value = c(vals_genDA), parameter = params, method = methods, metric = metric, n = n, m = m)
                   
                   res <- rbind(res, res_new)

                    }
                        
                  return(res)
              }

# Stop Distribution
stopCluster(cl)
close(pb)

# Save results to RData file
save(results.mixed.genDA, file="mixed_genDA.RData")