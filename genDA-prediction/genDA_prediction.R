

## Prediction Simulation ##

# Load Required Packages
library(ggplot2)
library(doSNOW)
library(tcltk)

# Set Up Parallel Process

TRIALS <-  250

cl <- makeSOCKcluster(45)
registerDoSNOW(cl)

pb <- txtProgressBar(max=TRIALS, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

# Run the simulation 

results.prediction <- foreach(i = 1:TRIALS,
                              .combine = rbind,
                              .packages = c("genDA", "MASS", "randomForest" ,"glmnet", "class", "e1071", "mvtnorm", "cvTools", "stats"),
                              .options.snow = opts) %dopar%{
                                
                                
                                ## Functions and necessary values
                                
                                n_vals <- c(100, 50, 50, 500)
                                m_vals <- c(50, 100, 1000, 1000)
                                sep_vals <- c(2, 1, 0.75, 0.5, 0.25, 0.1,0)
                                
                                covariance <- c("Equal", "Unequal")
                                
                                .vec2mat <- function(vy) {
                                  vy=as.numeric(vy)
                                  K <- length(unique(vy))
                                  n <- length(vy)
                                  mY <- matrix(0, nrow = n, ncol = K)
                                  for (i in 1:length(vy)) {
                                    val <- vy[i]
                                    mY[i, val] <- 1
                                  }
                                  return(mY)
                                }
                                
                                data_generation = function(n, m, k, i,sep){
                                  
                                  set.seed(i)
                                  
                                  ## Simulate Data ## 
                                  
                                  vphi <- rep(1, m)
                                  vtau <- rep(0, n)
                                  
                                  
                                  class = as.factor(sample(c(0,1), n, replace=TRUE))
                                  
                                  
                                  splits = ceiling(m/3)
                                  f1 <- rep("gaussian", splits)
                                  f2 <- rep("poisson", splits)
                                  f3 <- rep("binomial", (m-2*splits))
                                  
                                  f <- c(f1,f2,f3)
                                  f <- f[sample(m)]
                                  
                                  class <- sort(class)
                                  
                                  if(k == "Equal"){
                                    
                                    m1 <- c(0, 0)
                                    sigma <- diag(c(0.5,0.5))
                                    mU <- mvtnorm::rmvnorm(n, m1, sigma)
                                    
                                    
                                    mLambda <- matrix(0, m, 2)
                                    mLambda[,1] <- runif(m, -1, 1)
                                    mLambda[,2] <- runif(m, -1, 1)
                                    mL <- t(mLambda)
                                    
                                    mEta.latent <- mU%*%mL
                                  } else {
                                    
                                    mEta.latent <- c()
                                    for(c in 0:1){
                                      n_sub <- sum(class == c)
                                      
                                      m1 <- c(0, 0)
                                      sigma <- diag(c(0.5,0.5))
                                      mU <- mvtnorm::rmvnorm(n_sub, m1, sigma)
                                      
                                      mLambda <- matrix(0, m, 2)
                                      mLambda[,1] <- if(c==0){rep(c(1,1,1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1), m/25)} else if(c==1){rep(c(-1,-1,-1,-1,-1,1,1,1,1,1), m/10)} 
                                      mLambda[,2] <- if(c==0){rep(c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1), m/25)} else if(c==1){rep(c(0,0,0,0,0,1,1,1,1,1), m/10)} 
                                      
                                      
                                      mL <- t(mLambda)
                                    
                                      mEta.latent_sub <- mU%*%mL
                                      
                                      mEta.latent <- rbind(mEta.latent, mEta.latent_sub)
                                    }
                                  }
                                  
                                  mX <- as.matrix(.vec2mat(class)[,-1])
                                  mB <- matrix(sep, nrow = 1, ncol = m)
                                  
                                  vbeta0 <- rep(-sep, m)
                                  
                                  mEta.fixed  <- matrix(1,n,1)%*%vbeta0 +  mX%*%mB
                                  
                                  
                                  mEta <- mEta.fixed + mEta.latent
                                  
                                  
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
                                  
                                  samp <- sample(n)
                                  class <- class[samp]
                                  mY <- mY[samp, ]
                                  
                                  return(list(mY=mY, class=class, f=f))
                                }
                                
                                res <- data.frame(value = double(),
                                                  method = character(),
                                                  covariance = character(),
                                                  n = integer(),
                                                  m = integer(),
                                                  sep = double(), 
                                                  trial = integer(),
                                                  stringsAsFactors = FALSE)
                                
                                
                                V = 5
                                
                                ################## Begin Prediction Loop #####################
                                
                                t1 <- system.time({
                                  for(s in 1:length(n_vals)){
                                    
                                    for(l in 1:length(sep_vals)){
                                  
                                      for(k in covariance){
                                      set.seed(i)
                                      
                                      n <- n_vals[s]
                                      m <- m_vals[s]
                                      
                                      sep = sep_vals[l]
                                      
                                      data <- data_generation(n, m, k, i, sep)
                                      
                                      mY <- data$mY
                                      class = data$class
                                      f = data$f
                                      
                                      samp <- sample(n)
                                      
                                      cvSets <- cvFolds(n, V)
                                      errSet <- matrix(nrow=V, ncol=7, 0)
                                      
                                      for(j in 1:V){
                                        
                                        inds = which(cvSets$which==j)
                                        
                                        test.inds <- samp[inds]
                                        train.inds <- -samp[inds]
                                        
                                        class.test <- class[test.inds]
                                        class.train <- class[train.inds]
                                        mY.test <- mY[test.inds, ]
                                        mY.train <- mY[train.inds, ]
                                        
                                        data.train <- data.frame(class = class.train, mY.train)
                                        
                                        # genDA - common covariance
                                        
                                        res_genDA_LDA <- genDA(y = mY.train, class = class.train, num.lv =2, family = f, standard.errors = FALSE, common.covariance = TRUE)
                                        vals <- predict(res_genDA_LDA, newdata = mY.test)$class
                                        err_genDA_LDA <- sum(as.character(vals)!=as.character(class.test))
                                        
                                        # genDA - different covariance 
                                        
                                        res_genDA_QDA <- genDA(y = mY.train, class = class.train, num.lv =2, family = f, standard.errors = FALSE, common.covariance = FALSE)
                                        vals <- predict(res_genDA_QDA, newdata = mY.test)$class
                                        err_genDA_QDA <- sum(as.character(vals)!=as.character(class.test))
                                        
                                        # randomForest
                                        
                                        res_RF <- randomForest(mY.train, class.train)
                                        vals <- predict(res_RF, newdata = mY.test)
                                        err_RF <- sum(as.character(vals)!=as.character(class.test))
                                        
                                        
                                        # Test ordinary glm with no feature selection
                                        res_glm <- glm(class~., data = data.train,family=binomial)
                                        err_glm <- sum(round(predict(res_glm, data.frame(mY.test), type="response"))!=class.test)
                                        
                                        # glmnet - LASSO
                                        
                                        res_LC <- cv.glmnet(as.matrix(mY.train),as.factor(class.train),family="multinomial")
                                        vals <- predict(res_LC, newx = as.matrix(mY.test), s = "lambda.min", type = "class")
                                        err_LC <- sum(as.character(vals)!=as.character(class.test))
                                        
                                        
                                        # KNN (K=1)
                                        
                                        res_KNN <- knn(mY.train, mY.test, as.factor(class.train), k=1)
                                        err_KNN <- sum(as.character(res_KNN)!=as.character(class.test))
                                        
                                        
                                        # SVM (radial kernel)
                                        
                                        res_SVM <- svm(mY.train, as.factor(class.train), probability=FALSE)
                                        vals <-predict(res_SVM,mY.test, decision.values = TRUE, probability=FALSE)
                                        err_svm <- sum(as.character(vals)!=as.character(class.test))
                                        
                                        errSet[j,] <- c(err_genDA_LDA, err_genDA_QDA, err_RF, err_glm, err_LC, err_KNN, err_svm)
                                        
                                      }
                                      
                                      values <- apply(errSet,2, sum)/n
                                      methods <- c("genDA-common", "genDA-separate", "randomForest","glm", "LASSO", "KNN-1", "SVM")
                                      
                                      res_new <- data.frame(value = values, covariance = k,  n = n, m = m, method = methods, sep = sep_vals[l], trial = i)
                                      
                                      print(res_new[1:2,])
                                      res <- rbind(res, res_new)
                                      
                                    }
                                  }
                                }
                                })
                                
                                return(res)
                                
                              }

# Stop Distribution
stopCluster(cl)
close(pb)

# Save results to RData file
save(results.prediction, file="prediction.RData")

