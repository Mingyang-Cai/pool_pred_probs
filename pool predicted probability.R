# pool predicted probability MCAR
rm(list=ls())
library(mice, warn.conflicts = FALSE)
library(miceadds)
library(MASS)
library(magrittr)
library(dplyr)
library(furrr)
library(purrr)
library(tibble)
library(lme4)
library(readxl)
library(rstatix)
library(ggplot2)
library(caTools)
library(DPpack)
library(xgboost)
library(tidyr)
library(ggcorrplot)
library(ltm)
library(MuMIn)
library(performance)
library(specr)
library(ROCR)
library(caret)
library(car)
library(tree)
library(xgboost)
library(OptimalCutpoints)


# Create the correlation matrix for all variables
get_cormat <- function(cor, vars) {
  cormat <- diag(length(vars))
  cormat[cormat == 0] <- cor
  as.matrix(cormat)
}

# Make the regression coefficients given a prespecified effect size (specified in
# data simulation function)
get_coefs <- function(var_yhat, rel_weight, cormat) {
  sqrt(var_yhat / sum(rel_weight %*% t(rel_weight) * cormat)) * rel_weight
}

# Calculate an intercept based on a known prevalence (to vary class imbalance)
get_intercept <- function(betas, cormat, prevalence) {
  
  sigma <- c(sqrt(t(betas) %*% cormat %*% betas))
  
  integral <- function(x, b0) (1/sqrt(2 * pi)) * exp(-x^2 / 2) / (1 + exp(-b0 - sigma * x))
  
  f_integrate <- function(g, prevalence, g_lower, g_upper, ...) {
    prevalence - integrate(f = g, lower = g_lower, upper = g_upper, ...)$value
  }
  
  uniroot(f = f_integrate,
          interval = c(-99999999, 99999999),
          prevalence = prevalence,
          g = integral,
          g_lower = -Inf,
          g_upper = Inf)$root
}

#calculate predicted probability
pred <- function(formula, data){
  model <- glm(formula = formula, family="binomial", data = data)
  pred <- predict(model, data ,allow.new.levels =T) 
  pred <- exp(pred)/(exp(pred)+1)
  return(pred)
}

# Simulate the complete samples
sim_dat <- function(n, r2, prevalence = 0.5, rel_weight, cormat) {
  var_yhat <- (r2 * pi^2/3) / (1-r2)
  coefs <- get_coefs(var_yhat, rel_weight, cormat)
  b0 <- get_intercept(coefs, cormat, prevalence = prevalence)
  
  X <- MASS::mvrnorm(n, mu = rep(0, length(coefs)), Sigma = cormat)
  Y <- rbinom(n, 1, 1 / (1 + exp(-(b0 + X %*% coefs))))
  Y <- as.factor(Y)
  bind_cols(X = data.frame(X), Y = Y)
}


#predict probability bayesian
logreg.prob <- function(fit){
  fit.sum <- summary.glm(fit)
  beta <- fit.sum$coefficients[,"Estimate"]
  rv <- t(chol(fit.sum$cov.unscaled))
  beta.star <- beta + rv %*% rnorm(ncol(rv))
}

#predict probability frequentist
freq.predict <- function(glm, data){
  p <- predict(glm, newdata = data, type = "response")
  vec<-(runif(length(p)) <= p)
  vec[vec] <- 1 
  vec
}

#Generate incomplete data
set.seed(123)
mp <- 0.3  #missingness percentage
nsim <- 1
relative.strength <- rep(1, 4)
cormat <- get_cormat(0.3, relative.strength)
n    <- 1000
r2   <- c(0.5, 0.85)
prev <- c(0.2, 0.5)
mpattern <- matrix(c(1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0), nrow =3, byrow = TRUE)

plan(multisession)

sim <- future_map_dfr(1:nsim, 
                      ~expand_grid(N = n, R2 = r2, Prevalence = prev, 
                                   relative.strength = list(relative.strength),
                                   cormat = list(cormat)) %>%
                        mutate(dat = pmap(list(N, R2, Prevalence, relative.strength, cormat), 
                                          function(N, R2, Prevalence, relative.strength, cormat) {
                                            sim_dat(n = N, r2 = R2, prevalence = Prevalence,
                                                    rel_weight = relative.strength, cormat = cormat)
                                          }),
                               amp = map(dat, ~ampute(.x, prop = mp, mech = "MCAR", 
                                                      pattern = mpattern, 
                                                      freq = c(1/3, 1/3, 1/3))$amp)), # and specify other arguments here as well
                      .id = "NSIM",
                      .options = furrr_options(seed = TRUE),
                      .progress = TRUE)


evalution <- function(sim){
  output <- list()
  for (i in 1 : length(sim$dat)) {
    obs.data <- sim$dat[[i]][!is.na(sim$amp[[i]]$Y),]
    mis.data <- sim$dat[[i]][is.na(sim$amp[[i]]$Y),]
    glm1 <- glm(Y ~ X1 + X2 + X3 + X4, family = binomial, data = sim$dat[[i]])    
    true.1 <- replicate(n = 30, expr = freq.predict(glm1, mis.data) %>% as.factor %>%  
                confusionMatrix(sim$dat[[i]][is.na(sim$amp[[i]]$Y), "Y"]) %>% .$overall) %>% 
      rowMeans()
    glm2 <- glm(Y ~ X1 + X2 + X3 + X4, family = binomial, data = obs.data)
    true.2 <- replicate(n = 30, expr = freq.predict(glm2, mis.data) %>% as.factor %>%  
                confusionMatrix(sim$dat[[i]][is.na(sim$amp[[i]]$Y), "Y"]) %>% .$overall) %>% 
      rowMeans()
    cutoff.complete <- sim$dat[[i]] %>% mutate(pred.prob = pred(Y~ X1 + X2 + X3 + X4, .)) %>% 
      optimal.cutpoints(X = pred.prob  ~ Y, tag.healthy = 0, 
                        control = control.cutpoints(),
                        methods = "ROC01", data = ., ci.fit = FALSE, conf.level = 0.95, trace = FALSE) %>% 
      .$ROC01 %>% .$Global %>% .$optimal.cutoff %>% .$cutoff
    
    cutoff.obs<-  obs.data %>% mutate(pred.prob = pred(Y~ X1 + X2 + X3 + X4, .)) %>% 
      optimal.cutpoints(X = pred.prob  ~ Y, tag.healthy = 0, 
                        control = control.cutpoints(),
                        methods = "ROC01", data = ., ci.fit = FALSE, conf.level = 0.95, trace = FALSE) %>% 
      .$ROC01 %>% .$Global %>% .$optimal.cutoff %>% .$cutoff
    true.3 <- as.numeric(predict(glm1, newdata = mis.data, type = "response") > cutoff.complete) %>% as.factor() %>%  
      confusionMatrix(sim$dat[[i]][is.na(sim$amp[[i]]$Y), "Y"]) %>% 
      .$overall 
    true.4 <- as.numeric(predict(glm2, newdata = mis.data, type = "response") > cutoff.obs) %>% as.factor() %>%  
      confusionMatrix(sim$dat[[i]][is.na(sim$amp[[i]]$Y), "Y"]) %>% 
      .$overall   
    #imputation
    imp <- mice(sim$amp[[i]], m = 50, meth = c("", "", "", "pmm", "logreg"), print = FALSE)
    imp.data <-imp$imp$Y %>% as.data.frame 
    imp.data %<>% mutate_if(., is.factor, as.numeric)
    imp.data <- imp.data - 1
    pooled.prob.1 <- imp.data %>%apply(., 1, function(x) sum(x)/ ncol(imp$imp$Y))   
    #calculate cut-off points in each imputed dataset
    cut.off <- complete(imp, "all") %>% map(~.x %>% filter((is.na(sim$amp[[i]]$Y) == FALSE),) %>% 
                                                mutate(pred.prob = pred(Y~ X1 + X2 + X3 + X4, .))%>% 
                                                optimal.cutpoints(X = pred.prob  ~ Y, tag.healthy = 0, 
                                                                  control = control.cutpoints(),
                                                                  methods = "ROC01", data = ., ci.fit = FALSE, conf.level = 0.95, trace = FALSE)
    ) %>%
      map(~.x %>% .$ROC01 %>% .$Global %>% .$optimal.cutoff %>%
            .$cutoff) %>% unlist
    #impute with value 0 and 1, compare the pooled probability with some cut off values
    result.1 <- as.numeric(pooled.prob.1 > mean(cut.off)) %>% as.factor() %>%  
      confusionMatrix(sim$dat[[i]][is.na(sim$amp[[i]]$Y), "Y"])
    #impute with value 0 and 1, generate a bernoulli random variable with the pool probability
    result.2 <- replicate(30, expr = vapply(pooled.prob.1, function(x) rbinom(1, 1, x), as.integer(1L)) %>% 
                            as.factor() %>%  
                            confusionMatrix(sim$dat[[i]][is.na(sim$amp[[i]]$Y), "Y"]), simplify = FALSE)
    
    result.2 <- result.2 %>% map(~.x %>% .$overall) %>% unlist() %>% matrix(, ncol = 7, byrow = TRUE)  
    colnames(result.2) <- c("Accuracy", "Kappa", "AccuracyLower", "AccuracyUpper", 
                            "AccuracyNull", "AccuracyPValue", "McnemarPValue")
    #pool estimates of the model in each imputed dataset
    result.3 <- with(imp, glm(Y ~ X1 + X2 + X3 + X4, family = binomial)) %>% pool() %>% 
      .$pooled %>% .$estimate 
    #calculate the pooled probability with the aggregated model
    pooled.prob.2 <- sim$dat[[i]][is.na(sim$amp[[i]]$Y), c("X1","X2", "X3", "X4")] %>% bind_cols(Intercept = 1,.) %>% 
      as.matrix %*% as.matrix(result.3) %>% apply(., 2, function(x) 1 / (1 + exp(-x)))
    #compare the pooled probability with some cut off values
    result.3 <- as.numeric(pooled.prob.2 > mean(cut.off)) %>% as.factor() %>%  
      confusionMatrix(sim$dat[[i]][is.na(sim$amp[[i]]$Y), "Y"])
    #generate a bernoulli random variable with the pool probability
    result.4 <- replicate(30, expr = vapply(pooled.prob.2, function(x) rbinom(1, 1, x), as.integer(1L)) %>% 
                            as.factor() %>%  
                            confusionMatrix(sim$dat[[i]][is.na(sim$amp[[i]]$Y), "Y"]), simplify = FALSE)
    result.4 <- result.4 %>% map(~.x %>% .$overall) %>% unlist() %>% matrix(, ncol = 7, byrow = TRUE)  
    colnames(result.4) <- c("Accuracy", "Kappa", "AccuracyLower", "AccuracyUpper", 
                            "AccuracyNull", "AccuracyPValue", "McnemarPValue")
    #calculate predicted probability for incomplete cases in each imputed dataset
    para <- complete(imp, "all") %>% map(~.x %>% filter(!is.na(sim$amp[[i]]$Y)) %$% 
                                             glm(Y ~ X1 + X2 + X3 + X4, family = binomial) 
                                           %>% logreg.prob) 
    predict.data <- complete(imp, "all") %>% map(~.x %>% filter(is.na(sim$amp[[i]]$Y))%>% 
                                                     bind_cols(Intercept = 1,.) %>%
                                                     select(-Y)) 
    predict.probability <- matrix(0, ncol = length(para), nrow = sum(is.na(sim$amp[[i]]$Y)))
    for (j in 1 : length(para)) {
      predict.probability[, j] <- as.matrix(predict.data[[j]])%*%as.matrix(para[[j]]) %>% 
        apply(., 2, function(x) 1 / (1 + exp(-x)))
    }
    pooled.prob.3 <- apply(predict.probability, 1, mean)
    #predict probability in each imputed dataset, pool and compare with cut off 
    result.5 <- as.numeric(pooled.prob.3 > mean(cut.off)) %>% as.factor() %>%  
      confusionMatrix(sim$dat[[i]][is.na(sim$amp[[i]]$Y), "Y"])
    #predict probability in each imputed dataset, pool and generate random variable
    result.6 <- replicate(30, expr = vapply(pooled.prob.3, function(x) rbinom(1, 1, x), as.integer(1L)) %>% 
                            as.factor() %>%  
                            confusionMatrix(sim$dat[[i]][is.na(sim$amp[[i]]$Y), "Y"]), simplify = FALSE)
    result.6 <- result.6 %>% map(~.x %>% .$overall) %>% unlist() %>% matrix(, ncol = 7, byrow = TRUE)  
    colnames(result.6) <- c("Accuracy", "Kappa", "AccuracyLower", "AccuracyUpper", 
                            "AccuracyNull", "AccuracyPValue", "McnemarPValue")
    result.7 <- matrix(0, ncol = ncol(predict.probability), 
                       nrow = nrow(predict.probability))
    #compare probability with cut off in each imputed dataset and pool (majority vote)
    for (k in 1 : ncol(predict.probability)) {
      result.7[, k] <- as.numeric(predict.probability[, k] > cut.off[k])
    }
    
    result.7 <- result.7 %>% apply(., 1, function(x) sum(x)/ ncol(result.7)) %>% 
      is_greater_than(0.5) %>% as.numeric %>% as.factor() %>%  
      confusionMatrix(sim$dat[[i]][is.na(sim$amp[[i]]$Y), "Y"]) %>% .$overall
    output[[i]] <- list(true.1 = true.1[1], true.2 = true.2[1], 
         true.3 = true.3[1], true.4 = true.4[1], 
           result.1 = result.1$overall[1], 
           result.2 = result.2[,"Accuracy"] %>%mean(),
           result.3 = result.3$overall[1], 
           result.4 = result.4[,"Accuracy"] %>%mean(),
           result.5 = result.5$overall[1],
           result.6 = result.6[,"Accuracy"] %>%mean(),
           result.7 = result.7[1])
    
  }
  return(output)
}
result <- evalution(sim) %>% map (unlist)
sim$R2
sim$Prevalence
result


