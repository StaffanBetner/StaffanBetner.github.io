---
layout: post
title: Simultaneous confidence intervals for (derivatives of) GAM smooths
published: true
---

You need a *gam* object created with mgcv::gam(method = “REML”) in R. This procedure is the same that is implemented in Gavin Simpson’s package *gratia*, available from CRAN.

1. Extract corrected/unconditional (incorporating the smoothing parameter uncertainty) covariance matrix of the estimated parameters from the model, hereby denoted by $\hat{V}_c$. 

   ```{R}
   Vc <- vcov(gam_object, unconditional = TRUE) 
   # Removing the last argument might make it work for other regression objects as well, I have not tried
   ```

2. Extract the model matrix, $X_p$, from the estimated model with a dense evaluation grid, $X_i$ for the smooth(s) of interest. Other predictors can be set to a constant. 

   ```{r}
   Xp <- model.matrix(gam_object, newdata = dense_evaluation_grid)
   ```

3. Simulate many (e.g. 10 000) observations from $N(0,\hat{V}_c)$, denoted by $B_u$.

   ```{r}
   Bu <- MASS::mvrnorm(n = 10000, mu = rep(0, nrow(Vc)), Sigma = Vc)
   ```

      3.0 For derivatives:

      This is motivated by: $$\hat{f} (x) = \sum_{k=1}^{M} \hat{\beta}_k g_k(x)\rightarrow \hat{f’}(x) = \sum_{k=1}^{M} \hat{\beta}_k g’_k(x)$$

      3.1. Extract $X_p$ for $X_i + \epsilon$, where $\epsilon$ is close to zero (e.g. 0.00001). ($X^{(2)}_p$) Be sure to avoid to add $\epsilon$ to factor variables though. 

      ```{r}
      eps <- 0.00001
      Xp2 <- model.matrix(gam_object, newdata = dense_evaluation_grid + eps)
      ```

      3.2. Calculate approximation of first derivative by $ \frac{X^{(2)}_{p} - X_{p}}{\epsilon} $, and use as $X_p$ in consecutive steps. 
      
      ```{r}
      Xp <- (Xp2 - Xp) / eps
      ```

And then for each smooth of interest:

4. Set irrelevant columns in $X_p$ to 0, e.g. the intercept and other terms. If you want simultaneousness over several terms at the same time, you should put the uninteresting ones to 0. In some cases you want to include the intercept, which is typically a good idea if the estimated smooth is linear, or close to linear, or if you want the simultaneous confidence interval for $E[Y|X]$.

   ```{r}
   # Example with just one smooth, where we set the intercept in the model matrix to 0
   # The intercept is in the first column
   Xp[,1] <- 0
   ```

5. Calculate standard errors for the estimated smooth evaluated in $X_i$ by $diag(X_p\hat{V}_cX_p')^{1/2}$ or equivalently (and more efficient) in R: $rowSums(X_p\hat{V}_c\cdot X_p)^{1/2}$ and denote $S_e$.

   ```{r}
   Se <- sqrt(rowSums(Xp %*% (Vc*Xp)))
   ```

6. Calculate $abs(\frac{X_pB_u'}{S_x})$ (dimensions: rows($X_i$) times number of simulations).

   ```{r}
   absVal <- abs((Xp %*% t(Bu)) / Sx)
   ```

7. For each simulation (i.e. column) find maximum of those absolute values.

   ```{r}
   maximums <- apply(absVal, 2, max)
   ```

8. Use $1-\alpha$ quantile (type 8 in R) of the maximums as critical value for the simultaneous confidence intervals for that smooth ($m_{1-\alpha}$).

   ```{r}
   crit <- quantile(maximums, type = 8, prob = 0.95)
   ```

9. Predict the value of the estimated smooth, and calculate the simultaneous (with regard to the points in $X_i$) confidence interval for the smooth by $X_p\hat{\beta}\pm m_{1-\alpha}S_e=\hat{f(x)}\pm m_{1-\alpha}S_e$.

   ```{r}
   est <- (Xp %*% coef(gam_model))
   lower <- est - crit * Se
   upper <- est + crit * Se
   confidence_intervals <- data.frame(x = dense_evaluation_grid$x, 
                                      est = est, 
                                      lower = lower, 
                                      upper = upper)
   ```

10. Finally, you could plot it with ggplot2.

    ```{r}
    ggplot(data = confidence_intervals, 
           aes(x = x, 
           y = est, 
           ymin = lower, 
           ymax = upper)) + 
          geom_line() + 
          geom_ribbon(alpha = 0.3)
    ```
