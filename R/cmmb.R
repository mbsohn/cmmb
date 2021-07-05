#library(parallel); library(doParallel); library(ggplot2)
convert2comp <- function(x){
  return(x/sum(x));
}

get.init.lambda <- function(n.sample, n.feature, tol){
  f <- function(x, p) {(qnorm(1-x/p))^4 + 2*((qnorm(1-x/p))^2) - x;}
  k <- uniroot(f, lower=0, upper=n.feature-1, tol=tol, p=n.feature)$root
  lambda_0 <- sqrt(2/n.sample)*qnorm(1-k/n.feature)
  return(lambda_0)
}

comp_lasso_probit <- function(y, x, l_constraint, lam, mean_value, tol, max.iter){
  n <- nrow(x); p <- ncol(x)
  q <- 2*y-1; bet0 <- rep(0, p)
  gramC <- crossprod(l_constraint); diagC <- diag(gramC); dediagC <- gramC - diag(diagC)
  iter <- 0; ksi <- 0; ksi0 <- 1
  while (sum(abs(ksi-ksi0))>tol && iter<max.iter){
    ksi0 <- ksi; iter2 <- 0
    mu <- 1; bet <- rep(1,p)/p
    while (sum(abs(bet-bet0))>tol && iter2<1000){
      bet0 <- bet
      xbeta <- x %*% bet0
      phi <- dnorm(q*xbeta); PHI <- pnorm(q*xbeta)
      lmbd <- q*phi/PHI; LMBD <- diag(as.numeric(lmbd*(xbeta+lmbd)))
      u <- xbeta + solve(LMBD) %*% lmbd
      gramX <- crossprod(x, LMBD%*%x); diagX <- diag(gramX); dediagX <- gramX - diag(diagX)
      covXY <- crossprod(x, LMBD%*%u)
      term0 <- (covXY-dediagX%*%bet)/n - mu*(t(l_constraint)%*%ksi + dediagC%*%bet)
      term2 <- diagX/n + mu*diagC
      for(j in 1:p){
        term1 <- sign(term0[j])*max(0, abs(term0[j])-lam)
        bet[j] <- term1/term2[j]
        dif <- bet[j] - bet0[j]
        term0 <- term0 - dediagX[,j]*dif/n - dediagC[,j]*dif*mu
      }
      if(mean_value==TRUE){
        stp <- bet-bet0; t <- 1
        xbeta <- x %*% bet
        phi <- dnorm(q*xbeta); PHI <- pnorm(q*xbeta)
        lmbd <- q*phi/PHI
        ll_cur <- sum(log(PHI))
        eps <- crossprod(lmbd, x%*%stp); ll_new <- ll_cur + abs(eps)
        iter3 <- 1
        while ((ll_new+0.3*t*eps>ll_cur) & iter3<500){
          bet_new <- bet0 + t*stp
          xbeta_new <- x %*% bet_new
          PHI_new <- pnorm(q*xbeta_new)
          ll_new <- sum(log(PHI_new))
          t <- 0.9*t
          iter3 <- iter3 + 1
        }
        bet <- bet_new
      }
      iter2 <- iter2 + 1
    }
    dif2 <- l_constraint%*%bet
    ksi <- ksi + dif2
    term0 <- term0 - mu*t(l_constraint)%*%dif2
    iter <- iter + 1
  }
  return(list(beta=bet, u=u, LMBD=LMBD))
}

cmm_slr_probit <- function(y, X, contr, lam0, flag_mv, tol, max.iter){
  n <- nrow(X); p <- ncol(X); contr2 <- t(as.matrix(contr))
  if(is.null(lam0)) lam0 <- get.init.lambda(n, p, tol)
  sigma <- 1; sigma_2 <- 1; sigma_s <- 0.5; iter <- 1
  while (abs(sigma-sigma_s)>0.01 & iter<100){
    iter <- iter + 1
    sigma <- (sigma_s + sigma_2)/2
    lam <- sigma*lam0
    rslt <- comp_lasso_probit(y, X, contr2, lam, flag_mv, tol, max.iter)
    bet2 <- rslt$beta
    u <- rslt$u
    LMBD <- rslt$LMBD
    s <- sum(abs(bet2)>0.001); s <- min(s, n-1)
    sigma_s <- sqrt(t(u-X%*%bet2) %*% LMBD %*% (u-X%*%bet2))/sqrt(n-s-1)
    sigma_2 <- sigma
  }
  if(iter==100) print("Not converge!")
  sigma <- sigma_s
  bet <- bet2
  return(list(beta=bet, u=u, LMBD=LMBD, sigma=sigma, lambda0=lam0))
}

est.debias.B.probit <- function(y, M, t, x, lam0=NULL, flag_mv=FALSE, tol=1e-4, max.iter=5000){
  if(!is.vector(y)) y <- as.vector(y)
  if(!is.vector(t)) t <- as.vector(t)
  n <- length(y); k <- ncol(M)
  if(is.null(x)){
    Z <- cbind(1, log(M), t)
    n.tx <- 1
  } else{
    Z <- cbind(1, log(M), x, t)
    if(!is.matrix(x)) x <- as.matrix(x)
    n.tx <- 1 + ncol(x)
  }
  n.vrs <- 1 + k + n.tx
  contr <- c(0, rep(1/sqrt(k), k), rep(0, n.tx))
  Z.til <- Z %*% (diag(n.vrs)-tcrossprod(contr))
  est.param <- cmm_slr_probit(y, Z.til, contr, lam0, flag_mv, tol, max.iter)
  gam0 <- est.param$lambda0/1
  zbeta <- Z.til %*% est.param$beta; q <- 2*y-1
  Sig <- crossprod(Z.til, est.param$LMBD%*%Z.til)/n
  Sig2 <- Sig - diag(diag(Sig))
  Q <- diag(n.vrs) - tcrossprod(contr)
  M.til <- matrix(0, n.vrs, n.vrs)
  for(i in 1:n.vrs){
    gam <- gam0/2
    while(gam<0.5){
      gam <- gam*2
      mi <- rep(1,n.vrs)
      mi0 <- rep(0,n.vrs)
      iter <- 1
      while(sum(abs(mi-mi0))>tol & iter<max.iter){
        mi0 <- mi
        for(j in 1:n.vrs){
          v <- -Sig2[j,]%*%mi+Q[j,i]
          mi[j] <- sign(v)*max(0, abs(v)-gam)/Sig[j,j]
        }
        iter <- iter + 1
      }
      if(iter<max.iter) break
    }
    M.til[i,] <- mi
  }
  M.til <- Q%*%M.til
  debias.B <- est.param$beta + M.til%*%t(Z.til)%*%(est.param$LMBD%*%(est.param$u-Z.til%*%est.param$beta))/n
  cov.debias.B <- as.numeric(est.param$sigma^2)*M.til%*%Sig%*%t(M.til)/n
  return(list(debias.B=t(debias.B), cov.debias.B=cov.debias.B, bias.B=est.param$beta))
}

get_D <- function(k, u){
  D_of_u <- matrix(-u, nrow=(k-1), ncol=(k-1))
  diag(D_of_u) <- (k-1)*u
  return(D_of_u)
}

build_mat_Ds <- function(n.comp, n.covar){
  n.rc <- (n.comp-1)*(n.covar+2)
  mat_Ds <- matrix(0, nrow=n.rc, ncol=n.rc)
  return(mat_Ds)
}

est.comp.param <- function(M, t, x, w){
  n <- nrow(M); k <- ncol(M)
  if(is.null(x)){
    mat_Ds <- build_mat_Ds(k, 0)
    sum.wt.sq <- sum(w*t^2)
    sum.wt <- sum(w*t)
    mat_Ds[1:(k-1), 1:(k-1)] <- get_D(k, sum.wt.sq)
    mat_Ds[k:(2*(k-1)), 1:(k-1)] <- get_D(k, sum.wt)
    mat_Ds[k:(2*(k-1)), k:(2*(k-1))] <- get_D(k, sum(w))
    mf <- k*t(log(M)) %*% (w*t) - sum(log(M)*w*t)
    mg <- k*t(log(M)) %*% w - sum(log(M)*w)
    vec_ms <- c(mf[-k], mg[-k])
  } else{
    n.x <- ncol(x)
    mat_Ds <- build_mat_Ds(k, n.x)
    sum.wt.sq <- sum(w*t^2)
    sum.wt <- sum(w*t)
    mat_Ds[1:(k-1), 1:(k-1)] <- get_D(k, sum.wt.sq)
    mat_Ds[k:(2*(k-1)), 1:(k-1)] <- get_D(k, sum.wt)
    mat_Ds[k:(2*(k-1)), k:(2*(k-1))] <- get_D(k, sum(w))
    sum.wxt <- colSums(w*x*t)
    sum.wx <- colSums(w*x)
    for(i in 3:(n.x+2)){
      mat_Ds[((i-1)*(k-1)+1):(i*(k-1)), 1:(k-1)] <- get_D(k, sum.wxt[i-2])
      mat_Ds[((i-1)*(k-1)+1):(i*(k-1)), k:(2*(k-1))] <- get_D(k, sum.wx[i-2])
      sum.wxx <- colSums(w*x*x[,(i-2)])
      for(j in i:(n.x+2)){
        mat_Ds[((j-1)*(k-1)+1):(j*(k-1)), ((i-1)*(k-1)+1):(i*(k-1))] <- get_D(k, sum.wxx[j-2])
      }
    }
    mf <- k*t(log(M)) %*% (w*t) - sum(log(M)*w*t)
    mg <- k*t(log(M)) %*% w - sum(log(M)*w)
    m_psy <- k*t(log(M)) %*% (w*x) - tcrossprod(rep(1,k), colSums(t(log(M)*w) %*% x))
    vec_ms <- c(mf[-k], mg[-k], m_psy[-k,])
  }

  mat_Ds[upper.tri(mat_Ds)] <- t(mat_Ds)[upper.tri(mat_Ds)]
  ests.k_minus_1 <- matrix(exp(solve(mat_Ds, vec_ms)), nrow=(k-1))
  ests_k <- 1/(colSums(ests.k_minus_1)+1)
  ests.k_minus_1 <- ests.k_minus_1 * tcrossprod(rep(1,k-1), ests_k)
  comp.param <- rbind(ests.k_minus_1, ests_k)
  rownames(comp.param) <- 1:k

  if(is.null(x)){
    colnames(comp.param) <- c("a", "m0")
  } else{
    colnames(comp.param) <- c("a", "m0", paste("psi", 1:n.x, sep=""))
  }
  return(comp.param)
}

cmm_alt <- function(x){
  l.x <- length(x)
  return(log(x[-l.x]/x[l.x]))
}

est_var_U1 <- function(M, t, x, w, comp.param){
  alt_M <- t(apply(M, 1, cmm_alt))
  alt_param <- apply(comp.param, 2, cmm_alt)
  if(!is.matrix(alt_param)){
    alt_M <- t(alt_M)
    alt_param <- matrix(alt_param, nrow=1)
  }
  hat_alt_M <- tcrossprod(cbind(t,1,x), alt_param)
  alt_U1 <- alt_M - hat_alt_M
  return(cov(alt_U1))
}

est_cmm_probit <- function(Y, M, t, X, w=rep(1,length(Y))){
  rearranging_data <- function(treatment,x=X,n=length(Y)){
    tmp_dat <- cbind(rep(treatment,n),1,x)
    return(t(tmp_dat))
  }
  rslt_1 <- est.comp.param(M,t,X,w)
  rslt_2 <- est.debias.B.probit(Y, M, t, X)
  var_U1 <- est_var_U1(M,t,X,w,rslt_1)
  log_comp_param <- log(rslt_1)
  n.comp <- ncol(M)
  cln_b <- rslt_2$debias.B[2:(n.comp+1)]
  n.cln_param <- length(rslt_2$debias.B)
  if(is.null(X)){
    cln_param <- rslt_2$debias.B[c(n.cln_param,1)]
  } else{
    cln_param <- rslt_2$debias.B[c(n.cln_param,1,(ncol(M)+2):(n.cln_param-1))]
  }
  cde_val_1 <- t(cln_param) %*% rearranging_data(1) + t(cln_b) %*% log_comp_param %*% rearranging_data(0)
  cde_val_2 <- t(cln_param) %*% rearranging_data(0) + t(cln_b) %*% log_comp_param %*% rearranging_data(0)
  cide_val_1 <- t(cln_param) %*% rearranging_data(1) + t(cln_b) %*% log_comp_param %*% rearranging_data(1)
  cide_val_2 <- t(cln_param) %*% rearranging_data(1) + t(cln_b) %*% log_comp_param %*% rearranging_data(0)
  scaling_factor <- as.numeric(sqrt(t(cln_b[-n.comp]) %*% var_U1 %*% cln_b[-n.comp] + 1))
  ecde <- mean(pnorm(cde_val_1/scaling_factor) - pnorm(cde_val_2/scaling_factor))
  ecide <- mean(pnorm(cide_val_1/scaling_factor) - pnorm(cide_val_2/scaling_factor))
  cde_cide <- c(ecde, ecide)
  names(cde_cide) <- c("cde", "cide")
  lmparam <- as.numeric(rslt_2$debias.B)
  if(is.null(X)){
    names(lmparam) <- c("c0", paste0("b",1:n.comp), "c")
  } else{
    names(lmparam) <- c("c0", paste0("b",1:n.comp), paste0("g", 1:(n.cln_param-n.comp-2)), "c")
  }
  return(list(total=cde_cide, cwprod=log(n.comp*rslt_1[,1]) * cln_b,
              comparam=rslt_1, lmparam=lmparam, Sig1=var_U1))
}

cmmb <- function(Y, M, tr, X, n.cores=NULL, n.boot=2000, ci.method="empirical",
                  p.value=FALSE, ForSA=FALSE, max.rho=0.5, sig.level=0.05, FWER=FALSE,
                  w=rep(1,length(Y)), prec=1e-4, max.iter=1000){
  cmmb_rslt <- est_cmm_probit(Y,M,tr,X)
  cmmb_cde <- cmmb_rslt$total["cde"]
  cmmb_cide <- cmmb_rslt$total["cide"]
  if(is.null(n.cores)){
    n.cores <- detectCores(logical=TRUE)
  }
  registerDoParallel(n.cores)
  X.sa <- NULL
  cmm.nest <- foreach(j=1:n.boot) %dopar% {
    indx_1 <- which(Y>0); indx_0 <- which(Y==0)
    n.1 <- length(indx_1); n.0 <- length(indx_0)
    indx <- sort(c(sample(indx_1, n.1, replace=TRUE), sample(indx_0, n.0, replace=TRUE)))
    boot.Y <- Y[indx]; boot.M <- M[indx,]; boot.tr <- c(scale(tr[indx], scale=FALSE))
    if(is.null(X)){
      boot.X <- NULL
    } else{
      boot.X <- scale(X[indx,], scale=FALSE)
      X.sa <- cbind(X.sa, boot.X)
    }
    tmp_cmmb_est <- est_cmm_probit(boot.Y,boot.M,boot.tr,boot.X)
    tmp_cmmb_est$X.sa <- X.sa
    return(tmp_cmmb_est)
  }
  if(FWER==TRUE){
    n.test <- 2
  } else{
    n.test <- 1
  }
  boot_est_total_cmm <- unlist(lapply(cmm.nest, "[[", 1))
  boot_est_cde_cmm <- boot_est_total_cmm[names(boot_est_total_cmm)=="cde"]
  boot_est_cide_cmm <- boot_est_total_cmm[names(boot_est_total_cmm)=="cide"]
  cmmb_cwprod <- cmmb_rslt$cwprod
  boot_est_cwprod_cmm <- do.call(rbind, lapply(cmm.nest, "[[", 2))
  if(ForSA==FALSE){
    if(ci.method=="empirical"){
      del_cde_quant <- quantile(boot_est_cde_cmm-mean(boot_est_cde_cmm),
                                c(sig.level/(2*n.test), 1-sig.level/(2*n.test)))
      cde_lu <- cmmb_cde - rev(del_cde_quant)
      del_cide_quant <- quantile(boot_est_cide_cmm-mean(boot_est_cide_cmm),
                                 c(sig.level/(2*n.test), 1-sig.level/(2*n.test)))
      cide_lu <- cmmb_cide - rev(del_cide_quant)
      tmp_del_cwprod <- sweep(boot_est_cwprod_cmm, 2, colMeans(boot_est_cwprod_cmm), "-")
      del_cwprod_quant <- apply(tmp_del_cwprod, 2,
                                function(x) quantile(x, c(sig.level/(2*n.test), 1-sig.level/(2*n.test))))
      cwprod_lu <- cmmb_cwprod - t(apply(del_cwprod_quant,2,rev))
      if(p.value==TRUE){
        cde.null.dist <- boot_est_cde_cmm-mean(boot_est_cde_cmm)
        cide.null.dist <- boot_est_cide_cmm-mean(boot_est_cide_cmm)
        if(cmmb_cde<0){
          cde_pval <- min(1, mean(cde.null.dist<cmmb_cde)*2*n.test)
        } else{
          cde_pval <- min(1, mean(cde.null.dist>cmmb_cde)*2*n.test)
        }
        if(cmmb_cide<0){
          cide_pval <- min(1, mean(cide.null.dist<cmmb_cide)*2*n.test)
        } else{
          cide_pval <- min(1, mean(cide.null.dist>cmmb_cide)*2*n.test)
        }
        frslt <- list(total=rbind(c(cmmb_cde, cde_lu, cde_pval), c(cmmb_cide, cide_lu, cide_pval)),
                      cwprod=cbind(cmmb_cwprod, cwprod_lu))
        colnames(frslt$total) <- c("Estimate", "Lower Limit", "Upper Limit", "P value")
        rownames(frslt$total) <- c("DE","TIDE")
        colnames(frslt$cwprod) <- c("Estimate", "Lower Limit", "Upper Limit")
        rownames(frslt$cwprod) <- colnames(M)
      } else{
        frslt <- list(total=rbind(c(cmmb_cde, cde_lu), c(cmmb_cide, cide_lu)),
                      cwprod=cbind(cmmb_cwprod, cwprod_lu))
        colnames(frslt$total) <- c("Estimate", "Lower Limit", "Upper Limit")
        rownames(frslt$total) <- c("DE","TIDE")
        colnames(frslt$cwprod) <- c("Estimate", "Lower Limit", "Upper Limit")
        rownames(frslt$cwprod) <- colnames(M)
      }
    } else if(ci.method=="percentile"){
      del_cde_quant <- quantile(boot_est_cde_cmm, c(sig.level/(2*n.test), 1-sig.level/(2*n.test)))
      cde_lu <- del_cde_quant
      del_cide_quant <- quantile(boot_est_cide_cmm, c(sig.level/(2*n.test), 1-sig.level/(2*n.test)))
      cide_lu <- del_cide_quant
      del_cwprod_quant <- apply(boot_est_cwprod_cmm, 2,
                                function(x) quantile(x, c(sig.level/(2*n.test), 1-sig.level/(2*n.test))))
      cwprod_lu <- t(del_cwprod_quant)
      frslt <- list(total=rbind(c(cmmb_cde, cde_lu), c(cmmb_cide, cide_lu)),
                    cwprod=cbind(cmmb_cwprod, cwprod_lu))
      colnames(frslt$total) <- c("Estimate", "Lower Limit", "Upper Limit")
      rownames(frslt$total) <- c("DE","TIDE")
      colnames(frslt$cwprod) <- c("Estimate", "Lower Limit", "Upper Limit")
      rownames(frslt$cwprod) <- colnames(M)
    } else{
      stop("The ci.method argument must be either 'empirical' or 'percentile'!", call.=FALSE)
    }
  } else{
    del_cde_quant <- quantile(boot_est_cde_cmm-mean(boot_est_cde_cmm),
                              c(sig.level/(2*n.test), 1-sig.level/(2*n.test)))
    cde_lu <- cmmb_cde - rev(del_cde_quant)
    del_cide_quant <- quantile(boot_est_cide_cmm-mean(boot_est_cide_cmm),
                               c(sig.level/(2*n.test), 1-sig.level/(2*n.test)))
    cide_lu <- cmmb_cide - rev(del_cide_quant)
    tmp_del_cwprod <- sweep(boot_est_cwprod_cmm, 2, colMeans(boot_est_cwprod_cmm), "-")
    del_cwprod_quant <- apply(tmp_del_cwprod, 2,
                              function(x) quantile(x, c(sig.level/(2*n.test), 1-sig.level/(2*n.test))))
    cwprod_lu <- cmmb_cwprod - t(apply(del_cwprod_quant,2,rev))
    boot_est_comparam_cmm <- do.call(cbind, lapply(cmm.nest, "[[", 3))
    n.comparam <- length(unique(colnames(cmm.nest[[1]]$comparam)))
    boot_est_lmparam_cmm <- do.call(rbind, lapply(cmm.nest, "[[", 4))
    if(!is.null(X)){
      boot_Xs_cmm <- do.call(cbind, lapply(cmm.nest, "[[", 6))
      n.X <- ncol(X)
    } else{
      boot_Xs_cmm <- NULL
    }
    boot_est_altU1_cmm <- do.call(cbind, lapply(cmm.nest, "[[", 5))
    k <- nrow(cmm.nest[[1]]$comparam)

    if(abs(max.rho)>0.99) max.rho <- 0.99
    v.rho.est <- seq(0.01, max.rho, 0.01)
    n.rho.est <- length(v.rho.est)
    tmp.sa.b.rho.p <- tmp.sa.cide.rho.p <- NULL
    tmp.sa.b.rho.p <- as.matrix(t(cmmb_rslt$lmparam[2:k]),,drop=FALSE)
    tmp.sa.cide.rho.p <- cmmb_cide
    for(i in 1:n.rho.est){
      tmp.cide.rho.j <- cmmbsa(rho=v.rho.est[i],
                               comparam=cmmb_rslt$comparam,
                               lmparam=cmmb_rslt$lmparam,
                               Sig1=cmmb_rslt$Sig1,
                               b_rho=tmp.sa.b.rho.p[i,],
                               covar=X,
                               n=nrow(M), prec=prec, max.iter=max.iter)
      tmp.sa.cide.rho.p <- c(tmp.sa.cide.rho.p, tmp.cide.rho.j$cide)
      tmp.sa.b.rho.p <- rbind(tmp.sa.b.rho.p, tmp.cide.rho.j$b_rho)
    }
    tmp.sa.b.rho.n <- tmp.sa.cide.rho.n <- NULL
    tmp.sa.b.rho.n <- as.matrix(t(cmmb_rslt$lmparam[2:k]),,drop=FALSE)
    tmp.sa.cide.rho.n <- cmmb_cide
    v.rho.est <- rev(seq(-max.rho, -0.01, 0.01))
    n.rho.est <- length(v.rho.est)
    for(i in 1:n.rho.est){
      tmp.cide.rho.j <- cmmbsa(rho=v.rho.est[i],
                               comparam=cmmb_rslt$comparam,
                               lmparam=cmmb_rslt$lmparam,
                               Sig1=cmmb_rslt$Sig1,
                               b_rho=tmp.sa.b.rho.n[i,],
                               covar=X,
                               n=nrow(M), prec=prec, max.iter=max.iter)
      tmp.sa.cide.rho.n <- c(tmp.sa.cide.rho.n, tmp.cide.rho.j$cide)
      tmp.sa.b.rho.n <- rbind(tmp.sa.b.rho.n, tmp.cide.rho.j$b_rho)
    }
    sa.b.rho <- rbind(apply(tmp.sa.b.rho.n,2,rev), tmp.sa.b.rho.p[-1,,drop=FALSE])

    v.rho <- seq(-max.rho, max.rho, 0.01)
    n.rho <- length(v.rho)
    n_boots <- nrow(boot_est_lmparam_cmm)
    sa.rslt <- matrix(0, nrow=n.rho, ncol=4)
    colnames(sa.rslt) <- c("rho", "Estimate", "Lower Limit", "Upper Limit")
    sa.rslt[,"Estimate"] <- c(rev(tmp.sa.cide.rho.n), tmp.sa.cide.rho.p[-1])
    registerDoParallel(n.cores)
    for(i in 1:n.rho){
      tmp.sa.rslt <- NULL
      if(is.null(X)){
        tmp.sa.rslt <- foreach(j=1:n_boots, .combine=c) %dopar% {
          tmp.cide.rho.j <- cmmbsa(rho=v.rho[i],
                                   comparam=boot_est_comparam_cmm[,(1+n.comparam*(j-1)):(n.comparam*j)],
                                   lmparam=boot_est_lmparam_cmm[j,],
                                   Sig1=boot_est_altU1_cmm[,(1+(k-1)*(j-1)):((k-1)*j)],
                                   b_rho=sa.b.rho[i,],
                                   covar=NULL,
                                   n=nrow(M), prec=prec, max.iter=max.iter)
          tmp.cide.rho.j$cide
        }
      } else{
        tmp.sa.rslt <- foreach(j=1:n_boots, .combine=c) %dopar% {
          tmp.cide.rho.j <- cmmbsa(rho=v.rho[i],
                                   comparam=boot_est_comparam_cmm[,(1+n.comparam*(j-1)):(n.comparam*j)],
                                   lmparam=boot_est_lmparam_cmm[j,],
                                   Sig1=boot_est_altU1_cmm[,(1+(k-1)*(j-1)):((k-1)*j)],
                                   b_rho=sa.b.rho[i,],
                                   covar=boot_Xs_cmm[,(1+n.X*(j-1)):(n.X*j)],
                                   n=nrow(M), prec=prec, max.iter=max.iter)
          tmp.cide.rho.j$cide
        }
      }
      if(sum(is.na(tmp.sa.rslt))/n_boots > 0.1){
        cide_sa_lu <- c(NA, NA)
      } else{
        tmp_sa_del_cide_quant <- quantile(tmp.sa.rslt-mean(tmp.sa.rslt),
                                          c(sig.level/(2*n.test), 1-sig.level/(2*n.test)), na.rm=TRUE)
        cide_sa_lu <- sa.rslt[i,"Estimate"] - rev(tmp_sa_del_cide_quant)
      }
      tmp_sa_del_cide_quant <- quantile(tmp.sa.rslt-mean(tmp.sa.rslt),
                                        c(sig.level/(2*n.test), 1-sig.level/(2*n.test)), na.rm=TRUE)
      cide_sa_lu <- sa.rslt[i,"Estimate"] - rev(tmp_sa_del_cide_quant)
      sa.rslt[i, c("rho", "Lower Limit", "Upper Limit")] <- c(v.rho[i], cide_sa_lu)
    }
    sa.rslt <- sa.rslt[complete.cases(sa.rslt),]
    if(p.value==TRUE){
      cde.null.dist <- boot_est_cde_cmm-mean(boot_est_cde_cmm)
      cide.null.dist <- boot_est_cide_cmm-mean(boot_est_cide_cmm)
      if(cmmb_cde<0){
        cde_pval <- min(1, mean(cde.null.dist<cmmb_cde)*2*n.test)
      } else{
        cde_pval <- min(1, mean(cde.null.dist>cmmb_cde)*2*n.test)
      }
      if(cmmb_cide<0){
        cide_pval <- min(1, mean(cide.null.dist<cmmb_cide)*2*n.test)
      } else{
        cide_pval <- min(1, mean(cide.null.dist>cmmb_cide)*2*n.test)
      }
      frslt <- list(total=rbind(c(cmmb_cde, cde_lu, cde_pval), c(cmmb_cide, cide_lu, cide_pval)),
                    cwprod=cbind(cmmb_cwprod, cwprod_lu),
                    cide.rho=sa.rslt)
      colnames(frslt$total) <- c("Estimate", "Lower Limit", "Upper Limit", "P value")
      rownames(frslt$total) <- c("DE","TIDE")
      colnames(frslt$cwprod) <- c("Estimate", "Lower Limit", "Upper Limit")
      rownames(frslt$cwprod) <- colnames(M)
    } else{
      frslt <- list(total=rbind(c(cmmb_cde, cde_lu), c(cmmb_cide, cide_lu)),
                    cwprod=cbind(cmmb_cwprod, cwprod_lu),
                    cide.rho=sa.rslt)
      colnames(frslt$total) <- c("Estimate", "Lower Limit", "Upper Limit")
      rownames(frslt$total) <- c("DE","TIDE")
      colnames(frslt$cwprod) <- c("Estimate", "Lower Limit", "Upper Limit")
      rownames(frslt$cwprod) <- colnames(M)
    }
    l.cut <- u.cut <- NULL
    for(i in 1:10){
      if(abs(frslt$cide.rho[i+1,3])-abs(frslt$cide.rho[i,3]) > 0|
         abs(frslt$cide.rho[i+1,4])-abs(frslt$cide.rho[i,4]) > 0){
        l.cut <- c(l.cut, i)
      }
    }
    n.rho <- nrow(frslt$cide.rho)
    for(i in n.rho:(n.rho-9)){
      if(abs(frslt$cide.rho[i-1,3])-abs(frslt$cide.rho[i,3]) > 0|
         abs(frslt$cide.rho[i-1,4])-abs(frslt$cide.rho[i,4]) > 0){
        u.cut <- c(u.cut, i)
      }
    }
    if(length(u.cut)>0){
      frslt$cide.rho <- frslt$cide.rho[-c(min(u.cut):n.rho),]
    }
    if(length(l.cut)>0){
      frslt$cide.rho <- frslt$cide.rho[-c(1:max(l.cut)),]
    }
  }
  class(frslt) <- "cmmb"
  return(frslt)
}

cmmbsa <- function(rho, comparam, lmparam, Sig1, b_rho=NULL, covar=X, n=nrow(M), prec, max.iter){
  if(is.na(b_rho[1])) return(list(cide=NA, b_rho=NA))
  k <- nrow(comparam)
  alt_M <- log(sweep(comparam[-k,,drop=F], 2, comparam[k,], FUN="/"))
  cln_b <- lmparam[2:(k+1)]
  b_k <- cln_b[-k]
  b.alt.terms <- crossprod(alt_M, b_k)
  n.cln_param <- length(lmparam)
  if(is.null(covar)){
    cln_param <- lmparam[c(n.cln_param,1)]
  } else{
    cln_param <- lmparam[c(n.cln_param,1,(k+2):(n.cln_param-1))]
  }
  scaling_factor <- as.numeric(sqrt(crossprod(b_k, Sig1%*%b_k)+1))
  cln_param_mdl0 <- (cln_param + b.alt.terms)/scaling_factor
  cov01 <- (Sig1%*%b_k)/scaling_factor
  m <- length(b_k)
  if(is.null(b_rho)) b_rho <- b_k
  if(m==1){
    S <- sqrt(Sig1)
    A <- S^2*(cov01^2-S^2)
    B <- (cov01^2-S^2)*rho*S
    C <- (cov01^2-rho^2*S^2)
    sqrt.root.term <- round(B^2-A*C, 8)
    tmp.b_rho <- c((-B+sqrt(sqrt.root.term))/A, (-B-sqrt(sqrt.root.term))/A)
    del_tmp_b_rho <- abs(tmp.b_rho-b_rho)
    b_rho <- tmp.b_rho[which.min(del_tmp_b_rho)]
    scaling_factor_rho <- as.numeric(sqrt(crossprod(b_rho, Sig1%*%b_rho)+2*rho*crossprod(b_rho, S)+1))
    cide_rho_1 <- cbind(1, rep(1,n), covar)%*%cln_param_mdl0
    cide_rho_2 <- rep(crossprod(alt_M[,"a"], b_rho)/scaling_factor_rho, n)
    cide_rho <- mean(pnorm(cide_rho_1) - pnorm(cide_rho_1-cide_rho_2))
    return(list(cide=cide_rho, b_rho=b_rho))
  } else{
    S <- sqrt(diag(Sig1))
    del_b_rho <- 1; n.iter <- 1
    while(del_b_rho > prec & n.iter < max.iter){
      b_rho_old <- b_rho
      for(j in 1:m){
        A <- S[j]^2*(cov01[j]^2-S[j]^2)
        B <- (cov01[j]^2-S[j]^2)*(Sig1[j,-j]%*%b_rho[-j]+rho*S[j])
        Sig1_j <- Sig1[-j,-j,drop=FALSE]; diag(Sig1_j) <- 0
        c.term1 <- crossprod(b_rho[-j]^2, S[-j]^2) + crossprod(b_rho[-j], Sig1_j%*%b_rho[-j]) +
          2*rho*crossprod(b_rho[-j], S[-j]) + 1
        Sig1_j2 <- crossprod(Sig1[j,-j,drop=FALSE]); diag(Sig1_j2) <- 0
        c.term2 <- Sig1[j,-j]^2%*%b_rho[-j]^2 + rho^2*S[j]^2 + 2*rho*S[j]*Sig1[j,-j]%*%b_rho[-j] +
          crossprod(b_rho[-j], Sig1_j2%*%b_rho[-j])
        C <- cov01[j]^2*c.term1 - c.term2
        sqrt.root.term <- round(B^2-A*C, 8)
        if(sqrt.root.term < 0 | is.na(sqrt.root.term)){
          return(list(cide=NA, b_rho=NA))
        }
        tmp.b_rho <- c((-B+sqrt(sqrt.root.term))/A, (-B-sqrt(sqrt.root.term))/A)
        del_tmp_b_rho <- abs(tmp.b_rho-b_rho[j])
        b_rho[j] <- tmp.b_rho[which.min(del_tmp_b_rho)]
      }
      del_b_rho <- sum((b_rho-b_rho_old)^2)
      n.iter <- n.iter + 1
    }
    scaling_factor_rho <- as.numeric(sqrt(crossprod(b_rho, Sig1%*%b_rho)+2*rho*crossprod(b_rho, S)+1))
    cide_rho_1 <- cbind(1, rep(1,n), covar)%*%cln_param_mdl0
    cide_rho_2 <- rep(crossprod(alt_M[,"a"], b_rho)/scaling_factor_rho, n)
    cide_rho <- mean(pnorm(cide_rho_1) - pnorm(cide_rho_1-cide_rho_2))
    return(list(cide=cide_rho, b_rho=b_rho))
  }
}

plot_cmmb_sa <- function(cmmb.output){
  if(!is(cmmb.output, "cmmb")) stop("Input must be an output of cmmb()!", call.=FALSE)
  tide.sa <- cmmb.output$cide.rho
  if(is.null(tide.sa)) stop("Output does not contain sensitivity analysis results.\n Check if the argument, ForSA, was set to TRUE in cmmb()!", call.=FALSE)
  tide.sa <- as.data.frame(tide.sa)
  gg.sa <- ggplot(tide.sa, aes(rho, Estimate)) +
    geom_line() +
    geom_ribbon(aes(ymin=`Lower Limit`, ymax=`Upper Limit`), fill="grey", alpha=0.5) +
    labs(x=expression(rho), y="Expected Mediation Effect") +
    theme(axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10),
          legend.position="none")
  return(gg.sa)
}

plot_cw_ide <- function(cmmb.output){
  if(!is(cmmb.output, "cmmb")) stop("Input must be an output of cmmb()!", call.=FALSE)
  par(mar = c(7, 6, 3, 3) + 0.1, xpd=FALSE)
  cw.ide <- cmmb.output$cwprod
  n.comp <- nrow(cw.ide); err_w <- 0.05
  g <- plot(x=c(1:n.comp), y=cw.ide[,"Estimate"],
            ylim=c(min(cw.ide[,"Lower Limit"]), max(cw.ide[,"Upper Limit"])),
            xlab="", ylab="Product of Path Coefficients", xaxt="n",
            cex.axis=0.75, cex=0.75, cex.lab=0.9)
  for(i in 1:n.comp){
    segments(i, cw.ide[i,"Lower Limit"], i, cw.ide[i,"Upper Limit"])
    segments(i-err_w, cw.ide[i,"Lower Limit"], i+err_w, cw.ide[i,"Lower Limit"])
    segments(i-err_w, cw.ide[i,"Upper Limit"], i+err_w, cw.ide[i,"Upper Limit"])
  }
  abline(h=0, col="darkgray", lty=2)
  text(c(1:n.comp), par("usr")[3]-0.01, srt=-90, adj=0, labels=rownames(cw.ide), xpd=TRUE, cex=0.75)
}