library(parallel); library(brglm2)
gen.cmm.sim.data <- function(n.sample=50, n.comp=5, effect_mag=1, a=NULL, b=NULL){
 if(n.comp < 5) stop("Number of components must be >= 5!", call.=FALSE)
 n.non.sig <- n.comp - 4
 convert2comp <- function(x){
  return(x/sum(x))
 }
 if(is.null(a)){
  a <- c(20,10,5,2,rep(1,n.non.sig)); a <- convert2comp(a)
 }
 tr <- rbinom(n.sample, 1, 0.5)
 M0 <- rep(1, n.comp); M0 <- convert2comp(M0)
 X <- matrix(rnorm(n.sample,0,1), nrow=n.sample); colnames(X) <- "X"
 E1 <- matrix(rlnorm(n.sample*n.comp,0,1), nrow=n.sample)
 E1 <- t(apply(E1, 1, convert2comp))
 A2t <- t(outer(a,tr,FUN="^"))/rowSums(t(outer(a,tr,FUN="^")))
 tmp.M <- tcrossprod(rep(1,n.sample), M0) * A2t
 tmp.M <- tmp.M/rowSums(tmp.M)
 psy <- c(rep(1,length(a))); psy <- convert2comp(psy)
 Psy2X <- t(outer(psy,as.numeric(X),FUN="^"))/rowSums(t(outer(psy,as.numeric(X),FUN="^")))
 M <- (tmp.M*Psy2X*E1)/rowSums(tmp.M*Psy2X*E1)
 b <- c(0,0.5,-0.5,0.5,-0.5,rep(0,n.non.sig),rep(0,ncol(X)),1)*effect_mag
 colnames(M) <- paste0("taxon", 1:ncol(M))
 MX <- cbind(1, log(M), X, tr)
 Y <- rbinom(n.sample, 1, pnorm(MX%*%b))
 return(list(Y=Y, M=M, tr=tr, X=X))
}

est_pcs <- function(Y, M, tr, X){
 pc <- prcomp(log(M))
 pc.vars <- pc$sdev^2
 pc.vars <- pc.vars/sum(pc.vars)
 cum.prop <- cumsum(pc.vars)
 n.pc <- min(which(cum.prop>0.95)[1], round(length(Y)*0.1))
 res <- lm(pc$x[,1:n.pc] ~ tr + X)
 n.trx <- 2 + ncol(X)
 res.coeff <- lapply(summary(res), coefficients)
 est.pc.m.param <- matrix(0, nrow=n.pc, ncol=n.trx)
 pc.m.var <- NULL
 for(i in 1:n.pc){
  est.pc.m.param[i,] <- res.coeff[[i]][,1]
  pc.eps <- pc$x[,i] - cbind(1,tr,X)%*%est.pc.m.param[i,]
  pc.m.var <- c(pc.m.var, var(pc.eps))
 }
 Z.pc <- cbind(pc$x[,1:n.pc], X, tr)
 resy <- glm(Y ~ Z.pc, family = binomial(link = "probit"), method="brglmFit")
 pc.bet <- resy$coefficients[2:(n.pc+1)]
 n.pc.ln.param <- length(resy$coefficients)
 pc.ln.param <- resy$coefficients[c(1,n.pc.ln.param,(n.pc+2):(n.pc.ln.param-1))]
 ide.pc <- crossprod(est.pc.m.param[,2], pc.bet)
 de_ide <- c(resy$coefficients[n.pc+2], ide.pc)
 names(de_ide) <- c("cde", "cide")
 return(de_ide)
}

est_pcp <- function(Y, M, tr, X){
 pc.rearranging_data <- function(treatment,x=X,n=length(Y)){
  tmp_dat <- cbind(1,rep(treatment,n),x)
  return(t(tmp_dat))
 }
 pc <- prcomp(log(M))
 pc.vars <- pc$sdev^2
 pc.vars <- pc.vars/sum(pc.vars)
 cum.prop <- cumsum(pc.vars)
 n.pc <- min(which(cum.prop>0.95)[1], round(length(Y)*0.1))
 res <- lm(pc$x[,1:n.pc] ~ tr + X)
 n.trx <- 2 + ncol(X)
 res.coeff <- lapply(summary(res), coefficients)
 est.pc.m.param <- matrix(0, nrow=n.pc, ncol=n.trx)
 pc.m.var <- NULL
 for(i in 1:n.pc){
  est.pc.m.param[i,] <- res.coeff[[i]][,1]
  pc.eps <- pc$x[,i] - cbind(1,tr,X)%*%est.pc.m.param[i,]
  pc.m.var <- c(pc.m.var, var(pc.eps))
 }
 Z.pc <- cbind(pc$x[,1:n.pc], X, tr)
 resy <- glm(Y ~ Z.pc, family = binomial(link = "probit"), method="brglmFit")
 pc.bet <- resy$coefficients[2:(n.pc+1)]
 n.pc.ln.param <- length(resy$coefficients)
 pc.ln.param <- resy$coefficients[c(1,n.pc.ln.param,(n.pc+2):(n.pc.ln.param-1))]
 pc.cde_val_1 <- t(pc.ln.param) %*% pc.rearranging_data(1) + t(pc.bet) %*% est.pc.m.param %*% pc.rearranging_data(0)
 pc.cde_val_2 <- t(pc.ln.param) %*% pc.rearranging_data(0) + t(pc.bet) %*% est.pc.m.param %*% pc.rearranging_data(0)
 pc.cide_val_1 <- t(pc.ln.param) %*% pc.rearranging_data(1) + t(pc.bet) %*% est.pc.m.param %*% pc.rearranging_data(1)
 pc.cide_val_2 <- t(pc.ln.param) %*% pc.rearranging_data(1) + t(pc.bet) %*% est.pc.m.param %*% pc.rearranging_data(0)
 pc.scaling_factor <- as.numeric(sqrt(sum(pc.bet^2*pc.m.var) + 1))
 pc.ecde <- mean(pnorm(pc.cde_val_1/pc.scaling_factor) - pnorm(pc.cde_val_2/pc.scaling_factor))
 pc.ecide <- mean(pnorm(pc.cide_val_1/pc.scaling_factor) - pnorm(pc.cide_val_2/pc.scaling_factor))
 pc.cde_cide <- c(pc.ecde, pc.ecide)
 names(pc.cde_cide) <- c("cde", "cide")
 return(pc.cde_cide)
}

PCS <- function(Y, M, tr, X, n.boot=2000, sig.level=0.05, n.cores=NULL){
 pcs_boot <- function(...){
  cide_pcr_log <- cde_pcr_log <- NULL
  indx_1 <- which(Y>0); indx_0 <- which(Y==0)
  n.1 <- length(indx_1); n.0 <- length(indx_0)
  for(i in 1:n.boots){
   indx <- sort(c(sample(indx_1, n.1, replace=TRUE), sample(indx_0, n.0, replace=TRUE)))
   boot.Y <- Y[indx]; boot.M <- M[indx,]; boot.tr <- tr[indx]
   boot.X <- X[indx,,drop=FALSE]
   rslt_pcr_log <- est_pcs(boot.Y, boot.M, boot.tr, boot.X)
   cde_pcr_log <- c(cde_pcr_log, rslt_pcr_log[1])
   cide_pcr_log <- c(cide_pcr_log, rslt_pcr_log[2])
  }
  return(list(cde=cde_pcr_log, cide=cide_pcr_log))
 }
 pcs_rslt <- est_pcs(sim.dat$Y, sim.dat$M, sim.dat$tr, sim.dat$X)
 pcs_cde <- pcs_rslt["cde"]
 pcs_cide <- pcs_rslt["cide"]
 if(is.null(n.cores)){
  n.cores <- round(parallel::detectCores(logical=TRUE)/1.5)
 }
 n.boots <- ceiling(n.boot/n.cores)
 pcs.nest <- do.call(list, parallel::mclapply(seq_len(n.cores), pcs_boot, mc.cores=n.cores))
 boot_est_cde_pcs <- unlist(lapply(pcs.nest, "[[", 1))
 del_cde_quant <- quantile(boot_est_cde_pcs-mean(boot_est_cde_pcs),
                           c(sig.level/2, 1-sig.level/2))
 cde_lu <- pcs_cde - rev(del_cde_quant)
 boot_est_cide_pcs <- unlist(lapply(pcs.nest, "[[", 2))
 del_cide_quant <- quantile(boot_est_cide_pcs-mean(boot_est_cide_pcs),
                            c(sig.level/2, 1-sig.level/2))
 cide_lu <- pcs_cide - rev(del_cide_quant)
 frslt <- rbind(c(pcs_cde, cde_lu), c(pcs_cide, cide_lu))
 colnames(frslt) <- c("Estimate", "Lower Limit", "Upper Limit")
 rownames(frslt) <- c("DE","TIDE")
 return(frslt)
}

PCP <- function(Y, M, tr, X, n.boot=2000, sig.level=0.05, n.cores=NULL){
 pcp_boot <- function(...){
  cide_pcr_log <- cde_pcr_log <- NULL
  indx_1 <- which(Y>0); indx_0 <- which(Y==0)
  n.1 <- length(indx_1); n.0 <- length(indx_0)
  for(i in 1:n.boots){
   indx <- sort(c(sample(indx_1, n.1, replace=TRUE), sample(indx_0, n.0, replace=TRUE)))
   boot.Y <- Y[indx]; boot.M <- M[indx,]; boot.tr <- tr[indx]
   boot.X <- X[indx,,drop=FALSE]
   rslt_pcr_log <- est_pcp(boot.Y, boot.M, boot.tr, boot.X)
   cde_pcr_log <- c(cde_pcr_log, rslt_pcr_log[1])
   cide_pcr_log <- c(cide_pcr_log, rslt_pcr_log[2])
  }
  return(list(cde=cde_pcr_log, cide=cide_pcr_log))
 }
 pcp_rslt <- est_pcp(sim.dat$Y, sim.dat$M, sim.dat$tr, sim.dat$X)
 pcp_cde <- pcp_rslt["cde"]
 pcp_cide <- pcp_rslt["cide"]
 if(is.null(n.cores)){
  n.cores <- round(parallel::detectCores(logical=TRUE)/1.5)
 }
 n.boots <- ceiling(n.boot/n.cores)
 pcp.nest <- do.call(list, parallel::mclapply(seq_len(n.cores), pcp_boot, mc.cores=n.cores))
 boot_est_cde_pcp <- unlist(lapply(pcp.nest, "[[", 1))
 del_cde_quant <- quantile(boot_est_cde_pcp-mean(boot_est_cde_pcp),
                           c(sig.level/2, 1-sig.level/2))
 cde_lu <- pcp_cde - rev(del_cde_quant)
 boot_est_cide_pcp <- unlist(lapply(pcp.nest, "[[", 2))
 del_cide_quant <- quantile(boot_est_cide_pcp-mean(boot_est_cide_pcp),
                            c(sig.level/2, 1-sig.level/2))
 cide_lu <- pcp_cide - rev(del_cide_quant)
 frslt <- rbind(c(pcp_cde, cde_lu), c(pcp_cide, cide_lu))
 colnames(frslt) <- c("Estimate", "Lower Limit", "Upper Limit")
 rownames(frslt) <- c("DE","TIDE")
 return(frslt)
}