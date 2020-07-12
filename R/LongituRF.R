#' (S)MERF algorithm
#'
#' @param X [matrix]:
#' @param Y [vector]:
#' @param id [vector]:
#' @param Z [matrix]:
#' @param iter [numeric]:
#' @param mtry [numeric]:
#' @param ntree [numeric]:
#' @param time [vector]:
#' @param sto [character]:
#' @param delta [numeric]:
#'
#' @import randomForest
#' @import stats
#' @return
#' @export
#'
MERF <- function(X,Y,id,Z,iter,mtry,ntree, time, sto, delta = 0.001){
  q <- dim(Z)[2]
  nind <- length(unique(id))
  btilde <- matrix(0,nind,q) #### Pour la ligne i, on a les effets al?atoires de l'individu i
  sigmahat <- 1 #### init
  Btilde <- diag(rep(1,q)) ### init
  epsilonhat <- rep(0,length(Y))
  id_btilde <- unique(id)
  Tiime <- sort(unique(time))
  omega <- rep(0,length(Y))
  sigma2 <- 1
  Vrai <- NULL
  inc <- 1
  OOB <- NULL

  if (class(sto)=="character"){
    if (sto=="fbm"){
      id_omega <- matrix(0,nind,length(unique(time)))
      for (i in 1:length(unique(id))){
        w <- which(id ==id_btilde[i])
        time11 <- time[w]
        where <- NULL
        for (j in 1:length(time11)){
          where <- c(where,which(Tiime==time11[j]))
        }
        id_omega[i,where] <- 1
      }
      omega <- matrix(0,nind,length(unique(time)))
      omega2 <- rep(0,length(Y))
      h <- opti.FBM(X,Y,id,Z,iter, mtry,ntree,time)
      for (i in 1:iter){
        ystar <- rep(0,length(Y))
        for (k in 1:nind){ #### on retrace les effets al?atoires
          indiv <- which(id==unique(id)[k])
          ystar[indiv] <- Y[indiv]- Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
        }

        forest <- randomForest(X,ystar,mtry=mtry,ntree=ntree, importance = TRUE) ### on construit l'arbre
        fhat <- predict(forest) #### pr?diction avec l'arbre
        OOB[i] <- forest$mse[ntree]
        for (k in 1:nind){ ### calcul des effets al?atoires par individu
          indiv <- which(id==unique(id)[k])
          K <- cov.fbm(time[indiv], h)
          V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,, drop=FALSE])+diag(as.numeric(sigmahat),length(indiv),length(indiv))+ sigma2*K
          btilde[k,] <- Btilde%*%t(Z[indiv,, drop=FALSE])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
          omega[k,which(id_omega[k,]==1)] <- sigma2*K%*%solve(V)%*%(Y[indiv]-fhat[indiv])
          omega2[indiv] <- omega[k,which(id_omega[k,]==1)]
          epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[k,which(id_omega[k,]==1)]
        }
        sigm <- sigmahat
        B <- Btilde
        sigmahat <- sig.fbm(Y,sigmahat,id, Z, epsilonhat, Btilde, time, sigma2,h) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
        Btilde  <- bay.fbm(btilde,Btilde,Z,id,sigm, time, sigma2,h) #### MAJ des param?tres de la variance des effets al?atoires.
        ### MAJ de la volatilit? du processus stochastique
        sigma2 <- gam_fbm(Y,sigm,id,Z,B,time,sigma2,omega,id_omega,h)
        Vrai <- c(Vrai, logV.fbm(Y,fhat,Z[,,drop=FALSE],time,id,Btilde,sigma2,sigmahat,h))
        if (i>1) inc <- (Vrai[i-1]-Vrai[i])/Vrai[i-1]
        if (inc < delta) {
          print(cat("stopped after", i, "iterations."))
          return(list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat,sigma_sto=sigma2, id_btilde=unique(id), sto= sto, vraisemblance = Vrai,id=id, time =time, Hurst=h, OOB =OOB, omega=omega2))
        }
      }
      return(list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id),sigma_sto=sigma2,omega=omega2, sigma_sto =sigma2, time = time, sto= sto, Hurst =h, id=id, Vraisemblance=Vrai, OOB =OOB))
    }


    if (sto=="exp"){
      id_omega <- matrix(0,nind,length(unique(time)))
      for (i in 1:length(unique(id))){
        w <- which(id ==id_btilde[i])
        time11 <- time[w]
        where <- NULL
        for (j in 1:length(time11)){
          where <- c(where,which(Tiime==time11[j]))
        }
        id_omega[i,where] <- 1
      }
      omega <- matrix(0,nind,length(unique(time)))
      omega2 <- rep(0,length(Y))
      alpha <- opti.exp(X,Y,id,Z,iter, mtry,ntree,time)
      for (i in 1:iter){
        ystar <- rep(0,length(Y))
        for (k in 1:nind){ #### on retrace les effets al?atoires
          indiv <- which(id==unique(id)[k])
          ystar[indiv] <- Y[indiv]- Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
        }

        forest <- randomForest(X,ystar,mtry=mtry,ntree=ntree, importance = TRUE) ### on construit l'arbre
        fhat <- predict(forest) #### pr?diction avec l'arbre
        OOB[i] <- forest$mse[ntree]
        for (k in 1:nind){ ### calcul des effets al?atoires par individu
          indiv <- which(id==unique(id)[k])
          K <- cov.exp(time[indiv], alpha)
          V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,, drop=FALSE])+diag(as.numeric(sigmahat),length(indiv),length(indiv))+ sigma2*K
          btilde[k,] <- Btilde%*%t(Z[indiv,, drop=FALSE])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
          omega[k,which(id_omega[k,]==1)] <- sigma2*K%*%solve(V)%*%(Y[indiv]-fhat[indiv])
          omega2[indiv] <- omega[k,which(id_omega[k,]==1)]
          epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
        }
        sigm <- sigmahat
        B <- Btilde
        sigmahat <- sig.exp(Y,sigmahat,id, Z, epsilonhat, Btilde, time, sigma2,alpha) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
        Btilde  <- bay.exp(btilde,Btilde,Z,id,sigm, time, sigma2,alpha) #### MAJ des param?tres de la variance des effets al?atoires.
        ### MAJ de la volatilit? du processus stochastique
        sigma2 <- gam_exp(Y,sigm,id,Z,B,time,sigma2,omega,id_omega,alpha)
        Vrai <- c(Vrai,logV.exp(Y,fhat,Z[,,drop=FALSE],time,id,Btilde,sigma2,sigmahat,alpha))
        if (i>1) inc <- (Vrai[i-1]-Vrai[i])/Vrai[i-1]
        if (inc < delta) {
          print(cat("stopped after", i, "iterations."))
          return(list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, vraisemblance = Vrai,id=id, time=time, alpha = alpha, OOB =OOB, omega=omega2))
        }
      }
      return(list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), omega=omega2, sigma_sto =sigma2, time = time, sto= sto, alpha=alpha, id=id, Vraisemblance=Vrai, OOB =OOB))
    }

    if ( sto=="none"){
      for (i in 1:iter){
        ystar <- rep(NA,length(Y))
        for (k in 1:nind){ #### on retrace les effets al?atoires
          indiv <- which(id==unique(id)[k])
          ystar[indiv] <- Y[indiv]- Z[indiv,,drop=FALSE]%*%btilde[k,]
        }

        forest <- randomForest(X,ystar,mtry=mtry,ntree=ntree, importance = TRUE) ### on construit l'arbre
        fhat <- predict(forest)
        OOB[i] <- forest$mse[ntree]
        for (k in 1:nind){
          indiv <- which(id==unique(id)[k])
          V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,, drop=FALSE])+diag(as.numeric(sigmahat),length(indiv),length(indiv))
          btilde[k,] <- Btilde%*%t(Z[indiv,, drop=FALSE])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
          epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,, drop=FALSE]%*%btilde[k,]
        }

        sigm <- sigmahat
        sigmahat <- sig(sigma = sigmahat,id = id, Z = Z, epsilon = epsilonhat,Btilde = Btilde)
        Btilde  <- bay(bhat = btilde,Bhat = Btilde,Z = Z,id = id,sigmahat = sigm)
        Vrai <- c(Vrai, logV(Y,fhat,Z,time,id,Btilde,0,sigmahat,sto))
        if (i>1) inc <-abs((Vrai[i-1]-Vrai[i])/Vrai[i-1])
        if (inc < delta) {
          print(cat("stopped after", i, "iterations."))
          return(list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, vraisemblance = Vrai,id=id, time=time, OOB =OOB))
        }
      }
      return(list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, vraisemblance=Vrai,id=id, time=time, OOB =OOB))
    }
  }
  for (i in 1:iter){
    ystar <- rep(0,length(Y))
    for (k in 1:nind){
      indiv <- which(id==unique(id)[k])
      ystar[indiv] <- Y[indiv]- Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
    }

    forest <- randomForest(X,ystar,mtry=mtry,ntree=ntree, importance=TRUE)
    fhat <- predict(forest)
    OOB[i] <- forest$mse[ntree]
    for (k in 1:nind){
      indiv <- which(id==unique(id)[k])
      K <- sto_analysis(sto,time[indiv])
      V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,, drop=FALSE])+diag(as.numeric(sigmahat),length(indiv),length(indiv))+ sigma2*K
      btilde[k,] <- Btilde%*%t(Z[indiv,, drop=FALSE])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
      omega[indiv] <- sigma2*K%*%solve(V)%*%(Y[indiv]-fhat[indiv])
      epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
    }
    sigm <- sigmahat
    B <- Btilde
    sigmahat <- sig_sto(sigmahat,id, Z, epsilonhat, Btilde, time, sigma2,sto)
    Btilde  <- bay_sto(btilde,Btilde,Z,id,sigm, time, sigma2,sto)
    sigma2 <- gam_sto(sigm,id,Z,B,time,sigma2,sto,omega)
    Vrai <- c(Vrai, logV(Y,fhat,Z[,,drop=FALSE],time,id,Btilde,sigma2,sigmahat,sto))
    if (i>1) inc <- abs((Vrai[i-1]-Vrai[i])/Vrai[i-1])
    if (inc < delta) {
      print(cat("stopped after", i, "iterations."))
      return(list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), omega=omega, sigma_sto =sigma2, time = time, sto= sto,Vraisemblance=Vrai,id=id, OOB =OOB))
    }
  }
  return(list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id),omega=omega, sigma_sto =sigma2, time = time, sto= sto,Vraisemblance=Vrai,id=id, OOB =OOB))
}

#' Title
#'
#'
#' @import stats
#'
#' @keywords internal
bay_sto <- function(bhat,Bhat,Z,id, sigmahat, time, sigma2,sto){ #### actualisation des param?tres de B
  nind <- length(unique(id))
  q <- dim(Z)[2]
  Nombre <- length(id)
  D <- 0
  for (j in 1:nind){
    w <- which(id==unique(id)[j])
    K <- sto_analysis(sto,time[w])
    V <- Z[w,, drop=FALSE]%*%Bhat%*%t(Z[w,,drop=FALSE])+diag(as.numeric(sigmahat),length(w),length(w))+sigma2*K
    D <- D+ (bhat[j,]%*%t(bhat[j,]))+ (Bhat- Bhat%*%t(Z[w,, drop=FALSE])%*%solve(V)%*%Z[w,,drop=FALSE]%*%Bhat)
  }
  D <- D/nind
  return(D)
}

#' Title
#'
#' @import stats
#'
#' @keywords internal
sig_sto <- function(sigma,id,Z, epsilon, Btilde, time, sigma2,sto){ #### fonction d'actualisation du param?tre de la variance des erreurs
  nind <- length(unique(id))
  Nombre <- length(id)
  sigm <- 0
  for (j in 1:nind){
    w <- which(id==unique(id)[j])
    K <- sto_analysis(sto,time[w])
    V <- Z[w,, drop=FALSE]%*%Btilde%*%t(Z[w,, drop=FALSE])+diag(as.numeric(sigma),length(w),length(w))+sigma2*K
    sigm <- sigm + t(epsilon[w])%*%epsilon[w] + sigma*(length(w)-sigma*(sum(diag(solve(V)))))
  }
  sigm <- sigm/Nombre
  return(sigm)
}

#' Title
#'
#' @import stats
#'
#' @keywords internal
gam_sto <- function(sigma,id,Z, Btilde, time, sigma2,sto, omega){
  nind <- length(unique(id))
  Nombre <- length(id)
  gam <- 0
  for (k in 1:nind){
    indiv <- which(id==unique(id)[k])
    K <- sto_analysis(sto,time[indiv])
    V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,,drop=FALSE])+diag(as.numeric(sigma),length(indiv),length(indiv))+ sigma2*K
    Omeg <- omega[indiv]
    gam <-gam+ (t(Omeg)%*%solve(K)%*%Omeg) + sigma2*(length(indiv)-sigma2*sum(diag(solve(V)%*%K)))
  }
  return(as.numeric(gam)/Nombre)
}


#' Predict with a (S)MERF alorithm
#'
#' @param merf_sto :
#' @param X :
#' @param Z :
#' @param id :
#' @param time :
#'
#' @import stats
#' @import randomForest
#'
#' @return
#' @export
predict_merf <- function(merf_sto, X,Z,id,time){
  n <- length(unique(id))
  id_btilde <- merf_sto$id_btilde
  f <- predict(merf_sto$forest,X)
  Time <- merf_sto$time
  id_btilde <- merf_sto$id_btilde
  Ypred <- rep(0,length(id))
  id.app=merf_sto$id
  if (merf_sto$sto=="none"){
    for (i in 1:length(unique(id))){
      w <- which(id==unique(id)[i])
      k <- which(id_btilde==unique(id)[i])
      Ypred[w] <- f[w] + Z[w,, drop=FALSE]%*%merf_sto$random_effects[k,]
    }
    return(Ypred)
  }

  if (merf_sto$sto=="exp"){
    for (i in 1:length(unique(id))){
      w <- which(id==unique(id)[i])
      k <- which(id_btilde==unique(id)[i])
      om <- which(id.app==unique(id)[i])
      Ypred[w] <- f[w] + Z[w,, drop=FALSE]%*%merf_sto$random_effects[k,] + predict.exp(merf_sto$omega[om],Time[om],time[w], merf_sto$alpha)
    }
    return(Ypred)
  }

  if (merf_sto$sto=="fbm"){
    for (i in 1:length(unique(id))){
      w <- which(id==unique(id)[i])
      k <- which(id_btilde==unique(id)[i])
      om <- which(id.app==unique(id)[i])
      Ypred[w] <- f[w] + Z[w,, drop=FALSE]%*%merf_sto$random_effects[k,] + predict.fbm(merf_sto$omega[om],Time[om],time[w], merf_sto$Hurst)
    }
    return(Ypred)
  }

  for (i in 1:length(unique(id))){
    w <- which(id==unique(id)[i])
    k <- which(id_btilde==unique(id)[i])
    om <- which(id.app==unique(id)[i])
    Ypred[w] <- f[w] + Z[w,, drop=FALSE]%*%merf_sto$random_effects[k,] + predict.sto(merf_sto$omega[om],Time[om],time[w], merf_sto$sto)
  }
  return(Ypred)
}



#' blabla
#'
#' @import stats
#'
#' @keywords internal
predict.exp <- function(omega,time.app,time.test, alpha){
  pred <- rep(0,length(time.test))
  for (i in 1:length(time.test)){
    inf <- which(time.app<=time.test[i])
    sup <- which(time.app>time.test[i])
    if (length(inf)>0){
      if (length(sup)>0){
        time_inf <- max(time.app[inf])
        time_sup <- min(time.app[sup])
        pred[i] <- mean(c(omega[which(time.app==time_inf)],omega[which(time.app==time_sup)]))}
      if(length(sup)==0) {time_inf <- max(time.app[inf])
      pred[i] <- omega[which(time.app==time_inf)]*exp(-alpha*abs(time.test[i]-max(time.app)))
      }
    }
    if (length(sup)>0 & length(inf)==0){
      time_sup <- min(time.app[sup])
      pred[i] <- omega[which(time.app==time_sup)]*exp(-alpha*abs(time.test[i]-max(time.app)))
    }
    return(pred)
  }
}

#' Title
#'
#' @import stats
#'
#' @keywords internal
predict.fbm <- function(omega,time.app,time.test, h){
  pred <- rep(0,length(time.test))
  for (i in 1:length(time.test)){
    inf <- which(time.app<=time.test[i])
    sup <- which(time.app>time.test[i])
    if (length(inf)>0){
      if (length(sup)>0){
        time_inf <- max(time.app[inf])
        time_sup <- min(time.app[sup])
        pred[i] <- mean(c(omega[which(time.app==time_inf)],omega[which(time.app==time_sup)]))}
      if(length(sup)==0) {time_inf <- max(time.app[inf])
      pred[i] <- omega[which(time.app==time_inf)]*0.5*(time.test[i]^(2*h)+max(time.app)^(2*h)-abs(time.test[i]-max(time.app))^(2*h))/(max(time.app)^(2*h))
      }
    }
    if (length(sup)>0 & length(inf)==0){
      time_sup <- min(time.app[sup])
      pred[i] <- omega[which(time.app==time_sup)]*0.5*(time.test[i]^(2*h)+min(time.app)^(2*h)-abs(time.test[i]-min(time.app))^(2*h))/(min(time.app)^(2*h))
    }
    return(pred)
  }
}

#' Title
#'
#' @import stats
#'
#' @keywords internal
predict.sto <- function(omega,time.app,time.test, sto){
  pred <- rep(0,length(time.test))

  if (class(sto)=="function"){
    for (i in 1:length(time.test)){
      inf <- which(time.app<=time.test[i])
      sup <- which(time.app>time.test[i])
      if (length(inf)>0){
        if (length(sup)>0){
          time_inf <- max(time.app[inf])
          time_sup <- min(time.app[sup])
          pred[i] <- mean(c(omega[which(time.app==time_inf)],omega[which(time.app==time_sup)]))}
        if(length(sup)==0) {time_inf <- max(time.app[inf])
        pred[i] <- omega[which(time.app==time_inf)]*((sto(time.test[i],max(time.app)))/sto(max(time.app),max(time.app)))
        }
      }
      if (length(sup)>0 & length(inf)==0){
        time_sup <- min(time.app[sup])
        pred[i] <- omega[which(time.app==time_sup)]*((sto(time.test[i],min(time.app)))/sto(min(time.app),min(time.app)))
      }
    }
    return(pred)
  }

  else {
    for (i in 1:length(time.test)){
      inf <- which(time.app<=time.test[i])
      sup <- which(time.app>time.test[i])
      if (length(inf)>0){
        if (length(sup)>0){
          time_inf <- max(time.app[inf])
          time_sup <- min(time.app[sup])
          pred[i] <- mean(c(omega[which(time.app==time_inf)],omega[which(time.app==time_sup)]))}
        if(length(sup)==0) {time_inf <- max(time.app[inf])
        if (sto=="BM"){
          pred[i] <- omega[which(time.app==time_inf)]}
        if (sto=="OrnUhl"){
          pred[i] <- omega[which(time.app==time_inf)]*(exp(-abs(time.test[i]-max(time.app))/2))
        }
        if (sto=="BBridge"){
          pred[i] <- omega[which(time.app==time_inf)]*((1-time.test[i])/(1-max(time.app)^2))
        }
        }
      }
      if (length(sup)>0 & length(inf)==0){
        time_sup <- min(time.app[sup])
        if (sto=="BM"){
          pred[i] <- omega[which(time.app==time_sup)]*(time.test[i]/min(time.app))}
        if (sto=="OrnUhl"){
          pred[i] <- omega[which(time.app==time_sup)]*(exp(-abs(time.test[i]-min(time.app))/2))
        }
        if (sto=="BBridge"){
          pred[i] <- omega[which(time.app==time_sup)]*(time.test[i]/min(time.app))
        }
      }
    }}
  return(pred)
}

#' Title
#'
#' @import stats
#'
#' @keywords internal
sto_analysis <- function(sto, time){
  MAT <- matrix(0,length(time), length(time))

  if (class(sto)=="function"){
    for (i in 1:length(time)){
      for (j in 1:length(time)){
        MAT[i,j] <- sto(time[i], time[j])
      }
    }
    return(MAT)
  }

  if (sto=="BM"){
    for (i in 1:length(time)){
      for (j in 1:length(time)){
        MAT[i,j] <- min(time[i], time[j])
      }
    }
    return(MAT)
  }

  if (sto=="OrnUhl"){
    for (i in 1:length(time)){
      for (j in 1:length(time)){
        MAT[i,j] <- exp(-abs(time[i]-time[j])/2)
      }
    }
    return(MAT)
  }

  if (sto=="BBridge"){
    for (i in 1:length(time)){
      for (j in 1:length(time)){
        MAT[i,j] <- min(time[i], time[j]) - time[i]*time[j]
      }
    }
    return(MAT)
  }

}


#' Title
#'
#'
#' @import stats
#'
#' @keywords internal
sig <- function(sigma,id,Z, epsilon, Btilde){ #### fonction d'actualisation du param?tre de la variance des erreurs
  nind <- length(unique(id))
  Nombre <- length(id)
  sigm <- 0
  for (j in 1:nind){
    w <- which(id==unique(id)[j])
    V <- Z[w,, drop=FALSE]%*%Btilde%*%t(Z[w,, drop=FALSE])+diag(as.numeric(sigma),length(w),length(w))
    sigm <- sigm + t(epsilon[w])%*%epsilon[w] + sigma*(length(w)-sigma*(sum(diag(solve(V)))))
  }
  sigm <- sigm/Nombre
  return(sigm)
}

#' Title
#'
#' @import stats
#'
#' @keywords internal
bay <- function(bhat,Bhat,Z,id, sigmahat){ #### actualisation des param?tres de B
  nind <- length(unique(id))
  q <- dim(Z)[2]
  Nombre <- length(id)
  D <- 0
  for (j in 1:nind){
    w <- which(id==unique(id)[j])
    V <- Z[w,, drop=FALSE]%*%Bhat%*%t(Z[w,, drop=FALSE])+diag(as.numeric(sigmahat),length(w),length(w))
    D <- D+ (bhat[j,]%*%t(bhat[j,]))+ (Bhat- Bhat%*%t(Z[w,, drop=FALSE])%*%solve(V)%*%Z[w,, drop=FALSE]%*%Bhat)
  }
  D <- D/nind
  return(D)
}

#' Title
#'
#' @import stats
#'
#' @keywords internal
Moy <- function(id,Btilde,sigmahat,Phi,Y,Z){
  S1<- 0
  S2<- 0
  nind <- length(unique(id))
  for (i in 1:nind){
    w <- which(id==unique(id)[i])
    V <- Z[w,, drop=FALSE]%*%Btilde%*%t(Z[w,, drop=FALSE])+diag(as.numeric(sigmahat),length(w),length(w))
    S1 <- S1 + t(Phi[w,, drop=FALSE])%*%solve(V)%*%Phi[w,, drop=FALSE]
    S2 <- S2 + t(Phi[w,, drop = FALSE])%*%solve(V)%*%Y[w]
  }
  return(solve(S1)%*%S2)
}


#' (S)REEMforest alogorithm
#'
#' @param X :
#' @param Y :
#' @param id :
#' @param Z :
#' @param iter :
#' @param mtry :
#' @param ntree :
#' @param time :
#' @param sto :
#' @param delta :
#'
#' @import stats
#' @import randomForest
#'
#' @return
#' @export
#'
REEMforest <- function(X,Y,id,Z,iter,mtry,ntree, time, sto, delta = 0.001){
  q <- dim(Z)[2]
  nind <- length(unique(id))
  btilde <- matrix(0,nind,q) #### Pour la ligne i, on a les effets al?atoires de l'individu i
  sigmahat <- 1 #### init
  Btilde <- diag(rep(1,q)) ### init
  epsilonhat <- 0
  id_btilde <- unique(id)
  Tiime <- sort(unique(time))
  omega <- rep(0,length(Y))
  sigma2 <- 1
  Vrai <- NULL
  inc <- 1
  OOB <- NULL

  if (class(sto)=="character"){
    if (sto=="fbm"){

      id_omega <- matrix(0,nind,length(unique(time)))
      for (i in 1:length(unique(id))){
        w <- which(id ==id_btilde[i])
        time11 <- time[w]
        where <- NULL
        for (j in 1:length(time11)){
          where <- c(where,which(Tiime==time11[j]))
        }
        id_omega[i,where] <- 1
      }
      omega <- matrix(0,nind,length(unique(time)))
      omega2 <- rep(0,length(Y))
      h <- opti.FBMreem(X,Y,id,Z,iter, mtry,ntree,time)
      for (i in 1:iter){
        ystar <- rep(0,length(Y))
        for (k in 1:nind){ #### on retrace les effets al?atoires
          indiv <- which(id==unique(id)[k])
          ystar[indiv] <- Y[indiv]- Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
        }
        forest <- randomForest(X, ystar,mtry=mtry,ntree=ntree, importance = TRUE, keep.inbag=TRUE)
        f1 <- predict(forest,X,nodes=TRUE)
        OOB[i] <- forest$mse[ntree]
        trees <- attributes(f1)
        inbag <- forest$inbag
        matrice.pred <- matrix(NA,length(Y),ntree)

        for (k in 1:ntree){
          Phi <- matrix(0,length(Y),length(unique(trees$nodes[,k])))
          indii <- which(forest$forest$nodestatus[,k]==-1)
          for (l in 1:dim(Phi)[2]){
            w <- which(trees$nodes[,k]==indii[l])
            Phi[w,l] <- 1
          }
          oobags <- unique(which(inbag[,k]==0))
          beta <- Moy_fbm(id[-oobags],Btilde,sigmahat,Phi[-oobags,],ystar[-oobags],Z[-oobags,,drop=FALSE],h,time[-oobags], sigma2)
          forest$forest$nodepred[indii,k] <- beta
          matrice.pred[oobags,k] <- Phi[oobags,]%*%beta
        }

        fhat <- rep(NA,length(Y))
        for (k in 1:length(Y)){
          w <- which(is.na(matrice.pred[k,])==TRUE)
          fhat[k] <- mean(matrice.pred[k,-w])
        }

        for (k in 1:nind){ ### calcul des effets al?atoires par individu
          indiv <- which(id==unique(id)[k])
          K <- cov.fbm(time[indiv],h)
          V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,, drop=FALSE])+diag(as.numeric(sigmahat),length(indiv),length(indiv))+ sigma2*K
          btilde[k,] <- Btilde%*%t(Z[indiv,, drop=FALSE])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
          omega[k,which(id_omega[k,]==1)] <- sigma2*K%*%solve(V)%*%(Y[indiv]-fhat[indiv])
          omega2[indiv] <- omega[k,which(id_omega[k,]==1)]
          epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
        }
        sigm <- sigmahat
        B <- Btilde
        sigmahat <- sig.fbm(Y,sigmahat,id, Z, epsilonhat, Btilde, time, sigma2,h) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
        Btilde  <- bay.fbm(btilde,Btilde,Z,id,sigm, time, sigma2,h) #### MAJ des param?tres de la variance des effets al?atoires.
        ### MAJ de la volatilit? du processus stochastique
        sigma2 <- gam_fbm(Y,sigm,id,Z,B,time,sigma2,omega,id_omega,h)
        Vrai <- c(Vrai, logV.fbm(Y,fhat,Z[,,drop=FALSE],time,id,Btilde,sigma2,sigmahat,h))
        if (i>1) inc <- abs(Vrai[i-1]-Vrai[i])/abs(Vrai[i-1])
        if (inc< delta) {
          print(cat("stopped after", i, "iterations."))
          return(list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, vraisemblance = Vrai,id=id, time =time, Hurst=h, OOB =OOB, omega=omega2))
        }
      }
      return(list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto=sto,omega=omega2, sigma_sto =sigma2, time =time, sto= sto, Hurst =h, Vraisemblance=Vrai, OOB =OOB))
    }

    if ( sto=="none"){
      for (i in 1:iter){
        ystar <- rep(0,length(Y))
        for (k in 1:nind){ #### on retrace les effets al?atoires
          indiv <- which(id==unique(id)[k])
          ystar[indiv] <- Y[indiv]- Z[indiv,, drop=FALSE]%*%btilde[k,]
        }
        forest <- randomForest(X, ystar,mtry=mtry,ntree=ntree, importance = TRUE, keep.inbag=TRUE)
        f1 <- predict(forest,X,nodes=TRUE)
        trees <- attributes(f1)
        OOB[i] <- forest$mse[ntree]
        inbag <- forest$inbag
        matrice.pred <- matrix(NA,length(Y),ntree)


        for (k in 1:ntree){
          Phi <- matrix(0,length(Y),length(unique(trees$nodes[,k])))
          indii <- which(forest$forest$nodestatus[,k]==-1)
          for (l in 1:dim(Phi)[2]){
            w <- which(trees$nodes[,k]==indii[l])
            Phi[w,l] <- 1
          }
          oobags <- unique(which(inbag[,k]==0))
          beta <- Moy(id[-oobags],Btilde,sigmahat,Phi[-oobags,],ystar[-oobags],Z[-oobags,,drop=FALSE])
          forest$forest$nodepred[indii,k] <- beta
          matrice.pred[oobags,k] <- Phi[oobags,]%*%beta
        }

        fhat <- rep(NA,length(Y))
        for (k in 1:length(Y)){
          w <- which(is.na(matrice.pred[k,])==TRUE)
          fhat[k] <- mean(matrice.pred[k,-w])
        }

        for (k in 1:nind){ ### calcul des effets al?atoires par individu
          indiv <- which(id==unique(id)[k])
          V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,, drop=FALSE])+diag(as.numeric(sigmahat),length(indiv),length(indiv))
          btilde[k,] <- Btilde%*%t(Z[indiv,, drop=FALSE])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
          epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,, drop=FALSE]%*%btilde[k,]
        }

        sigm <- sigmahat
        sigmahat <- sig(sigmahat,id, Z, epsilonhat, Btilde) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
        Btilde  <- bay(btilde,Btilde,Z,id,sigm) #### MAJ des param?tres de la variance des effets al?atoires.
        Vrai <- c(Vrai, logV(Y,fhat,Z,time,id,Btilde,0,sigmahat,sto))
        if (i>1) inc <- abs((Vrai[i-1]-Vrai[i])/Vrai[i-1])
        if (inc< delta) {
          print(cat("stopped after", i, "iterations."))
          return(list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, vraisemblance = Vrai,id=id, time =time, OOB =OOB))
        }
      }
      return(list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, id = id , time = time , Vraisemblance=Vrai, OOB =OOB))
    }
  }
  for (i in 1:iter){

    ystar <- rep(0,length(Y))
    for (k in 1:nind){ #### on retrace les effets al?atoires
      indiv <- which(id==unique(id)[k])
      ystar[indiv] <- Y[indiv]- Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
    }

    forest <- randomForest(X, ystar,mtry=mtry,ntree=ntree, importance = TRUE, keep.inbag=TRUE)
    f1 <- predict(forest,X,nodes=TRUE)
    OOB[i] <- forest$mse[ntree]
    trees <- attributes(f1)
    inbag <- forest$inbag
    matrice.pred <- matrix(NA,length(Y),ntree)

    for (k in 1:ntree){
      Phi <- matrix(0,length(Y),length(unique(trees$nodes[,k])))
      indii <- which(forest$forest$nodestatus[,k]==-1)
      for (l in 1:dim(Phi)[2]){
        w <- which(trees$nodes[,k]==indii[l])
        Phi[w,l] <- 1
      }
      oobags <- unique(which(inbag[,k]==0))
      beta <- Moy_sto(id[-oobags],Btilde,sigmahat,Phi[-oobags,, drop=FALSE],ystar[-oobags],Z[-oobags,,drop=FALSE], sto, time[-oobags], sigma2)
      forest$forest$nodepred[indii,k] <- beta
      matrice.pred[oobags,k] <- Phi[oobags,]%*%beta
    }

    fhat <- rep(NA,length(Y))
    for (k in 1:length(Y)){
      w <- which(is.na(matrice.pred[k,])==TRUE)
      fhat[k] <- mean(matrice.pred[k,-w])
    }

    for (k in 1:nind){ ### calcul des effets al?atoires par individu
      indiv <- which(id==unique(id)[k])
      K <- sto_analysis(sto,time[indiv])
      V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,, drop=FALSE])+diag(as.numeric(sigmahat),length(indiv),length(indiv))+ sigma2*K
      btilde[k,] <- Btilde%*%t(Z[indiv,, drop=FALSE])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
      omega[indiv] <- sigma2*K%*%solve(V)%*%(Y[indiv]-fhat[indiv])
      epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
    }
    sigm <- sigmahat
    B <- Btilde
    sigmahat <- sig_sto(sigmahat,id, Z, epsilonhat, Btilde, time, sigma2,sto) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
    Btilde  <- bay_sto(btilde,Btilde,Z,id,sigm, time, sigma2,sto) #### MAJ des param?tres de la variance des effets al?atoires.
    ### MAJ de la volatilit? du processus stochastique
    sigma2 <- gam_sto(sigm,id,Z,B,time,sigma2,sto,omega)
    Vrai <- c(Vrai, logV(Y,fhat,Z[,,drop=FALSE],time,id,Btilde,sigma2,sigmahat,sto))
    if (i>1) {inc <- abs((Vrai[i-1]-Vrai[i])/Vrai[i-1])
    if (Vrai[i]<Vrai[i-1]) {reemfouille <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), omega=omega, sigma_sto =sigma2, time = time, sto= sto,Vraisemblance=Vrai,id=id, OOB =OOB)}
    }
    if (inc< delta) {
      print(cat("stopped after", i, "iterations."))
      return(reemfouille)
    }
  }
  return(list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id),omega=omega, sigma_sto =sigma2, time = time, sto= sto, id=id, OOB =OOB, Vraisemblance=Vrai))
}

#' Title
#'
#' @import stats
#'
#' @keywords internal
Moy_sto <- function(id,Btilde,sigmahat,Phi,Y,Z, sto, time, sigma2){
  S1<- 0
  S2<- 0
  nind <- length(unique(id))
  for (i in 1:nind){
    w <- which(id==unique(id)[i])
    K <- sto_analysis(sto,time[w])
    V <- Z[w,, drop=FALSE]%*%Btilde%*%t(Z[w,, drop=FALSE])+diag(as.numeric(sigmahat),length(w),length(w))+ sigma2*K
    S1 <- S1 + t(Phi[w,, drop=FALSE])%*%solve(V)%*%Phi[w,, drop=FALSE]
    S2 <- S2 + t(Phi[w,, drop=FALSE])%*%solve(V)%*%Y[w]
  }
  return(solve(S1)%*%S2)
}

#' Title
#'
#' @import stats
#'
#' @keywords internal
Moy_exp <- function(id,Btilde,sigmahat,Phi,Y,Z, alpha, time, sigma2){
  S1<- 0
  S2<- 0
  nind <- length(unique(id))
  for (i in 1:nind){
    w <- which(id==unique(id)[i])
    K <- cov.exp(time[w], alpha)
    V <- Z[w,, drop=FALSE]%*%Btilde%*%t(Z[w,, drop=FALSE])+diag(as.numeric(sigmahat),length(w),length(w))+ sigma2*K
    S1 <- S1 + t(Phi[w,, drop=FALSE])%*%solve(V)%*%Phi[w,, drop=FALSE]
    S2 <- S2 + t(Phi[w,, drop=FALSE])%*%solve(V)%*%Y[w]
  }
  return(solve(S1)%*%S2)
}


#' (S)MERT algorithm
#'
#' @param X :
#' @param Y :
#' @param id :
#' @param Z :
#' @param iter :
#' @param time :
#' @param sto :
#' @param delta :
#'
#' @import stats
#' @import rpart
#'
#' @return
#' @export
MERT <- function(X,Y,id,Z,iter,time, sto, delta = 0.001){
  q <- dim(Z)[2]
  nind <- length(unique(id))
  btilde <- matrix(0,nind,q) #### Pour la ligne i, on a les effets al?atoires de l'individu i
  sigmahat <- 1 #### init
  Btilde <- diag(rep(1,q)) ### init
  epsilonhat <- 0
  id_btilde <- unique(id)
  Tiime <- sort(unique(time))
  omega <- rep(0,length(Y))
  sigma2 <- 1
  inc <- 1
  Vrai <- NULL
  id_omega=sto

  if (class(sto)=="character"){
    if ( sto=="none"){
      for (i in 1:iter){
        ystar <- rep(0,length(Y))
        for (k in 1:nind){ #### on retrace les effets al?atoires
          indiv <- which(id==unique(id)[k])
          ystar[indiv] <- Y[indiv]- Z[indiv,, drop=FALSE]%*%btilde[k,]
        }

        tree <- rpart(ystar~.,as.data.frame(X)) ### on construit l'arbre
        fhat <- predict(tree, as.data.frame(X)) #### pr?diction avec l'arbre
        for (k in 1:nind){ ### calcul des effets al?atoires par individu
          indiv <- which(id==unique(id)[k])
          V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,, drop=FALSE])+diag(as.numeric(sigmahat),length(indiv),length(indiv))
          btilde[k,] <- Btilde%*%t(Z[indiv,, drop=FALSE])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
          epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,, drop=FALSE]%*%btilde[k,]
        }
        sigm <- sigmahat
        sigmahat <- sig(sigmahat,id, Z, epsilonhat, Btilde) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
        Btilde  <- bay(btilde,Btilde,Z,id,sigm) #### MAJ des param?tres de la variance des effets al?atoires.
        Vrai <- c(Vrai, logV(Y,fhat,Z[,,drop=FALSE],time,id,Btilde,0,sigmahat,"none"))
        if (i>1) inc <- (Vrai[i-1]-Vrai[i])/Vrai[i-1]
        if (inc< delta) {
          print(cat("stopped after", i, "iterations."))
          return(list(forest=tree,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, Vraisemblance = Vrai, id =id, time=time))
        }
      }
      return(list(forest=tree,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, Vraisemblance=Vrai, id=id, time=time))
    }
  }
  for (i in 1:iter){
    ystar <- rep(0,length(Y))
    for (k in 1:nind){ #### on retrace les effets al?atoires
      indiv <- which(id==unique(id)[k])
      ystar[indiv] <- Y[indiv]- Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
    }

    tree <- rpart(ystar~.,X) ### on construit l'arbre
    fhat <- predict(tree, X) #### pr?diction avec l'arbre
    for (k in 1:nind){ ### calcul des effets al?atoires par individu
      indiv <- which(id==unique(id)[k])
      K <- sto_analysis(sto,time[indiv])
      V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,, drop=FALSE])+diag(as.numeric(sigmahat),length(indiv),length(indiv))+ sigma2*K
      btilde[k,] <- Btilde%*%t(Z[indiv,, drop=FALSE])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
      omega[indiv] <- sigma2*K%*%solve(V)%*%(Y[indiv]-fhat[indiv]-Z[indiv,,drop=FALSE]%*%btilde[k,])
      epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
    }
    #### pr?diction du processus stochastique:
    sigm <- sigmahat
    B <- Btilde
    sigmahat <- sig_sto(sigmahat,id, Z, epsilonhat, Btilde, time, sigma2,sto) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
    Btilde  <- bay_sto(btilde,Btilde,Z,id,sigm, time, sigma2,sto) #### MAJ des param?tres de la variance des effets al?atoires.
    ### MAJ de la volatilit? du processus stochastique
    sigma2 <- gam_sto(sigm,id,Z,B,time,sigma2,sto,omega)
    Vrai <- c(Vrai, logV(Y,fhat,Z[,,drop=FALSE],time,id,Btilde,sigma2,sigmahat,sto))
    if (i>1) inc <- (Vrai[i-1]-Vrai[i])/Vrai[i-1]
    if (inc< delta) {
      print(cat("stopped after", i, "iterations."))
      return(list(forest=tree,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat,id_omega=id_omega, id_btilde=unique(id),omega=omega, sigma_sto =sigma2, time = time, sto= sto,Vraisemblance=Vrai, id = id))
    }
  }
  return(list(forest=tree,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id),id_omega=id_omega,omega=omega, sigma_sto =sigma2, time = time, sto= sto, Vraisemblance=Vrai, id=id))
}


#' (S)REEMtree algorithm
#'
#' @param X :
#' @param Y :
#' @param id :
#' @param Z :
#' @param iter :
#' @param time :
#' @param sto :
#' @param delta :
#'
#' @import stats
#' @import rpart
#'
#' @return
#' @export
REEMtree <- function(X,Y,id,Z,iter, time, sto, delta = 0.001){
  q <- dim(Z)[2]
  nind <- length(unique(id))
  btilde <- matrix(0,nind,q) #### Pour la ligne i, on a les effets al?atoires de l'individu i
  sigmahat <- 1 #### init
  Btilde <- diag(rep(1,q)) ### init
  epsilonhat <- 0
  id_btilde <- unique(id)
  Tiime <- sort(unique(time))
  omega <- rep(0,length(Y))
  sigma2 <- 1
  Vrai <- NULL
  inc <- 1
  id_omega=sto

  if (class(sto)=="character"){
    if ( sto=="none"){
      for (i in 1:iter){
        ystar <- rep(0,length(Y))
        for (k in 1:nind){ #### on retrace les effets al?atoires
          indiv <- which(id==unique(id)[k])
          ystar[indiv] <- Y[indiv]- Z[indiv,, drop=FALSE]%*%btilde[k,]
        }

        tree <- rpart(ystar~.,X)
        Phi <- matrix(0,length(Y), length(unique(tree$where)))
        feuilles <- predict(tree,X)
        leaf <- unique(feuilles)
        for (p in 1:length(leaf)){
          w <- which(feuilles==leaf[p])
          Phi[unique(w),p] <- 1
        }

        beta <- Moy(id,Btilde,sigmahat,Phi,Y,Z) ### fit des feuilles

        for (k in 1:length(unique(tree$where))){
          ou <- which(tree$frame[,5]==leaf[k])
          lee <- which(tree$frame[,1]=="<leaf>")
          w <- intersect(ou,lee)
          tree$frame[w,5] <- beta[k]
        }
        fhat <- predict(tree, as.data.frame(X))
        for (k in 1:nind){ ### calcul des effets al?atoires par individu
          indiv <- which(id==unique(id)[k])
          V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,, drop=FALSE])+diag(as.numeric(sigmahat),length(indiv),length(indiv))
          btilde[k,] <- Btilde%*%t(Z[indiv,, drop=FALSE])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
          epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,, drop=FALSE]%*%btilde[k,]
        }
        sigm <- sigmahat
        sigmahat <- sig(sigmahat,id, Z, epsilonhat, Btilde) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
        Btilde  <- bay(btilde,Btilde,Z,id,sigm) #### MAJ des param?tres de la variance des effets al?atoires.
        Vrai <- c(Vrai, logV(Y,fhat,Z[,,drop=FALSE],time,id,Btilde,0,sigmahat,sto))
        if (i>1) inc <- (Vrai[i-1]-Vrai[i])/Vrai[i-1]
        if (inc <  delta) {
          print(cat("stopped after", i, "iterations."))
          return(list(forest=tree,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, vraisemblance = Vrai,id=id, time=time))
        }
      }
      return(list(forest=tree,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, Vraisemblance=Vrai, time =time, id=id ))
    }
  }
  for (i in 1:iter){
    ystar <- rep(0,length(Y))
    for (k in 1:nind){ #### on retrace les effets al?atoires
      indiv <- which(id==unique(id)[k])
      ystar[indiv] <- Y[indiv]- Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
    }

    tree <- rpart(ystar~.,X)
    Phi <- matrix(0,length(Y), length(unique(tree$where)))
    feuilles <- predict(tree,X)
    leaf <- unique(feuilles)
    for (p in 1:length(leaf)){
      w <- which(feuilles==leaf[p])
      Phi[unique(w),p] <- 1
    }

    beta <- Moy_sto(id,Btilde,sigmahat,Phi,Y,Z,sto,time,sigma2) ### fit des feuilles

    for (k in 1:length(unique(tree$where))){
      ou <- which(tree$frame[,5]==leaf[k])
      lee <- which(tree$frame[,1]=="<leaf>")
      w <- intersect(ou,lee)
      tree$frame[w,5] <- beta[k]
    }
    fhat <- predict(tree, X)
    for (k in 1:nind){ ### calcul des effets al?atoires par individu
      indiv <- which(id==unique(id)[k])
      K <- sto_analysis(sto,time[indiv])
      V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,, drop=FALSE])+diag(as.numeric(sigmahat),length(indiv),length(indiv))+ sigma2*K
      btilde[k,] <- Btilde%*%t(Z[indiv,, drop=FALSE])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
      omega[indiv] <- sigma2*K%*%solve(V)%*%(Y[indiv]-fhat[indiv])
      epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
    }
    #### pr?diction du processus stochastique:
    sigm <- sigmahat
    B <- Btilde
    sigmahat <- sig_sto(sigmahat,id, Z, epsilonhat, Btilde, time, sigma2,sto) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
    Btilde  <- bay_sto(btilde,Btilde,Z,id,sigm, time, sigma2,sto) #### MAJ des param?tres de la variance des effets al?atoires.
    sigma2 <- gam_sto(sigm,id,Z,B,time,sigma2,sto,omega)
    Vrai <- c(Vrai, logV(Y,fhat,Z[,,drop=FALSE],time,id,Btilde,sigma2,sigmahat,sto))
    if (i>1) inc <- (Vrai[i-1]-Vrai[i])/Vrai[i-1]
    if (inc< delta) {
      print(cat("stopped after", i, "iterations."))
      return(list(forest=tree,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), id_omega=id_omega, omega=omega, sigma_sto =sigma2, time = time, sto= sto,Vraisemblance=Vrai,id=id))
    }
  }
  return(list(forest=tree,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), id_omega=id_omega,omega=omega, sigma_sto =sigma2, time = time, sto= sto, Vraisemblance=Vrai, id=id))
}

#' Title
#'
#'
#' @import stats
#'
#' @keywords internal
Moy_fbm <- function(id,Btilde,sigmahat,Phi,Y,Z, H, time, sigma2){
  S1<- 0
  S2<- 0
  nind <- length(unique(id))
  for (i in 1:nind){
    w <- which(id==unique(id)[i])
    K <- cov.fbm(time[w],H)
    V <- Z[w,, drop=FALSE]%*%Btilde%*%t(Z[w,, drop=FALSE])+diag(as.numeric(sigmahat),length(w),length(w))+ sigma2*K
    S1 <- S1 + t(Phi[w,, drop=FALSE])%*%solve(V)%*%Phi[w,, drop=FALSE]
    S2 <- S2 + t(Phi[w,, drop=FALSE])%*%solve(V)%*%Y[w]
  }
  M <- solve(S1)%*%S2
}

g <- function(s,t) { 1/abs(s-t)}
#### faire le mod?le lin?aire mixte stochastique ::::

#' (S)LMEM algorithm
#'
#' @param X [matrix]:
#' @param Y [vector]:
#' @param id [vector]:
#' @param Z [matrix]:
#' @param iter [numeric]:
#' @param time [vector]:
#' @param sto [numeric]:
#' @param delta [numeric]:
#'
#' @import stats
#'
#' @return
#' @export
Linmix <- function(X,Y,id,Z,iter, time, sto, delta = 0.001){
  q <- dim(Z)[2]
  nind <- length(unique(id))
  btilde <- matrix(0,nind,q) #### Pour la ligne i, on a les effets al?atoires de l'individu i
  sigmahat <- 1 #### init
  if (q==1) Btilde <- 1
  if (q>1) Btilde <- diag(rep(1,q)) ### init
  epsilonhat <- 0
  id_btilde <- unique(id)
  Tiime <- sort(unique(time))
  omega <- rep(0,length(Y))
  sigma2 <- 1
  inc <- 1
  Vrai <- NULL

  if (class(sto)=="character"){
    if ( sto=="none"){
      for (i in 1:iter){
        ystar <- rep(0,length(Y))
        for (k in 1:nind){ #### on retrace les effets al?atoires
          indiv <- which(id==unique(id)[k])
          ystar[indiv] <- Y[indiv]- Z[indiv,, drop=FALSE]%*%btilde[k,]
        }

        beta <- Moy(id,Btilde,sigmahat,X,Y,Z) ### fit des feuilles
        fhat <- X%*%beta
        for (k in 1:nind){ ### calcul des effets al?atoires par individu
          indiv <- which(id==unique(id)[k])
          V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,, drop=FALSE])+diag(as.numeric(sigmahat),length(indiv),length(indiv))
          btilde[k,] <- Btilde%*%t(Z[indiv,,drop=FALSE])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
        }

        for (k in 1:nind){
          indiv <- which(id==unique(id)[k])
          epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,, drop=FALSE]%*%btilde[k,]
        }
        sigm <- sigmahat
        sigmahat <- sig(sigmahat,id, Z, epsilonhat, Btilde) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
        Btilde  <- bay(btilde,Btilde,Z,id,sigm) #### MAJ des param?tres de la variance des effets al?atoires.

      }
      return(list(beta=beta,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto))
    }
  }
  for (i in 1:iter){
    ystar <- rep(0,length(Y))
    for (k in 1:nind){ #### on retrace les effets al?atoires
      indiv <- which(id==unique(id)[k])
      ystar[indiv] <- Y[indiv]- Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
    }


    beta <- Moy_sto(id,Btilde,sigmahat,X,Y,Z,sto,time,sigma2) ### fit des feuilles
    fhat <- X%*%beta
    for (k in 1:nind){
      indiv <- which(id==unique(id)[k])
      K <- sto_analysis(sto,time[indiv])
      V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,, drop=FALSE])+diag(as.numeric(sigmahat),length(indiv),length(indiv))+ sigma2*K
      btilde[k,] <- Btilde%*%t(Z[indiv,, drop=FALSE])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
    }
    #### pr?diction du processus stochastique:
    for (k in 1:length(unique(id))){
      indiv <- which(id==unique(id)[k])
      K <- sto_analysis(sto,time[indiv])
      V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,, drop=FALSE])+diag(as.numeric(sigmahat),length(indiv),length(indiv))+sigma2*K
      omega[indiv] <- sigma2*K%*%solve(V)%*%(Y[indiv]-fhat[indiv])
    }

    for (k in 1:nind){
      indiv <- which(id==unique(id)[k])
      epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
    }
    sigm <- sigmahat
    B <- Btilde
    sigmahat <- sig_sto(sigmahat,id, Z, epsilonhat, Btilde, time, sigma2,sto) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
    Btilde  <- bay_sto(btilde,Btilde,Z,id,sigm, time, sigma2,sto) #### MAJ des param?tres de la variance des effets al?atoires.
    sigma2 <- gam_sto(sigm,id,Z,B,time,sigma2,sto,omega)
    Vrai <- c(Vrai, logV(Y,fhat,Z[,,drop=FALSE],time,id,Btilde,sigma2,sigmahat,sto))
    if (i>1) inc <- (Vrai[i-1]-Vrai[i])/Vrai[i-1]
    if (inc< delta) {
      print(cat("stopped after", i, "iterations."))
      return(list(beta=beta,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id),omega=omega, sigma_sto =sigma2, time = Tiime, sto= sto,Vraisemblance=Vrai))
    }
  }
  return(list(beta=beta,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id),omega=omega, sigma_sto =sigma2, time = Tiime, sto= sto, Vraisemblance=Vrai))
}


#' Predict a fitted (S)LMEM algorithm
#'
#' @param lin :
#' @param X :
#' @param Z :
#' @param id :
#' @param time :
#'
#' @import stats
#'
#' @return
#' @export
#'
predict_lin <- function(lin, X,Z,id,time){
  n <- length(unique(id))
  Y <- rep(0,length(id))
  id_btilde <- lin$id_btilde
  f <- X%*%lin$beta
  Time <- lin$time
  id_btilde <- lin$id_btilde
  Ypred <- rep(0,length(id))
  if (lin$sto=="none"){
    for (i in 1:length(unique(id))){
      w <- which(id==unique(id)[i])
      k <- which(id_btilde==unique(id)[i])
      Ypred[w] <- f[w] + Z[w,, drop=FALSE]%*%lin$random_effects[k,]
    }
    return(Ypred)
  }

  if (lin$sto=="exp"){
    for (i in 1:length(unique(id))){
      w <- which(id==unique(id)[i])
      k <- which(id_btilde==unique(id)[i])
      om <- which(lin$id_omega[k,]==1)
      Ypred[w] <- f[w] + Z[w,, drop=FALSE]%*%lin$random_effects[k,] + predict.exp(lin$omega[k,om],Time[om],time[w], lin$alpha)
    }
    return(Ypred)
  }

  if (lin$sto=="fbm"){
    for (i in 1:length(unique(id))){
      w <- which(id==unique(id)[i])
      k <- which(id_btilde==unique(id)[i])
      om <- which(lin$id_omega[k,]==1)
      Ypred[w] <- f[w] + Z[w,, drop=FALSE]%*%lin$random_effects[k,] + predict.fbm(lin$omega[k,om],Time[om],time[w], lin$Hurst)
    }
    return(Ypred)
  }

  for (i in 1:length(unique(id))){
    w <- which(id==unique(id)[i])
    k <- which(id_btilde==unique(id)[i])
    om <- which(lin$id_omega[k,]==1)
    Ypred[w] <- f[w] + Z[w,, drop=FALSE]%*%lin$random_effects[k,] + predict.sto(lin$omega[k,om],Time[om],time[w], lin$sto)
  }
  return(Ypred)
}

#' Title
#'
#'
#' @import stats
#'
#' @keywords internal
logV <- function(Y,f,Z,time,id,B,gamma,sigma, sto){
  Vraisem <- 0
  if (sto=="none"){
    for (i in 1:length(unique(id))){
      w <- which(id==unique(id)[i])
      V <- Z[w,,drop=FALSE]%*%B%*%t(Z[w,,drop=FALSE])+diag(as.numeric(sigma),length(w),length(w))
      Vraisem <- Vraisem + log(det(V))+ t(Y[w]-f[w])%*%solve(V)%*%(Y[w]-f[w])
    }
    return(Vraisem)
  }
  for (i in 1:length(unique(id))){
    w <- which(id==unique(id)[i])
    K <- sto_analysis(sto,time[w])
    V <- Z[w,,drop=FALSE]%*%B%*%t(Z[w,,drop=FALSE])+gamma*K+ diag(as.numeric(sigma),length(w),length(w))
    Vraisem <- Vraisem + log(det(V))+ t(Y[w]-f[w])%*%solve(V)%*%(Y[w]-f[w])
  }
  return(Vraisem)
}



#' Title
#'
#' @param time
#' @param H
#'
#' @keywords internal
cov.fbm <- function(time,H){
  K <- matrix(0,length(time),length(time))
  for (i in 1:length(time)){
    for (j in 1:length(time)){
      K[i,j] <- 0.5*(time[i]^(2*H)+time[j]^(2*H)-abs(time[i]-time[j])^(2*H))
    }
  }
  return(K)
}

#' Title
#'
#' @param X
#' @param Y
#' @param id
#' @param Z
#' @param iter
#' @param mtry
#' @param ntree
#' @param time
#'
#' @import stats
#'
#' @keywords internal
opti.FBM <- function(X,Y,id,Z,iter,mtry,ntree,time){
  print("Do you want to enter an ensemble for the Hurst parameter ? (1/0)")
  resp <- scan(nmax=1)
  if (resp ==1){
    print("please enter your ensemble (vector):")
    H <- scan()
    opti <- NULL}

  if (resp==0) {H <- seq(0.1,0.9,0.1)}
  opti <- NULL
  for (h in H){
    q <- dim(Z)[2]
    nind <- length(unique(id))
    btilde <- matrix(0,nind,q) #### Pour la ligne i, on a les effets al?atoires de l'individu i
    sigmahat <- 1 #### init
    Btilde <- diag(rep(1,q)) ### init
    epsilonhat <- 0
    id_btilde <- unique(id)
    Tiime <- sort(unique(time))
    omega <- matrix(0,nind,length(unique(time)))
    mu = rep(0,length(id_btilde))
    sigma2 <- 1
    id_omega <- matrix(0,nind,length(unique(time)))
    for (i in 1:length(unique(id))){
      w <- which(id ==id_btilde[i])
      time11 <- time[w]
      where <- NULL
      for (j in 1:length(time11)){
        where <- c(where,which(Tiime==time11[j]))
      }
      id_omega[i,where] <- 1
    }
    for (i in 1:iter){
      ystar <- rep(0,length(Y))
      for (k in 1:nind){ #### on retrace les effets al?atoires
        indiv <- which(id==unique(id)[k])
        ystar[indiv] <- Y[indiv]- Z[indiv,]%*%btilde[k,]- omega[k,which(id_omega[k,]==1)]
      }

      forest <- randomForest(as.data.frame(X),ystar,mtry=mtry,ntree=ntree) ### on construit l'arbre
      fhat <- predict(forest, as.data.frame(X)) #### pr?diction avec l'arbre
      for (k in 1:nind){ ### calcul des effets al?atoires par individu
        indiv <- which(id==unique(id)[k])
        K <- cov.fbm(time[indiv], h)
        V <- Z[indiv,]%*%Btilde%*%t(Z[indiv,])+diag(rep(sigmahat,length(Y[indiv])))+ sigma2*K
        btilde[k,] <- Btilde%*%t(Z[indiv,])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
      }
      #### pr?diction du processus stochastique:
      for (k in 1:length(id_btilde)){
        indiv <- which(id==unique(id)[k])
        K <- cov.fbm(time[indiv], h)
        V <- Z[indiv,]%*%Btilde%*%t(Z[indiv,])+diag(rep(sigmahat,length(Y[indiv])))+sigma2*K
        omega[k,which(id_omega[k,]==1)] <- sigma2*K%*%solve(V)%*%(Y[indiv]-fhat[indiv])
      }

      for (k in 1:nind){
        indiv <- which(id==unique(id)[k])
        epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,]%*%btilde[k,]- omega[k,which(id_omega[k,]==1)]
      }
      sigm <- sigmahat
      B <- Btilde
      sigmahat <- sig.fbm(Y,sigmahat,id, Z, epsilonhat, Btilde, time, sigma2,h) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
      Btilde  <- bay.fbm(btilde,Btilde,Z,id,sigm, time, sigma2,h) #### MAJ des param?tres de la variance des effets al?atoires.
      sigma2 <- gam_fbm(Y,sigm,id,Z,B,time,sigma2,omega,id_omega, h)
    }
    opti <- c(opti,logV.fbm(Y, fhat, Z, time, id, Btilde, sigma2,sigmahat,h))
  }
  return(H[which.max(opti)])
}


#' Title
#'
#' @param bhat
#' @param Bhat
#' @param Z
#' @param id
#' @param sigmahat
#' @param time
#' @param sigma2
#' @param h
#'
#' @import stats
#'
#' @keywords internal
bay.fbm <- function(bhat,Bhat,Z,id, sigmahat, time, sigma2,h){ #### actualisation des param?tres de B
  nind <- length(unique(id))
  q <- dim(Z)[2]
  Nombre <- length(id)
  D <- 0
  for (j in 1:nind){
    w <- which(id==unique(id)[j])
    K <- cov.fbm(time[w], h)
    V <- Z[w,]%*%Bhat%*%t(Z[w,])+diag(rep(sigmahat,length(w)))+sigma2*K
    D <- D+ (bhat[j,]%*%t(bhat[j,]))+ (Bhat- Bhat%*%t(Z[w,])%*%solve(V)%*%Z[w,]%*%Bhat)
  }
  D <- D/nind
  return(D)
}

#' Title
#'
#' @param sigma
#' @param id
#' @param Z
#' @param epsilon
#' @param Btilde
#' @param time
#' @param sigma2
#' @param h
#'
#' @import stats
#'
#' @keywords internal
sig.fbm <- function(Y,sigma,id,Z, epsilon, Btilde, time, sigma2,h){ #### fonction d'actualisation du param?tre de la variance des erreurs
  nind <- length(unique(id))
  Nombre <- length(id)
  sigm <- 0
  for (j in 1:nind){
    w <- which(id==unique(id)[j])
    K <- cov.fbm(time[w], h)
    V <- Z[w,]%*%Btilde%*%t(Z[w,])+diag(rep(sigma,length(Y[w])))+sigma2*K
    sigm <- sigm + t(epsilon[w])%*%epsilon[w] + sigma*(length(w)-sigma*(sum(diag(solve(V)))))
  }
  sigm <- sigm/Nombre
  return(sigm)
}

#' Title
#'
#' @param sigm
#' @param id
#' @param Z
#' @param B
#' @param time
#' @param sigma2
#' @param omega
#' @param id_omega
#' @param h
#'
#' @import stats
#'
#' @keywords internal
gam_fbm <- function(Y,sigm,id,Z,Btilde,time,sigma2,omega,id_omega,h){
  nind <- length(unique(id))
  Nombre <- length(id)
  gam <- 0
  for (k in 1:nind){
    indiv <- which(id==unique(id)[k])
    K <- cov.fbm(time[indiv],h)
    V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,,drop=FALSE])+diag(rep(sigm,length(Y[indiv])))+ sigma2*K
    Omeg <- omega[k,which(id_omega[k,]==1)]
    gam <-gam+ (t(Omeg)%*%solve(K)%*%Omeg) + sigma2*(length(indiv)-sigma2*sum(diag(solve(V)%*%K)))
  }
  return(as.numeric(gam)/Nombre)
}

#' Title
#'
#' @param sigm
#' @param id
#' @param Z
#' @param B
#' @param time
#' @param sigma2
#' @param omega
#' @param id_omega
#' @param alpha
#'
#' @import stats
#'
#' @keywords internal
gam_exp <- function(Y,sigm,id,Z,Btilde,time,sigma2,omega,id_omega, alpha){
  nind <- length(unique(id))
  Nombre <- length(id)
  gam <- 0
  for (k in 1:nind){
    indiv <- which(id==unique(id)[k])
    K <- cov.exp(time[indiv],alpha)
    V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,,drop=FALSE])+diag(rep(sigm,length(Y[indiv])))+ sigma2*K
    Omeg <- omega[k,which(id_omega[k,]==1)]
    gam <-gam+ (t(Omeg)%*%solve(K)%*%Omeg) + sigma2*(length(indiv)-sigma2*sum(diag(solve(V)%*%K)))
  }
  return(as.numeric(gam)/Nombre)
}

#' Title
#'
#' @param time
#' @param alpha
#'
#' @import stats
#'
#' @keywords internal
cov.exp <- function(time,alpha){
  K <- matrix(0,length(time),length(time))
  for (i in 1:length(time)){
    for (j in 1:length(time)){
      K[i,j] <- exp(-(alpha*abs(time[i]-time[j])))
    }
  }
  return(K)
}

#' Title
#'
#' @param bhat
#' @param Bhat
#' @param Z
#' @param id
#' @param sigmahat
#' @param time
#' @param sigma2
#' @param alpha
#'
#'  @import stats
#'
#' @keywords internal
bay.exp <- function(bhat,Bhat,Z,id, sigmahat, time, sigma2,alpha){ #### actualisation des param?tres de B
  nind <- length(unique(id))
  q <- dim(Z)[2]
  Nombre <- length(id)
  D <- 0
  for (j in 1:nind){
    w <- which(id==unique(id)[j])
    K <- cov.exp(time[w], alpha)
    V <- Z[w,]%*%Bhat%*%t(Z[w,])+diag(rep(sigmahat,length(w)))+sigma2*K
    D <- D+ (bhat[j,]%*%t(bhat[j,]))+ (Bhat- Bhat%*%t(Z[w,])%*%solve(V)%*%Z[w,]%*%Bhat)
  }
  D <- D/nind
  return(D)
}

#' Title
#'
#' @param sigma
#' @param id
#' @param Z
#' @param epsilon
#' @param Btilde
#' @param time
#' @param sigma2
#' @param alpha
#'
#' @import stats
#'
#' @keywords internal
sig.exp <- function(Y,sigma,id,Z, epsilon, Btilde, time, sigma2,alpha){ #### fonction d'actualisation du param?tre de la variance des erreurs
  nind <- length(unique(id))
  Nombre <- length(id)
  sigm <- 0
  for (j in 1:nind){
    w <- which(id==unique(id)[j])
    K <- cov.exp(time[w], alpha)
    V <- Z[w,]%*%Btilde%*%t(Z[w,])+diag(rep(sigma,length(Y[w])))+sigma2*K
    sigm <- sigm + t(epsilon[w])%*%epsilon[w] + sigma*(length(w)-sigma*(sum(diag(solve(V)))))
  }
  sigm <- sigm/Nombre
  return(sigm)
}

#' Title
#'
#' @param X
#' @param Y
#' @param id
#' @param Z
#' @param iter
#' @param mtry
#' @param ntree
#' @param time
#'
#' @import stats
#'
#' @keywords internal
opti.exp <- function(X,Y,id,Z,iter,mtry,ntree,time){
  print("Do you want to enter a set for the alpha parameter ? (1/0)")
  resp <- scan(nmax=1)
  if (resp ==1){
    print("please enter your set (vector):")
    alpha <- scan()
    opti <- NULL}

  if (resp==0) {alpha <- seq(0.1,0.9,0.05)}
  opti <- NULL
  for (al in alpha){
    q <- dim(Z)[2]
    nind <- length(unique(id))
    btilde <- matrix(0,nind,q) #### Pour la ligne i, on a les effets al?atoires de l'individu i
    sigmahat <- 1 #### init
    Btilde <- diag(rep(1,q)) ### init
    epsilonhat <- 0
    id_btilde <- unique(id)
    Tiime <- sort(unique(time))
    omega <- matrix(0,nind,length(unique(time)))
    sigma2 <- 1
    id_omega <- matrix(0,nind,length(unique(time)))
    for (i in 1:length(unique(id))){
      w <- which(id ==id_btilde[i])
      time11 <- time[w]
      where <- NULL
      for (j in 1:length(time11)){
        where <- c(where,which(Tiime==time11[j]))
      }
      id_omega[i,where] <- 1
    }
    for (i in 1:iter){
      ystar <- rep(0,length(Y))
      for (k in 1:nind){ #### on retrace les effets al?atoires
        indiv <- which(id==unique(id)[k])
        ystar[indiv] <- Y[indiv]- Z[indiv,]%*%btilde[k,]- omega[k,which(id_omega[k,]==1)]
      }

      forest <- randomForest(as.data.frame(X),ystar,mtry=mtry,ntree=ntree) ### on construit l'arbre
      fhat <- predict(forest, as.data.frame(X)) #### pr?diction avec l'arbre
      for (k in 1:nind){ ### calcul des effets al?atoires par individu
        indiv <- which(id==unique(id)[k])
        K <- cov.exp(time[indiv], al)
        V <- Z[indiv,]%*%Btilde%*%t(Z[indiv,])+diag(rep(sigmahat,length(Y[indiv])))+ sigma2*K
        btilde[k,] <- Btilde%*%t(Z[indiv,])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
      }
      #### pr?diction du processus stochastique:
      for (k in 1:length(id_btilde)){
        indiv <- which(id==unique(id)[k])
        K <- cov.exp(time[indiv], al)
        V <- Z[indiv,]%*%Btilde%*%t(Z[indiv,])+diag(rep(sigmahat,length(Y[indiv])))+sigma2*K
        omega[k,which(id_omega[k,]==1)] <- sigma2*K%*%solve(V)%*%(Y[indiv]-fhat[indiv])
      }

      for (k in 1:nind){
        indiv <- which(id==unique(id)[k])
        epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,]%*%btilde[k,]- omega[k,which(id_omega[k,]==1)]
      }
      sigm <- sigmahat
      B <- Btilde
      sigmahat <- sig.exp(Y,sigmahat,id, Z, epsilonhat, Btilde, time, sigma2,al) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
      Btilde  <- bay.exp(btilde,Btilde,Z,id,sigm, time, sigma2,al) #### MAJ des param?tres de la variance des effets al?atoires.
      ### MAJ de la volatilit? du processus stochastique
      sigma2 <- gam_exp(Y,sigm,id,Z,B,time,sigma2,omega,id_omega,al)
    }
    opti <- c(opti,sigmahat)
  }
  return(alpha[which.min(opti)])
}


#' Title
#'
#' @param X
#' @param Y
#' @param id
#' @param Z
#' @param iter
#' @param mtry
#' @param ntree
#' @param time
#'
#' @import stats
#'
#' @keywords internal
opti.FBMreem <- function(X,Y,id,Z,iter,mtry,ntree,time){
  print("Do you want to enter a set for the Hurst parameter ? (1/0)")
  resp <- scan(nmax=1)
  if (resp ==1){
    print("please enter your ensemble (vector):")
    H <- scan()
    opti <- NULL}

  if (resp==0) {H <- seq(0.1,0.9,0.1)}
  opti <- NULL
  for (h in H){
    q <- dim(Z)[2]
    nind <- length(unique(id))
    btilde <- matrix(0,nind,q) #### Pour la ligne i, on a les effets al?atoires de l'individu i
    sigmahat <- 1 #### init
    Btilde <- diag(rep(1,q)) ### init
    epsilonhat <- 0
    id_btilde <- unique(id)
    Tiime <- sort(unique(time))
    omega <- matrix(0,nind,length(unique(time)))
    mu = rep(0,length(id_btilde))
    sigma2 <- 1
    id_omega <- matrix(0,nind,length(unique(time)))
    for (i in 1:length(unique(id))){
      w <- which(id ==id_btilde[i])
      time11 <- time[w]
      where <- NULL
      for (j in 1:length(time11)){
        where <- c(where,which(Tiime==time11[j]))
      }
      id_omega[i,where] <- 1
    }
    for (i in 1:iter){
      ystar <- rep(0,length(Y))
      for (k in 1:nind){ #### on retrace les effets al?atoires
        indiv <- which(id==unique(id)[k])
        ystar[indiv] <- Y[indiv]- Z[indiv,]%*%btilde[k,]- omega[k,which(id_omega[k,]==1)]
      }
      forest <- randomForest(as.data.frame(X), ystar,mtry=mtry,ntree=ntree)
      f1 <- predict(forest, as.data.frame(X), nodes=TRUE)
      trees <- attributes(f1)
      K <- ntree

      for (k in 1:K){
        Phi <- matrix(0,length(Y),length(unique(trees$nodes[,k])))
        indii <- which(forest$forest$nodestatus[,k]==-1)
        for (l in 1:dim(Phi)[2]){
          w <- which(trees$nodes[,k]==indii[l])
          Phi[w,l] <- 1
        }
        beta <- Moy_fbm(id,Btilde,sigmahat,Phi,Y,Z, h ,time, sigma2)
        forest$forest$nodepred[indii,k] <- beta
      }

      fhat <- predict(forest,as.data.frame(X))

      for (k in 1:nind){ ### calcul des effets al?atoires par individu
        indiv <- which(id==unique(id)[k])
        K <- cov.fbm(time[indiv],h)
        V <- Z[indiv,]%*%Btilde%*%t(Z[indiv,])+diag(rep(sigmahat,length(Y[indiv])))+ sigma2*K
        btilde[k,] <- Btilde%*%t(Z[indiv,])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
      }
      #### pr?diction du processus stochastique:
      for (k in 1:length(id_btilde)){
        indiv <- which(id==unique(id)[k])
        K <- cov.fbm(time[indiv], h)
        V <- Z[indiv,]%*%Btilde%*%t(Z[indiv,])+diag(rep(sigmahat,length(Y[indiv])))+sigma2*K
        omega[k,which(id_omega[k,]==1)] <- sigma2*K%*%solve(V)%*%(Y[indiv]-fhat[indiv])
      }

      for (k in 1:nind){
        indiv <- which(id==unique(id)[k])
        epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,]%*%btilde[k,]- omega[k,which(id_omega[k,]==1)]
      }
      sigm <- sigmahat
      B <- Btilde
      sigmahat <- sig.fbm(Y,sigmahat,id, Z, epsilonhat, Btilde, time, sigma2,h) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
      Btilde  <- bay.fbm(btilde,Btilde,Z,id,sigm, time, sigma2,h) #### MAJ des param?tres de la variance des effets al?atoires.
      ### MAJ de la volatilit? du processus stochastique
      sigma2 <- gam_fbm(Y,sigm,id,Z,B,time,sigma2,omega,id_omega, h)
    }
    opti <- c(opti,sigmahat)
  }
  return(H[which.min(opti)])
}


#' Title
#'
#' @param Y
#' @param Z
#' @param id
#' @param fhat
#' @param sigmahat
#' @param sigma2
#' @param H
#' @param B
#'
#' @import stats
#'
#' @keywords internal
logV.fbm <- function(Y, fhat, Z, time, id, B, sigma2,sigmahat,H){
  logl <- 0
  for (i in 1:length(unique(id))){
    indiv <- which(id==unique(id)[i])
    K <- cov.fbm(time[indiv], H)
    V <- Z[indiv,]%*%B%*%t(Z[indiv,])+ sigma2*K+ as.numeric(sigmahat)*diag(rep(1,length(indiv)))
    logl <- logl + log(det(V)) + t(Y[indiv]-fhat[indiv])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
  }
  return(-logl)
}

#' Title
#'
#' @param Y
#' @param Z
#' @param id
#' @param fhat
#' @param sigmahat
#' @param sigma2
#' @param H
#' @param B
#'
#' @import stats
#'
#' @keywords internal
logV.exp <- function(Y, fhat, Z, time, id, B, sigma2,sigmahat,H){
  logl <- 0
  for (i in 1:length(unique(id))){
    indiv <- which(id==unique(id)[i])
    K <- cov.exp(time[indiv], H)
    V <- Z[indiv,]%*%B%*%t(Z[indiv,])+ sigma2*K+ as.numeric(sigmahat)*diag(rep(1,length(indiv)))
    logl <- logl + log(det(V)) + t(Y[indiv]-fhat[indiv])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
  }
  return(-logl)
}
