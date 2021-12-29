#' (S)MERF algorithm
#'
#'
#' (S)MERF is an adaptation of the random forest regression method to longitudinal data introduced by Hajjem et. al. (2014) <doi:10.1080/00949655.2012.741599>.
#' The model has been improved by Capitaine et. al. (2020) <doi:10.1177/0962280220946080> with the addition of a stochastic process.
#' The algorithm will estimate the parameters of the following semi-parametric stochastic mixed-effects model: \deqn{Y_i(t)=f(X_i(t))+Z_i(t)\beta_i + \omega_i(t)+\epsilon_i}
#' with \eqn{Y_i(t)} the output at time \eqn{t} for the \eqn{i}th individual; \eqn{X_i(t)} the input predictors (fixed effects) at time \eqn{t} for the \eqn{i}th individual;
#' \eqn{Z_i(t)} are the random effects at time \eqn{t} for the \eqn{i}th individual; \eqn{\omega_i(t)} is the stochastic process at time \eqn{t} for the \eqn{i}th individual
#'  which model the serial correlations of the output measurements; \eqn{\epsilon_i} is the residual error.
#'
#' @param X [matrix]: A \code{N}x\code{p} matrix containing the \code{p} predictors of the fixed effects, column codes for a predictor.
#' @param Y [vector]: A vector containing the output trajectories.
#' @param id [vector]: Is the vector of the identifiers for the different trajectories.
#' @param Z [matrix]: A \code{N}x\code{q} matrix containing the \code{q} predictor of the random effects.
#' @param iter [numeric]: Maximal number of iterations of the algorithm. The default is set to \code{iter=100}
#' @param mtry [numeric]: Number of variables randomly sampled as candidates at each split. The default value is \code{p/3}.
#' @param ntree [numeric]: Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times. The default value is \code{ntree=500}.
#' @param time [vector]: Is the vector of the measurement times associated with the trajectories in \code{Y},\code{Z} and \code{X}.
#' @param sto [character]: Defines the covariance function of the stochastic process, can be either \code{"none"} for no stochastic process, \code{"BM"} for Brownian motion, \code{OrnUhl} for standard Ornstein-Uhlenbeck process, \code{BBridge} for Brownian Bridge, \code{fbm} for Fractional Brownian motion; can also be a function defined by the user.
#' @param delta [numeric]: The algorithm stops when the difference in log likelihood between two iterations is smaller than \code{delta}. The default value is set to O.O01
#' @param conditional [logical]: Determines if the random forest algorithm is the one implemented by Breiman (2001) and available in \code{\link[randomForest]{randomForest}}, which is the default \code{FALSE}, or if it is the implemented by Strobl et al. (2007) using conditional inference trees in \code{\link[party]{cforest}}.
#'
#' @import randomForest
#' @import party
#' @import stats
#'
#' @return A fitted (S)MERF model which is a list of the following elements: \itemize{
#' \item \code{forest:} Random forest obtained at the last iteration.
#' \item \code{random_effects :} predictions of random effects for different trajectories.
#' \item \code{id_btilde:} Identifiers of individuals associated with the predictions \code{random_effects}.
#' \item \code{var_random_effects: } Estimation of the variance covariance matrix of random effects.
#' \item \code{sigma_sto: } Estimation of the volatility parameter of the stochastic process.
#' \item \code{sigma: } Estimation of the residual variance parameter.
#' \item \code{time: } The vector of the measurement times associated with the trajectories in \code{Y},\code{Z} and \code{X}.
#' \item \code{sto: } Stochastic process used in the model.
#' \item \code{Vraisemblance:} Log-likelihood of the different iterations.
#' \item \code{id: } Vector of the identifiers for the different trajectories.
#' \item \code{OOB: } OOB error of the fitted random forest at each iteration, only when \code{conditional=FALSE}.
#' }
#'
#' @export
#'
#'
#' @references
#' Ahlem Hajjem, François Bellavance, and Denis Larocque (2014). Mixed-effects random forest for clustered data. Journal of Statistical Computation and Simulation, 84(6), 1313–1328. \doi{10.1080/00949655.2012.741599}
#'
#' Louis Capitaine, Robin Genuer, and Rodolphe Thiébaut (2020). Random forests for high-dimensional longitudinal data. Statistical Methods in Medical Research, 096228022094608. \doi{10.1177/0962280220946080}
#'
#' Leo Breiman (2001). Random Forests. Machine Learning, 45(1), 5–32.
#'
#' Carolin Strobl, Anne-Laure Boulesteix, Achim Zeileis and Torsten Hothorn (2007). Bias in Random Forest Variable Importance Measures: Illustrations, Sources and a Solution. BMC Bioinformatics, 8(25).
#'
#'
#' @examples
#' set.seed(123)
#' data <- DataLongGenerator(n=20) # Generate the data composed by n=20 individuals.
#' # Train a SMERF model on the generated data. Should take ~ 50 seconds
#' # The data are generated with a Brownian motion,
#' # so we use the parameter sto="BM" to specify a Brownian motion as stochastic process
#' smerf <- MERF(X=data$X,Y=data$Y,Z=data$Z,id=data$id,time=data$time,mtry=2,ntree=500,sto="BM")
#' smerf$forest # is the fitted random forest (obtained at the last iteration).
#' smerf$random_effects # are the predicted random effects for each individual.
#' smerf$omega # are the predicted stochastic processes.
#' plot(smerf$Vraisemblance) # evolution of the log-likelihood.
#' smerf$OOB # OOB error at each iteration.
#' csmerf <- MERF(data$X,data$Y,Z=data$Z,id=data$id,time=data$time,mtry=2,sto="BM",conditional=TRUE)
#' csmerf$forest # is the fitted random forest (obtained at the last iteration).
#' csmerf$random_effects # are the predicted random effects for each individual.
#' csmerf$omega # are the predicted stochastic processes.
#' plot(csmerf$Vraisemblance) # evolution of the log-likelihood.
#'
#'
MERF <- function(X,Y,id,Z,iter=100,mtry=ceiling(ncol(X)/3),ntree=500, time, sto, delta = 0.001, conditional = FALSE){
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

	if(!is.logical(conditional) || length(conditional) != 1 || is.na(conditional)){
		conditional <- FALSE
	}
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
			h <- opti.FBM(X,Y,id,Z,iter, mtry,ntree,time, conditional)
			for (i in 1:iter){
				ystar <- rep(0,length(Y))
				for (k in 1:nind){ #### on retrace les effets al?atoires
					indiv <- which(id==unique(id)[k])
					ystar[indiv] <- Y[indiv]- Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
				}

				if(!conditional){
					forest <- randomForest(as.data.frame(X),ystar,mtry=mtry,ntree=ntree, importance = TRUE) ### on construit l'arbre
					fhat <- predict(forest) #### pr?diction avec l'arbre
					OOB[i] <- forest$mse[ntree]
				} else{
					citdata <- cbind(ystar, as.data.frame(X))
					forest <- cforest(ystar ~ ., data = citdata, controls = cforest_unbiased(mtry = mtry, ntree = ntree))
					fhat <- predict(forest, OOB = TRUE, type = "response") #### pr?diction avec l'arbre
				}
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
					print(paste0("stopped after ", i, " iterations."))
					if(!conditional){
						sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat,sigma_sto=sigma2, id_btilde=unique(id), sto= sto, vraisemblance = Vrai,id=id, time =time, Hurst=h, OOB =OOB, omega=omega2)
					} else {
						sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat,sigma_sto=sigma2, id_btilde=unique(id), sto= sto, vraisemblance = Vrai,id=id, time =time, Hurst=h, omega=omega2)
					}
					class(sortie)<-"longituRF"
					return(sortie)
				}
			}
			if(!conditional){
				sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id),sigma_sto=sigma2,omega=omega2, sigma_sto =sigma2, time = time, sto= sto, Hurst =h, id=id, Vraisemblance=Vrai, OOB =OOB)
			} else {
				sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id),sigma_sto=sigma2,omega=omega2, sigma_sto =sigma2, time = time, sto= sto, Hurst =h, id=id, Vraisemblance=Vrai)
			}
			class(sortie) <- "longituRF"
			return(sortie)
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

				if(!conditional){
					forest <- randomForest(as.data.frame(X),ystar,mtry=mtry,ntree=ntree, importance = TRUE) ### on construit l'arbre
					fhat <- predict(forest) #### pr?diction avec l'arbre
					OOB[i] <- forest$mse[ntree]
				} else{
					citdata <- cbind(ystar, as.data.frame(X))
					forest <- cforest(ystar ~ ., data = citdata, controls = cforest_unbiased(mtry = mtry, ntree = ntree))
					fhat <- predict(forest, OOB = TRUE, type = "response") #### pr?diction avec l'arbre
				}
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
					print(paste0("stopped after ", i, " iterations."))
					if(!conditional){
						sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, vraisemblance = Vrai,id=id, time=time, alpha = alpha, OOB =OOB, omega=omega2)
					} else {
						sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, vraisemblance = Vrai,id=id, time=time, alpha = alpha, omega=omega2)
					}
					class(sortie) <- "longituRF"
					return(sortie)
				}
			}
			if(!conditional){
				sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), omega=omega2, sigma_sto =sigma2, time = time, sto= sto, alpha=alpha, id=id, Vraisemblance=Vrai, OOB =OOB)
			}else{
				sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), omega=omega2, sigma_sto =sigma2, time = time, sto= sto, alpha=alpha, id=id, Vraisemblance=Vrai)
			}
			class(sortie) <- "longituRF"
			return(sortie)
		}

		if(sto == "none"){
			for(i in 1:iter){
				ystar <- rep(NA,length(Y))
				for (k in 1:nind){ #### on retrace les effets al?atoires
					indiv <- which(id==unique(id)[k])
					ystar[indiv] <- Y[indiv]- Z[indiv,,drop=FALSE]%*%btilde[k,]
				}

				if(!conditional){
					forest <- randomForest(as.data.frame(X),ystar,mtry=mtry,ntree=ntree, importance = TRUE) ### on construit l'arbre
					fhat <- predict(forest) #### pr?diction avec l'arbre
					OOB[i] <- forest$mse[ntree]
				}else{
					citdata <- cbind(ystar, as.data.frame(X))
					forest <- cforest(ystar ~ ., data = citdata, controls = cforest_unbiased(mtry = mtry, ntree = ntree))
					fhat <- predict(forest, OOB = TRUE, type = "response") #### pr?diction avec l'arbre
				}
				for(k in 1:nind){
					indiv <- which(id==unique(id)[k])
					V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,, drop=FALSE])+diag(as.numeric(sigmahat),length(indiv),length(indiv))
					btilde[k,] <- Btilde%*%t(Z[indiv,, drop=FALSE])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
					epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,, drop=FALSE]%*%btilde[k,]
				}

				sigm <- sigmahat
				sigmahat <- sig(sigma = sigmahat,id = id, Z = Z, epsilon = epsilonhat,Btilde = Btilde)
				Btilde  <- bay(bhat = btilde,Bhat = Btilde,Z = Z,id = id,sigmahat = sigm)
				Vrai <- c(Vrai, logV(Y,fhat,Z,time,id,Btilde,0,sigmahat,sto))
				if(i > 1) inc <- abs((Vrai[i-1]-Vrai[i])/Vrai[i-1])
				if(inc < delta) {
					print(paste0("stopped after ", i, " iterations."))
					if(!conditional){
						sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, vraisemblance = Vrai,id=id, time=time, OOB =OOB)
					}else {
						sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, vraisemblance = Vrai,id=id, time=time)
					}
					class(sortie) <- "longituRF"
					return(sortie)
				}
			}
			if(!conditional){
				sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, vraisemblance=Vrai,id=id, time=time, OOB =OOB)
			}else{
				sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, vraisemblance=Vrai,id=id, time=time)
			}
			class(sortie) <- "longituRF"
			return(sortie)
		}
	}
	for(i in 1:iter){
		ystar <- rep(0,length(Y))
		for(k in 1:nind){
			indiv <- which(id==unique(id)[k])
			ystar[indiv] <- Y[indiv]- Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
		}

		if(!conditional){
			forest <- randomForest(as.data.frame(X),ystar,mtry=mtry,ntree=ntree, importance = TRUE) ### on construit l'arbre
			fhat <- predict(forest) #### pr?diction avec l'arbre
			OOB[i] <- forest$mse[ntree]
		}else{
			citdata <- cbind(ystar, as.data.frame(X))
			forest <- cforest(ystar ~ ., data = citdata, controls = cforest_unbiased(mtry = mtry, ntree = ntree))
			fhat <- predict(forest, OOB = TRUE, type = "response") #### pr?diction avec l'arbre
		}
		for(k in 1:nind){
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
		if(i > 1) inc <- abs((Vrai[i-1]-Vrai[i])/Vrai[i-1])
		if(inc < delta) {
			print(paste0("stopped after ", i, " iterations."))
			if(!conditional){
				sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), omega=omega, sigma_sto =sigma2, time = time, sto= sto,Vraisemblance=Vrai,id=id, OOB =OOB)
			} else{
				sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), omega=omega, sigma_sto =sigma2, time = time, sto= sto,Vraisemblance=Vrai,id=id)
			}
			class(sortie) <- "longituRF"
			return(sortie)
		}
	}
	if(!conditional){
		sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id),omega=omega, sigma_sto =sigma2, time = time, sto= sto,Vraisemblance=Vrai,id=id, OOB =OOB)
	}else{
		sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id),omega=omega, sigma_sto =sigma2, time = time, sto= sto,Vraisemblance=Vrai,id=id)
	}
	class(sortie) <- "longituRF"
	return(sortie)
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


#' predict with longitudinal trees and random forests.
#'
#' @param object : a \code{longituRF} output of (S)MERF; (S)REEMforest; (S)MERT or (S)REEMtree function.
#' @param X [matrix]: matrix of the fixed effects for the new observations to be predicted.
#' @param Z [matrix]: matrix of the random effects for the new observations to be predicted.
#' @param id [vector]: vector of the identifiers of the new observations to be predicted.
#' @param time [vector]: vector of the time measurements of the new observations to be predicted.
#' @param ... : low levels arguments.
#'
#' @import stats
#' @import randomForest
#' @import party
#'
#' @return vector of the predicted output for the new observations.
#'
#' @export
#'
#' @examples \donttest{
#' set.seed(123)
#' data <- DataLongGenerator(n=20) # Generate the data composed by n=20 individuals.
#' REEMF <- REEMforest(X=data$X,Y=data$Y,Z=data$Z,id=data$id,time=data$time,mtry=2,ntree=500,sto="BM")
#' # Then we predict on the learning sample :
#' pred.REEMF <- predict(REEMF, X=data$X,Z=data$Z,id=data$id, time=data$time)
#' # Let's have a look at the predictions
#' # the predictions are in red while the real output trajectories are in blue:
#' par(mfrow=c(4,5),mar=c(2,2,2,2))
#' for (i in unique(data$id)){
#'   w <- which(data$id==i)
#'   plot(data$time[w],data$Y[w],type="l",col="blue")
#'   lines(data$time[w],pred.REEMF[w], col="red")
#' }
#' # Train error :
#' mean((pred.REEMF-data$Y)^2)
#'
#' # The same function can be used with a fitted SMERF model:
#' smerf <-MERF(X=data$X,Y=data$Y,Z=data$Z,id=data$id,time=data$time,mtry=2,ntree=500,sto="BM")
#' pred.smerf <- predict(smerf, X=data$X,Z=data$Z,id=data$id, time=data$time)
#' # Train error :
#' mean((pred.smerf-data$Y)^2)
#' # This function can be used even on a MERF model (when no stochastic process is specified)
#' merf <-MERF(X=data$X,Y=data$Y,Z=data$Z,id=data$id,time=data$time,mtry=2,ntree=500,sto="none")
#' pred.merf <- predict(merf, X=data$X,Z=data$Z,id=data$id, time=data$time)
#' # Train error :
#' mean((pred.merf-data$Y)^2)
#'
#' }
predict.longituRF <- function(object, X,Z,id,time,...){
	n <- length(unique(id))
	id_btilde <- object$id_btilde
	if(class(object$forest) == "BinaryTree"){
		f <- predict(object$forest, newdata = as.data.frame(X))
		if(exists("beta", object)){
			leaf <- unique(f)
			nnodes <- length(object$beta)
			if(length(leaf) != nnodes){
				stop("Distinct predicted values are not in 1:1 correspondence with tree leaves")
			}
			aux <- f
			for(p in seq_len(nnodes)){
				f[which(aux == leaf[p])] <- object$beta[p]
			}
			rm(aux)
		}
	} else {
		f <- predict(object$forest, newdata = as.data.frame(X))
	}

	Time <- object$time
	id_btilde <- object$id_btilde
	Ypred <- rep(0,length(id))
	id.app=object$id
	if (object$sto=="none"){
		for (i in 1:length(unique(id))){
			w <- which(id==unique(id)[i])
			k <- which(id_btilde==unique(id)[i])
			Ypred[w] <- f[w] + Z[w,, drop=FALSE]%*%object$random_effects[k,]
		}
		return(Ypred)
	}

	if (object$sto=="exp"){
		for (i in 1:length(unique(id))){
			w <- which(id==unique(id)[i])
			k <- which(id_btilde==unique(id)[i])
			om <- which(id.app==unique(id)[i])
			Ypred[w] <- f[w] + Z[w,, drop=FALSE]%*%object$random_effects[k,] + predict.exp(object$omega[om],Time[om],time[w], object$alpha)
		}
		return(Ypred)
	}

	if (object$sto=="fbm"){
		for (i in 1:length(unique(id))){
			w <- which(id==unique(id)[i])
			k <- which(id_btilde==unique(id)[i])
			om <- which(id.app==unique(id)[i])
			Ypred[w] <- f[w] + Z[w,, drop=FALSE]%*%object$random_effects[k,] + predict.fbm(object$omega[om],Time[om],time[w], object$Hurst)
		}
		return(Ypred)
	}

	for (i in 1:length(unique(id))){
		w <- which(id==unique(id)[i])
		k <- which(id_btilde==unique(id)[i])
		om <- which(id.app==unique(id)[i])
		Ypred[w] <- f[w] + Z[w,, drop=FALSE]%*%object$random_effects[k,] + predict.sto(object$omega[om],Time[om],time[w], object$sto)
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

	if(class(sto)=="function"){
		for (i in 1:length(time.test)){
			inf <- which(time.app<=time.test[i])
			sup <- which(time.app>time.test[i])
			if(length(inf)>0){
				if(length(sup)>0){
					time_inf <- max(time.app[inf])
					time_sup <- min(time.app[sup])
					pred[i] <- mean(c(omega[which(time.app==time_inf)],omega[which(time.app==time_sup)]))}
				if(length(sup)==0) {time_inf <- max(time.app[inf])
				pred[i] <- omega[which(time.app==time_inf)]*((sto(time.test[i],max(time.app)))/sto(max(time.app),max(time.app)))
				}
			}
			if(length(sup)>0 & length(inf)==0){
				time_sup <- min(time.app[sup])
				pred[i] <- omega[which(time.app==time_sup)]*((sto(time.test[i],min(time.app)))/sto(min(time.app),min(time.app)))
			}
		}
		return(pred)
	}else{
		for(i in 1:length(time.test)){
			inf <- which(time.app<=time.test[i])
			sup <- which(time.app>time.test[i])
			if(length(inf)>0){
				if(length(sup)>0){
					time_inf <- max(time.app[inf])
					time_sup <- min(time.app[sup])
					pred[i] <- mean(c(omega[which(time.app==time_inf)],omega[which(time.app==time_sup)]))}
				if(length(sup)==0) {
					time_inf <- max(time.app[inf])
					if(sto=="BM"){
						pred[i] <- omega[which(time.app==time_inf)]}
					if(sto=="OrnUhl"){
						pred[i] <- omega[which(time.app==time_inf)]*(exp(-abs(time.test[i]-max(time.app))/2))
					}
					if(sto=="BBridge"){
						pred[i] <- omega[which(time.app==time_inf)]*((1-time.test[i])/(1-max(time.app)^2))
					}
				}
			}
			if(length(sup)>0 & length(inf)==0){
				time_sup <- min(time.app[sup])
				if(sto=="BM"){
					pred[i] <- omega[which(time.app==time_sup)]*(time.test[i]/min(time.app))}
				if(sto=="OrnUhl"){
					pred[i] <- omega[which(time.app==time_sup)]*(exp(-abs(time.test[i]-min(time.app))/2))
				}
				if(sto=="BBridge"){
					pred[i] <- omega[which(time.app==time_sup)]*(time.test[i]/min(time.app))
				}
			}
		}
	}
	return(pred)
}

#' Title
#'
#' @import stats
#'
#' @keywords internal
sto_analysis <- function(sto, time){
	MAT <- matrix(0,length(time), length(time))

	if(class(sto)=="function"){
		for (i in 1:length(time)){
			for (j in 1:length(time)){
				MAT[i,j] <- sto(time[i], time[j])
			}
		}
		return(MAT)
	}

	if(sto=="BM"){
		for (i in 1:length(time)){
			for (j in 1:length(time)){
				MAT[i,j] <- min(time[i], time[j])
			}
		}
		return(MAT)
	}

	if(sto=="OrnUhl"){
		for (i in 1:length(time)){
			for (j in 1:length(time)){
				MAT[i,j] <- exp(-abs(time[i]-time[j])/2)
			}
		}
		return(MAT)
	}

	if(sto=="BBridge"){
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
	for(j in 1:nind){
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
	for(j in 1:nind){
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
	for(i in 1:nind){
		w <- which(id==unique(id)[i])
		V <- Z[w,, drop=FALSE]%*%Btilde%*%t(Z[w,, drop=FALSE])+diag(as.numeric(sigmahat),length(w),length(w))
		S1 <- S1 + t(Phi[w,, drop=FALSE])%*%solve(V)%*%Phi[w,, drop=FALSE] # square matrix of dimension nnodes
		S2 <- S2 + t(Phi[w,, drop = FALSE])%*%solve(V)%*%Y[w] # nnodes rows x 1 colum # nnodes rows x 1 columnn
	}
	return(solve(S1)%*%S2)
}


#' (S)REEMforest algorithm
#'
#'
#'
#' (S)REEMforest is an adaptation of the random forest regression method to longitudinal data introduced by Capitaine et. al. (2020) <doi:10.1177/0962280220946080>.
#' The algorithm will estimate the parameters of the following semi-parametric stochastic mixed-effects model: \deqn{Y_i(t)=f(X_i(t))+Z_i(t)\beta_i + \omega_i(t)+\epsilon_i}
#' with \eqn{Y_i(t)} the output at time \eqn{t} for the \eqn{i}th individual; \eqn{X_i(t)} the input predictors (fixed effects) at time \eqn{t} for the \eqn{i}th individual;
#' \eqn{Z_i(t)} are the random effects at time \eqn{t} for the \eqn{i}th individual; \eqn{\omega_i(t)} is the stochastic process at time \eqn{t} for the \eqn{i}th individual
#'  which model the serial correlations of the output measurements; \eqn{\epsilon_i} is the residual error.
#'
#'
#' @param X [matrix]: A \code{N}x\code{p} matrix containing the \code{p} predictors of the fixed effects, column codes for a predictor.
#' @param Y [vector]: A vector containing the output trajectories.
#' @param id [vector]: Is the vector of the identifiers for the different trajectories.
#' @param Z [matrix]: A \code{N}x\code{q} matrix containing the \code{q} predictor of the random effects.
#' @param iter [numeric]: Maximal number of iterations of the algorithm. The default is set to \code{iter=100}
#' @param mtry [numeric]: Number of variables randomly sampled as candidates at each split. The default value is \code{p/3}.
#' @param ntree [numeric]: Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times. The default value is \code{ntree=500}.
#' @param time [time]: Is the vector of the measurement times associated with the trajectories in \code{Y},\code{Z} and \code{X}.
#' @param sto [character]: Defines the covariance function of the stochastic process, can be either \code{"none"} for no stochastic process, \code{"BM"} for Brownian motion, \code{OrnUhl} for standard Ornstein-Uhlenbeck process, \code{BBridge} for Brownian Bridge, \code{fbm} for Fractional Brownian motion; can also be a function defined by the user.
#' @param delta [numeric]: The algorithm stops when the difference in log likelihood between two iterations is smaller than \code{delta}. The default value is set to O.O01
#' @param conditional [logical]: Determines if the random forest algorithm is the one implemented by Breiman (2001) and available in \code{\link[randomForest]{randomForest}}, which is the default \code{FALSE}, or if it is the implemented by Strobl et al. (2007) using conditional inference trees in \code{\link[party]{cforest}}.
#'
#' @import stats
#' @import randomForest
#' @import party
#'
#' @return A fitted (S)REEMforest model which is a list of the following elements: \itemize{
#' \item \code{forest:} Random forest obtained at the last iteration.
#' \item \code{random_effects :} predictions of random effects for different trajectories.
#' \item \code{id_btilde:} Identifiers of individuals associated with the predictions \code{random_effects}.
#' \item \code{var_random_effects: } Estimation of the variance covariance matrix of random effects.
#' \item \code{sigma_sto: } Estimation of the volatility parameter of the stochastic process.
#' \item \code{sigma: } Estimation of the residual variance parameter.
#' \item \code{time: } The vector of the measurement times associated with the trajectories in \code{Y},\code{Z} and \code{X}.
#' \item \code{sto: } Stochastic process used in the model.
#' \item \code{Vraisemblance:} Log-likelihood of the different iterations.
#' \item \code{id: } Vector of the identifiers for the different trajectories.
#' \item \code{OOB: } OOB error of the fitted random forest at each iteration, only when \code{conditional=FALSE}.
#' }
#' @export
#'
#' @references
#'
#' Louis Capitaine, Robin Genuer, and Rodolphe Thiébaut (2020). Random forests for high-dimensional longitudinal data. Statistical Methods in Medical Research, 096228022094608. \doi{10.1177/0962280220946080}
#'
#' Leo Breiman (2001). Random Forests. Machine Learning, 45(1), 5–32.
#'
#' Carolin Strobl, Anne-Laure Boulesteix, Achim Zeileis and Torsten Hothorn (2007). Bias in Random Forest Variable Importance Measures: Illustrations, Sources and a Solution. BMC Bioinformatics, 8(25).
#'
#'
#' @examples \donttest{
#' set.seed(123)
#' data <- DataLongGenerator(n=20) # Generate the data composed by n=20 individuals.
#' # Train a SREEMforest model on the generated data. Should take ~ 50 secondes
#' # The data are generated with a Brownian motion
#' #  so we use the parameter sto="BM" to specify a Brownian motion as stochastic process
#' SREEMF <- REEMforest(X=data$X,Y=data$Y,Z=data$Z,id=data$id,time=data$time,mtry=2,ntree=500,sto="BM")
#' SREEMF$forest # is the fitted random forest (obtained at the last iteration).
#' SREEMF$random_effects # are the predicted random effects for each individual.
#' SREEMF$omega # are the predicted stochastic processes.
#' plot(SREEMF$Vraisemblance) #evolution of the log-likelihood.
#' SREEMF$OOB # OOB error at each iteration.
#' cSREEMF <- REEMforest(X=data$X,Y=data$Y,Z=data$Z,id=data$id,time=data$time,mtry=2,ntree=500,sto="BM",conditional=TRUE)
#' cSREEMF$forest # is the fitted random forest (obtained at the last iteration).
#' cSREEMF$random_effects # are the predicted random effects for each individual.
#' cSREEMF$omega # are the predicted stochastic processes.
#' plot(cSREEMF$Vraisemblance) #evolution of the log-likelihood.
#' }
#'
REEMforest <- function(X,Y,id,Z,iter=100,mtry=ceiling(ncol(X)/3),ntree=500, time, sto, delta = 0.001, conditional = FALSE){
	q <- dim(Z)[2]
	nind <- length(unique(id))
	btilde <- matrix(0,nind,q) #### Pour la ligne i, on a les effets aléatoires de l'individu i
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

	if(!is.logical(conditional) || length(conditional) != 1 || is.na(conditional)){
		conditional <- FALSE
	}

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
			h <- opti.FBMreem(X,Y,id,Z,iter, mtry,ntree,time, conditional)
			for (i in 1:iter){
				ystar <- rep(0,length(Y))
				for (k in 1:nind){ #### on retrace les effets al?atoires
					indiv <- which(id==unique(id)[k])
					ystar[indiv] <- Y[indiv]- Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
				}

				if(!conditional){
					forest <- randomForest(as.data.frame(X), ystar,mtry=mtry,ntree=ntree, importance = TRUE, keep.inbag=TRUE)
					f1 <- predict(forest, as.data.frame(X),nodes=TRUE)
					# f1 == apply(predict(forest, as.data.frame(X), predict.all = TRUE)$individual, 1, mean)
					trees <- attr(f1, "nodes")
					inbag <- forest$inbag
					OOB[i] <- forest$mse[ntree]
				}else if(conditional){
					citdata <- cbind(ystar, as.data.frame(X))
					forest <- cforest(ystar ~ ., data = citdata, controls = cforest_unbiased(mtry = mtry, ntree = ntree))
					f1 <- predict(forest, OOB = TRUE, type = "response")
					# ?cforest:
					# The aggregation scheme works by averaging observation weights extracted from each of the ntree trees and NOT by averaging predictions directly as in randomForest.
					# See Hothorn et al. (2004) for a description.
					trees <- as.matrix(as.data.frame(forest@where,
									 row.names = row.names(citdata),
									 col.names = paste0("V", seq_len(ntree))))
					inbag <- as.matrix(as.data.frame(forest@weights,
									 row.names = row.names(citdata),
									 col.names = paste0("V", seq_len(ntree))))
				}

				matrice.pred <- matrix(NA,length(Y),ntree)
				# this matrix keeps the predictions for each observation in Y through all the trees where such observation is IN-BAG

				beta <- vector(mode = "list", length = ntree)
				for(k in 1:ntree){
					indii <- unique(trees[,k])
					nnodes <- length(indii)
					Phi <- matrix(0,length(Y), nnodes)

					for(l in seq_len(nnodes)){
						w <- which(trees[,k]==indii[l])
						Phi[w,l] <- 1
					}
					oobags <- unique(which(inbag[,k]==0))
					beta[[k]] <- Moy_fbm(id[-oobags],Btilde,sigmahat,Phi[-oobags,],ystar[-oobags],Z[-oobags,,drop=FALSE],h,time[-oobags], sigma2)
					matrice.pred[oobags,k] <- Phi[oobags,]%*%beta[[k]]
				}

				fhat <- rep(NA,length(Y))
				for(k in 1:length(Y)){
					w <- which(is.na(matrice.pred[k,]))
					fhat[k] <- mean(matrice.pred[k,-w])
					# this is the average of predictions of observation k through those trees which have k IN-BAG
					# making fhat <- predict(forest, as.data.frame(X)) after modifying leaves with beta's
					# would be for each observation k the average of the predictions through ALL the trees
				}

				for(k in 1:nind){ ### calcul des effets al?atoires par individu
					indiv <- which(id==unique(id)[k])
					K <- cov.fbm(time[indiv],h)
					V <- Z[indiv,, drop=FALSE] %*% Btilde %*% t(Z[indiv,, drop=FALSE]) + diag(as.numeric(sigmahat),length(indiv),length(indiv)) + sigma2 * K
					btilde[k,] <- Btilde %*% t(Z[indiv,, drop=FALSE]) %*% solve(V) %*% (Y[indiv] - fhat[indiv])
					omega[k,which(id_omega[k,]==1)] <- sigma2 * K %*% solve(V) %*% (Y[indiv] - fhat[indiv])
					omega2[indiv] <- omega[k,which(id_omega[k,]==1)]
					epsilonhat[indiv] <- Y[indiv] - fhat[indiv] - Z[indiv,, drop=FALSE] %*% btilde[k,] - omega[indiv]
				}
				sigm <- sigmahat
				B <- Btilde
				sigmahat <- sig.fbm(Y,sigmahat,id, Z, epsilonhat, Btilde, time, sigma2,h) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
				Btilde  <- bay.fbm(btilde,Btilde,Z,id,sigm, time, sigma2,h) #### MAJ des param?tres de la variance des effets al?atoires.
				### MAJ de la volatilit? du processus stochastique
				sigma2 <- gam_fbm(Y,sigm,id,Z,B,time,sigma2,omega,id_omega,h)
				Vrai <- c(Vrai, logV.fbm(Y,fhat,Z[,,drop=FALSE],time,id,Btilde,sigma2,sigmahat,h))
				if(i > 1) inc <- abs(Vrai[i-1]-Vrai[i]) / abs(Vrai[i-1])
				if(inc < delta){
					print(paste0("stopped after ", i, " iterations."))
					if(!conditional){
						for(k in 1:ntree){
							indii <- unique(trees[,k])
							forest$forest$nodepred[indii,k] <- beta[[k]]
						}
						sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, vraisemblance = Vrai,id=id, time =time, Hurst=h, OOB =OOB, omega=omega2)
					}else{
						sortie <- list(forest=forest,beta=beta,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, vraisemblance = Vrai,id=id, time =time, Hurst=h, omega=omega2)
					}
					class(sortie) <- "longituRF"
					return(sortie)
				}
			}
			if(!conditional){
				for(k in 1:ntree){
					indii <- unique(trees[,k])
					forest$forest$nodepred[indii,k] <- beta[[k]]
				}
				sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto=sto,omega=omega2, sigma_sto =sigma2, time =time, sto= sto, Hurst =h, Vraisemblance=Vrai, OOB =OOB)
			}else{
				sortie <- list(forest=forest,beta=beta,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto=sto,omega=omega2, sigma_sto =sigma2, time =time, sto= sto, Hurst =h, Vraisemblance=Vrai)
			}
			class(sortie) <- "longituRF"
			return(sortie)
		}

		if(sto=="none"){
			for(i in 1:iter){
				ystar <- rep(0,length(Y))
				for(k in 1:nind){ #### on retrace les effets al?atoires
					indiv <- which(id==unique(id)[k])
					ystar[indiv] <- Y[indiv]- Z[indiv,, drop=FALSE]%*%btilde[k,]
				}
				if(!conditional){
					forest <- randomForest(as.data.frame(X), ystar,mtry=mtry,ntree=ntree, importance = TRUE, keep.inbag=TRUE)
					f1 <- predict(forest,as.data.frame(X),nodes=TRUE)
					trees <- attr(f1, "nodes")
					inbag <- forest$inbag
					OOB[i] <- forest$mse[ntree]
				}else if(conditional){
					citdata <- cbind(ystar, as.data.frame(X))
					forest <- cforest(ystar ~ ., data = citdata, controls = cforest_unbiased(mtry = mtry, ntree = ntree))
					f1 <- predict(forest, OOB = TRUE, type = "response")
					trees <- as.matrix(as.data.frame(forest@where,
									 row.names = row.names(citdata),
									 col.names = paste0("V", seq_len(ntree))))
					inbag <- as.matrix(as.data.frame(forest@weights,
									 row.names = row.names(citdata),
									 col.names = paste0("V", seq_len(ntree))))
				}
				matrice.pred <- matrix(NA,length(Y),ntree)


				beta <- vector(mode = "list", length = ntree)
				for(k in 1:ntree){
					indii <- unique(trees[,k])
					nnodes <- length(indii)
					Phi <- matrix(0,length(Y),nnodes)
					for(l in seq_len(nnodes)){
						w <- which(trees[,k]==indii[l])
						Phi[w,l] <- 1
					}
					oobags <- unique(which(inbag[,k]==0))
					beta[[k]] <- Moy(id[-oobags],Btilde,sigmahat,Phi[-oobags,],ystar[-oobags],Z[-oobags,,drop=FALSE])
					matrice.pred[oobags,k] <- Phi[oobags,]%*%beta[[k]]
				}

				fhat <- rep(NA,length(Y))
				for(k in 1:length(Y)){
					w <- which(is.na(matrice.pred[k,]))
					fhat[k] <- mean(matrice.pred[k,-w])
				}

				for(k in 1:nind){ ### calcul des effets al?atoires par individu
					indiv <- which(id==unique(id)[k])
					V <- Z[indiv,, drop=FALSE] %*% Btilde %*% t(Z[indiv,, drop=FALSE]) + diag(as.numeric(sigmahat),length(indiv),length(indiv))
					btilde[k,] <- Btilde %*% t(Z[indiv,, drop=FALSE]) %*% solve(V) %*% (Y[indiv] - fhat[indiv])
					epsilonhat[indiv] <- Y[indiv] - fhat[indiv] - Z[indiv,, drop=FALSE] %*% btilde[k,]
				}

				sigm <- sigmahat
				sigmahat <- sig(sigmahat,id, Z, epsilonhat, Btilde) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
				Btilde  <- bay(btilde,Btilde,Z,id,sigm) #### MAJ des param?tres de la variance des effets al?atoires.
				Vrai <- c(Vrai, logV(Y,fhat,Z,time,id,Btilde,0,sigmahat,sto))
				if(i > 1) inc <- abs((Vrai[i-1]-Vrai[i])/Vrai[i-1])
				if(inc < delta) {
					print(paste0("stopped after ", i, " iterations."))
					if(!conditional){
						for(k in seq_len(ntree)){
							indii <- unique(trees[,k])
							forest$forest$nodepred[indii,k] <- beta[[k]]
						}
						sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, vraisemblance = Vrai,id=id, time =time, OOB =OOB)
					}else{
						sortie <- list(forest=forest,beta=beta,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, vraisemblance = Vrai,id=id, time =time)
					}
					class(sortie) <- "longituRF"
					return(sortie)
				}
			}
			if(!conditional){
				for(k in seq_len(ntree)){
					indii <- unique(trees[,k])
					forest$forest$nodepred[indii,k] <- beta[[k]]
				}
				sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, id = id , time = time , Vraisemblance=Vrai, OOB =OOB)
			}else{
				sortie <- list(forest=forest,beta=beta,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, id = id , time = time , Vraisemblance=Vrai)
			}
			class(sortie) <- "longituRF"
			return(sortie)
		}
	}
	for (i in 1:iter){

		ystar <- rep(0,length(Y))
		for (k in 1:nind){ #### on retrace les effets al?atoires
			indiv <- which(id==unique(id)[k])
			ystar[indiv] <- Y[indiv]- Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
		}

		if(!conditional){
			forest <- randomForest(as.data.frame(X), ystar,mtry=mtry,ntree=ntree, importance = TRUE, keep.inbag=TRUE)
			f1 <- predict(forest,as.data.frame(X),nodes=TRUE)
			trees <- attr(f1, "nodes")
			inbag <- forest$inbag
			OOB[i] <- forest$mse[ntree]
		}else{
					citdata <- cbind(ystar, as.data.frame(X))
					forest <- cforest(ystar ~ ., data = citdata, controls = cforest_unbiased(mtry = mtry, ntree = ntree))
					f1 <- predict(forest, OOB = TRUE, type = "response")
					trees <- as.matrix(as.data.frame(forest@where,
									 row.names = row.names(citdata),
									 col.names = paste0("V", seq_len(ntree))))
					inbag <- as.matrix(as.data.frame(forest@weights,
									 row.names = row.names(citdata),
									 col.names = paste0("V", seq_len(ntree))))
		}
		matrice.pred <- matrix(NA,length(Y),ntree)

		beta <- vector(mode = "list", length = ntree)
		for(k in 1:ntree){
			indii <- unique(trees[,k])
			nnodes <- length(indii)
			Phi <- matrix(0,length(Y), nnodes)
			for(l in seq_len(nnodes)){
				w <- which(trees[,k]==indii[l])
				Phi[w,l] <- 1
			}
			oobags <- unique(which(inbag[,k]==0))
			beta[[k]] <- Moy_sto(id[-oobags],Btilde,sigmahat,Phi[-oobags,, drop=FALSE],ystar[-oobags],Z[-oobags,,drop=FALSE], sto, time[-oobags], sigma2)
			matrice.pred[oobags,k] <- Phi[oobags,]%*%beta[[k]]
		}

		fhat <- rep(NA,length(Y))
		for(k in 1:length(Y)){
			w <- which(is.na(matrice.pred[k,])==TRUE)
			fhat[k] <- mean(matrice.pred[k,-w])
		}

		for(k in 1:nind){ ### calcul des effets al?atoires par individu
			indiv <- which(id==unique(id)[k])
			K <- sto_analysis(sto,time[indiv])
			V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,, drop=FALSE])+diag(as.numeric(sigmahat),length(indiv),length(indiv))+ sigma2*K
			btilde[k,] <- Btilde %*% t(Z[indiv,, drop=FALSE]) %*% solve(V) %*% (Y[indiv] - fhat[indiv])
			omega[indiv] <- sigma2 * K %*% solve(V) %*% (Y[indiv]-fhat[indiv])
			epsilonhat[indiv] <- Y[indiv] - fhat[indiv] - Z[indiv,, drop=FALSE] %*% btilde[k,] - omega[indiv]
		}
		sigm <- sigmahat
		B <- Btilde
		sigmahat <- sig_sto(sigmahat,id, Z, epsilonhat, Btilde, time, sigma2,sto) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
		Btilde  <- bay_sto(btilde,Btilde,Z,id,sigm, time, sigma2,sto) #### MAJ des param?tres de la variance des effets al?atoires.
		### MAJ de la volatilit? du processus stochastique
		sigma2 <- gam_sto(sigm,id,Z,B,time,sigma2,sto,omega)
		Vrai <- c(Vrai, logV(Y,fhat,Z[,,drop=FALSE],time,id,Btilde,sigma2,sigmahat,sto))
		if(i > 1) {
			inc <- abs((Vrai[i-1]-Vrai[i]) / Vrai[i-1])
			if(Vrai[i] < Vrai[i-1]) {
				if(!conditional){
					for(k in seq_len(ntree)){
						indii <- unique(trees[,k])
						forest$forest$nodepred[indii,k] <- beta[[k]]
					}
					reemfouille <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), omega=omega, sigma_sto =sigma2, time = time, sto= sto,Vraisemblance=Vrai,id=id, OOB =OOB)
				}else{
					reemfouille <- list(forest=forest,beta=beta,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), omega=omega, sigma_sto =sigma2, time = time, sto= sto,Vraisemblance=Vrai,id=id)
				}
			}
		}
		if(inc < delta) {
			print(paste0("stopped after ", i, " iterations."))
			class(reemfouille) <- "longituRF"
			return(reemfouille)
		}
	}
	if(!conditional){
		for(k in seq_len(ntree)){
			indii <- unique(trees[,k])
			forest$forest$nodepred[indii,k] <- beta[[k]]
		}
		sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id),omega=omega, sigma_sto =sigma2, time = time, sto= sto, id=id, OOB =OOB, Vraisemblance=Vrai)
	}else{
		sortie <- list(forest=forest,beta=beta,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id),omega=omega, sigma_sto =sigma2, time = time, sto= sto, id=id, OOB =OOB, Vraisemblance=Vrai)
	}
	class(sortie) <- "longituRF"
	return(sortie)
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
#'
#' (S)MERT is an adaptation of the regression trees method to longitudinal data introduced by Hajjem et. al. (2011) <doi:10.1016/j.spl.2010.12.003>.
#' The model has been improved by Capitaine et. al. (2020) <doi:10.1177/0962280220946080> with the addition of a stochastic process.
#' The algorithm will estimate the parameters of the following semi-parametric stochastic mixed-effects model: \deqn{Y_i(t)=f(X_i(t))+Z_i(t)\beta_i + \omega_i(t)+\epsilon_i}
#' with \eqn{Y_i(t)} the output at time \eqn{t} for the \eqn{i}th individual; \eqn{X_i(t)} the input predictors (fixed effects) at time \eqn{t} for the \eqn{i}th individual;
#' \eqn{Z_i(t)} are the random effects at time \eqn{t} for the \eqn{i}th individual; \eqn{\omega_i(t)} is the stochastic process at time \eqn{t} for the \eqn{i}th individual
#'  which model the serial correlations of the output measurements; \eqn{\epsilon_i} is the residual error.
#'
#' @param X [matrix]: A \code{N}x\code{p} matrix containing the \code{p} predictors of the fixed effects, column codes for a predictor.
#' @param Y [vector]: A vector containing the output trajectories.
#' @param id [vector]: Is the vector of the identifiers for the different trajectories.
#' @param Z [matrix]: A \code{N}x\code{q} matrix containing the \code{q} predictor of the random effects.
#' @param iter [numeric]: Maximal number of iterations of the algorithm. The default is set to \code{iter=100}
#' @param time [vector]: Is the vector of the measurement times associated with the trajectories in \code{Y},\code{Z} and \code{X}.
#' @param sto [character]: Defines the covariance function of the stochastic process, can be either \code{"none"} for no stochastic process, \code{"BM"} for Brownian motion, \code{OrnUhl} for standard Ornstein-Uhlenbeck process, \code{BBridge} for Brownian Bridge, \code{fbm} for Fractional Brownian motion; can also be a function defined by the user.
#' @param delta [numeric]: The algorithm stops when the difference in log likelihood between two iterations is smaller than \code{delta}. The default value is set to O.O01
#' @param conditional [logical]: Determines if the (S)MERT algorithm uses conditional inference trees as implemented by Hothorn et al. (2006) in \code{\link[party]{ctree}} (\code{TRUE}) or the usual regression trees as implemented in \code{[rpart]{rpart}}, being this case (\code{TRUE}) the default.
#'
#'
#'
#' @import stats
#' @import rpart
#' @import party
#'
#'
#'
#' @return A fitted (S)MERT model which is a list of the following elements: \itemize{
#' \item \code{forest:} Tree obtained at the last iteration.
#' \item \code{random_effects :} predictions of random effects for different trajectories.
#' \item \code{id_btilde:} Identifiers of individuals associated with the predictions \code{random_effects}.
#' \item \code{var_random_effects: } Estimation of the variance covariance matrix of random effects.
#' \item \code{sigma_sto: } Estimation of the volatility parameter of the stochastic process.
#' \item \code{sigma: } Estimation of the residual variance parameter.
#' \item \code{time: } The vector of the measurement times associated with the trajectories in \code{Y},\code{Z} and \code{X}.
#' \item \code{sto: } Stochastic process used in the model.
#' \item \code{Vraisemblance:} Log-likelihood of the different iterations.
#' \item \code{id: } Vector of the identifiers for the different trajectories.
#' }
#'
#' @export
#'
#'
#' @references
#' Ahlem Hajjem, François Bellavance, and Denis Larocque (2011). Mixed effects regression trees for clustered data. Statistics & Probability Letters, 81(4), 451–459. \doi{10.1016/j.spl.2010.12.003}
#'
#' Louis Capitaine, Robin Genuer, and Rodolphe Thiébaut (2020). Random forests for high-dimensional longitudinal data. Statistical Methods in Medical Research, 096228022094608. \doi{10.1177/0962280220946080}
#'
#' Breiman L., Friedman J. H., Olshen R. A., and Stone, C. J. (1984) Classification and Regression Trees. Wadsworth.
#'
#' Torsten Hothorn, Kurt Hornik and Achim Zeileis (2006). Unbiased Recursive Partitioning: A Conditional Inference Framework. Journal of Computational and Graphical Statistics, 15(3), 651--674.
#'
#'
#' @examples
#' set.seed(123)
#' data <- DataLongGenerator(n=20) # Generate the data composed by n=20 individuals.
#' # Train a SMERF model on the generated data. Should take ~ 50 secondes
#' # The data are generated with a Brownian motion,
#' # so we use the parameter sto="BM" to specify a Brownian motion as stochastic process
#' smert <- MERT(X=data$X,Y=data$Y,Z=data$Z,id=data$id,time=data$time,sto="BM")
#' smert$forest # is the fitted regression tree (obtained at the last iteration).
#' smert$random_effects # are the predicted random effects for each individual.
#' smert$omega # are the predicted stochastic processes.
#' plot(smert$Vraisemblance) #evolution of the log-likelihood.
#' csmert <- MERT(X=data$X,Y=data$Y,Z=data$Z,id=data$id,time=data$time,sto="BM",conditional=TRUE)
#' csmert$forest # is the fitted regression tree (obtained at the last iteration).
#' csmert$random_effects # are the predicted random effects for each individual.
#' csmert$omega # are the predicted stochastic processes.
#' plot(csmert$Vraisemblance) #evolution of the log-likelihood.
#'
#'
MERT <- function(X,Y,id,Z,iter=100,time, sto, delta = 0.001, conditional = FALSE){
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

	if(!is.logical(conditional) || length(conditional) != 1 || is.na(conditional)){
		conditional <- FALSE
	}

	if(class(sto)=="character"){
		if(sto=="none"){
			for(i in 1:iter){
				ystar <- rep(0,length(Y))
				for(k in 1:nind){ #### on retrace les effets al?atoires
					indiv <- which(id==unique(id)[k])
					ystar[indiv] <- Y[indiv]- Z[indiv,, drop=FALSE]%*%btilde[k,]
				}

				if(!conditional){
					tree <- rpart(ystar~.,as.data.frame(X)) ### on construit l'arbre
				}else{
					citdata <- cbind(ystar, as.data.frame(X))
					tree <- ctree(ystar~., citdata)
				}
				fhat <- predict(tree, as.data.frame(X)) #### pr?diction avec l'arbre
				for(k in 1:nind){ ### calcul des effets al?atoires par individu
					indiv <- which(id==unique(id)[k])
					V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,, drop=FALSE])+diag(as.numeric(sigmahat),length(indiv),length(indiv))
					btilde[k,] <- Btilde%*%t(Z[indiv,, drop=FALSE])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
					epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,, drop=FALSE]%*%btilde[k,]
				}
				sigm <- sigmahat
				sigmahat <- sig(sigmahat,id, Z, epsilonhat, Btilde) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
				Btilde  <- bay(btilde,Btilde,Z,id,sigm) #### MAJ des param?tres de la variance des effets al?atoires.
				Vrai <- c(Vrai, logV(Y,fhat,Z[,,drop=FALSE],time,id,Btilde,0,sigmahat,"none"))
				if(i > 1) inc <- (Vrai[i-1]-Vrai[i])/Vrai[i-1]
				if(inc< delta) {
					print(paste0("stopped after ", i, " iterations."))
					sortie <- list(forest=tree,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, Vraisemblance = Vrai, id =id, time=time)
					class(sortie) <- "longituRF"
					return(sortie)
				}
			}
			sortie <- list(forest=tree,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, Vraisemblance=Vrai, id=id, time=time)
			class(sortie) <- "longituRF"
			return(sortie)
		}
	}
	for(i in 1:iter){
		ystar <- rep(0,length(Y))
		for(k in 1:nind){ #### on retrace les effets al?atoires
			indiv <- which(id==unique(id)[k])
			ystar[indiv] <- Y[indiv]- Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
		}

		if(!conditional){
			tree <- rpart(ystar~.,as.data.frame(X)) ### on construit l'arbre
		}else{
			citdata <- cbind(ystar, as.data.frame(X))
			tree <- ctree(ystar~., citdata)
		}
		fhat <- predict(tree, as.data.frame(X)) #### pr?diction avec l'arbre
		for(k in 1:nind){ ### calcul des effets al?atoires par individu
			indiv <- which(id==unique(id)[k])
			K <- sto_analysis(sto,time[indiv])
			V <- Z[indiv,, drop=FALSE] %*% Btilde %*% t(Z[indiv,, drop=FALSE]) + diag(as.numeric(sigmahat),length(indiv),length(indiv)) + sigma2 * K
			btilde[k,] <- Btilde %*% t(Z[indiv,, drop=FALSE]) %*% solve(V) %*% (Y[indiv] - fhat[indiv])
			omega[indiv] <- sigma2 * K %*% solve(V) %*% (Y[indiv] - fhat[indiv] - Z[indiv,,drop=FALSE] %*% btilde[k,])
			epsilonhat[indiv] <- Y[indiv] - fhat[indiv] - Z[indiv,, drop=FALSE] %*% btilde[k,] - omega[indiv]
		}
		#### pr?diction du processus stochastique:
		sigm <- sigmahat
		B <- Btilde
		sigmahat <- sig_sto(sigmahat,id, Z, epsilonhat, Btilde, time, sigma2,sto) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
		Btilde  <- bay_sto(btilde,Btilde,Z,id,sigm, time, sigma2,sto) #### MAJ des param?tres de la variance des effets al?atoires.
		### MAJ de la volatilit? du processus stochastique
		sigma2 <- gam_sto(sigm,id,Z,B,time,sigma2,sto,omega)
		Vrai <- c(Vrai, logV(Y,fhat,Z[,,drop=FALSE],time,id,Btilde,sigma2,sigmahat,sto))
		if(i > 1) inc <- (Vrai[i-1]-Vrai[i])/Vrai[i-1]
		if(inc < delta) {
			print(paste0("stopped after ", i, " iterations."))
			sortie <- list(forest=tree,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat,id_omega=id_omega, id_btilde=unique(id),omega=omega, sigma_sto =sigma2, time = time, sto= sto,Vraisemblance=Vrai, id = id)
			class(sortie) <- "longituRF"
			return(sortie)
		}
	}
	sortie <- list(forest=tree,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id),id_omega=id_omega,omega=omega, sigma_sto =sigma2, time = time, sto= sto, Vraisemblance=Vrai, id=id)
	class(sortie) <- "longituRF"
	return(sortie)
}

# After
# names(beta) <- unique(tree@where)
# the next function (used as feed_leaves(tree@tree))  would allow to modify the terminal leaves of the conditional inference trees

# feed_leaves <- function(ctobject, xstr = ""){
#   xstr <- paste0(xstr, sub("ctobject", "", deparse(substitute(ctobject))))
#   if(ctobject$terminal){
#     xstr <- paste0(xstr, "$prediction")
#     eval(parse(text = paste(xstr, "<<-", beta[[as.character(ctobject$nodeID)]])))
#   } else {
#     Recall(ctobject$left, xstr = xstr)
#     Recall(ctobject$right, xstr = xstr)
#   }
# }

# but it has none consequence on the prediction functions/methods predict/predict_response,
# so it only changes the not easily accessible prediction values on the terminal nodes.
# Therefore, it seems to me that changing just the values associated with terminal nodes in tree object is a bit weird.
# It is also applicable to rpart trees and to tree$frame$yval values, so I believe those trees should also be kept
# as they are computed by rpart, returning the beta values to use with the predict function as for party::ctree

#' (S)REEMtree algorithm
#'
#'
#' (S)REEMtree is an adaptation of the random-effects regression trees method to longitudinal data introduced by Sela and Simonoff. (2012) <doi:10.1007/s10994-011-5258-3>.
#' The algorithm will estimate the parameters of the following semi-parametric stochastic mixed-effects model: \deqn{Y_i(t)=f(X_i(t))+Z_i(t)\beta_i + \omega_i(t)+\epsilon_i}
#' with \eqn{Y_i(t)} the output at time \eqn{t} for the \eqn{i}th individual; \eqn{X_i(t)} the input predictors (fixed effects) at time \eqn{t} for the \eqn{i}th individual;
#' \eqn{Z_i(t)} are the random effects at time \eqn{t} for the \eqn{i}th individual; \eqn{\omega_i(t)} is the stochastic process at time \eqn{t} for the \eqn{i}th individual
#'  which model the serial correlations of the output measurements; \eqn{\epsilon_i} is the residual error.
#'
#'
#' @param X [matrix]: A \code{N}x\code{p} matrix containing the \code{p} predictors of the fixed effects, column codes for a predictor.
#' @param Y [vector]: A vector containing the output trajectories.
#' @param id [vector]: Is the vector of the identifiers for the different trajectories.
#' @param Z [matrix]: A \code{N}x\code{q} matrix containing the \code{q} predictor of the random effects.
#' @param iter [numeric]: Maximal number of iterations of the algorithm. The default is set to \code{iter=100}
#' @param time [vector]: Is the vector of the measurement times associated with the trajectories in \code{Y},\code{Z} and \code{X}.
#' @param sto [character]: Defines the covariance function of the stochastic process, can be either \code{"none"} for no stochastic process, \code{"BM"} for Brownian motion, \code{OrnUhl} for standard Ornstein-Uhlenbeck process, \code{BBridge} for Brownian Bridge, \code{fbm} for Fractional Brownian motion; can also be a function defined by the user.
#' @param delta [numeric]: The algorithm stops when the difference in log likelihood between two iterations is smaller than \code{delta}. The default value is set to O.O01
#' @param conditional [logical]: Determines if the (S)REEMtree algorithm uses conditional inference trees as implemented by Hothorn et al. (2006) in \code{\link[party]{ctree}} (\code{TRUE}) or the usual regression trees as implemented in \code{[rpart]{rpart}}, being this case (\code{TRUE}) the default.
#'
#'
#' @import stats
#' @import rpart
#' @import party
#'
#'
#'
#' @return A fitted (S)REEMtree model which is a list of the following elements: \itemize{
#' \item \code{forest:} Tree obtained at the last iteration.
#' \item \code{beta:} predicted response values at the terminal nodes of the tree, only when \code{conditional==TRUE}.
#' \item \code{random_effects :} predictions of random effects for different trajectories.
#' \item \code{id_btilde:} Identifiers of individuals associated with the predictions \code{random_effects}.
#' \item \code{var_random_effects: } Estimation of the variance covariance matrix of random effects.
#' \item \code{sigma_sto: } Estimation of the volatility parameter of the stochastic process.
#' \item \code{sigma: } Estimation of the residual variance parameter.
#' \item \code{time: } The vector of the measurement times associated with the trajectories in \code{Y},\code{Z} and \code{X}.
#' \item \code{sto: } Stochastic process used in the model.
#' \item \code{Vraisemblance:} Log-likelihood of the different iterations.
#' \item \code{id: } Vector of the identifiers for the different trajectories.
#' }
#'
#' @export
#'
#'
#' @references
#'
#' Rebecca J. Sela and Jeffrey S. Simonoff (2012). RE-EM trees: a data mining approach for longitudinal and clustered data. Machine Learning, 86(2), 169–207. \doi{10.1007/s10994-011-5258-3}
#'
#' Wei Fu and Jeffrey S. Simonoff (2015). Unbiased regression trees for longitudinal and clustered data. Computational Statistics & Data Analysis, 88, 53–74. \doi{10.1016/j.csda.2015.02.004}
#'
#' Louis Capitaine, Robin Genuer, and Rodolphe Thiébaut (2020). Random forests for high-dimensional longitudinal data. Statistical Methods in Medical Research, 096228022094608. \doi{10.1177/0962280220946080}
#'
#' Breiman L., Friedman J. H., Olshen R. A., and Stone, C. J. (1984) Classification and Regression Trees. Wadsworth.
#'
#' Torsten Hothorn, Kurt Hornik and Achim Zeileis (2006). Unbiased Recursive Partitioning: A Conditional Inference Framework. Journal of Computational and Graphical Statistics, 15(3), 651--674.
#'
#'
#' @examples
#' set.seed(123)
#' data <- DataLongGenerator(n=20) # Generate the data composed by n=20 individuals.
#' # Train a SREEMtree model on the generated data.
#' # The data are generated with a Brownian motion,
#' # so we use the parameter sto="BM" to specify a Brownian motion as stochastic process
#' sreemt <- REEMtree(X=data$X,Y=data$Y,Z=data$Z,id=data$id,time=data$time,
#' sto="BM", delta=0.0001)
#' sreemt$forest # is the fitted random forest (obtained at the last iteration).
#' sreemt$random_effects # are the predicted random effects for each individual.
#' sreemt$omega # are the predicted stochastic processes.
#' plot(sreemt$Vraisemblance) #evolution of the log-likelihood.
#' csreemt <- REEMtree(X=data$X,Y=data$Y,Z=data$Z,id=data$id,time=data$time,
#' sto="BM", delta=0.0001, conditional=TRUE)
#' csreemt$forest # is the fitted random forest (obtained at the last iteration).
#' csreemt$random_effects # are the predicted random effects for each individual.
#' csreemt$omega # are the predicted stochastic processes.
#' plot(csreemt$Vraisemblance) #evolution of the log-likelihood.
#'
REEMtree <- function(X,Y,id,Z,iter=10, time, sto, delta = 0.001, conditional = FALSE){
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

	if(!is.logical(conditional) || length(conditional) != 1 || is.na(conditional)){
		conditional <- FALSE
	}

	if(class(sto)=="character"){
		if(sto=="none"){
			for(i in 1:iter){
				ystar <- rep(0,length(Y))
				for(k in 1:nind){ #### on retrace les effets al?atoires
					indiv <- which(id==unique(id)[k])
					ystar[indiv] <- Y[indiv]- Z[indiv,, drop=FALSE]%*%btilde[k,]
				}

				if(!conditional){
					tree <- rpart(ystar~.,as.data.frame(X))
				}else{
					citdata <- cbind(ystar, as.data.frame(X))
					tree <- ctree(ystar~., citdata)
				}
				feuilles <- predict(tree, as.data.frame(X))
				leaf <- unique(feuilles)
				nnodes <- length(leaf)
				Phi <- matrix(0,length(Y), nnodes)

				for(p in 1:nnodes){
					w <- which(feuilles==leaf[p])
					Phi[unique(w),p] <- 1
				}

				beta <- Moy(id,Btilde,sigmahat,Phi,Y,Z) ### fit des feuilles

				# Replace the predicted response at each terminal node of the tree with the estimated population level predicted response from the LMM fit
				fhat <- feuilles
				for(p in seq_len(nnodes)){
					fhat[which(feuilles==leaf[p])] <- beta[p]
				}


				# through next steps update btilde and other parameters
				for(k in 1:nind){ ### calcul des effets al?atoires par individu
					indiv <- which(id==unique(id)[k])
					V <- Z[indiv,, drop=FALSE] %*% Btilde %*% t(Z[indiv,, drop=FALSE]) + diag(as.numeric(sigmahat),length(indiv),length(indiv))
					btilde[k,] <- Btilde %*% t(Z[indiv,, drop=FALSE]) %*% solve(V) %*% (Y[indiv] - fhat[indiv])
					epsilonhat[indiv] <- Y[indiv] - fhat[indiv] - Z[indiv,, drop=FALSE] %*% btilde[k,]
				}
				sigm <- sigmahat
				sigmahat <- sig(sigmahat,id, Z, epsilonhat, Btilde) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
				Btilde  <- bay(btilde,Btilde,Z,id,sigm) #### MAJ des param?tres de la variance des effets al?atoires.
				# compute and check the log-likelihood
				Vrai <- c(Vrai, logV(Y,fhat,Z[,,drop=FALSE],time,id,Btilde,0,sigmahat,sto))
				if(i > 1) inc <- (Vrai[i-1]-Vrai[i])/Vrai[i-1]
				if(inc <  delta) {
					print(paste0("stopped after ", i, " iterations."))
					if(!conditional){
						lee <- which(tree$frame[,"var"]=="<leaf>")
						for(k in 1:nnodes){
							ou <- which(tree$frame[,"yval"]==leaf[k])
							w <- intersect(ou,lee)
							tree$frame[w,"yval"] <- beta[k]
						}
						sortie <- list(forest=tree,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, vraisemblance = Vrai,id=id, time=time)
					}else{
						# the tree produces a non-correct prediction, it needs to be corrected with beta values
						sortie <- list(forest=tree,beta=beta,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, vraisemblance = Vrai,id=id, time=time)
					}
					class(sortie) <- "longituRF"
					return(sortie)
				}
			}
			if(!conditional){
				lee <- which(tree$frame[,"var"]=="<leaf>")
				for(k in 1:nnodes){
					ou <- which(tree$frame[,"yval"]==leaf[k])
					w <- intersect(ou,lee)
					tree$frame[w,"yval"] <- beta[k]
				}
				sortie <- list(forest=tree,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, Vraisemblance=Vrai, time =time, id=id )
			}else{
				sortie <- list(forest=tree,beta=beta,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, Vraisemblance=Vrai, time =time, id=id )
			}
			class(sortie) <- "longituRF"
			return(sortie)
		}
	}
	for(i in 1:iter){
		ystar <- rep(0,length(Y))
		for(k in 1:nind){ #### on retrace les effets al?atoires
			indiv <- which(id==unique(id)[k])
			ystar[indiv] <- Y[indiv]- Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[indiv]
		}

		if(!conditional){
			tree <- rpart(ystar~.,as.data.frame(X))
		}else{
			citdata <- cbind(ystar, as.data.frame(X))
			tree <- ctree(ystar~., citdata)
		}
		feuilles <- predict(tree,as.data.frame(X))
		leaf <- unique(feuilles)
		nnodes <- length(leaf) # predicted values at terminal nodes
		Phi <- matrix(0,length(Y), nnodes)

		for(p in seq_len(nnodes)){
			w <- which(feuilles==leaf[p])
			Phi[unique(w),p] <- 1
		}

		beta <- Moy_sto(id,Btilde,sigmahat,Phi,Y,Z,sto,time,sigma2) ### fit des feuilles

		fhat <- feuilles
		for(p in seq_len(nnodes)){
			fhat[which(feuilles==leaf[p])] <- beta[p]
		}

		for(k in 1:nind){ ### calcul des effets al?atoires par individu
			indiv <- which(id==unique(id)[k])
			K <- sto_analysis(sto,time[indiv])
			V <- Z[indiv,, drop=FALSE] %*% Btilde %*% t(Z[indiv,, drop=FALSE]) + diag(as.numeric(sigmahat),length(indiv),length(indiv)) + sigma2 * K
			btilde[k,] <- Btilde %*% t(Z[indiv,, drop=FALSE]) %*% solve(V) %*% (Y[indiv] - fhat[indiv])
			omega[indiv] <- sigma2 * K %*% solve(V) %*% (Y[indiv] - fhat[indiv])
			epsilonhat[indiv] <- Y[indiv] - fhat[indiv] - Z[indiv,, drop=FALSE] %*% btilde[k,] - omega[indiv]
		}
		#### pr?diction du processus stochastique:
		sigm <- sigmahat
		B <- Btilde
		sigmahat <- sig_sto(sigmahat,id, Z, epsilonhat, Btilde, time, sigma2,sto) ##### MAJ de la variance des erreurs ! ici que doit se trouver le probl?me !
		Btilde  <- bay_sto(btilde,Btilde,Z,id,sigm, time, sigma2,sto) #### MAJ des param?tres de la variance des effets al?atoires.
		sigma2 <- gam_sto(sigm,id,Z,B,time,sigma2,sto,omega)
		Vrai <- c(Vrai, logV(Y,fhat,Z[,,drop=FALSE],time,id,Btilde,sigma2,sigmahat,sto))
		if(i > 1) inc <- (Vrai[i-1]-Vrai[i])/Vrai[i-1]
		if(inc < delta) {
			print(paste0("stopped after ", i, " iterations."))
			if(!conditional){
				lee <- which(tree$frame[,"var"]=="<leaf>")
				for(k in 1:nnodes){
					ou <- which(tree$frame[,"yval"]==leaf[k])
					w <- intersect(ou,lee)
					tree$frame[w,"yval"] <- beta[k]
				}

				sortie <- list(forest=tree,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), id_omega=id_omega, omega=omega, sigma_sto =sigma2, time = time, sto= sto,Vraisemblance=Vrai,id=id)
			}else{
				sortie <- list(forest=tree,beta=beta,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), id_omega=id_omega, omega=omega, sigma_sto =sigma2, time = time, sto= sto,Vraisemblance=Vrai,id=id)
			}
			class(sortie) <- "longituRF"
			return(sortie)
		}
	}
	if(!conditional){
		lee <- which(tree$frame[,"var"]=="<leaf>")
		for(k in 1:nnodes){
			ou <- which(tree$frame[,"yval"]==leaf[k])
			w <- intersect(ou,lee)
			tree$frame[w,"yval"] <- beta[k]
		}
		sortie <- list(forest=tree,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), id_omega=id_omega,omega=omega, sigma_sto =sigma2, time = time, sto= sto, Vraisemblance=Vrai, id=id)
	} else if (conditional) {
		sortie <- list(forest=tree,beta=beta,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), id_omega=id_omega,omega=omega, sigma_sto =sigma2, time = time, sto= sto, Vraisemblance=Vrai, id=id)
	}
	class(sortie) <- "longituRF"
	return(sortie)
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
			V2 = log(det(V))+ t(Y[w]-f[w])%*%solve(V)%*%(Y[w]-f[w])
			if (V2<Inf){
				Vraisem <- Vraisem + V2
			}
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
#' @param time [time]
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
#' @param X [matrix]
#' @param Y [vector]
#' @param id [vector]
#' @param Z [matrix]
#' @param iter [numeric]
#' @param mtry [numeric]
#' @param ntree [numeric]
#' @param time [time]
#' @param conditional [logical]
#'
#' @import stats
#'
#' @keywords internal
opti.FBM <- function(X,Y,id,Z,iter,mtry,ntree,time,conditional){
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
				ystar[indiv] <- Y[indiv]- Z[indiv,, drop=FALSE]%*%btilde[k,]- omega[k,which(id_omega[k,]==1)]
			}

			if(!conditional){
				forest <- randomForest(as.data.frame(X),ystar,mtry=mtry,ntree=ntree) ### on construit l'arbre
				fhat <- predict(forest, as.data.frame(X)) #### pr?diction avec l'arbre
			} else if(conditional){
				citdata <- cbind(ystar, as.data.frame(X))
				forest <- cforest(ystar ~ ., data = citdata, controls = cforest_unbiased(mtry = mtry, ntree = ntree))
				fhat <- predict(forest, OOB = TRUE, type = "response") #### pr?diction avec l'arbre
			}
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
#' @param X [matrix]
#' @param Y [vector]
#' @param id [vector]
#' @param Z [matrix]
#' @param iter [numeric]
#' @param mtry [numeric]
#' @param ntree [numeric]
#' @param time [time]
#' @param conditional [logical]
#'
#' @import stats
#'
#' @keywords internal
opti.FBMreem <- function(X,Y,id,Z,iter,mtry,ntree,time,conditional){
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
				ystar[indiv] <- Y[indiv]- Z[indiv,,drop=FALSE]%*%btilde[k,]- omega[k,which(id_omega[k,]==1)]
			}
			if(!conditional){
				forest <- randomForest(as.data.frame(X), ystar,mtry=mtry,ntree=ntree, keep.inbag=TRUE)
				f1 <- predict(forest, as.data.frame(X), nodes=TRUE)
				trees <- attributes(f1)$nodes
				inbag <- forest$inbag
			} else if (conditional){
				citdata <- cbind(ystar, as.data.frame(X))
				forest <- cforest(ystar ~ ., data = citdata, controls = cforest_unbiased(mtry = mtry, ntree = ntree))
				f1 <- predict(forest, OOB = TRUE, type = "response")
				trees <- as.matrix(as.data.frame(forest@where,
								 row.names = row.names(citdata),
								 col.names = paste0("V", seq_len(ntree))))
				inbag <- as.matrix(as.data.frame(forest@weights,
								 row.names = row.names(citdata),
								 col.names = paste0("V", seq_len(ntree))))
			}
			K <- ntree

			matrice.pred <- matrix(NA,length(Y),ntree)
			for (k in 1:K){
				indii <- unique(trees[,k])
				nnodes <- length(indii)
				Phi <- matrix(0,length(Y),nnodes)
				for (l in seq_len(nnodes)){
					w <- which(trees[,k]==indii[l])
					Phi[w,l] <- 1
				}
				oobags <- unique(which(inbag[,k]==0))
				beta <- Moy_fbm(id[-oobags],Btilde,sigmahat,Phi[-oobags,],ystar[-oobags],Z[-oobags,,drop=FALSE], h ,time[-oobags], sigma2)
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
				V <- Z[indiv,,drop=FALSE]%*%Btilde%*%t(Z[indiv,,drop=FALSE])+diag(rep(sigmahat,length(Y[indiv])))+ sigma2*K
				btilde[k,] <- Btilde%*%t(Z[indiv,])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
			}
			#### pr?diction du processus stochastique:
			for (k in 1:length(id_btilde)){
				indiv <- which(id==unique(id)[k])
				K <- cov.fbm(time[indiv], h)
				V <- Z[indiv,,drop=FALSE]%*%Btilde%*%t(Z[indiv,,drop=FALSE])+diag(rep(sigmahat,length(Y[indiv])))+sigma2*K
				omega[k,which(id_omega[k,]==1)] <- sigma2*K%*%solve(V)%*%(Y[indiv]-fhat[indiv])
			}

			for (k in 1:nind){
				indiv <- which(id==unique(id)[k])
				epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,,drop=FALSE]%*%btilde[k,]- omega[k,which(id_omega[k,]==1)]
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



#' Longitudinal data generator
#'
#'
#' Simulate longitudinal data according to the semi-parametric stochastic mixed-effects model given by: \deqn{Y_i(t)=f(X_i(t))+Z_i(t)\beta_i + \omega_i(t)+\epsilon_i}
#' with \eqn{Y_i(t)} the output at time \eqn{t} for the \eqn{i}th individual; \eqn{X_i(t)} the input predictors (fixed effects) at time \eqn{t} for the \eqn{i}th individual;
#' \eqn{Z_i(t)} are the random effects at time \eqn{t} for the \eqn{i}th individual; \eqn{\omega_i(t)} is a Brownian motion with volatility \eqn{\gamma^2=0.8} at time \eqn{t} for the \eqn{i}th individual; \eqn{\epsilon_i} is the residual error with
#' variance \eqn{\sigma^2=0.5}.
#' The data are simulated according to the simulations in low dimensional in the low dimensional scheme of the paper <doi:10.1177/0962280220946080>
#'
#' @param n [numeric]: Number of individuals. The default value is \code{n=50}.
#' @param p [numeric]: Number of predictors. The default value is \code{p=6}.
#' @param G [numeric]: Number of groups of predictors with temporal behavior, generates \code{p-G} input variables with no temporal behavior.
#'
#' @import mvtnorm
#' @import latex2exp
#'
#' @return a list of the following elements: \itemize{
#' \item \code{Y:} vector of the output trajectories.
#' \item \code{X :} matrix of the fixed-effects predictors.
#' \item \code{Z:} matrix of the random-effects predictors.
#' \item \code{id: } vector of the identifiers for each individual.
#' \item \code{time: } vector the the time measurements for each individual.
#' }
#'
#' @export
#'
#' @examples
#' oldpar <- par()
#' oldopt <- options()
#' data <- DataLongGenerator(n=17, p=6,G=6) # Generate the data
#' # Let's see the output :
#' w <- which(data$id==1)
#' plot(data$time[w],data$Y[w],type="l",ylim=c(min(data$Y),max(data$Y)), col="grey")
#' for (i in unique(data$id)){
#'   w <- which(data$id==i)
#'   lines(data$time[w],data$Y[w], col='grey')
#' }
#' # Let's see the fixed effects predictors:
#' par(mfrow=c(2,3), mar=c(2,3,3,2))
#' for (i in 1:ncol(data$X)){
#'   w <- which(data$id==1)
#'   plot(data$time[w],data$X[w,i], col="grey",ylim=c(min(data$X[,i]),
#'   max(data$X[,i])),xlim=c(1,max(data$time)),main=latex2exp::TeX(paste0("$X^{(",i,")}$")))
#'   for (k in unique(data$id)){
#'     w <- which(data$id==k)
#'     lines(data$time[w],data$X[w,i], col="grey")
#'   }
#' }
#' par(oldpar)
#' options(oldopt)
#'
DataLongGenerator <- function(n=50,p=6,G=6){

	mes <-floor(4*runif(n)+8)
	time <- NULL
	id <- NULL
	nb2 <- c(1:n)
	for (i in 1:n){
		time <- c(time, seq(1,mes[i], by=1))
		id <- c(id, rep(nb2[i], length(seq(1,mes[i], by=1))))
	}

	bruit <- floor(0*p)
	bruit <- bruit+ (p-bruit)%%G
	nices <- NULL
	for (i in 1:G){
		nices <- c(nices,rep(i,(p-bruit)/G))
	}

	comportements <- matrix(0,length(time),G)
	comportements[,1] <- 2.44+0.04*(time-((time-6)^2)/(time/3))
	comportements[,2] <- 0.5*time-0.1*(time-5)^2
	comportements[,3] <- 0.25*time-0.05*(time-6)^2
	comportements[,4] <- cos((time-1)/3)
	comportements[,5] <- 0.1*time + sin(0.6*time+1.3)
	comportements[,6] <- -0.1*time^2


	X <- matrix(0,length(time), p)
	for (i in 1:(p-bruit)){
		X[,i] <- comportements[,nices[i]] + rnorm(length(time),0 ,0.2)
	}

	for (j in 1:n){
		w <- which(id==j)
		X[w,1:(p-bruit)] <- X[w,1:(p-bruit)] + rnorm(1,0,0.1)
	}

	for (i in (p-bruit):p){
		X[,i] <- rnorm(length(time),0, 3)
	}

	f <- 1.3*X[,1]^2 + 2*sqrt(abs(X[,which(nices==2)[1]]))

	sigma <- cbind(c(0.5,0.6),c(0.6,3))
	Btilde<- matrix(0,length(unique(id)),2)
	for (i in 1:length(unique(id))){
		Btilde[i,] <- rmvnorm(1, mean=rep(0,2),sigma=sigma)
	}

	Z <- as.matrix(cbind(rep(1,length(f)),2*runif(length(f))))

	effets  <- NULL
	for (i in 1:length(unique(id))){
		w <- which(id==unique(id)[i])
		effets <- c(effets, Z[w,, drop=FALSE]%*%Btilde[i,])
	}
	##### simulation de mouvemments brownien
	gam <- 0.8
	BM <- NULL
	m <- length(unique(id))
	for (i in 1:m){
		w <- which(id==unique(id)[i])
		W <- rep(0,length(w))
		t <- time[w]
		for (j in 2:length(w)){
			W[j] <- W[j-1]+sqrt(gam*(t[j]-t[j-1]))*rnorm(1,0,1)
		}
		BM <- c(BM,W)
	}

	sigma2 <- 0.5
	Y <- f + effets +rnorm(length(f),0,sigma2)+BM
	return(list(Y=Y,X=X,Z=Z,id=id, time=time))
}


#' Stability score function for (S)MERF and (S)REEMforest methods
#'
#' Computes the stability scores for (S)MERF and (S)REEMforest methods.
#'
#' @param X [matrix]: A \code{N}x\code{p} matrix containing the \code{p} predictors of the fixed effects, column codes for a predictor.
#' @param Y [vector]: A vector containing the output trajectories.
#' @param id [vector]: Is the vector of the identifiers for the different trajectories.
#' @param Z [matrix]: A \code{N}x\code{q} matrix containing the \code{q} predictor of the random effects.
#' @param time [vector]: Is the vector of the measurement times associated with the trajectories in \code{Y},\code{Z} and \code{X}.
#' @param mtry [numeric]: Number of variables ramdomly picked to split each node.
#' @param ntree [numeric]: Number of trees in the RF.
#' @param sto [string]: Defines the covariance function of the stochastic process, can be either \code{"none"} for no stochastic process, \code{"BM"} for Brownian motion, \code{OrnUhl} for standard Ornstein-Uhlenbeck process, \code{BBridge} for Brownian Bridge, \code{fbm} for Fractional Brownian motion; can also be a function defined by the user.
#' @param method [string]: Defines the method to be used, can be either "MERF" or "REEMforest".
#' @param eta [numeric]: The size of the neighborhood for the stability score. Can be a vector, in this case, returns the stability scores corresponding to all the values of the vector.
#' @param nvars [numeric]: The number of variables to consider among the most impotant variables. Can be a vector, in this case, the function returns the stability scores corresponding to all the values of the vector.
#' @param cforest [logical]: Determines if the random forest algorithm is the one implemented by Breiman (2001) and available in \code{\link[randomForest]{randomForest}}, which is the default \code{FALSE}, or if it is the implemented by Strobl et al. (2007) using conditional inference trees in \code{\link[party]{cforest}}.
#' @param conditionalVI [logical]: A logical determining whether unconditional or conditional computation of the importance is performed. Only applicable when \code{cforest==TRUE}.
#'
#' @export
#'
#' @return A matrix with all the stability scores corresponding to the eta and nvars values. The $i$th row corresponds to the $i$th value of eta while the $i$th column corresponds to the $i$ value of nvars.
#
Stability_Score <- function(X,Y,Z,id,time,mtry,ntree, sto="BM",method="MERF", eta = c(1:ncol(X)),nvars=c(1:ncol(X)), cforest = FALSE, conditionalVI = FALSE){

	if(!is.logical(cforest) || length(cforest) != 1 || is.na(cforest)){
		cforest <- FALSE
	}
	if(!is.logical(conditionalVI) || length(conditionalVI) != 1 || is.na(conditionalVI)){
		conditionalVI <- FALSE
	}

	if (method=="REEMforest"){
		sortie1 <- REEMforest(X=X,Y=Y,Z=Z,id=id,time=time,mtry=mtry,ntree=ntree,sto=sto)
		sortie2 <- REEMforest(X=X,Y=Y,Z=Z,id=id,time=time,mtry=mtry,ntree=ntree,sto=sto)
	}

	else {
		sortie1 <- MERF(X=X,Y=Y,Z=Z,id=id,time=time,mtry=mtry,ntree=ntree,sto=sto,conditional=cforest)
		sortie2 <- MERF(X=X,Y=Y,Z=Z,id=id,time=time,mtry=mtry,ntree=ntree,sto=sto,conditional=cforest)
	}

	if(!cforest){
		imp1 <- sort(sortie1$forest$importance[,1], decreasing = TRUE, index.return=TRUE)
		imp2 <- sort(sortie2$forest$importance[,1], decreasing = TRUE, index.return=TRUE)
	} else {
		imp1 <- sort(varimp(sortie1$forest, conditional=conditionalVI), decreasing = TRUE, index.return=TRUE)
		imp2 <- sort(varimp(sortie2$forest, conditional=conditionalVI), decreasing = TRUE, index.return=TRUE)
	}

	ss <- matrix(NA,length(eta),length(nvars))
	for (i in 1:length(eta)){
		for (k in 1:length(nvars)){
			nind = rep(0,nvars[k])
			for (l in 1:nvars[k]){
				variable = imp1$ix[l]
				variable2_interv = imp2$ix[c(max(1,l-eta[i]):min(nvars[k],l+eta[i]))]
				nind[l] = 1*length(which(variable2_interv == variable)>0)
			}
			ss[i,k] = mean(nind)
		}
	}
	SS <- as.data.frame(ss)
	colnames(SS) = nvars
	rownames(SS) = eta
	return(SS)
}


