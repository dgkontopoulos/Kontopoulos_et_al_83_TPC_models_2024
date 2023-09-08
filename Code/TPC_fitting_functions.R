# This script provides functions for fitting the 83 TPC models.
#
# These functions need only a data frame as input, with columns of
# 'temp' (temperature) and 'trait_value'.

#############################################
# 2ND ORDER POLYNOMIAL MODEL (3 parameters) #
#                                           #
# Parameters: a, b, d                       #
#############################################
fit_2nd_order_polynomial_3_pars <- function(dataset)
{
	
	# Set the starting value of a to the minimum trait value.
	a_start <- min(dataset$trait_value)
	
	# Set the starting value of b arbitrarily to 0.1.
	b_start <- 0.1
	
	# Set the starting value of d arbitrarily to 0.003.
	d_start <- 0.003
	
	function_to_be_fitted <- function(a, b, d, temp)
	{
		return(
			log(
				a + b * temp + d * temp^2
			)
		)
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, b, d, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,	b = 0.5 * b_start,
				d = 0.5 * d_start
			),
			start_upper = c(
				a = 1.5 * a_start,	b = 1.5 * b_start,
				d = 1.5 * d_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(-Inf, -Inf, -Inf),
			upper = c(Inf, Inf, Inf)
		)
	)
	
	return(fit)
}

#############################################
# 3RD ORDER POLYNOMIAL MODEL (4 parameters) #
#                                           #
# Parameters: a, b, d, e                    #
#############################################
fit_3rd_order_polynomial_4_pars <- function(dataset)
{
	
	# Set the starting value of a to the minimum trait value.
	a_start <- min(dataset$trait_value)
	
	# Set the starting value of b arbitrarily to 0.1.
	b_start <- 0.1
	
	# Set the starting value of d arbitrarily to 0.003.
	d_start <- 0.003
	
	# Set the starting value of e arbitrarily to 0.0000002.
	e_start <- 0.0000002
	
	function_to_be_fitted <- function(a, b, d, e, temp)
	{
		return(
			log(
				a + b * temp + d * temp^2 + e * temp^3
			)
		)
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, b, d, e, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,	b = 0.5 * b_start,
				d = 0.5 * d_start,	e = 0.5 * e_start
			),
			start_upper = c(
				a = 1.5 * a_start,	b = 1.5 * b_start,
				d = 1.5 * d_start,	e = 1.5 * e_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(-Inf, -Inf, -Inf, -Inf),
			upper = c(Inf, Inf, Inf, Inf)
		)
	)
	
	return(fit)
}

#############################################
# 4TH ORDER POLYNOMIAL MODEL (5 parameters) #
#                                           #
# Parameters: a, b, d, e, f                 #
#############################################
fit_4th_order_polynomial_5_pars <- function(dataset)
{
	
	# Set the starting value of a to the minimum trait value.
	a_start <- min(dataset$trait_value)
	
	# Set the starting value of b arbitrarily to 0.1.
	b_start <- 0.1
	
	# Set the starting value of d arbitrarily to 0.003.
	d_start <- 0.003
	
	# Set the starting value of e arbitrarily to 0.0000002.
	e_start <- 0.0000002
	
	# Set the starting value of f arbitrarily to 0.0000000001.
	f_start <- 0.0000000001
	
	function_to_be_fitted <- function(a, b, d, e, f, temp)
	{
		result <- a + b * temp + d * temp^2 + e * temp^3 + f * temp^4
		
		if ( any(result <= 0 ) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, b, d, e, f, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,	b = 0.5 * b_start,
				d = 0.5 * d_start,	e = 0.5 * e_start,
				f = 0.5 * f_start
			),
			start_upper = c(
				a = 1.5 * a_start,	b = 1.5 * b_start,
				d = 1.5 * d_start,	e = 1.5 * e_start,
				f = 1.5 * f_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(-Inf, -Inf, -Inf, -Inf, -Inf),
			upper = c(Inf, Inf, Inf, Inf, Inf)
		)
	)
	
	return(fit)
}

#############################################
# 5TH ORDER POLYNOMIAL MODEL (6 parameters) #
#                                           #
# Parameters: a, b, d, e, f, g              #
#############################################
fit_5th_order_polynomial_6_pars <- function(dataset)
{
	
	# Set the starting value of a to the minimum trait value.
	a_start <- min(dataset$trait_value)
	
	# Set the starting value of b arbitrarily to 0.1.
	b_start <- 0.1
	
	# Set the starting value of d arbitrarily to 0.003.
	d_start <- 0.003
	
	# Set the starting value of e arbitrarily to 0.0000002.
	e_start <- 0.0000002
	
	# Set the starting value of f arbitrarily to 0.0000000001.
	f_start <- 0.0000000001
	
	# Set the starting value of g arbitrarily to 0.0000000001.
	g_start <- 0.0000000001
	
	function_to_be_fitted <- function(a, b, d, e, f, g, temp)
	{
		result <- a + b * temp + d * temp^2 + e * temp^3 + f * temp^4 + g * temp^5
		
		if ( any(result <= 0 ) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, b, d, e, f, g, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,	b = 0.5 * b_start,
				d = 0.5 * d_start,	e = 0.5 * e_start,
				f = 0.5 * f_start,	g = 0.5 * g_start
			),
			start_upper = c(
				a = 1.5 * a_start,	b = 1.5 * b_start,
				d = 1.5 * d_start,	e = 1.5 * e_start,
				f = 1.5 * f_start,	g = 1.5 * g_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf),
			upper = c(Inf, Inf, Inf, Inf, Inf, Inf)
		)
	)
	
	return(fit)
}

#############################################
# ANALYTIS - ALLAHYARI MODEL (5 parameters) #
# (Allahyari, PhD thesis, 2005)             #
#                                           #
# Parameters: a, b, d, T_min, T_max         #
#############################################
fit_Analytis_Allahyari_5_pars <- function(dataset)
{

	# Set the starting value of a arbitrarily to 1.
	a_start <- 1

	# Set the starting value of b arbitrarily to 1.
	b_start <- 1
	
	# Set the starting value of d arbitrarily to 0.5.
	d_start <- 0.5
	
	# Set the starting value of T_min 1 degree below the minimum 
	# temperature in the dataset.
	T_min_start <- min(dataset$temp) - 1

	# Set the starting value of T_max 1 degree above the maximum 
	# temperature in the dataset.
	T_max_start <- max(dataset$temp) + 1
		
	function_to_be_fitted <- function(a, b, d, T_min, T_max, temp)
	{
		if ( T_min >= T_max || any(temp <= T_min) || any(temp >= T_max) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			delta <- (temp - T_min)/(T_max - T_min)
			
			return(
				log(
					a * delta^b * (1 - delta^d)
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, b, d, T_min, T_max, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,			b = 0.5 * b_start,
				d = 0.5 * d_start,			T_min = 0.5 * T_min_start,
				T_max = 0.5 * T_max_start
			),
			start_upper = c(
				a = 1.5 * a_start,			b = 1.5 * b_start,
				d = 1.5 * d_start,			T_min = 1.5 * T_min_start,
				T_max = 1.5 * T_max_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0, -20, 0),
			upper = c(Inf, Inf, Inf, 150, 150)
		)
	)
	
	return(fit)
}

################################################
# ANALYTIS - KONTODIMAS MODEL (3 parameters)   #
# (Kontodimas et al., Environ. Entomol., 2004) #
#                                              #
# Parameters: a, T_min, T_max                  #
################################################
fit_Analytis_Kontodimas_3_pars <- function(dataset)
{

	# Set the starting value of a arbitrarily to 1.
	a_start <- 1
	
	# Set the starting value of T_min 1 degree below the minimum 
	# temperature in the dataset.
	T_min_start <- min(dataset$temp) - 1

	# Set the starting value of T_max 1 degree above the maximum 
	# temperature in the dataset.
	T_max_start <- max(dataset$temp) + 1
		
	function_to_be_fitted <- function(a, T_min, T_max, temp)
	{
		if ( T_min >= T_max || any(temp <= T_min) || any(temp >= T_max) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(
				log(
					a * ( (temp - T_min)^2 ) * (T_max - temp)
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, T_min, T_max, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,		T_min = 0.5 * T_min_start,
				T_max = 0.5 * T_max_start
			),
			start_upper = c(
				B_0 = 1.5 * a_start,		T_min = 1.5 * T_min_start,
				T_max = 1.5 * T_max_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, -20, 0),
			upper = c(Inf, 150, 150)
		)
	)
	
	return(fit)
}

###################################################
# ASBURY - ANGILLETTA MODEL (6 parameters)        #
# (Asbury & Angilletta, Am. Nat., 2010)           #
#                                                 #
# Parameters: alpha, delta, epsilon, zeta, eta, E #
###################################################
fit_Asbury_Angilletta_6_pars <- function(dataset)
{
	
	# Set the Boltzmann constant.
	k <- 8.617 * 10^-5
	
	# Set the starting value of alpha arbitrarily to 1°C (274.15 K).
	alpha_start <- 274.15
	
	# Set the starting value of delta arbitrarily to 0.1.
	delta_start <- 0.1
	
	# Set the starting value of epsilon arbitrarily to 0.7.
	epsilon_start <- 0.7
	
	# Set the starting value of zeta arbitrarily to 100.
	zeta_start <- 100
	
	# Set the starting value of eta arbitrarily to 1.
	eta_start <- 1
	
	T_pk <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])

	# If the dataset has measurements before the thermal optimum, 
	# a starting value for E can be set as roughly the slope of the 
	# following regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * 1/(k * (T + 273.15))
	#
	# Otherwise, just set E to 0.6 eV.
	
	dataset_before_peak <- dataset[dataset$temp < T_pk,]
	if ( 
		nrow(dataset_before_peak) < 3 || 
		length(unique(dataset_before_peak$temp)) < 3 || 
		length(unique(dataset_before_peak$trait_value)) < 3 
	)
	{
		E_start <- 0.6
	} else
	{
		y_vals <- log(dataset_before_peak$trait_value)
		x_vals <- 1/(k * (dataset_before_peak$temp + 273.15))
		
		# Take the absolute value.
		E_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
		
		if ( E_start >= 10 )
		{
			E_start <- 0.6
		}
	}
		
	function_to_be_fitted <- function(alpha, delta, epsilon, zeta, eta, E, temp)
	{
		temp <- temp + 273.15

		# Set the Boltzmann constant.
		k <- 8.617 * 10^-5
		
		if ( any( (temp - alpha) > zeta ) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(
				log(
					eta * ( ( ( (temp-alpha)/zeta ) ^ (epsilon/delta-1) * (1 - (temp-alpha) / zeta) ^ ( (1-epsilon) / delta - 1) ) / beta( epsilon / delta, (1 - epsilon) / delta)) * exp(-E/(k * temp))
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				alpha, delta, epsilon, zeta, eta, E, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				alpha = 0.5 * alpha_start,		delta = 0.5 * delta_start,
				epsilon = 0.5 * epsilon_start,	zeta = 0.5 * zeta_start,
				eta = 0.5 * eta_start,			E = 0.5 * E_start
			),
			start_upper = c(
				alpha = 1.5 * alpha_start,		delta = 1.5 * delta_start,
				epsilon = 1.5 * epsilon_start,	zeta = 1.5 * zeta_start,
				eta = 1.5 * eta_start,			E = 1.5 * E_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0, 0, 0, 0),
			upper = c(Inf, Inf, 1, Inf, Inf, 10)
		)
	)
	
	return(fit)
}

#######################################
# ASHRAFI I MODEL (3 parameters)      #
# (Ashrafi et al., Evol. Appl., 2018) #
#                                     #
# Parameters: a, b, d                 #
#######################################
fit_Ashrafi_I_3_pars <- function(dataset)
{
	
	# Set the starting value of a to the minimum trait value.
	a_start <- min(dataset$trait_value)
	
	# Set the starting value of b arbitrarily to 0.1.
	b_start <- 0.1
	
	# Set the starting value of d arbitrarily to 0.003.
	d_start <- 0.003
	
	function_to_be_fitted <- function(a, b, d, temp)
	{
		return(
			log(
				a + b * temp^2 * log(temp) + d * temp^3
			)
		)
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, b, d, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,	b = 0.5 * b_start,
				d = 0.5 * d_start
			),
			start_upper = c(
				a = 1.5 * a_start,	b = 1.5 * b_start,
				d = 1.5 * d_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(-Inf, -Inf, -Inf),
			upper = c(Inf, Inf, Inf)
		)
	)
	
	return(fit)
}

#######################################
# ASHRAFI II MODEL (3 parameters)     #
# (Ashrafi et al., Evol. Appl., 2018) #
#                                     #
# Parameters: a, b, d                 #
#######################################
fit_Ashrafi_II_3_pars <- function(dataset)
{
	
	# Set the starting value of a to the minimum trait value.
	a_start <- min(dataset$trait_value)
	
	# Set the starting value of b arbitrarily to 0.1.
	b_start <- 0.1
	
	# Set the starting value of d arbitrarily to 0.003.
	d_start <- 0.003
	
	function_to_be_fitted <- function(a, b, d, temp)
	{
		return(
			log(
				a + b * temp^(3/2) + d * temp^2
			)
		)
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, b, d, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,	b = 0.5 * b_start,
				d = 0.5 * d_start
			),
			start_upper = c(
				a = 1.5 * a_start,	b = 1.5 * b_start,
				d = 1.5 * d_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(-Inf, -Inf, -Inf),
			upper = c(Inf, Inf, Inf)
		)
	)
	
	return(fit)
}

#######################################
# ASHRAFI III MODEL (3 parameters)    #
# (Ashrafi et al., Evol. Appl., 2018) #
#                                     #
# Parameters: a, b, d                 #
#######################################
fit_Ashrafi_III_3_pars <- function(dataset)
{
	
	# Set the starting value of a to the minimum trait value.
	a_start <- min(dataset$trait_value)
	
	# Set the starting value of b arbitrarily to 0.1.
	b_start <- 0.1
	
	# Set the starting value of d arbitrarily to 0.003.
	d_start <- 0.003
	
	function_to_be_fitted <- function(a, b, d, temp)
	{
		result <- 1 / ( a + b * exp(temp) + d * exp(-temp) )
		if ( any(result <= 0) || b == 0 || d == 0 )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, b, d, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,	b = 0.5 * b_start,
				d = 0.5 * d_start
			),
			start_upper = c(
				a = 1.5 * a_start,	b = 1.5 * b_start,
				d = 1.5 * d_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(-Inf, 0, 0),
			upper = c(Inf, Inf, Inf)
		)
	)
	
	return(fit)
}

#######################################
# ASHRAFI IV MODEL (4 parameters)     #
# (Ashrafi et al., Evol. Appl., 2018) #
#                                     #
# Parameters: a, b, d, e              #
#######################################
fit_Ashrafi_IV_4_pars <- function(dataset)
{
	
	# Set the starting value of a to the minimum trait value.
	a_start <- min(dataset$trait_value)
	
	# Set the starting value of b arbitrarily to 0.1.
	b_start <- 0.1
	
	# Set the starting value of d arbitrarily to 0.003.
	d_start <- 0.003
	
	# Set the starting value of e arbitrarily to 0.5.
	e_start <- 0.5
	
	function_to_be_fitted <- function(a, b, d, e, temp)
	{
		temp <- temp + 273.15
		
		return(
			log(
				a + b * temp + d * log(temp)^2 + e * sqrt(temp)
			)
		)
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, b, d, e, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,	b = 0.5 * b_start,
				d = 0.5 * d_start,	e = 0.5 * e_start
			),
			start_upper = c(
				a = 1.5 * a_start,	b = 1.5 * b_start,
				d = 1.5 * d_start,	e = 0.5 * e_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(-Inf, -Inf, -Inf, -Inf),
			upper = c(Inf, Inf, Inf, Inf)
		)
	)
	
	return(fit)
}

#######################################
# ASHRAFI V MODEL (4 parameters)      #
# (Ashrafi et al., Evol. Appl., 2018) #
#                                     #
# Parameters: a, b, d, e              #
#######################################
fit_Ashrafi_V_4_pars <- function(dataset)
{
	
	# Set the starting value of a to the minimum trait value.
	a_start <- min(dataset$trait_value)
	
	# Set the starting value of b arbitrarily to 0.1.
	b_start <- 0.1
	
	# Set the starting value of d arbitrarily to 0.003.
	d_start <- 0.003
	
	# Set the starting value of e arbitrarily to 0.5.
	e_start <- 0.5
	
	function_to_be_fitted <- function(a, b, d, e, temp)
	{
		temp <- temp + 273.15
		
		return(
			log(
				a + b * log(temp)^2 + d * log(temp) + ((e * log(temp))/temp)
			)
		)
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, b, d, e, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,	b = 0.5 * b_start,
				d = 0.5 * d_start,	e = 0.5 * e_start
			),
			start_upper = c(
				a = 1.5 * a_start,	b = 1.5 * b_start,
				d = 1.5 * d_start,	e = 0.5 * e_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(-Inf, -Inf, -Inf, -Inf),
			upper = c(Inf, Inf, Inf, Inf)
		)
	)
	
	return(fit)
}

##############################
# ATKIN MODEL (3 parameters) #
# (Atkin et al., 2005)       #
#                            #
# Parameters: B_0, a, b      #
##############################
fit_Atkin_3_pars <- function(dataset)
{
	
	# Set the minimum trait measurement at the rise of the TPC as the 
	# starting value for B_0.
	B_0_start <- min(dataset$trait_value[dataset$temp == min(dataset$temp)])
	
	# Set the starting value of a arbitrarily to 3.09.
	a_start <- 3.09
	
	# Set the starting value of b arbitrarily to 0.043.
	b_start <- 0.043
	
	function_to_be_fitted <- function(B_0, a, b, temp)
	{
		result <- B_0 * (a - b * temp)^(temp/10)
		
		if ( any(is.na(result)) || any(result <= 0) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_0, a, b, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_0 = 0.5 * B_0_start,	a = 0.5 * a_start,
				b = 0.5 * b_start
			),
			start_upper = c(
				B_0 = 1.5 * B_0_start,	a = 1.5 * a_start,
				b = 1.5 * b_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0),
			upper = c(Inf, Inf, Inf)
		)
	)
	
	return(fit)
}

########################################
# BILINEAR MODEL (4 parameters)        #
# (Olsen et al., Agron. J., 1993)      #
#                                      #
# Parameters: B_pk, T_min, T_max, T_pk #
########################################
fit_Bilinear_4_pars <- function(dataset)
{

	# Set the starting value of B_pk to the maximum trait value in the 
	# dataset.
	B_pk_start <- max(dataset$trait_value)

	# Set the starting value of T_min 1 degree below the minimum 
	# temperature in the dataset.
	T_min_start <- min(dataset$temp) - 1

	# Set the starting value of T_max 1 degree above the maximum 
	# temperature in the dataset.
	T_max_start <- max(dataset$temp) + 1
	
	# Set the starting value of T_pk to the temperature at which the 
	# trait reaches its maximum value.
	T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
		
	function_to_be_fitted <- function(B_pk, T_min, T_max, T_pk, temp)
	{
		if ( T_min >= T_max || any(temp < T_min) || any(temp > T_max) || T_min >= T_pk || T_max <= T_pk )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			result <- rep(NA, length(temp))
			
			temps_below_T_pk <- which(temp <= T_pk)
			temps_above_T_pk <- which(temp > T_pk)
		
			result[temps_below_T_pk] <-  B_pk * (temp[temps_below_T_pk] - T_min) / (T_pk - T_min)
			result[temps_above_T_pk] <-  B_pk * (T_max - temp[temps_above_T_pk]) / (T_max - T_pk)
			
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_pk, T_min, T_max, T_pk, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_pk = 0.5 * B_pk_start,	T_min = 0.5 * T_min_start,
				T_max = 0.5 * T_max_start,	T_pk = 0.5 * T_pk_start
			),
			start_upper = c(
				B_pk = 1.5 * B_pk_start,	T_min = 1.5 * T_min_start,
				T_max = 1.5 * T_max_start,	T_pk = 1.5 * T_pk_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, -20, 0, 0),
			upper = c(Inf, 150, 150, 150)
		)
	)
	
	return(fit)
}

##############################################
# BOATMAN MODEL (5 parameters)               #
# (Boatman et al., PLoS One, 2017)           #
#                                            #
# Parameters: B_pk, T_min, T_max, theta, phi #
##############################################
fit_Boatman_5_pars <- function(dataset)
{
	
	# Set the starting value of B_pk to the maximum trait value in the 
	# dataset.
	B_pk_start <- max(dataset$trait_value)
	
	# Set the starting value of T_min 1 degree below the minimum 
	# temperature in the dataset.
	T_min_start <- min(dataset$temp) - 1

	# Set the starting value of T_max 1 degree above the maximum 
	# temperature in the dataset.
	T_max_start <- max(dataset$temp) + 1
	
	# Set the starting value of theta arbitrarily to 1.5.
	theta_start <- 1.5
	
	# Set the starting value of phi arbitrarily to 0.35.
	phi_start <- 0.35

	function_to_be_fitted <- function(B_pk, T_min, T_max, theta, phi, temp)
	{
		result <- B_pk * ( sin( pi * ( (temp - T_min) / (T_max - T_min) ) ^ theta ) ) ^ phi
		
		if ( any(result > B_pk) || any(result <= 0) || any(temp >= T_max) || any(temp <= T_min) || T_min >= T_max )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_pk, T_min, T_max, theta, phi, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_pk = 0.5 * B_pk_start,	T_min = 0.5 * T_min_start,
				T_max = 0.5 * T_max_start,	theta = 0.5 * theta_start,
				phi = 0.5 * phi_start
			),
			start_upper = c(
				B_pk = 1.5 * B_pk_start,	T_min = 1.5 * T_min_start,
				T_max = 1.5 * T_max_start,	theta = 1.5 * theta_start,
				phi = 1.5 * phi_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, -20, 0, 0, 0),
			upper = c(Inf, 150, 150, Inf, Inf)
		)
	)
	
	return(fit)
}

############################################
# BRIERE I MODEL (3 parameters)            #
# (Briere et al., Environ. Entomol., 1999) #
#                                          #
# Parameters: a, T_min, T_max              #
############################################
fit_Briere_I_3_pars <- function(dataset)
{

	# Set the starting value of a arbitrarily to 1.
	a_start <- 1
	
	# Set the starting value of T_min 1 degree below the minimum 
	# temperature in the dataset.
	T_min_start <- min(dataset$temp) - 1

	# Set the starting value of T_max 1 degree above the maximum 
	# temperature in the dataset.
	T_max_start <- max(dataset$temp) + 1
	
	function_to_be_fitted <- function(a, T_min, T_max, temp)
	{
		if ( T_min >= T_max || any(temp < T_min) || any(temp > T_max) )
		{
			return(rep(1e10, length(temp)))
		} else
		{	
			return(
				log(
					a * temp * ( temp - T_min ) * sqrt( T_max - temp )
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, T_min, T_max, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,		T_min = 0.5 * T_min_start,
				T_max = 0.5 * T_max_start
			),
			start_upper = c(
				a = 1.5 * a_start,		T_min = 1.5 * T_min_start,
				T_max = 1.5 * T_max_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, -20, 0),
			upper = c(Inf, 150, 150)
		)
	)
	
	return(fit)
}

############################################
# BRIERE II MODEL (4 parameters)           #
# (Briere et al., Environ. Entomol., 1999) #
#                                          #
# Parameters: a, T_min, T_max, m           #
############################################
fit_Briere_II_4_pars <- function(dataset)
{

	# Set the starting value of a arbitrarily to 1.
	a_start <- 1
	
	# Set the starting value of T_min 1 degree below the minimum 
	# temperature in the dataset.
	T_min_start <- min(dataset$temp) - 1

	# Set the starting value of T_max 1 degree above the maximum 
	# temperature in the dataset.
	T_max_start <- max(dataset$temp) + 1
	
	# Set the starting value of m arbitrarily to 2.
	m_start <- 2
	
	function_to_be_fitted <- function(a, T_min, T_max, m, temp)
	{
		if ( T_min >= T_max || any(temp < T_min) || any(temp > T_max) )
		{
			return(rep(1e10, length(temp)))
		} else
		{	
			return(
				log(
					a * temp * ( temp - T_min ) * ( T_max - temp )^(1/m)
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, T_min, T_max, m, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,			T_min = 0.5 * T_min_start,
				T_max = 0.5 * T_max_start,	m = 0.5 * m_start
			),
			start_upper = c(
				a = 1.5 * a_start,			T_min = 1.5 * T_min_start,
				T_max = 1.5 * T_max_start,	m = 1.5 * m_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, -20, 0, 0),
			upper = c(Inf, 150, 150, Inf)
		)
	)
	
	return(fit)
}

#############################################
# CARDINAL TEMPERATURE MODEL (4 parameters) #
# (Lobry et al., Binary, 1991)              #
#                                           #
# Parameters: B_pk, T_pk, T_min, T_max      #
#############################################
fit_Cardinal_Temperature_4_pars <- function(dataset)
{

	# Set the starting value of B_pk to the maximum trait value in the 
	# dataset.
	B_pk_start <- max(dataset$trait_value)

	# Set the starting value of T_pk to the temperature at which the 
	# trait reaches its maximum value.
	T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])

	# Set the starting value of T_min to the minimum temperature in the 
	# dataset.
	T_min_start <- min(dataset$temp)

	# Set the starting value of T_max to the maximum temperature in the 
	# dataset.
	T_max_start <- max(dataset$temp)
			
	function_to_be_fitted <- function(B_pk, T_pk, T_min, T_max, temp)
	{
		result <- B_pk * ( 1 - ( ((temp - T_pk)^2) / ( ((temp - T_pk)^2) + temp * (T_max + T_min - temp) - T_max * T_min ) ) )
		
		if ( T_min >= T_max || T_pk <= T_min || T_pk >= T_max || any(temp < T_min) || any(temp > T_max) || any(result <= 0) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_pk, T_pk, T_min, T_max, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_pk = 0.5 * B_pk_start,	T_pk = 0.5 * T_pk_start,
				T_min = 0.5 * T_min_start,	T_max = 0.5 * T_max_start
			),
			start_upper = c(
				B_pk = 1.5 * B_pk_start,	T_pk = 1.5 * T_pk_start,
				T_min = 1.5 * T_min_start,	T_max = 1.5 * T_max_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, -20, 0),
			upper = c(Inf, 150, 150, 150)
		)
	)
	
	return(fit)
}

#############################################################
# CARDINAL TEMPERATURE MODEL WITH INFLECTION (4 parameters) #
# (Rosso et al., J. Theor. Biol., 1993)                     #
#                                                           #
# Parameters: B_pk, T_pk, T_min, T_max                      #
#############################################################
fit_Cardinal_Temperature_with_inflection_4_pars <- function(dataset)
{

	# Set the starting value of B_pk to the maximum trait value in the 
	# dataset.
	B_pk_start <- max(dataset$trait_value)

	# Set the starting value of T_pk to the temperature at which the 
	# trait reaches its maximum value.
	T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])

	# Set the starting value of T_min to the minimum temperature in the 
	# dataset.
	T_min_start <- min(dataset$temp)

	# Set the starting value of T_max to the maximum temperature in the 
	# dataset.
	T_max_start <- max(dataset$temp)
			
	function_to_be_fitted <- function(B_pk, T_pk, T_min, T_max, temp)
	{
		result <- B_pk * ((temp - T_max) * ((temp - T_min)^2))/( (T_pk - T_min) * ( (T_pk - T_min) * (temp - T_pk) - (T_pk - T_max) * (T_pk + T_min - 2 * temp)  ) )
		
		if ( any(result > B_pk) || T_min >= T_max || T_pk <= T_min || T_pk >= T_max || any(temp < T_min) || any(temp > T_max) || any(result <= 0) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_pk, T_pk, T_min, T_max, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_pk = 0.5 * B_pk_start,	T_pk = 0.5 * T_pk_start,
				T_min = 0.5 * T_min_start,	T_max = 0.5 * T_max_start
			),
			start_upper = c(
				B_pk = 1.5 * B_pk_start,	T_pk = 1.5 * T_pk_start,
				T_min = 1.5 * T_min_start,	T_max = 1.5 * T_max_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, -20, 0),
			upper = c(Inf, 150, 150, 150)
		)
	)
	
	return(fit)
}

##################################################
# DENT-LIKE MODEL (5 parameters)                 #
# (Soltani et al., Field Crops Res., 2006)       #
#                                                #
# Parameters: B_pk, T_min, T_max, T_pk_l, T_pk_u #
##################################################
fit_Dent_like_5_pars <- function(dataset)
{

	# Set the starting value of B_pk to the maximum trait value in the 
	# dataset.
	B_pk_start <- max(dataset$trait_value)

	# Set the starting value of T_min 1 degree below the minimum 
	# temperature in the dataset.
	T_min_start <- min(dataset$temp) - 1

	# Set the starting value of T_max 1 degree above the maximum 
	# temperature in the dataset.
	T_max_start <- max(dataset$temp) + 1
	
	# Set the starting value of T_pk_u to the temperature at which the 
	# trait reaches its maximum value.
	T_pk_u_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
	
	# Set the starting value of T_pk_l one degree below T_pk_u.
	T_pk_l_start <- T_pk_u_start - 1
		
	function_to_be_fitted <- function(B_pk, T_min, T_max, T_pk_l, T_pk_u, temp)
	{
		if ( T_min >= T_max || any(temp < T_min) || any(temp > T_max) || T_pk_l >= T_pk_u || T_min >= T_pk_l || T_max <= T_pk_u )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			result <- rep(NA, length(temp))
			
			temps_below_T_pk_l <- which(temp < T_pk_l)
			temps_between_T_pk_l_and_T_pk_u <- which(temp >= T_pk_l & temp <= T_pk_u)
			temps_above_T_pk_u <- which(temp > T_pk_u)
		
			result[temps_below_T_pk_l] <-  B_pk * (temp[temps_below_T_pk_l] - T_min) / (T_pk_l - T_min)
			result[temps_between_T_pk_l_and_T_pk_u] <- B_pk
			result[temps_above_T_pk_u] <-  B_pk * (T_max - temp[temps_above_T_pk_u]) / (T_max - T_pk_u)
			
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_pk, T_min, T_max, T_pk_l, T_pk_u, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_pk = 0.5 * B_pk_start,		T_min = 0.5 * T_min_start,
				T_max = 0.5 * T_max_start,		T_pk_l = 0.5 * T_pk_l_start,
				T_pk_u = 0.5 * T_pk_u_start
			),
			start_upper = c(
				B_pk = 1.5 * B_pk_start,	T_min = 1.5 * T_min_start,
				T_max = 1.5 * T_max_start,	T_pk_l = 1.5 * T_pk_l_start,
				T_pk_u = 1.5 * T_pk_u_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, -20, 0, 0, 0),
			upper = c(Inf, 150, 150, 150, 150)
		)
	)
	
	return(fit)
}

###########################################
# DOUBLE GAUSSIAN MODEL (4 parameters)    #
# (Phillips et al., J. Evol. Biol., 2014) #
#                                         #
# Parameters: B_pk, T_pk, a, b            #
###########################################
fit_double_Gaussian_4_pars <- function(dataset)
{
	
	# Set the starting value of B_pk to the maximum trait value in the 
	# dataset.
	B_pk_start <- max(dataset$trait_value)
	
	# Set the starting value of T_pk to the temperature at which the 
	# trait reaches its maximum value.
	T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
	
	# Set the starting value of a arbitrarily to 1.
	a_start <- 1
	
	# Set the starting value of b arbitrarily to 0.4.
	b_start <- 0.4
	
	function_to_be_fitted <- function(B_pk, T_pk, a, b, temp)
	{
		result <- rep(NA, length(temp))
		
		temps_below_T_pk <- which(temp < T_pk)
		temps_above_T_pk <- which(temp >= T_pk)
		
		result[temps_below_T_pk] <-  B_pk ^ ( exp(1) - ( ( - (temp[temps_below_T_pk] - T_pk)^2 )/( 2 * a^2 ) ) )
		result[temps_above_T_pk] <-  B_pk ^ ( exp(1) - ( ( - (temp[temps_above_T_pk] - T_pk)^2 )/( 2 * (a * b)^2 ) ) )
		
		if ( any( result > B_pk^(exp(1)) ) || any( result <= 0 ) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_pk, T_pk, a, b, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_pk = 0.5 * B_pk_start,	T_pk = 0.5 * T_pk_start,
				a = 0.5 * a_start,			b = 0.5 * b_start
			),
			start_upper = c(
				B_pk = 1.5 * B_pk_start,	T_pk = 1.5 * T_pk_start,
				a = 1.5 * a_start,			b = 1.5 * b_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0, 0),
			upper = c(Inf, 150, Inf, 1)
		)
	)
	
	return(fit)
}

####################################################
# ENZYME-ASSISTED ARRHENIUS MODEL (5 parameters)   #
# (DeLong et al., Ecol. Evol., 2017)               #
#                                                  #
# Parameters: a, E_b, E_DH, T_m, E_DCp             #
####################################################
fit_Enzyme_Assisted_Arrhenius_5_pars <- function(dataset)
{
	
	# Set the starting value of a arbitrarily to 1.
	a_start <- 1

	# Set the starting value of E_b arbitrarily to 0.01.
	E_b_start <- 0.01
	
	# Set the starting value of E_DH arbitrarily to 2.
	E_DH_start <- 2
	
	# Set the starting value of E_DCp arbitrarily to 0.1.
	E_DCp_start <- 0.1
	
	# Set the starting value of T_m 1 degree above the highest 
	# temperature in the dataset.
	T_m_start <- max(dataset$temp) + 1
		
	function_to_be_fitted <- function(a, E_b, E_DH, E_DCp, T_m, temp)
	{
		temp <- temp + 273.15
		T_m <- T_m + 273.15

		# Set the Boltzmann constant.
		k <- 8.617 * 10^-5
		
		if ( any(temp > T_m) )
		{
			return(rep(1e10, length(temp)))
		} else
		{		
			return(
				log(
					a * exp( (- (E_b - (E_DH * (1 - (temp/T_m) ) + E_DCp * (temp - T_m - temp * log(temp/T_m) ) ) ) ) / (k * temp) ) 
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, E_b, E_DH, E_DCp, T_m, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,			E_b = 0.5 * E_b_start,
				E_DH = 0.5 * E_DH_start,	E_DCp = 0.5 * E_DCp_start,
				T_m = 0.5 * T_m_start
			),
			start_upper = c(
				a = 1.5 * a_start,			E_b = 1.5 * E_b_start,
				E_DH = 1.5 * E_DH_start,	E_DCp = 1.5 * E_DCp_start,
				T_m = 1.5 * T_m_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, -Inf, 0, 0),
			upper = c(Inf, Inf, Inf, Inf, 150)
		)
	)
	
	return(fit)
}

############################################
# EUBANK MODEL (3 parameters)              #
# (Eubank et al., Environ. Entomol., 1973) #
#                                          #
# Parameters: a, T_pk, b                   #
############################################
fit_Eubank_3_pars <- function(dataset)
{
	# Set the starting value of a arbitrarily to 300.
	a_start <- 300
	
	# Set the starting value of b arbitrarily to 50.
	b_start <- 50

	# Set the starting value of T_pk to the temperature at which the 
	# trait reaches its maximum value.
	T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
	
	function_to_be_fitted <- function(a, T_pk, b, temp)
	{
		return(
			log(
				a / ((temp - T_pk)^2 + b)
			)
		)
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, T_pk, b, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,	T_pk = 0.5 * T_pk_start,
				b = 0.5 * b_start
			),
			start_upper = c(
				a = 1.5 * a_start,	T_pk = 1.5 * T_pk_start,
				b = 1.5 * b_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0),
			upper = c(Inf,150, Inf)
		)
	)
	
	return(fit)
}

########################################################
# EXPONENTIALLY MODIFIED GAUSSIAN MODEL (5 parameters) #
# (Woods et al., Methods Ecol. Evol., 2018)            #
#                                                      #
# Parameters: a, b, d, e, f                            #
########################################################
fit_exponentially_modified_Gaussian_5_pars <- function(dataset)
{
	
	# Set the starting value of a arbitrarily to 1.
	a_start <- 1
	
	# Set the starting value of b arbitrarily to 5.
	b_start <- 5
	
	# Set the starting value of d arbitrarily to 3.
	d_start <- 3
	
	# Set the starting value of e arbitrarily to 2.
	e_start <- 2
	
	# Set the starting value of f to the minimum trait value in the 
	# dataset.
	f_start <- min(dataset$trait_value)

	function_to_be_fitted <- function(a, b, d, e, f, temp)
	{
		result <- ( ( a * d * sqrt(2 * pi) ) / ( 2 * e ) ) * exp( ( ( b - temp )/e ) + ( d^2/( 2 * e^2 ) ) ) * ( ( e/abs(e) ) - erf( ( (b - temp)/( sqrt(2) * d ) ) + ( d / (sqrt(2) * e) ) ) ) + f
	
		if ( !( min(result) %in% c(result[1], result[length(result)]) ) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, b, d, e, f, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,	b = 0.5 * b_start,
				d = 0.5 * d_start,	e = 0.5 * e_start,
				f = 0.5 * f_start
			),
			start_upper = c(
				a = 1.5 * a_start,	b = 1.5 * b_start,
				d = 1.5 * d_start,	e = 1.5 * e_start,
				f = 1.5 * f_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(-Inf, -Inf, -Inf, -Inf, -Inf),
			upper = c(Inf, Inf, Inf, Inf, Inf)
		)
	)
	
	return(fit)
}

########################################
# EXTENDED BRIERE MODEL (5 parameters) #
# (Cruz-Loya et al., bioRxiv, 2020)    #
#                                      #
# Parameters: a, T_min, T_max, b, d    #
########################################
fit_extended_Briere_5_pars <- function(dataset)
{

	# Set the starting value of a arbitrarily to 1.
	a_start <- 1
	
	# Set the starting value of T_min 1 degree below the minimum 
	# temperature in the dataset.
	T_min_start <- min(dataset$temp) - 1

	# Set the starting value of T_max 1 degree above the maximum 
	# temperature in the dataset.
	T_max_start <- max(dataset$temp) + 1
	
	# Set the starting value of b arbitrarily to 1.
	b_start <- 1
	
	# Set the starting value of d arbitrarily to 1.5
	d_start <- 1.5
	
	function_to_be_fitted <- function(a, T_min, T_max, b, d, temp)
	{
		if ( T_min >= T_max || any(temp < T_min) || any(temp > T_max) )
		{
			return(rep(1e10, length(temp)))
		} else
		{	
			return(
				log(
					a * temp * ((temp - T_min)^b) * ((T_max - temp)^d)
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, T_min, T_max, b, d, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,			T_min = 0.5 * T_min_start,
				T_max = 0.5 * T_max_start,	b = 0.5 * b_start,
				d = 0.5 * d_start
			),
			start_upper = c(
				a = 1.5 * a_start,			T_min = 1.5 * T_min_start,
				T_max = 1.5 * T_max_start,	b = 1.5 * b_start,
				d = 1.5 * d_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, -20, 0, 0, 0),
			upper = c(Inf, 150, 150, Inf, Inf)
		)
	)
	
	return(fit)
}

#################################################
# EXTENDED JOHNSON - LEWIN MODEL (5 parameters) #
# (Guenther, Ecol. Appl., 1997)                 #
#                                               #
# Parameters: B_0, E, T_pk, E_D, a              #
#################################################
fit_extended_Johnson_Lewin_5_pars <- function(dataset)
{
	
	# Set the Boltzmann constant.
	k <- 8.617 * 10^-5
	
	# Set the minimum trait measurement at the rise of the TPC as the 
	# starting value for B_0.
	B_0_start <- min(dataset$trait_value[dataset$temp == min(dataset$temp)])

	# Set the starting value of T_pk to the temperature at which the 
	# trait reaches its maximum value.
	T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])

	# If the dataset has measurements before the thermal optimum, 
	# a starting value for E can be set as roughly the slope of the 
	# following regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * 1/(k * (T + 273.15))
	#
	# Otherwise, just set E to 0.6 eV.
	
	dataset_before_peak <- dataset[dataset$temp < T_pk_start,]
	if ( 
		nrow(dataset_before_peak) < 3 || 
		length(unique(dataset_before_peak$temp)) < 3 || 
		length(unique(dataset_before_peak$trait_value)) < 3 
	)
	{
		E_start <- 0.6
	} else
	{
		y_vals <- log(dataset_before_peak$trait_value)
		x_vals <- 1/(k * (dataset_before_peak$temp + 273.15))
		
		# Take the absolute value.
		E_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
		
		if ( E_start >= 10 )
		{
			E_start <- 0.6
		}
	}
	
	# If the dataset has measurements after the thermal optimum, 
	# a starting value for E_D can be set as the slope of the following 
	# regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * 1/(k * (T + 273.15))
	#
	# Otherwise, just set E_D to 3 eV.
	dataset_after_peak <- dataset[dataset$temp > T_pk_start,]
	if ( 
		nrow(dataset_after_peak) < 3 || 
		length(unique(dataset_after_peak$temp)) < 3 || 
		length(unique(dataset_after_peak$trait_value)) < 3 
	)
	{
		E_D_start <- 3
	} else
	{
		y_vals <- log(dataset_after_peak$trait_value)
		x_vals <- 1/(k * (dataset_after_peak$temp + 273.15))
		
		# Take the absolute value.
		E_D_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
		
		if ( E_D_start >= 50 )
		{
			E_D_start <- 3
		}
	}
	
	if ( E_start >= E_D_start )
	{
		E_start <- 0.9 * E_D_start
	}
	
	# Set the starting value of a arbitrarily to 1.
	a_start <- 1
	
	function_to_be_fitted <- function(B_0, E, T_pk, E_D, a, temp)
	{
		temp <- temp + 273.15

		# Set the Boltzmann constant.
		k <- 8.617 * 10^-5
				
		if ( E == E_D )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(
				log(
					B_0 * temp * exp((-E/k) * (1/temp)) / ( a + ( E/( E_D - E ) ) * exp( (E_D/k) * (1/T_pk - 1/temp) ) )
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_0, E, T_pk, E_D, a, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_0 = 0.5 * B_0_start,				E = 0.5 * E_start,
				T_pk = 0.5 * T_pk_start + 273.15,	E_D = 0.5 * E_D_start,
				a = 0.5 * a_start
			),
			start_upper = c(
				B_0 = 1.5 * B_0_start,				E = 1.5 * E_start,
				T_pk = 1.5 * T_pk_start + 273.15,	E_D = 1.5 * E_D_start,
				a = 1.5 * a_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 273.15, 0, 0),
			upper = c(Inf, 10, 273.15 + 150, 50, Inf)
		)
	)
	
	return(fit)
}

######################################################
# EXTENDED SHARPE - SCHOOLFIELD MODEL (7 parameters) #
# (Guenther, Ecol. Appl., 1997)                      #
#                                                    #
# Parameters: B_0, DH_A, DH_L, T_L50, DH_H, T_H50, a #
######################################################
fit_extended_Sharpe_Schoolfield_7_pars <- function(dataset)
{
	
	# Set the universal gas constant.
	R <- 1.987
	
	# Set the minimum trait measurement at the rise of the TPC as the 
	# starting value for B_0.
	B_0_start <- min(dataset$trait_value[dataset$temp == min(dataset$temp)])
	
	# If the dataset has measurements before the thermal optimum, 
	# a starting value for DH_A can be set as the slope of the following 
	# regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * 1/(R * (T + 273.15))
	#
	# Otherwise, just set DH_A to 15000.
	T_pk <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
	
	dataset_before_peak <- dataset[dataset$temp < T_pk,]
	if ( 
		nrow(dataset_before_peak) < 3 || 
		length(unique(dataset_before_peak$temp)) < 3 || 
		length(unique(dataset_before_peak$trait_value)) < 3 
	)
	{
		DH_A_start <- 15000
	} else
	{
		y_vals <- log(dataset_before_peak$trait_value)
		x_vals <- 1/(R * (dataset_before_peak$temp + 273.15))
		
		# Take the absolute value.
		DH_A_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
	}
	
	# If the dataset has measurements after the thermal optimum, 
	# a starting value for DH_H can be set as roughly the slope of the 
	# following regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * 1/(R * (T + 273.15))
	#
	# Otherwise, just set DH_H to 50000.
	dataset_after_peak <- dataset[dataset$temp > T_pk,]
	if ( 
		nrow(dataset_after_peak) < 3 || 
		length(unique(dataset_after_peak$temp)) < 3 || 
		length(unique(dataset_after_peak$trait_value)) < 3 
	)
	{
		DH_H_start <- 50000
	} else
	{
		y_vals <- log(dataset_after_peak$trait_value)
		x_vals <- 1/(R * (dataset_after_peak$temp + 273.15))
		
		# Take the absolute value.
		DH_H_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
	}
	
	# Set the starting value of DH_L arbitrarily to 70000.
	DH_L_start <- 70000
	
	# T_H50 should be close to the thermal optimum (T_pk), so set the 
	# starting value of T_H50 = T_pk.
	T_H50_start <- T_pk + 273.15
	
	# Set the starting value of T_L50 arbitrarily to the minimum 
	# temperature in the dataset.
	T_L50_start <- min(dataset$temp) + 273.15
	
	# Set the starting value of a arbitrarily to 1.
	a_start <- 1
	
	function_to_be_fitted <- function(B_0, DH_A, DH_L, T_L50, DH_H, T_H50, a, temp)
	{
		temp <- temp + 273.15

		# Set the universal gas constant.
		R <- 1.987
		
		# Set the reference temperature to 0°C.
		T_ref <- 273.15
			
		result <- (B_0 * (temp/T_ref) * exp( (DH_A/R) * ((1/T_ref) - (1/temp)) ) ) / ( a + exp( (-DH_L/R) * ((1/T_L50) - (1/temp)) ) + exp( (DH_H/R) * ((1/T_H50) - (1/temp)) ) )
				
		if ( T_L50 >= T_H50 || is.na(result) || any(result <= 0) || DH_A == 0 || a == 0 )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_0, DH_A, DH_L, T_L50, DH_H, T_H50, a, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_0 = 0.5 * B_0_start,		DH_A = 0.5 * DH_A_start,
				DH_L = 0.5 * DH_L_start,	T_L50 = 0.5 * T_L50_start,
				DH_H = 0.5 * DH_H_start,	T_H50 = 0.5 * T_H50_start,
				a = 0.5 * a_start
			),
			start_upper = c(
				B_0 = 1.5 * B_0_start,		DH_A = 1.5 * DH_A_start,
				DH_L = 1.5 * DH_L_start,	T_L50 = 1.5 * T_L50_start,
				DH_H = 1.5 * DH_H_start,	T_H50 = 1.5 * T_H50_start,
				a = 1.5 * a_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0, -20 + 273.15, 0, 0 + 273.15, 0),
			upper = c(Inf, 300000, 300000, 150 + 273.15, 300000, 150 + 273.15, Inf)
		)
	)
	
	return(fit)
}

#####################################################
# FINSTAD - JONSSON MODEL (4 parameters)           #
# (Finstad & Jonsson, Mar. Ecol. Prog. Ser., 2012) #
#                                                   #
# Parameters: a, T_min, T_max, b                    #
#####################################################
fit_Finstad_Jonsson_4_pars <- function(dataset)
{

	# Set the starting value of a arbitrarily to 0.4
	a_start <- 0.4
	
	# Set the starting value of T_min 1 degree below the minimum 
	# temperature in the dataset.
	T_min_start <- min(dataset$temp) - 1

	# Set the starting value of T_max 1 degree above the maximum 
	# temperature in the dataset.
	T_max_start <- max(dataset$temp) + 1
		
	# Set the starting value of b arbitrarily to 0.25.
	b_start <- 0.25
		
	function_to_be_fitted <- function(a, T_min, T_max, b, temp)
	{
		if ( T_min >= T_max || any(temp < T_min) || any(temp > T_max) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(
				log(
					a * (temp - T_min) * (1 - exp(b * (temp - T_max)))
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, T_min, T_max, b, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,			T_min = 0.5 * T_min_start,
				T_max = 0.5 * T_max_start,	b = 0.5 * b_start
			),
			start_upper = c(
				a = 1.5 * a_start,			T_min = 1.5 * T_min_start,
				T_max = 1.5 * T_max_start,	b = 1.5 * b_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, -20, 0, 0),
			upper = c(1, 150, 150, 1)
		)
	)
	
	return(fit)
}

#######################################
# GAUSSIAN MODEL (3 parameters)       #
# (Angilletta, J. Therm. Biol., 2006) #
#                                     #
# Parameters: B_pk, T_pk, a           #
#######################################
fit_Gaussian_3_pars <- function(dataset)
{
	
	# Set the starting value of B_pk to the maximum trait value in the 
	# dataset.
	B_pk_start <- max(dataset$trait_value)
	
	# Set the starting value of T_pk to the temperature at which the 
	# trait reaches its maximum value.
	T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
	
	# Set the starting value of a arbitrarily to 90.
	a_start <- 90

	function_to_be_fitted <- function(B_pk, T_pk, a, temp)
	{
		return(
			log(
				B_pk * exp( - 0.5 * ( abs( temp - T_pk ) / a )^2 )
			)
		)
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_pk, T_pk, a, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_pk = 0.5 * B_pk_start,	T_pk = 0.5 * T_pk_start,
				a = 0.5 * a_start
			),
			start_upper = c(
				B_pk = 1.5 * B_pk_start,	T_pk = 1.5 * T_pk_start,
				a = 1.5 * a_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0),
			upper = c(Inf, 150, Inf)
		)
	)
	
	return(fit)
}

############################################
# GAUSSIAN - GOMPERTZ MODEL (5 parameters) #
# (Frazier et al., Am. Nat., 2006)         #
#                                          #
# Parameters: d, alpha, beta, T_pk, theta  #
############################################
fit_Gaussian_Gompertz_5_pars <- function(dataset)
{
	
	# Set the starting value of theta arbitrarily to 6, similarly to 
	# the Frazier et al. paper.
	theta_start <- 6

	# Set the starting value of T_pk to the temperature at which the 
	# trait reaches its maximum value.
	T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
	
	# Set the starting value of d to the maximum trait value in the 
	# dataset.
	d_start <- max(dataset$trait_value)

	# If the dataset has measurements before the thermal optimum, 
	# a starting value for alpha can be set as the slope of the following 
	# regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * T
	#
	# Otherwise, just arbitrarily set alpha to 0.6.
	
	dataset_before_peak <- dataset[dataset$temp < T_pk_start,]
	if ( 
		nrow(dataset_before_peak) < 3 || 
		length(unique(dataset_before_peak$temp)) < 3 || 
		length(unique(dataset_before_peak$trait_value)) < 3 
	)
	{
		alpha_start <- 0.6
	} else
	{
		y_vals <- log(dataset_before_peak$trait_value)
		x_vals <- dataset_before_peak$temp
		
		# Take the absolute value.
		alpha_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
	}
	
	# If the dataset has measurements after the thermal optimum, 
	# a starting value for beta can be set as roughly the slope of 
	# the following regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * T
	#
	# Otherwise, just arbitrarily set beta to 3.
	dataset_after_peak <- dataset[dataset$temp > T_pk_start,]
	if ( 
		nrow(dataset_after_peak) < 3 || 
		length(unique(dataset_after_peak$temp)) < 3 || 
		length(unique(dataset_after_peak$trait_value)) < 3 
	)
	{
		beta_start <- 3
	} else
	{
		y_vals <- log(dataset_after_peak$trait_value)
		x_vals <- dataset_after_peak$temp
		
		# Take the absolute value.
		beta_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
	}
	
	function_to_be_fitted <- function(d, alpha, beta, T_pk, theta, temp)
	{
		result <- d * exp(-exp( beta * (temp - T_pk) - theta ) - alpha * ( temp - T_pk )^2  )
		
		if ( any( is.nan(result) ) || any( result <= 0 ) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				d, alpha, beta, T_pk, theta, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				d = 0.5 * d_start,			alpha = 0.5 * alpha_start,
				beta = 0.5 * beta_start,	T_pk = 0.5 * T_pk_start,
				theta = 0.5 * theta_start
			),
			start_upper = c(
				d = 1.5 * d_start,			alpha = 1.5 * alpha_start,
				beta = 1.5 * beta_start,	T_pk = 1.5 * T_pk_start,
				theta = 1.5 * theta_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0, 0, -Inf),
			upper = c(Inf, Inf, Inf, 150, Inf)
		)
	)
	
	return(fit)
}

####################################
# HINSHELWOOD MODEL (4 parameters) #
# (Hinshelwood, 1946)              #
#                                  #
# Parameters: a, E_1, b, E_2       #
####################################
fit_Hinshelwood_4_pars <- function(dataset)
{
	
	# Set the universal gas constant.
	R <- 0.001987
	
	# Set the starting value of a to the minimum trait value at the 
	# low temperature end.
	a_start <- min(dataset$trait_value[dataset$temp == min(dataset$temp)])
	
	# Set the starting value of b to the minimum trait value at the 
	# high temperature end.
	b_start <- min(dataset$trait_value[dataset$temp == max(dataset$temp)])

	# If the dataset has measurements before the thermal optimum, 
	# a starting value for E_1 can be set as roughly the slope of the 
	# following regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * 1/(R * (T + 273.15)))
	#
	# Otherwise, just set E_1 to 0.6.
	T_pk <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])

	dataset_before_peak <- dataset[dataset$temp < T_pk,]
	if ( 
		nrow(dataset_before_peak) < 3 || 
		length(unique(dataset_before_peak$temp)) < 3 || 
		length(unique(dataset_before_peak$trait_value)) < 3 
	)
	{
		E_1_start <- 0.6
	} else
	{
		y_vals <- log(dataset_before_peak$trait_value)
		x_vals <- 1/(R * (dataset_before_peak$temp + 273.15))
		
		# Take the absolute value.
		E_1_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
	}
	
	# If the dataset has measurements after the thermal optimum, 
	# a starting value for E_2 can be set as the slope of the following 
	# regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * 1/(R * (T + 273.15))
	#
	# Otherwise, just set E_2 to 3.
	dataset_after_peak <- dataset[dataset$temp > T_pk,]
	if ( 
		nrow(dataset_after_peak) < 3 || 
		length(unique(dataset_after_peak$temp)) < 3 || 
		length(unique(dataset_after_peak$trait_value)) < 3 
	)
	{
		E_2_start <- 3
	} else
	{
		y_vals <- log(dataset_after_peak$trait_value)
		x_vals <- 1/(R * (dataset_after_peak$temp + 273.15))
		
		# Take the absolute value.
		E_2_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
	}
	
	function_to_be_fitted <- function(a, E_1, b, E_2, temp)
	{
		temp <- temp + 273.15

		# Set the universal gas constant.
		R <- 0.001987
		
		result <- a * exp( -E_1 / (R * temp) ) - b * exp( -E_2 / (R * temp) )

		if ( any(result <= 0 ) || E_1 == 0 || E_2 == 0 )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, E_1, b, E_2, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,	E_1 = 0.5 * E_1_start,
				b = 0.5 * b_start,	E_2 = 0.5 * E_2_start
			),
			start_upper = c(
				a = 1.5 * a_start,	E_1 = 1.5 * E_1_start,
				b = 1.5 * b_start,	E_2 = 1.5 * E_2_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0, 0),
			upper = c(Inf, Inf, Inf, Inf)
		)
	)
	
	return(fit)
}

#########################################
# HOBBS MODEL (4 parameters)            #
# (Hobbs et al., ACS Chem. Biol., 2013) #
#                                       #
# Parameters: B_0, DH, DC_p, DS         #
#########################################
fit_Hobbs_4_pars <- function(dataset)
{

	# Set the minimum trait measurement at the rise of the TPC as the 
	# starting value for B_0.
	B_0_start <- min(dataset$trait_value[dataset$temp == min(dataset$temp)])
	
	# Set the starting value of DH arbitrarily to 30000.
	DH_start <- 30000
	
	# Set the starting value of DC_p arbitrarily to 1000.
	DC_p_start <- 1000
	
	# Set the starting value of DS arbitrarily to 17.
	DS_start <- 17
	
	function_to_be_fitted <- function(B_0, DH, DC_p, DS, temp)
	{		
		temp <- temp + 273.15
		
		# Set the reference temperature to 0°C.
		T_ref <- 273.15
		
		# Set the universal gas constant.
		R <- 8.314
		
		# Set the Boltzmann constant.
		k <- 1.38 * 1e-23
		
		# Set the Planck constant.
		h <- 6.626 * 1e-34
		
		result <- B_0 * ( (k * temp) / h ) * exp( ( - ( DH - DC_p * (temp - T_ref) ) / (R * temp) ) + ( (DS - DC_p * log(temp/T_ref)) / R ) ) 

		if ( any(is.nan(result)) || any(result <= 0 ) || B_0 == 0 || DH == 0 || DC_p == 0 || DS == 0 )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_0, DH, DC_p, DS, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_0 = 0.5 * B_0_start,		DH = 0.5 * DH_start,
				DC_p = 0.5 * DC_p_start,	DS = 0.5 * DS_start
			),
			start_upper = c(
				B_0 = 1.5 * B_0_start,		DH = 1.5 * DH_start,
				DC_p = 1.5 * DC_p_start,	DS = 1.5 * DS_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0, 0),
			upper = c(Inf, Inf, Inf, Inf)
		)
	)
	
	return(fit)
}

#########################################
# HUEY - STEVENSON MODEL (5 parameters) #
# (Huey & Stevenson, Amer. Zool., 1979) #
#                                       #
# Parameters: a, K_l, T_min, K_u, T_max #
#########################################
fit_Huey_Stevenson_5_pars <- function(dataset)
{

	# Set the starting value of T_min to the minimum temperature in the 
	# dataset.
	T_min_start <- min(dataset$temp)

	# Set the starting value of T_max to the maximum temperature in the 
	# dataset.
	T_max_start <- max(dataset$temp)

	# Set the starting value of K_l arbitrarily to 0.1.
	K_l_start <- 0.1
	
	# Set the starting value of K_u arbitrarily to 0.7.
	K_u_start <- 0.7
	
	# Set the starting value of a arbitrarily to 1.
	a_start <- 1

	function_to_be_fitted <- function(a, K_l, T_min, K_u, T_max, temp)
	{
		
		if ( T_min >= T_max || any(temp < T_min) || any(temp > T_max) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(
				log(
					a * ( 1 - exp(-K_l * (temp - T_min) ) ) * ( 1 - exp( K_u * ( temp - T_max ) )  )
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, K_l, T_min, K_u, T_max, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,			K_l = 0.5 * K_l_start,
				T_min = 0.5 * T_min_start,	K_u = 0.5 * K_u_start,
				T_max = 0.5 * T_max_start
			),
			start_upper = c(
				a = 1.5 * a_start,			K_l = 1.5 * K_l_start,
				T_min = 1.5 * T_min_start,	K_u = 1.5 * K_u_start,
				T_max = 1.5 * T_max_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, -20, 0, 0),
			upper = c(Inf, 1, 150, 1, 150)
		)
	)
	
	return(fit)
}

#############################################
# JANISCH I MODEL (3 parameters)            #
# (Janisch, Pflüger's Arch. Physiol., 1925) #
#                                           #
# Parameters: m, a, T_pk                    #
#############################################
fit_Janisch_I_3_pars <- function(dataset)
{
	
	# Set the starting value of m arbitrarily to 3.
	m_start <- 3
	
	# Set the starting value of a arbitrarily to 0.3.
	a_start <- 0.3
	
	# Set the starting value of T_pk to the temperature at which the 
	# trait reaches its maximum value.
	T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])

	function_to_be_fitted <- function(m, a, T_pk, temp)
	{
		return(
			log(
				1 / ( (m/2) * ( a^(temp - T_pk) + a^( - (temp - T_pk) ) ) )
			)
		)
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				m, a, T_pk, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				m = 0.5 * m_start,		a = 0.5 * a_start,
				T_pk = 0.5 * T_pk_start
			),
			start_upper = c(
				m = 1.5 * m_start,		a = 1.5 * a_start,
				T_pk = 1.5 * T_pk_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0),
			upper = c(Inf, Inf, 150)
		)
	)
	
	return(fit)
}

#############################################
# JANISCH II MODEL (4 parameters)           #
# (Janisch, Pflüger's Arch. Physiol., 1925) #
#                                           #
# Parameters: m, a, b, T_pk                 #
#############################################
fit_Janisch_II_4_pars <- function(dataset)
{
	
	# Set the starting value of m arbitrarily to 3.
	m_start <- 3
	
	# Set the starting value of a arbitrarily to 0.3.
	a_start <- 0.3
	
	# Set the starting value of T_pk to the temperature at which the 
	# trait reaches its maximum value.
	T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])

	function_to_be_fitted <- function(m, a, b, T_pk, temp)
	{
		return(
			log(
				1 / ( (m/2) * ( a^(temp - T_pk) + b^( - (temp - T_pk) ) ) )
			)
		)
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				m, a, b, T_pk, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				m = 0.5 * m_start,	a = 0.5 * a_start,
				b = 0.5 * b_start,	T_pk = 0.5 * T_pk_start	
			),
			start_upper = c(
				m = 1.5 * m_start,	a = 1.5 * a_start,
				b = 1.5 * b_start,	T_pk = 1.5 * T_pk_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0, 0),
			upper = c(Inf, Inf, Inf, 150)
		)
	)
	
	return(fit)
}

############################################
# JÖHNK MODEL (5 parameters)               #
# (Jöhnk et al., Glob. Chang. Biol., 2008) #
#                                          #
# Parameters: B_pk, a, T_pk, b, d          #
############################################
fit_Johnk_5_pars <- function(dataset)
{
	
	# Set the starting value of B_pk to the maximum trait value in the 
	# dataset.
	B_pk_start <- max(dataset$trait_value)
	
	# Set the starting value of a arbitrarily to 5.
	a_start <- 5
	
	# Set the starting value of T_pk to the temperature at which the 
	# trait reaches its maximum value.
	T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
		
	# Set the starting value of b arbitrarily to 1.1.
	b_start <- 1.1
	
	# Set the starting value of d arbitrarily to 1.2.
	d_start <- 1.2

	function_to_be_fitted <- function(B_pk, a, T_pk, b, d, temp)
	{
		result <- B_pk * ( 1 + a * ( ( b ^ ( temp - T_pk ) - 1 ) - ( log(b)/log(d) ) * ( d ^ ( temp - T_pk ) - 1 ) ) )

		if ( any( result > B_pk ) || any ( result <= 0 ) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_pk, a, T_pk, b, d, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_pk = 0.5 * B_pk_start,	a = 0.5 * a_start,
				T_pk = 0.5 * T_pk_start,	b = 0.5 * b_start,
				d = 0.5 * d_start
			),
			start_upper = c(
				B_pk = 1.5 * B_pk_start,	a = 1.5 * a_start,
				T_pk = 1.5 * T_pk_start,	b = 1.5 * b_start,
				d = 1.5 * d_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0, 1, 1),
			upper = c(Inf, 100, 150, 10, 10)
		)
	)
	
	return(fit)
}

####################################################
# JOHNSON - LEWIN MODEL (4 parameters)             #
# (Johnson & Lewin, J. Cell. Comp. Physiol., 1946) #
#                                                  #
# Parameters: B_0, E, T_pk, E_D                    #
####################################################
fit_Johnson_Lewin_4_pars <- function(dataset)
{
	
	# Set the Boltzmann constant.
	k <- 8.617 * 10^-5
	
	# Set the minimum trait measurement at the rise of the TPC as the 
	# starting value for B_0.
	B_0_start <- min(dataset$trait_value[dataset$temp == min(dataset$temp)])

	# Set the starting value of T_pk to the temperature at which the 
	# trait reaches its maximum value.
	T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])

	# If the dataset has measurements before the thermal optimum, 
	# a starting value for E can be set as roughly the slope of the 
	# following regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * 1/(k * (T + 273.15))
	#
	# Otherwise, just set E to 0.6 eV.
	
	dataset_before_peak <- dataset[dataset$temp < T_pk_start,]
	if ( 
		nrow(dataset_before_peak) < 3 || 
		length(unique(dataset_before_peak$temp)) < 3 || 
		length(unique(dataset_before_peak$trait_value)) < 3 
	)
	{
		E_start <- 0.6
	} else
	{
		y_vals <- log(dataset_before_peak$trait_value)
		x_vals <- 1/(k * (dataset_before_peak$temp + 273.15))
		
		# Take the absolute value.
		E_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
		
		if ( E_start >= 10 )
		{
			E_start <- 0.6
		}
	}
	
	# If the dataset has measurements after the thermal optimum, 
	# a starting value for E_D can be set as the slope of the following 
	# regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * 1/(k * (T + 273.15))
	#
	# Otherwise, just set E_D to 3 eV.
	dataset_after_peak <- dataset[dataset$temp > T_pk_start,]
	if ( 
		nrow(dataset_after_peak) < 3 || 
		length(unique(dataset_after_peak$temp)) < 3 || 
		length(unique(dataset_after_peak$trait_value)) < 3 
	)
	{
		E_D_start <- 3
	} else
	{
		y_vals <- log(dataset_after_peak$trait_value)
		x_vals <- 1/(k * (dataset_after_peak$temp + 273.15))
		
		# Take the absolute value.
		E_D_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
		
		if ( E_D_start >= 50 )
		{
			E_D_start <- 3
		}
	}
	
	if ( E_start >= E_D_start )
	{
		E_start <- 0.9 * E_D_start
	}
	
	function_to_be_fitted <- function(B_0, E, T_pk, E_D, temp)
	{
		temp <- temp + 273.15

		# Set the Boltzmann constant.
		k <- 8.617 * 10^-5
				
		if ( E == E_D )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(
				log(
					B_0 * temp * exp((-E/k) * (1/temp)) / ( 1 + ( E/( E_D - E ) ) * exp( (E_D/k) * (1/T_pk - 1/temp) ) )
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_0, E, T_pk, E_D, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_0 = 0.5 * B_0_start,				E = 0.5 * E_start,
				T_pk = 0.5 * T_pk_start + 273.15,	E_D = 0.5 * E_D_start
			),
			start_upper = c(
				B_0 = 1.5 * B_0_start,				E = 1.5 * E_start,
				T_pk = 1.5 * T_pk_start + 273.15,	E_D = 1.5 * E_D_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 273.15, 0),
			upper = c(Inf, 10, 273.15 + 150, 50)
		)
	)
	
	return(fit)
}

#####################################
# KUMARASWAMY MODEL (5 parameters)  #
# (Tittes et al., Am. Nat., 2019)   #
#                                   #
# Parameters: a, b, T_min, T_max, d #
#####################################
fit_Kumaraswamy_5_pars <- function(dataset)
{

	# Set the starting value of a arbitrarily to 1.
	a_start <- 1
	
	# Set the starting value of b arbitrarily to 5.
	b_start <- 5
	
	# Set the starting value of d arbitrarily to 4.
	d_start <- 4

	# Set the starting value of T_min 1 degree below the lowest 
	# temperature in the dataset.
	T_min_start <- min(dataset$temp) - 1
	
	# Set the starting value of T_max 1 degree above the highest 
	# temperature in the dataset.
	T_max_start <- max(dataset$temp) + 1
		
	function_to_be_fitted <- function(a, b, T_min, T_max, d, temp)
	{
		temp_factor <- ( (temp - T_min) / (T_max - T_min) )
		
		result <- d * a * b * ( temp_factor^(a - 1) ) * ( ( 1 - temp_factor^a )^(b - 1) )
				
		if ( any(temp <= T_min) || any(temp >= T_max) || any(is.nan(result)) || any(result <= 0) )
		{
			return(rep(1e10, length(temp)))
		} else
		{	
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, b, T_min, T_max, d, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,			b = 0.5 * b_start,
				T_min = 0.5 * T_min_start,	T_max = 0.5 * T_max_start,
				d = 0.5 * d_start
			),
			start_upper = c(
				a = 1.5 * a_start,			b = 1.5 * b_start,
				T_min = 1.5 * T_min_start,	T_max = 1.5 * T_max_start,
				d = 1.5 * d_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, -20, 0, 0),
			upper = c(10, 10, 150, 150, Inf)
		)
	)
	
	return(fit)
}

############################################
# LACTIN I MODEL (3 parameters)            #
# (Lactin et al., Environ. Entomol., 1995) #
#                                          #
# Parameters: rho, T_max, DT               #
############################################
fit_Lactin_I_3_pars <- function(dataset)
{

	# Set the starting value of T_max to the maximum temperature in the 
	# dataset.
	T_max_start <- max(dataset$temp)

	# Calculate the value of the thermal optimum (T_pk).
	T_pk <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])

	# Set the starting value of DT as 0.8 times the difference between 
	# T_max_start and T_pk.
	DT_start <- 0.8 * (T_max_start - T_pk)

	# If the dataset has measurements before the thermal optimum, 
	# a starting value for rho can be set as the slope of the following 
	# regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * T
	#
	# Otherwise, just arbitrarily set rho to 0.6.
		
	dataset_before_peak <- dataset[dataset$temp < T_pk,]
	if ( 
		nrow(dataset_before_peak) < 3 || 
		length(unique(dataset_before_peak$temp)) < 3 || 
		length(unique(dataset_before_peak$trait_value)) < 3 
	)
	{
		rho_start <- 0.6
	} else
	{
		y_vals <- log(dataset_before_peak$trait_value)
		x_vals <- dataset_before_peak$temp
		
		# Take the absolute value.
		rho_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
	}
		
	function_to_be_fitted <- function(rho, T_max, DT, temp)
	{
		if ( any(temp > T_max) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(
				log(
					exp(rho * temp) - exp( rho * T_max - ( (T_max - temp) / DT ) )
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				rho, T_max, DT, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				rho = 0.5 * rho_start,	T_max = 0.5 * T_max_start,
				DT = 0.5 * DT_start
			),
			start_upper = c(
				rho = 1.5 * rho_start,	T_max = 1.5 * T_max_start,
				DT = 1.5 * DT_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0),
			upper = c(Inf, 150, 150)
		)
	)
	
	return(fit)
}

############################################
# LACTIN II MODEL (4 parameters)           #
# (Lactin et al., Environ. Entomol., 1995) #
#                                          #
# Parameters: rho, T_max, DT, lambda       #
############################################
fit_Lactin_II_4_pars <- function(dataset)
{

	# Set the starting value of lambda arbitrarily to 0.3.
	lambda_start <- 0.3

	# Set the starting value of T_max to the maximum temperature in the 
	# dataset.
	T_max_start <- max(dataset$temp)

	# Calculate the value of the thermal optimum (T_pk).
	T_pk <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])

	# Set the starting value of DT as 0.8 times the difference between 
	# T_max_start and T_pk.
	DT_start <- 0.8 * (T_max_start - T_pk)

	# If the dataset has measurements before the thermal optimum, 
	# a starting value for rho can be set as the slope of the following 
	# regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * T
	#
	# Otherwise, just arbitrarily set rho to 0.6.
		
	dataset_before_peak <- dataset[dataset$temp < T_pk,]
	if ( 
		nrow(dataset_before_peak) < 3 || 
		length(unique(dataset_before_peak$temp)) < 3 || 
		length(unique(dataset_before_peak$trait_value)) < 3 
	)
	{
		rho_start <- 0.6
	} else
	{
		y_vals <- log(dataset_before_peak$trait_value)
		x_vals <- dataset_before_peak$temp
		
		# Take the absolute value.
		rho_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
	}
		
	function_to_be_fitted <- function(rho, T_max, DT, lambda, temp)
	{
		if ( any(temp > T_max) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(
				log(
					exp(rho * temp) - exp( rho * T_max - ( (T_max - temp) / DT ) ) + lambda
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				rho, T_max, DT, lambda, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				rho = 0.5 * rho_start,	T_max = 0.5 * T_max_start,
				DT = 0.5 * DT_start,	lambda = 0.5 * lambda_start
			),
			start_upper = c(
				rho = 1.5 * rho_start,	T_max = 1.5 * T_max_start,
				DT = 1.5 * DT_start,	lambda = 1.5 * lambda_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0, -Inf),
			upper = c(Inf, 150, 150, Inf)
		)
	)
	
	return(fit)
}

##########################################
# LINEAR - LOGISTIC MODEL (5 parameters) #
# (Hui et al., Oecologia, 2019)          #
#                                        #
# Parameters: a, T_min, b, T_max, d      #
##########################################
fit_linear_logistic_5_pars <- function(dataset)
{
	
	# Set the starting value of a arbitrarily to 0.09.
	a_start <- 0.09
	
	# Set the starting value of T_min 1 degree below the lowest 
	# temperature in the dataset.
	T_min_start <- min(dataset$temp) - 1
	
	# Set the starting value of T_max 1 degree above the highest 
	# temperature in the dataset.
	T_max_start <- max(dataset$temp) + 1

	# Set the starting value of b arbitrarily to 0.85.
	b_start <- 0.85
	
	# Set the starting value of d arbitrarily to 0.73.
	d_start <- 0.73
	
	function_to_be_fitted <- function(a, T_min, b, T_max, d, temp)
	{		
		if ( any(temp <= T_min) || any(temp >= T_max) || ( T_min >= T_max ) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(
				log(
					a * (temp - T_min) * ( (1 - exp(-b * (T_max - temp))) / (1 + exp(-b * (d - temp))) )
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, T_min, b, T_max, d, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,	T_min = 0.5 * T_min_start,
				b = 0.5 * b_start,	T_max = 0.5 * T_max_start,
				d = 0.5 * d_start
			),
			start_upper = c(
				a = 1.5 * a_start,	T_min = 1.5 * T_min_start,
				b = 1.5 * b_start,	T_max = 1.5 * T_max_start,
				d = 1.5 * d_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, -20, 0, 0, 0),
			upper = c(Inf, 150, Inf, 150, Inf)
		)
	)
	
	return(fit)
}

###########################################
# LOGAN I MODEL (4 parameters)            #
# (Logan et al., Environ. Entomol., 1976) #
#                                         #
# Parameters: psi, rho, T_max, DT         #
###########################################
fit_Logan_I_4_pars <- function(dataset)
{

	# Set the minimum trait measurement as the starting value for psi.
	psi_start <- min(dataset$trait_value[dataset$temp == min(dataset$temp)])

	# Set the starting value of T_max to the maximum temperature in the 
	# dataset.
	T_max_start <- max(dataset$temp)

	# Calculate the value of the thermal optimum (T_pk).
	T_pk <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])

	# Set the starting value of DT as 0.8 times the difference between 
	# T_max_start and T_pk.
	DT_start <- 0.8 * (T_max_start - T_pk)

	# If the dataset has measurements before the thermal optimum, 
	# a starting value for rho can be set as the slope of the following 
	# regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * T
	#
	# Otherwise, just arbitrarily set rho to 0.6.
		
	dataset_before_peak <- dataset[dataset$temp < T_pk,]
	if ( 
		nrow(dataset_before_peak) < 3 || 
		length(unique(dataset_before_peak$temp)) < 3 || 
		length(unique(dataset_before_peak$trait_value)) < 3 
	)
	{
		rho_start <- 0.6
	} else
	{
		y_vals <- log(dataset_before_peak$trait_value)
		x_vals <- dataset_before_peak$temp
		
		# Take the absolute value.
		rho_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
	}
		
	function_to_be_fitted <- function(psi, rho, T_max, DT, temp)
	{
		if ( any(temp > T_max) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(
				log(
					psi * ( exp( rho * temp ) - exp( rho * T_max - ( (T_max - temp) / DT ) ) )
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				psi, rho, T_max, DT, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				psi = 0.5 * psi_start,		rho = 0.5 * rho_start,
				T_max = 0.5 * T_max_start,	DT = 0.5 * DT_start
			),
			start_upper = c(
				psi = 1.5 * psi_start,		rho = 1.5 * rho_start,
				T_max = 1.5 * T_max_start,	DT = 1.5 * DT_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0, 0),
			upper = c(Inf, Inf, 150, 150)
		)
	)
	
	return(fit)
}

###########################################
# LOGAN II MODEL (5 parameters)           #
# (Logan et al., Environ. Entomol., 1976) #
#                                         #
# Parameters: alpha, k, rho, T_max, DT    #
###########################################
fit_Logan_II_5_pars <- function(dataset)
{

	# Set the starting value of alpha arbitrarily to 0.4.
	alpha_start <- 0.4
	
	# Set the starting value of k arbitrarily to 0.3
	k_start <- 0.3

	# Set the starting value of T_max to the maximum temperature in the 
	# dataset + 10.
	T_max_start <- max(dataset$temp) + 10

	# Calculate the value of the thermal optimum (T_pk).
	T_pk <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])

	# Set the starting value of DT as 0.8 times the difference between 
	# T_max_start and T_pk.
	DT_start <- 0.8 * (T_max_start - T_pk)

	# If the dataset has measurements before the thermal optimum, 
	# a starting value for rho can be set as the slope of the following 
	# regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * T
	#
	# Otherwise, just arbitrarily set rho to 0.6.
		
	dataset_before_peak <- dataset[dataset$temp < T_pk,]
	if ( 
		nrow(dataset_before_peak) < 3 || 
		length(unique(dataset_before_peak$temp)) < 3 || 
		length(unique(dataset_before_peak$trait_value)) < 3 
	)
	{
		rho_start <- 0.6
	} else
	{
		y_vals <- log(dataset_before_peak$trait_value)
		x_vals <- dataset_before_peak$temp
		
		# Take the absolute value.
		rho_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
	}
	
	function_to_be_fitted <- function(alpha, k, rho, T_max, DT, temp)
	{
		
		if ( any(temp > T_max) || any( ( 1 / ( 1 + k * exp( -rho * ( temp ) ) ) ) <= exp( - ( ( T_max - temp ) / DT ) ) ) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(
				log(
					alpha * ( ( 1 / ( 1 + k * exp( -rho * ( temp ) ) ) ) - exp( - ( ( T_max - temp ) / DT ) ) )
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				alpha, k, rho, T_max, DT, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				alpha = 0.5 * alpha_start,	k = 0.5 * k_start,
				rho = 0.5 * rho_start,		T_max = 0.5 * T_max_start,
				DT = 0.5 * DT_start
			),
			start_upper = c(
				alpha = 1.5 * alpha_start,	k = 1.5 * k_start,
				rho = 1.5 * rho_start,		T_max = 1.5 * T_max_start,
				DT = 1.5 * DT_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0, 0, 0),
			upper = c(Inf, Inf, Inf, 150, 150)
		)
	)
	
	return(fit)
}

##############################################
# LOGAN III MODEL (4 parameters)             #
# (Hilbert & Logan, Environ. Entomol., 1983) #
#                                            #
# Parameters: psi, D, T_max, DT              #
##############################################
fit_Logan_III_4_pars <- function(dataset)
{

	# Set the starting value of psi arbitrarily to 30.
	psi_start <- 30

	# Set the starting value of T_max to the maximum temperature in the 
	# dataset.
	T_max_start <- max(dataset$temp)

	# Calculate the value of the thermal optimum (T_pk).
	T_pk <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])

	# Set the starting value of DT as 0.8 times the difference between 
	# T_max_start and T_pk.
	DT_start <- 0.8 * (T_max_start - T_pk)
	
	# Set the starting value of D as twice the difference between the 
	# maximum and minimum temperatures in the dataset.
	D_start <- 2 * (max(dataset$temp) - min(dataset$temp))
		
	function_to_be_fitted <- function(psi, D, T_max, DT, temp)
	{
		if ( any(temp > T_max) || any( ( psi * ( temp^2 / ( temp^2 + D^2 ) ) ) <= exp( - ( ( T_max - temp ) / DT ) ) ) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(
				log(
					psi * ( temp^2 / ( temp^2 + D^2 ) ) - exp( - ( (T_max - temp) / DT ) )
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				psi, D, T_max, DT, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				psi = 0.5 * psi_start,	D = 0.5 * D_start,
				T_max = 0.5 * T_max_start,	DT = 0.5 * DT_start
			),
			start_upper = c(
				psi = 1.5 * psi_start,	D = 1.5 * D_start,
				T_max = 1.5 * T_max_start,	DT = 1.5 * DT_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0, 0),
			upper = c(Inf, Inf, 150, 150)
		)
	)
	
	return(fit)
}

###############################################
# MITCHELL - ANGILLETTA MODEL (3 parameters)  #
# (Mitchell & Angilletta, Funct. Ecol., 2009) #
#                                             #
# Parameters: a, b, T_pk                      #
###############################################
fit_Mitchell_Angilletta_3_pars <- function(dataset)
{
	
	# Set the starting value of a arbitrarily to 1000.
	a_start <- 1000
	
	# Set the starting value of b to the difference between the maximum 
	# and minimum temperature in the dataset.
	b_start <- max(dataset$temp) - min(dataset$temp)
	
	# Set the starting value of T_pk to the temperature at which the 
	# trait reaches its maximum value.
	T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
		
	function_to_be_fitted <- function(a, b, T_pk, temp)
	{
		result <- ( 1 / (2 * b) ) * ( 1 + cos( ( (temp - T_pk) / b ) * pi ) ) * a
		
		# Make sure that all values are positive and that the resulting 
		# curve is not multimodal.
		if ( any(result <= 0 ) || any( ((temp - T_pk)/b) > 2 ) || any( ((temp - T_pk)/b) < -2 ) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, b, T_pk, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,		b = 0.5 * b_start,
				T_pk = 0.5 * T_pk_start
			),
			start_upper = c(
				a = 1.5 * a_start,		b = 1.5 * b_start,
				T_pk = 1.5 * T_pk_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0),
			upper = c(Inf, 150, 150)
		)
	)
	
	return(fit)
}

##############################################
# MODIFIED BILINEAR MODEL (6 parameters)     #
# (Torabi et al., Int. J. Plant Prod., 2020) #
#                                            #
# Parameters: B_pk, T_min, T_max, T_pk, a, b #
##############################################
fit_modified_Bilinear_6_pars <- function(dataset)
{

	# Set the starting value of B_pk to the maximum trait value in the 
	# dataset.
	B_pk_start <- max(dataset$trait_value)

	# Set the starting value of T_min 1 degree below the minimum 
	# temperature in the dataset.
	T_min_start <- min(dataset$temp) - 1

	# Set the starting value of T_max 1 degree above the maximum 
	# temperature in the dataset.
	T_max_start <- max(dataset$temp) + 1
	
	# Set the starting value of T_pk to the temperature at which the 
	# trait reaches its maximum value.
	T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
	
	# Set the starting value of a arbitrarily to 1.5.
	a_start <- 1.5
	
	# Set the starting value of b arbitrarily to 0.5.
	b_start <- 0.5
		
	function_to_be_fitted <- function(B_pk, T_min, T_max, T_pk, a, b, temp)
	{
		if ( T_min >= T_max || any(temp < T_min) || any(temp > T_max) || T_min >= T_pk || T_max <= T_pk )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			result <- rep(NA, length(temp))
			
			temps_below_T_pk <- which(temp <= T_pk)
			temps_above_T_pk <- which(temp > T_pk)
		
			result[temps_below_T_pk] <-  B_pk * ((temp[temps_below_T_pk] - T_min) / (T_pk - T_min))^a
			result[temps_above_T_pk] <-  B_pk * ((T_max - temp[temps_above_T_pk]) / (T_max - T_pk))^b
			
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_pk, T_min, T_max, T_pk, a, b, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_pk = 0.5 * B_pk_start,	T_min = 0.5 * T_min_start,
				T_max = 0.5 * T_max_start,	T_pk = 0.5 * T_pk_start,
				a = 0.5 * a_start,			b = 0.5 * b_start
			),
			start_upper = c(
				B_pk = 1.5 * B_pk_start,	T_min = 1.5 * T_min_start,
				T_max = 1.5 * T_max_start,	T_pk = 1.5 * T_pk_start,
				a = 1.5 * a_start,			b = 1.5 * b_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, -20, 0, 0, 0, 0),
			upper = c(Inf, 150, 150, 150, Inf, Inf)
		)
	)
	
	return(fit)
}

########################################################
# MODIFIED DENT-LIKE MODEL (7 parameters)              #
# (Torabi et al., Int. J. Plant Prod., 2020)           #
#                                                      #
# Parameters: B_pk, T_min, T_max, T_pk_l, T_pk_u, a, b #
########################################################
fit_modified_Dent_like_7_pars <- function(dataset)
{

	# Set the starting value of B_pk to the maximum trait value in the 
	# dataset.
	B_pk_start <- max(dataset$trait_value)

	# Set the starting value of T_min 1 degree below the minimum 
	# temperature in the dataset.
	T_min_start <- min(dataset$temp) - 1

	# Set the starting value of T_max 1 degree above the maximum 
	# temperature in the dataset.
	T_max_start <- max(dataset$temp) + 1
	
	# Set the starting value of T_pk_u to the temperature at which the 
	# trait reaches its maximum value.
	T_pk_u_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
	
	# Set the starting value of T_pk_l one degree below T_pk_u.
	T_pk_l_start <- T_pk_u_start - 1
	
	# Set the starting value of a arbitrarily to 1.5.
	a_start <- 1.5
	
	# Set the starting value of b arbitrarily to 0.5.
	b_start <- 0.5
		
	function_to_be_fitted <- function(B_pk, T_min, T_max, T_pk_l, T_pk_u, a, b, temp)
	{
		if ( T_min >= T_max || any(temp < T_min) || any(temp > T_max) || T_pk_l >= T_pk_u || T_min >= T_pk_l || T_max <= T_pk_u )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			result <- rep(NA, length(temp))
			
			temps_below_T_pk_l <- which(temp < T_pk_l)
			temps_between_T_pk_l_and_T_pk_u <- which(temp >= T_pk_l & temp <= T_pk_u)
			temps_above_T_pk_u <- which(temp > T_pk_u)
		
			result[temps_below_T_pk_l] <-  B_pk * ((temp[temps_below_T_pk_l] - T_min) / (T_pk_l - T_min))^a
			result[temps_between_T_pk_l_and_T_pk_u] <- B_pk
			result[temps_above_T_pk_u] <-  B_pk * ((T_max - temp[temps_above_T_pk_u]) / (T_max - T_pk_u))^b
			
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_pk, T_min, T_max, T_pk_l, T_pk_u, a, b, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_pk = 0.5 * B_pk_start,		T_min = 0.5 * T_min_start,
				T_max = 0.5 * T_max_start,		T_pk_l = 0.5 * T_pk_l_start,
				T_pk_u = 0.5 * T_pk_u_start,	a = 0.5 * a_start,
				b = 0.5 * b_start
			),
			start_upper = c(
				B_pk = 1.5 * B_pk_start,		T_min = 1.5 * T_min_start,
				T_max = 1.5 * T_max_start,		T_pk_l = 1.5 * T_pk_l_start,
				T_pk_u = 1.5 * T_pk_u_start,	a = 1.5 * a_start,
				b = 1.5 * b_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, -20, 0, 0, 0, 0, 0),
			upper = c(Inf, 150, 150, 150, 150, Inf, Inf)
		)
	)
	
	return(fit)
}

#############################################
# MODIFIED DEUTSCH MODEL (4 parameters)     #
# (Krenek et al., Eur. J. Protistol., 2011) #
#                                           #
# Parameters: T_pk, sigma, T_max, B_pk      #
#############################################
fit_modified_Deutsch_4_pars <- function(dataset)
{

	# Set the starting value of T_pk to the temperature at which the 
	# trait reaches its maximum value.
	T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
	
	# Set the starting value of T_max to the maximum temperature in the 
	# dataset.
	T_max_start <- max(dataset$temp)
	
	# Set the starting value of sigma arbitrarily to 3.
	sigma_start <- 3
	
	# Set the starting value of B_pk to the maximum trait value in the 
	# dataset.
	B_pk_start <- max(dataset$trait_value)
	
	function_to_be_fitted <- function(T_pk, sigma, T_max, B_pk, temp)
	{
		if ( T_pk >= T_max || any(temp > T_max) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			result <- rep(NA, length(temp))
		
			temps_below_T_pk <- which(temp <= T_pk)
			temps_above_T_pk <- which(temp > T_pk)
		
			result[temps_below_T_pk] <-  B_pk * exp( - ( ( temp[temps_below_T_pk] - T_pk )/( 2 * sigma ) )^2 )
			result[temps_above_T_pk] <-  B_pk - B_pk * ( (temp[temps_above_T_pk] - T_pk) / (T_pk - T_max) )^2
		
			if ( any( result > B_pk ) )
			{
				return(rep(1e10, length(temp)))
			} else
			{
				return(log(result))
			}
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				T_pk, sigma, T_max, B_pk, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				T_pk = 0.5 * T_pk_start,	sigma = 0.5 * sigma_start,
				T_max = 0.5 * T_max_start,	B_pk = 0.5 * B_pk_start
			),
			start_upper = c(
				T_pk = 1.5 * T_pk_start,	sigma = 1.5 * sigma_start,
				T_max = 1.5 * T_max_start,	B_pk = 1.5 * B_pk_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0, 0),
			upper = c(150, 150, 150, Inf)
		)
	)
	
	return(fit)
}

##########################################
# MODIFIED GAUSSIAN MODEL (4 parameters) #
# (Angilletta, J. Therm. Biol., 2006)    #
#                                        #
# Parameters: B_pk, T_pk, a, b           #
##########################################
fit_modified_Gaussian_4_pars <- function(dataset)
{
	
	# Set the starting value of B_pk to the maximum trait value in the 
	# dataset.
	B_pk_start <- max(dataset$trait_value)
	
	# Set the starting value of T_pk to the temperature at which the 
	# trait reaches its maximum value.
	T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
	
	# Set the starting value of a arbitrarily to 90.
	a_start <- 90
	
	# Set the starting value of b arbitrarily to 2.
	b_start <- 2
	
	function_to_be_fitted <- function(B_pk, T_pk, a, b, temp)
	{
		return(
			log(
				B_pk * exp( - 0.5 * ( abs( temp - T_pk ) / a )^b )
			)
		)
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_pk, T_pk, a, b, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_pk = 0.5 * B_pk_start,	T_pk = 0.5 * T_pk_start,
				a = 0.5 * a_start,			b = 0.5 * b_start
			),
			start_upper = c(
				B_pk = 1.5 * B_pk_start,	T_pk = 1.5 * T_pk_start,
				a = 1.5 * a_start,			b = 1.5 * b_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0, 0),
			upper = c(Inf, 150, Inf, Inf)
		)
	)
	
	return(fit)
}

#########################################
# NEWBERY MODEL (4 parameters)          #
# (Newbery et al., Plant Pathol., 2020) #
#                                       #
# Parameters: a, b, d, e                #
#########################################
fit_Newbery_4_pars <- function(dataset)
{

	# Set the starting value of a arbitrarily to 1.
	a_start <- 1
	
	# Set the starting value of b arbitrarily to 0.4.
	b_start <- 0.4
	
	# Set the starting value of d arbitrarily to 1.
	d_start <- 1
	
	# Set the starting value of e arbitrarily to 0.0012.
	e_start <- 0.0012
	
	function_to_be_fitted <- function(a, b, d, e, temp)
	{
		result <- a + b * temp + d * (1 - exp(e * temp^2))
		
		if ( any(is.na(result)) || any(result <= 0) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, b, d, e, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,	b = 0.5 * b_start,
				d = 0.5 * d_start,	e = 0.5 * e_start
			),
			start_upper = c(
				a = 1.5 * a_start,	b = 1.5 * b_start,
				d = 1.5 * d_start,	e = 1.5 * e_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(-Inf, -Inf, -Inf, -Inf),
			upper = c(Inf, Inf, Inf, Inf)
		)
	)
	
	return(fit)
}

#######################################
# O'NEILL MODEL (4 parameters)        #
# (O'Neill et al., 1972)              #
#                                     #
# Parameters: Q_10, T_max, T_pk, B_pk #
#######################################
fit_ONeill_4_pars <- function(dataset)
{

	# Set the starting value of Q10 arbitrarily to 2.
	Q_10_start <- 2
	
	# Set the starting value of T_max to the highest temperature in the 
	# dataset.
	T_max_start <- max(dataset$temp)

	# Set the starting value of T_pk to the temperature at which the 
	# trait reaches its maximum value.
	T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
	
	# Set the starting value of B_pk to the maximum trait value in the 
	# dataset.
	B_pk_start <- max(dataset$trait_value)
	
	function_to_be_fitted <- function(Q_10, T_max, T_pk, B_pk, temp)
	{
		if ( (T_max <= T_pk) || any(temp > T_max) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			W <- (Q_10 - 1) * (T_max - T_pk)
			X <- ( W^2 * ( 1 + (1 + (40/W))^(1/2) )^2 ) / 400
			V <- (T_max - temp) / (T_max - T_pk)
			
			result <- B_pk * V^X * exp(X * (1 - V))
			
			if ( any( result > B_pk ) )
			{
				return(rep(1e10, length(temp)))
			} else
			{
				return(log(result))
			}
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				Q_10, T_max, T_pk, B_pk, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				Q_10 = 0.5 * Q_10_start,	T_max = 0.5 * T_max_start,
				T_pk = 0.5 * T_pk_start,	B_pk = 0.5 * B_pk_start
			),
			start_upper = c(
				Q_10 = 1.5 * Q_10_start,	T_max = 1.5 * T_max_start,
				T_pk = 1.5 * T_pk_start,	B_pk = 1.5 * B_pk_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0, 0),
			upper = c(Inf, 150, 150, Inf)
		)
	)
	
	return(fit)
}

###########################################
# RATKOWSKY MODEL (4 parameters)          #
# (Ratkowsky et al., J. Bacteriol., 1983) #
#                                         #
# Parameters: a, T_min, b, T_max          #
###########################################
fit_Ratkowsky_4_pars <- function(dataset)
{

	# Set the starting value of T_min to the lowest temperature in the 
	# dataset.
	T_min_start <- min(dataset$temp) + 273.15
	
	# Set the starting value of T_max to the highest temperature in the 
	# dataset.
	T_max_start <- max(dataset$temp) + 273.15

	T_pk <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
	
	# If the dataset has measurements before the thermal optimum, 
	# a starting value for a can be set as the slope of the following 
	# regression for that subset of the data:
	#
	# sqrt(trait_value) ~ intercept + slope * (T + 273.15)
	#
	# Otherwise, just arbitrarily set a to 0.6.
	
	dataset_before_peak <- dataset[dataset$temp < T_pk,]
	if ( 
		nrow(dataset_before_peak) < 3 || 
		length(unique(dataset_before_peak$temp)) < 3 || 
		length(unique(dataset_before_peak$trait_value)) < 3 
	)
	{
		a_start <- 0.6
	} else
	{
		y_vals <- sqrt(dataset_before_peak$trait_value)
		x_vals <- dataset_before_peak$temp + 273.15
		
		# Take the absolute value.
		a_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
	}
	
	# If the dataset has measurements after the thermal optimum, 
	# a starting value for b can be set as roughly the slope of 
	# the following regression for that subset of the data:
	#
	# log(sqrt(trait_value)) ~ intercept + slope * (T + 273.15)
	#
	# Otherwise, just arbitrarily set b to 3.
	dataset_after_peak <- dataset[dataset$temp > T_pk,]
	if ( 
		nrow(dataset_after_peak) < 3 || 
		length(unique(dataset_after_peak$temp)) < 3 || 
		length(unique(dataset_after_peak$trait_value)) < 3 
	)
	{
		b_start <- 3
	} else
	{
		y_vals <- log(sqrt(dataset_after_peak$trait_value))
		x_vals <- dataset_after_peak$temp + 273.15
		
		# Take the absolute value.
		b_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
	}
	
	function_to_be_fitted <- function(a, T_min, b, T_max, temp)
	{
		temp <- temp + 273.15
		
		if ( any(temp <= T_min) || any(temp >= T_max) || ( T_min >= T_max ) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(
				log(
					(a * (temp - T_min)*(1 - exp(b * (temp - T_max))))^2
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, T_min, b, T_max, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,	T_min = 0.5 * T_min_start,
				b = 0.5 * b_start,	T_max = 0.5 * T_max_start
			),
			start_upper = c(
				a = 1.5 * a_start,	T_min = 1.5 * T_min_start,
				b = 1.5 * b_start,	T_max = 1.5 * T_max_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 273.15 - 20, 0, 273.15),
			upper = c(Inf, 273.15 + 150, Inf, 273.15 + 150)
		)
	)
	
	return(fit)
}

#######################################################
# REGNIERE MODEL (6 parameters)                       #
# (Regniere et al., J. Insect Physiol., 2012)         #
#                                                     #
# Parameters: psi, rho, T_min, T_max, DT_low, DT_high #
#######################################################
fit_Regniere_6_pars <- function(dataset)
{

	# Set the minimum trait measurement as the starting value for psi.
	psi_start <- min(dataset$trait_value[dataset$temp == min(dataset$temp)])

	# Set the starting value of T_min 1 degree below the minimum 
	# temperature in the dataset.
	T_min_start <- min(dataset$temp) - 1
	
	# Set the starting value of T_max 1 degree above the maximum 
	# temperature in the dataset.
	T_max_start <- max(dataset$temp) + 1

	# Calculate the value of the thermal optimum (T_pk).
	T_pk <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])

	# Set the starting value of DT_high as 0.8 times the difference 
	# between T_max_start and T_pk.
	DT_high_start <- 0.8 * (T_max_start - T_pk)
	
	# Set the starting value of DT_low arbitrarily to 0.2.
	DT_low_start <- 0.2

	# If the dataset has measurements before the thermal optimum, 
	# a starting value for rho can be set as the slope of the following 
	# regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * T
	#
	# Otherwise, just arbitrarily set rho to 0.6.
		
	dataset_before_peak <- dataset[dataset$temp < T_pk,]
	if ( 
		nrow(dataset_before_peak) < 3 || 
		length(unique(dataset_before_peak$temp)) < 3 || 
		length(unique(dataset_before_peak$trait_value)) < 3 
	)
	{
		rho_start <- 0.6
	} else
	{
		y_vals <- log(dataset_before_peak$trait_value)
		x_vals <- dataset_before_peak$temp
		
		# Take the absolute value.
		rho_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
	}
		
	function_to_be_fitted <- function(psi, rho, T_min, T_max, DT_low, DT_high, temp)
	{
		result <- psi * ( exp(rho * (temp - T_min)) - ( (T_max - temp) / (T_max - T_min) ) * exp( -rho * (temp - T_min)/DT_low ) - ( (temp - T_min) / (T_max - T_min) ) * exp( rho * (T_max - T_min) - ( (T_max - temp) / DT_high ) ) ) 
		T_pk <- temp[result == max(result)]

		if ( any(temp >= T_max) || any(temp <= T_min) || T_min >= T_max || DT_low > (T_pk - T_min) || DT_high > (T_max - T_pk) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				psi, rho, T_min, T_max, DT_low, DT_high, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				psi = 0.5 * psi_start,			rho = 0.5 * rho_start,
				T_min = 0.5 * T_min_start,		T_max = 0.5 * T_max_start,
				DT_low = 0.5 * DT_low_start,	DT_high = 0.5 * DT_high_start
			),
			start_upper = c(
				psi = 1.5 * psi_start,			rho = 1.5 * rho_start,
				T_min = 1.5 * T_min_start,		T_max = 1.5 * T_max_start,
				DT_low = 1.5 * DT_low_start,	DT_high = 1.5 * DT_high_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, -20, 0, 0, 0),
			upper = c(Inf, Inf, 150, 150, 150, 150)
		)
	)
	
	return(fit)
}

############################################################################
# REZENDE - BOZINOVIC MODEL (4 parameters)                                 #
# (Rezende & Bozinovic, Philos. Trans. R. Soc. Lond., B, Biol. Sci., 2019) #
#                                                                          #
# Parameters: B_0, Q_10, T_th, d                                           #
############################################################################
fit_Rezende_Bozinovic_4_pars <- function(dataset)
{

	# Set the minimum trait measurement at the rise of the TPC as the 
	# starting value for B_0.
	B_0_start <- min(dataset$trait_value[dataset$temp == min(dataset$temp)])
	
	# Set the starting value of Q10 arbitrarily to 2.
	Q_10_start <- 2
	
	# Set the starting value of d arbitrarily to 0.003.
	d_start <- 0.003
	
	# Get a starting value for T_th based on the relationship among 
	# T_th, T_pk, Q_10, and d.
	T_pk <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])

	T_th_start <- T_pk + 10/log(Q_10_start) - sqrt((1/d_start) + (10/log(Q_10_start))^2)
	
	function_to_be_fitted <- function(B_0, Q_10, T_th, d, temp)
	{
		result <- rep(NA, length(temp))
		
		temps_below_T_th <- which(temp <= T_th)
		temps_above_T_th <- which(temp > T_th)
		
		result[temps_below_T_th] <-  B_0 * exp(temp[temps_below_T_th] * log(Q_10)/10)
		result[temps_above_T_th] <-  B_0 * exp(temp[temps_above_T_th] * log(Q_10)/10) * ( 1 - d * (temp[temps_above_T_th] - T_th)^2 )
		
		return(log(result))
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_0, Q_10, T_th, d, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_0 = 0.5 * B_0_start,		Q_10 = 0.5 * Q_10_start,
				T_th = 0.5 * T_th_start,	d = 0.5 * d_start
			),
			start_upper = c(
				B_0 = 1.5 * B_0_start,		Q_10 = 1.5 * Q_10_start,
				T_th = 1.5 * T_th_start,	d = 1.5 * d_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0, 0),
			upper = c(Inf, Inf, 150, Inf)
		)
	)
	
	return(fit)
}

##############################################
# RICE CLOCK MODEL (6 parameters)            #
# (Gao et al., Agric. For. Meteorol., 1992)  #
#                                            #
# Parameters: B_pk, a, b, T_min, T_max, T_pk #
##############################################
fit_Rice_Clock_6_pars <- function(dataset)
{

	# Set the starting value of B_pk to the maximum trait value in the 
	# dataset.
	B_pk_start <- max(dataset$trait_value)

	# Set the starting value of a arbitrarily to 1.
	a_start <- 1
	
	# Set the starting value of b arbitrarily to 0.5.
	b_start <- 0.5
	
	# Set the starting value of T_min 1 degree below the minimum 
	# temperature in the dataset.
	T_min_start <- min(dataset$temp) - 1

	# Set the starting value of T_max 1 degree above the maximum 
	# temperature in the dataset.
	T_max_start <- max(dataset$temp) + 1
	
	# Set the starting value of T_pk to the temperature at which the 
	# trait reaches its maximum value.
	T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
		
	function_to_be_fitted <- function(B_pk, a, b, T_min, T_max, T_pk, temp)
	{
		if ( T_min >= T_max || any(temp <= T_min) || any(temp >= T_max) || T_min >= T_pk || T_max <= T_pk )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			second_part <- (( (temp - T_min) / (T_pk - T_min) )^a) * (( (T_max - temp) / (T_max - T_pk) )^b)
			
			second_part[second_part > 1] <- 1
			
			return(
				log(
					B_pk * second_part
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_pk, a, b, T_min, T_max, T_pk, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_pk = 0.5 * B_pk_start,	a = 0.5 * a_start,
				b = 0.5 * b_start,			T_min = 0.5 * T_min_start,
				T_max = 0.5 * T_max_start,	T_pk = 0.5 * T_pk_start
			),
			start_upper = c(
				B_pk = 1.5 * B_pk_start,	a = 1.5 * a_start,
				b = 1.5 * b_start,			T_min = 1.5 * T_min_start,
				T_max = 1.5 * T_max_start,	T_pk = 1.5 * T_pk_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0, -20, 0, 0),
			upper = c(Inf, Inf, Inf, 150, 150, 150)
		)
	)
	
	return(fit)
}

################################
# RITCHIE MODEL (4 parameters) #
# (Ritchie, Sci. Rep., 2018)   #
#                              #
# Parameters: d_0, E_D, DE, w  #
################################
fit_Ritchie_4_pars <- function(dataset)
{

	# Set the minimum trait measurement at the rise of the TPC as the 
	# starting value for d_0.
	d_0_start <- min(dataset$trait_value[dataset$temp == min(dataset$temp)])

	# Set the starting value of E_D arbitrarily to 20.
	E_D_start <- 20
	
	# Set the starting value of DE arbitrarily to 30.
	DE_start <- 30
	
	# Set the starting value of w arbitrarily to 11.
	w_start <- 11

	function_to_be_fitted <- function(d_0, E_D, DE, w, temp)
	{
		temp <- temp + 273.15
		
		R <- 8.318 * 1e-3
		
		if ( any( ( DE/(R * temp) ) <= w ) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(
				log(
					R * d_0 * exp( (-E_D) / (R * temp) ) * ((DE/(R * temp)) - w)
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				d_0, E_D, DE, w, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				d_0 = 0.5 * d_0_start,	E_D = 0.5 * E_D_start,
				DE = 0.5 * DE_start,	w = 0.5 * w_start
			),
			start_upper = c(
				d_0 = 1.5 * d_0_start,	E_D = 1.5 * E_D_start,
				DE = 1.5 * DE_start,	w = 1.5 * w_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0, -Inf),
			upper = c(Inf, Inf, Inf, Inf)
		)
	)
	
	return(fit)
}

#############################################
# ROSS - RATKOWSKY MODEL (5 parameters)     #
# (Ratkowsky et al., J. Theor. Biol., 2005) #
#                                           #
# Parameters: a, DH_A, n, DH, DC_p          #
#############################################
fit_Ross_Ratkowsky_5_pars <- function(dataset)
{
	
	# Set the universal gas constant.
	R <- 8.314
	
	# Set the starting value of a arbitrarily to 1.
	a_start <- 1
	
	# Set the starting value of DH_A arbitrarily to 50000.
	DH_A_start <- 50000
	
	# Set the starting value of n arbitrarily to 300.
	n_start <- 300
	
	# Set the starting value of DH arbitrarily to 5000.
	DH_start <- 5000
	
	# Set the starting value of DC_p arbitrarily to 60.
	DC_p_start <- 60
	
	function_to_be_fitted <- function(a, DH_A, n, DH, DC_p, temp)
	{
		temp <- temp + 273.15

		# Set the universal gas constant.
		R <- 8.314
		
		D <- 1 + exp( -n * ( (DH - temp * 18.1 + DC_p * ( (temp - 373.6) - temp * log(temp/385.2) ) ) / (R * temp) ) )
		result <- ( a * temp * exp((-DH_A/(R * temp)) ) ) / D

		if ( any(is.nan(result)) || any(result <= 0 ) || a == 0 )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, DH_A, n, DH, DC_p, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,		DH_A = 0.5 * DH_A_start,
				n = 0.5 * n_start,		DH = 0.5 * DH_start,
				DC_p = 0.5 * DC_p_start
			),
			start_upper = c(
				a = 1.5 * a_start,		DH_A = 1.5 * DH_A_start,
				n = 1.5 * n_start,		DH = 1.5 * DH_start,
				DC_p = 1.5 * DC_p_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 10000, 10, 3000, 37),
			upper = c(Inf, 120000, 20000, 7000, 118)
		)
	)
	
	return(fit)
}

####################################
# RUIZ MODEL (4 parameters)        #
# (Ruiz et al., bioRxiv, 2019)     #
#                                  #
# Parameters: B_0, DB_max, a, T_pk #
####################################
fit_Ruiz_4_pars <- function(dataset)
{

	# Set the minimum trait measurement at the rise of the TPC as the 
	# starting value for B_0.
	B_0_start <- min(dataset$trait_value[dataset$temp == min(dataset$temp)])
	
	# Set the difference between the minimum and the maximum trait values 
	# as the starting value for DB_max.
	DB_max_start <- max(dataset$trait_value) - B_0_start
	
	# Set the starting value of T_pk to the temperature at which the 
	# trait reaches its maximum value.
	T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
		
	# Set the starting value of a arbitrarily to 0.25.
	a_start <- 0.25

	function_to_be_fitted <- function(B_0, DB_max, a, T_pk, temp)
	{
		return(
			log(
				B_0 + DB_max * exp(-a * (temp - T_pk)^2)
			)
		)
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_0, DB_max, a, T_pk, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_0 = 0.5 * B_0_start,	DB_max = 0.5 * DB_max_start,
				a = 0.5 * a_start,		T_pk = 0.5 * T_pk_start
			),
			start_upper = c(
				B_0 = 1.5 * B_0_start,	DB_max = 1.5 * DB_max_start,
				a = 1.5 * a_start,		T_pk = 1.5 * T_pk_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0, 0),
			upper = c(Inf, Inf, Inf, 150)
		)
	)
	
	return(fit)
}

###################################################
# SHARPE - SCHOOLFIELD MODEL (6 parameters)       #
# (Schoolfield et al., J. Theor. Biol., 1981)     #
#                                                 #
# Parameters: B_0, DH_A, DH_L, T_L50, DH_H, T_H50 #
###################################################
fit_Sharpe_Schoolfield_6_pars <- function(dataset)
{
	
	# Set the universal gas constant.
	R <- 1.987
	
	# Set the minimum trait measurement at the rise of the TPC as the 
	# starting value for B_0.
	B_0_start <- min(dataset$trait_value[dataset$temp == min(dataset$temp)])
	
	# If the dataset has measurements before the thermal optimum, 
	# a starting value for DH_A can be set as the slope of the following 
	# regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * 1/(R * (T + 273.15))
	#
	# Otherwise, just set DH_A to 15000.
	T_pk <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
	
	dataset_before_peak <- dataset[dataset$temp < T_pk,]
	if ( 
		nrow(dataset_before_peak) < 3 || 
		length(unique(dataset_before_peak$temp)) < 3 || 
		length(unique(dataset_before_peak$trait_value)) < 3 
	)
	{
		DH_A_start <- 15000
	} else
	{
		y_vals <- log(dataset_before_peak$trait_value)
		x_vals <- 1/(R * (dataset_before_peak$temp + 273.15))
		
		# Take the absolute value.
		DH_A_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
	}
	
	# If the dataset has measurements after the thermal optimum, 
	# a starting value for DH_H can be set as roughly the slope of the 
	# following regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * 1/(R * (T + 273.15))
	#
	# Otherwise, just set DH_H to 50000.
	dataset_after_peak <- dataset[dataset$temp > T_pk,]
	if ( 
		nrow(dataset_after_peak) < 3 || 
		length(unique(dataset_after_peak$temp)) < 3 || 
		length(unique(dataset_after_peak$trait_value)) < 3 
	)
	{
		DH_H_start <- 50000
	} else
	{
		y_vals <- log(dataset_after_peak$trait_value)
		x_vals <- 1/(R * (dataset_after_peak$temp + 273.15))
		
		# Take the absolute value.
		DH_H_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
	}
	
	# Set the starting value of DH_L arbitrarily to 70000.
	DH_L_start <- 70000
	
	# T_H50 should be close to the thermal optimum (T_pk), so set the 
	# starting value of T_H50 = T_pk.
	T_H50_start <- T_pk + 273.15
	
	# Set the starting value of T_L50 arbitrarily to the minimum 
	# temperature in the dataset.
	T_L50_start <- min(dataset$temp) + 273.15
	
	function_to_be_fitted <- function(B_0, DH_A, DH_L, T_L50, DH_H, T_H50, temp)
	{
		temp <- temp + 273.15

		# Set the universal gas constant.
		R <- 1.987
		
		# Set the reference temperature to 0°C.
		T_ref <- 273.15
			
		result <- (B_0 * (temp/T_ref) * exp( (DH_A/R) * ((1/T_ref) - (1/temp)) ) ) / ( 1 + exp( (-DH_L/R) * ((1/T_L50) - (1/temp)) ) + exp( (DH_H/R) * ((1/T_H50) - (1/temp)) ) )
				
		if ( T_L50 >= T_H50 || is.na(result) || any(result <= 0) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_0, DH_A, DH_L, T_L50, DH_H, T_H50, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_0 = 0.5 * B_0_start,		DH_A = 0.5 * DH_A_start,
				DH_L = 0.5 * DH_L_start,	T_L50 = 0.5 * T_L50_start,
				DH_H = 0.5 * DH_H_start,	T_H50 = 0.5 * T_H50_start
			),
			start_upper = c(
				B_0 = 1.5 * B_0_start,		DH_A = 1.5 * DH_A_start,
				DH_L = 1.5 * DH_L_start,	T_L50 = 1.5 * T_L50_start,
				DH_H = 1.5 * DH_H_start,	T_H50 = 1.5 * T_H50_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0, -20 + 273.15, 0, 0 + 273.15),
			upper = c(Inf, 300000, 300000, 150 + 273.15, 300000, 150 + 273.15)
		)
	)
	
	return(fit)
}

#######################################################
# SIMPLIFIED ASBURY - ANGILLETTA MODEL (4 parameters) #
# (Asbury & Angilletta, Am. Nat., 2010)               #
#                                                     #
# Parameters: alpha, delta, epsilon, zeta             #
#######################################################
fit_simplified_Asbury_Angilletta_4_pars <- function(dataset)
{
	
	# Set the starting value of alpha arbitrarily to 1°C (274.15 K).
	alpha_start <- 274.15
	
	# Set the starting value of delta arbitrarily to 0.1.
	delta_start <- 0.1
	
	# Set the starting value of epsilon arbitrarily to 0.7.
	epsilon_start <- 0.7
	
	# Set the starting value of zeta arbitrarily to 100.
	zeta_start <- 100
		
	function_to_be_fitted <- function(alpha, delta, epsilon, zeta, temp)
	{
		temp <- temp + 273.15
		
		if ( any( (temp - alpha) > zeta ) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(
				log(
					( ( (temp-alpha)/zeta ) ^ (epsilon/delta-1) * (1 - (temp-alpha) / zeta) ^ ( (1-epsilon) / delta - 1) ) / beta( epsilon / delta, (1 - epsilon) / delta) 
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				alpha, delta, epsilon, zeta, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				alpha = 0.5 * alpha_start,		delta = 0.5 * delta_start,
				epsilon = 0.5 * epsilon_start,	zeta = 0.5 * zeta_start
			),
			start_upper = c(
				alpha = 1.5 * alpha_start,		delta = 1.5 * delta_start,
				epsilon = 1.5 * epsilon_start,	zeta = 1.5 * zeta_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0, 0),
			upper = c(Inf, Inf, 1, Inf)
		)
	)
	
	return(fit)
}

##########################################################
# SIMPLIFIED BETA TYPE MODEL (3 parameters)              #
# (Damos & Savopoulou-Soultani, J. Econ. Entomol., 2008) #
#                                                        #
# Parameters: rho, alpha, beta                           #
##########################################################
fit_simplified_beta_type_3_pars <- function(dataset)
{

	# Set the starting value of rho arbitrarily to 0.0012.
	rho_start <- 0.0012

	# Set the starting value of alpha arbitrarily to 4.
	alpha_start <- 4

	# Set the starting value of beta arbitrarily to 5.
	beta_start <- 5
		
	function_to_be_fitted <- function(rho, alpha, beta, temp)
	{		
		return(
			log(
				rho * ( alpha - (temp/10) ) * ( temp/10 ) ^ beta
			)
		)
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				rho, alpha, beta, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				rho = 0.5 * rho_start,	alpha = 0.5 * alpha_start,
				beta = 0.5 * beta_start
			),
			start_upper = c(
				rho = 1.5 * rho_start,	alpha = 1.5 * alpha_start,
				beta = 1.5 * beta_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, -Inf),
			upper = c(Inf, Inf, Inf)
		)
	)
	
	return(fit)
}

############################################
# SIMPLIFIED BRIERE I MODEL (3 parameters) #
# (Briere et al., Environ. Entomol., 1999) #
#                                          #
# Parameters: a, T_min, T_max              #
############################################
fit_simplified_Briere_I_3_pars <- function(dataset)
{
	# Set the starting value of a arbitrarily to 1.
	a_start <- 1
	
	# Set the starting value of T_min 1 degree below the minimum 
	# temperature in the dataset.
	T_min_start <- min(dataset$temp) - 1

	# Set the starting value of T_max 1 degree above the maximum 
	# temperature in the dataset.
	T_max_start <- max(dataset$temp) + 1
	
	function_to_be_fitted <- function(a, T_min, T_max, temp)
	{
		if ( T_min >= T_max || any(temp < T_min) || any(temp > T_max) )
		{
			return(rep(1e10, length(temp)))
		} else
		{	
			return(
				log(
					a * ( temp - T_min ) * sqrt( T_max - temp )
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, T_min, T_max, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,			T_min = 0.5 * T_min_start,
				T_max = 0.5 * T_max_start
			),
			start_upper = c(
				a = 1.5 * a_start,			T_min = 1.5 * T_min_start,
				T_max = 1.5 * T_max_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, -20, 0),
			upper = c(Inf, 150, 150)
		)
	)
	
	return(fit)
}

#############################################
# SIMPLIFIED BRIERE II MODEL (4 parameters) #
# (Briere et al., Environ. Entomol., 1999)  #
#                                           #
# Parameters: a, T_min, T_max, m            #
#############################################
fit_simplified_Briere_II_4_pars <- function(dataset)
{

	# Set the starting value of a arbitrarily to 1.
	a_start <- 1
	
	# Set the starting value of T_min 1 degree below the minimum 
	# temperature in the dataset.
	T_min_start <- min(dataset$temp) - 1

	# Set the starting value of T_max 1 degree above the maximum 
	# temperature in the dataset.
	T_max_start <- max(dataset$temp) + 1
	
	# Set the starting value of m arbitrarily to 2.
	m_start <- 2
	
	function_to_be_fitted <- function(a, T_min, T_max, m, temp)
	{
		if ( T_min >= T_max || any(temp < T_min) || any(temp > T_max) )
		{
			return(rep(1e10, length(temp)))
		} else
		{		
			return(
				log(
					a * ( temp - T_min ) * ( T_max - temp )^(1/m)
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, T_min, T_max, m, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,			T_min = 0.5 * T_min_start,
				T_max = 0.5 * T_max_start,	m = 0.5 * m_start
			),
			start_upper = c(
				a = 1.5 * a_start,			T_min = 1.5 * T_min_start,
				T_max = 1.5 * T_max_start,	m = 1.5 * m_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, -20, 0, 0),
			upper = c(Inf, 150, 150, Inf)
		)
	)
	
	return(fit)
}

###################################################
# SIMPLIFIED EXTENDED BRIERE MODEL (5 parameters) #
# (Cruz-Loya et al., bioRxiv, 2020)               #
#                                                 #
# Parameters: a, T_min, T_max, b, d               #
###################################################
fit_simplified_extended_Briere_5_pars <- function(dataset)
{

	# Set the starting value of a arbitrarily to 1.
	a_start <- 1
	
	# Set the starting value of T_min 1 degree below the minimum 
	# temperature in the dataset.
	T_min_start <- min(dataset$temp) - 1

	# Set the starting value of T_max 1 degree above the maximum 
	# temperature in the dataset.
	T_max_start <- max(dataset$temp) + 1
	
	# Set the starting value of b arbitrarily to 1.
	b_start <- 1
	
	# Set the starting value of d arbitrarily to 1.5
	d_start <- 1.5
	
	function_to_be_fitted <- function(a, T_min, T_max, b, d, temp)
	{
		if ( T_min >= T_max || any(temp < T_min) || any(temp > T_max) )
		{
			return(rep(1e10, length(temp)))
		} else
		{	
			return(
				log(
					a * ((temp - T_min)^b) * ((T_max - temp)^d)
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, T_min, T_max, b, d, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,			T_min = 0.5 * T_min_start,
				T_max = 0.5 * T_max_start,	b = 0.5 * b_start,
				d = 0.5 * d_start
			),
			start_upper = c(
				a = 1.5 * a_start,			T_min = 1.5 * T_min_start,
				T_max = 1.5 * T_max_start,	b = 1.5 * b_start,
				d = 1.5 * d_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, -20, 0, 0, 0),
			upper = c(Inf, 150, 150, Inf, Inf)
		)
	)
	
	return(fit)
}

############################################################
# SIMPLIFIED EXTENDED JOHNSON - LEWIN MODEL (5 parameters) #
# (Guenther, Ecol. Appl., 1997)                            #
#                                                          #
# Parameters: B_0, E, T_pk, E_D, a                         #
############################################################
fit_simplified_extended_Johnson_Lewin_5_pars <- function(dataset)
{
	
	# Set the Boltzmann constant.
	k <- 8.617 * 10^-5
	
	# Set the minimum trait measurement at the rise of the TPC as the 
	# starting value for B_0.
	B_0_start <- min(dataset$trait_value[dataset$temp == min(dataset$temp)])

	# Set the starting value of T_pk to the temperature at which the 
	# trait reaches its maximum value.
	T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])

	# If the dataset has measurements before the thermal optimum, 
	# a starting value for E can be set as roughly the slope of the 
	# following regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * 1/(k * (T + 273.15))
	#
	# Otherwise, just set E to 0.6 eV.
	
	dataset_before_peak <- dataset[dataset$temp < T_pk_start,]
	if ( 
		nrow(dataset_before_peak) < 3 || 
		length(unique(dataset_before_peak$temp)) < 3 || 
		length(unique(dataset_before_peak$trait_value)) < 3 
	)
	{
		E_start <- 0.6
	} else
	{
		y_vals <- log(dataset_before_peak$trait_value)
		x_vals <- 1/(k * (dataset_before_peak$temp + 273.15))
		
		# Take the absolute value.
		E_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
		
		if ( E_start >= 10 )
		{
			E_start <- 0.6
		}
	}
	
	# If the dataset has measurements after the thermal optimum, 
	# a starting value for E_D can be set as the slope of the following 
	# regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * 1/(k * (T + 273.15))
	#
	# Otherwise, just set E_D to 3 eV.
	dataset_after_peak <- dataset[dataset$temp > T_pk_start,]
	if ( 
		nrow(dataset_after_peak) < 3 || 
		length(unique(dataset_after_peak$temp)) < 3 || 
		length(unique(dataset_after_peak$trait_value)) < 3 
	)
	{
		E_D_start <- 3
	} else
	{
		y_vals <- log(dataset_after_peak$trait_value)
		x_vals <- 1/(k * (dataset_after_peak$temp + 273.15))
		
		# Take the absolute value.
		E_D_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
		
		if ( E_D_start >= 50 )
		{
			E_D_start <- 3
		}
	}
	
	if ( E_start >= E_D_start )
	{
		E_start <- 0.9 * E_D_start
	}
	
	# Set the starting value of a arbitrarily to 1.
	a_start <- 1
	
	function_to_be_fitted <- function(B_0, E, T_pk, E_D, a, temp)
	{
		temp <- temp + 273.15

		# Set the Boltzmann constant.
		k <- 8.617 * 10^-5
				
		if ( E == E_D )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(
				log(
					B_0 * exp((-E/k) * (1/temp)) / ( a + ( E/( E_D - E ) ) * exp( (E_D/k) * (1/T_pk - 1/temp) ) )
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_0, E, T_pk, E_D, a, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_0 = 0.5 * B_0_start,				E = 0.5 * E_start,
				T_pk = 0.5 * T_pk_start + 273.15,	E_D = 0.5 * E_D_start,
				a = 0.5 * a_start
			),
			start_upper = c(
				B_0 = 1.5 * B_0_start,				E = 1.5 * E_start,
				T_pk = 1.5 * T_pk_start + 273.15,	E_D = 1.5 * E_D_start,
				a = 1.5 * a_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 273.15, 0, 0),
			upper = c(Inf, 10, 273.15 + 150, 50, Inf)
		)
	)
	
	return(fit)
}

#################################################################
# SIMPLIFIED EXTENDED SHARPE - SCHOOLFIELD MODEL (7 parameters) #
# (Guenther, Ecol. Appl., 1997)                                 #
#                                                               #
# Parameters: B_0, DH_A, DH_L, T_L50, DH_H, T_H50, a            #
#################################################################
fit_simplified_extended_Sharpe_Schoolfield_7_pars <- function(dataset)
{
	
	# Set the universal gas constant.
	R <- 1.987
	
	# Set the minimum trait measurement at the rise of the TPC as the 
	# starting value for B_0.
	B_0_start <- min(dataset$trait_value[dataset$temp == min(dataset$temp)])
	
	# If the dataset has measurements before the thermal optimum, 
	# a starting value for DH_A can be set as the slope of the following 
	# regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * 1/(R * (T + 273.15))
	#
	# Otherwise, just set DH_A to 15000.
	T_pk <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
	
	dataset_before_peak <- dataset[dataset$temp < T_pk,]
	if ( 
		nrow(dataset_before_peak) < 3 || 
		length(unique(dataset_before_peak$temp)) < 3 || 
		length(unique(dataset_before_peak$trait_value)) < 3 
	)
	{
		DH_A_start <- 15000
	} else
	{
		y_vals <- log(dataset_before_peak$trait_value)
		x_vals <- 1/(R * (dataset_before_peak$temp + 273.15))
		
		# Take the absolute value.
		DH_A_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
	}
	
	# If the dataset has measurements after the thermal optimum, 
	# a starting value for DH_H can be set as roughly the slope of the 
	# following regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * 1/(R * (T + 273.15))
	#
	# Otherwise, just set DH_H to 50000.
	dataset_after_peak <- dataset[dataset$temp > T_pk,]
	if ( 
		nrow(dataset_after_peak) < 3 || 
		length(unique(dataset_after_peak$temp)) < 3 || 
		length(unique(dataset_after_peak$trait_value)) < 3 
	)
	{
		DH_H_start <- 50000
	} else
	{
		y_vals <- log(dataset_after_peak$trait_value)
		x_vals <- 1/(R * (dataset_after_peak$temp + 273.15))
		
		# Take the absolute value.
		DH_H_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
	}
	
	# Set the starting value of DH_L arbitrarily to 70000.
	DH_L_start <- 70000
	
	# T_H50 should be close to the thermal optimum (T_pk), so set the 
	# starting value of T_H50 = T_pk.
	T_H50_start <- T_pk + 273.15
	
	# Set the starting value of T_L50 arbitrarily to the minimum 
	# temperature in the dataset.
	T_L50_start <- min(dataset$temp) + 273.15
	
	# Set the starting value of a arbitrarily to 1.
	a_start <- 1
	
	function_to_be_fitted <- function(B_0, DH_A, DH_L, T_L50, DH_H, T_H50, a, temp)
	{
		temp <- temp + 273.15

		# Set the universal gas constant.
		R <- 1.987
		
		# Set the reference temperature to 0°C.
		T_ref <- 273.15
			
		result <- (B_0 * exp( (DH_A/R) * ((1/T_ref) - (1/temp)) ) ) / ( a + exp( (-DH_L/R) * ((1/T_L50) - (1/temp)) ) + exp( (DH_H/R) * ((1/T_H50) - (1/temp)) ) )
				
		if ( T_L50 >= T_H50 || is.na(result) || any(result <= 0) || DH_A == 0 || a == 0 )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_0, DH_A, DH_L, T_L50, DH_H, T_H50, a, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_0 = 0.5 * B_0_start,		DH_A = 0.5 * DH_A_start,
				DH_L = 0.5 * DH_L_start,	T_L50 = 0.5 * T_L50_start,
				DH_H = 0.5 * DH_H_start,	T_H50 = 0.5 * T_H50_start,
				a = 0.5 * a_start
			),
			start_upper = c(
				B_0 = 1.5 * B_0_start,		DH_A = 1.5 * DH_A_start,
				DH_L = 1.5 * DH_L_start,	T_L50 = 1.5 * T_L50_start,
				DH_H = 1.5 * DH_H_start,	T_H50 = 1.5 * T_H50_start,
				a = 1.5 * a_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0, -20 + 273.15, 0, 0 + 273.15, 0),
			upper = c(Inf, 300000, 300000, 150 + 273.15, 300000, 150 + 273.15, Inf)
		)
	)
	
	return(fit)
}

####################################################
# SIMPLIFIED JOHNSON - LEWIN MODEL (4 parameters)  #
# (Johnson & Lewin, J. Cell. Comp. Physiol., 1946) #
#                                                  #
# Parameters: B_0, E, T_pk, E_D                    #
####################################################
fit_simplified_Johnson_Lewin_4_pars <- function(dataset)
{
	
	# Set the Boltzmann constant.
	k <- 8.617 * 10^-5
	
	# Set the minimum trait measurement at the rise of the TPC as the 
	# starting value for B_0.
	B_0_start <- min(dataset$trait_value[dataset$temp == min(dataset$temp)])

	# Set the starting value of T_pk to the temperature at which the 
	# trait reaches its maximum value.
	T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])

	# If the dataset has measurements before the thermal optimum, 
	# a starting value for E can be set as roughly the slope of the 
	# following regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * 1/(k * (T + 273.15))
	#
	# Otherwise, just set E to 0.6 eV.
	
	dataset_before_peak <- dataset[dataset$temp < T_pk_start,]
	if ( 
		nrow(dataset_before_peak) < 3 || 
		length(unique(dataset_before_peak$temp)) < 3 || 
		length(unique(dataset_before_peak$trait_value)) < 3 
	)
	{
		E_start <- 0.6
	} else
	{
		y_vals <- log(dataset_before_peak$trait_value)
		x_vals <- 1/(k * (dataset_before_peak$temp + 273.15))
		
		# Take the absolute value.
		E_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
		
		if ( E_start >= 10 )
		{
			E_start <- 0.6
		}
	}
	
	# If the dataset has measurements after the thermal optimum, 
	# a starting value for E_D can be set as the slope of the following 
	# regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * 1/(k * (T + 273.15))
	#
	# Otherwise, just set E_D to 3 eV.
	dataset_after_peak <- dataset[dataset$temp > T_pk_start,]
	if ( 
		nrow(dataset_after_peak) < 3 || 
		length(unique(dataset_after_peak$temp)) < 3 || 
		length(unique(dataset_after_peak$trait_value)) < 3 
	)
	{
		E_D_start <- 3
	} else
	{
		y_vals <- log(dataset_after_peak$trait_value)
		x_vals <- 1/(k * (dataset_after_peak$temp + 273.15))
		
		# Take the absolute value.
		E_D_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
		
		if ( E_D_start >= 50 )
		{
			E_D_start <- 3
		}
	}
	
	if ( E_start >= E_D_start )
	{
		E_start <- 0.9 * E_D_start
	}
	
	function_to_be_fitted <- function(B_0, E, T_pk, E_D, temp)
	{
		temp <- temp + 273.15

		# Set the Boltzmann constant.
		k <- 8.617 * 10^-5
				
		if ( E == E_D )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(
				log(
					B_0 * exp((-E/k) * (1/temp)) / ( 1 + ( E/( E_D - E ) ) * exp( (E_D/k) * (1/T_pk - 1/temp) ) )
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_0, E, T_pk, E_D, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_0 = 0.5 * B_0_start,				E = 0.5 * E_start,
				T_pk = 0.5 * T_pk_start + 273.15,	E_D = 0.5 * E_D_start
			),
			start_upper = c(
				B_0 = 1.5 * B_0_start,				E = 1.5 * E_start,
				T_pk = 1.5 * T_pk_start + 273.15,	E_D = 1.5 * E_D_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 273.15, 0),
			upper = c(Inf, 10, 273.15 + 150, 50)
		)
	)
	
	return(fit)
}

########################################################
# SIMPLIFIED SHARPE - SCHOOLFIELD MODEL (6 parameters) #
# (Schoolfield et al., J. Theor. Biol., 1981)          #
#                                                      #
# Parameters: B_0, DH_A, DH_L, T_L50, DH_H, T_H50      #
########################################################
fit_simplified_Sharpe_Schoolfield_6_pars <- function(dataset)
{
	
	# Set the universal gas constant.
	R <- 1.987
	
	# Set the minimum trait measurement at the rise of the TPC as the 
	# starting value for B_0.
	B_0_start <- min(dataset$trait_value[dataset$temp == min(dataset$temp)])
	
	# If the dataset has measurements before the thermal optimum, 
	# a starting value for DH_A can be set as the slope of the following 
	# regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * 1/(R * (T + 273.15))
	#
	# Otherwise, just set DH_A to 15000.
	T_pk <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
	
	dataset_before_peak <- dataset[dataset$temp < T_pk,]
	if ( 
		nrow(dataset_before_peak) < 3 || 
		length(unique(dataset_before_peak$temp)) < 3 || 
		length(unique(dataset_before_peak$trait_value)) < 3 
	)
	{
		DH_A_start <- 15000
	} else
	{
		y_vals <- log(dataset_before_peak$trait_value)
		x_vals <- 1/(R * (dataset_before_peak$temp + 273.15))
		
		# Take the absolute value.
		DH_A_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
	}
	
	# If the dataset has measurements after the thermal optimum, 
	# a starting value for DH_H can be set as roughly the slope of the 
	# following regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * 1/(R * (T + 273.15))
	#
	# Otherwise, just set DH_H to 50000.
	dataset_after_peak <- dataset[dataset$temp > T_pk,]
	if ( 
		nrow(dataset_after_peak) < 3 || 
		length(unique(dataset_after_peak$temp)) < 3 || 
		length(unique(dataset_after_peak$trait_value)) < 3 
	)
	{
		DH_H_start <- 50000
	} else
	{
		y_vals <- log(dataset_after_peak$trait_value)
		x_vals <- 1/(R * (dataset_after_peak$temp + 273.15))
		
		# Take the absolute value.
		DH_H_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
	}
	
	# Set the starting value of DH_L arbitrarily to 70000.
	DH_L_start <- 70000
	
	# T_H50 should be close to the thermal optimum (T_pk), so set the 
	# starting value of T_H50 = T_pk.
	T_H50_start <- T_pk + 273.15
	
	# Set the starting value of T_L50 arbitrarily to the minimum 
	# temperature in the dataset.
	T_L50_start <- min(dataset$temp) + 273.15
	
	function_to_be_fitted <- function(B_0, DH_A, DH_L, T_L50, DH_H, T_H50, temp)
	{
		temp <- temp + 273.15

		# Set the universal gas constant.
		R <- 1.987
		
		# Set the reference temperature to 0°C.
		T_ref <- 273.15
			
		result <- (B_0 * exp( (DH_A/R) * ((1/T_ref) - (1/temp)) ) ) / ( 1 + exp( (-DH_L/R) * ((1/T_L50) - (1/temp)) ) + exp( (DH_H/R) * ((1/T_H50) - (1/temp)) ) )
				
		if ( T_L50 >= T_H50 || is.na(result) || any(result <= 0) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_0, DH_A, DH_L, T_L50, DH_H, T_H50, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_0 = 0.5 * B_0_start,		DH_A = 0.5 * DH_A_start,
				DH_L = 0.5 * DH_L_start,	T_L50 = 0.5 * T_L50_start,
				DH_H = 0.5 * DH_H_start,	T_H50 = 0.5 * T_H50_start
			),
			start_upper = c(
				B_0 = 1.5 * B_0_start,		DH_A = 1.5 * DH_A_start,
				DH_L = 1.5 * DH_L_start,	T_L50 = 1.5 * T_L50_start,
				DH_H = 1.5 * DH_H_start,	T_H50 = 1.5 * T_H50_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0, -20 + 273.15, 0, 0 + 273.15),
			upper = c(Inf, 300000, 300000, 150 + 273.15, 300000, 150 + 273.15)
		)
	)
	
	return(fit)
}

#########################################
# SKEW-NORMAL MODEL (4 parameters)      #
# (Urban et al., Proc. R. Soc. B, 2012) #
#                                       #
# Parameters: a, b, lambda, sigma       #
#########################################
fit_skew_normal_4_pars <- function(dataset)
{
	
	# Set the starting value of a arbitrarily to 1.
	a_start <- 1
	
	# b should be above the temperature at which the trait reaches its 
	# maximum value (T_pk). Thus, set the starting value of b to T_pk.
	b_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
	
	# Set the starting value of lambda arbitrarily to 2.7.
	lambda_start <- 2.7
	
	# Set the starting value of sigma to the difference between the 
	# maximum and minimum temperature in the dataset.
	sigma_start <- max(dataset$temp) - min(dataset$temp)
	
	function_to_be_fitted <- function(a, b, lambda, sigma, temp)
	{
		return(
			log(
				a * ( exp( ( - (temp - b)^2 ) / (sigma^2) ) ) * ( 1 + erf( ( -lambda * (temp - b) ) / sigma ) )
			)
		)
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, b, lambda, sigma, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,				b = 0.5 * b_start,
				lambda = 0.5 * lambda_start,	sigma = 0.5 * sigma_start
			),
			start_upper = c(
				a = 1.5 * a_start,				b = 1.5 * b_start,
				lambda = 1.5 * lambda_start,	sigma = 1.5 * sigma_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, -Inf, 0),
			upper = c(Inf, 150, Inf, 150)
		)
	)
	
	return(fit)
}

############################################
# STEVENSON MODEL (6 parameters)           #
# (Stevenson et al., Physiol. Zool., 1985) #
#                                          #
# Parameters: S, K1, K2, T_min, K3, T_max  #
############################################
fit_Stevenson_6_pars <- function(dataset)
{

	# Set the starting value of S arbitrarily to 1.
	S_start <- 1

	# Set the starting value of T_min 1 degree below the minimum 
	# temperature in the dataset.
	T_min_start <- min(dataset$temp) - 1

	# Set the starting value of T_max 1 degree above the maximum 
	# temperature in the dataset.
	T_max_start <- max(dataset$temp) + 1

	# Set the starting value of K1 arbitrarily to 6.
	K1_start <- 6
	
	# Set the starting value of K2 arbitrarily to 0.4.
	K2_start <- 0.4
	
	# Set the starting value of K3 arbitrarily to 0.7.
	K3_start <- 0.7
			
	function_to_be_fitted <- function(S, K1, K2, T_min, K3, T_max, temp)
	{		
		if ( T_min >= T_max || any(temp < T_min) || any(temp > T_max) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(
				log(
					S * ( 1 / ( 1 + K1 * exp( (-K2) * (temp - T_min) ) ) ) * ( 1 - exp(K3 * (temp - T_max)) )
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				S, K1, K2, T_min, K3, T_max, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				S = 0.5 * S_start,			K1 = 0.5 * K1_start,
				K2 = 0.5 * K2_start,		T_min = 0.5 * T_min_start,	
				K3 = 0.5 * K3_start,		T_max = 0.5 * T_max_start
			),
			start_upper = c(
				S = 1.5 * S_start,			K1 = 1.5 * K1_start,
				K2 = 1.5 * K2_start,		T_min = 1.5 * T_min_start,	
				K3 = 1.5 * K3_start,		T_max = 1.5 * T_max_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0, -20, 0, 0),
			upper = c(Inf, Inf, 1, 150, 1, 150)
		)
	)
	
	return(fit)
}

#####################################
# STINNER MODEL (4 parameters)      #
# (Stinner et al., Can. Ent., 1974) #
#                                   #
# Parameters: a, b, B_pk, T_pk      #
#####################################
fit_Stinner_4_pars <- function(dataset)
{

	# Set the starting value of a arbitrarily to 0.2.
	a_start <- 0.2
	
	# Set the starting value of b arbitrarily to 0.1.
	b_start <- 0.1

	# Set the starting value of T_pk to the temperature at which the 
	# trait reaches its maximum value.
	T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
	
	# Set the starting value of B_pk to the maximum trait value in the 
	# dataset.
	B_pk_start <- max(dataset$trait_value)
	
	function_to_be_fitted <- function(a, b, B_pk, T_pk, temp)
	{
		result <- rep(NA, length(temp))
		
		temps_below_T_pk <- which(temp <= T_pk)
		temps_above_T_pk <- which(temp > T_pk)
		
		result[temps_below_T_pk] <-  B_pk * (1 + exp(a + b * T_pk) ) / ( 1 + exp( a + b * temp[temps_below_T_pk] ) )
		result[temps_above_T_pk] <-  B_pk * (1 + exp(a + b * T_pk) ) / ( 1 + exp( a + b * ( 2 * T_pk - temp[temps_above_T_pk] ) ) )
		
		return(log(result))
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, b, B_pk, T_pk, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,			b = 0.5 * b_start,
				B_pk = 0.5 * B_pk_start, 	T_pk = 0.5 * T_pk_start
			),
			start_upper = c(
				a = 1.5 * a_start,			b = 1.5 * b_start,
				B_pk = 1.5 * B_pk_start,	T_pk = 1.5 * T_pk_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(-Inf, -Inf, 0, 0),
			upper = c(Inf, Inf, Inf, 150)
		)
	)
	
	return(fit)
}

########################################
# TAYLOR - SEXTON MODEL (3 parameters) #
# (Taylor & Sexton, Ecology, 1972)     #
#                                      #
# Parameters: B_pk, T_min, T_pk        #
########################################
fit_Taylor_Sexton_3_pars <- function(dataset)
{
	
	# Set the starting value of B_pk to the maximum trait value in the 
	# dataset.
	B_pk_start <- max(dataset$trait_value)
	
	# Set the starting value of T_min to the minimum temperature in the 
	# dataset.
	T_min_start <- min(dataset$temp)
	
	# Set the starting value of T_pk to the temperature at which the 
	# trait reaches its maximum value.
	T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])

	function_to_be_fitted <- function(B_pk, T_min, T_pk, temp)
	{
		if ( T_pk <= T_min || any(temp < T_min) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(
				log(
					B_pk * ( - (temp - T_min)^4 + 2 * ((temp - T_min)^2) * (T_pk - T_min)^2 ) / (T_pk - T_min)^4
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_pk, T_min, T_pk, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_pk = 0.5 * B_pk_start,	T_min = 0.5 * T_min_start,
				T_pk = 0.5 * T_pk_start
			),
			start_upper = c(
				B_pk = 1.5 * B_pk_start,	T_min = 1.5 * T_min_start,
				T_pk = 1.5 * T_pk_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, -20, 0),
			upper = c(Inf, 150, 150)
		)
	)
	
	return(fit)
}

##################################
# THOMAS I MODEL (4 parameters)  #
# (Thomas et al., Science, 2012) #
#                                #
# Parameters: a, b, w, z         #
##################################
fit_Thomas_I_4_pars <- function(dataset)
{
	
	# Set the starting value of w to the difference between the maximum 
	# and minimum temperature in the dataset.
	w_start <- max(dataset$temp) - min(dataset$temp)
	
	# Set the starting value of a to the maximum trait value in the 
	# dataset.
	a_start <- max(dataset$trait_value)
	
	# Set the starting value of b arbitrarily to 0.3.
	b_start <- 0.3
	
	# z should be close to the thermal optimum. So, set the starting 
	# value of z to the temperature with the maximum trait value.
	z_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])

	function_to_be_fitted <- function(a, b, z, w, temp)
	{
		if ( any( abs( temp - z ) >= ( w / 2 ) ) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(
				log(
					a * exp(b * temp) * (1 - ( ( temp - z ) / (w/2) )^2 ) 
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, b, z, w, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,	b = 0.5 * b_start,
				z = 0.5 * z_start,	w = 0.5 * w_start
			),
			start_upper = c(
				a = 1.5 * a_start,	b = 1.5 * b_start,
				z = 1.5 * z_start,	w = 1.5 * w_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, -Inf, 0, 0),
			upper = c(Inf, Inf, 150, 150)
		)
	)
	
	return(fit)
}

#############################################
# THOMAS II MODEL (5 parameters)            #
# (Thomas et al., Glob. Chang. Biol., 2017) #
#                                           #
# Parameters: a, b, d, e, f                 #
#############################################
fit_Thomas_II_5_pars <- function(dataset)
{
	
	# Set the starting value of a arbitrarily to 1.
	a_start <- 1

	# If the dataset has measurements before the thermal optimum, 
	# a starting value for b can be set as the slope of the following 
	# regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * T
	#
	# Otherwise, just arbitrarily set b to 0.6.
	T_pk <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
	
	dataset_before_peak <- dataset[dataset$temp < T_pk,]
	if ( 
		nrow(dataset_before_peak) < 3 || 
		length(unique(dataset_before_peak$temp)) < 3 || 
		length(unique(dataset_before_peak$trait_value)) < 3 
	)
	{
		b_start <- 0.6
	} else
	{
		y_vals <- log(dataset_before_peak$trait_value)
		x_vals <- dataset_before_peak$temp
		
		# Take the absolute value.
		b_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
	}
	
	# Set the starting value of d arbitrarily to 1.
	d_start <- 1
	
	# Set the starting value of e arbitrarily to 0.4.
	e_start <- 0.2
	
	# If the dataset has measurements after the thermal optimum, 
	# a starting value for f can be set as roughly the slope of 
	# the following regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * T
	#
	# Otherwise, just arbitrarily set f to 0.1.
	dataset_after_peak <- dataset[dataset$temp > T_pk,]
	if ( 
		nrow(dataset_after_peak) < 3 || 
		length(unique(dataset_after_peak$temp)) < 3 || 
		length(unique(dataset_after_peak$trait_value)) < 3 
	)
	{
		f_start <- 0.1
	} else
	{
		y_vals <- log(dataset_after_peak$trait_value)
		x_vals <- dataset_after_peak$temp
		
		# Take the absolute value.
		f_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
	}
	
	function_to_be_fitted <- function(a, b, d, e, f, temp)
	{
		result <- a * exp(b * temp) - (d + e * exp(f * temp) )

		if ( any(is.na(result)) || any(result <= 0) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, b, d, e, f, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,	b = 0.5 * b_start,
				d = 0.5 * d_start,	e = 0.5 * e_start,
				f = 0.5 * f_start
			),
			start_upper = c(
				a = 1.5 * a_start,	b = 1.5 * b_start,
				d = 1.5 * d_start,	e = 1.5 * e_start,
				f = 1.5 * f_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, -Inf, 0, 0),
			upper = c(Inf, Inf, Inf, Inf, Inf)
		)
	)
	
	return(fit)
}

####################################################
# THORNTON - LESSEM MODEL (6 parameters)           #
# (Thornton & Lessem, Trans. Am. Fish. Soc., 1978) #
#                                                  #
# Parameters: a, K1, T_1, T_2, T_3, K4             #
####################################################
fit_Thornton_Lessem_6_pars <- function(dataset)
{
	
	# Set the starting value of a arbitrarily to 1.
	a_start <- 1

	# Set the starting value of T_1 to the minimum temperature in the 
	# dataset.
	T_1_start <- min(dataset$temp)

	# Set the starting value of T_3 to the maximum temperature in the 
	# dataset.
	T_3_start <- max(dataset$temp)

	# Set the starting value of K1 arbitrarily to 0.1.
	K1_start <- 0.1
	
	# Set the starting value of K4 arbitrarily to 0.7.
	K4_start <- 0.7
	
	# Set the starting value of T_2 to the temperature at which the 
	# trait reaches its maximum value.
	T_2_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
		
	function_to_be_fitted <- function(a, K1, T_2, T_1, K4, T_3, temp)
	{
		gamma_1 <- ( 1 / ( T_2 - T_1 ) ) * log( ( 0.98 * (1 - K1) ) / ( K1 * (1 - 0.98) ) )
		gamma_2 <- ( 1 / ( T_3 - T_2 ) ) * log( ( 0.98 * (1 - K4) ) / ( K4 * (1 - 0.98) ) )
		
		result <- a * ( (K1 * exp(gamma_1 * (temp - T_1) ) ) / (1 + K1 * ( exp( gamma_1 * (temp - T_1) ) - 1 ) ) ) * ( ( K4 * exp( gamma_2 * (T_3 - temp) ) ) / ( 1 + K4 * ( exp( gamma_2 * (T_3 - temp) ) - 1 ) ) )

		if (
			T_1 >= T_3 || T_2 <= T_1 || T_2 >= T_3 
		)
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, K1, T_2, T_1, K4, T_3, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,
				K1 = 0.5 * K1_start,		T_2 = 0.5 * T_2_start,
				T_1 = 0.5 * T_1_start,	K4 = 0.5 * K4_start,
				T_3 = 0.5 * T_3_start
			),
			start_upper = c(
				a = 1.5 * a_start,
				K1 = 1.5 * K1_start,		T_2 = 1.5 * T_2_start,
				T_1 = 1.5 * T_1_start,	K4 = 1.5 * K4_start,
				T_3 = 1.5 * T_3_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0,0, 0, -20, 0, 0),
			upper = c(Inf, 1, 150, 150, 1, 150)
		)
	)
	
	return(fit)
}

#############################################################################
# TOMLINSON - MENZ MODEL (4 parameters)                                     #
# (Tomlinson & Menz, Comp. Biochem. Physiol. A Mol. Integr. Physiol., 2015) #
#                                                                           #
# Parameters: B_0, a, b, d                                                  #
#############################################################################
fit_Tomlinson_Menz_4_pars <- function(dataset)
{

	# Set the minimum trait measurement at the rise of the TPC as the 
	# starting value for B_0.
	B_0_start <- min(dataset$trait_value[dataset$temp == min(dataset$temp)])

	T_pk <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])

	# If the dataset has measurements before the thermal optimum, 
	# a starting value for a can be set as roughly the slope of the 
	# following regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * T
	#
	# Otherwise, just set a to 0.6.
	
	dataset_before_peak <- dataset[dataset$temp < T_pk,]
	if ( 
		nrow(dataset_before_peak) < 3 || 
		length(unique(dataset_before_peak$temp)) < 3 || 
		length(unique(dataset_before_peak$trait_value)) < 3 
	)
	{
		a_start <- 0.6
	} else
	{
		y_vals <- log(dataset_before_peak$trait_value)
		x_vals <- dataset_before_peak$temp
		
		# Take the absolute value.
		a_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
	}

	# Set the starting value of b arbitrarily to the lowest temperature 
	# in the dataset.
	b_start <- min(dataset$temp)

	# d should be close to the thermal optimum. Thus, set the starting 
	# value of d to T_pk.
	d_start <- T_pk
		
	function_to_be_fitted <- function(B_0, a, b, d, temp)
	{
		result <- B_0 * ( exp(a * temp) - exp(b - temp) - exp(temp - d) )
		
		if ( any(result <= 0) || b >= d )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_0, a, b, d, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_0 = 0.5 * B_0_start,	a = 0.5 * a_start,
				b = 0.5 * b_start,		d = 0.5 * d_start
			),
			start_upper = c(
				B_0 = 1.5 * B_0_start,	a = 1.5 * a_start,
				b = 1.5 * b_start,		d = 1.5 * d_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, -20, 0),
			upper = c(Inf, Inf, 150, 150)
		)
	)
	
	return(fit)
}

####################################################
# TOMLINSON - PHILLIPS MODEL (3 parameters)        #
# (Tomlinson & Phillips, J. Insect Physiol., 2015) #
#                                                  #
# Parameters: B_0, a, b                            #
####################################################
fit_Tomlinson_Phillips_3_pars <- function(dataset)
{

	# Set the minimum trait measurement at the rise of the TPC as the 
	# starting value for B_0.
	B_0_start <- min(dataset$trait_value[dataset$temp == min(dataset$temp)])

	T_pk <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])

	# If the dataset has measurements before the thermal optimum, 
	# a starting value for a can be set as roughly the slope of the 
	# following regression for that subset of the data:
	#
	# ln(trait_value) ~ intercept + slope * T
	#
	# Otherwise, just set a to 0.6.
	
	dataset_before_peak <- dataset[dataset$temp < T_pk,]
	if ( 
		nrow(dataset_before_peak) < 3 || 
		length(unique(dataset_before_peak$temp)) < 3 || 
		length(unique(dataset_before_peak$trait_value)) < 3 
	)
	{
		a_start <- 0.6
	} else
	{
		y_vals <- log(dataset_before_peak$trait_value)
		x_vals <- dataset_before_peak$temp
		
		# Take the absolute value.
		a_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
	}

	# b should be close to the thermal optimum. Thus, set the starting 
	# value of b to T_pk.
	b_start <- T_pk
	
	function_to_be_fitted <- function(B_0, a, b, temp)
	{
		result <- B_0 * ( exp(a * temp) - exp(temp - b) )
		
		if ( any(result <= 0) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_0, a, b, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_0 = 0.5 * B_0_start,	a = 0.5 * a_start,
				b = 0.5 * b_start
			),
			start_upper = c(
				B_0 = 1.5 * B_0_start,	a = 1.5 * a_start,
				b = 1.5 * b_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0),
			upper = c(Inf, Inf, 150)
		)
	)
	
	return(fit)
}

##########################################
# VAN'T HOFF MODEL (4 parameters)        #
# (Portner et al., Biogeosciences, 2010) #
#                                        #
# Parameters: a, b, d, e                 #
##########################################
fit_Vant_Hoff_4_pars <- function(dataset)
{

	# Set the starting value of a arbitrarily to 1.
	a_start <- 1
	
	# Set the starting value of b arbitrarily to 0.6.
	b_start <- 0.6
	
	# Set the starting value of d arbitrarily to 1.
	d_start <- 1
	
	# Set the starting value of e arbitrarily to 2.
	e_start <- 2
		
	function_to_be_fitted <- function(a, b, d, e, temp)
	{
		temp <- temp + 273.15
		
		if ( b == 0 || d == 0 || e == 0 )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(
				log(
					a * exp(-b/temp) * (temp^d) * exp(e*temp) 
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, b, d, e, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,	b = 0.5 * b_start,
				d = 0.5 * d_start,	e = 0.5 * e_start
			),
			start_upper = c(
				a = 1.5 * a_start,	b = 1.5 * b_start,
				d = 1.5 * d_start,	e = 1.5 * e_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, -Inf, -Inf, -Inf),
			upper = c(Inf, Inf, Inf, Inf)
		)
	)
	
	return(fit)
}

######################################
# WANG - ENGEL MODEL (4 parameters)  #
# (Wang & Engel, Agric. Syst., 1998) #
#                                    #
# Parameters: B_pk, T_min, T_pk, a   #
######################################
fit_Wang_Engel_4_pars <- function(dataset)
{

	# Set the starting value of B_pk to the maximum trait value in the 
	# dataset.
	B_pk_start <- max(dataset$trait_value)

	# Set the starting value of T_min to the minimum temperature in the 
	# dataset.
	T_min_start <- min(dataset$temp)

	# Set the starting value of T_pk to the temperature at which the 
	# trait reaches its maximum value.
	T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])

	# Set the starting value of a arbitrarily to 2.
	a_start <- 2

	function_to_be_fitted <- function(B_pk, T_min, T_pk, a, temp)
	{
		T_max <- (2^(1/a)) * (T_pk - T_min) + T_min	
		
		if ( T_min >= T_max || T_pk <= T_min || T_pk >= T_max || any(temp < T_min) || any(temp > T_max) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(
				log(
					B_pk * ( ( 2 * ((temp - T_min)^a) * ((T_pk - T_min)^a) - ((temp - T_min)^(2*a)) ) / ( (T_pk - T_min)^(2*a) ) )
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_pk, T_min, T_pk, a, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_pk = 0.5 * B_pk_start,	T_min = 0.5 * T_min_start,
				T_pk = 0.5 * T_pk_start,	a = 0.5 * a_start
			),
			start_upper = c(
				B_pk = 1.5 * B_pk_start,	T_min = 1.5 * T_min_start,
				T_pk = 1.5 * T_pk_start,	a = 1.5 * a_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, -20, 0, 0),
			upper = c(Inf, 150, 150, Inf)
		)
	)
	
	return(fit)
}

#################################################
# WANG - LAN - DING MODEL (7 parameters)        #
# (Wang et al., Acta Ecol. Sin., 1982)          #
#                                               #
# Parameters: a, K_l, T_min, K_u, T_max, b, T_0 #
#################################################
fit_Wang_Lan_Ding_7_pars <- function(dataset)
{

	# Set the starting value of T_min to the minimum temperature in the 
	# dataset.
	T_min_start <- min(dataset$temp)

	# Set the starting value of T_max to the maximum temperature in the 
	# dataset.
	T_max_start <- max(dataset$temp)

	# Set the starting value of K_l arbitrarily to 0.1.
	K_l_start <- 0.1
	
	# Set the starting value of K_u arbitrarily to 0.7.
	K_u_start <- 0.7
	
	# Set the starting value of a arbitrarily to 1.
	a_start <- 1
	
	# Set the starting value of b arbitrarily to 0.5.
	b_start <- 0.5
	
	# T_0 should be close to the thermal optimum. Thus, set the starting 
	# value of T_0 to the temperature at which the trait reaches its 
	# maximum value.
	T_0_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])

	function_to_be_fitted <- function(a, K_l, T_min, K_u, T_max, b, T_0, temp)
	{
		
		if ( T_min >= T_max || T_0 <= T_min || T_0 >= T_max || any(temp < T_min) || any(temp > T_max) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(
				log(
					a * ( 1 - exp(-K_l * (temp - T_min) ) ) * ( 1 - exp( K_u * ( temp - T_max ) )  ) / ( 1 + exp( - b * (temp - T_0) ) )
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				a, K_l, T_min, K_u, T_max, b, T_0, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				a = 0.5 * a_start,			K_l = 0.5 * K_l_start,
				T_min = 0.5 * T_min_start,	K_u = 0.5 * K_u_start,
				T_max = 0.5 * T_max_start,	b = 0.5 * b_start,
				T_0 = 0.5 * T_0_start
			),
			start_upper = c(
				a = 1.5 * a_start,			K_l = 1.5 * K_l_start,
				T_min = 1.5 * T_min_start,	K_u = 1.5 * K_u_start,
				T_max = 1.5 * T_max_start,	b = 1.5 * b_start,
				T_0 = 1.5 * T_0_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, -20, 0, 0, 0, 0),
			upper = c(Inf, 1, 150, 1, 150, 1, 150)
		)
	)
	
	return(fit)
}

#########################################
# WARREN - DREYER MODEL (3 parameters)  #
# (Warren & Dreyer, J. Exp. Bot., 2006) #
#                                       #
# Parameters: B_pk, T_pk, a             #
#########################################
fit_Warren_Dreyer_3_pars <- function(dataset)
{
	
	# Set the starting value of B_pk to the maximum trait value in the 
	# dataset.
	B_pk_start <- max(dataset$trait_value)
		
	# Set the starting value of T_pk to the temperature at which the 
	# trait reaches its maximum value.
	T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
	
	# Set the starting value of a arbitrarily to 1.
	a_start <- 1
	
	function_to_be_fitted <- function(B_pk, T_pk, a, temp)
	{
		return(
			log(
				B_pk * exp( - 0.5 * ( log(temp/T_pk) / a )^2 )
			)
		)
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_pk, T_pk, a, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_pk = 0.5 * B_pk_start,	T_pk = 0.5 * T_pk_start,
				a = 0.5 * a_start
			),
			start_upper = c(
				B_pk = 1.5 * B_pk_start,	T_pk = 1.5 * T_pk_start,
				a = 1.5 * a_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, 0),
			upper = c(Inf, 150, Inf)
		)
	)
	
	return(fit)
}

#####################################
# WEIBULL MODEL (4 parameters)      #
# (Shi & Ge, J. Therm. Biol., 2010) #
#                                   #
# Parameters: B_pk, T_pk, d, e      #
#####################################
fit_Weibull_4_pars <- function(dataset)
{
	
	# Set the starting value of B_pk to the maximum trait value in the 
	# dataset.
	B_pk_start <- max(dataset$trait_value)
	
	# Set the starting value of T_pk to the temperature at which the 
	# trait reaches its maximum value.
	T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
	
	# Set the starting value of d arbitrarily to 90.
	d_start <- 90
	
	# Set the starting value of e arbitrarily to 2.
	e_start <- 2

	function_to_be_fitted <- function(B_pk, T_pk, d, e, temp)
	{
		result <- B_pk * (( (e - 1)/e )^( (1 - e)/e )) * ( ( ( (temp - T_pk)/d ) + ( (e - 1)/e )^(1/e) ) ^ (e - 1) ) * exp( - ( ( ( temp -T_pk ) / d ) + ( ( e - 1 )/e )^(1/e) )^e  + ( (e - 1)/e ) )  
		
		if ( any(is.na(result)) || any(result <= 0 ) )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(log(result))
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_pk, T_pk, d, e, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_pk = 0.5 * B_pk_start,	T_pk = 0.5 * T_pk_start,
				d = 0.5 * d_start,			e = 0.5 * e_start
			),
			start_upper = c(
				B_pk = 1.5 * B_pk_start,	T_pk = 1.5 * T_pk_start,
				d = 1.5 * d_start,			e = 1.5 * e_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, 0, -Inf, -Inf),
			upper = c(Inf, 150, Inf, Inf)
		)
	)
	
	return(fit)
}

########################################
# YAN - HUNT MODEL (4 parameters)      #
# (Yan & Hunt, Ann. Bot., 1999)        #
#                                      #
# Parameters: B_pk, T_min, T_max, T_pk #
########################################
fit_Yan_Hunt_4_pars <- function(dataset)
{

	# Set the starting value of B_pk to the maximum trait value in the 
	# dataset.
	B_pk_start <- max(dataset$trait_value)
	
	# Set the starting value of T_min 1 degree below the minimum 
	# temperature in the dataset.
	T_min_start <- min(dataset$temp) - 1

	# Set the starting value of T_max 1 degree above the maximum 
	# temperature in the dataset.
	T_max_start <- max(dataset$temp) + 1
		
	# Set the starting value of T_pk to the temperature at which the 
	# trait reaches its maximum value.
	T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
		
	function_to_be_fitted <- function(B_pk, T_min, T_max, T_pk, temp)
	{
		if ( T_min >= T_max || any(temp <= T_min) || any(temp >= T_max) || T_min >= T_pk || T_max <= T_pk )
		{
			return(rep(1e10, length(temp)))
		} else
		{
			return(
				log(
					B_pk * ( (T_max - temp) / (T_max - T_pk) ) * ( (temp - T_min) / (T_pk - T_min) )^( (T_pk - T_min) / (T_max - T_pk) )
				)
			)
		}
	}
	
	fit <- NULL
	
	try(
		fit <- nls_multstart(
			log(trait_value) ~ function_to_be_fitted(
				B_pk, T_min, T_max, T_pk, temp = temp
			),
			data = dataset,
			iter = 1000,
			start_lower = c(
				B_pk = 0.5 * B_pk_start,	T_min = 0.5 * T_min_start,
				T_max = 0.5 * T_max_start,	T_pk = 0.5 * T_pk_start
			),
			start_upper = c(
				B_pk = 1.5 * B_pk_start,	T_min = 1.5 * T_min_start,
				T_max = 1.5 * T_max_start,	T_pk = 1.5 * T_pk_start
			),
			supp_errors = 'Y',
			convergence_count = FALSE,
			control = nls.lm.control(
				ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
				maxfev = 100000
			),
			lower = c(0, -20, 0, 0),
			upper = c(Inf, 150, 150, 150)
		)
	)
	
	return(fit)
}
