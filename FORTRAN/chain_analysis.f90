module chain_analysis_mod
use types_mod
use maths_mod
use global_vars_mod, only : one_sigma, two_sigma, neural, PI, what_to_do, database
use mcmc_mod, only : read_filter
use ann_mod, only : neural_eval, lininterpol_eval
implicit none

contains

!------------------------------------------------------------
! Read the chains
!------------------------------------------------------------
	subroutine read_chains_multinest(fich_chain, chain_analysis)
	type(markov_chain_type) :: chain_analysis
	character(len=200) :: fich_chain
	integer :: nparams, chain_length, i, j
	logical :: exists

		i = 0
		
		inquire(file=trim(adjustl(fich_chain))//'post_equal_weights.dat', exist=exists)
		if (.not.exists) then
			write(*,*) 'File not found : ', trim(adjustl(fich_chain))//'post_equal_weights.dat'
			stop
		endif
		
		open(unit=12,file=trim(adjustl(fich_chain))//'post_equal_weights.dat',action='read',status='old')

! This works only in Intel Fortran
! 		do while (.not.eof(12))
! 			read(12,*)
! 			i = i + 1
! 		enddo
! 		close(12)

		do while (.true.)
   		read (12, *, end=999)
   		i = i + 1
		enddo
999   continue
		close(12)
		
		chain_analysis%nlength = i
		print *, 'Length of posterior samples : ', chain_analysis%nlength
		print *, 'Number of parameters : ', chain_analysis%nparam
		
		allocate(chain_analysis%chain(chain_analysis%nparam+1,chain_analysis%nlength))
		open(unit=12,file=trim(adjustl(fich_chain))//'post_equal_weights.dat',action='read',&
			status='old')
		do i = 1, chain_analysis%nlength
			read(12,*) (chain_analysis%chain(j,i),j=1,chain_analysis%nparam+1)
		enddo
		close(12)

! Read evidence
		inquire(file=trim(adjustl(fich_chain))//'stats.dat', exist=exists)
		if (.not.exists) then
			write(*,*) 'File not found : ', trim(adjustl(fich_chain))//'stats.dat'
			stop
		endif
		
		open(unit=12,file=trim(adjustl(fich_chain))//'stats.dat',action='read',status='old')
		read(12,FMT='(16X,E20.12)') chain_analysis%evidence
		close(12)
		
! Read the MAP parameters
		open(unit=12,file=trim(adjustl(fich_chain))//'.MAP',action='read',status='old')
		read(12,*) chain_analysis%best_parameters
		close(12)
		
	end subroutine read_chains_multinest

!------------------------------------------------------------
! Remove burn-in
!------------------------------------------------------------
	subroutine remove_burnin(chain_analysis)
	type(markov_chain_type) :: chain_analysis
	real(kind=KIND_FLOAT), pointer :: chain_new(:,:), temp(:)
	real(kind=KIND_FLOAT) :: percent
	integer :: burn_in
				
		write(*,*) 'Burn-in of ', chain_analysis%burnin, '%'
		
		burn_in = chain_analysis%niter_max * chain_analysis%burnin / 100.d0
		if (burn_in > chain_analysis%niter_max) burn_in = chain_analysis%niter_max
						
		write(*,*) 'Leaving a chain of length : ', chain_analysis%niter_max-burn_in

! Change length of chain without burn-in
		chain_analysis%niter_max = chain_analysis%niter_max-burn_in

! Re-allocate new memory
		allocate(chain_new(chain_analysis%nparam,chain_analysis%niter_max))
		chain_new = chain_analysis%chain(:,burn_in:)
		
		deallocate(chain_analysis%chain)
		allocate(chain_analysis%chain(chain_analysis%nparam,chain_analysis%niter_max))
		chain_analysis%chain = chain_new
				
		deallocate(chain_new)

! Do the same with the acceptance rate and the log-likelihood
		allocate(temp(chain_analysis%niter_max))
		temp = chain_analysis%posterior(burn_in:)		
		deallocate(chain_analysis%posterior)
		allocate(chain_analysis%posterior(chain_analysis%niter_max))
		chain_analysis%posterior = temp

		temp = chain_analysis%acceptance_rate(burn_in:)
		deallocate(chain_analysis%acceptance_rate)
		allocate(chain_analysis%acceptance_rate(chain_analysis%niter_max))
		chain_analysis%acceptance_rate = temp
				
		deallocate(temp)
			
	end subroutine remove_burnin

!------------------------------------------------------------
! Get the optimal binning for building a histogram of a variable
! Ref: Freedman, D, & Diaconis, P. (1981). "On the histogram as a density
!      estimator: L2 theory". Zeitschrift für Wahrscheinlichkeitstheorie und
!      verwandte Gebiete 57 (4): 453–476
!  BIN=2*IQR(var)/n^(1/3)   with n the number of elements of the array 'var'
!------------------------------------------------------------
	function optbin(var)
	real(kind=KIND_FLOAT) :: optbin, var(:)
	integer :: n, quart1, quart3
	integer, allocatable :: ind(:)
	
		n = size(var)
		allocate(ind(n))
		
		call qsortd(var,ind,n)
		
		quart1 = n*0.25
		quart3 = n*0.75
		
		optbin = 2.0*(var(ind(quart3)) - var(ind(quart1))) / (n*1.0)**(1.0/3.0)
		
	end function optbin

!------------------------------------------------------------
! 1D histograms
!------------------------------------------------------------
	subroutine oned_histogram(chain, x, yGauss, yStep)
	real(kind=KIND_FLOAT) :: chain(:)
	real(kind=KIND_FLOAT) :: step, xmin, xmax
	real(kind=KIND_FLOAT), pointer :: x(:), yGauss(:), yStep(:)
	real(kind=KIND_FLOAT), allocatable :: xx(:)
	integer :: i, n, nn, ncut
	real(kind=KIND_FLOAT) :: wei, norm, error_norm, sig
	
! Bins
		step = optbin(chain)

		xmin = minval(chain)
		xmax = maxval(chain)
		n = (xmax-xmin) / step
		step = (xmax-xmin) / float(n)
		
! Variables
		allocate(x(n))
		allocate(yGauss(n))
		allocate(yStep(n))
		
		do i = 1, n
			x(i) = xmin + step * (i-1.0)
		enddo
		
! Variables for the normalization
		if (n <= 50) then
			nn = 50
			allocate(xx(nn))
			do i = 1, nn
				xx(i) = xmin + (xmax-xmin)/nn * i
			enddo
		else
			allocate(xx(n))
			xx = x
		endif
				
! Doing the histogram
		sig = 1.2d0
		do i = 1, n
			ncut = count(chain >= x(i)-step/2.d0 .and. chain < x(i)+step/2.0)
			
! Smoothing kernel
			wei = sum(exp(-0.50 * (chain-x(i))**2 / (sig*step)**2.0))
			norm = int_tabulated(xx, exp(-0.50 * (xx-x(i))**2 / (sig*step)**2.0))
			yGauss(i) = wei / norm
			yStep(i) = ncut
		enddo
		
		deallocate(xx)
		
		yGauss = yGauss / maxval(yGauss)
		yStep = yStep / maxval(yStep)
		
	end subroutine oned_histogram

!------------------------------------------------------------
! Confidence intervals
! INPUT
!   xdata : x position of the histogram
!   ydata : y position of the histogram
! OUTPUT
!   est : estimated value of the parameter
! errup : value of x at the upper confidence interval
! errdown : value of x at the lower confidence interval
! conf : desired confidence level
!------------------------------------------------------------
	subroutine conf_limits_ccf(xdata, ydata, conf, est, errup, errdown)
	real(kind=KIND_FLOAT) :: xdata(:), ydata(:)
	real(kind=KIND_FLOAT) :: conf, est, errup, errdown
	real(kind=KIND_FLOAT) :: xpdf_low, xpdf_up, xpdf_mid, xmin, xmax, lower, upper, norm
	real(kind=KIND_FLOAT), allocatable :: xx(:), yy(:), x(:), y(:), xF(:), F(:)
	integer, allocatable :: ind(:)
	integer :: n, npto, loc(1), i
	
		n = size(xdata)
		
		xpdf_low = 0.50 - conf/200.
		xpdf_mid = 0.5
		xpdf_up  = 0.5 + conf/200.
		
		if (xpdf_low < 0.0 .or. xpdf_up > 1) then 
			print *,'wrong value for CONF'
			stop
		endif
				
		allocate(xx(n))
		allocate(yy(n))
		allocate(ind(n))
		
		xx = xdata
		yy = ydata / maxval(ydata)
						
! Sorting
		call qsortd(xx,ind,n)
		xx = xx(ind)
		yy = yy(ind)
		deallocate(ind)
						
! Interpolation by splines
		npto = 300
		xmin = xx(1)
		xmax = xx(n)
		
		if (npto > n) then
			allocate(x(npto))
			allocate(y(npto))
			
			do i = 1, npto
				x(i) = xmin + (xmax-xmin)/npto * i
			enddo
			call spline(xx,yy,x,y)
			n = npto
		else
			allocate(x(n))
			allocate(y(n))
			x = xx
			y = yy
		endif
		
! Peak normalization
		y = y / maxval(y)
		
! Computing the cumulative distribution function
		norm = int_tabulated(x,y)
		
		allocate(xF(0:n-1))
		allocate(F(0:n-1))
		xF = x(1:n)
				
		do i = 5, n
			F(i-1) = int_tabulated(x(1:i-1),y(1:i-1)) / norm
		enddo
				
		loc = minloc(abs(F - xpdf_low))		
		lower = xF(loc(1))
		
		loc = minloc(abs(F - xpdf_mid))		
		est = xF(loc(1))
		
		loc = minloc(abs(F - xpdf_up))
		upper = xF(loc(1))
		
		errup = upper - est
		errdown = est - lower

		deallocate(x,y)
		deallocate(xF,F)
		deallocate(xx,yy)
	
	end subroutine conf_limits_ccf

!------------------------------------------------------------
! Analyze the chains
!------------------------------------------------------------
	subroutine analyze_chains(chain_analysis, prior)
	type(markov_chain_type) :: chain_analysis	
	type(prior_type) :: prior
	character(len=120) :: str_parameter
	real(kind=KIND_FLOAT), pointer :: x(:), yGauss(:), yStep(:)
	integer :: i, n
	real(kind=KIND_FLOAT) :: est, errup1s, errdown1s, errup2s, errdown2s

! Remove burn-in
		call remove_burnin(chain_analysis)

! Save chain to a file
		open(unit=12,file=trim(adjustl(chain_analysis%filename))//'.chain',&
			action='write',status='replace',	form='unformatted')
		write(12) chain_analysis%nparam, chain_analysis%niter_max
		write(12) chain_analysis%chain
		write(12) chain_analysis%posterior
		write(12) chain_analysis%acceptance_rate
		close(12)

! Save also marginalized posterior distributions
		open(unit=12,file=trim(adjustl(chain_analysis%filename))//'.hist1D',&
			action='write',status='replace',&
			form='unformatted')
		write(12) chain_analysis%nparam

! And confidence intervals
		open(unit=13,file=trim(adjustl(chain_analysis%filename))//'.confidence',&
			action='write',status='replace')
					
		do i = 1, chain_analysis%nparam				
			if (prior%typ(i) /= 'D') then
			
				call oned_histogram(chain_analysis%chain(i,:), x, yGauss, yStep)

! 1sigma confidence intervals
				call conf_limits_ccf(x, yStep, one_sigma, est, errup1s, errdown1s)
				call conf_limits_ccf(x, yStep, two_sigma, est, errup2s, errdown2s)
			else
				allocate(x(10))
				x = prior%mu(i)
				allocate(yGauss(10))
				yGauss = 1.0
				allocate(yStep(10))
				yStep = 1.0
				est = prior%mu(i)
				errdown1s = 0.d0
				errup1s = 0.d0
				errdown2s = 0.d0
				errup2s = 0.d0
			endif
			
			chain_analysis%most_probable(i) = est
			write(13,*) est, errdown1s, errup1s, errdown2s, errup2s
			str_parameter = chain_analysis%param_names(i)
			write(*,FMT='(A,A,F9.4)') trim(adjustl(str_parameter)), ' : E(x)=', est
			write(*,FMT='(A,F9.4,A,F9.4)') '1-sigma : -', errdown1s, '  +', errup1s
			write(*,FMT='(A,F9.4,A,F9.4)') '2-sigma : -', errdown2s, '  +', errup2s
			n = size(x)
			write(12) i, n
			write(12) x, yGauss, yStep
			deallocate(x)
			deallocate(yGauss)
			deallocate(yStep)
		enddo

		write(13,*) chain_analysis%best_parameters
		
		close(12)
		close(13)
		
	end subroutine analyze_chains

!------------------------------------------------------------
! Bayesian information gain
!------------------------------------------------------------
	subroutine information_gain(chain_analysis, prior)
	type(markov_chain_type) :: chain_analysis
	type(prior_type) :: prior
	integer :: Ntot_filters, i, j, k, loop, n_remaining_filters, M, N, ind
	logical :: use_filter
	real(kind=KIND_FLOAT) :: temp, noise, pars(9), posterior, minobs, maxobs, lambda0, dj
	type(filter_type), pointer :: filter_new(:)
	character(len=15), allocatable :: filter_names(:)
	character(len=15) :: str
	real(kind=KIND_FLOAT), allocatable :: inf_gain(:), djs(:), logdjs(:), predictive(:), model_eval(:), differential_entropy(:)
	integer :: which, spectroscopy_from, spectroscopy_to
	
! Read list of filters
		open(unit=12, file='FILTERS/filters_eu.dat',action='read',status='old')
		read(12,*) Ntot_filters
		allocate(filter_names(Ntot_filters))
		
		n_remaining_filters = Ntot_filters - chain_analysis%observation%npoints
		allocate(filter_new(chain_analysis%observation%npoints+1))

! Read the name of all filters
		do i = 1, Ntot_filters
			read(12,*) str, temp
			str = trim(adjustl(str))			
			filter_names(i) = str
		enddo
		close(12)

! Fill the new filter structure with the filters in the SED
		do i = 1, chain_analysis%observation%npoints
			temp = read_filter(chain_analysis%observation%filter(i)%name, filter_new(i))
		enddo

		allocate(inf_gain(n_remaining_filters))

		M = 100
		N = 200

		allocate(djs(M))
		allocate(logdjs(M))
		allocate(predictive(M))
		allocate(differential_entropy(M))
		allocate(model_eval(chain_analysis%observation%npoints+1))
		ind = 1

		open(unit=12,file=trim(adjustl(chain_analysis%filename))//'_predictive.dat',action='write',status='replace')
		open(unit=13,file=trim(adjustl(chain_analysis%filename))//'_information_gain.dat',action='write',status='replace')
		write(12,*) M

		minobs = 1.d10
		maxobs = -1.d10
		
! Compute range of variation of the SEDs for the later integration
! These will be used to compute relatively good estimations of the lower and
! upper limits for the integral over observations
		do i = 1, chain_analysis%nlength
			pars = chain_analysis%chain(1:9,i)

! Points obtained with new filter
 			if (what_to_do > 0) then
				call neural_eval(pars, neural%SED_highres, prior%include_agn, prior%reddening_law)
				neural%SED_highres = neural%SED_highres / 1.e10 * 10.d0**pars(7)
				if (maxval(neural%SED_highres) > maxobs) maxobs = maxval(neural%SED_highres)
				if (minval(neural%SED_highres) < minobs) minobs = minval(neural%SED_highres)
			else
				call lininterpol_eval(pars, database%SED_highres, prior%include_agn, prior%reddening_law)
				database%SED_highres = database%SED_highres / 1.e10 * 10.d0**pars(7)
				if (maxval(database%SED_highres) > maxobs) maxobs = maxval(database%SED_highres)
				if (minval(database%SED_highres) < minobs) minobs = minval(database%SED_highres)
			endif

		enddo

 		minobs = log10(minobs)-0.5
 		maxobs = log10(maxobs)+0.5

! Read the rest of filters
		do loop = 1, Ntot_filters

! Verify that this is not already included in the SED so that we only use new filters
			use_filter = .TRUE.
			do j = 1, chain_analysis%observation%npoints
				if (index(chain_analysis%observation%filter(j)%name, filter_names(loop)) /= 0) then
					use_filter = .FALSE.
				endif
			enddo

! If we have to use the filter
			if (use_filter) then

! Read it
				temp = read_filter(filter_names(loop), filter_new(chain_analysis%observation%npoints+1))

				predictive = 0.d0

! Compute the values of the observation for the quadrature int(o,p(o|D,If)*log(p(o|D,If)))
				do j = 1, M
				
					logdjs(j) = minobs + (j-1.d0)/(M-1.d0) * (maxobs-minobs)
					djs(j) = 10.d0**(logdjs(j))

! Noise. We use the following typical values:
! 10% for lambda <= 5 micron
! 20% for 5<lambda<=25
! 30% for 25<lambda<=200
! 50% for lambda>200
					lambda0 = filter_new(chain_analysis%observation%npoints+1)%central
					if (lambda0 <= 5.d0) then
						noise = 0.1d0 * djs(j)
					endif
					if (lambda0 > 5.d0 .and. lambda0 <= 25.d0) then
						noise = 0.2d0 * djs(j)
					endif
					if (lambda0 > 25.d0 .and. lambda0 <= 200.d0) then
						noise = 0.3d0 * djs(j)
					endif
					if (lambda0 > 200.d0) then
						noise = 0.3d0 * djs(j)
					endif

					do i = 1, N
						
						which = int((chain_analysis%nlength-1) * randomu()) + 1

! Synthesize SED in new filter set
						pars = chain_analysis%chain(1:9,which)

! Points obtained with new filter
 						k = chain_analysis%observation%npoints+1
 						
 						if (what_to_do > 0) then
							call neural_eval(pars, neural%SED_highres, prior%include_agn, prior%reddening_law)
							call lin_interpol(neural%lambda*(one+pars(9)), neural%SED_highres, &
								filter_new(k)%lambda, filter_new(k)%SED_interpol)
						else
							call lininterpol_eval(pars, database%SED_highres, prior%include_agn, prior%reddening_law)
							call lin_interpol(database%lambda*(one+pars(9)), database%SED_highres, &
								filter_new(k)%lambda, filter_new(k)%SED_interpol)
						endif
 
  						model_eval(k) = &
  							int_tabulated(filter_new(k)%lambda, filter_new(k)%SED_interpol*filter_new(k)%transmission) /&
  							int_tabulated(filter_new(k)%lambda, filter_new(k)%transmission)

						model_eval = model_eval / 1.e10 * 10.d0**pars(7)

! Compute predictive distribution as a Monte Carlo estimation. We sample
! from the posterior and compute the likelihood of the potential value of the new observation given
! by djs(j)
						posterior = (djs(j) - model_eval(k))**2 / noise**2						
						posterior = 1.d0 / (sqrt(2.d0*PI)*noise) * exp(-0.5d0 * posterior)						
						predictive(j) = predictive(j) + posterior
						
					enddo

					predictive(j) = predictive(j) / (1.d0*N)
					
				enddo

! We have to compute the entropy of the predictive distribution
				differential_entropy = predictive * log(predictive)


				predictive = predictive / (log(10.d0) * int_tabulated(logdjs, djs*predictive))

				do j = 1, M
					if (predictive(j) <= 0.d0) differential_entropy(j) = 0.d0
					write(12,*) djs(j), predictive(j)
! 					print *, djs(j), predictive(j)
				enddo
				write(12,*)

! Compute the integral. We use the integral in log-scale to be more precise				
				inf_gain(ind) = -log(10.d0) * int_tabulated(logdjs, djs*differential_entropy) - 0.5d0*log(2.d0*PI*exp(1.d0)*noise**2)
				
				write(*,FMT='(I2,A,I2,3X,A,A,3X,A,F9.5,3X,A,F10.5)') ind, '/', n_remaining_filters, 'Filter= ', filter_names(loop),&
					'Wavelength=', filter_new(chain_analysis%observation%npoints+1)%central, 'EU=', inf_gain(ind)
				write(13,*) filter_names(loop), filter_new(chain_analysis%observation%npoints+1)%central, inf_gain(ind)
				ind = ind + 1
			endif

		enddo


!--------------------------------------------------------------
! Now the spectroscopy between 8 and 13 micron
! We do this using a Montecarlo approach, since the direct integral has to be
! carried in too many dimensions
!--------------------------------------------------------------


! ! First locate the index for 8 and 13 micron
! 100		if (what_to_do > 0) then
! 			spectroscopy_from = close_to(neural%lambda, 8.0)
! 			spectroscopy_to = close_to(neural%lambda, 13.0)
! 		else
! 			spectroscopy_from = close_to(database%lambda, 8.0)
! 			spectroscopy_to = close_to(database%lambda, 13.0)
! 		endif
! 
! 		print *, spectroscopy_from, spectroscopy_to
! 
! 		M = 1000
! 		N = 250
! 
! 		deallocate(djs)
! 		deallocate(predictive)
! 
! 		allocate(djs(M))
! 		allocate(predictive(M))
! 		
! 		do j = 1, M
! 
! 			do i = 1, N
! 
! 				which = int((chain_analysis%nlength-1) * randomu()) + 1
! 
! ! Synthesize SED in new filter set
! 				pars = chain_analysis%chain(1:9,which)
! 
! ! Points obtained with new filter
! 				k = chain_analysis%observation%npoints+1
! 
! 				if (what_to_do > 0) then
! 					call neural_eval(pars, neural%SED_highres, prior%include_agn, prior%reddening_law)
! 					neural%SED_highres = neural%SED_highres / 1.e10 * 10.d0**pars(7)
! 					posterior = 0.d0
! 					do k = spectroscopy_from, spectroscopy_to
! 						djs(k) = noise * randomn() + neural%SED_highres(k)
! 						noise = 0.5d0 * djs(k)
! 						posterior = posterior + (djs(k) - neural%SED_highres(k))**2 / noise**2
! 					enddo
! 				else
! 					call lininterpol_eval(pars, database%SED_highres, prior%include_agn, prior%reddening_law)
! 					database%SED_highres = database%SED_highres / 1.e10 * 10.d0**pars(7)
! 					posterior = 0.d0
! 					do k = spectroscopy_from, spectroscopy_to
! 						djs(k) = noise * randomn() + database%SED_highres(k)
! 						noise = 0.5d0 * djs(k)
! 						posterior = posterior + (djs(k) - database%SED_highres(k))**2 / noise**2
! 					enddo
! 				endif
! 
! ! Compute predictive distribution as a Monte Carlo estimation. We sample
! ! from the posterior and compute the likelihood of the potential value of the new observation given
! ! by djs(j)
! 
! 				posterior = exp(-0.5d0 * posterior)
! 				predictive(j) = predictive(j) + posterior
! 
! 			enddo
! 
! 			predictive(j) = predictive(j) / (1.d0*N)
! 
! 		enddo
! 
! ! Compute the integral. We use the integral in log-scale to be more precise
! 		write(*,FMT='(A,F10.5)') 'Spectro EU=', -sum(log(predictive)) / (1.d0*M) - &
! 			0.5d0*log((2.d0*PI*exp(1.d0)*noise**2)**(spectroscopy_to-spectroscopy_from))
! 		write(13,*) 'Spectro   ', -sum(log(predictive)) / (1.d0*M) - 0.5d0*log((2.d0*PI*exp(1.d0)*noise**2)**(spectroscopy_to-spectroscopy_from))
		
		close(12)
		close(13)
		
	end subroutine information_gain

!------------------------------------------------------------
! Analyze the chains
!------------------------------------------------------------
	subroutine analyze_chains_multinest(chain_analysis, prior, suggest)
	type(markov_chain_type) :: chain_analysis	
	type(prior_type) :: prior
	logical :: suggest
	real(kind=KIND_FLOAT), pointer :: x(:), yGauss(:), yStep(:)
	integer :: i, j, n
	character(len=120) :: str_parameter
	real(kind=KIND_FLOAT) :: est, errup1s, errup2s, errdown1s, errdown2s, chi2_mean, chi2_max, pars(9), pctg
	integer, pointer :: indx(:)
	real(kind=KIND_FLOAT), allocatable :: chain_function(:)
	
		call read_chains_multinest(chain_analysis%filename, chain_analysis)
		
		open(unit=12,file=trim(adjustl(chain_analysis%filename))//'.hist1D',&
			action='write',status='replace')
		write(12,*) chain_analysis%nparam
		
		open(unit=13,file=trim(adjustl(chain_analysis%filename))//'.confidence',&
			action='write',status='replace')

! Do histograms for all inverted variables
		do i = 1, chain_analysis%nparam

			if (prior%typ(i) /= 'D') then
				call oned_histogram(chain_analysis%chain(i,:), x, yGauss, yStep)

! 1sigma confidence intervals			
				call conf_limits_ccf(x, yStep, one_sigma, est, errup1s, errdown1s)
				call conf_limits_ccf(x, yStep, two_sigma, est, errup2s, errdown2s)
			else
				allocate(x(10))
				x = prior%mu(i)
				allocate(yGauss(10))
				yGauss = 1.0
				allocate(yStep(10))
				yStep = 1.0
				est = prior%mu(i)
				errup1s = 0.d0
				errup2s = 0.d0
				errdown1s = 0.d0
				errdown2s = 0.d0
			endif				
			
			chain_analysis%most_probable(i) = est
			write(13,*) est, chain_analysis%best_parameters(i), errdown1s, errup1s, errdown2s, errup2s			
			str_parameter = chain_analysis%param_names(i)
			write(*,FMT='(A,A,F9.4,A7,F9.4)') trim(adjustl(str_parameter)), ' : E(x)=', est, ' - MAP=', &
				chain_analysis%best_parameters(i)
			write(*,FMT='(A,F9.4,A,F9.4)') '1-sigma : -', errdown1s, '  +', errup1s
			write(*,FMT='(A,F9.4,A,F9.4)') '2-sigma : -', errdown2s, '  +', errup2s
			n = size(x)
			write(12,*) i, n
			do j = 1, n
				write(12,*) x(j), yGauss(j), yStep(j)
			enddo
			deallocate(x)
			deallocate(yGauss)
			deallocate(yStep)
			
		enddo		

! Calculate average ln L for calculating Kullback-Leibler distance
		chain_analysis%avg_lnL = sum(chain_analysis%chain(chain_analysis%nparam+1,:)) / chain_analysis%nlength
		
		write(13,*) 'Evidence, <ln L> and Kullback-Leibler divergence'
		write(13,*) chain_analysis%evidence, chain_analysis%avg_lnL, &
			-chain_analysis%evidence + chain_analysis%avg_lnL

		write(*,*) 'Evidence, <ln L> and Kullback-Leibler divergence'
		write(*,*) chain_analysis%evidence, chain_analysis%avg_lnL, &
			-chain_analysis%evidence + chain_analysis%avg_lnL
			
		close(12)
		close(13)
		
		open(unit=12,file=trim(adjustl(chain_analysis%filename))//'.SED_samples',&
			action='write',status='replace',form='unformatted')
			
! The two last are the MAP SED and the estimated one
		write(12) chain_analysis%nlength+4
		
! Save file with evaluated SEDs at the posterior samples
		do i = 1, chain_analysis%nlength
			pars = chain_analysis%chain(1:9,i)

! Points obtained with new filter
 			if (what_to_do > 0) then
				call neural_eval(pars, neural%SED_highres, prior%include_agn, prior%reddening_law)
				neural%SED_highres = neural%SED_highres / 1.e10 * 10.d0**pars(7)
				write(12) neural%SED_highres
			else
				call lininterpol_eval(pars, database%SED_highres, prior%include_agn, prior%reddening_law)
				database%SED_highres = database%SED_highres / 1.e10 * 10.d0**pars(7)
				write(12) database%SED_highres
			endif
			
		enddo
		
! The SED for the value of the median parameters
! with and without extinction (if applied, of course).
! If no extinction is used, the two synthesis will be the same
		if (what_to_do > 0) then

! Full SED
			call neural_eval(chain_analysis%most_probable, neural%SED_highres, prior%include_agn, prior%reddening_law)
			neural%SED_highres = neural%SED_highres / 1.e10 * 10.d0**chain_analysis%most_probable(7)
			write(12) neural%SED_highres

! SED with (potential) AGN and no extinction
			call neural_eval(chain_analysis%most_probable, neural%SED_highres, prior%include_agn, 0)
			neural%SED_highres = neural%SED_highres / 1.e10 * 10.d0**chain_analysis%most_probable(7)
			write(12) neural%SED_highres

! SED without (potential) AGN and with extinction
			call neural_eval(chain_analysis%most_probable, neural%SED_highres, 0, prior%reddening_law)
			neural%SED_highres = neural%SED_highres / 1.e10 * 10.d0**chain_analysis%most_probable(7)
			write(12) neural%SED_highres

! SED without AGN and no extinction
			call neural_eval(chain_analysis%most_probable, neural%SED_highres, 0, 0)
			neural%SED_highres = neural%SED_highres / 1.e10 * 10.d0**chain_analysis%most_probable(7)
			write(12) neural%SED_highres
			
		else
			call lininterpol_eval(chain_analysis%most_probable, database%SED_highres, prior%include_agn, prior%reddening_law)
			database%SED_highres = database%SED_highres / 1.e10 * 10.d0**chain_analysis%most_probable(7)
			write(12) database%SED_highres
			
			call lininterpol_eval(chain_analysis%most_probable, database%SED_highres, prior%include_agn, 0)
			database%SED_highres = database%SED_highres / 1.e10 * 10.d0**chain_analysis%most_probable(7)
			write(12) database%SED_highres

			call lininterpol_eval(chain_analysis%most_probable, database%SED_highres, 0, prior%reddening_law)
			database%SED_highres = database%SED_highres / 1.e10 * 10.d0**chain_analysis%most_probable(7)
			write(12) database%SED_highres

			call lininterpol_eval(chain_analysis%most_probable, database%SED_highres, 0, 0)
			database%SED_highres = database%SED_highres / 1.e10 * 10.d0**chain_analysis%most_probable(7)
			write(12) database%SED_highres
		endif

! The SED for the value of the MAP parameters
		if (what_to_do > 0) then
			call neural_eval(chain_analysis%best_parameters, neural%SED_highres, prior%include_agn, prior%reddening_law)
			neural%SED_highres = neural%SED_highres / 1.e10 * 10.d0**chain_analysis%best_parameters(7)
			write(12) neural%SED_highres
			
			call neural_eval(chain_analysis%best_parameters, neural%SED_highres, prior%include_agn, 0)
			neural%SED_highres = neural%SED_highres / 1.e10 * 10.d0**chain_analysis%best_parameters(7)
			write(12) neural%SED_highres

			call neural_eval(chain_analysis%best_parameters, neural%SED_highres, 0, prior%reddening_law)
			neural%SED_highres = neural%SED_highres / 1.e10 * 10.d0**chain_analysis%best_parameters(7)
			write(12) neural%SED_highres

			call neural_eval(chain_analysis%best_parameters, neural%SED_highres, 0, 0)
			neural%SED_highres = neural%SED_highres / 1.e10 * 10.d0**chain_analysis%best_parameters(7)
			write(12) neural%SED_highres
		else
			call lininterpol_eval(chain_analysis%best_parameters, database%SED_highres, prior%include_agn, prior%reddening_law)
			database%SED_highres = database%SED_highres / 1.e10 * 10.d0**chain_analysis%best_parameters(7)
			write(12) database%SED_highres
			
			call lininterpol_eval(chain_analysis%best_parameters, database%SED_highres, prior%include_agn, 0)
			database%SED_highres = database%SED_highres / 1.e10 * 10.d0**chain_analysis%best_parameters(7)
			write(12) database%SED_highres

			call lininterpol_eval(chain_analysis%best_parameters, database%SED_highres, 0, prior%reddening_law)
			database%SED_highres = database%SED_highres / 1.e10 * 10.d0**chain_analysis%best_parameters(7)
			write(12) database%SED_highres

			call lininterpol_eval(chain_analysis%best_parameters, database%SED_highres, 0, 0)
			database%SED_highres = database%SED_highres / 1.e10 * 10.d0**chain_analysis%best_parameters(7)
			write(12) database%SED_highres
		endif

		close(12)

! Carry out estimation of best next filter
		if (suggest) call information_gain(chain_analysis, prior)
		
	end subroutine analyze_chains_multinest

end module chain_analysis_mod
