module mcmc_mod
use types_mod
use maths_mod
use ann_mod, only : neural_eval_filter, lininterpol_eval_filter
use global_vars_mod, only : PI, verbose_level, what_to_do, version_bc
implicit none
contains

!-------------------------------------------------------------------
! Initialize Markov Chain
!-------------------------------------------------------------------
	subroutine initialize_chain(chain, nparam, nitermax, typ)
	type(markov_chain_type) :: chain
	integer :: nparam, nitermax	
	character(len=4) :: typ
	
		chain%nparam = nparam
		chain%niter_max = nitermax
		
		allocate(chain%proposed(nparam))
		allocate(chain%best_parameters(nparam))
		allocate(chain%last(nparam))
		allocate(chain%most_probable(nparam))
		
		allocate(chain%chain(nparam,nitermax))
		chain%chain = 0.0
		allocate(chain%mean(nparam))
		chain%mean = 0.0
		allocate(chain%mean_old(nparam))
		chain%mean_old = 0.0
		allocate(chain%covariance(nparam,nparam))
		chain%covariance = 0.0
		
		allocate(chain%posterior(nitermax))
		allocate(chain%acceptance_rate(nitermax))
		
		chain%typ = typ
		chain%alpha = 1.d0

		allocate(chain%param_names(nparam))
		
	end subroutine initialize_chain

!-------------------------------------------------------------------
! Read one of the filters and return the central wavelength
!-------------------------------------------------------------------
	function read_filter(filter_name, filter_info)
	real(kind=KIND_FLOAT) :: read_filter
	character(len=15) :: filter_name
	type(filter_type) :: filter_info
	integer, allocatable :: ind(:)
	real(kind=KIND_FLOAT), allocatable :: trans(:)
	integer :: i

		open(unit=15,file='FILTERS/'//trim(adjustl(filter_name))//'.res',action='read',status='old')
		read(15,*) filter_info%Nlambdas, filter_info%normalization

		allocate(filter_info%lambda(filter_info%Nlambdas))		
		allocate(filter_info%transmission(filter_info%Nlambdas))
		allocate(trans(filter_info%Nlambdas))
		allocate(filter_info%SED_interpol(filter_info%Nlambdas))

		do i = 1, filter_info%Nlambdas
			read(15,*) filter_info%lambda(i), trans(i)
			
! Avoid negative and very small values in the filters (putting them to zero does not
! affect much the final result)
			if (filter_info%transmission(i) < 1.e-5) then
				filter_info%transmission(i) = 0.0
			endif
			
		enddo

! Reorder the filter in case the wavelength is in decreasing order
		allocate(ind(filter_info%Nlambdas))
		call qsortd(filter_info%lambda,ind,filter_info%Nlambdas)
		do i = 1, filter_info%Nlambdas
			filter_info%transmission(i) = trans(ind(i))
		enddo
		deallocate(ind)

! Normalize
		filter_info%transmission = filter_info%transmission / filter_info%normalization

! Locate central wavelength of the filter
		filter_info%central = int_tabulated(filter_info%lambda, filter_info%lambda*filter_info%transmission) /&
			int_tabulated(filter_info%lambda, filter_info%transmission)
		
		close(15)

		filter_info%name = filter_name

		read_filter = filter_info%central
		
	end function read_filter

	
!-------------------------------------------------------------------
! Initialize observations and filters used
!-------------------------------------------------------------------
	subroutine initialize_observations(observation, obsfile)
	type(observations_type) :: observation
	character(len=200) :: obsfile, filter
	integer :: i

		verbose_level = 2
		if (verbose_level > 1) &
			write(*,*) 'Reading observations file : ', trim(adjustl(obsfile))
		open(unit=13,file=trim(adjustl(obsfile)),action='read',status='old')
		read(13,*) observation%description

! Version with spectroscopy
		if (version_bc < 0) then
			read(13,*) observation%npoints, observation%npoints_spectrum
! Total number of points (including filters and spectrum)
			observation%npoints_total = observation%npoints + observation%npoints_spectrum
			if (verbose_level > 1) then
				write(*,*) 'Number of filter points : ', observation%npoints
				write(*,*) 'Number of spectrum points : ', observation%npoints_spectrum
				write(*,*) 'Number of total points : ', observation%npoints+observation%npoints_spectrum
			endif
			allocate(observation%obs_x(observation%npoints_total))
			allocate(observation%obs_y(observation%npoints_total))
			allocate(observation%obs_sigma(observation%npoints_total))
			allocate(observation%model_eval(observation%npoints_total))
			allocate(observation%filter(observation%npoints))

		else
			read(13,*) observation%npoints
			if (verbose_level > 1) then
				write(*,*) 'Number of filter points : ', observation%npoints
			endif
			allocate(observation%obs_x(observation%npoints))
			allocate(observation%obs_y(observation%npoints))
			allocate(observation%obs_sigma(observation%npoints))
			allocate(observation%model_eval(observation%npoints))
			allocate(observation%filter(observation%npoints))
		endif
		
! Read the SED that is given by filters
		do i = 1, observation%npoints
			read(13,*) filter, observation%obs_y(i), observation%obs_sigma(i)
			
			observation%obs_x(i) = read_filter(filter,observation%filter(i))			
			if (verbose_level > 1) &
			write(*,*) 'Filter ', trim(adjustl(filter)), &
				' : N_lambda = ', observation%filter(i)%Nlambdas,&
				' / Lambda_c = ', observation%filter(i)%central
			
! Identify points with only upper limits. They are put as negative values
! followed by the confidence limit
			if (observation%obs_y(i) < 0.0) then
			
! The value of sigma is such that, given a confidence level (0.68, 0.95, etc.)
! the integral between 0 and the tabulated value is 0.68, 0.95, etc. times
! the total area of the gaussian
				observation%obs_sigma(i) = &
					abs(observation%obs_y(i)) / inverse_normal(observation%obs_sigma(i))
					
! Center the likelihood for these points at zero
				observation%obs_y(i) = 0.d0
				
				if (verbose_level > 1) then
					write(*,FMT='(A,F10.4)') '       Upper limit -> sigma set to ',&
						observation%obs_sigma(i)
				endif
				
			else
				if (verbose_level > 1) then
					write(*,FMT='(A,F10.4)') '       Normal -> sigma is ',	observation%obs_sigma(i)
				endif
			endif
			
		enddo

! Version with spectroscopy
		if (version_bc < 0 .and. observation%npoints_spectrum /= 0) then
			read(13,*)
! Read the rest of the SED that is given by a spectrum
			do i = observation%npoints+1, observation%npoints_total
				read(13,*) observation%obs_x(i), observation%obs_y(i), observation%obs_sigma(i)
			enddo
		endif

		close(13)
		
	end subroutine initialize_observations

!-------------------------------------------------------------------
! Initialize priors
!-------------------------------------------------------------------
	subroutine initialize_prior(chain, prior, ranges_from, ranges_to, typ, mu, sigma)
	type(markov_chain_type) :: chain
	type(prior_type) :: prior
	real(kind=KIND_FLOAT) :: ranges_from(:), ranges_to(:)
	character(len=1) :: typ(:)	
	real(kind=KIND_FLOAT) :: mu(:), sigma(:)
	integer :: i
		
		prior%nparam = chain%nparam
		allocate(prior%typ(chain%nparam))
		allocate(prior%mu(chain%nparam))
		allocate(prior%sigma(chain%nparam))
		
! Limits for the parameters
		allocate(prior%minim(prior%nparam))
		prior%minim = ranges_from
		allocate(prior%maxim(prior%nparam))
		prior%maxim = ranges_to		
				
! Fill information and if the prior is gaussian, fill in
! the width and position of the gaussian
		do i = 1, prior%nparam
			prior%typ(i) = typ(i)
			
! Gaussian type
			if (typ(i) == 'G') then
				prior%mu(i) = mu(i)
				prior%sigma(i) = sigma(i)
			endif
			
! Dirac type
			if (typ(i) == 'D') then
				prior%mu(i) = mu(i)
				prior%minim(i) = mu(i) - 1.d-3
				prior%maxim(i) = mu(i) + 1.d-3
			endif			
		enddo

	end subroutine initialize_prior

!-------------------------------------------------------------------
! Initialize problem
!-------------------------------------------------------------------
	subroutine initialize_problem(inputfile, chain, prior, read_prior)
	type(markov_chain_type) :: chain
	type(prior_type) :: prior
	real(kind=KIND_FLOAT) :: temp
	logical :: read_prior
	character(len=200) :: inputfile, typ_chain, outfile, obsfile
	integer :: nitermax, nparam, i
	real(kind=KIND_FLOAT), allocatable :: ranges_from(:), ranges_to(:), mu(:), sigma(:)
	character(len=1), allocatable :: typ(:)
	character(len=20), allocatable :: param_names(:)
	
	
		open(unit=12,file=inputfile,action='read',status='old')		
		call lb(12,2)
		read(12,*) what_to_do
		call lb(12,1)
		read(12,*) nparam
		call lb(12,1)
		read(12,*) nitermax
		call lb(12,1)
		read(12,*) chain%burnin
		call lb(12,1)
		read(12,*) typ_chain

		if (typ_chain == 'MCMC') then
			call lb(12,1)
			read(12,*) obsfile
			call initialize_observations(chain%observation, obsfile)
		else
			call lb(12,2)
		endif
		
		call lb(12,1)
		read(12,*) outfile
		chain%filename = outfile
		
		call initialize_chain(chain, nparam, nitermax, typ_chain)

! Use AGN+SED or only SED
		call lb(12,1)
		read(12,*) prior%include_agn

! Which reddening law to use
		call lb(12,1)
		read(12,*) prior%reddening_law

! Read Chiar & Tielens law if used
		if (prior%reddening_law == 6) then
			write(*,*) 'Reading Chiar & Tielens extinction curve'
			open(unit=14,file='FILTERS/extinction_chiar_tielens2006.dat',action='read',status='old')
			call lb(14, 14)
			allocate(prior%lambda_pixie(258), prior%extinction_pixie(258))
			do i = 1, 258
				read(14,*) prior%lambda_pixie(i), temp, prior%extinction_pixie(i)
			enddo
			close(14)			
		endif
		
		allocate(ranges_from(nparam))
		allocate(ranges_to(nparam))
		allocate(typ(nparam))
		allocate(mu(nparam))
		allocate(sigma(nparam))
		allocate(param_names(nparam))
		
		call lb(12,3)
		do i = 1, nparam
			read(12,*) param_names(i), ranges_from(i), ranges_to(i), typ(i), mu(i), sigma(i)
		enddo

! Names given to the parameters for later identification
		chain%param_names = param_names
		
		if (read_prior) then
			call initialize_prior(chain, prior, ranges_from, ranges_to, typ, mu, sigma)
		endif
		
		deallocate(ranges_from)
		deallocate(ranges_to)
		deallocate(typ)
		deallocate(mu)
		deallocate(sigma)
		
		close(12)
	end subroutine initialize_problem

!-------------------------------------------------------------------
! Update the mean and the covariance
!-------------------------------------------------------------------
	subroutine update_statistics(loop, chain)
	integer :: loop
	type(markov_chain_type) :: chain
	integer :: i, j
	
! Update the mean
		chain%mean_old = chain%mean
		chain%mean = chain%mean_old + (chain%last - chain%mean_old) / (loop+1.d0)
		
! Update the covariance
		do i = 1, chain%nparam
			do j = 1, chain%nparam
				chain%covariance(i,j) = (loop-1.d0)/loop * chain%covariance(i,j) + &
					(chain%last(i)-chain%mean_old(i))*(chain%last(j)-chain%mean_old(j)) / &
					(loop+1.d0)**2 + (chain%last(i)-chain%mean(i))*&
					(chain%last(j)-chain%mean(j)) / (loop*1.d0)
			enddo
		enddo
	
	end subroutine update_statistics

!-------------------------------------------------------------------
! Multi-dimensional gaussian proposal distribution function
!-------------------------------------------------------------------
	function proposal(last, covariance, alpha)
	real(kind=KIND_FLOAT) :: last(:), covariance(:,:), alpha
	real(kind=KIND_FLOAT) :: proposal(size(last))

		proposal = mrandomn(last,alpha*covariance)
		
	end function proposal
	
!-------------------------------------------------------------------
! Calculate prior distribution
!-------------------------------------------------------------------
	function calculate_prior(last, prior)
	real(kind=KIND_FLOAT) :: last(:), calculate_prior
	type(prior_type) :: prior
	integer :: i
	
		calculate_prior = 1.0
		do i = 1, prior%nparam

! Verify Dirac delta priors
! In this case, set last(i) to the value of mu
			if (prior%typ(i) == 'D') then
				last(i) = prior%mu(i)
			endif
			
! If the value of any parameter is outside the available range
! return 0 for the prior
			if (last(i) > prior%maxim(i) .or. last(i) < prior%minim(i)) then
				calculate_prior = 0.0
				return
			endif
			
! Now verify possible Gaussian priors
			if (prior%typ(i) == 'G') then
				calculate_prior = calculate_prior * &
					exp(-(last(i)-prior%mu(i))**2 / (2.0*prior%sigma(i)**2))
			endif
		enddo
				
	end function calculate_prior

!-------------------------------------------------------------------
! Calculate posterior distribution
!-------------------------------------------------------------------
	function calculate_posterior(last, prior, observations, typ)
	real(kind=KIND_FLOAT) :: last(:), calculate_posterior
	type(prior_type) :: prior
	character(len=4) :: typ
	type(observations_type) :: observations
	real(kind=KIND_FLOAT) :: P_prior, logP_likelihood, logP_posterior, chi2, chi2_phot, chi2_spec

! Calculate prior distribution
		P_prior = calculate_prior(last, prior)

		if (P_prior == 0.0) then
			calculate_posterior = 0.0
			return
		else
			select case(typ)

! Standard MCMC with gaussian errors
				case('MCMC')

! Evaluate the direct problem					
! If what_to_do > 0, use the neural network
! If what_to_do < 0, use the linear interpolation in the database
					if (what_to_do > 0) then
						call neural_eval_filter(last, observations, prior%include_agn, prior%reddening_law)
					else
						call lininterpol_eval_filter(last, observations, prior%include_agn, prior%reddening_law)
					endif

! Weight differently both datasets according to Lahav et al. MNRAS 315, L45 (2000)
					if (abs(what_to_do) == 4) then
! chi^2 = Na * ln(chi_B^2) + Nb * ln(chi_B^2)

! First photometric points
						chi2_phot = sum((observations%model_eval(1:observations%npoints)-&
							observations%obs_y(1:observations%npoints))**2 / &
							(observations%obs_sigma(1:observations%npoints)**2))

! Then spectroscopic points
						chi2_spec = sum((observations%model_eval(observations%npoints+1:observations%npoints_total)-&
							observations%obs_y(observations%npoints+1:observations%npoints_total))**2 / &
							(observations%obs_sigma(observations%npoints+1:observations%npoints_total)**2))

						chi2 = observations%npoints * log(chi2_phot) + &
							(observations%npoints_total - observations%npoints) * log(chi2_spec)
					else
! Calculate the standard chi^2
						chi2 = sum((observations%model_eval-observations%obs_y)**2 / &
							(observations%obs_sigma**2))
					endif

! Calculate likelihood and posterior
					logP_likelihood = -0.5 * chi2
					logP_posterior = logP_likelihood + log(P_prior)
					
				case('SAMP')
				
				case('PRIO')
					logP_posterior = log(P_prior)
			end select
		endif
		
		calculate_posterior = logP_posterior
		
	end function calculate_posterior

!-------------------------------------------------------------------
! Markov Chain Monte Carlo
!-------------------------------------------------------------------
	subroutine mcmc(chain, prior, step_from, step_to, initialize)
	type(markov_chain_type) :: chain
	type(prior_type) :: prior
	type(observations_type) :: observations	
	integer :: step_from, step_to
	logical :: initialize
	integer :: i, loop, j, init_model
	real(kind=KIND_FLOAT) :: fproposed, flast, fmax, r, alpha, ran
	real(kind=KIND_FLOAT) :: acceptance_rate_global, acceptance_rate_local, alpha_theory
	real(kind=KIND_FLOAT), allocatable :: list_accepted(:), lhs_pos(:,:)
	character(len=200) :: format_output
	character(len=2) :: str

		allocate(list_accepted(200))
		verbose_level = 2
								
		if (verbose_level > 1) &
			write(*,*) 'Running MCMC for ', chain%nparam,' parameters...'
		alpha_theory = 2.4**2 / chain%nparam
		
		write(str,FMT='(I2)') chain%nparam
		format_output = '(4X,A7,'//trim(adjustl(str))//'(F9.4,1X))'
		
! Initialize chain
		if (initialize) then
		
! Starting point is selected randomnly from the available space
			do i = 1, chain%nparam
				chain%proposed(i) = randomu() * (prior%maxim(i) - prior%minim(i)) + &
					prior%minim(i)
				chain%last(i) = chain%proposed(i)
			enddo
			
			chain%mean = chain%proposed
			
! Covariance matrix is set initially to a diagonal equal to 10% the
! available space of parameters
			do i = 1, chain%nparam
				chain%covariance(i,i) = 0.1 * (prior%maxim(i) - prior%minim(i))
			enddo
			
! Evaluate the posterior at this point
			flast = calculate_posterior(chain%proposed, prior, chain%observation, &
				chain%typ)
			
			fmax = flast

			chain%accepted_models = 0.0

			print *, 'LHS initialization'
! Sample in a latin hypercube sampling of the space of parameters and start
! as close to the maximum as possible to allow efficient posterior sampling
			allocate(lhs_pos(chain%nparam,1000))

! Fill the array with uniform random numbers
			do i = 1, 1000
				do j = 1, chain%nparam
					lhs_pos(j,i) = randomu()
				enddo
			enddo

! Latinize the sampling
			call dtable_latinize (chain%nparam, 1000, lhs_pos)

! Evaluate posterior in each point and save the one with the largest for initialization
			fmax = -1.d30
			do i = 1, 1000
				do j = 1, chain%nparam
					chain%proposed(j) = lhs_pos(j,i) * (prior%maxim(j) - prior%minim(j)) + &
						prior%minim(j)				
				enddo
				flast = calculate_posterior(chain%proposed, prior, chain%observation, &
					chain%typ)
				if (flast > fmax) then
					fmax = flast
					init_model = i
				endif
			enddo

! Select the initial model
			do j = 1, chain%nparam
				chain%proposed(j) = lhs_pos(j,init_model) * (prior%maxim(j) - prior%minim(j)) + &
					prior%minim(j)
			enddo
			deallocate(lhs_pos)
			
		else
			flast = chain%flast
			fmax = chain%fmax
		endif
		
! Start/continue the chain
		do loop = step_from, step_to
		
! Propose a model until it is inside the range available for the parameters			
			chain%proposed = proposal(chain%last, chain%covariance, chain%alpha)			
			fproposed = calculate_posterior(chain%proposed, prior, chain%observation, chain%typ)
			
			list_accepted(mod(loop,200)+1) = 0
			
			if (fproposed /= 0.0) then
				r = exp(fproposed - flast)
! 				r = fproposed / flast
				alpha = min(1.0,r)
				ran = randomu()
         	if (ran < alpha) then
	         	chain%last = chain%proposed
         		flast = fproposed
         		chain%flast = flast
         	
         		if (flast > fmax) then
	         		fmax = flast
	         		chain%fmax = flast
	         		chain%best_parameters = chain%last
         		endif
	         	
         		chain%accepted_models = chain%accepted_models + 1
         		list_accepted(mod(loop,200)+1) = 1
         	
         	endif
         endif
         
         chain%chain(:,loop) = chain%last
         chain%acceptance_rate(loop) = acceptance_rate_local
         chain%posterior(loop) = flast
         
! Update statistical information
         if (loop > 100) then
         	call update_statistics(loop, chain)
         endif
         
         if (loop > 100 .and. loop / 100 == loop / 100.0) then
				acceptance_rate_global = chain%accepted_models / (1.0*loop) * 100.0
				acceptance_rate_local = sum(list_accepted) / 200.0 * 100.0

				if (verbose_level > 1) then
					write(*,FMT='(I6,2X,A,F6.2,A,F6.2,A,F9.3,1X,F9.3,A,E13.4)'),&
						loop, 'Acc. rate: ', &
						acceptance_rate_global, '-', acceptance_rate_local, &
						'% -- (alp,alp/th): ', chain%alpha, &
						chain%alpha/alpha_theory, ' -- logL(max): ', fmax
						
					write(*,FMT=format_output) 'Avg:   ', chain%mean
					write(*,FMT=format_output) 'Sigma:   ', sqrt(diagon(chain%covariance))
				endif
			
			
! Control the value of alpha by forcing the acceptance rate
				if (acceptance_rate_local > 30.d0 .and. chain%alpha < 1.d0) then
					chain%alpha = chain%alpha * 1.04d0			
				endif
				if (acceptance_rate_local < 20.d0 .and. chain%alpha > 0.1) then
					chain%alpha = chain%alpha * 0.96d0
				endif

				if (acceptance_rate_global < 1.d-2 .and. loop > step_to / 3) then
					print *, 'Sampler not working properly. Rerun.'
					stop
				endif
				
			endif

		enddo

		deallocate(list_accepted)

		open(unit=18,file='end.info',action='write',status='replace')
		write(18,*) 'MCMC finished'
		close(18)
				
	end subroutine mcmc

end module mcmc_mod
