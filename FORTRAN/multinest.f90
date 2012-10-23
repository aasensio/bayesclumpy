! Parameters for the multinest algorithm
module multinest_mod
use global_vars_mod
use Nested
use mcmc_mod, only : calculate_posterior
use chain_analysis_mod
implicit none

!dimensionality
	integer :: sdim
      
      
!priors on the parameters
!uniform priors (-6,6) are used for all dimensions & are set in main.f90
   real(kind=8), allocatable :: prior_range(:,:) !spriorran(sdim,2)

   integer, allocatable :: nest_pWrap(:) ! Wraparound parameters
      
! Parameters for Nested Sampler
	
!whether to do multimodal sampling
	logical :: nest_mmodal
 	parameter(nest_mmodal=.true.)
	
!max no. of live points
   integer :: nest_nlive
	parameter(nest_nlive=500)
      
!tot no. of parameters, should be sdim in most cases but if you need to
!store some additional parameters with the actual parameters then
!you need to pass them through the likelihood routine
	integer :: nest_nPar

	logical :: nest_ceff
	parameter(nest_ceff = .FALSE.)
      
!seed for nested sampler, -ve means take it from sys clock
	integer :: nest_rseed
	parameter(nest_rseed=-1)
      
!evidence tolerance factor
   real(kind=8) :: nest_tol
   parameter(nest_tol=0.5)
      
!enlargement factor reduction parameter
	real(kind=8) :: nest_efr
   parameter(nest_efr=0.8d0)
      
!root for saving posterior files
	character*100 nest_root	
	
!no. of iterations after which the ouput files should be updated
	integer nest_updInt
	parameter(nest_updInt=100)

!null evidence (set it to very high negative no. if null evidence is unknown)
	real*8 nest_Ztol
	parameter(nest_Ztol=-1.d90)
      
!max modes expected, for memory allocation
  	integer nest_maxModes
  	parameter(nest_maxModes=20)
      
!no. of parameters to cluster (for mode detection)
  	integer nest_nClsPar  	
      
!whether to resume from a previous run
  	logical nest_resume
  	parameter(nest_resume=.false.)
      
!feedback on the sampling progress?
  	logical nest_fb
  	parameter(nest_fb=.true.)

	real(kind=8) :: maxlhood
	
contains

!----------------------------------------------------------------------
	subroutine nest_Sample
   integer :: nclusters,context !total number of clusters found
   integer :: maxNode !variables used by the posterior routine
   
   	call nestRun(nest_mmodal,nest_ceff,nest_nlive,nest_tol,nest_efr,sdim,nest_nPar, &
   		nest_nClsPar,nest_maxModes,nest_updInt,nest_Ztol,nest_root,nest_rseed, nest_pWrap, &
   		nest_fb,nest_resume,getLogLike,dumper,context)

	end subroutine nest_Sample
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! Dumper is called after updInt*10 iterations
!----------------------------------------------------------------------
	subroutine dumper(nSamples, nlive, nPar, physLive, posterior, paramConstr, maxLogLike, logZ)

	integer :: nSamples                                ! number of samples in posterior array
	integer :: nlive                                   ! number of live points
	integer :: nPar                                    ! number of parameters saved (physical plus derived)
	real(kind=8), pointer :: physLive(:,:)      ! array containing the last set of live points
	real(kind=8), pointer :: posterior(:,:)     ! array with the posterior distribution
	real(kind=8), pointer :: paramConstr(:)     ! array with mean, sigmas, maxlike & MAP parameters
	real(kind=8) :: maxLogLike                     ! max loglikelihood value
	real(kind=8) :: logZ                           ! log evidence

	end subroutine dumper

!----------------------------------------------------------------------
! Wrapper around Likelihood Function
! Cube(1:n_dim) has nonphysical parameters
! scale Cube(1:n_dim) & return the scaled parameters in Cube(1:n_dim) &
! additional parameters in Cube(n_dim+1:nPar)
! return the log-likelihood in lnew
!----------------------------------------------------------------------
	subroutine getLogLike(Cube,n_dim,nPar,lnew,context)

   integer :: n_dim,nPar,context
   real(kind=8) :: lnew,Cube(nPar)
   
   !call your loglike function here   
   	call slikelihood(Cube,lnew)

	end subroutine getLogLike

!----------------------------------------------------------------------
! Evaluate likelihood
!----------------------------------------------------------------------
	subroutine slikelihood(Cube,slhood)              
	real(kind=8) :: Cube(nest_nPar),slhood
	real(kind=8) :: temp(sdim),dist,loclik, chi2, ff_total
	integer :: i,j
	logical :: physical
	

		Cube(1:sdim) = (prior_range(1:sdim,2)-prior_range(1:sdim,1))*Cube(1:sdim) + prior_range(1:sdim,1)

		chain%proposed = Cube

		slhood = calculate_posterior(chain%proposed, prior, chain%observation, chain%typ)

		if (slhood > maxlhood) then
			maxlhood = slhood
			chain%best_parameters = chain%proposed
			
! Save the best parameters for later use
			open(unit=16,file=trim(adjustl(nest_root))//'.MAP',action='write',status='replace')
			write(16,*) chain%best_parameters
			close(16)
		endif

		return
			
	end subroutine slikelihood


!----------------------------------------------------------------------
! Do MULTINEST MCMC
!----------------------------------------------------------------------
	subroutine do_multinest(chain, prior, step_from, step_to, initialize)
	type(markov_chain_type) :: chain
	type(prior_type) :: prior
	type(observations_type) :: observations	
	integer :: step_from, step_to
	logical :: initialize
	
	real(kind=8) :: chi2
	logical :: physical
	integer :: i

		nest_nPar = chain%nparam
		sdim = nest_npar
		nest_nClsPar = nest_npar

		nest_root = chain%filename

		maxlhood = -1.d100

! Allocate memory for the prior ranges
		allocate(prior_range(nest_nPar,2))

! For periodic boundary conditions
		allocate(nest_pWrap(sdim))
		nest_pWrap = 0

! Set them equal to the prior ranges read before
		prior_range(:,1) = prior%minim(:)
		prior_range(:,2) = prior%maxim(:)


! Sample
		call nest_Sample

! ! Save MAP profile
! 		stat%parameters_proposed = stat%parameters_MAP
!  		call fill_model(inv,stat,proposed)
!  
!  		call synthesize(proposed,linea,Observation,Stokes_Syn(1))
!  		
!  		open(unit=12,file=trim(adjustl(fich_chain))//'.stokes_map',&
!  			action='write',status='replace')
!  		do i = 1, Stokes_Syn(1)%nlambda
!  			write(12,FMT='(F10.4,2X,4(E12.4))') Stokes_Syn(1)%lambda(i),&
!  				Stokes_Syn(1)%stokes(1,i), Stokes_Syn(1)%stokes(2,i), &
!  				Stokes_Syn(1)%stokes(3,i), Stokes_Syn(1)%stokes(4,i)
!  		enddo
!  		close(12)

! Analyze chains
! 		call analyze_chains_multinest(chain_analysis,fich_chain,stat%parameters_MAP)
! 
! 		write(*,*) 'MAP values'
! 		write(*,*) stat%parameters_MAP
! 		write(*,*)
! 
!  		print *, 'Maximum likelihood value = ', maxlhood 		

		deallocate(prior_range)

		open(unit=18,file='end.info',action='write',status='replace')
		write(18,*) 'Multinest finished'
		close(18)
		
	end subroutine do_multinest

end module multinest_mod
