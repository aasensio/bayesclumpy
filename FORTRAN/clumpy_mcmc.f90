program mcmc_ann
use ann_mod
use mcmc_mod
use global_vars_mod
use chain_analysis_mod
use multinest_mod, only : do_multinest
use maths_mod, only : detect_version
implicit none

real(kind=4) :: dum
character(len=200) :: inputfile
integer :: nargs

! Detect if spectroscopy is allowed in this version
	version_bc = detect_version()	
	
	call init_random_seed()
	
	call initialize_neural_networks()
		
	inputfile = 'chain.cfg'

! Verify if a configuration file is passed as an argument
! If not, use the default 'config' file
	nargs = iargc()
   if (nargs == 1) then
   	call getarg(1,inputfile)
   endif

   call initialize_problem(inputfile, chain, prior, .TRUE.)
   
! If we want to do inference with the full database and local linear interpolation,
! read the database
   if (what_to_do < 0) then
		call initialize_database()
	endif
   
! Carry out MCMC

	open(unit=18,file='end.info',action='write',status='replace')
	write(18,*) 'Sampling'
	close(18)
	select case(abs(what_to_do))
		case(1)			
			call mcmc(chain, prior, 1, chain%niter_max, .TRUE.)
			call analyze_chains(chain, prior)
		case(2)
 			call do_multinest(chain, prior, 1, chain%niter_max, .TRUE.)
 			call analyze_chains_multinest(chain, prior, .FALSE.)
 		case(3)
 			call analyze_chains_multinest(chain, prior, .TRUE.)
 		case(4)
 			call do_multinest(chain, prior, 1, chain%niter_max, .TRUE.)
 			call analyze_chains_multinest(chain, prior, .FALSE.)
	end select
				
end program mcmc_ann