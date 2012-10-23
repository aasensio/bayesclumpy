module global_vars_mod
use types_mod
implicit none
					
	type(neural_network_type) :: neural
	type(markov_chain_type) :: chain
	type(prior_type) :: prior
	type(full_database_type) :: database
			
	integer :: verbose_level, what_to_do
	real(kind=KIND_FLOAT) :: version_bc

	real(kind=KIND_FLOAT), parameter :: one_sigma = 68.268955e0, &
		two_sigma = 95.449972e0, PI=3.14159265359d0
	
end module global_vars_mod