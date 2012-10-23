module types_mod
implicit none

! Precision of the calculations
	integer, parameter :: KIND_FLOAT=4
	real(kind=KIND_FLOAT), parameter :: very_large = 1.e30, one = 1.e0
	
	type neural_network_PCA_type
		integer :: ninput, nhidden, noutput, nunknowns, PCAsize
      real(kind=KIND_FLOAT), pointer :: input(:), hidden(:), bias(:), output(:)
      real(kind=KIND_FLOAT), pointer :: hidden_output(:,:), input_hidden(:,:)
      real(kind=KIND_FLOAT), pointer :: input_norm(:), output_norm(:)
      real(kind=KIND_FLOAT), pointer :: input_mean(:), output_mean(:)
      real(kind=KIND_FLOAT), pointer :: PCAvector(:)
	end type neural_network_PCA_type
	
	type neural_network_type
		integer :: Nnets
		real(kind=KIND_FLOAT), pointer :: lambda(:), meanSpec(:), SED_highres(:)
		type(neural_network_PCA_type), pointer :: net(:)
	end type neural_network_type

	type filter_type
		integer :: Nlambdas
		real(kind=KIND_FLOAT) :: central, normalization
		character(len=15) :: name
		real(kind=KIND_FLOAT), pointer :: lambda(:), transmission(:), SED_interpol(:)
	end type filter_type
	
	type observations_type
		integer :: npoints, npoints_spectrum, npoints_total
		character(len=200) :: description
		real(kind=KIND_FLOAT), pointer :: obs_x(:), obs_y(:), obs_sigma(:), model_eval(:)
		type(filter_type), pointer :: filter(:)
	end type observations_type
	
	type prior_type
		integer :: nparam, include_agn, reddening_law
		character(len=1), pointer :: typ(:)
		real(kind=KIND_FLOAT), pointer :: sigma(:), mu(:), minim(:), maxim(:), lambda_pixie(:), extinction_pixie(:)
	end type prior_type
	
	type markov_chain_type
		character(len=4) :: typ
		integer :: nparam, accepted_models, initial, niter_max, nlength
		real(kind=KIND_FLOAT), pointer :: proposed(:), last(:), most_probable(:)
		real(kind=KIND_FLOAT), pointer :: chain(:,:), mean(:), mean_old(:), covariance(:,:)
		real(kind=KIND_FLOAT), pointer :: posterior(:), acceptance_rate(:), best_parameters(:)
		real(kind=KIND_FLOAT) :: flast, fmax, alpha, burnin, avg_lnL, evidence
		character(len=200) :: filename
		character(len=20), pointer :: param_names(:)
		type(observations_type) :: observation
	end type markov_chain_type
	
	type parallel_tempering_type
		integer :: nchains, avg_steps, niter_max
		real(kind=KIND_FLOAT), pointer :: betas(:)
		type(markov_chain_type), pointer :: chain(:)
	end type parallel_tempering_type
	
	type full_database_type
		integer :: nY, nsig, nN0, nq, ntauv, ni, nlambda, db_id
		integer :: ny_id, nsig_id, nn0_id, nq_id, ntauv_id, ni_id, npca_id, nlam_id
		integer :: lam_id, coefs_id, base_id, meansed_id
		integer :: y_id, sig_id, n0_id, q_id, tauv_id, i_id, npca, nlam
      real(kind=KIND_FLOAT), pointer :: Y(:), sig(:), N0(:), q(:), tauv(:), i(:)
      real(kind=KIND_FLOAT), pointer :: coefs(:,:,:,:,:,:,:), base(:,:)
      real(kind=KIND_FLOAT), pointer :: meanSed(:), lambda(:), wrk(:,:), params(:,:), SED_highres(:)
      integer, pointer :: indi(:), ee(:,:)
	end type full_database_type
	
end module types_mod
