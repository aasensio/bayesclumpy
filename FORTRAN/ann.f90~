module ann_mod
use types_mod
use global_vars_mod
use maths_mod, only : int_tabulated, lin_interpol, hunt
implicit none

contains

!-------------------------------------------------------------------
! Return the AGN spectrum
!-------------------------------------------------------------------
	subroutine agn_shape(lambda, agn)
	real(kind=KIND_FLOAT) :: lambda(:), agn(:)
	real(kind=KIND_FLOAT), parameter :: lambdah = 0.01, lambdau = 0.1, lambdaRJ = 1.0, p = 0.5, cte = 0.2784

		where(lambda <= lambdah)
			agn = cte * lambda**1.2 / lambdah**1.2
		endwhere

		where(lambda > lambdah .and. lambda <= lambdau)
			agn = cte
		endwhere

		where(lambda > lambdau .and. lambda <= lambdaRJ)
			agn = cte * lambda**(-p) / lambdau**(-p)
		endwhere

		where(lambda > lambdaRJ)
			agn = cte * lambdaRJ**(-p) / lambdau**(-p) * lambda**(-3) / lambdaRJ**(-3)
		endwhere
		
	end subroutine agn_shape

! ---------------------------------------------------------
! Extinction laws
! ---------------------------------------------------------
	function extinction_curve(lambda, model, Av)
	real(kind=KIND_FLOAT) :: lambda(:), extinction_curve(size(lambda)), Av
	integer :: model, n
	real(kind=KIND_FLOAT) :: invlambda0, gamm, C1, C2, C3, C4, slope, zero, Rv
	real(kind=KIND_FLOAT), allocatable :: invlambda(:)
	
		n = size(lambda)
		
		select case(model)

! No extinction
			case(0)
				extinction_curve = 0.0

! Allen (1976) Milky Way
			case(1)
				extinction_curve = 0.0
				
! Seaton (1979) Milky Way
			case(2)
				allocate(invlambda(n))
				invlambda0 = 4.595
				gamm = 1.051
				C1 = -0.38
				C2 = 0.74
				C3 = 3.96
				C4 = 0.26
				invlambda = 1.d0 / lambda
				where(invlambda >= 5.9)
					extinction_curve = C1 + C2*invlambda + C3 / ((invlambda - invlambda0**2/invlambda)**2 + gamm**2) + &
						C4*(0.539*(invlambda-5.9)**2 + 0.0564*(invlambda-5.9)**3)
				elsewhere
					extinction_curve = C1 + C2*invlambda + C3 / ((invlambda - invlambda0**2/invlambda)**2 + gamm**2)
				endwhere

				extinction_curve = extinction_curve * Av				

! Fitzpatrick (1986) Large Magellanic Cloud
			case(3)
				allocate(invlambda(n))
				invlambda0 = 4.608
				gamm = 0.994
				C1 = -0.69
				C2 = 0.89
				C3 = 2.55
				C4 = 0.50
				invlambda = 1.d0 / lambda
				where(invlambda >= 5.9)
					extinction_curve = C1 + C2*invlambda + C3 / ((invlambda - invlambda0**2/invlambda)**2 + gamm**2) + &
						C4*(0.539*(invlambda-5.9)**2 + 0.0564*(invlambda-5.9)**3)
				elsewhere
					extinction_curve = C1 + C2*invlambda + C3 / ((invlambda - invlambda0**2/invlambda)**2 + gamm**2)
				endwhere

				extinction_curve = extinction_curve * Av				

! Prevot et al. (1984) Small Magellanic Cloud
			case(4)
				extinction_curve = 1.0
				
! Calzetti et al. (2000) Starburst galaxies
			case(5)
				Rv = 4.05
				
! From 0.12 micron to 0.63 micron
				where(lambda >= 0.12 .and. lambda < 0.63)
					extinction_curve = 2.659*(-2.156 + 1.509/lambda - 0.198/lambda**2 + 0.011/lambda**3) + Rv
				endwhere

! From 0.63 micron to 2.20 micron
				where(lambda >= 0.63 .and. lambda < 2.20)
					extinction_curve = 2.659*(-1.857 + 1.040/lambda) + Rv
				endwhere

! Below 0.12 micron
				C1 = 2.659*(-2.156 + 1.509/0.12 - 0.198/0.12**2 + 0.011/0.12**3) + Rv
				C2 = 2.659*(-2.156 + 1.509/0.11 - 0.198/0.11**2 + 0.011/0.11**3) + Rv
				slope = (C1-C2) / (0.12-0.11)
				zero = C1 - slope * 0.12
				where(lambda < 0.12)
					extinction_curve = slope * lambda + zero
				endwhere

! Above 2.20 micron
				C1 = 2.659*(-1.857 + 1.040/2.19) + Rv
				C2 = 2.659*(-1.857 + 1.040/2.20) + Rv
				slope = (C1-C2) / (2.19-2.20)
				zero = C1 - slope * 2.19
				where(lambda > 2.20)
					extinction_curve = slope * lambda + zero
				endwhere

! Avoid negative values
				where(extinction_curve < 0)
					extinction_curve = 0.0
				endwhere

				extinction_curve = extinction_curve * Av  / Rv
					
! Chiar & Tielens (2006) extinction curve
			case(6)
				call lin_interpol(prior%lambda_pixie, prior%extinction_pixie, lambda, extinction_curve)
! A_K = A_V * 0.09, from Whittet (2003)
				extinction_curve = extinction_curve * 0.09d0 * Av
		end select
		
		extinction_curve = 10.0**(-0.4*extinction_curve)
		
	end function extinction_curve


!-------------------------------------------------------------------
! Activation function for the neural networks
!-------------------------------------------------------------------
	function sigma(x)
	real(kind=KIND_FLOAT) :: x, sigma
		sigma = tanh(x)
	end function sigma
	
!-------------------------------------------------------------------
! Evaluate neural network
!-------------------------------------------------------------------	
	subroutine eval_neuralnetwork(indnet, pars)
! 	type(neural_network_type) :: neural
	real(kind=KIND_FLOAT) :: pars(:)
	integer :: i, j, indnet
	
		neural%net(indnet)%hidden = 0.e0
		neural%net(indnet)%output = 0.e0

! Hidden neurons
		do i = 1, neural%net(indnet)%nhidden
			do j = 1, neural%net(indnet)%ninput
				neural%net(indnet)%hidden(i) = neural%net(indnet)%hidden(i) + &
					neural%net(indnet)%input_hidden(j,i) * pars(j)
			enddo			
			neural%net(indnet)%hidden(i) = &
				sigma(neural%net(indnet)%hidden(i) + neural%net(indnet)%bias(i))
		enddo
		
! Output neurons
		do i = 1, neural%net(indnet)%noutput
			do j = 1, neural%net(indnet)%nhidden
				neural%net(indnet)%output(i) = neural%net(indnet)%output(i) + &
					neural%net(indnet)%hidden_output(j,i) * neural%net(indnet)%hidden(j)
			enddo
		enddo
		
	end subroutine eval_neuralnetwork

!------------------------------------
! Initialize the neural networks
!------------------------------------
subroutine initialize_neural_networks()
! type(neural_network_type) :: neural
integer :: i, j, k, loop, Nnets, npca, ind, n
real(kind=KIND_FLOAT), allocatable :: evec(:,:), weights(:)
character(len=2) :: str
	
	write(*,*) 'Initializing neural networks...'
	
! Read number of neural networks and number of input/outputs
	open(unit=12,file='NETWORKS/neural_topologies.dat',action='read',status='old')
	read(12,*) Nnets	
	allocate(neural%net(Nnets))
	neural%Nnets = Nnets
	do i = 1, Nnets
		read(12,*) ind, neural%net(i)%ninput, neural%net(i)%noutput
		allocate(neural%net(i)%input_norm(neural%net(i)%ninput))
		allocate(neural%net(i)%input_mean(neural%net(i)%ninput))
		allocate(neural%net(i)%output_norm(neural%net(i)%noutput))
		allocate(neural%net(i)%output_mean(neural%net(i)%noutput))
	enddo
	close(12)
		
! Reading PCA eigenvectors
	open(unit=12,file='NETWORKS/PCA_eigenvectors.bin',status='old',form='unformatted')
	read(12) npca
	write(*,*) 'Reading ', npca, ' PCA eigenvectors...'
	allocate(evec(npca,npca))
	read(12) evec
	close(12)
	do i = 1, Nnets
		neural%net(i)%PCAsize = npca
		allocate(neural%net(i)%PCAvector(npca))
		neural%net(i)%PCAvector = evec(i,:)
	enddo
	deallocate(evec)
	
! Read normalization factors for all the networks
	write(*,*) 'Reading neural network normalizations...'
	do i = 1, Nnets
		write(str,FMT='(I2)') i-1		
		open(unit=12,file='NETWORKS/normalization'//trim(adjustl(str))//'.net',&
			action='read',status='old')
		do j = 1, neural%net(i)%ninput
			read(12,*) neural%net(i)%input_mean(j), neural%net(i)%input_norm(j)
		enddo
		do j = 1, neural%net(i)%noutput
			read(12,*) neural%net(i)%output_mean(j), neural%net(i)%output_norm(j)
		enddo
		close(12)
	enddo
	
! Read mean spectrum
	write(*,*) 'Reading mean SED...'
	open(unit=12,file='NETWORKS/meanSED.bin',status='old',form='unformatted')
	allocate(neural%lambda(npca))
	allocate(neural%meanSpec(npca))
	allocate(neural%SED_highres(npca))
	
	read(12) neural%lambda
	read(12) neural%meanSpec
	close(12)
	
! Read neural networks' weights
	do i = 1, Nnets
		write(str,FMT='(I2)') i-1
		open(unit=12,file='NETWORKS/PCA'//trim(adjustl(str))//'.net',&
			action='read',status='old')
			
		read(12,*) neural%net(i)%ninput, neural%net(i)%nhidden, neural%net(i)%noutput, &
			ind
			
! Total number of weights+bias
		neural%net(i)%nunknowns = neural%net(i)%ninput * neural%net(i)%nhidden + &
			neural%net(i)%noutput * neural%net(i)%nhidden + neural%net(i)%nhidden
			
! Allocate memory
		allocate(neural%net(i)%input_hidden(neural%net(i)%ninput,neural%net(i)%nhidden))
		allocate(neural%net(i)%hidden_output(neural%net(i)%nhidden,neural%net(i)%noutput))
		allocate(neural%net(i)%bias(neural%net(i)%nhidden))
		allocate(neural%net(i)%input(neural%net(i)%ninput))
		allocate(neural%net(i)%hidden(neural%net(i)%nhidden))
		allocate(neural%net(i)%output(neural%net(i)%noutput))
		
! Read weights
		allocate(weights(neural%net(i)%nunknowns))
		do j = 1, neural%net(i)%nunknowns
			read(12,*) weights(j)
		enddo
		
! Assing weights and bias accordingly
		loop = 1
      do k = 1, neural%net(i)%nhidden
      	do j = 1, neural%net(i)%ninput
         	neural%net(i)%input_hidden(j,k) = weights(loop)
            loop = loop + 1
         enddo
      enddo

		do k = 1, neural%net(i)%nhidden
      	do j = 1, neural%net(i)%noutput
         	neural%net(i)%hidden_output(k,j) = weights(loop)
            loop = loop + 1
         enddo
      enddo

      do k = 1, neural%net(i)%nhidden
      	neural%net(i)%bias(k) = weights(loop)
         loop = loop + 1
      enddo
		close(12)
		
		deallocate(weights)
	enddo
	
	write(*,*) 'Read weights for ', Nnets, ' neural networks...'
	
end subroutine initialize_neural_networks

!-------------------------------------------------------------------
! Linear interpolation in database
! Note that the parameters in the code are in a different order to
! those on the database
! DB: i,tauv,q,N0,sigma,Y
! code: sigma,Y,N0,q,tauv,i
! and this modification has to be done before calling this routine
! In this routine, the parameters are passed as:
! Y,sigma,N0,q,tauv,i
!-------------------------------------------------------------------	
	subroutine lininterpol_db(pars, coefs)
	real(kind=KIND_FLOAT) :: pars(:), coefs(:), delta
	integer :: i, j, ndim, near(6), indices(6), ind
		
		ndim = 6

! Find the indices of the hypercube around the desired value
		call hunt(database%Y,database%nY,pars(1),near(1))
		call hunt(database%sig,database%nsig,pars(2),near(2))
		call hunt(database%n0,database%nn0,pars(3),near(3))
		call hunt(database%q,database%nq,pars(4),near(4))
		call hunt(database%tauv,database%ntauv,pars(5),near(5))
		call hunt(database%i,database%ni,pars(6),near(6))
		
		
! Extract the values of the function that will be used
		do i = 1, 2**ndim
			indices = near + database%ee(:,i)
			database%wrk(:,i) = database%coefs(:,indices(6),indices(5),indices(4),indices(3),indices(2),indices(1))
		enddo
		
! Do the actual linear interpolation
		do i = 1, ndim
			ind = database%indi(i)
			
			delta = -(pars(ind) - database%params(ind,near(ind))) / &
				(database%params(ind,near(ind)) - database%params(ind,near(ind)+1))
				
			do j = 1, 2**(ndim-i)
				
				database%wrk(:,j) = (database%wrk(:,2*j) - database%wrk(:,2*j-1)) * delta + database%wrk(:,2*j-1)
				
			enddo
			
		enddo
		
		coefs = database%wrk(:,1)		
		
	end subroutine lininterpol_db

	
!------------------------------------------------------------
! Check if the NetCDF action ended correctly
!------------------------------------------------------------
! 	subroutine check(status)
!    integer, intent ( in) :: status
! 
! 		if (status /= nf90_noerr) then
! 			print *, trim(nf90_strerror(status))
! 			stop
! 		endif
! 	end subroutine check
	
! !------------------------------------
! ! Read the full database
! !------------------------------------
! subroutine initialize_database_netcdf()
! character(len=320) :: name_var
! integer :: ndim, i, j, k, nr, c1, c2, nmax
! real(kind=KIND_FLOAT) :: dbpars(6), pars(6), coefs(20), output(119)
! 	
! 	write(*,*) 'Initializing database...'
! 
! ! Open the file in read-only access
! 	call check( nf90_open('DATABASE/compressed_database.nc', NF90_NOWRITE, database%db_id) )
! 	
! ! Get the varid of the data variable, based on their names
! 	call check( nf90_inq_dimid(database%db_id, "ny", database%ny_id) )
! 	call check( nf90_inquire_dimension(database%db_id, database%ny_id, name_var, database%ny) )
! 	
! 	call check( nf90_inq_dimid(database%db_id, "nN0", database%nn0_id) )
! 	call check( nf90_inquire_dimension(database%db_id, database%nn0_id, name_var, database%nn0) )
! 	
! 	call check( nf90_inq_dimid(database%db_id, "nq", database%nq_id) )
! 	call check( nf90_inquire_dimension(database%db_id, database%nq_id, name_var, database%nq) )
! 	
! 	call check( nf90_inq_dimid(database%db_id, "ntauv", database%ntauv_id) )
! 	call check( nf90_inquire_dimension(database%db_id, database%ntauv_id, name_var, database%ntauv) )
! 	
! 	call check( nf90_inq_dimid(database%db_id, "nsig", database%nsig_id) )
! 	call check( nf90_inquire_dimension(database%db_id, database%nsig_id, name_var, database%nsig) )
! 	
! 	call check( nf90_inq_dimid(database%db_id, "ni", database%ni_id) )
! 	call check( nf90_inquire_dimension(database%db_id, database%ni_id, name_var, database%ni) )
! 	
! 	call check( nf90_inq_dimid(database%db_id, "npca", database%npca_id) )
! 	call check( nf90_inquire_dimension(database%db_id, database%npca_id, name_var, database%npca) )
! 	
! 	call check( nf90_inq_dimid(database%db_id, "nlam", database%nlam_id) )
! 	call check( nf90_inquire_dimension(database%db_id, database%nlam_id, name_var, database%nlam) )
! 	
! 	call check( nf90_inq_varid(database%db_id, "lambda", database%lam_id) )
! 	call check( nf90_inq_varid(database%db_id, "base", database%base_id) )
! 	call check( nf90_inq_varid(database%db_id, "coefs", database%coefs_id) )
! 	call check( nf90_inq_varid(database%db_id, "meansed", database%meansed_id) )
! 	call check( nf90_inq_varid(database%db_id, "y", database%y_id) )
! 	call check( nf90_inq_varid(database%db_id, "sig", database%sig_id) )
! 	call check( nf90_inq_varid(database%db_id, "n0", database%n0_id) )
! 	call check( nf90_inq_varid(database%db_id, "q", database%q_id) )
! 	call check( nf90_inq_varid(database%db_id, "tauv", database%tauv_id) )
! 	call check( nf90_inq_varid(database%db_id, "i", database%i_id) )
! 		
! ! Read the database
! 	allocate(database%lambda(database%nlam))
! 	allocate(database%base(database%npca,database%nlam))
! 	allocate(database%coefs(database%npca,database%ni,database%ntauv,database%nq,database%nn0,database%nsig,database%ny))
! 	allocate(database%meansed(database%nlam))
! 	allocate(database%y(database%ny))
! 	allocate(database%sig(database%nsig))
! 	allocate(database%n0(database%nn0))
! 	allocate(database%q(database%nq))
! 	allocate(database%tauv(database%ntauv))
! 	allocate(database%i(database%ni))
! 	allocate(database%wrk(database%npca,2**6))
! 	allocate(database%SED_highres(database%nlam))
! 	
! 	call check( nf90_get_var(database%db_id, database%lam_id, database%lambda) )
! 	call check( nf90_get_var(database%db_id, database%base_id, database%base) )
! 	call check( nf90_get_var(database%db_id, database%coefs_id, database%coefs) )
! 	call check( nf90_get_var(database%db_id, database%meansed_id, database%meansed) )
! 	call check( nf90_get_var(database%db_id, database%y_id, database%y) )
! 	call check( nf90_get_var(database%db_id, database%sig_id, database%sig) )
! 	call check( nf90_get_var(database%db_id, database%n0_id, database%n0) )
! 	call check( nf90_get_var(database%db_id, database%q_id, database%q) )
! 	call check( nf90_get_var(database%db_id, database%tauv_id, database%tauv) )
! 	call check( nf90_get_var(database%db_id, database%i_id, database%i) )
! 	
! 	call check( nf90_close(database%db_id) )
! 	
! 	print *, 'Done'
! 	
! ! Do some precomputations for accelerating the linear interpolation routines on the database
! 	ndim = 6
! 	allocate(database%indi(ndim))
! 	do i = 1, ndim
! 		database%indi(i) = ndim - i + 1
! 	enddo
! 	
! 	allocate(database%ee(ndim,2**ndim))
! 	do i = 1, ndim
! 		nr = 2**(i-1)
! 		c1 = 1
! 		do j = 1, nr
! 			c2 = 1
! 			do k = 1, 2**ndim/nr
! 				database%ee(ndim-database%indi(i)+1,c1) = (c2-1) / 2**(ndim-i)
! 				c1 = c1 + 1
! 				c2 = c2 + 1
! 			enddo
! 		enddo
! 	enddo
! 	
! ! Make array of parameters
! 	nmax = maxval( (/database%nY,database%nsig,database%nn0,database%nq,database%ntauv,database%ni/) )
! 	allocate(database%params(6,nmax))
! 	
! 	database%params = 1.d10
! 	
! 	database%params(1,1:database%nY) = database%Y
! 	database%params(2,1:database%nsig) = database%sig
! 	database%params(3,1:database%nn0) = database%n0
! 	database%params(4,1:database%nq) = database%q
! 	database%params(5,1:database%ntauv) = database%tauv
! 	database%params(6,1:database%ni) = database%i
! 		
! end subroutine initialize_database_netcdf

!------------------------------------
! Read the full database
!------------------------------------
subroutine initialize_database()
character(len=320) :: name_var
integer :: ndim, i, j, k, nr, c1, c2, nmax
real(kind=KIND_FLOAT) :: dbpars(6), pars(6), coefs(20), output(119)

	write(*,*) 'Initializing database...'

! Open the file in read-only access
	open(unit=12,file='DATABASE/compressed_database.bin', action='read', status='old', form='unformatted')

	read(12) database%nY, database%nn0, database%nq, database%ntauv, database%nsig, database%ni, database%npca, database%nlam

! Read the database
 	allocate(database%lambda(database%nlam))
 	allocate(database%base(database%npca,database%nlam))
 	allocate(database%coefs(database%npca,database%ni,database%ntauv,database%nq,database%nn0,database%nsig,database%ny))
 	allocate(database%meansed(database%nlam))	
	allocate(database%y(database%ny))
	allocate(database%sig(database%nsig))
	allocate(database%n0(database%nn0))
	allocate(database%q(database%nq))
	allocate(database%tauv(database%ntauv))
	allocate(database%i(database%ni))
 	allocate(database%wrk(database%npca,2**6))
 	allocate(database%SED_highres(database%nlam))
	
	read(12) database%Y
	read(12) database%sig
	read(12) database%n0	
	read(12) database%q	
	read(12) database%tauv	
	read(12) database%i	

 	read(12) database%lambda
 	read(12) database%base
 	read(12) database%coefs
 	read(12) database%meanSED

	print *, 'Done'

! Do some precomputations for accelerating the linear interpolation routines on the database
	ndim = 6
	allocate(database%indi(ndim))
	do i = 1, ndim
		database%indi(i) = ndim - i + 1
	enddo

	allocate(database%ee(ndim,2**ndim))
	do i = 1, ndim
		nr = 2**(i-1)
		c1 = 1
		do j = 1, nr
			c2 = 1
			do k = 1, 2**ndim/nr
				database%ee(ndim-database%indi(i)+1,c1) = (c2-1) / 2**(ndim-i)
				c1 = c1 + 1
				c2 = c2 + 1
			enddo
		enddo
	enddo

! Make array of parameters
	nmax = maxval( (/database%nY,database%nsig,database%nn0,database%nq,database%ntauv,database%ni/) )
	allocate(database%params(6,nmax))

	database%params = 1.d10

	database%params(1,1:database%nY) = database%Y
	database%params(2,1:database%nsig) = database%sig
	database%params(3,1:database%nn0) = database%n0
	database%params(4,1:database%nq) = database%q
	database%params(5,1:database%ntauv) = database%tauv
	database%params(6,1:database%ni) = database%i

end subroutine initialize_database


!------------------------------------
! Evaluate neural network
!------------------------------------
subroutine neural_eval(pars, output, include_agn, reddening_law)
! type(neural_network_type) :: neural
integer :: loopnet, j, include_agn, reddening_law
real(kind=KIND_FLOAT) :: pars(:), output(:), PCA_coeff
real(kind=KIND_FLOAT), allocatable :: pars_normalized(:), agn(:)

	output = 0.0
	allocate(pars_normalized(neural%net(1)%ninput))

	allocate(agn(size(neural%lambda)))
	
! Normalize input
	do loopnet = 1, neural%Nnets
		do j = 1, neural%net(loopnet)%ninput
			pars_normalized(j) = (pars(j) - neural%net(loopnet)%input_mean(j)) / &
				neural%net(loopnet)%input_norm(j)
		enddo
		
! Evaluate neural network
		call eval_neuralnetwork(loopnet, pars_normalized)
		
! Apply inverse normalization
		PCA_coeff = neural%net(loopnet)%output(1) * neural%net(loopnet)%output_norm(1) + &
			neural%net(loopnet)%output_mean(1)

! Add the PCA coefficient times the eigenvector
		output = output + PCA_coeff * neural%net(loopnet)%PCAvector
		
	enddo

! Add the mean spectrum and transform to flux
	output = 10.e0**(output + neural%meanSpec)

! Add the AGN contribution if desired
	if (include_agn == 1) then
		call agn_shape(neural%lambda, agn)
		output = output + agn
	endif
	
! Finally, transform to Jansky
	output = (neural%lambda*1.e-4 / 2.99792458d10) / 1.e-26 * output

! If extinction is taken into account
	output = output * extinction_curve(neural%lambda, reddening_law, pars(8))
				
	deallocate(pars_normalized)
	deallocate(agn)
	
end subroutine neural_eval

!------------------------------------
! Evaluate neural network
!------------------------------------
subroutine lininterpol_eval(pars, output, include_agn, reddening_law)
integer :: loopnet, j, include_agn, reddening_law
real(kind=KIND_FLOAT) :: pars(:), output(:)
real(kind=KIND_FLOAT), allocatable :: PCAcoefs(:), dbpars(:), agn(:)

	allocate(PCAcoefs(database%npca))
	allocate(dbpars(6))

	allocate(agn(size(neural%lambda)))

! Rotate Y and sigma, which are not in the same order as in the DB
	dbpars = pars(1:6)
	dbpars(1) = pars(2)
	dbpars(2) = pars(1)

	output = 0.0

	call lininterpol_db(dbpars, PCAcoefs)

	do j = 1, database%npca
		output = output + PCAcoefs(j) * database%base(j,:)
	enddo
	
! Add the mean spectrum and transform to flux
	output = 10.e0**(output + database%meanSED)	
	
! Add the AGN contribution if desired
	if (include_agn == 1) then
		call agn_shape(neural%lambda, agn)
		output = output + agn
	endif

	
! Finally, transform to Jansky
	output = (database%lambda*1.e-4 / 2.99792458d10) / 1.e-26 * output

! If extinction is taken into account
	output = output * extinction_curve(database%lambda, reddening_law, pars(8))

	deallocate(PCAcoefs)
	deallocate(dbpars)
	deallocate(agn)
	
end subroutine lininterpol_eval

!------------------------------------
! Evaluate neural network
!------------------------------------
subroutine neural_eval_filter(pars, obs, include_agn, reddening_law)
type(observations_type) :: obs
real(kind=KIND_FLOAT) :: pars(:)
integer :: i, include_agn, reddening_law

	call neural_eval(pars, neural%SED_highres, include_agn, reddening_law)

! Points obtained with filters
	do i = 1, obs%npoints

! Reinterpolate the SED to the wavelength axis of the filter
! This gives a better sampling since the SED is usually smooth and the filter extent
! in wavelength is much reduced than the SED

! Take into account the redshift in this case
		call lin_interpol(neural%lambda*(one+pars(9)), neural%SED_highres, &
			obs%filter(i)%lambda, obs%filter(i)%SED_interpol)

		obs%model_eval(i) = &
			int_tabulated(obs%filter(i)%lambda, obs%filter(i)%SED_interpol*obs%filter(i)%transmission) /&
			int_tabulated(obs%filter(i)%lambda, obs%filter(i)%transmission)

	enddo

	if (version_bc < 0) then
! Points obtained with spectrograph

! Reinterpolate the SED to the wavelength axis of the observation
! This gives a better sampling since the SED is usually smooth and the filter extent
! in wavelength is much reduced than the SED
		call lin_interpol(neural%lambda*(one+pars(9)), neural%SED_highres, obs%obs_x(obs%npoints+1:obs%npoints_total),&
			obs%model_eval(obs%npoints+1:obs%npoints_total))
	endif
					 
	obs%model_eval = obs%model_eval / 1.e10 * 10.d0**pars(7)

end subroutine neural_eval_filter

!------------------------------------
! Evaluate neural network
!------------------------------------
subroutine lininterpol_eval_filter(pars, obs, include_agn, reddening_law)
type(observations_type) :: obs
real(kind=KIND_FLOAT) :: pars(:)
integer :: i, include_agn, reddening_law

	call lininterpol_eval(pars, database%SED_highres, include_agn, reddening_law)

! Points obtained with filters
	do i = 1, obs%npoints

! Reinterpolate the SED to the wavelength axis of the filter
! This gives a better sampling since the SED is usually smooth and the filter extent
! in wavelength is much reduced than the SED

! Take into account the redshift in this case
		call lin_interpol(database%lambda*(one+pars(9)), database%SED_highres, &
			obs%filter(i)%lambda, obs%filter(i)%SED_interpol)

		obs%model_eval(i) = &
			int_tabulated(obs%filter(i)%lambda, obs%filter(i)%SED_interpol*obs%filter(i)%transmission) /&
			int_tabulated(obs%filter(i)%lambda, obs%filter(i)%transmission)

	enddo

	if (version_bc < 0) then
! Points obtained with spectrograph

! Reinterpolate the SED to the wavelength axis of the observation
! This gives a better sampling since the SED is usually smooth and the filter extent
! in wavelength is much reduced than the SED
		call lin_interpol(database%lambda*(one+pars(9)), database%SED_highres, obs%obs_x(obs%npoints+1:obs%npoints_total),&
			obs%model_eval(obs%npoints+1:obs%npoints_total))
	endif
					 
	obs%model_eval = obs%model_eval / 1.e10 * 10.d0**pars(7)

end subroutine lininterpol_eval_filter


end module ann_mod
