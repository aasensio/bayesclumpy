module maths_mod
use types_mod, only : KIND_FLOAT, very_large
implicit none
contains

!-------------------------------------------------------------------
! Detect version
!-------------------------------------------------------------------
	function detect_version()
	real(kind=KIND_FLOAT) :: detect_version
		open(unit=15,file='BayesClumpy.version',action='read',status='old', err=101)
		read(15,*) detect_version
		close(15)
		if (detect_version < 0) then
			write(*,FMT='(A,F3.1,A1)') ' BayesClumpy v', abs(detect_version), 's'
		else
			write(*,FMT='(A,F3.1)') ' BayesClumpy v', abs(detect_version)
		endif
		return
		
101   detect_version = 2.0
		write(*,FMT='(A,F3.1)') ' BayesClumpy v', abs(detect_version)		
		
	end function detect_version

!-------------------------------------------------------------------
! Read n lines at the selected unit
!-------------------------------------------------------------------
	subroutine lb(u, n)
	integer :: u, n
   integer :: i
   	do i = 1, n
      	read(u,*)
      enddo
	end subroutine lb
	
!-------------------------------------------------------------------
! Returns which is the index of the array wave_total closer to wave
!-------------------------------------------------------------------
	function close_to(wave_total, wave)
	real(kind=KIND_FLOAT) :: wave_total(:), wave
	real(kind=KIND_FLOAT), allocatable :: diff(:)
	integer :: n, i, location(1)
	integer :: close_to

		n = size(wave_total)
	   allocate(diff(n))
	   diff = abs(wave_total-wave)
	   location = minloc(diff)
	   deallocate(diff)
	   close_to = location(1)

	end function close_to
!-------------------------------------------------------------
! Carry out the Cholesky decomposition of a symmetric matrix
!-------------------------------------------------------------
	subroutine cholesky(a,n,p)
	integer :: n
	real(kind=KIND_FLOAT) :: a(n,n), p(n)
	integer :: i, j, k
	real(kind=KIND_FLOAT) :: sum
	
		do i = 1, n
			do j = i, n
				sum = a(i,j)
				do k = i-1, 1, -1
					sum = sum - a(i,k)*a(j,k)
				enddo
				if (i == j) then
					if (sum == 0.d0) then
						print *, 'Cholesky decomposition failed...'	
					endif
					p(i) = sqrt(sum)
				else
					a(j,i) = sum / p(i)
				endif
			enddo
		enddo
						
	end subroutine cholesky
	
!-------------------------------------------------------------
! Initialize the random number generator
!-------------------------------------------------------------
	subroutine init_random_seed()
	integer :: i, n, clock
   integer, dimension(:), allocatable :: seed
          
   	call random_seed(size = n)
      allocate(seed(n))
          
      call system_clock(count=clock)
          
      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      call random_seed(put = seed)
          
      deallocate(seed)
   end subroutine init_random_seed
          
!-------------------------------------------------------------
! Generates a random number following an uniform distribution in the interval [0,1]
! Call it with idum<0 for initialization
!-------------------------------------------------------------
	function randomu()
	real(kind=KIND_FLOAT) :: randomu
	
		call random_number(randomu)
							
	end function randomu
		
!-------------------------------------------------------------
! Generates a random number following an normal distribution with zero mean
! and unit variance
!-------------------------------------------------------------
	function randomn()

	real(kind=KIND_FLOAT) :: randomn

	real(kind=KIND_FLOAT) :: u, sum
	real(kind=KIND_FLOAT), save :: v, sln
	logical, save :: second = .false.
	real(kind=KIND_FLOAT), parameter :: one = 1.0, vsmall = tiny( one )

! if second, use the second random number generated on last call
	if (second) then

		second = .false.
  		randomn = v*sln
	else
! first call; generate a pair of random normals

  		second = .true.
  		do
    		call random_number( u )
    		call random_number( v )
    		u = scale( u, 1 ) - one
    		v = scale( v, 1 ) - one
    		sum = u*u + v*v + vsmall         ! vsmall added to prevent log(zero) / zero
    		if(sum < one) exit
  		end do
  		sln = sqrt(- scale( log(sum), 1 ) / sum)
  		randomn = u*sln
	end if

	return
	end function randomn

! 	function randomn()
! 	real(kind=KIND_FLOAT) :: randomn, gasdev, fac, gset, rsq, v1, v2, ran1
! 	integer :: iset
! 	save :: iset, gset		
! 		if (iset == 0) then
! 			rsq = 0.0
! 			do while(rsq > 1 .or. rsq == 0.0)
! 				v1 = 2.0*randomu()-1.0
! 				v2 = 2.0*randomu()-1.0
! 				rsq = v1**2+v2**2
! 			enddo
! 			fac = sqrt(-2.0*log(rsq)/rsq)
! 			gset = v1*fac
! 			randomn = v2*fac
! 			iset = 1
! 		else
! 			randomn = gset
! 			iset = 0
! 		endif
! 	end function randomn
	
!-------------------------------------------------------------
! Generates a multivariate normal random number with a given
! mean and covariance matrix
!-------------------------------------------------------------
	function mrandomn(rmean,covar)
	integer :: idum, n, i, j
	real(kind=KIND_FLOAT) :: rmean(:), covar(:,:), mrandomn(size(rmean))
	real(kind=KIND_FLOAT) :: chol(size(rmean),size(rmean)), p(size(rmean)), eps(size(rmean))
	
		n = size(rmean)
				
		chol = covar
		
		do i = 1, n
			eps(i) = randomn()
			chol(i,i) = chol(i,i) + 1.d-7    ! Some regularization
		enddo
								
		call cholesky(chol,n,p)
										
		do j = 1, n			
			do i = j, n
				chol(j,i) = 0.d0
			enddo
			chol(j,j) = p(j)
		enddo
								
		mrandomn = matmul(chol,eps) + rmean			
				
		
	end function mrandomn
	
!-------------------------------------------------------------
! Return the diagonal of a matrix
!-------------------------------------------------------------
	function diagon(mat)
	real(kind=KIND_FLOAT) :: mat(:,:)
	real(kind=KIND_FLOAT) :: diagon(size(mat(1,:)))
	integer :: i
	
		do i = 1, size(mat(1,:))
			diagon(i) = mat(i,i)
		enddo
		
	end function diagon

!-------------------------------------------------------------
! CUBINT approximates an integral using cubic interpolation of data.
!  Parameters:
!
!    Input, real FTAB(NTAB), contains the tabulated function
!    values, FTAB(I) = F(XTAB(I)).
!
!    Input, real XTAB(NTAB), contains the points at which the
!    function was tabulated.  XTAB should contain distinct
!    values, given in ascending order.
!
!    Input, integer NTAB, the number of tabulated points.
!    NTAB must be at least 4.
!
!    Input, integer IA, the entry of XTAB at which integration
!    is to begin.  IA must be no less than 1 and no greater
!    than NTAB.
!
!    Input, integer IB, the entry of XTAB at which integration
!    is to end.  IB must be no less than 1 and no greater than
!    NTAB.
!
!    Output, real RESULT, the approximate value of the
!    integral from XTAB(IA) to XTAB(IB) of the function.
!
!    Output, real ERROR, an estimate of the error in
!    integration.
!-------------------------------------------------------------
	subroutine cubint ( ftab, xtab, ntab, ia, ib, result, error )

  	integer ntab
!
  	real(kind=KIND_FLOAT) :: c, d1, d2, d3, error, ftab(ntab), h1, h2, h3, h4
  	integer :: i, ia, ib, ind, it, j, k
  	real(kind=KIND_FLOAT) r1, r2, r3, r4, result, s, term, xtab(ntab)
!
  	result = 0.0E+00
  	error = 0.0E+00
 
  	if ( ia == ib ) then
    	return
  	end if
 
  	if ( ntab < 4 ) then
    	write ( *, '(a)' ) ' '
    	write ( *, '(a)' ) 'CUBINT - Fatal error!'
    	write ( *, '(a,i6)' ) '  NTAB must be at least 4, but input NTAB = ',ntab
    	stop
  	endif
 
  	if ( ia < 1 ) then
    	write ( *, '(a)' ) ' '
    	write ( *, '(a)' ) 'CUBINT - Fatal error!'
    	write ( *, '(a,i6)' ) '  IA must be at least 1, but input IA = ',ia
    	stop
  	endif
 
  	if ( ia > ntab ) then
   	write ( *, '(a)' ) ' '
    	write ( *, '(a)' ) 'CUBINT - Fatal error!'
    	write ( *, '(a,i6)' ) '  IA must be <= NTAB, but input IA=',ia
    	stop
  	endif
 
  	if ( ib < 1 ) then
    	write ( *, '(a)' ) ' '
    	write ( *, '(a)' ) 'CUBINT - Fatal error!'
    	write ( *, '(a,i6)' ) '  IB must be at least 1, but input IB = ',ib
    	stop
  	endif
 
  	if ( ib > ntab ) then
    	write ( *, '(a)' ) ' '
    	write ( *, '(a)' ) 'CUBINT - Fatal error!'
    	write ( *, '(a,i6)' ) '  IB must be <= NTAB, but input IB=',ib
    	stop
  	endif
!
!  Temporarily switch IA and IB, and store minus sign in IND
!  so that, while integration is carried out from low X's
!  to high ones, the sense of the integral is preserved.
!
  	if ( ia > ib ) then
    	ind = -1
    	it = ib
    	ib = ia
    	ia = it
  	else
    	ind = 1
  	endif
 
  	s = 0.0E+00
  	c = 0.0E+00
  	r4 = 0.0E+00
  	j = ntab-2
  	if ( ia < ntab-1 .or. ntab == 4 ) then
    	j=max(3,ia)
  	endif

  	k = 4
  	if ( ib > 2 .or. ntab == 4 ) then
    	k=min(ntab,ib+2)-1
  	endif
 
  	do i = j, k
 
    	if ( i <= j ) then
 
      	h2 = xtab(j-1)-xtab(j-2)
      	d3 = (ftab(j-1)-ftab(j-2)) / h2
      	h3 = xtab(j)-xtab(j-1)
      	d1 = (ftab(j)-ftab(j-1)) / h3
      	h1 = h2+h3
      	d2 = (d1-d3)/h1
      	h4 = xtab(j+1)-xtab(j)
      	r1 = (ftab(j+1)-ftab(j)) / h4
      	r2 = (r1-d1) / (h4+h3)
      	h1 = h1+h4
      	r3 = (r2-d2) / h1
 
      	if ( ia <= 1 ) then
        		result = h2 * (ftab(1)+h2*(0.5*d3-h2*(d2/6.0-(h2+h3+h3)*r3/12.)))
        		s = -h2**3 * (h2*(3.0*h2+5.0*h4)+10.0*h3*h1)/60.0
      	endif
 
    	else
 
	      h4 = xtab(i+1)-xtab(i)
      	r1 = (ftab(i+1)-ftab(i))/h4
      	r4 = h4+h3
      	r2 = (r1-d1)/r4
      	r4 = r4+h2
      	r3 = (r2-d2)/r4
      	r4 = (r3-d3)/(r4+h1)
 
    	endif
 
    	if ( i > ia .and. i <= ib ) then
 
      	term = h3*((ftab(i)+ftab(i-1))*0.5-h3*h3*(d2+r2+(h2-h4)*r3) / 12.0 )
      	result = result+term
      	c = h3**3*(2.0E+00 *h3*h3+5.*(h3*(h4+h2) + 2.0 * h2 * h4 ) ) / 120.0E+00
      	error = error+(c+s)*r4
 
      	if ( i /= j ) then
        		s = c
      	else
        		s = s+c+c
      	endif
 
    	else
 
	      error = error+r4*s
 
    	endif
 
    	if ( i >= k ) then
 
      	if ( ib >= ntab ) then
	        	term = h4*(ftab(ntab) - h4*(0.5*r1+h4*(r2/6.0 +(h3+h3+h4)*r3/12.)))
        		result = result + term
        		error = error - h4**3 * r4 * &
          		( h4 * ( 3.0 * h4 + 5.0 * h2 ) &
          		+ 10.0 * h3 * ( h2 + h3 + h4 ) ) / 60.0E+00
      	endif
 
      	if ( ib >= ntab-1 ) error=error+s*r4
    	else
	      h1 = h2
      	h2 = h3
      	h3 = h4
      	d1 = r1
      	d2 = r2
      	d3 = r3
    	endif
 
  	enddo
!
!  Restore original values of IA and IB, reverse signs
!  of RESULT and ERROR, to account for integration
!  that proceeded from high X to low X.
!
  	if ( ind /= 1 ) then
    	it = ib
    	ib = ia
    	ia = it
    	result = -result
    	error = -error
  	endif
 
  	return
	end subroutine

! ---------------------------------------------------------
!	Given an array xx(1:n), and given a value x, returns a value jlo such that x is between
!	xx(jlo) and xx(jlo+1). xx(1:n) must be monotonic, either increasing or decreasing.
!	jlo=0 or jlo=n is returned to indicate that x is out of range. jlo on input is taken as
!	the initial guess for jlo on output.
! ---------------------------------------------------------
	subroutine hunt(xx,n,x,jlo)
	integer :: jlo,n
	real(kind=KIND_FLOAT) :: x,xx(n)
	integer :: inc,jhi,jm
	logical :: ascnd

		ascnd=xx(n).ge.xx(1)
		if (jlo.le.0.or.jlo.gt.n) then
			jlo=0
			jhi=n+1
			goto 3
		endif
		inc=1
		if (x.ge.xx(jlo).eqv.ascnd) then
1     	jhi=jlo+inc
			if (jhi.gt.n) then
				jhi=n+1
			else if (x.ge.xx(jhi).eqv.ascnd) then
				jlo=jhi
				inc=inc+inc
				goto 1
			endif
		else
			jhi=jlo
2     	jlo=jhi-inc
			if (jlo.lt.1) then
				jlo=0
			else if (x.lt.xx(jlo).eqv.ascnd) then
				jhi=jlo
				inc=inc+inc
				goto 2
			endif
		endif
3 		if (jhi-jlo.eq.1) then
			if(x.eq.xx(n)) jlo=n-1
			if(x.eq.xx(1)) jlo=1
			return
		endif
		jm = (jhi+jlo)/2
		if (x.ge.xx(jm).eqv.ascnd) then
			jlo=jm
		else
			jhi=jm
		endif
		goto 3
	end subroutine hunt

! ---------------------------------------------------------
! Given xa(:) and ya(:) which tabulate a function, returns the interpolation using
! linear interpolation of vector x(:) in y(:)
! ---------------------------------------------------------
	subroutine lin_interpol(xa,ya,x,y)
   real(kind=KIND_FLOAT), INTENT(IN) :: xa(:), ya(:), x(:)
   real(kind=KIND_FLOAT), INTENT(INOUT) :: y(:)
   integer :: i, n, na
   integer :: loc, jlo

   	n = size(x)
   	na = size(xa)

      do i = 1, n
      	call hunt(xa,na,x(i),jlo)
         if (jlo == 0) then
         	y(i) = ya(1)
         else if (jlo == na) then
         	y(i) = ya(na)
         else
            y(i) = (ya(jlo+1)-ya(jlo))/(xa(jlo+1)-xa(jlo)) * (x(i)-xa(jlo)) + ya(jlo)
         endif
      enddo

   end subroutine lin_interpol


!----------------------------------------------------------------
! This function integrates a tabulated function
!----------------------------------------------------------------		
	function int_tabulated(x, f)
	real(kind=KIND_FLOAT) :: x(:), f(:), int_tabulated, res, error_res
	integer :: n
		n = size(x)
		call cubint (f, x, n, 1, n, res, error_res)
		int_tabulated = res
	end function int_tabulated

	subroutine qsortd(x,ind,n)
	
	implicit none
	integer, parameter  :: dp = 4
	
	real (kind=KIND_FLOAT), intent(in)  :: x(:)
	integer, intent(out)   :: ind(:)
	integer, intent(in)    :: n
	
	!***************************************************************************
	
	!                                                         robert renka
	!                                                 oak ridge natl. lab.
	
	!   this subroutine uses an order n*log(n) quick sort to sort a real (dp)
	! array x into increasing order.  the algorithm is as follows.  ind is
	! initialized to the ordered sequence of indices 1,...,n, and all interchanges
	! are applied to ind.  x is divided into two portions by picking a central
	! element t.  the first and last elements are compared with t, and
	! interchanges are applied as necessary so that the three values are in
	! ascending order.  interchanges are then applied so that all elements
	! greater than t are in the upper portion of the array and all elements
	! less than t are in the lower portion.  the upper and lower indices of one
	! of the portions are saved in local arrays, and the process is repeated
	! iteratively on the other portion.  when a portion is completely sorted,
	! the process begins again by retrieving the indices bounding another
	! unsorted portion.
	
	! input parameters -   n - length of the array x.
	
	!                      x - vector of length n to be sorted.
	
	!                    ind - vector of length >= n.
	
	! n and x are not altered by this routine.
	
	! output parameter - ind - sequence of indices 1,...,n permuted in the same
	!                          fashion as x would be.  thus, the ordering on
	!                          x is defined by y(i) = x(ind(i)).
	
	!*********************************************************************
	
	! note -- iu and il must be dimensioned >= log(n) where log has base 2.
	
	!*********************************************************************
	
	integer   :: iu(21), il(21)
	integer   :: m, i, j, k, l, ij, it, itt, indx
	real(kind=KIND_FLOAT)      :: r
	real(kind=KIND_FLOAT) :: t
	
	! local parameters -
	
	! iu,il =  temporary storage for the upper and lower
	!            indices of portions of the array x
	! m =      index for iu and il
	! i,j =    lower and upper indices of a portion of x
	! k,l =    indices in the range i,...,j
	! ij =     randomly chosen index between i and j
	! it,itt = temporary storage for interchanges in ind
	! indx =   temporary index for x
	! r =      pseudo random number for generating ij
	! t =      central element of x
	
	if (n <= 0) return
	
	! initialize ind, m, i, j, and r
	
	do  i = 1, n
	ind(i) = i
	end do
	
	m = 1
	i = 1
	j = n
	r = .375
	
	! top of loop
	
	20 if (i >= j) go to 70
	if (r <= .5898437) then
	r = r + .0390625
	else
	r = r - .21875
	end if
	
	! initialize k
	
	30 k = i
	
	! select a central element of x and save it in t
	
	ij = i + r*(j-i)
	it = ind(ij)
	t = x(it)
	
	! if the first element of the array is greater than t,
	!   interchange it with t
	
	indx = ind(i)
	if (x(indx) > t) then
	ind(ij) = indx
	ind(i) = it
	it = indx
	t = x(it)
	end if
	
	! initialize l
	
	l = j
	
	! if the last element of the array is less than t,
	!   interchange it with t
	
	indx = ind(j)
	if (x(indx) >= t) go to 50
	ind(ij) = indx
	ind(j) = it
	it = indx
	t = x(it)
	
	! if the first element of the array is greater than t,
	!   interchange it with t
	
	indx = ind(i)
	if (x(indx) <= t) go to 50
	ind(ij) = indx
	ind(i) = it
	it = indx
	t = x(it)
	go to 50
	
	! interchange elements k and l
	
	40 itt = ind(l)
	ind(l) = ind(k)
	ind(k) = itt
	
	! find an element in the upper part of the array which is
	!   not larger than t
	
	50 l = l - 1
	indx = ind(l)
	if (x(indx) > t) go to 50
	
	! find an element in the lower part of the array whcih is not smaller than t
	
	60 k = k + 1
	indx = ind(k)
	if (x(indx) < t) go to 60
	
	! if k <= l, interchange elements k and l
	
	if (k <= l) go to 40
	
	! save the upper and lower subscripts of the portion of the
	!   array yet to be sorted
	
	if (l-i > j-k) then
	il(m) = i
	iu(m) = l
	i = k
	m = m + 1
	go to 80
	end if
	
	il(m) = k
	iu(m) = j
	j = l
	m = m + 1
	go to 80
	
	! begin again on another unsorted portion of the array
	
	70 m = m - 1
	if (m == 0) return
	i = il(m)
	j = iu(m)
	
	80 if (j-i >= 11) go to 30
	if (i == 1) go to 20
	i = i - 1
	
	! sort elements i+1,...,j.  note that 1 <= i < j and j-i < 11.
	
	90 i = i + 1
	if (i == j) go to 70
	indx = ind(i+1)
	t = x(indx)
	it = indx
	indx = ind(i)
	if (x(indx) <= t) go to 90
	k = i
	
	100 ind(k+1) = ind(k)
	k = k - 1
	indx = ind(k)
	if (t < x(indx)) go to 100
	
	ind(k+1) = it
	go to 90
	end subroutine qsortd

! ---------------------------------------------------------
! Given x(:) and y(:) which tabulate a function and the derivative at the boundary points
! this function returns the second derivative of the spline at each point
! ---------------------------------------------------------
		subroutine splin1(x,y,yp1,ypn,y2)
		real(kind=KIND_FLOAT), INTENT(IN) :: x(:), y(:), yp1, ypn
		real(kind=KIND_FLOAT), INTENT(INOUT) :: y2(size(x))
		integer :: n, i, k
		real(kind=KIND_FLOAT) :: p, qn, sig, un, u(size(x))

			n = size(x)
			
			if (yp1 > .99e30) then
				y2(1) = 0.e0
				u(1) = 0.e0
			else
				y2(1) = -0.5e0
				u(1) = (3.e0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
			endif

			do i = 2, n-1
				sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))				
				p = sig * y2(i-1)+2.e0
				y2(i) = (sig-1.e0)/p
				u(i) = (6.e0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/&
					(x(i+1)-x(i-1))-sig*u(i-1))/p
			enddo
			if (ypn > .99e30) then
				qn = 0.e0
				un = 0.e0
			else
				qn = 0.5e0
				un = (3.e0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
			endif
			
			y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.e0)

			do k = n-1, 1, -1
				y2(k) = y2(k)*y2(k+1)+u(k)
			enddo

		end subroutine splin1

! ---------------------------------------------------------
! Given xa(:) and ya(:) which tabulate a function, returns the interpolation using
! splines of vector x(:) in y(:)
! ---------------------------------------------------------
		subroutine spline(xa,ya,x,y)
		real(kind=KIND_FLOAT), INTENT(INOUT) :: y(:)
		real(kind=KIND_FLOAT), INTENT(IN) :: xa(:), ya(:), x(:)
		real(kind=KIND_FLOAT) :: y2a(size(xa))
		integer :: n_x, n, i, k, khi, klo
		real(kind=KIND_FLOAT) :: a, b, h, extrap
			
			n = size(xa)
			n_x = size(x)
			call splin1(xa,ya,very_large,very_large,y2a)

			do i = 1, n_x					

! Downward extrapolation 
				if (x(i) < xa(1)) then
!					y(i) = ya(1)
					y(i) = ya(1) + (ya(1)-ya(2))/(xa(1)-xa(2)) * (xa(1) - x(i))
				else 

! Upward extrapolation
				if (x(i) > xa(n)) then
!					y(i) = ya(n)
					y(i) = ya(n) + (ya(n)-ya(n-1)) / (xa(n)-xa(n-1)) * (x(i) - xa(n))
				else
! In range
						klo = 1
						khi = n
1						if(khi-klo > 1) then
							k = (khi+klo)/2
							if (xa(k) > x(i)) then
								khi = k
							else
								klo = k
							endif					
							go to 1
						endif

						h = xa(khi)-xa(klo)

						if (h == 0.d0) then
							print *, 'bad xa input in spline'
							stop
						endif
						a = (xa(khi)-x(i))/h
						b = (x(i)-xa(klo))/h

						y(i) = a*ya(klo)+b*ya(khi)+((a**3.e0-a)*y2a(klo)+(b**3.e0-b)*y2a(khi))*(h**2.e0)/6.e0
					endif
				endif
			enddo

		end subroutine spline
		
! ---------------------------------------------------------
! This routine returns the inverse error function
! ---------------------------------------------------------
		function dierfc(y)
      implicit real(kind=KIND_FLOAT) (a - h, o - z)
      parameter (&
         qa = 9.16461398268964d-01, &
         qb = 2.31729200323405d-01, &
         qc = 4.88826640273108d-01, &
         qd = 1.24610454613712d-01, &
         q0 = 4.99999303439796d-01, &
         q1 = 1.16065025341614d-01, &
         q2 = 1.50689047360223d-01, &
         q3 = 2.69999308670029d-01, &
         q4 = -7.28846765585675d-02)
      parameter (&
         pa = 3.97886080735226000d+00, &
         pb = 1.20782237635245222d-01, &
         p0 = 2.44044510593190935d-01, &
         p1 = 4.34397492331430115d-01, &
         p2 = 6.86265948274097816d-01, &
         p3 = 9.56464974744799006d-01, &
         p4 = 1.16374581931560831d+00, &
         p5 = 1.21448730779995237d+00, &
         p6 = 1.05375024970847138d+00, &
         p7 = 7.13657635868730364d-01, &
         p8 = 3.16847638520135944d-01, &
         p9 = 1.47297938331485121d-02, &
         p10 = -1.05872177941595488d-01,& 
         p11 = -7.43424357241784861d-02)
      parameter (&
         p12 = 2.20995927012179067d-03, &
         p13 = 3.46494207789099922d-02, &
         p14 = 1.42961988697898018d-02, &
         p15 = -1.18598117047771104d-02, &
         p16 = -1.12749169332504870d-02, &
         p17 = 3.39721910367775861d-03, &
         p18 = 6.85649426074558612d-03, &
         p19 = -7.71708358954120939d-04, &
         p20 = -3.51287146129100025d-03, &
         p21 = 1.05739299623423047d-04, &
         p22 = 1.12648096188977922d-03)
      z = y
      if (y .gt. 1) z = 2 - y
      w = qa - log(z)
      u = sqrt(w)
      s = (qc + log(u)) / w
      t = 1 / (u + qb)
      x = u * (1 - s * (0.5d0 + s * qd)) - &
         ((((q4 * t + q3) * t + q2) * t + q1) * t + q0) * t
      t = pa / (pa + x)
      u = t - 0.5d0
      s = (((((((((p22 * u + p21) * u + p20) * u + &
         p19) * u + p18) * u + p17) * u + p16) * u + &
         p15) * u + p14) * u + p13) * u + p12
      s = ((((((((((((s * u + p11) * u + p10) * u + &
         p9) * u + p8) * u + p7) * u + p6) * u + p5) * u + &
         p4) * u + p3) * u + p2) * u + p1) * u + p0) * t - &
         z * exp(x * x - pb)
      x = x + s * (1 + x * s)
      if (y .gt. 1) x = -x
      dierfc = x
      end function dierfc
      
! ---------------------------------------------------------
! This routine returns the number of standard deviations needed for
! getting A area under the gaussian curve
! ---------------------------------------------------------
		function inverse_normal(A)
		real(kind=KIND_FLOAT) :: inverse_normal, A
			inverse_normal = sqrt(2.)*dierfc(1.0-A)
			return
		end function inverse_normal

		subroutine dtable_latinize ( dim_num, n, table )

!******************************************************************************
!
!! DTABLE_LATINIZE "Latinizes" a real table dataset.
!
!  Discussion:
!
!    It is assumed, though not necessary, that the input dataset
!    has points that lie in the unit hypercube.
!
!    In any case, the output dataset will have this property.
!
!  Modified:
!
!    30 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer N, the number of cells.
!
!    Input/output, real ( kind = 8 ) TABLE(DIM_NUM,N).  On input, the dataset to
!    be "Latinized".  On output, the Latinized dataset.
!
  implicit none

  integer dim_num
  integer n

  integer i
  integer indx(n)
  integer j
  real ( kind=KIND_FLOAT ) table(dim_num,n)

  do i = 1, dim_num
    call r8vec_sort_heap_index_a ( n, table(i,1:n), indx )
    do j = 1, n
      table(i,indx(j)) = real ( 2 * j - 1, kind = KIND_FLOAT ) &
                       / real ( 2 * n    , kind = KIND_FLOAT )
    end do
  end do

  return
end subroutine dtable_latinize

subroutine r8vec_sort_heap_index_a ( n, a, indx )

!*****************************************************************************80
!
!! R8VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of a real vector.
!
!  Discussion:
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(INDX(I)), I = 1 to N is sorted,
!
!    or explicitly, by the call
!
!      call R8VEC_PERMUTE ( N, A, INDX )
!
!    after which A(I), I = 1 to N is sorted.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), an array to be index-sorted.
!
!    Output, integer INDX(N), the sort index.  The
!    I-th element of the sorted array is A(INDX(I)).
!
  implicit none

  integer n

  real ( kind=KIND_FLOAT ) a(n)
  real ( kind=KIND_FLOAT ) aval
  integer i
  integer indx(n)
  integer indxt
  integer ir
  integer j
  integer l

  if ( n < 1 ) then
    return
  end if

  if ( n == 1 ) then
    indx(1) = 1
    return
  end if

  do i = 1, n
    indx(i) = i
  end do

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      aval = a(indxt)

    else

      indxt = indx(ir)
      aval = a(indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if ( a(indx(j)) < a(indx(j+1)) ) then
          j = j + 1
        end if
      end if

      if ( aval < a(indx(j)) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
	end subroutine r8vec_sort_heap_index_a


end module maths_mod
