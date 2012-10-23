;-------------------------------------------------------------
;+
; NAME:
;       FACTOR
; PURPOSE:
;       Find prime factors of a given number.
; CATEGORY:
; CALLING SEQUENCE:
;       factor, x, p, n
; INPUTS:
;       x = Number to factor (>1).       in
; KEYWORD PARAMETERS:
;       Keywords:
;         /QUIET  means do not print factors.
;         /DEBUG  Means list steps as they happen.
;         /TRY    Go beyond 20000 primes.
; OUTPUTS:
;       p = Array of prime numbers.      out
;       n = Count of each element of p.  out
; COMMON BLOCKS:
; NOTES:
;       Note: see also prime, numfactors, print_fact.
; MODIFICATION HISTORY:
;       R. Sterner.  4 Oct, 1988.
;       RES 25 Oct, 1990 --- converted to IDL V2.
;       R. Sterner, 1999 Jun 30 --- Improved (faster, bigger).
;       R. Sterner, 1999 Jul  7 --- Bigger values (used unsigned).
;       R. Sterner, 1999 Jul  9 --- Tried to make backward compatable.
;       R. Sterner, 2000 Jan 06 --- Fixed to ignore non-positive numbers.
;       Johns Hopkins University Applied Physics Laboratory.
;
; Copyright (C) 1988, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;-
;-------------------------------------------------------------
; NAME:
;       SPC
; PURPOSE:
;       Return a string with the specified number of spaces (or other char).
; CATEGORY:
; CALLING SEQUENCE:
;       s = spc(n, [text])
; INPUTS:
;       n = number of spaces (= string length).   in 
;        text = optional text string.              in
;          # spaces returned is n-strlen(strtrim(text,2))
; KEYWORD PARAMETERS:
;       Keywords:
;         CHARACTER=ch  Character other than a space.
;           Ex: CHAR='-'.
;         /NOTRIM means do not do a strtrim on text.
; OUTPUTS:
;       s = resulting string.                     out
; COMMON BLOCKS:
; NOTES:
;       Note: Number of requested spaces is reduced by the
;         length of given string.  Useful for text formatting.
; MODIFICATION HISTORY:
;       Written by R. Sterner, 16 Dec, 1984.
;       RES --- rewritten 14 Jan, 1986.
;       R. Sterner, 27 Jun, 1990 --- added text.
;       R. Sterner, 1994 Sep  7 --- Allowed text arrays.
;       R. Sterner, 1999 Jul  2 --- Added /NOTRIM keyword.
;       Johns Hopkins University Applied Physics Laboratory.
;
; Copyright (C) 1984, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;-------------------------------------------------------------
 
	function spc,n, text, character=char, notrim=notrim, help=hlp
 
	if (n_params(0) lt 1) or keyword_set(hlp) then begin
	  print,' Return a string with the specified number of spaces (or '+$
	    'other char).' 
	  print,' s = spc(n, [text])' 
	  print, '  n = number of spaces (= string length).   in '
	  print,'   text = optional text string.              in'
	  print,'     # spaces returned is n-strlen(strtrim(text,2))'
	  print,'   s = resulting string.                     out' 
	  print,' Keywords:'
	  print,'   CHARACTER=ch  Character other than a space.'
	  print,"     Ex: CHAR='-'."
	  print,'   /NOTRIM means do not do a strtrim on text.'
	  print,' Note: Number of requested spaces is reduced by the'
	  print,'   length of given string.  Useful for text formatting.'
	  return, -1
	endif
 
	if n_params(0) eq 1 then begin
	  n2 = n
	endif else begin
	  if keyword_set(notrim) then $
	    ntxt=strlen(text) else ntxt=strlen(strtrim(text,2))
;	  n2 = n - strlen(strtrim(text,2))
	  n2 = n - ntxt
	endelse
 
	ascii = 32B
	if n_elements(char) ne 0 then ascii = (byte(char))[0]
 
	num = n_elements(n2)
	out = strarr(num)
	for i = 0, num-1 do begin
	  if n2[i] le 0 then out[i] = '' else $
	    out[i] = string(bytarr(n2[i]) + ascii)
	endfor
 
	if n_elements(out) eq 1 then out=out[0]
	return, out
 
	end


;-------------------------------------------------------------
; NAME:
;       PRINT_FACT
; PURPOSE:
;       Print prime factors found by the factor routine.
; CATEGORY:
; CALLING SEQUENCE:
;       print_fact, p, n
; INPUTS:
;       p = prime factors.          in
;       n = number of each factor.  in
; KEYWORD PARAMETERS:
; OUTPUTS:
; COMMON BLOCKS:
; NOTES:
; MODIFICATION HISTORY:
;       R. Sterner  4 Oct, 1988.
;       RES 25 Oct, 1990 --- converted to IDL V2.
;       R. Sterner, 26 Feb, 1991 --- Renamed from print_factors.pro
;       R. Sterner, 1999 Jun 30 --- Better output format.
;       R. Sterner, 1999 Jul  7 --- Bigger values (used unsigned).
;       R. Sterner, 1999 Jul  9 --- Made backward compatable.
;
; Copyright (C) 1988, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;-------------------------------------------------------------
 
	pro print_fact, p, n, help=hlp
 
	if (n_params(0) lt 2) or keyword_set(hlp) then begin
	  print,' Print prime factors found by the factor routine.'
	  print,' print_fact, p, n'
	  print,'   p = prime factors.          in'
	  print,'   n = number of each factor.  in'
	  return
	endif
 
	;-------  Drop unused primes  ---------------
	w = where(n gt 0)	; Find only primes used.
	p2 = p[w]
	n2 = n[w]
 
	;-------  Use largest available integer type  --------------
	flag = !version.release ge 5.2
	if flag eq 1 then begin
	  err=execute('t=1ULL')		; Use 64 bit int (hide from old IDL).
	endif else begin
	  t = 1L			; Use long int (best available in old).
	endelse
 
	;-------  Compute number from it's prime factors.  ----------
	for i = 0, n_elements(p2)-1 do t = t * p2[i]^n2[i]
 
	;-------  Prepare output  -----------------------
	a = strtrim(t,2)+' = '			; Start factors string.
	b = ''					; Start exponents string.
	last = n_elements(p2)-1			; Last factors index.
	for i=0, last do begin
	  a = a + strtrim(p2[i],2)		; Insert next factor.
	  lena = strlen(a)			; Length of factor string.
	  nxtb = strtrim(n2[i],2)		; Next exponent.
	  if nxtb eq '1' then nxtb=' '		; Weed out 1s.
	  b = b+spc(lena,b,/notrim)+nxtb	; Insert next exponent.
	  if i ne last then a=a+' x '		; Not last, add x.
	endfor
 
	;------  Print exponents and factors  -----------
	print,' '
	print,' '+b
	print,' '+a
 
	return
	end


 
	pro factor, x, p, n, quiet=quiet, debug=debug, try=try, help=hlp
 
	if (n_params(0) lt 1) or keyword_set(hlp) then begin
	  print,' Find prime factors of a given number.'
	  print,' factor, x, p, n'
	  print,'   x = Number to factor (>1).       in'
	  print,'   p = Array of prime numbers.      out'
	  print,'   n = Count of each element of p.  out'
	  print,' Keywords:'
	  print,'   /QUIET  means do not print factors.'
	  print,'   /DEBUG  Means list steps as they happen.'
	  print,'   /TRY    Go beyond 20000 primes.'
	  print,' Note: see also prime, numfactors, print_fact.'
	  return
	endif
 
	if x le 0 then return
 
	flag = !version.release ge 5.2
 
	s = sqrt(x)			; Only need primes up to sqrt(x).
	g = long(50 + 0.13457*s)	; Upper limit of # primes up to s.
	np = 50				; Start with np (50) primes.
	p = prime(np)			; Find np primes.
	n = intarr(n_elements(p))	; Divisor count.
 
	if flag eq 1 then $		; Working number.
	  err=execute('t=ulong64(x)') $	; Use best integer available.
	  else t=long(x)		; Best pre-5.2 integer.
	i = 0L				; Index of test prime.
 
loop:	pt = p[i]			; Pull test prime.
	if keyword_set(debug) then $
	  print,' Trying '+strtrim(pt,2)+' into '+strtrim(t,2)
	if flag eq 1 then $
	  err=execute('t2=ulong64(t/pt)') $
	  else t2=long(t/pt)
	if t eq t2*pt then begin	; Check if it divides.
	  if keyword_set(debug) then $
	    print,'   Was a factor.  Now do '+strtrim(t2,2)
	  n[i] = n[i] + 1		; Yes, count it.
	  t = t2			; Result after division.
	  if t2 eq 1 then goto, done	; Check if done.
	  goto, loop			; Continue.
	endif else begin
	  i = i + 1			; Try next prime.
	  if i ge np then begin
	    s = sqrt(t)			; Only need primes up to sqrt(x).
	    g = long(50 + 0.13457*s)	; Upper limit of # primes up to s.
	    if g le np then goto, last	; Must be done.
	    np = (np+50)<g		; Want 50 more primes.
	    if (np gt 20000) and (not keyword_set(try)) then begin
	      print,' Too hard.  Tried '+strtrim(np-50,2)+' primes.'
	      print,' Trying to crack '+strtrim(t,2)
	      print,' To go farther use keyword /TRY.'
	      return
	    endif
	    if keyword_set(debug) then $
	      print,' Finding more primes: '+strtrim(np,2)+ $
	      '.  Max needed = '+strtrim(g,2)
	    p = prime(np)		; Find primes.
	    n = [n,intarr(50)]		; Make room for more factors.
	  endif
	  if i ge g then goto, last	; Nothing up to sqrt works.
	  goto, loop			; Continue.
	endelse
 
last:	p = [p,t]			; Residue was > sqrt, must be prime.
	n = [n,1]			; Must occur only once. (else < sqrt).
 
done:	w = where(n gt 0)
	n = n[w]			; Trim excess off tables.
	p = p[w]
 
	if not keyword_set(quiet) then print_fact, p, n
 
	return
	end

;-------------------------------------------------------------
;+
; NAME:
;     PRIME
; PURPOSE:
;     Return an array with the specified number of prime numbers.
; EXPLANATATION:
;     This procedure is similar to PRIMES in the standard IDL distribution,
;     but stores results in a common block, and so is much faster 
;
; CALLING SEQUENCE:
;       p = prime(n)
; INPUTS:
;       n = desired number of primes, scalar positive integer
; OUTPUTS:
;       p = resulting array of primes, vector of positive integers
; COMMON BLOCKS:
;       prime_com
; NOTES:
;       Note: Primes that have been found in previous calls are
;         remembered and are not regenerated.
; MODIFICATION HISTORY:
;       R. Sterner  17 Oct, 1985.
;       R. Sterner,  5 Feb, 1993 --- fixed a bug that missed a few primes.
;       Converted to IDL V5          March 1999
;
; Copyright (C) 1985, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;-
;-------------------------------------------------------------
 
	function prime,n, help=hlp
 
	common prime_com, max, pmax
 
	if (n_params(0) lt 1) or keyword_set(hlp) then begin
	  print,' Return an array with the specified number of prime numbers.'
	  print,' p = prime(n)'
	  print,'   n = desired number of primes.    in'
	  print,'   p = resulting array of primes.   out'
	  print,' Note: Primes that have been found in previous calls are'
	  print,'   remembered and are not regenerated.'
	  return, -1
	endif
 
	if n_elements(max) eq 0 then max = 0	; Make MAX defined.
	if n le max then return, pmax[0:n-1]	; Enough primes in memory.
	p = lonarr(n)				; Need to find primes.
	if max eq 0 then begin			; Have none now. Start with 8.
	  p[0] = [2,3,5,7,11,13,17,19]
	  if n le 8 then return, p[0:n-1]	; Need 8 or less.
	  i = 8					; Need more than 8.
	  t = 19L				; Search start value.
	endif else begin			; Start with old primes.
	  p[0] = pmax				; Move old primes into big arr.
	  i = max				; Current prime count.
	  t = p[max-1]				; Biggest prime so far.
	endelse
 
loop:	if i eq n then begin			; Have enough primes.
	  max = n				; Count.
	  pmax = p				; Array of primes.
	  return, p				; Return primes.
	endif
loop2:	t = t + 2				; Next test value, t.
	it = 1					; Start testing with 1st prime.
loop3:	pr = p[it]				; Pick next test prime.
	pr2 = pr*pr				; Square it.
	if pr2 gt t then begin			; Selected prime > sqrt(t)?
	  i = i + 1				; Yes, count
	  p[i-1] = t				; and store new prime.
	  goto, loop				; Go check if done.
	endif
	if pr2 eq t then goto, loop2		; Test number, t, was a square.
	if (t mod pr) eq 0 then goto, loop2	; Curr prime divides t.
	it = it + 1				; Check next prime.
	goto, loop3
	end

function convolve, image, psf, FT_PSF=psf_FT, FT_IMAGE=imFT, NO_FT=noft, $
			CORRELATE=correlate, AUTO_CORRELATION=auto
;+
; NAME:
;	CONVOLVE
; PURPOSE:
;	Convolution of an image with a Point Spread Function (PSF)
; EXPLANATION:
;	The default is to compute the convolution using a product of 
;	Fourier transforms (for speed).
;
; CALLING SEQUENCE:
;
;	imconv = convolve( image1, psf, FT_PSF = psf_FT )
;  or:
;	correl = convolve( image1, image2, /CORREL )
;  or:
;	correl = convolve( image, /AUTO )
;
; INPUTS:
;	image = 2-D array (matrix) to be convolved with psf
;	psf = the Point Spread Function, (size < or = to size of image).
;
; OPTIONAL INPUT KEYWORDS:
;
;	FT_PSF = passes out/in the Fourier transform of the PSF,
;		(so that it can be re-used the next time function is called).
;	FT_IMAGE = passes out/in the Fourier transform of image.
;
;	/CORRELATE uses the conjugate of the Fourier transform of PSF,
;		to compute the cross-correlation of image and PSF,
;		(equivalent to IDL function convol() with NO rotation of PSF)
;
;	/AUTO_CORR computes the auto-correlation function of image using FFT.
;
;	/NO_FT overrides the use of FFT, using IDL function convol() instead.
;		(then PSF is rotated by 180 degrees to give same result)
; METHOD:
;	When using FFT, PSF is centered & expanded to size of image.
; HISTORY:
;	written, Frank Varosi, NASA/GSFC 1992.
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
	sp = size( psf_FT )  &  sif = size( imFT )
	sim = size( image )  &  sc = sim/2  &  npix = N_elements( image )

	if (sim[0] NE 2) OR keyword_set( noft ) then begin
		if keyword_set( auto ) then begin
			message,"auto-correlation only for images with FFT",/INF
			return, image
		  endif else if keyword_set( correlate ) then $
				return, convol( image, psf ) $
			else	return, convol( image, rotate( psf, 2 ) )
	   endif

	if (sif[0] NE 2) OR (sif[sif[0]+1] NE 6) OR $
	   (sif[1] NE sim[1]) OR (sif[2] NE sim[2]) then imFT = FFT( image,-1 )

	if keyword_set( auto ) then $
	 return, shift( npix*float( FFT( imFT*conj( imFT ),1 ) ), sc[1],sc[2] )

	if (sp[0] NE 2) OR (sp[sp[0]+1] NE 6) OR $
	   (sp[1] NE sim[1]) OR (sp[2] NE sim[2]) then begin
		sp = size( psf )
		if (sp[0] NE 2) then begin
			message,"must supply PSF matrix (2nd arg.)",/INFO
			return, image
		   endif
		Loc = ( sc - sp/2 ) > 0		;center PSF in new array,
		s = (sp/2 - sc) > 0	   ;handle all cases: smaller or bigger
		L = (s + sim-1) < (sp-1)
		psf_FT = complexarr( sim[1], sim[2] )
		psf_FT[ Loc[1], Loc[2] ] = psf[ s[1]:L[1], s[2]:L[2] ]
		psf_FT = FFT( psf_FT, -1, /OVERWRITE )
	   endif

	if keyword_set( correlate ) then $
		conv = npix * float( FFT( imFT * conj( psf_FT ), 1 ) ) $
	  else	conv = npix * float( FFT( imFT * psf_FT, 1 ) )

	sc = sc + (sim MOD 2)	;shift correction for odd size images.

return, shift( conv, sc[1], sc[2] )
end

function psf_gaussian, parameters, NPIXEL=npix, NDIMENSION=ndim, FWHM=fwhm,  $
                        DOUBLE = double, CENTROID=cntrd, ST_DEV=st_dev,  $
                        XY_CORREL=xy_corr, NORMALIZE=normalize
;+
; NAME:
;       PSF_GAUSSIAN
;
; PURPOSE:
;       Create a 1-d, 2-d, or 3-d Gaussian with specified FWHM, center 
; EXPLANATION:
;       Return a point spread function having Gaussian profiles,
;       as either a 1D vector, a 2D image, or 3D volumetric-data.
;
; CALLING SEQUENCE:
;       psf = psf_Gaussian( NPIXEL=, FWHM= , CENTROID = 
;                     [ /DOUBLE, /NORMALIZE, ST_DEV=,  NDIMEN= ] ) 
; or:
;       psf = psf_Gaussian( parameters, NPIXEL = ,NDIMEN = )
;
; REQUIRED INPUT KEYWORD:
;       NPIXEL = number pixels for each dimension, specify as an array,
;               or just one number to make all sizes equal.
;
; OPTIONAL KEYWORDS:
;       CENTROID = floating scalar or vector giving position of  PSF center.    
;               default is exact center of requested vector/image/volume.
;               The number of elements in CENTROID should equal the number of
;               dimensions.    **The definition of Centroid was changed in
;               March 2002, and now an integer defines the center of a pixel.**
;
;       /DOUBLE  = If set, then the output array is computed in double precision
;               the default is to return a floating point array.
;
;       FWHM = the desired Full-Width Half-Max (pixels) in each dimension,
;               specify as an array, or single number to make all the same.
;
;       NDIMEN = integer dimension of result: either 1 (vector), 2 (image), or 
;                3 (volume), default = 2 (an image result).
;
;       /NORMALIZE causes resulting PSF to be normalized so Total( psf ) = 1.
;
;       ST_DEV = optional way to specify width by standard deviation param.
;                Ignored if FWHM is specified.
;
;       XY_CORREL = scalar between 0 and 1 specifying correlation coefficient
;               Use this keyword, for example, to specify an elliptical 
;               Gaussian oriented at an angle to the X,Y axis.   Only valid
;               for 2-dimensional case.
;
;
; INPUTS (optional):
;
;       parameters = an NDIMEN by 3 array giving for each dimension:
;                       [ maxval, center, st_dev ],  overrides other keywords.
;
; EXAMPLE:
;       (1) Create a 31 x 31 array containing a normalized centered Gaussian 
;       with an X FWHM = 4.3 and a Y FWHM = 3.6
;
;       IDL> array = PSF_GAUSSIAN( Npixel=31, FWHM=[4.3,3.6], /NORMAL )
;
;       (2) Create a 50 pixel 1-d Gaussian vector with a maximum of 12, 
;          centered at  pixel 23 with a sigma of 19.2
;
;       IDL> psf = psf_gaussian([12,23,19.2],npixel=50)
; EXTERNAL CALLS:
;       function Gaussian()
; NOTES:
;       To improve speed, floating underflow exceptions are suppressed (using 
;       the MASK=32  keyword of CHECK_MATH() rather than being flagged.
;
; HISTORY:
;       Written, Frank Varosi NASA/GSFC 1991.
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Suppress underflow messages, add DOUBLE keyword. **Modified centroid
;       definition so integer position is pixel center** W. Landsman March 2002
;       Allow use of the ST_DEV (not STDEV) keyword W. Landsman Nov. 2002
;-
        On_error,2

        if (N_params() LT 1 ) and $
            not (keyword_set( FWHM) or keyword_set(ST_DEV)) then begin
                print,'Syntax - psf = PSF_GAUSSIAN( parameters, NPIXEL = )'
                print, $
       'or       psf = PSF_GAUSSIAN( FWHM = ,ST_DEV = ,NPIXEL = ,[CENTROID = ])'
                return, -1
        endif

        sp = size( parameters )
        if sp[0] EQ 1 then begin               ;Vector supplied?
                ndim = 1
                factor = parameters[0]
                cntrd = parameters[1]
                st_dev = parameters[2] 
         endif  else  if (sp[0] GE 1) then begin    ;Ndimen x 3 array supplied?
                 ndim = sp[1]
                 factor = total( parameters[*,0] )/float( ndim )
                cntrd = parameters[*,1]
                st_dev = parameters[*,2]
           endif

        double = keyword_set(double)
        if double then idltype = 5 else idltype = 4
        if N_elements( ndim ) NE 1 then ndim=2
        ndim = ndim>1

        if N_elements( npix ) LE 0 then begin
                message,"must specify size of result with NPIX=",/INFO
                return,(-1)
          endif else if N_elements( npix ) LT ndim then $
                        npix = replicate( npix[0], ndim )

        if (N_elements( cntrd ) LT ndim) AND (N_elements( cntrd ) GT 0) then $
                        cntrd = replicate( cntrd[0], ndim )

        if N_elements( cntrd ) LE 0 then cntrd=(npix-1)/2. 
        if N_elements( fwhm ) GT 0 then begin 
               st_dev = fwhm/( 2.0d* sqrt( 2.0d* aLog(2.0d) ) )
               if not double then st_dev  = float(st_dev)
        endif 

        if N_elements( st_dev ) LE 0 then begin
                message,"must specify ST_DEV= or FWHM=",/INFO
                return,(-1)
          endif

        if N_elements( st_dev ) LT ndim then $
                        st_dev = replicate( st_dev[0], ndim )

        CASE ndim OF

        1: BEGIN
                x = findgen( npix[0] ) - cntrd[0]
                psf = gaussian( x, [1,0,st_dev] )
             END

        2: BEGIN
                psf = make_array( DIM=npix[0:ndim-1], TYPE = idltype )
                x = make_array( npix[0], /INDEX, TYPE=idltype ) - cntrd[0]
                y = make_array( npix[1], /INDEX, TYPE=idltype ) - cntrd[1]

                if N_elements( xy_corr ) EQ 1 then begin
                        sigfac = 1 / (2. * st_dev^2 )
                        y2 = sigfac[1] * y^2
                        x1 = sigfac[0] * x
                        yc = y * ( xy_corr/(st_dev[0]*st_dev[1]) )
                        for j=0,npix[1]-1 do begin
                                zz = x * (yc[j] + x1) + y2[j]
                                w = where( zz LT 86, nw )
                                if (nw GT 0) then psf[w,j] = exp( -zz[w] )
                          endfor
                  endif else begin
                        psfx = gaussian( x, [ 1, 0, st_dev[0] ], DOUBLE=double )
                        psfy = gaussian( y, [ 1, 0, st_dev[1] ], DOUBLE=double )
                        error = check_math(/print, MASK=32)
                        save_except = !EXCEPT & !EXCEPT = 0
                        for j=0,npix[1]-1 do psf[0,j] = psfx * psfy[j]
                        error = check_math(MASK=32)    ;Clear floating underflow
                        !EXCEPT = save_except  
                   endelse
             END

        3: BEGIN
                psf = make_array( DIM=npix[0:ndim-1], TYPE = idltype )
                x = make_array( npix[0], /INDEX, TYPE=idltype ) - cntrd[0]
                y = make_array( npix[1], /INDEX, TYPE=idltype ) - cntrd[1]
                z = make_array( npix[2], /INDEX, TYPE=idltype ) - cntrd[2]
                psfx = gaussian( x, [ 1, 0, st_dev[0] ], DOUBLE = double )
                psfy = gaussian( y, [ 1, 0, st_dev[1] ], DOUBLE = double)
                psfz = gaussian( z, [ 1, 0, st_dev[2] ], DOUBLE = double )
                error = check_math(MASK=32,/PRINT)
                save_except = !EXCEPT & !EXCEPT = 0
                for k=0,npix[2]-1 do begin
                    for j=0,npix[1]-1 do psf[0,j,k] = psfx * psfy[j] * psfz[k]
                 endfor
                 error = check_math(MASK=32)
                 !EXCEPT = save_except  
             END

        ENDCASE

        if keyword_set( normalize ) then return, psf/total( psf )

        if N_elements( factor ) EQ 1 then begin
                if (factor NE 1) then return,factor*psf else return,psf
           endif else return, psf
end

function filter_image, image, SMOOTH=width_smooth, ITERATE_SMOOTH=iterate, $
                              MEDIAN=width_median, ALL_PIXELS=all_pixels, $
                              FWHM_GAUSSIAN=fwhm, NO_FT_CONVOL=no_ft, PSF=psf
;+
; NAME:
;       FILTER_IMAGE
;
; PURPOSE:
;       Identical to MEDIAN or SMOOTH but handle edges and allow iterations.
; EXPLANATION:
;       Computes the average and/or median of pixels in moving box,
;       replacing center pixel with the computed average and/or median,
;       (using the IDL SMOOTH() or MEDIAN() functions).
;       The main reason for using this function is the options to
;       also process the pixels at edges and corners of image, and,
;       to apply iterative smoothing simulating convolution with Gaussian,
;       and/or to convolve image with a Gaussian kernel.
;
; CALLING SEQUENCE:
;       Result = filter_image( image, SMOOTH=width, MEDIAN = width, /ALL_PIXELS
;                               /ITERATE, FWHM =,  /NO_FT_CONVOL)
;
; INPUT:
;       image = 2-D array (matrix)
;
; OPTIONAL INPUT KEYWORDS:
;       SMOOTH = scalar (odd) integer specifying the width of a square box 
;               for moving average, in # pixels.  /SMOOTH  means use box 
;               width = 3 pixels for smoothing.
;
;        MEDIAN = scalar (usually odd) integer specifying the width of square 
;               moving box for median filter, in # pixels.   /MEDIAN  means use
;               box width = 3 pixels for median filter.
;   
;       /ALL_PIXELS causes the edges of image to be filtered as well.   This
;               is accomplished by reflecting pixels adjacent to edges outward
;               (similar to the /EDGE_WRAP keyword in CONVOL).
;               Note that this is a different algorithm from the /EDGE_TRUCATE 
;               keyword to SMOOTH or CONVOL, which duplicates the nearest pixel.   
;
;       /ITERATE means apply smooth(image,3) iteratively for a count of
;               (box_width-1)/2 times (=radius), when box_width >= 5.
;               This is equivalent to convolution with a Gaussian PSF
;               of FWHM = 2 * sqrt( radius ) as radius gets large.
;               Note that /ALL_PIXELS is automatically applied,
;               giving better results in the iteration limit.
;               (also, MEDIAN keyword is ignored when /ITER is specified).
;
;       FWHM_GAUSSIAN = Full-width half-max of Gaussian to convolve with image. 
;                       FWHM can be a single number (circular beam),
;                       or 2 numbers giving axes of elliptical beam.
;
;       /NO_FT_CONVOL causes the convolution to be computed directly,
;               with intrinsic IDL CONVOL function.   The default is to use 
;               FFT when factors of size are all LE 13.   Note that 
;               external function convolve.pro handles both cases)
;
; OPTIONAL INPUT/OUTPUT KEYWORD:
;     PSF = Array containing the PSF used during the convolution.   This 
;           keyword is only active if the FWHM_GAUSSIAN keyword is also 
;           specified.     If PSF is undefined on input, then upon output it
;           contains the Gaussian convolution specified by the FWHM_GAUSSIAN
;           keyword.    If the PSF array is defined on input then it is used 
;           as the convolution kernel,  the value of the  FWHM_GAUSSIAN keyword
;           is ignored.      Typically, on a first call set PSF to an undefined
;           variable, which can be reused for subsequent calls to prevent 
;           recalculation of the Gaussian PSF.
; RESULT:
;       Function returns the smoothed, median filtered, or convolved image.
;       If both SMOOTH and MEDIAN are specified, median filter is applied first.
;
; EXAMPLES:
;       To apply 3x3 moving median filter and
;       then 3x3 moving average, both applied to all pixels:
;
;               Result = filter_image( image, /SMOOTH, /MEDIAN, /ALL )
;
;       To iteratively apply 3x3 moving average filter for 4 = (9-1)/2 times,
;       thus approximating convolution with Gaussian of FWHM = 2*sqrt(4) = 4 :
;
;               Result = filter_image( image, SMOOTH=9, /ITER )
;
;       To convolve all pixels with Gaussian of FWHM = 3.7 x 5.2 pixels:
;
;               Result = filter_image( image, FWHM=[3.7,5.2], /ALL )
;
; EXTERNAL CALLS:
;       function psf_gaussian
;       function convolve
;       pro factor
;       function prime          ;all these called only if FWHM is specified
;
; PROCEDURE:
;       If both /ALL_PIXELS (or /ITERATE)  keywords are set then
;       create a larger image by reflecting the edges outward, then call the 
;       IDL MEDIAN() or SMOOTH() function on the larger image, and just return 
;       the central part (the original size image).
;
;       NAN values are recognized during calls to MEDIAN() or SMOOTH(), but 
;       not for convolution with a Gaussian (FWHM keyword supplied). 
; HISTORY:
;       Written, 1991, Frank Varosi, NASA/GSFC.
;       FV, 1992, added /ITERATE option.
;       FV, 1993, added FWHM_GAUSSIAN= option.
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Use /EVEN call to median, recognize NAN values in SMOOTH 
;                  W. Landsman   June 2001
;       Added PSF keyword,   Bjorn Heijligers/WL, September 2001
;-

  if N_params() LT 1 then begin
      print,'Syntax - Result = filter_image( image, SMOOTH=width, /ALL_PIXELS'
      print,'                 MEDIAN= width, ITERATE, FWHM=,  /NO_FT_CONVOL'
      return, -1
  endif

        sim = size( image )
        Lx = sim[1]-1
        Ly = sim[2]-1

        if (sim[0] NE 2) OR (sim[4] LE 4) then begin
                message,"input must be an image (a matrix)",/INFO
                return,image
           endif

        if keyword_set( iterate ) then begin
                if N_elements( width_smooth ) NE 1 then return,image
                if (width_smooth LT 1) then return,image
                imf = image
                nit = (width_smooth>3)/2
                for i=1,nit do  imf = filter_image( imf, /SMOOTH, /ALL )
                return,imf
           endif

        box_wid = 0
        if keyword_set( width_smooth ) then box_wid = width_smooth > 3
        if keyword_set( width_median ) then box_wid = (width_median > box_wid)>3

        if keyword_set( fwhm ) then begin
                npix = ( 3 * fwhm[ 0: ( (N_elements( fwhm )-1) < 1 ) ] ) > 3
                npix = 2 * fix( npix/2 ) + 1    ;make # pixels odd.
                box_wid = box_wid > max( [npix] )
           endif

        if (box_wid LT 3) then return, image

        if keyword_set(all_pixels) then begin
                
                box_wid = fix( box_wid )
                radius = (box_wid/2) > 1
                Lxr = Lx+radius
                Lyr = Ly+radius
                rr = 2*radius
                imf = fltarr( sim[1]+rr, sim[2]+rr )
                imf[radius,radius] = image              ; reflect edges outward
                                                        ; to make larger image.
                imf[  0,0] = rotate( imf[radius:rr,*], 5 )      ;Left
                imf[Lxr,0] = rotate( imf[Lx:Lxr,*], 5 )         ;right
                imf[0,  0] = rotate( imf[*,radius:rr], 7 )      ;bottom
                imf[0,Lyr] = rotate( imf[*,Ly:Lyr], 7 )         ;top

          endif else begin

                radius=0
                imf = image
           endelse

        if keyword_set( width_median ) then $
                       imf = median(/even, imf, width_median>3 ) 
                            
        if keyword_set( width_smooth ) then $
              imf = smooth( imf, width_smooth>3, /NAN )

        if keyword_set( fwhm ) then begin

                if N_elements( no_ft ) NE 1 then begin
                        sim = size( imf )
                        factor,sim[1],pfx,nfx,/quiet
                        factor,sim[2],pfy,nfy,/quiet
                        no_ft = max( [pfx,pfy] ) GT 13
                   endif

                if N_elements(PSF) EQ 0 then $
                          psf=psf_gaussian( NP=npix,FWHM=fwhm,/NORM )
                imf = convolve( imf,  NO_FT=no_ft, psf) 
          endif

    if radius GT 0 then $
                return, imf[ radius:(Lx+radius), radius:(Ly+radius) ] $
           else return, imf
end


;+
; NAME:
;        TVFRAME
; PURPOSE:
;        Display an image with axis.
; CATEGORY:
;        Image dispaly
; CALLING SEQUENCE:
;        TVBOX,data
; INPUTS:
;        data = Two dimensional array of numerical type
; KEYWORD PARAMETERS:
;     XRANGE   : Array with at least 2 elements giving the range for
;                the annotation of the x-axis.
;                Only min(XRANGE) and max(XRANGE) are used.
;                Default is annotation with pixel indices.
;                (should not be confused with the standard plotting keword)
;     YRANGE   : Same as XRANGE but for y-axis.
;     POSITION : (output) position parameter for the axis in device units,
;                may be used in subsequent calls to PLOT or CONTOUR.
;                Example : TVBOX,a,/aspect,position=p
;                          CONTOUR,a,nlev=10,xsty=5,ysty=5,pos=p,/dev,/noeras
;                          (no keyword like /ASPECT exists for CONTOUR)
;     /sample  : If set, nearest neighbourhood method is used in
;                resizing the data array.
;                Default : averaging or linear interpolation.
;                (Corresponds to the SAMPLE-keyword of the REBIN function)
;     /center  : If set, tickmarks are placed at the center of each pixel.
;                ( example: tvframe,randomu(seed,20,20),/sam [,/center] )
;     /aspect  : If set, the aspect ratio of the data array is preserved
;                (quadratic pixel shapes).
;     /noscale : If set, data is not scaled. (much like TV compared to TVSCL)
;     /bar     : If set, an intensity bar is displayed.
;     BRANGE   : Range for annotation of the intensity bar.
;                Default is the range of data.
;     BTITLE   : Title of the intensity bar.
;     BTICKS,BMINOR : Control the number of tickmarks for the intensity bar.
; Standard plotting keywords :
;      TITLE,SUBTITLE,XTITLE,YTITLE,TICKLEN,CHARSIZE,XSTYLE,YSTYLE,
;      XTICKS,YTICKS,XMINOR,YMINOR
;                Normally they are just passed through,
;                with the following exceptions :
;                the default value for TICKLEN is -.01 instead of .02 .
;                bit 0 of X/Y-STYLE is allways set (exact axis range).
;                bit 1 of X/Y-STYLE is  never  set (no extended axis range).
; OUTPUTS:
;        POSITION of the plot frame.
; COMMON BLOCKS:
;        none.
; SIDE EFFECTS:
;        none.
; RESTRICTIONS:
;        none.
; PROCEDURE:
;        The image is scaled to the size of the axis-box via the
;        RESIZE-function, and then displayed using TV or TVSCL.
; MODIFICATION HISTORY:
;        Written, A. Welz, Univ. Wuerzburg, Germany, Feb. 1991
;        Extended (/bar,xticks,..) by A.W. June 1992
;	 	   A. Asensio Ramos November 2002 - Included the use of a top in the used color table
;	 	   A. Asensio Ramos June 2009 - Included btickformat
;-      
;
pro tvframe , a    $
    , sample=sample , center=center , aspect=aspect , noscale=noscale    $
    , POSITION=POSITION    $
    , XRANGE=XRANGI , YRANGE=YRANGI   $
    , TITLE=TITLE , XTITLE=XTITLE , YTITLE=YTITLE , SUBTITLE=SUBTITLE   $
    , TICKLEN=TICKLEN , CHARSIZE=CHARSIZE   $
    , XTICKS=XTICKS , YTICKS=YTICKS , XMINOR=XMINOR , YMINOR=YMINOR $
    , XSTYLE=XSTYLI , YSTYLE=YSTYLI   $
    , bar=bar , BTITLE=BTITLE , BRANGE=BRANGE  $
    , BTICKS=BTICKS , BMINOR=BMINOR, TOP=TOP, BOTTOM=BOTTOM, BBOTTOM=BBOTTOM, XCHARSIZE=XCHARSIZE, YCHARSIZE=YCHARSIZE $
	 , BCHARSIZE=BCHARSIZE,XTICKFORMAT=XTICKFORMAT, YTICKFORMAT=YTICKFORMAT $
	 , XLOG=XLOG, YLOG=YLOG, WINPOS=WINPOS, XTICKV=XTICKV, NOERASE=NOERASE $
	 , XTICKNAME=XTICKNAME, INVISIBLEBAR=INVISIBLEBAR, BTICKFORMAT=BTICKFORMAT

on_error,2
a = reform(a)
sa=size(reform(a)) & if sa(0) ne 2 then goto,errout
mina=min(a,max=maxa)
;
; set keyword parameters to default values if not present
if n_elements(XRANGI) eq 0 then XRANGI=[0,sa(1)-1]
if n_elements(YRANGI) eq 0 then YRANGI=[0,sa(2)-1]
if n_elements(   TITLE) eq 0 then    TITLE=''
if n_elements(SUBTITLE) eq 0 then SUBTITLE=''
if n_elements(  XTITLE) eq 0 then   XTITLE=''
if n_elements(  YTITLE) eq 0 then   YTITLE=''
if n_elements(TICKLEN) eq 0 then TICKLEN=-.01
if n_elements(XTICKS) eq 0 then XTICKS=0
if n_elements(XMINOR) eq 0 then XMINOR=0
if n_elements(YTICKS) eq 0 then YTICKS=0
if n_elements(YMINOR) eq 0 then YMINOR=0
if n_elements(CHARSIZE) eq 0 then CHARSIZE=1.0
if n_elements(XCHARSIZE) eq 0 then XCHARSIZE=0
if n_elements(YCHARSIZE) eq 0 then YCHARSIZE=0
if n_elements(BCHARSIZE) eq 0 then BCHARSIZE=CHARSIZE
XSTYLE=1  &  YSTYLE=1
if n_elements(XSTYLI) eq 1 then XSTYLE=( 1 or XSTYLI ) and 29
if n_elements(YSTYLI) eq 1 then YSTYLE=( 1 or YSTYLI ) and 29
if n_elements(BTITLE) eq 0 then BTITLE=''
if n_elements(BRANGE) eq 0 then BRANGE=float([mina,maxa])
if n_elements(BTICKS) eq 0 then BTICKS=0
if n_elements(BMINOR) eq 0 then BMINOR=0
if n_elements(NOERASE) eq 0 then NOERASE=0
if n_elements(BOTTOM) eq 0 then BOTTOM=mina
if n_elements(BBOTTOM) eq 0 then BBOTTOM=0
if n_elements(TOP) eq 0 then mcol=255
;
XRANGE=float(minmax(XRANGI))
YRANGE=float(minmax(YRANGI))
;
if keyword_set(center) then begin
    xunit=0.5*(XRANGE(1)-XRANGE(0))/float(sa(1)-1)
    yunit=0.5*(YRANGE(1)-YRANGE(0))/float(sa(2)-1)
    XRANGE(0)=XRANGE(0)-xunit  &  XRANGE(1)=XRANGE(1)+xunit
    YRANGE(0)=YRANGE(0)-yunit  &  YRANGE(1)=YRANGE(1)+yunit
endif else begin
; CHANGE: Andres Asensio Ramos (24/02/03)
    xunit=(XRANGE(1)-XRANGE(0))/float(sa(1)-1)
    yunit=(YRANGE(1)-YRANGE(0))/float(sa(2)-1)
    XRANGE(1)=XRANGE(1);+xunit
    YRANGE(1)=YRANGE(1);+yunit
endelse
;
if (keyword_set(WINPOS)) then begin
plot,XRANGE,YRANGE,/nodata,xstyle=xstyle or 4,ystyle=ystyle or 4    $
     ,  TITLE=' ',XTITLE=XTITLE,YTITLE=YTITLE    $
     ,  SUBTITLE=SUBTITLE      $
     ,  TICKLEN=TICKLEN,CHARSIZE=CHARSIZE   $
     ,  color=!p.background $
	  ,  XCHARSIZE=XCHARSIZE, YCHARSIZE=YCHARSIZE, XTICKFORMAT=XTICKFORMAT, YTICKFORMAT=YTICKFORMAT $
	  ,  XLOG=XLOG, YLOG=YLOG, POSITION=WINPOS, NOERASE=NOERASE;, XTICKV=XTICKV, XTICKS=XTICKS
endif else begin
plot,XRANGE,YRANGE,/nodata,xstyle=xstyle or 4,ystyle=ystyle or 4    $
     ,  TITLE=' ',XTITLE=XTITLE,YTITLE=YTITLE    $
     ,  SUBTITLE=SUBTITLE      $
     ,  TICKLEN=TICKLEN,CHARSIZE=CHARSIZE   $
     ,  color=!p.background $
	  ,  XCHARSIZE=XCHARSIZE, YCHARSIZE=YCHARSIZE, XTICKFORMAT=XTICKFORMAT, YTICKFORMAT=YTICKFORMAT $
	  ,  XLOG=XLOG, YLOG=YLOG, XTICKV=XTICKV,NOERASE=NOERASE,XTICKNAME=XTICKNAME
endelse
	  

;
px = !x.window * !d.x_vsize     ;Position of frame in device units
py = !y.window * !d.y_vsize
sx = px(1)-px(0)                ;Size of frame in device units
sy = py(1)-py(0)
if keyword_set(bar) or keyword_set(invisiblebar) then sx = sx/1.25
if keyword_set(aspect) then begin
     f = float(sa(1))/sa(2)*sy/sx
     if f ge 1. then sy=sy/f else sx=sx*f
     sx=fix(sx)
endif
POSITION = [px(0),py(0),px(0)+sx,py(0)+sy]

; Also add the possibility to add an invisible bar
if keyword_set(bar) or keyword_set(invisiblebar) then begin
   bx    = fix(px(0)+sx*1.04)
   by    = fix(py(0))
   bsx   = fix(sx*0.08)
   bsy   = fix(sy)
   barpos= [bx,by,bx+bsx,by+bsy]
endif
;
; If we have a top color, reduce the number of colors used
mcol = ( !D.N_COLORS - 1) > 0       ; changed 29/09/09
mcol = 255
if (n_elements(top) eq 1) then begin
	 if (top lt !D.N_COLORS) then begin
	 	  mcol=( top - 1) > 0
	 endif
endif
;
if (!d.flags and 1) ne 0 then begin

      ;  scalable pixels	
; CHANGE: Andres Asensio Ramos (24/02/03)
;     if keyword_set(sample) then b=a else b=resize(a,256,256)
	  if keyword_set(sample) then b=a else b=congrid(a,256,256,/interp,/minus_one)
     if keyword_set(noscale) then begin
         tv, bbottom>b<mcol ,px(0),py(0),xsize=sx,ysize=sy,/device
     endif else begin
     		renorm = (1.d0*b-mina) / (maxa-mina) * (mcol-bbottom)+bbottom
         tv, renorm ,px(0),py(0),xsize=sx,ysize=sy,/device,top=mcol
     endelse
     if keyword_set(bar) then begin
        barim=findgen(1,256)/255*(maxa-mina)+mina
        renorm = (1.d0*barim-mina) / (maxa-mina) * (mcol-bbottom)+bbottom
        if keyword_set(noscale) then begin
            tv, bbottom>barim<mcol ,bx,by,xsize=bsx,ysize=bsy,/device
        endif else begin
            tv, renorm,bx,by,xsize=bsx,ysize=bsy,/device,top=mcol
        endelse
     endif

endif else begin

      ;  not scalable pixels

     if sx*sy gt 10e6 then begin
        print,' do you really want to allocate ',sx*sy   $
             ,' words of memory ? [y/n]'
        answer='n'   &  read,answer
        if answer ne 'y' then return
     endif
; CHANGE: Andres Asensio Ramos (24/02/03)
;     if keyword_set(sample) then b=resize(a,sx,sy,/sample)  $
;        else b=resize(a,sx,sy)		
     if keyword_set(sample) then b=congrid(a,sx,sy,/minus_one)  $
        else b=congrid(a,sx,sy,/interp,/minus_one)		  
     if keyword_set(noscale) then begin
        tv, 0>b<mcol ,px(0),py(0),/device
     endif else begin
     	  renorm = (1.d0*b-mina) / (maxa-mina) * (mcol-bbottom)+bbottom 
        tv, renorm,px(0),py(0),/device,top=mcol
     endelse
     if keyword_set(bar) then begin
        barim=findgen(1,bsy)/(bsy-1)*(maxa-mina)+mina
        barim=rebin(barim,bsx,bsy,/sample)
        renorm = (1.d0*barim-mina) / (maxa-mina) * (mcol-bbottom)+bbottom        
        if keyword_set(noscale) then begin
           tv, 0>barim<mcol ,bx,by,/device
        endif else begin
           tv, renorm,bx,by,/device,top=mcol
        endelse
     endif

endelse
;
if keyword_set(bar) then begin
   plot,[0,1],BRANGE,/nodata,/noerase,pos=barpos,/device,xsty=5,ysty=5
   plots,[bx+bsx,bx,bx,bx+bsx],[by,by,by+bsy,by+bsy],/device
   axis,yaxis=1,bx+bsx,/device,yrange=BRANGE,ystyle=1   $
     ,  ytitle=BTITLE $
     ,  TICKLEN=-.15,CHARSIZE=BCHARSIZE,YTICKS=BTICKS,YMINOR=BMINOR,ytickformat=btickformat
endif
;
plot,XRANGE,YRANGE,/nodata,/noerase,xstyle=XSTYLE,ystyle=YSTYLE    $
     ,  POSITION=POSITION  ,/device    $
     ,  TITLE=TITLE,XTITLE=XTITLE,YTITLE=YTITLE    $
     ,  SUBTITLE=SUBTITLE      $
     ,  TICKLEN=TICKLEN,CHARSIZE=CHARSIZE   $
     ,  XTICKS=XTICKS,XMINOR=XMINOR , YTICKS=YTICKS,YMINOR=YMINOR $
	  ,  XCHARSIZE=XCHARSIZE, YCHARSIZE=YCHARSIZE, XTICKFORMAT=XTICKFORMAT, YTICKFORMAT=YTICKFORMAT $
	  ,  XLOG=XLOG, YLOG=YLOG, XTICKV=XTICKV
;
return
;
errout: print,' TVFRAME : data must be 2-dimensional !'
return
end

pro abre_ps, nombre, TODO=todo, ENCAPSULATED=encapsulated, COLOR=color, LANDSCAPE=landscape,$
	FONT=font, MITAD=mitad, SMALLFONT=smallfont, CENTER=center, XOFFSET=XOFFSET
	
	if (not keyword_set(XOFFSET)) then begin
		xoffset=0.5
	endif
	if (keyword_set(FONT)) then begin
	 	get_lun,u
	 	a = !p.charsize
	   openw,u,'~/idl/fuente.dat'
	 	printf,u,a
	 	close,u	
		free_lun,u
	 	!p.charsize=1.4
	endif
	if (keyword_set(SMALLFONT)) then begin
	 	get_lun,u
	 	a = !p.charsize
	   openw,u,'~/idl/fuente.dat'
	 	printf,u,a
	 	close,u	
		free_lun,u
	 	!p.charsize=1.0
	endif
	a = pswindow(/cm)
	set_plot,'ps',/copy
	encap = 0
	col = 0
	landsca = 0
	if (keyword_set(ENCAPSULATED)) then encap=1
	if (keyword_set(COLOR)) then col = 1
	if (keyword_set(LANDSCAPE)) then landsca = 1
	if (not keyword_set(TODO)) then begin
	 	  if (keyword_set(MITAD)) then begin
		   	device, filename=nombre, ENCAPSULATED=encap, COLOR=col, $
	 	  			LANDSCAPE=landsca, bits_per_pixel=8,xoffset=XOFFSET,yoffset=14.5,$
	 	  			xsize=20.5-XOFFSET,ysize=14	 	  			
		  endif else if (keyword_set(CENTER)) then begin		  
		   	device, filename=nombre, ENCAPSULATED=encap, COLOR=col, $
	 	  			LANDSCAPE=landsca, bits_per_pixel=8, $
					XSIZE=a.xsize, YSIZE=a.ysize, XOFFSET=a.xoffset, $
    	   		 YOFFSET=a.yoffset
		  endif else begin
		   	device, filename=nombre, ENCAPSULATED=encap, COLOR=col, $
	 	  			LANDSCAPE=landsca, bits_per_pixel=8
		  endelse
	 endif else begin
		device, filename=nombre,xoffset=XOFFSET,yoffset=1.5,xsize=19.5,ysize=26, ENCAPSULATED=encap,$
				COLOR=col, LANDSCAPE=landsca, bits_per_pixel=8
	 endelse
end

function histog, a, nbins=nbins, binsize=binsize, plot=plot, optimbin=optimbin, min=min, max=max, _extra=_extra
	if (not keyword_set(min)) then begin
		mina = min(a)
	endif else begin
		mina = min
	endelse
	if (not keyword_set(max)) then begin
		maxa = max(a)
	endif else begin
		maxa = max
	endelse
	if (keyword_set(optimbin)) then begin
		binopt = optbin(a)
		h = histogram(a,binsize=binopt,min=mina,max=maxa)
		n = n_elements(h)
		x = findgen(n) * binopt + mina
		if (keyword_set(plot)) then begin
	 		plot,x,h, psym=10, _extra=_extra
	 	endif
		return,[[x],[h]]
	endif
		
	if (keyword_set(nbins)) then begin
		h = histogram(a,nbins=nbins,min=mina,max=maxa)	 
	 	x = findgen(nbins)/(nbins-1.d0) * (maxa-mina) + mina
	 	if (keyword_set(plot)) then begin
	 		plot,x,h, psym=10, _extra=_extra
	 	endif
	 	return,[[x],[h]]
	endif
	 
	if (keyword_set(binsize)) then begin
		h = histogram(a,binsize=binsize,min=mina,max=maxa)
		n = n_elements(h)
		x = findgen(n) * binsize + mina
		if (keyword_set(plot)) then begin
	 		plot,x,h, psym=10, _extra=_extra
	 	endif
		return,[[x],[h]]
	endif
	 	 
end



pro cierra_ps, FONT=font

; Reset the postcript window	

	device,/inches,xoffset=3./4,yoffset=5,xsize=7,ysize=5
	device,/close
	set_plot,'x'
	if (keyword_set(FONT)) then begin
	 	get_lun,u
	    openr,u,'~/idl/fuente.dat'
	 	readf,u,a
	 	close,u
		free_lun,u
	 	!p.charsize=a
	endif
end

pro posterior_2d, file

; Read the Markov chains

	nlength = file_lines(file+'post_equal_weights.dat')

	one_sigma = (1.d0-2.d0*(1.d0-gauss_pdf(1.)))*100.
	two_sigma = (1.d0-2.d0*(1.d0-gauss_pdf(2.)))*100.

	samples = fltarr(10,nlength)
	openr,2,file+'post_equal_weights.dat'
   readf,2,samples
   close,2
	
	params = ['!7r!3','Y','N','q','!7s!3!dv!n','i','shift','A_v','z']
	params_names = ['sigma','Y','N','q','tauv','i','shift','extinction','redshift']

	!p.multi = [0,7,7]

; Do the plots
	for i = 0, 6 do begin
		for j = 0, 6 do begin

			if (i eq j) then begin
				h = histog(samples[i,*], nbin=40)
				plot, h[*,0], h[*,1] / max(h[*,1]), psym=10
			endif else begin
			
				chainx = reform(samples[i,*])
				chainy = reform(samples[j,*])

				deltax = max(chainx)-min(chainx)
				deltay = max(chainy)-min(chainy)
					
				h = hist_2d(chainx,chainy,bin1=deltax/25,bin2=deltay/25,$
					min1=min(chainx),max1=max(chainx),min2=min(chainy),max2=max(chainy))

				nx = n_elements(h[*,0])
				ny = n_elements(h[0,*])
					
				x = findgen(nx)/(nx-1.d0)*deltax + min(chainx)
				y = findgen(ny)/(ny-1.d0)*deltay + min(chainy)

				h_smooth = filter_image(h,/iterate,smooth=5)

				tvframe, h_smooth/max(h_smooth)*100., xran=[min(x),max(x)], xsty=1, yran=[min(y),max(y)], ysty=1,$
					xtit=params[i],ytit=params[j], charsize=1.1, $
					btit='p('+strtrim(string(params[i]),2)+','+$
					strtrim(string(params[j]),2)+')/p!dmax!n [%]'
				contour, h_smooth/max(h_smooth)*100., x, y, levels=[100.d0-two_sigma,100.d0-one_sigma],/overpl,$
					c_col=[255,255],c_thick=[3,3],c_line=[2,0],c_annotation=['95%','68%'],c_labels=[1,1]
			endelse
		endfor		
   endfor

   !p.multi=0

   stop
end

