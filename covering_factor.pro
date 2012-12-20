@bayesclumpy_routines
; This example shows how to estimate the posterior for the covering factor from the samples of the posterior
; Just pass the root of the file used in Bayesclumpy and the distance to the source in Mpc
; You have to select which method you used for interpolation (LINEAR/NEURAL)
pro covering_factor, file, dist, neural=neural, linear=linear

	if (n_elements(neural)+n_elements(linear) eq 0) then begin
		print, 'You have to select which interpolation method you used (/linear or /neural in the call to read_sed_samples)'
		stop
	endif

	if (n_elements(linear) ne 0) then begin
		restore,'database.idl'
		lambda = database.lambda
		nlam = database.nlam
	endif

	if (n_elements(neural) ne 0) then begin
		restore,'neural.idl'
		lambda = neural.lambda
		nlam = n_elements(lambda)
	endif
	
; Define the parameters
	params = ['sigma','Y','N','q','tauv','angle','shift','extinction','redshift']
	nparam = n_elements(params)

; Read the Markov chains	
	openr,2,file+'post_equal_weights.dat',error=err
	nlength = file_lines(file+'post_equal_weights.dat')
	chain = fltarr(nparam+1,nlength)
	readf,2,chain
	logposterior = reform(chain[nparam,*])
	chain = chain[0:nparam-1,*]
	close,2

; Calculate AGN luminosity
	shif = chain[6,*]
	t1 = 1.0d0*10.d0^(-10.d0+shif)
	t2 = t1*4.*!pi*dist^2 ; dd in Mpc
	Lagn = t2*10^(0.98)*1.d48

	Fbol = fltarr(nlength)

	loop = 0

	seds = fltarr(nlength,nlam)

; Read the SEDs
	n = 0L
	openr,2,file+'.SED_samples',/f77
	readu,2,n
	temp = fltarr(nlam)
	for j = 0L, nlength-1 do begin
		readu,2,temp
		seds[j,*] = temp
	endfor
	close,2

	for i = 0, nlength-1 do begin
		SED_noAGN = seds[i,*]
		Fbol[i] = tsum(lambda*1.d-4,sed_noagn*3.d10/(lambda*1.d-4)^2) * 1.d-26
	endfor

	d = dist * 3.08506d18 * 1.d6
	lbol = Fbol * 4 * !dpi * d^2

	r = Lbol / Lagn

	h = histog(r, nbin=40)
	plot, h[*,0], h[*,1] / max(h[*,1]), xtit='Covering factor', ytit='Normalized posterior', psym=10
	perc = percentile(r, [50.d0 - 50.d0*erf(1.d0/sqrt(2.d0)), 50.d0, 50.d0 + 50.d0*erf(1.d0/sqrt(2.d0))])
	verx, perc[0], line=1
	verx, perc[1], line=0
	verx, perc[2], line=1
	print, 'Apparent covering factor'
	print, 'Distance [Mpc] = ', dist
	print, 'Median = ', perc[1]
	print, '-1sigma = ', perc[1]-perc[0]
	print, '+1sigma = ', perc[2]-perc[1]

end


pro test_covering
	covering_factor, 'MARKOVCHAINS/circinus', 4.d0, /linear
end