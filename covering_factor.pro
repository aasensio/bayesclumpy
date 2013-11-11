@bayesclumpy_routines
@read_samples
@read_seds
; This example shows how to estimate the posterior for the covering factor from the samples of the posterior
; Just pass the root of the file used in Bayesclumpy and the distance to the source in Mpc
; It also computes the total covering factor f2
pro covering_factor, file, dist

	restore,'database.idl'
	lambda = database.lambda
	nlam = database.nlam
	
; Define the parameters
	params = ['sigma','Y','N','q','tauv','angle','shift','extinction','redshift']
	nparam = n_elements(params)

; Read the Markov chains
	read_samples, file, chain
	
	nlength = n_elements(chain[0,*])
	
; Read the SED samples
	read_sed_samples, file, lambda, seds, SED_median=SED_median, SED_noextinction_median=SED_noextinction_median,$
		SED_noagn_median=SED_noagn_median, SED_noagn_noextinction_median=SED_noagn_noextinction_median,$
		SED_MAP=SED_MAP, SED_noextinction_MAP=SED_noextinction_MAP,$
		SED_noagn_MAP=SED_noagn_MAP, SED_noagn_noextinction_MAP=SED_noagn_noextinction_MAP, oneSigmaUp=oneSigmaUp, oneSigmaDown=oneSigmaDown
		
; Calculate AGN luminosity
	shif = chain[6,*]
	t1 = 1.0d0*10.d0^(-10.d0+shif)
	t2 = t1*4.*!pi*dist^2 ; dd in Mpc
	Lagn = t2*10^(0.98)*1.d48

	Fbol = fltarr(nlength)

	loop = 0

	for i = 0, nlength-1 do begin
		SED_noAGN = seds[i,*]
		Fbol[i] = tsum(lambda*1.d-4,sed_noagn*3.d10/(lambda*1.d-4)^2) * 1.d-26
	endfor

	d = dist * 3.08506d18 * 1.d6
	lbol = Fbol * 4 * !dpi * d^2

	r = Lbol / Lagn

	!p.multi = [0,1,2]
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
	
	f2 = fltarr(nlength)
        
	beta = findgen(200) / 199.d0 * 90.d0
	mui = cos(beta* !DPI / 180.d0)
	N0 = chain[2,*]
	sigma = chain[0,*]

	for i = 0, nlength-1 do begin
		mui = cos(beta * !DPI / 180.d0)
		Nt = N0[i] * exp(-beta^2 / sigma[i]^2)
		f2[i] = 1.d0 - int_tabulated(mui, exp(-Nt), /sort)
	endfor
	
	h = histog(f2, nbin=40)
	plot, h[*,0], h[*,1] / max(h[*,1]), xtit='Total covering factor', ytit='Normalized posterior', psym=10
	
	!p.multi = 0
	
	stop

end


pro test_covering
	covering_factor, 'MARKOVCHAINS/circinus', 4.d0
end