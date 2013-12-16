; This example file shows how to read the marginal posteriors from the output of Bayesclumpy
; Just pass the root of the file used in Bayesclumpy

pro read_marginal, file

; Define the parameters
	params = ['sigma','Y','N','q','tauv','angle','shift','extinction','redshift']
	nparam = n_elements(params)

	est = fltarr(9)
	errup = fltarr(9)
	errdown = fltarr(9)
	temp = fltarr(6,9)

; Read the estimated parameters
	openr,2,file+'.confidence',error=err
	readf,2,temp
	close,2
	est = temp[0,*]
	map = temp[1,*]
	errdown = temp[2,*]
	errup = temp[3,*]

; Read the marginalized posteriors
	openr,2,file+'.hist1D',error=err

; Do the plots
	nparam = 0L
	readf,2,nparam
	loop = 0L
	n = 0L
	!p.multi = [0,3,3]
	mode = est*0.d0
	for i = 0, 8 do begin
		readf,2,loop,n
		x = fltarr(n)
		yStep = fltarr(n)
		yGauss = fltarr(n)
		temp = fltarr(3,n)
		readf,2,temp
		x = temp[0,*]
		yGauss = temp[1,*]
		yStep = temp[2,*]
		plot, x, yStep, psym=10, xtit=params[i],charsize=1.4
		verx, est[i]
		verx, map[i], line=1
		verx, est[i]-errdown[i], line=2
		verx, est[i]+errup[i], line=2
		
; Compute the mode of the histogram using a parabolic fit to the maximum to get sub-pixel precision
		maximum = max(yStep, maxLoc)
		left = maxLoc - 2
		if (left lt 0) then left = 0
		right = maxLoc + 2
		if (right gt n_elements(yStep)) then right = n_elements(yStep)-1
		res = poly_fit(x[left:right], yStep[left:right], 2)
		mode[i] = -res[1] / (2.d0*res[2])
	endfor
	close,2
	!p.multi = 0

end

pro test_read_marginal
	read_marginal, 'MARKOVCHAINS/circinus'
end