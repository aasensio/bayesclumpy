; Return the perc percetile of an array
function percentile, data, perc
	sidx = sort(data)
	ndata = n_elements(data)
	return, data[sidx[perc*ndata / 100]]
end

; Read SEDs obtained from samples of the MCMC
; It returns the samples in "seds", optionally returning the
; median and the MAP SEDs
; You have to select which method you used for interpolation (LINEAR/NEURAL)
pro read_sed_samples, file, lambda, seds, SED_median=SED_median, SED_noextinction_median=SED_noextinction_median,$
	SED_noagn_median=SED_noagn_median, SED_noagn_noextinction_median=SED_noagn_noextinction_median,$
	SED_MAP=SED_MAP, SED_noextinction_MAP=SED_noextinction_MAP,$
	SED_noagn_MAP=SED_noagn_MAP, SED_noagn_noextinction_MAP=SED_noagn_noextinction_MAP,$
	oneSigmaUp=oneSigmaUp, oneSigmaDown=oneSigmaDown

	restore,'database.idl'
	lambda = database.lambda
	nlam = database.nlam
	
	n = 0L
	openr,2,file+'.SED_samples';,/f77
	readu,2,n
	n = n - 4
	temp = fltarr(nlam)
	seds = dblarr(n,nlam)
	for j = 0L, n-1 do begin		
		readu,2,temp		
		seds[j,*] = temp
	endfor

; The four last couple of SEDs are the one of the median parameters and the MAP
	readu,2,temp
	SED_median = temp
	readu,2,temp
	SED_noextinction_median = temp
	readu,2,temp
	SED_noagn_median = temp
	readu,2,temp
	SED_noagn_noextinction_median = temp

	readu,2,temp
	SED_MAP = temp
	readu,2,temp
	SED_noextinction_MAP = temp
	readu,2,temp
	SED_noagn_MAP = temp
	readu,2,temp
	SED_noagn_noextinction_MAP = temp
	
	oneSigmaUp = fltarr(nlam)
	oneSigmaDown = fltarr(nlam)
	
	for i = 0, nlam-1 do begin
 		oneSigmaUp[i] = percentile(seds[*,i], 50+68.d0/2.d0)
 		oneSigmaDown[i] = percentile(seds[*,i], 50-68.d0/2.d0)
	endfor

	close,2
	
end

; Test program for reading the SEDs and plotting them
pro test_read_sed_samples

	r=[0,255,234,0  ,0  ,255,0  ,255,255,0  ,175,255]
   g=[0,255,0  ,140,216,235,235,0  ,148,126,0  ,67 ]
   b=[0,255,0  ,234,0  ,0  ,228,200,0  ,0  ,201,67 ]
	tvlct, r, g, b
        
	read_sed_samples, 'MARKOVCHAINS/circinus', l, seds, SED_median=SED_median, SED_noextinction_median=SED_noextinction_median,$
		SED_noagn_median=SED_noagn_median, SED_noagn_noextinction_median=SED_noagn_noextinction_median,$
		SED_MAP=SED_MAP, SED_noextinction_MAP=SED_noextinction_MAP,$
		SED_noagn_MAP=SED_noagn_MAP, SED_noagn_noextinction_MAP=SED_noagn_noextinction_MAP, oneSigmaUp=oneSigmaUp, oneSigmaDown=oneSigmaDown

	nseds = n_elements(seds[*,0])
	plot, l, seds[0,*], /nodata, /xlog, /ylog, xran=[0.1,100]
	for i = 0, nseds-1 do begin
		oplot, l, seds[i,*]
	endfor

	oplot, l, SED_median, col=2, thick=3
	oplot, l, SED_MAP, col=3, thick=3
	
	oplot, l, oneSigmaUp, col=4
	oplot, l, oneSigmaDown, col=4

	stop
end