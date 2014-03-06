@bayesclumpy_ann
; Return the perc percetile of an array
function percentile, data, perc
	sidx = sort(data)
	ndata = n_elements(data)
	return, data[sidx[perc*ndata / 100]]
end

; Read the observed SED
pro readObservedSed, file, x, flux, error, agn_name
	agn_name = ''
	if (file_test(file)) then begin
		openr,2,file
		readf,2,agn_name

; Spectroscopy?
		readf,2,nlines,nspec
		flux = fltarr(nlines+nspec)			
		error = fltarr(nlines+nspec)
		
		str = ''
		if (nlines ne 0) then begin
			filt = strarr(nlines)
			for i = 0, nlines-1 do begin
				readf,2,str
				res = strsplit(str,/extract)
				filt[i] = res[0]
				flux[i] = float(res[1])
				error[i] = float(res[2])
			endfor
		endif

; Spectroscopy?
		if (nspec gt 0) then begin
			temp = ''
			readf,2,temp

			x = fltarr(nspec)
			for i = nlines, nlines+nspec-1 do begin
				readf,2,t1,t2,t3
				x[i-nlines] = t1
				flux[i] = t2
				error[i] = t3
			endfor
		endif
		close,2
		
; Read corresponding filters if necessary
		if (nlines ne 0) then begin
			filters = read_filters(filt)			
			if (nspec gt 0) then begin
				x = [filters.info.central,x]
			endif else begin
				x = [filters.info.central]
			endelse
		endif
		
	endif else begin
		res = dialog_message('File with observations does not exist.')
	endelse
end

; Read SEDs obtained from samples of the MCMC
; It returns the samples in "seds", optionally returning the
; median and the MAP SEDs
; You have to select which method you used for interpolation (LINEAR/NEURAL)
pro read_sed_samples, file, lambda, seds, seds_noAGN=seds_noAGN, SED_median=SED_median, SED_noextinction_median=SED_noextinction_median,$
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
	seds_noAGN = dblarr(n,nlam)
	for j = 0L, n-1 do begin
		readu,2,temp		
		seds[j,*] = temp
	endfor
	
	for j = 0L, n-1 do begin
		readu,2,temp		
		seds_noAGN[j,*] = temp
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
        
    readObservedSed, 'OBSERVATIONS/circinus.cat', x, flux, error, agn_name
	read_sed_samples, 'MARKOVCHAINS/circinus', l, seds, SED_median=SED_median, SED_noextinction_median=SED_noextinction_median,$
		SED_noagn_median=SED_noagn_median, SED_noagn_noextinction_median=SED_noagn_noextinction_median,$
		SED_MAP=SED_MAP, SED_noextinction_MAP=SED_noextinction_MAP,$
		SED_noagn_MAP=SED_noagn_MAP, SED_noagn_noextinction_MAP=SED_noagn_noextinction_MAP, oneSigmaUp=oneSigmaUp, oneSigmaDown=oneSigmaDown
		
	plot, x, abs(flux), psym=8, /ylog, /xlog, xran=[0.1,100], xsty=1,$
		tit=agn_name, xtit='Wavelength [!7l!6m]',ytit='Flux',yran=[min(flux)*0.1,max(flux)*10.0]
		
	ind = where(flux gt 0)
	if (ind[0] ne -1) then begin
		errplot, x[ind], flux[ind]-error[ind], flux[ind]+error[ind];, col=5
	endif

; Points with upper limits
	ind = where(flux lt 0)
	if (ind[0] ne -1) then begin
		for i = 0, n_elements(ind)-1 do begin
			arrow, x[ind], abs(flux[ind]), $
				x[ind], flux[ind]-0.5*abs(flux[ind]), /data;, col=5
		endfor
	endif

	nseds = n_elements(seds[*,0])	
	for i = 0, nseds-1 do begin
; 		oplot, l, seds[i,*]
	endfor

	oplot, l, SED_median, col=2, thick=3
	oplot, l, SED_MAP, col=3, thick=3
	
	oplot, l, oneSigmaUp, line=2, thick=3
	oplot, l, oneSigmaDown, line=2, thick=3

	stop
end