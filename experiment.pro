; Run the inference for a given galaxy
pro experiment_doinference, file, nobs

	str = strarr(31)
	openr,2,'chain_original.cfg'
	readf,2,str
	close,2

	name = file+'_'+strtrim(string(nobs),2)
	sy_type = 2
	redshift = 0.0014


	str[12] = "'OBSERVATIONS/NEXTFILTER/"+name+".cat'"

	str[14] = "'MARKOVCHAINS/NEXTFILTER/"+name+"'"

; First do inference
	str[2] = '-2'

; If Sy1, then use AGN and extinction
	if (sy_type eq 1) then begin
		str[4] = '9'
		str[16] = '1'
		str[18] = '6'
		str[29] = "'extinction'  0.00000   5.00000  'U'  0.00000  0.00000"
	endif else begin
		str[4] = '9'
		str[16] = '0'
		str[18] = '0'
		str[29] = "'extinction'  0.00000   5.00000  'D'  0.00000  0.00000"
	endelse

	str[30] = "'redshift'  0.00000   6.00000  'D'  "+strtrim(string(redshift),2)+"  0.00000"

	openw,2,'chain.cfg'
	for k = 0, 30 do printf,2,str[k]
	close,2
	spawn,'./clumpy_mcmc'


; If Sy1, then use AGN and extinction
	if (sy_type eq 1) then begin
		str[4] = '9'
		str[16] = '1'
		str[18] = '6'
		str[29] = "'extinction'  0.00000   5.00000  'U'  0.00000  0.00000"
	endif else begin
		str[4] = '9'
		str[16] = '0'
		str[18] = '0'
		str[29] = "'extinction'  0.00000   5.00000  'D'  0.00000  0.00000"
	endelse

	str[30] = "'redshift'  0.00000   6.00000  'D'  "+strtrim(string(redshift),2)+"  0.00000"

; Then suggest new filter
	str[2] = '-3'

	openw,2,'chain.cfg'
	for k = 0, 30 do printf,2,str[k]
	close,2

	spawn,'./clumpy_mcmc'

end

; Read the expected utility
pro read_eu, file, filt_name, lambda, information, xtit
	nfilt = file_lines(file+'_information_gain.dat')
	filt_name = strarr(nfilt)
	lambda = dblarr(nfilt)
	information = dblarr(nfilt)
	openr,2,file+'_information_gain.dat'
	t = ''
	for i = 0, nfilt-1 do begin
		readf,2,t
		res = strsplit(t,' ',/extract)
	   filt_name[i] = res[0]
	   lambda[i] = double(res[1])
	   information[i] = double(res[2])
	endfor
	close,2

	colors = intarr(nfilt)
	colors[where(lambda lt 6.d0)] = 0
	colors[where(lambda ge 6.d0 and lambda lt 25)] = 2
	colors[where(lambda ge 25.d0 and lambda lt 200)] = 3
	colors[where(lambda ge 200.d0)] = 4

	information = information - min(information)
	plot, lambda, information / max(information), /xlog, psym=8, yran=[0.001,1.2], ysty=1, xran=[0.5, 999.99], xsty=1, $
		charsize=1.5, /nodata, xtit=xtit
	for i = 0, nfilt-1 do plots, lambda[i], information[i] / max(information), col=colors[i], psym=8	
	verx, 4.d0, line=1

end

pro readSedSamples, file, lambda, seds

	nlength = file_lines(file+'post_equal_weights.dat')

	seds = fltarr(nlength,n_elements(lambda))

	read_sed_samples, file, lambda, seds
	
; 	n = 0L
; 	openr,2,file+'.SED_samples',/f77
; 	readu,2,n
; 	temp = fltarr(n_elements(lambda))
; 	for j = 0L, nlength-1 do begin
; 		readu,2,temp
; 		seds[j,*] = temp
; 	endfor
; 	close,2

; Select only a subsample of 200 from it
	ind = randperm(nlength)

	seds = seds[ind[0:199],*]

end

pro central_lambda, lambda0, name_filter
	openr,2,'FILTERS/normalizations.dat'
	readf,2,nfilters
	str = ''
	lambda0 = fltarr(nfilters)
	name_filter = strarr(nfilters)
	for i = 0, nfilters-1 do begin
		readf,2,str
		res = strsplit(str,' ',/extract)
		res2 = strsplit(res[0],"'",/extract)

		file = 'FILTERS/'+res2+'.res'
		openr,3,file
		readf,3,ndat,normal
		trans = dblarr(2,ndat)
		readf,3,trans
		close,3

		lambda0[i] = tsum(trans[0,*],trans[1,*]*trans[0,*]) / tsum(trans[0,*],trans[1,*])
		name_filter[i] = res2

	endfor
	close,2
end


pro read_sed, file, lambda, filt, flux, error
	central_lambda, lambda0, name_filter

	agn_name = ''
	str = ''
	openr,2,file
	readf,2,agn_name
	readf,2,nlines,nspec
	flux = dblarr(nlines)
	error = dblarr(nlines)
	lambda = dblarr(nlines)
	filt = strarr(nlines)
	for i = 0, nlines-1 do begin
		readf,2,str
		res = strsplit(str,/extract)
		filt[i] = strsplit(res[0],"'",/extract)
		flux[i] = double(res[1])
		error[i] = double(res[2])
		lambda[i] = lambda0[where(name_filter eq filt[i])]
	endfor
	close,2
end

pro runall
	cwpal

; 	file = 'circinus'
	file = 'ngc3081'
	
	restore,'database.idl'

	read_sed, 'OBSERVATIONS/NEXTFILTER/'+file+'.cat', lambda_sed, filt_available, flux_available, error_available

	for k = 1, n_elements(flux_available) do begin

; Carry out the inference
   	 experiment_doinference, file, k

; Read the outputs
		readSedSamples, 'MARKOVCHAINS/NEXTFILTER/'+file+'_'+strtrim(string(k),2), database.lambda, seds
		read_sed, 'OBSERVATIONS/NEXTFILTER/'+file+'_'+strtrim(string(k),2)+'.cat', lambda_sed, filt, flux, error

; Plot the SEDs and the expected utility
		!p.multi = [0,2,1]
		plot, lambda_sed, flux, psym=8, /ylog, /xlog, xran=[0.7,1000], xsty=1, charsize=1.5, yran=[1.d-2,1.d5], ysty=1,/nodata
		for i = 0, n_elements(seds[*,0])-1 do oplot, database.lambda, seds[i,*], col=150
		oplot, lambda_sed, flux, psym=8

		read_eu, 'MARKOVCHAINS/NEXTFILTER/'+file+'_'+strtrim(string(k),2), filt_name, lambda, information

; Locate the next best available filter

; Try to find the filter with a wavelength difference below 4 micron to the best filter
; If none is available, then go down until finding the available filter with the
; larges expected utility
		maxinf = max(information, loc)
		lambda0 = lambda[maxinf]
		lambda_dif = abs(lambda-lambda0)

; First verify if the filter with the largest EU is available
; 		res = where(filt_available eq filt_name[loc], count)
; 
; 		if (count eq 0) then begin
; 
; 			res = -1
; 		
; ; Try to find the filter with a wavelength difference below 4 micron to the best filter
; 			ind = where(lambda_dif gt 0 and lambda_dif lt 4.d0, counttot)
; 			if (counttot ne 0) then begin
; 				found = 0
; 				i = 0
; 				while (found eq 0 and i lt counttot) do begin
; 					out = where(filt_available eq filt_name[ind[i]], count)
; 					if (count ne 0) then begin
; 						res = out
; 						found = 1
; 					endif
; 					i = i + 1
; 				endwhile
; 			endif
; 
; ; If none is available, then go down until finding the available filter with the
; ; larges expected utility
; 			if (res eq -1) then begin
; 				ind = reverse(sort(information))
; 				loop = 0
; 				while (where(filt_available eq filt_name[ind[loop]]) eq -1) do begin
; 					loop = loop + 1
; 				endwhile
; 				res = where(filt_available eq filt_name[ind[loop]])
; 			endif
; 		endif

		ind = reverse(sort(information))
		loop = 0
		while (where(filt_available eq filt_name[ind[loop]]) eq -1) do begin
			loop = loop + 1
		endwhile
		res = where(filt_available eq filt_name[ind[loop]])

; Generate the new SED if we have not used all observations
		if (k ne n_elements(flux)) then begin
			openw,2,'OBSERVATIONS/NEXTFILTER/'+file+'_'+strtrim(string(k+1),2)+'.cat'
			printf,2,"'"+file+"'"
			printf,2,k+1,0
			for i = 0, k-1 do begin
				printf,2,"'"+filt[i]+"'", flux[i], error[i]
			endfor
			printf,2,FORMAT='(A,F,F)',"'"+filt_available[res]+"'",flux_available[res],error_available[res]
			close,2

			print, 'Adding ', filt_available[res]
		endif

	endfor

	!p.multi = 0

	stop
end

pro analyze_all, postcript=postcript

	!x.thick = 4
	!y.thick = 4
	cwpal
	
	file = 'ngc3081'
; 	file = 'circinus'
	
	restore,'database.idl'
	read_sed, 'OBSERVATIONS/NEXTFILTER/'+file+'.cat', lambda_sed, filt_available, flux_available, error_available

	normaliz = flux_available[cerca(lambda_sed,10.d0)]

	if (keyword_set(postcript)) then begin
		abre_ps,'sampling_evolution_'+file+'.eps',/encaps,/perfect,/color
	endif else begin
		erase
		window, 0, xsize=800, ysize=1000
	endelse

	multiplot2, [6,7], ygap=0.005, xgap=0.005

	nf = n_elements(flux_available)

	for k = 1, nf do begin
	
		readSedSamples, 'MARKOVCHAINS/NEXTFILTER/'+file+'_'+strtrim(string(k),2), database.lambda, seds
		read_sed, 'OBSERVATIONS/NEXTFILTER/'+file+'_'+strtrim(string(k),2)+'.cat', lambda_sed, filt, flux, error

		xtit = ''
		if (k eq nf) then begin
			xtit = textoidl('\lambda')+' ['+textoidl('\mu')+'m]'
		endif
		
; Plot the SEDs and the expected utility		
		plot, lambda_sed, flux / normaliz, psym=8, /ylog, /xlog, xran=[0.7,999.99], xsty=1, charsize=1.5, yran=[1.d-4,99.99], ysty=1,/nodata,xtit=xtit
		for i = 0, n_elements(seds[*,0])-1 do oplot, database.lambda, seds[i,*] / normaliz, col=150
		oplot, lambda_sed, flux / normaliz, psym=8
		errplot, lambda_sed, (flux-error) / normaliz, (flux+error) / normaliz, thick=3

		multiplot2


		read_eu, 'MARKOVCHAINS/NEXTFILTER/'+file+'_'+strtrim(string(k),2), filt_name, lambda, information, xtit

		multiplot2

		openr,2,'MARKOVCHAINS/NEXTFILTER/'+file+'_'+strtrim(string(k),2)+'.hist1D',error=err
		nparam = 0L
		readf,2,nparam
		loop = 0L
		n = 0L		
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

; sigma
			if (i eq 0) then begin
				if (k eq nf) then begin
					plot, x, yStep / max(yStep), psym=10, charsize=1.4, xran=[15,max(x)], xsty=1, xtickv=[15,30,45,60,max(x)], xticks=3,$
						xtit=textoidl('\sigma'), thick=4, yran=[0,1.2]
				endif else begin
					plot, x, yStep / max(yStep), psym=10, charsize=1.4, xran=[15,max(x)], xsty=1, xtickv=[15,30,45,60,max(x)], xticks=3, thick=4, yran=[0,1.2]
				endelse				
				multiplot2
			endif


; Y
			if (i eq 1) then begin
				if (k eq nf) then begin
					plot, x, yStep / max(yStep), psym=10, charsize=1.4, xran=[5,max(x)], xsty=1, xtickv=[20,40,60,80], xticks=3, xtit='Y', $
						thick=4, yran=[0,1.2]
				endif else begin
					plot, x, yStep / max(yStep), psym=10, charsize=1.4, xran=[5,max(x)], xsty=1, xtickv=[20,40,60,80], xticks=3, thick=4, yran=[0,1.2]
				endelse
				multiplot2
			endif
			
; N0
			if (i eq 2) then begin
				if (k eq nf) then begin
					plot, x, yStep / max(yStep), psym=10, charsize=1.4, xran=[1,max(x)], xsty=1, xtickv=[1.001,5,10,15], xticks=3, xtit='N!d0!n', $
						thick=4, yran=[0,1.2]
				endif else begin
					plot, x, yStep / max(yStep), psym=10, charsize=1.4, xran=[1,max(x)], xsty=1, xtickv=[1.001,5,10,15], xticks=3, thick=4, yran=[0,1.2]
				endelse
				multiplot2, /doyaxis
			endif

; i
			if (i eq 5) then begin
				if (k eq nf) then begin
					plot, x, yStep / max(yStep), psym=10, charsize=1.4, xran=[0.001,max(x)], xsty=1, xtit='i', thick=4, yran=[0,1.2], ysty=8,$
						ytickname=replicate(' ',30)
				endif else begin
					plot, x, yStep / max(yStep), psym=10, charsize=1.4, xran=[0.001,max(x)], xsty=1, thick=4, yran=[0,1.2], ysty=8,$
						ytickname=replicate(' ',30)
				endelse
				axis, yaxis=1, yran=[0,1.1999], ysty=1
				multiplot2
			endif
	
		endfor
		close,2
		
	endfor

	if (keyword_set(postcript)) then begin
		cierra_ps,/bounds
	endif
	
	multiplot2, /default
	stop

end