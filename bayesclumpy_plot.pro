;-----------------------------------------
; Return an extiction law
; lambda in micron
;-----------------------------------------
function extinction, lambda, model, Av

	k_lambda = fltarr(n_elements(lambda))
	case (model) of

; Seaton (1979) MW
		2 : begin
				Rv = 3.1d0
				invlambda0 = 4.595
				gamm = 1.051
				C1 = -0.38
				C2 = 0.74
				C3 = 3.96
				C4 = 0.26
				invlambda = 1.d0 / lambda

				k_lambda = C1 + C2*invlambda + C3 / ((invlambda - invlambda0^2/invlambda)^2 + gamm^2) + $
					C4*(0.539*(invlambda-5.9)^2 + 0.0564*(invlambda-5.9)^3)

				ind = where(invlambda lt 5.9)
				k_lambda[ind] = C1 + C2*invlambda[ind] + $
					C3 / ((invlambda[ind] - invlambda0^2/invlambda[ind])^2 + gamm^2)

				A_lambda = k_lambda * Av

				return, 10.d0^(-0.4*A_lambda)

			 end

; Fitzpatrick (1986) LMC
		3 : begin
				Rv = 3.1d0
				invlambda0 = 4.608
				gamm = 0.994
				C1 = -0.69
				C2 = 0.89
				C3 = 2.55
				C4 = 0.50
				invlambda = 1.d0 / lambda

				k_lambda = C1 + C2*invlambda + C3 / ((invlambda - invlambda0^2/invlambda)^2 + gamm^2) + $
					C4*(0.539*(invlambda-5.9)^2 + 0.0564*(invlambda-5.9)^3)

				ind = where(invlambda lt 5.9)
				k_lambda[ind] = C1 + C2*invlambda[ind] + $
					C3 / ((invlambda[ind] - invlambda0^2/invlambda[ind])^2 + gamm^2)

				A_lambda = k_lambda * Av

				return, 10.d0^(-0.4*A_lambda)

			 end
				
; Calzetti et al. (2000)
		5 : begin
				Rv = 4.05d0

; From 0.12 micron to 0.63 micron
				ind = where(lambda ge 0.12 and lambda lt 0.63)
				if (ind[0] ne -1) then begin
					k_lambda[ind] = 2.659*(-2.156 + 1.509/lambda[ind] - 0.198/lambda[ind]^2 + $
						0.011/lambda[ind]^3) + Rv
				endif

; From 0.63 micron to 2.20 micron
				ind = where(lambda ge 0.63 and lambda lt 2.20)
				if (ind[0] ne -1) then begin
					k_lambda[ind] = 2.659*(-1.857 + 1.040/lambda[ind]) + Rv
				endif

; Below 0.63 micron
				v1 = 2.659*(-2.156 + 1.509/0.12 - 0.198/0.12^2 + 0.011/0.12^3) + Rv
				v2 = 2.659*(-2.156 + 1.509/0.11 - 0.198/0.11^2 + 0.011/0.11^3) + Rv
				slope = (v1-v2) / (0.12-0.11)
				zero = v1 - slope * 0.12
				ind = where(lambda lt 0.12)
				if (ind[0] ne -1) then begin
					k_lambda[ind] = slope * lambda[ind] + zero
				endif

; Above 2.2 micron
				v1 = 2.659*(-1.857 + 1.040/2.19) + Rv
				v2 = 2.659*(-1.857 + 1.040/2.20) + Rv
				slope = (v1-v2) / (2.19-2.20)
				zero = v1 - slope * 2.19
				ind = where(lambda gt 2.20)
				if (ind[0] ne -1) then begin
					k_lambda[ind] = slope * lambda[ind] + zero
				endif
				ind = where(k_lambda lt 0.0)
				if (ind[0] ne -1) then begin
					k_lambda[ind] = 0.0
				endif
				
				A_lambda = k_lambda / Rv * Av
				
				return, 10.d0^(-0.4*A_lambda)
			end
; Chiar & Tielens (2000)
		6 : begin
				openr,2,'FILTERS/extinction_chiar_tielens2006.dat'
				temp = ''
				for i = 0, 13 do readf,2,temp
				temp = fltarr(3,258)
				readf,2,temp
				close,2
				lambda_pixie = reform(temp[0,*])
				extinction_pixie = reform(temp[2,*])
				ext = interpol(extinction_pixie, lambda_pixie, lambda)
				A_lambda = ext * 0.09d0 * Av
				ind = where(lambda gt 30.d0, count)
				if (count ne 0) then A_lambda[ind] = 0.d0
				return, 10.d0^(-0.4*A_lambda)
			 end
	endcase
end


;-----------------------------------------
; Plot chains
;-----------------------------------------
pro plot_observed_SED, state

	cwpal	

	plot, (*state.obs_x), abs((*state.obs_y)), psym=8, /ylog, /xlog, xran=[0.1,100], xsty=1,$
		tit=state.agn_name, xtit='Wavelength [!7l!6m]',ytit='Flux',yran=[min((*state.obs_y))*0.1,max((*state.obs_y))*10000]
; 	oplot, (*state.obs_x), abs((*state.obs_y)), psym=8, col=5

; Points with standard gaussian errors
	ind = where((*state.obs_y) gt 0)
	if (ind[0] ne -1) then begin
		errplot, (*state.obs_x)[ind], (*state.obs_y)[ind]-(*state.obs_sigma)[ind], $
			(*state.obs_y)[ind]+(*state.obs_sigma)[ind];, col=5
	endif

; Points with upper limits
	ind = where((*state.obs_y) lt 0)
	if (ind[0] ne -1) then begin
		for i = 0, n_elements(ind)-1 do begin
			arrow, (*state.obs_x)[ind], abs((*state.obs_y)[ind]), $
				(*state.obs_x)[ind], (*state.obs_y)[ind]-0.5*abs((*state.obs_y)[ind]), /data;, col=5
		endfor
	endif	
	
end

;-----------------------------------------
; Plot chains
;-----------------------------------------
pro plot_chains, file, param_names, postcript=postscript, type=type, agn_name=agn_name

	consistent = verify_consistency(type)
	if (consistent eq 0) then begin
		res = dialog_message('Posterior samples not compatible with sampler.'+string(10B)+$
				'Change sampler of run again the inference',/error)		
		return
	endif
	
; Standard Markov chains
	if (type eq 0) then begin
; Read the Markov chains
		openr,2,file+'.chain',/f77,error=err
		if (err ne 0) then begin
			res = dialog_message('Error opening posterior samples.'+string(10B)+$
				'You should run again the inference',/error)			
			return
		endif else begin
			nparam = 0L
			nlength = 0L
			readu,2,nparam,nlength

			chain = fltarr(nparam,nlength)
			readu,2,chain
			logposterior = fltarr(nlength)
			acceptance = fltarr(nlength)
			readu,2,logposterior
			readu,2,acceptance
			close,2
		endelse
	endif

; Multinest
	if (type eq 1) then begin
; Read the Markov chains
		nparam = n_elements(param_names)
		openr,2,file+'post_equal_weights.dat',error=err
		if (err ne 0) then begin
			res = dialog_message('Error opening posterior samples.'+string(10B)+$
				'You should run again the inference',/error)
			return
		endif else begin
			nlength = file_lines(file+'post_equal_weights.dat')
			chain = fltarr(nparam+1,nlength)		
			readf,2,chain
			logposterior = reform(chain[nparam,*])
			chain = chain[0:nparam-1,*]
			close,2
		endelse
	endif	

; Do some plots
	!p.multi = [0,3,3]
	for i = 0, 8 do begin
		plot, chain[i,*], psym=3, ytit=param_names[i],charsize=1.4
	endfor
	!p.multi=0

; Do the same plots in POSTCRIPT files
	!p.multi = [0,3,3]	
	open_ps,'PLOTS/chains_'+strtrim(clean(agn_name))+'.ps'
	for i = 0, 8 do begin
		plot, chain[i,*], psym=3, ytit=param_names[i],charsize=1.4
	endfor
	!p.multi=0
	close_ps
	print, 'PLOTS/chains_'+strtrim(clean(agn_name))+'.ps created'

end

;-----------------------------------------
; Plot marginalized
;-----------------------------------------
pro plot_marginalized, file, param_names, state, postcript=postscript, agn_name=agn_name

	consistent = verify_consistency(1)
	if (consistent eq 0) then begin
		res = dialog_message('Posterior samples not compatible with sampler.'+string(10B)+$
				'Change sampler of run again the inference',/error)
		return
	endif
	
; Standard Markov chains
; 	if (state.mcmc_method eq 0) then begin
; 		est = fltarr(9)
; 		errup = fltarr(9)
; 		errdown = fltarr(9)
; 		temp = fltarr(5,9)
; 		
; 	; Read the estimated parameters
; 		openr,2,file+'.confidence',error=err
; 		if (err ne 0) then begin
; 			res = dialog_message('Error opening summary file.'+string(10B)+$
; 				'You should run again the inference',/error)
; 			return
; 		endif else begin
; 			readf,2,temp
; 			close,2
; 			est = temp[0,*]
; 			errdown = temp[1,*]
; 			errup = temp[2,*]
; 		endelse
; 		
; 	; Read the marginalized posteriors
; 		openr,2,file+'.hist1D',/f77,error=err
; 		if (err ne 0) then begin
; 			res = dialog_message('Error opening marginal posteriors.'+string(10B)+$
; 				'You should run again the inference',/error)
; 			return
; 		endif else begin
; 
; ; Do the plots
; 			nparam = 0L   
; 			readu,2,nparam
; 			loop = 0L
; 			n = 0L
; 			!p.multi = [0,3,3]
; 			for i = 0, 8 do begin
; 				readu,2,loop,n		
; 				x = fltarr(n)
; 				yStep = fltarr(n)
; 				yGauss = fltarr(n)
; 				readu,2,x,yGauss, yStep
; 				plot, x, yStep, psym=10, xtit=param_names[i],charsize=1.4, $
; 					xran=[state.ranges_from[i],state.ranges_to[i]],xsty=1
; 				verx, est[i]
; 				verx, est[i]-errdown[i], line=2
; 				verx, est[i]+errup[i], line=2
; 			endfor	
; 			close,2
; 			!p.multi = 0
; 
; ; Do the plots in POSTCRIPT
; 			open_ps,'PLOTS/marginal_'+strtrim(agn_name)+'.ps'
; 			nparam = 0L
; 			openr,2,file+'.hist1D',/f77,error=err
; 			readu,2,nparam
; 			loop = 0L
; 			n = 0L
; 			!p.multi = [0,3,3]
; 			for i = 0, 8 do begin
; 				readu,2,loop,n
; 				x = fltarr(n)
; 				yStep = fltarr(n)
; 				yGauss = fltarr(n)
; 				readu,2,x,yGauss, yStep
; 				plot, x, yStep, psym=10, xtit=param_names[i],charsize=1.4, $
; 					xran=[state.ranges_from[i],state.ranges_to[i]],xsty=1
; 				verx, est[i]
; 				verx, est[i]-errdown[i], line=2
; 				verx, est[i]+errup[i], line=2
; 			endfor
; 			close,2
; 			!p.multi = 0
; 			close_ps
; 		endelse
; 	endif

; Multinest
; 	if (state.mcmc_method eq 1) then begin
		est = fltarr(9)
		errup = fltarr(9)
		errdown = fltarr(9)
		temp = fltarr(6,9)
		
	; Read the estimated parameters
		openr,2,file+'.confidence',error=err
		if (err ne 0) then begin
			res = dialog_message('Error opening summary file.'+string(10B)+$
				'You should run again the inference',/error)
			return
		endif else begin
			readf,2,temp
			close,2
			est = temp[0,*]
			map = temp[1,*]
			errdown = temp[2,*]
			errup = temp[3,*]
		endelse
		
; Read the marginalized posteriors
		openr,2,file+'.hist1D',error=err
		if (err ne 0) then begin
			res = dialog_message('Error opening marginal posteriors.'+string(10B)+$
				'You should run again the inference',/error)
			return
		endif else begin
		
; Do the plots
			nparam = 0L
			readf,2,nparam
			loop = 0L
			n = 0L
			!p.multi = [0,3,3]
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
				plot, x, yStep, psym=10, xtit=param_names[i],charsize=1.4;, $
	; 				xran=[state.ranges_from[i],state.ranges_to[i]],xsty=1
				verx, est[i]
				verx, map[i], line=1
				verx, est[i]-errdown[i], line=2
				verx, est[i]+errup[i], line=2
			endfor	
			close,2
			!p.multi = 0

; Do the plots in POSTCRIPT
			open_ps,'PLOTS/marginal_'+strtrim(clean(agn_name))+'.ps'
			nparam = 0L
			openr,2,file+'.hist1D',error=err
			readf,2,nparam
			loop = 0L
			n = 0L
			!p.multi = [0,3,3]
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
				plot, x, yStep, psym=10, xtit=param_names[i],charsize=1.4;, $
	; 				xran=[state.ranges_from[i],state.ranges_to[i]],xsty=1
				verx, est[i]
				verx, map[i], line=1
				verx, est[i]-errdown[i], line=2
				verx, est[i]+errup[i], line=2
			endfor	
			close,2
			!p.multi = 0
			close_ps
			print, 'PLOTS/marginal_'+strtrim(clean(agn_name))+'.ps created'
		endelse
		
; 	endif
	
end

;-----------------------------------------
; Return MAP values
;-----------------------------------------
pro map_values, file, est, errup, errdown, best, type=type, error=err

	consistent = verify_consistency(type)
	if (consistent eq 0) then begin
		res = dialog_message('Posterior samples not compatible with sampler.'+string(10B)+$
				'Change sampler of run again the inference',/error)
		err = 1
		return
	endif
	
	if (type eq 0) then begin

		est = fltarr(9)
		errup = fltarr(9)
		errdown = fltarr(9)
		best = fltarr(9)
		temp = fltarr(5,9)
		
	; Read the estimated parameters
		openr,2,file+'.confidence',error=err
		if (err ne 0) then begin
			res = dialog_message('Error opening summary file.'+string(10B)+$
				'You should run again the inference',/error)
			return
		endif else begin
			readf,2,temp
			readf,2,best
			close,2
			est = temp[0,*]
			errdown = temp[1,*]
			errup = temp[2,*]
		endelse
	endif

	if (type eq 1) then begin
		est = fltarr(9)
		errup = fltarr(9)
		errdown = fltarr(9)
		temp = fltarr(6,9)
		
	; Read the estimated parameters
		openr,2,file+'.confidence',error=err
		if (err ne 0) then begin
			res = dialog_message('Error opening summary file.'+string(10B)+$
				'You should run again the inference',/error)
			return
		endif else begin
			readf,2,temp
			close,2
			est = temp[0,*]
			best = temp[1,*]
			errdown = temp[2,*]
			errup = temp[3,*]
		endelse
	endif
end

;-----------------------------------------
; Plot compatible models
;-----------------------------------------
pro plot_models, file, param_names, neural, est, errup, errdown, best, state, print_flux=print_flux, seds
	cwpal

	loop = 0

; Read the Markov chain	
	nparam = n_elements(param_names)
	openr,2,file+'post_equal_weights.dat',error=err
	if (err ne 0) then begin
		res = dialog_message('Error opening posterior samples.'+string(10B)+$
			'You should run again the inference',/error)
		return
	endif else begin
		nlength = file_lines(file+'post_equal_weights.dat')
		chain = fltarr(nparam+1,nlength)
		readf,2,chain
		logposterior = reform(chain[nparam,*])
		chain = chain[0:nparam-1,*]
		close,2
	endelse
	
	seds = fltarr(nlength,n_elements((*state.lambda)))
		
	n = 0L
	openr,2,file+'.SED_samples',/f77
	readu,2,n
	temp = fltarr(n_elements((*state.lambda)))
	for j = 0L, nlength-1 do begin
		readu,2,temp
		seds[j,*] = temp
	endfor
		
; The four last couple of SEDs are the one of the median parameters and the MAP
	readu,2,temp
	SED_median = temp
	readu,2,temp
	SED_median_noextinction = temp
	readu,2,temp
	SED_median_noagn = temp
	readu,2,temp
	SED_median_noagn_noextinction = temp
		
	readu,2,temp
	SED_MAP = temp
	readu,2,temp
	SED_MAP_noextinction = temp
	readu,2,temp
	SED_MAP_noagn = temp
	readu,2,temp
	SED_MAP_noagn_noextinction = temp
		
	close,2
	
; 	if (keyword_set(calcseds)) then begin
; 		seds = fltarr(nlength,n_elements(neural.lambda))
; 
; 		progressBar = Obj_New("SHOWPROGRESS", TITLE='Calculating', MESSAGE='Computing models',$
; 			XSIZE=200)
;    	progressBar->Start
;    
; 		loop = 0
; 		i = 0
; 		while (loop le nlength-1) do begin
; 
; 			percent = loop / (nlength-1.d0) * 100.d0
;       	progressBar->Update, percent
;       
; 			if (state.reddening_law ne 0) then begin
; 				extinction_curve = extinction(neural.lambda, state.reddening_law, chain[7,loop])
; 			endif else begin
; 				extinction_curve = replicate(1.d0,n_elements(neural.lambda))
; 			endelse
; 			
; ; 			delta = (chain[0:6,loop] - est)
; ; 			ind = where(abs(delta) lt errup or abs(delta) lt errdown, count)
; 			
; ; 			if (count eq 1) then begin
; 				SED = neural_SED(neural, chain[0,loop], chain[1,loop], chain[2,loop], chain[3,loop], chain[4,loop], chain[5,loop], $
; 					include_agn=state.agn_plus_sed, /jansky)
; 				seds[i,*] = SED/1.d10 * 10.d0^(chain[6,loop]) * extinction_curve
; 				i = i + 1
; ; 			endif
; 			
; 			loop = loop + 1
; 		endwhile
; 
; 		seds = seds[0:i-1,*]
; 
; 		progressBar->Destroy
;    	Obj_Destroy, progressBar
; 	endif


; 	if (state.reddening_law ne 0) then begin
; 		extinction_curve_MAP = extinction((*state.lambda), state.reddening_law, best[7])
; 	endif else begin
; 		extinction_curve_MAP = replicate(1.d0,n_elements((*state.lambda)))
; 	endelse
; 		
; 	SED_MAP = neural_SED(neural, best[0], best[1], best[2], best[3], $
; 		best[4], best[5], include_agn=state.agn_plus_sed, /jansky) / 1.d10 * 10.d0^best[6]

; Plot the MAP SED
	if (state.reddening_law ne 0) then begin
		oplot, (*state.lambda)*(1.d0+best[8]), SED_MAP_noextinction, thick=3, col=4
	endif
	oplot, (*state.lambda)*(1.d0+best[8]), SED_MAP, thick=3, col=2

; Compute the bolometric fluxes
	SED_Fbol = tsum((*state.lambda)*1.d-4,SED_MAP_noextinction*3.d10/((*state.lambda)*1.d-4)^2) * 1.d-26
	if (keyword_set(print_flux)) then print, 'MAP bolometric flux [erg*cm^-2*s^-1] : ', SED_Fbol
	SED_Fbol = tsum((*state.lambda)*1.d-4,SED_MAP*3.d10/((*state.lambda)*1.d-4)^2) * 1.d-26
	if (keyword_set(print_flux)) then print, 'MAP bolometric flux with extinction [erg*cm^-2*s^-1] : ', SED_Fbol

	if (state.agn_plus_sed eq 1) then begin
		SED_Fbol = tsum((*state.lambda)*1.d-4,SED_MAP_noagn_noextinction*3.d10/((*state.lambda)*1.d-4)^2) * 1.d-26
		if (keyword_set(print_flux)) then print, 'MAP bolometric flux [erg*cm^-2*s^-1] (without AGN): ', SED_Fbol
 		SED_Fbol = tsum((*state.lambda)*1.d-4,SED_MAP_noagn*3.d10/((*state.lambda)*1.d-4)^2) * 1.d-26
 		if (keyword_set(print_flux)) then print, 'MAP bolometric flux with extinction [erg*cm^-2*s^-1] (without AGN): ', SED_Fbol
	endif

; Compute the bolometric fluxes
	SED_Fbol = tsum((*state.lambda)*1.d-4,SED_median_noextinction*3.d10/((*state.lambda)*1.d-4)^2) * 1.d-26
	if (keyword_set(print_flux)) then print, 'Median bolometric flux [erg*cm^-2*s^-1] : ', SED_Fbol
	SED_Fbol = tsum((*state.lambda)*1.d-4,SED_median*3.d10/((*state.lambda)*1.d-4)^2) * 1.d-26
	if (keyword_set(print_flux)) then print, 'Median bolometric flux with extinction [erg*cm^-2*s^-1] : ', SED_Fbol

	if (state.agn_plus_sed eq 1) then begin
		SED_Fbol = tsum((*state.lambda)*1.d-4,SED_median_noagn_noextinction*3.d10/((*state.lambda)*1.d-4)^2) * 1.d-26
		if (keyword_set(print_flux)) then print, 'Median bolometric flux [erg*cm^-2*s^-1] (without AGN): ', SED_Fbol
 		SED_Fbol = tsum((*state.lambda)*1.d-4,SED_median_noagn*3.d10/((*state.lambda)*1.d-4)^2) * 1.d-26
 		if (keyword_set(print_flux)) then print, 'Median bolometric flux with extinction [erg*cm^-2*s^-1] (without AGN): ', SED_Fbol
	endif
	
	oplot, (*state.lambda)*(1.d0+est[8]), SED_median, thick=3, col=3

	top = fltarr(n_elements((*state.lambda)))
	bot = fltarr(n_elements((*state.lambda)))
	
	for i = 0, n_elements((*state.lambda))-1 do begin
 		top[i] = percentile(seds[*,i], 50+95.d0/2.d0)
 		bot[i] = percentile(seds[*,i], 50-95.d0/2.d0)
	endfor

	oplot, (*state.lambda)*(1.d0+est[8]), top, line=2, thick=3
	oplot, (*state.lambda)*(1.d0+est[8]), bot, line=2, thick=3

	if (keyword_set(calcseds)) then begin
		openw,2,file+'.synthesis',width=132
		printf,2,'   Wavelength        MAP           Median       Lower_bound     Upper_bound'
		for i = 0, n_elements((*state.lambda))-1 do begin
			printf,2,(*state.lambda)[i], SED_MAP[i], SED_median[i], bot[i], top[i]
		endfor
		close,2
	endif
end

;-----------------------------------------
; Calculates the ratio Ltorus/Lagn
;-----------------------------------------
pro Ltorus_Lagn, file, param_names, dist, postcript=postscript, agn_name=agn_name
	
; Read the Markov chains
	nparam = n_elements(param_names)
	openr,2,file+'post_equal_weights.dat',error=err
	if (err ne 0) then begin
		res = dialog_message('Error opening posterior samples.'+string(10B)+$
			'You should run again the inference',/error)
		return
	endif else begin
		nlength = file_lines(file+'post_equal_weights.dat')
		chain = fltarr(nparam+1,nlength)		
		readf,2,chain
		logposterior = reform(chain[nparam,*])
		chain = chain[0:nparam-1,*]
		close,2
	endelse

; Calculate AGN luminosity
	shif = chain[6,*]
	t1 = 1.0d0*10.d0^(-10.d0+shif)
	t2 = t1*4.*!pi*dist^2 ; dd in Mpc
	Lagn = t2*10^(0.98)*1.d48
	
; Read information about the neural networks
	restore,'neural.idl'
	
	Fbol = fltarr(nlength)
	
	progressBar = Obj_New("SHOWPROGRESS", TITLE='Calculating', MESSAGE='Computing covering factor',$
		XSIZE=200)
   progressBar->Start
   
	for i = 0, nlength-1 do begin
		percent = i / (nlength-1.d0) * 100.d0
      progressBar->Update, percent
      
		SED_noAGN = neural_SED(neural, chain[0,i], chain[1,i], chain[2,i], chain[3,i], $
			chain[4,i], chain[5,i], include_agn=0, /jansky) / 1.d10 * 10.d0^chain[6,i]
	
		Fbol[i] = tsum(neural.lambda*1.d-4,sed_noagn*3.d10/(neural.lambda*1.d-4)^2) * 1.d-26
	endfor
	progressBar->Destroy
   Obj_Destroy, progressBar
   
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

;-----------------------------------------
; Returns the SED luminosity for a set of parameters
;-----------------------------------------
pro sed_luminosity, sigma, Y, N, q, tauv, angle, shif, extinction, reddening_law
	restore,'neural.idl'
	
	SED_MAP_noAGN = neural_SED(neural, sigma, Y, N, q, tauv, angle, include_agn=0, /jansky) / $
		1.d10 * 10.d0^shif

	if (reddening_law ne 0) then begin
		extinction_curve_MAP = extinction(neural.lambda, reddening_law, extinction)
	endif else begin
		extinction_curve_MAP = replicate(1.d0,n_elements(neural.lambda))
	endelse
	
	SED_Fbol = tsum(neural.lambda*1.d-4,SED_MAP_noAGN*3.d10/(neural.lambda*1.d-4)^2) * 1.d-26
		print, 'MAP bolometric flux [erg*cm^-2*s^-1] (without AGN): ', SED_Fbol
	SED_Fbol = tsum(neural.lambda*1.d-4,SED_MAP_noAGN*extinction_curve_MAP*3.d10/(neural.lambda*1.d-4)^2) * 1.d-26
		print, 'MAP bolometric flux with extinction [erg*cm^-2*s^-1] (without AGN): ', SED_Fbol
end

;-----------------------------------------
; Plot a given model in the forward modeling tab
;-----------------------------------------
pro plot_model_database, state
	restore,'neural.idl'
	cwpal
	
	if (state.sed_or_derivative eq 1) then begin
		res = neural_SED_derivative(neural, state.synth[0], state.synth[1], state.synth[2], $
			state.synth[3], state.synth[4], state.synth[5], state.derivative_which_par,$
			include_agn=state.agn_plus_sed, /jansky,out_agn=out_agn) / 1.d10 * 10.d0^state.synth[6]
		tity = 'Derivative flux'
		logy = 0
		rany = [min(res),max(res)]
	endif else begin
		res = neural_SED(neural, state.synth[0], state.synth[1], $
			state.synth[2], state.synth[3], state.synth[4], state.synth[5], include_agn=state.agn_plus_sed, /jansky,$
			out_agn=out_agn) / 1.d10 * 10.d0^state.synth[6]
		tity = 'Flux [mJy]'
		logy = 1
		rany = [1.d-12,1.d2]
	endelse
	

; Introduce the effect of extinction
	if (state.reddening_law ne 0) then begin
		extinction_curve = extinction(neural.lambda, state.reddening_law, state.synth[7])
	endif else begin
		extinction_curve = replicate(1.d0,n_elements(neural.lambda))
	endelse

	redshift = 1.d0 + state.synth[8]

	probAGN = exp(-state.synth[2]*exp(-((90.0-state.synth[5])/state.synth[0])^state.synth[3]))
	widget_control, state.plotSEDWidget, GET_VALUE=which_window
	wset, which_window
	plot, neural.lambda*redshift, res*extinction_curve, xtit='Wavelength [!7l!6m]',ytit=tity,$
		xsty=1, /xlog, ylog=logy, yran=rany;,tit='p(AGN)='+strtrim(string(probAGN),2)
	if (state.agn_plus_sed eq 1) then begin
		oplot, neural.lambda*redshift, out_agn*extinction_curve, line=2
	endif

	if (state.agn_plus_sed eq 1) then begin
		legend,['SED+AGN','AGN'],line=[0,2],/top,/right
	endif else begin
		legend,['SED'],line=[0],/top,/right
	endelse
				
	fil = file_info(state.obsfile)
	if (fil.exists eq 1) then begin
		read_observations, state.obsfile, state
		cwpal
		oplot, (*state.obs_x), (*state.obs_y), col=2, psym=8
		errplot, (*state.obs_x), (*state.obs_y)-(*state.obs_sigma), (*state.obs_y)+(*state.obs_sigma)

		chi2 = total( ((*state.obs_y) - res)^2 / (*state.obs_sigma)^2)
		logL = -0.5 * chi2
; 		xyouts, 0.2,0.25,'log L='+strtrim(string(logL),2),/normal
; 		xyouts, 0.2,0.20,'!7v!6!u2!n='+strtrim(string(chi2),2),/normal
		loadct,0,/silent
	endif
end

;-----------------------------------------
; Plot a given response surface
;-----------------------------------------
pro plot_response_function, state
	restore,'neural.idl'
	cwpal

	nl = n_elements(neural.lambda)

	res = fltarr(nl,100)

; 	loglambda = findgen(nl) / (nl-1.d0) * (alog10(max(neural.lambda))-alog10(min(neural.lambda))) + alog10(min(neural.lambda))
	loglambda = findgen(nl) / (nl-1.d0) * 2.d0

	param = state.synth

	from = state.ranges_from[state.derivative_which_par]
	to = state.ranges_to[state.derivative_which_par]
		
	for i = 0, 99 do begin
		par = from + i / 100.d0 * (to-from)
		param[state.derivative_which_par] = par		
		temp = neural_SED_derivative(neural, param[0], param[1], param[2], $
			param[3], param[4], param[5], state.derivative_which_par,$
			include_agn=state.agn_plus_sed, /jansky,out_agn=out_agn) / 1.d10 * 10.d0^param[6]		
		res[*,i] = interpol(temp, alog10(neural.lambda), loglambda)
	endfor

	loadct,4, /silent
	tvframe, res, xran=[min(loglambda),max(loglambda)], yran=[from,to], /bar, xtit='log !7k!6', $
		ytit='Derivative '+state.params[state.derivative_which_par]
	loadct,0, /silent

end

;-----------------------------------------
; Read the spectrum and the errors from the file
;-----------------------------------------
pro read_observations, file, state

	agn_name = ''
	if (file_test(file)) then begin
		openr,2,file
		readf,2,agn_name

; Spectroscopy?
		if (state.version lt 0) then begin
			readf,2,nlines,nspec
			flux = fltarr(nlines+nspec)			
			error = fltarr(nlines+nspec)
		endif else begin
			readf,2,nlines
			nspec = 0
			flux = fltarr(nlines)
			filt = strarr(nlines)
			error = fltarr(nlines)
		endelse
		state.agn_name = agn_name
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
		if (state.version lt 0 and nspec gt 0) then begin
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
			state.filters = ptr_new(filters)
			if (state.version lt 0 and nspec gt 0) then begin
				x = [filters.info.central,x]
			endif else begin
				x = [filters.info.central]
			endelse
		endif
		
	; Fill the observation arrays
		state.obs_x = ptr_new(x)
		state.obs_y = ptr_new(flux)
		state.obs_sigma = ptr_new(error)
	endif else begin
		res = dialog_message('File with observations does not exist.')
	endelse
	
end