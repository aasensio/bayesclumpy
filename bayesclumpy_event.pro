;-----------------------------------------
; Event handler
;-----------------------------------------
pro bayesclumpy_Event, event
	 widget_control, Event.id, GET_UVALUE=Action
	 widget_control, Event.top, GET_UVALUE=state

	 handler = Event.Top
	 
	 case Action of
	 	'TABS': begin
                end
		'OBS_FILE' : begin
      	widget_control, Event.id, GET_VALUE=value
      	state.obsfile = value
			widget_control, handler, SET_UVALUE=state
			widget_control, state.fileText, SET_VALUE=value

; Read the spectrum and the errors from the file
			read_observations, value, state			
		end
		
		'OBS_FILE_DIALOG' : begin
			value = dialog_pickfile(/read)
			state.obsfile = value
			widget_control, handler, SET_UVALUE=state
			widget_control, state.fileText, SET_VALUE=value

; Read the SED and the errors from the file
			read_observations, value, state
			
		 end

		 'PLOT_SED' : begin		 	
			widget_control, state.plotrightWidget, GET_VALUE=which_window
			wset, which_window
			plot_observed_SED, state
		 end

		 'PLOT_MAP' : begin
			widget_control, state.plotrightWidget, GET_VALUE=which_window
			wset, which_window
			map_values, 'MARKOVCHAINS/'+state.outfile, est, errup, errdown, best, type=1, error=err	
			if (err eq 0) then begin
				restore,'neural.idl'

; Plot the best fits
				wset, which_window
				plot_observed_SED, state
				plot_models, 'MARKOVCHAINS/'+state.outfile, state.params, neural, est, errup, errdown, best, state, /print_flux, seds

; Plot the best fits in POSTCRIPT
				open_ps,'PLOTS/MAP_'+strtrim(clean(state.agn_name))+'.ps', /color
				plot_observed_SED, state				
				plot_models, 'MARKOVCHAINS/'+state.outfile, state.params, neural, est, errup, errdown, best, state, seds
				close_ps
				print, 'Figure PLOTS/MAP_'+strtrim(clean(state.agn_name))+'.ps created'
			endif
		 end

		 'DO_INFERENCE' : begin
		 
; Verify that an observation file is present
			if (not ptr_valid(state.obs_x)) then begin
				res = dialog_message('You have to define the observations first',/error)
				return
			endif

; Verify that, if any Gaussian prior is used, the correct data is entered (sigma and mu for
; each Gaussian)
			ind = where(state.prior eq 'G')			
			if (ind[0] ne -1) then begin
				for i = 0, n_elements(ind)-1 do begin
					if (state.gaussian_prior_info[0,ind[i]] eq '' or $
						state.gaussian_prior_info[1,ind[i]] eq '') then begin
							res = dialog_message('Gaussian prior info for '+state.params[ind[i]]+' not defined',/error)
							return
					endif
				endfor
			endif
			
			ind = where(state.prior eq 'D')
			if (ind[0] ne -1) then begin
				for i = 0, n_elements(ind)-1 do begin
					if (state.gaussian_prior_info[0,ind[i]] lt state.ranges_from[ind[i]] and $
					 state.gaussian_prior_info[0,ind[i]] gt state.ranges_to[ind[i]]) then begin
						res = dialog_message('Dirac prior info for '+state.params[ind[i]]+' not defined',/error)
						return
					endif
				endfor
			endif
			
			ind = where(state.ranges_to - state.ranges_from lt 0)			
			if (ind[0] ne -1) then begin
				res = dialog_message('Limit of parameters {'+state.params[ind]+'} not set correctly',/error)
			endif else begin

				str = strarr(31)
				openr,2,'chain_original.cfg'
				readf,2,str
				close,2

				if (state.avg_method eq 1) then begin
					if (state.neural_interpolation eq 0) then begin
						str[2] = '4'
					endif else begin
						str[2] = '-4'
					endelse
				endif else begin
					if (state.neural_interpolation eq 0) then begin
						str[2] = '2'
					endif else begin
						str[2] = '-2'
					endelse
				endelse
				
				str[4] = '9'
				str[6] = strtrim(string(state.niter_mcmc),2)
				str[8] = strtrim(string(state.burnin),2)				
				str[10] = "'MCMC'"
				str[12] = "'"+state.obsfile+"'"
				str[14] = "'MARKOVCHAINS/"+state.outfile+"'"
				str[16] = strtrim(string(state.agn_plus_sed),2)
				str[18] = strtrim(string(state.reddening_law),2)

				for i = 0, 8 do begin
					str[i+22] = "'"+state.params[i]+"'  "+strtrim(string(state.ranges_from[i]),2)+$
						"   "+strtrim(string(state.ranges_to[i]),2)+"  '"+state.prior[i]+"'  "+$
						strtrim(string(state.gaussian_prior_info[0,i]),2)+"  "+$
						strtrim(string(state.gaussian_prior_info[1,i]),2)
				endfor

				openw,2,'chain.cfg'
				for i = 0, 30 do begin
					printf,2,str[i]
				endfor				
				close,2

; Verify that the code exists and run it
				res = file_test('clumpy_mcmc')
				if (res eq 1) then begin
					spawn,'./clumpy_mcmc'
; Re-run until sampling is correct
					while (verify_finished() eq 0) do begin
						spawn,'./clumpy_mcmc'
					endwhile
				endif else begin
					res = dialog_message('BayesCLUMPY executable is not present.'+string(10B)+$
						'Compile the code and try again',/error)
					return
				endelse
				

				widget_control, state.plotleftWidget, GET_VALUE=which_window
				wset, which_window
				plot_marginalized, 'MARKOVCHAINS/'+state.outfile, state.params, state, $
					agn_name=state.agn_name

				widget_control, state.plotrightWidget, GET_VALUE=which_window
				wset, which_window				
				plot_chains, 'MARKOVCHAINS/'+state.outfile, state.params, type=1, $
					agn_name=state.agn_name
								
			endelse
		 end

		 'SUGGEST_OBSERVATION' : begin
		 
; Verify that an observation file is present
			if (not ptr_valid(state.obs_x)) then begin
				res = dialog_message('You have to define the observations first',/error)
				return
			endif

; Verify that, if any Gaussian prior is used, the correct data is entered (sigma and mu for
; each Gaussian)
			ind = where(state.prior eq 'G')			
			if (ind[0] ne -1) then begin
				for i = 0, n_elements(ind)-1 do begin
					if (state.gaussian_prior_info[0,ind[i]] eq '' or $
						state.gaussian_prior_info[1,ind[i]] eq '') then begin
							res = dialog_message('Gaussian prior info for '+state.params[ind[i]]+' not defined',/error)
							return
					endif
				endfor
			endif
			
			ind = where(state.prior eq 'D')
			if (ind[0] ne -1) then begin
				for i = 0, n_elements(ind)-1 do begin
					if (state.gaussian_prior_info[0,ind[i]] lt state.ranges_from[ind[i]] and $
					 state.gaussian_prior_info[0,ind[i]] gt state.ranges_to[ind[i]]) then begin
						res = dialog_message('Dirac prior info for '+state.params[ind[i]]+' not defined',/error)
						return
					endif
				endfor
			endif
			
			ind = where(state.ranges_to - state.ranges_from lt 0)			
			if (ind[0] ne -1) then begin
				res = dialog_message('Limit of parameters {'+state.params[ind]+'} not set correctly',/error)
			endif else begin

				str = strarr(31)
				openr,2,'chain_original.cfg'
				readf,2,str
				close,2

				if (state.neural_interpolation eq 0) then begin
					str[2] = '3'
				endif else begin
					str[2] = '-3'
				endelse
				
				str[4] = '9'
				str[6] = strtrim(string(state.niter_mcmc),2)
				str[8] = strtrim(string(state.burnin),2)				
				str[10] = "'MCMC'"
				str[12] = "'"+state.obsfile+"'"
				str[14] = "'MARKOVCHAINS/"+state.outfile+"'"
				str[16] = strtrim(string(state.agn_plus_sed),2)
				str[18] = strtrim(string(state.reddening_law),2)

				for i = 0, 8 do begin
					str[i+22] = "'"+state.params[i]+"'  "+strtrim(string(state.ranges_from[i]),2)+$
						"   "+strtrim(string(state.ranges_to[i]),2)+"  '"+state.prior[i]+"'  "+$
						strtrim(string(state.gaussian_prior_info[0,i]),2)+"  "+$
						strtrim(string(state.gaussian_prior_info[1,i]),2)
				endfor

				openw,2,'chain.cfg'
				for i = 0, 30 do begin
					printf,2,str[i]
				endfor				
				close,2

; Verify that the code exists and run it
				res = file_test('clumpy_mcmc')
				if (res eq 1) then begin
					spawn,'./clumpy_mcmc'
				endif else begin
					res = dialog_message('BayesCLUMPY executable is not present.'+string(10B)+$
						'Compile the code and try again',/error)
					return
				endelse
				
				widget_control, state.plotleftWidget, GET_VALUE=which_window_left
				widget_control, state.plotrightWidget, GET_VALUE=which_window_right				
				plot_expected, 'MARKOVCHAINS/'+state.outfile, which_window_left, which_window_right
								
			endelse
		 end

		 'PLOT_INFERENCE' : begin
; 		 	restore,'markov_chain.idl'
		 	
		 	widget_control, state.plotleftWidget, GET_VALUE=which_window
			wset, which_window			
			plot_marginalized, 'MARKOVCHAINS/'+state.outfile, state.params, state, $
					agn_name=state.agn_name

			widget_control, state.plotrightWidget, GET_VALUE=which_window
			wset, which_window
			plot_chains, 'MARKOVCHAINS/'+state.outfile, state.params, type=1, $
					agn_name=state.agn_name
		 end
		 
		 'CALCULATE_COVFACTOR' : begin
			widget_control, state.plotrightWidget, GET_VALUE=which_window
			wset, which_window
			Ltorus_Lagn, state, 'MARKOVCHAINS/'+state.outfile, state.params, state.distance, agn_name=state.agn_name
		 end

		 'SUGGEST_OBSERVATION' : begin
		 	 widget_control, state.plotleftWidget, GET_VALUE=which_window
			 wset, which_window
		 	 experiment_design, 'MARKOVCHAINS/'+state.outfile, state, type=1
		 end

		 'OUTFILE' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.outfile = value
		 end
		 'DISTANCE' : begin
			widget_control, Event.id, GET_VALUE=value
		 	state.distance = value
		 end
		 'BURNIN_MCMC' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.burnin = float(value)
		 end

		 'sigmaminSlider' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.ranges_from[0] = value		 	
		 end
		 'sigmamaxSlider' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.ranges_to[0] = value
		 end
		 'YminSlider' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.ranges_from[1] = value
		 end
		 'YmaxSlider' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.ranges_to[1] = value
		 end
		 'NminSlider' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.ranges_from[2] = value
		 end
		 'NmaxSlider' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.ranges_to[2] = value
		 end
		 'qminSlider' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.ranges_from[3] = value
		 end
		 'qmaxSlider' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.ranges_to[3] = value
		 end
		 'tauvminSlider' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.ranges_from[4] = value
		 end
		 'tauvmaxSlider' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.ranges_to[4] = value
		 end
		 'angleminSlider' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.ranges_from[5] = value
		 end
		 'anglemaxSlider' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.ranges_to[5] = value
		 end
		 'shiftminSlider' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.ranges_from[6] = value
		 end
		 'shiftmaxSlider' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.ranges_to[6] = value
		 end
		'extminSlider' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.ranges_from[7] = value
		 end
		 'extmaxSlider' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.ranges_to[7] = value
		 end
		 'zminSlider' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.ranges_from[8] = value
		 end
		 'zmaxSlider' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.ranges_to[8] = value
		 end
		 
		 'NITER_MCMC' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.niter_mcmc = float(value)
		 end
		 
; Synthesis
		 'sigmaSynthSlider' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.synth[0] = value
		 	plot_model_database, state
		 end
		 'YSynthSlider' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.synth[1] = value
		 	plot_model_database, state
		 end
		 'NSynthSlider' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.synth[2] = value
		 	plot_model_database, state
		 end
		 'qSynthSlider' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.synth[3] = value
		 	plot_model_database, state
		 end
		 'tauvSynthSlider' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.synth[4] = value
		 	plot_model_database, state
		 end
		 'angleSynthSlider' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.synth[5] = value
		 	plot_model_database, state
		 end
		 'shiftSynthSlider' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.synth[6] = value
		 	plot_model_database, state
		 end
		 'extSynthSlider' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.synth[7] = value
		 	plot_model_database, state
		 end
		 'zSynthSlider' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.synth[8] = value
		 	plot_model_database, state
		 end		 
		 'PLOT_SED_SYNTH' : begin
		 	plot_model_database, state
		 end
		 'PLOT_DISTRIBUTION' : begin
		 	state.sed_or_derivative = 0
		 	plot_model_database, state
		 end
		 'PLOT_DERIVATIVE' : begin
		 	state.sed_or_derivative = 1
		 	plot_model_database, state
		 end
		 'SIGMA_DERIVATIVE' : begin
		 	state.derivative_which_par = 0
		 	plot_model_database, state
		 end
		 'Y_DERIVATIVE' : begin
		 	state.derivative_which_par = 1
		 	plot_model_database, state
		 end
		 'N_DERIVATIVE' : begin
		 	state.derivative_which_par = 2
		 	plot_model_database, state
		 end
		 'q_DERIVATIVE' : begin
		 	state.derivative_which_par = 3
		 	plot_model_database, state
		 end
		 'tauv_DERIVATIVE' : begin
		 	state.derivative_which_par = 4
		 	plot_model_database, state
		 end
		 'ANGLE_DERIVATIVE' : begin
		 	state.derivative_which_par = 5
		 	plot_model_database, state
		 end
		 'RESPONSE_FUNCTION' : begin
		 	plot_response_function, state
		 end

; Priors
		 'UNIFORM_sigma' : begin
		 	state.prior[0] = 'U'
			widget_control, state.widthsigma, SENSITIVE=0
			widget_control, state.possigma, SENSITIVE=0
		 end
		 'GAUSSIAN_sigma' : begin
		 	state.prior[0] = 'G'
			widget_control, state.widthsigma, SENSITIVE=1
			widget_control, state.possigma, SENSITIVE=1
		 end
		 'DIRAC_sigma' : begin
		 	state.prior[0] = 'D'
		 	widget_control, state.widthsigma, SENSITIVE=0
			widget_control, state.possigma, SENSITIVE=1
		 end
		 'sigma_POS' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.gaussian_prior_info[0,0] = value
		 end
		 'sigma_WIDTH' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.gaussian_prior_info[1,0] = value
		 end
		 
		 'UNIFORM_Y' : begin
		 	state.prior[1] = 'U'
			widget_control, state.widthY, SENSITIVE=0
			widget_control, state.posY, SENSITIVE=0
		 end
		 'GAUSSIAN_Y' : begin
		 	state.prior[1] = 'G'
			widget_control, state.widthY, SENSITIVE=1
			widget_control, state.posY, SENSITIVE=1
		 end
		 'DIRAC_Y' : begin
		 	state.prior[1] = 'D'
		 	widget_control, state.widthY, SENSITIVE=0
			widget_control, state.posY, SENSITIVE=1
		 end
		 'Y_POS' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.gaussian_prior_info[0,1] = value
		 end
		 'Y_WIDTH' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.gaussian_prior_info[1,1] = value
		 end
		 
		 'UNIFORM_N' : begin
		 	state.prior[2] = 'U'
			widget_control, state.widthN, SENSITIVE=0
			widget_control, state.posN, SENSITIVE=0
		 end
		 'GAUSSIAN_N' : begin
		 	state.prior[2] = 'G'
			widget_control, state.widthN, SENSITIVE=1
			widget_control, state.posN, SENSITIVE=1
		 end
		 'DIRAC_N' : begin
		 	state.prior[2] = 'D'
			widget_control, state.widthN, SENSITIVE=0
			widget_control, state.posN, SENSITIVE=1
		 end
		 'N_POS' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.gaussian_prior_info[0,2] = value
		 end
		 'N_WIDTH' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.gaussian_prior_info[1,2] = value
		 end
		 
		 'UNIFORM_q' : begin
		 	state.prior[3] = 'U'
			widget_control, state.widthq, SENSITIVE=0
			widget_control, state.posq, SENSITIVE=0
		 end
		 'GAUSSIAN_q' : begin
		 	state.prior[3] = 'G'
			widget_control, state.widthq, SENSITIVE=1
			widget_control, state.posq, SENSITIVE=1
		 end
		 'DIRAC_q' : begin
		 	state.prior[3] = 'D'
			widget_control, state.widthq, SENSITIVE=0
			widget_control, state.posq, SENSITIVE=1
		 end
		 'q_POS' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.gaussian_prior_info[0,3] = value
		 end
		 'q_WIDTH' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.gaussian_prior_info[1,3] = value
		 end
		 
		 'UNIFORM_tauv' : begin
		 	state.prior[4] = 'U'
			widget_control, state.widthtauv, SENSITIVE=0
			widget_control, state.postauv, SENSITIVE=0
		 end
		 'GAUSSIAN_tauv' : begin
		 	state.prior[4] = 'G'
			widget_control, state.widthtauv, SENSITIVE=1
			widget_control, state.postauv, SENSITIVE=1
		 end
		 'DIRAC_tauv' : begin
		 	state.prior[4] = 'D'
			widget_control, state.widthtauv, SENSITIVE=0
			widget_control, state.postauv, SENSITIVE=1
		 end
		 'tauv_POS' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.gaussian_prior_info[0,4] = value
		 end
		 'tauv_WIDTH' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.gaussian_prior_info[1,4] = value
		 end
		 
		 'UNIFORM_angle' : begin
		 	state.prior[5] = 'U'
			widget_control, state.widthangle, SENSITIVE=0
			widget_control, state.posangle, SENSITIVE=0
		 end
		 'GAUSSIAN_angle' : begin
		 	state.prior[5] = 'G'
			widget_control, state.widthangle, SENSITIVE=1
			widget_control, state.posangle, SENSITIVE=1
		 end
		 'DIRAC_angle' : begin
		 	state.prior[5] = 'D'
			widget_control, state.widthangle, SENSITIVE=0
			widget_control, state.posangle, SENSITIVE=1
		 end
		 'angle_POS' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.gaussian_prior_info[0,5] = value
		 end
		 'angle_WIDTH' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.gaussian_prior_info[1,5] = value
		 end
		 
		 'UNIFORM_shift' : begin
		 	state.prior[6] = 'U'
			widget_control, state.widthshift, SENSITIVE=0
			widget_control, state.posshift, SENSITIVE=0
		 end
		 'GAUSSIAN_shift' : begin
		 	state.prior[6] = 'G'
			widget_control, state.widthshift, SENSITIVE=1
			widget_control, state.posshift, SENSITIVE=1
		 end
		 'DIRAC_shift' : begin
		 	state.prior[6] = 'D'
			widget_control, state.widthshift, SENSITIVE=0
			widget_control, state.posshift, SENSITIVE=1
		 end
		 'shift_POS' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.gaussian_prior_info[0,6] = value
		 end
		 'shift_WIDTH' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.gaussian_prior_info[1,6] = value
		 end

		 'UNIFORM_ext' : begin
		 	state.prior[7] = 'U'
			widget_control, state.widthext, SENSITIVE=0
			widget_control, state.posext, SENSITIVE=0
		 end
		 'GAUSSIAN_ext' : begin
		 	state.prior[7] = 'G'
			widget_control, state.widthext, SENSITIVE=1
			widget_control, state.posext, SENSITIVE=1
		 end
		 'DIRAC_ext' : begin
		 	state.prior[7] = 'D'
			widget_control, state.widthext, SENSITIVE=0
			widget_control, state.posext, SENSITIVE=1
		 end
		 'ext_POS' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.gaussian_prior_info[0,7] = value
		 end
		 'ext_WIDTH' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.gaussian_prior_info[1,7] = value
		 end

		 'UNIFORM_z' : begin
		 	state.prior[8] = 'U'
			widget_control, state.widthz, SENSITIVE=0
			widget_control, state.posz, SENSITIVE=0
		 end
		 'GAUSSIAN_z' : begin
		 	state.prior[8] = 'G'
			widget_control, state.widthz, SENSITIVE=1
			widget_control, state.posz, SENSITIVE=1
		 end
		 'DIRAC_z' : begin
		 	state.prior[8] = 'D'
			widget_control, state.widthz, SENSITIVE=0
			widget_control, state.posz, SENSITIVE=1
		 end
		 'z_POS' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.gaussian_prior_info[0,8] = value
		 end
		 'z_WIDTH' : begin
		 	widget_control, Event.id, GET_VALUE=value
		 	state.gaussian_prior_info[1,8] = value
		 end
		 
		 'REDDENING_LIST' : begin
		 	state.reddening_law = Event.index
; If no reddening is chosen, use a Delta dirac prior for Av
		 	if (state.reddening_law eq 0) then begin
		 		state.prior[7] = 'D'
		 		state.gaussian_prior_info[0,7] = 0.0
		 		state.gaussian_prior_info[1,7] = 0.0
		 		widget_control, state.diracextButton, /SET_BUTTON
		 	endif
		 	if (state.reddening_law ne 0) then begin
		 		state.prior[7] = 'U'
		 		widget_control, state.uniformextButton, /SET_BUTTON
		 	endif
		 end
		 'YES_AGN' : begin
		 	state.agn_plus_sed = 1
		 	plot_model_database, state
		 end
		 'NO_AGN' : begin
		 	state.agn_plus_sed = 0
		 	plot_model_database, state
		 end
		 'YES_AGN_INF' : begin
		 	state.agn_plus_sed = 1
		 end
		 'NO_AGN_INF' : begin
		 	state.agn_plus_sed = 0
		 end

		 'PHOTSPEC_STD' : begin
		 	state.avg_method = 0
		 end
		 'PHOTSPEC_BAY' : begin
		 	state.avg_method = 1
		 end
		 'NEURAL_INTERPOLATION': begin
			state.neural_interpolation = 0
			restore,'neural.idl'
			state.lambda = ptr_new(neural.lambda)
		 end
		 'LINEAR_INTERPOLATION': begin
			state.neural_interpolation = 1
			restore,'database.idl'
			state.lambda = ptr_new(database.lambda)
		 end
		 		 
	 endcase

	 widget_control, handler, SET_UVALUE=state
 	 save, state, filename='state.idl'

end
