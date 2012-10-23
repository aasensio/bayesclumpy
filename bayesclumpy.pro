@bayesclumpy_ann
@bayesclumpy_routines
@bayesclumpy_plot
@bayesclumpy_init
@bayesclumpy_event
@bayesclumpy_experiment
@showprogress__define

;-----------------------------------------
; Main routine
;-----------------------------------------
pro bayesclumpy, reset_state=reset_state
	!p.multi = 0
	device, decomposed=0, retain=2, true_color=24
	loadct,0,/silent
	close,/all
	state = bayesclumpy_init(reset_state=reset_state)
	 
	widget_control, state.baseWidget, /REALIZE
	 
	xmanager, 'bayesclumpy', state.baseWidget, EVENT_HANDLER='BayesClumpy_Event'
end
