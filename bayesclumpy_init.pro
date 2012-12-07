;-----------------------------------------
; Initialization routine. Draws the widget
;-----------------------------------------
function bayesclumpy_init, reset_state=reset_state

	spec = ''
	if (file_test('BayesClumpy.version')) then begin
		openr,2,'BayesClumpy.version'
		readf,2,version_bc
		close,2
		spec = 's'
	endif else begin
		version_bc = 2.0
	endelse
	version = 'BayesCLUMPY v'+strtrim(string(abs(version_bc),FORMAT='(I1)'),2)+spec
	print, 'Starting '+ version
	if (keyword_set(reset_state) or file_test('state.idl') eq 0) then begin
		state = {baseWidget: 0L, obsfile: '', fileText: 0L, obs_x: ptr_new(), obs_y: ptr_new(), $
			obs_sigma: ptr_new(), plotleftWidget: 0L, plotrightWidget: 0L, $
			ranges_from: fltarr(9), ranges_to: fltarr(9), niter_mcmc: 25000L, burnin: 40.0, $
			posSigma: 0L, widthSigma: 0L, posY: 0L, widthY: 0L, posN: 0L, widthN: 0L,$
			posq: 0L, widthq: 0L, postauv: 0L, widthtauv: 0L, posangle: 0L, widthangle: 0L, posext: 0L,$
			widthext: 0L, posshift: 0L, widthshift: 0L, prior: replicate('U',9), gaussian_prior_info: fltarr(2,9),$
			agn_name: '', posz: 0L, widthz: 0L, $
			params: ['sigma','Y','N','q','tauv','angle','shift','extinction','redshift'], outfile: 'test', $
			plotSEDWidget: 0L, $
			synth: fltarr(9), reddening_law: 0, agn_plus_sed: 0, uniformextButton: 0L, gaussianextButton: 0L,$
			diracextButton: 0L, sed_or_derivative: 0, derivative_which_par: 0, avg_method: 0,$
			version: version_bc, filters: ptr_new(), distance: 0.d0, neural_interpolation: 0, lambda: ptr_new()}
		state.synth = [15,5.0,1.0,0.0,10.0,0.0,0.0,0.0,0.0]
	endif else begin
		restore,'state.idl'
	endelse

	slider_xsize = 110

	list_reddenings = ['No reddening','Allen (1976)','Seaton (1979) MW',$
		'Fitzpatrick (1986) LMC',$
		'Pr�vot et al. (1984) SMC','Calzetti et al. (2000) SB','Chiar & Tielens (2000)']

 	symbo = findgen(17) * (!pi*2/16.)
	usersym, cos(symbo), sin(symbo), /fill

	state.baseWidget = widget_base(TITLE=version)

; If not initialized...
	if (total(abs(state.ranges_from)) eq 0) then begin
		state.ranges_from = [15,5,1,0,5,0,-2,0,0]
		state.ranges_to = [70,30,15,3,150,90,2,10,6]
	endif

; Initialize data for the neural network
; The wavelength axis is returned in neural.lambda
	neural = initialize_neural_networks()
	save,neural,file='neural.idl'
	
; Initialize some things of the full database.
	database = initialize_database()	
	save,database,file='database.idl'
	
; Generate the full wavelength axis
	if (state.neural_interpolation eq 0) then begin
		state.lambda = ptr_new(neural.lambda)
	endif else begin
		state.lambda = ptr_new(database.lambda)
	endelse
	
	device, get_screen_size=winsize
	ysiz = 350
	if (winsize[1] le 800) then begin
		ysiz = 300
	endif
	xsiz = 500
	if (winsize[1] le 800) then begin
		xsiz = 550
	endif
	
	tabs = widget_tab(state.baseWidget, VALUE='tabs',UVALUE='TABS')
	
	tab2 = widget_base(tabs, TITLE='Inference', /COLUMN, UVALUE='TAB2')
	tab1 = widget_base(tabs, TITLE='Synthesis', /COLUMN, UVALUE='TAB1')
	
; **********************
; SYNTHESIS TAB
; **********************

; Plotting windows
	horizplotBase = widget_base(tab1, /ROW)
	state.plotSEDWidget = widget_draw(horizplotBase, XSIZE=550, YSIZE=350, /FRAME)

; Parameters
	horizparamBase = widget_base(tab1, /COLUMN)
	set1paramBase = widget_base(horizparamBase, /ROW)
	set2paramBase = widget_base(horizparamBase, /ROW)

; sigma
	sigmaBase = widget_base(set1paramBase, /COLUMN, FRAME=1, TITLE='sigma')
	sigmaSlider = cw_fslider(sigmaBase,TITLE='sigma',UVALUE='sigmaSynthSlider',$
	 	  XSIZE=150,MINIMUM=15,MAXIMUM=70,VALUE=state.synth[0], /EDIT)
	
; Y
	YBase = widget_base(set1paramBase, /COLUMN, FRAME=1, TITLE='Y')
	YSlider = cw_fslider(YBase,TITLE='Y',UVALUE='YSynthSlider',$
	 	  XSIZE=150,MINIMUM=5,MAXIMUM=100,VALUE=state.synth[1], /EDIT)
	
; N
	NBase = widget_base(set1paramBase, /COLUMN, FRAME=1, TITLE='N')
	NSlider = cw_fslider(NBase,TITLE='N',UVALUE='NSynthSlider',$
	 	  XSIZE=150,MINIMUM=1,MAXIMUM=15,VALUE=state.synth[2], /EDIT)
	 	  
; q
	qBase = widget_base(set1paramBase, /COLUMN, FRAME=1, TITLE='q')
	qSlider = cw_fslider(qBase,TITLE='q',UVALUE='qSynthSlider',$
	 	  XSIZE=150,MINIMUM=0,MAXIMUM=3,VALUE=state.synth[3], /EDIT)
	 	  
; tauv
	tauvBase = widget_base(set2paramBase, /COLUMN, FRAME=1, TITLE='tauv')
	tauvSlider = cw_fslider(tauvBase,TITLE='tauv',UVALUE='tauvSynthSlider',$
	 	  XSIZE=150,MINIMUM=5,MAXIMUM=150,VALUE=state.synth[4], /EDIT)
	 	  
; angle
	angleBase = widget_base(set2paramBase, /COLUMN, FRAME=1, TITLE='angle')
	angleSlider = cw_fslider(angleBase,TITLE='angle',UVALUE='angleSynthSlider',$
	 	  XSIZE=150,MINIMUM=0,MAXIMUM=90,VALUE=state.synth[5], /EDIT)

; shift
	shiftBase = widget_base(set2paramBase, /COLUMN, FRAME=1, TITLE='Sigma')
	shiftSlider = cw_fslider(shiftBase,TITLE='shift',UVALUE='shiftSynthSlider',$
	 	  XSIZE=150,MINIMUM=-4,MAXIMUM=4,VALUE=state.synth[6], /EDIT)

; extinction
	extBase = widget_base(set2paramBase, /COLUMN, FRAME=1, TITLE='A_V')
	extSlider = cw_fslider(extBase,TITLE='extinction',UVALUE='extSynthSlider',$
	 	  XSIZE=150,MINIMUM=0,MAXIMUM=10,VALUE=state.synth[7], /EDIT)

; redshift
	zBase = widget_base(set2paramBase, /COLUMN, FRAME=1, TITLE='z')
	zSlider = cw_fslider(zBase,TITLE='z',UVALUE='zSynthSlider',$
	 	  XSIZE=150,MINIMUM=0,MAXIMUM=6,VALUE=state.synth[8])
	 	  
	plotspecButton = widget_button(set2paramBase, VALUE='Plot SED',$
		UVALUE='PLOT_SED_SYNTH')
	plotspecButton = widget_button(set2paramBase, VALUE='Plot Response',$
		UVALUE='RESPONSE_FUNCTION')

	optionsBase = widget_base(tab1, /ROW)
	redBase = widget_base(optionsBase, /COLUMN)
	redLabel = widget_label(redBase, VALUE='Reddening law : ')
   listList = widget_combobox(redBase, VALUE=list_reddenings, $
        UVALUE='REDDENING_LIST',SENSITIVE=1)
   widget_control, listList, SET_COMBOBOX_SELECT=state.reddening_law


; Widget jujutsu to create fancy label frames.
   IF !VERSION.OS_FAMILY EQ 'WINDOWS' THEN fancyFont = 'Times*16*Italic*Bold'
   bulletinBoardBase = Widget_Base(optionsBase)
   thisLabel = ' AGN contribution '
   label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)
   labelGeometry = Widget_Info(label, /GEOMETRY)
   labelYSize =  labelGeometry.ysize
   agnBase = Widget_Base(bulletinBoardBase, COLUMN=1, /FRAME, $
        YOFFSET=labelYSize/2, YPAD=10, XPAD=10, /EXCLUSIVE)
   label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)
; Widget jujutsu to create fancy label frames.

   yesButton = widget_button(agnBase, VALUE='SED+AGN', UVALUE='YES_AGN')
   noButton = widget_button(agnBase, VALUE='Only SED', UVALUE='NO_AGN')
   case (state.agn_plus_sed) of
   	0 : widget_control, noButton, /SET_BUTTON
   	1 : widget_control, yesButton, /SET_BUTTON
   endcase


; Widget jujutsu to create fancy label frames.
   IF !VERSION.OS_FAMILY EQ 'WINDOWS' THEN fancyFont = 'Times*16*Italic*Bold'
   bulletinBoardBase = Widget_Base(optionsBase)
   thisLabel = ' SED/Derivative '
   label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)
   labelGeometry = Widget_Info(label, /GEOMETRY)
   labelYSize =  labelGeometry.ysize
   agnBase = Widget_Base(bulletinBoardBase, COLUMN=1, /FRAME, $
        YOFFSET=labelYSize/2, YPAD=10, XPAD=10, /EXCLUSIVE)
   label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)
; Widget jujutsu to create fancy label frames.

   yesButton = widget_button(agnBase, VALUE='SED', UVALUE='PLOT_DISTRIBUTION')
   noButton = widget_button(agnBase, VALUE='Derivative', UVALUE='PLOT_DERIVATIVE')
   case (state.sed_or_derivative) of
   	0 : widget_control, yesButton, /SET_BUTTON
   	1 : widget_control, noButton, /SET_BUTTON
   endcase

; Widget jujutsu to create fancy label frames.
   IF !VERSION.OS_FAMILY EQ 'WINDOWS' THEN fancyFont = 'Times*16*Italic*Bold'
   bulletinBoardBase = Widget_Base(optionsBase)
   thisLabel = ' Derivative with respect to'
   label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)
   labelGeometry = Widget_Info(label, /GEOMETRY)
   labelYSize =  labelGeometry.ysize
   agnBase = Widget_Base(bulletinBoardBase, COLUMN=1, /FRAME, $
        YOFFSET=labelYSize/2, YPAD=10, XPAD=10, /EXCLUSIVE)
   label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)
; Widget jujutsu to create fancy label frames.

   sigmaButton = widget_button(agnBase, VALUE='Sigma', UVALUE='SIGMA_DERIVATIVE')
   YButton = widget_button(agnBase, VALUE='Y', UVALUE='Y_DERIVATIVE')
   NButton = widget_button(agnBase, VALUE='N', UVALUE='N_DERIVATIVE')
   qButton = widget_button(agnBase, VALUE='q', UVALUE='q_DERIVATIVE')
   tauvButton = widget_button(agnBase, VALUE='tauv', UVALUE='tauv_DERIVATIVE')
   angleButton = widget_button(agnBase, VALUE='angle', UVALUE='ANGLE_DERIVATIVE')
   case (state.derivative_which_par) of
   	0 : widget_control, sigmaButton, /SET_BUTTON
   	1 : widget_control, YButton, /SET_BUTTON
   	2 : widget_control, NButton, /SET_BUTTON
   	3 : widget_control, qButton, /SET_BUTTON
   	4 : widget_control, tauvButton, /SET_BUTTON
   	5 : widget_control, angleButton, /SET_BUTTON
   endcase

	 
; **********************
; INFERENCE TAB
; **********************
	bayesclumpyBase = widget_base(tab2, /COLUMN, /BASE_ALIGN_CENTER)

; Plotting windows
	horizplotBase = widget_base(bayesclumpyBase, /ROW)
	state.plotleftWidget = widget_draw(horizplotBase, XSIZE=xsiz, YSIZE=ysiz, /FRAME)
	state.plotrightWidget = widget_draw(horizplotBase, XSIZE=xsiz, YSIZE=ysiz, /FRAME)

; Parameters
	horizparamBase = widget_base(bayesclumpyBase, /ROW)

; sigma
	sigmaBase = widget_base(horizparamBase, /COLUMN, FRAME=1, TITLE='Sigma')
	sigmaminSlider = widget_slider(sigmaBase,TITLE='sigma [min]',UVALUE='sigmaminSlider',$
	 	  XSIZE=slider_xsize,MINIMUM=15,MAXIMUM=70,VALUE=state.ranges_from[0])
	sigmamaxSlider = widget_slider(sigmaBase,TITLE='sigma [max]',UVALUE='sigmamaxSlider',$
	 	  XSIZE=slider_xsize,MINIMUM=15,MAXIMUM=70,VALUE=state.ranges_to[0])
	priorBase = widget_base(sigmaBase, /COLUMN, /EXCLUSIVE)
	uniformButton = widget_button(priorBase, VALUE='Uniform', UVALUE='UNIFORM_sigma')
	gaussianButton = widget_button(priorBase, VALUE='Gaussian', UVALUE='GAUSSIAN_sigma')
	diracButton = widget_button(priorBase, VALUE='Dirac', UVALUE='DIRAC_sigma')
	case (state.prior[0]) of
		'U': begin
				widget_control, uniformButton, /SET_BUTTON
				sens = [0,0]
			  end
		'G': begin
				widget_control, gaussianButton, /SET_BUTTON
				sens = [1,1]
			  end
		'D': begin
				widget_control, diracButton, /SET_BUTTON
				sens = [1,0]
			  end
	endcase
	labelposSigma = widget_label(sigmaBase, VALUE='Center value:')
	state.posSigma = widget_text(sigmaBase, VALUE=strtrim(string(state.gaussian_prior_info[0,0]),2),$
		UVALUE='sigma_POS',/EDITABLE,XSIZE=8,YSIZE=1,SENSITIVE=sens[0])
	labelwidthSigma = widget_label(sigmaBase, VALUE='Width:')
	state.widthSigma = widget_text(sigmaBase, VALUE=strtrim(string(state.gaussian_prior_info[1,0]),2),$
		UVALUE='sigma_WIDTH',/EDITABLE,XSIZE=8,YSIZE=1,SENSITIVE=sens[1])
	
; Y
	YBase = widget_base(horizparamBase, /COLUMN, FRAME=1, TITLE='Sigma')
	YminSlider = widget_slider(YBase,TITLE='Y [min]',UVALUE='YminSlider',$
	 	  XSIZE=slider_xsize,MINIMUM=5,MAXIMUM=100,VALUE=state.ranges_from[1])
	YmaxSlider = widget_slider(YBase,TITLE='Y [max]',UVALUE='YmaxSlider',$
	 	  XSIZE=slider_xsize,MINIMUM=5,MAXIMUM=100,VALUE=state.ranges_to[1])
	priorBase = widget_base(YBase, /COLUMN, /EXCLUSIVE)
	uniformButton = widget_button(priorBase, VALUE='Uniform', UVALUE='UNIFORM_Y')
	gaussianButton = widget_button(priorBase, VALUE='Gaussian', UVALUE='GAUSSIAN_Y')
	diracButton = widget_button(priorBase, VALUE='Dirac', UVALUE='DIRAC_Y')
	case (state.prior[1]) of
		'U': begin
				widget_control, uniformButton, /SET_BUTTON
				sens = [0,0]
			  end
		'G': begin
				widget_control, gaussianButton, /SET_BUTTON
				sens = [1,1]
			  end
		'D': begin
				widget_control, diracButton, /SET_BUTTON
				sens = [1,0]
			  end
	endcase
	labelposY = widget_label(YBase, VALUE='Center value:')
	state.posY = widget_text(YBase, VALUE=strtrim(string(state.gaussian_prior_info[0,1]),2),$
		UVALUE='Y_POS',/EDITABLE,XSIZE=8,YSIZE=1,SENSITIVE=sens[0])
	labelwidthY = widget_label(YBase, VALUE='Width:')
	state.widthY = widget_text(YBase, VALUE=strtrim(string(state.gaussian_prior_info[1,1]),2),$
		UVALUE='Y_WIDTH',/EDITABLE,XSIZE=8,YSIZE=1,SENSITIVE=sens[1])

; N
	NBase = widget_base(horizparamBase, /COLUMN, FRAME=1, TITLE='Sigma')
	NminSlider = widget_slider(NBase,TITLE='N [min]',UVALUE='NminSlider',$
	 	  XSIZE=slider_xsize,MINIMUM=1,MAXIMUM=15,VALUE=state.ranges_from[2])
	NmaxSlider = widget_slider(NBase,TITLE='N [max]',UVALUE='NmaxSlider',$
	 	  XSIZE=slider_xsize,MINIMUM=1,MAXIMUM=15,VALUE=state.ranges_to[2])
	priorBase = widget_base(NBase, /COLUMN, /EXCLUSIVE)
	uniformButton = widget_button(priorBase, VALUE='Uniform', UVALUE='UNIFORM_N')
	gaussianButton = widget_button(priorBase, VALUE='Gaussian', UVALUE='GAUSSIAN_N')
	diracButton = widget_button(priorBase, VALUE='Dirac', UVALUE='DIRAC_N')
	case (state.prior[2]) of
		'U': begin
				widget_control, uniformButton, /SET_BUTTON
				sens = [0,0]
			  end
		'G': begin
				widget_control, gaussianButton, /SET_BUTTON
				sens = [1,1]
			  end
		'D': begin
				widget_control, diracButton, /SET_BUTTON
				sens = [1,0]
			  end
	endcase
	labelposN = widget_label(NBase, VALUE='Center value:')
	state.posN = widget_text(NBase, VALUE=strtrim(string(state.gaussian_prior_info[0,2]),2),$
		UVALUE='N_POS',/EDITABLE,XSIZE=8,YSIZE=1,SENSITIVE=sens[0])
	labelwidthN = widget_label(NBase, VALUE='Width:')
	state.widthN = widget_text(NBase, VALUE=strtrim(string(state.gaussian_prior_info[1,2]),2),$
		UVALUE='N_WIDTH',/EDITABLE,XSIZE=8,YSIZE=1,SENSITIVE=sens[1])

; q
	qBase = widget_base(horizparamBase, /COLUMN, FRAME=1, TITLE='Sigma')
	qminSlider = widget_slider(qBase,TITLE='q [min]',UVALUE='qminSlider',$
	 	  XSIZE=slider_xsize,MINIMUM=0,MAXIMUM=3,VALUE=state.ranges_from[3])
	qmaxSlider = widget_slider(qBase,TITLE='q [max]',UVALUE='qmaxSlider',$
	 	  XSIZE=slider_xsize,MINIMUM=0,MAXIMUM=3,VALUE=state.ranges_to[3])
	priorBase = widget_base(qBase, /COLUMN, /EXCLUSIVE)
	uniformButton = widget_button(priorBase, VALUE='Uniform', UVALUE='UNIFORM_q')
	gaussianButton = widget_button(priorBase, VALUE='Gaussian', UVALUE='GAUSSIAN_q')
	diracButton = widget_button(priorBase, VALUE='Dirac', UVALUE='DIRAC_q')
	case (state.prior[3]) of
		'U': begin
				widget_control, uniformButton, /SET_BUTTON
				sens = [0,0]
			  end
		'G': begin
				widget_control, gaussianButton, /SET_BUTTON
				sens = [1,1]
			  end
		'D': begin
				widget_control, diracButton, /SET_BUTTON
				sens = [1,0]
			  end
	endcase
	labelposq = widget_label(qBase, VALUE='Center value:')
	state.posq = widget_text(qBase, VALUE=strtrim(string(state.gaussian_prior_info[0,3]),2),$
		UVALUE='q_POS',/EDITABLE,XSIZE=8,YSIZE=1,SENSITIVE=sens[0])
	labelwidthq = widget_label(qBase, VALUE='Width:')
	state.widthq = widget_text(qBase, VALUE=strtrim(string(state.gaussian_prior_info[1,3]),2),$
		UVALUE='q_WIDTH',/EDITABLE,XSIZE=8,YSIZE=1,SENSITIVE=sens[1])

; tauv
	tauvBase = widget_base(horizparamBase, /COLUMN, FRAME=1, TITLE='Sigma')
	tauvminSlider = widget_slider(tauvBase,TITLE='tauv [min]',UVALUE='tauvminSlider',$
	 	  XSIZE=slider_xsize,MINIMUM=5,MAXIMUM=150,VALUE=state.ranges_from[4])
	tauvmaxSlider = widget_slider(tauvBase,TITLE='tauv [max]',UVALUE='tauvmaxSlider',$
	 	  XSIZE=slider_xsize,MINIMUM=5,MAXIMUM=150,VALUE=state.ranges_to[4])
	priorBase = widget_base(tauvBase, /COLUMN, /EXCLUSIVE)
	uniformButton = widget_button(priorBase, VALUE='Uniform', UVALUE='UNIFORM_tauv')
	gaussianButton = widget_button(priorBase, VALUE='Gaussian', UVALUE='GAUSSIAN_tauv')
	diracButton = widget_button(priorBase, VALUE='Dirac', UVALUE='DIRAC_tauv')
	case (state.prior[4]) of
		'U': begin
				widget_control, uniformButton, /SET_BUTTON
				sens = [0,0]
			  end
		'G': begin
				widget_control, gaussianButton, /SET_BUTTON
				sens = [1,1]
			  end
		'D': begin
				widget_control, diracButton, /SET_BUTTON
				sens = [1,0]
			  end
	endcase
	labelpostauv = widget_label(tauvBase, VALUE='Center value:')
	state.postauv = widget_text(tauvBase, VALUE=strtrim(string(state.gaussian_prior_info[0,4]),2),$
		UVALUE='tauv_POS',/EDITABLE,XSIZE=8,YSIZE=1,SENSITIVE=sens[0])
	labelwidthtauv = widget_label(tauvBase, VALUE='Width:')
	state.widthtauv = widget_text(tauvBase, VALUE=strtrim(string(state.gaussian_prior_info[1,4]),2),$
		UVALUE='tauv_WIDTH',/EDITABLE,XSIZE=8,YSIZE=1,SENSITIVE=sens[1])

; angle
	angleBase = widget_base(horizparamBase, /COLUMN, FRAME=1, TITLE='Sigma')
	angleminSlider = widget_slider(angleBase,TITLE='angle [min]',UVALUE='angleminSlider',$
	 	  XSIZE=slider_xsize,MINIMUM=0,MAXIMUM=90,VALUE=state.ranges_from[5])
	anglemaxSlider = widget_slider(angleBase,TITLE='angle [max]',UVALUE='anglemaxSlider',$
	 	  XSIZE=slider_xsize,MINIMUM=0,MAXIMUM=90,VALUE=state.ranges_to[5])
	priorBase = widget_base(angleBase, /COLUMN, /EXCLUSIVE)
	uniformButton = widget_button(priorBase, VALUE='Uniform', UVALUE='UNIFORM_angle')
	gaussianButton = widget_button(priorBase, VALUE='Gaussian', UVALUE='GAUSSIAN_angle')
	diracButton = widget_button(priorBase, VALUE='Dirac', UVALUE='DIRAC_angle')
	case (state.prior[5]) of
		'U': begin
				widget_control, uniformButton, /SET_BUTTON
				sens = [0,0]
			  end
		'G': begin
				widget_control, gaussianButton, /SET_BUTTON
				sens = [1,1]
			  end
		'D': begin
				widget_control, diracButton, /SET_BUTTON
				sens = [1,0]
			  end
	endcase
	labelposangle = widget_label(angleBase, VALUE='Center value:')
	state.posangle = widget_text(angleBase, VALUE=strtrim(string(state.gaussian_prior_info[0,5]),2),$
		UVALUE='angle_POS',/EDITABLE,XSIZE=8,YSIZE=1,SENSITIVE=sens[0])
	labelwidthangle = widget_label(angleBase, VALUE='Width:')
	state.widthangle = widget_text(angleBase, VALUE=strtrim(string(state.gaussian_prior_info[1,5]),2),$
		UVALUE='angle_WIDTH',/EDITABLE,XSIZE=8,YSIZE=1,SENSITIVE=sens[1])

; shift
	shiftBase = widget_base(horizparamBase, /COLUMN, FRAME=1, TITLE='Sigma')
	shiftminSlider = widget_slider(shiftBase,TITLE='shift [min]',UVALUE='shiftminSlider',$
	 	  XSIZE=slider_xsize,MINIMUM=-4,MAXIMUM=4,VALUE=state.ranges_from[6])
	shiftmaxSlider = widget_slider(shiftBase,TITLE='shift [max]',UVALUE='shiftmaxSlider',$
	 	  XSIZE=slider_xsize,MINIMUM=-4,MAXIMUM=4,VALUE=state.ranges_to[6])
	priorBase = widget_base(shiftBase, /COLUMN, /EXCLUSIVE)
	uniformButton = widget_button(priorBase, VALUE='Uniform', UVALUE='UNIFORM_shift')
	gaussianButton = widget_button(priorBase, VALUE='Gaussian', UVALUE='GAUSSIAN_shift')
	diracButton = widget_button(priorBase, VALUE='Dirac', UVALUE='DIRAC_shift')
	case (state.prior[6]) of
		'U': begin
				widget_control, uniformButton, /SET_BUTTON
				sens = [0,0]
			  end
		'G': begin
				widget_control, gaussianButton, /SET_BUTTON
				sens = [1,1]
			  end
		'D': begin
				widget_control, diracButton, /SET_BUTTON
				sens = [1,0]
			  end
	endcase
	labelposshift = widget_label(shiftBase, VALUE='Center value:')
	state.posshift = widget_text(shiftBase, VALUE=strtrim(string(state.gaussian_prior_info[0,6]),2),$
		UVALUE='shift_POS',/EDITABLE,XSIZE=8,YSIZE=1,SENSITIVE=sens[0])
	labelwidthshift = widget_label(shiftBase, VALUE='Width:')
	state.widthshift = widget_text(shiftBase, VALUE=strtrim(string(state.gaussian_prior_info[1,6]),2),$
		UVALUE='shift_WIDTH',/EDITABLE,XSIZE=8,YSIZE=1,SENSITIVE=sens[1])

; Extinction
	extBase = widget_base(horizparamBase, /COLUMN, FRAME=1, TITLE='A_V')
	extminSlider = widget_slider(extBase,TITLE='A_V [min]',UVALUE='extminSlider',$
	 	  XSIZE=slider_xsize,MINIMUM=0,MAXIMUM=40,VALUE=state.ranges_from[7])
	extmaxSlider = widget_slider(extBase,TITLE='A_V [max]',UVALUE='extmaxSlider',$
	 	  XSIZE=slider_xsize,MINIMUM=0,MAXIMUM=40,VALUE=state.ranges_to[7])
	priorBase = widget_base(extBase, /COLUMN, /EXCLUSIVE)
	state.uniformextButton = widget_button(priorBase, VALUE='Uniform', UVALUE='UNIFORM_ext')
	state.gaussianextButton = widget_button(priorBase, VALUE='Gaussian', UVALUE='GAUSSIAN_ext')
	state.diracextButton = widget_button(priorBase, VALUE='Dirac', UVALUE='DIRAC_ext')
	case (state.prior[7]) of
		'U': begin
				widget_control, state.uniformextButton, /SET_BUTTON
				sens = [0,0]
			  end
		'G': begin
				widget_control, state.gaussianextButton, /SET_BUTTON
				sens = [1,1]
			  end
		'D': begin
				widget_control, state.diracextButton, /SET_BUTTON
				sens = [1,0]
			  end
	endcase
	labelposext = widget_label(extBase, VALUE='Center value:')
	state.posext = widget_text(extBase, VALUE=strtrim(string(state.gaussian_prior_info[0,7]),2),$
		UVALUE='ext_POS',/EDITABLE,XSIZE=8,YSIZE=1,SENSITIVE=sens[0])
	labelwidthext = widget_label(extBase, VALUE='Width:')
	state.widthext = widget_text(extBase, VALUE=strtrim(string(state.gaussian_prior_info[1,7]),2),$
		UVALUE='ext_WIDTH',/EDITABLE,XSIZE=8,YSIZE=1,SENSITIVE=sens[1])

; Redshift
	extBase = widget_base(horizparamBase, /COLUMN, FRAME=1, TITLE='z')
	extminSlider = widget_slider(extBase,TITLE='z [min]',UVALUE='zminSlider',$
	 	  XSIZE=slider_xsize,MINIMUM=0,MAXIMUM=6,VALUE=state.ranges_from[8])
	extmaxSlider = widget_slider(extBase,TITLE='z [max]',UVALUE='zmaxSlider',$
	 	  XSIZE=slider_xsize,MINIMUM=0,MAXIMUM=6,VALUE=state.ranges_to[8])
	priorBase = widget_base(extBase, /COLUMN, /EXCLUSIVE)
	uniformzButton = widget_button(priorBase, VALUE='Uniform', UVALUE='UNIFORM_z')
	gaussianzButton = widget_button(priorBase, VALUE='Gaussian', UVALUE='GAUSSIAN_z')
	diraczButton = widget_button(priorBase, VALUE='Dirac', UVALUE='DIRAC_z')
	case (state.prior[8]) of
		'U': begin
				widget_control, uniformzButton, /SET_BUTTON
				sens = [0,0]
			  end
		'G': begin
				widget_control, gaussianzButton, /SET_BUTTON
				sens = [1,1]
			  end
		'D': begin
				widget_control, diraczButton, /SET_BUTTON
				sens = [1,0]
			  end
	endcase
	labelposext = widget_label(extBase, VALUE='Center value:')
	state.posz = widget_text(extBase, VALUE=strtrim(string(state.gaussian_prior_info[0,8]),2),$
		UVALUE='z_POS',/EDITABLE,XSIZE=8,YSIZE=1,SENSITIVE=sens[0])
	labelwidthext = widget_label(extBase, VALUE='Width:')
	state.widthz = widget_text(extBase, VALUE=strtrim(string(state.gaussian_prior_info[1,8]),2),$
		UVALUE='z_WIDTH',/EDITABLE,XSIZE=8,YSIZE=1,SENSITIVE=sens[1])

; MCMC parameters
	mcmcBase = widget_base(horizparamBase, /COLUMN)
; 	labelniter = widget_label(mcmcBase, VALUE='Number of iterations:')
; 	niter = widget_text(mcmcBase, VALUE=strtrim(string(state.niter_mcmc),2),$
; 		UVALUE='NITER_MCMC',/EDITABLE,XSIZE=6,YSIZE=1)
; 	labelniter = widget_label(mcmcBase, VALUE='Burn-in [%]:')
; 	burnin = widget_text(mcmcBase, VALUE=strtrim(string(state.burnin),2),$
; 		UVALUE='BURNIN_MCMC',/EDITABLE,XSIZE=6,YSIZE=1)
	labelniter = widget_label(mcmcBase, VALUE='Output file:')
	outfile = widget_text(mcmcBase, VALUE=strtrim(string(state.outfile),2),$
		UVALUE='OUTFILE',/EDITABLE,XSIZE=6,YSIZE=1)
	
	labelniter = widget_label(mcmcBase, VALUE='Distance [Mpc]')
	outfile = widget_text(mcmcBase, VALUE=strtrim(string(state.distance),2),$
		UVALUE='DISTANCE',/EDITABLE,XSIZE=6,YSIZE=1)

	redLabel = widget_label(mcmcBase, VALUE='Reddening law : ')
   listList = widget_combobox(mcmcBase, VALUE=list_reddenings, $
        UVALUE='REDDENING_LIST',SENSITIVE=1)
   widget_control, listList, SET_COMBOBOX_SELECT=state.reddening_law

	mcmcBase2 = widget_base(mcmcBase, /ROW)
	
; Widget jujutsu to create fancy label frames.
   IF !VERSION.OS_FAMILY EQ 'WINDOWS' THEN fancyFont = 'Times*16*Italic*Bold'
   bulletinBoardBase = Widget_Base(mcmcBase2)
   thisLabel = ' AGN contribution '
   label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)
   labelGeometry = Widget_Info(label, /GEOMETRY)
   labelYSize =  labelGeometry.ysize
   agnBase = Widget_Base(bulletinBoardBase, COLUMN=1, /FRAME, $
        YOFFSET=labelYSize/2, YPAD=10, XPAD=10, /EXCLUSIVE)
   label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)
; Widget jujutsu to create fancy label frames.

   yesButton = widget_button(agnBase, VALUE='SED+AGN', UVALUE='YES_AGN_INF')
   noButton = widget_button(agnBase, VALUE='Only SED', UVALUE='NO_AGN_INF')
   case (state.agn_plus_sed) of
   	0 : widget_control, noButton, /SET_BUTTON
   	1 : widget_control, yesButton, /SET_BUTTON
   endcase

; Widget jujutsu to create fancy label frames.
   IF !VERSION.OS_FAMILY EQ 'WINDOWS' THEN fancyFont = 'Times*16*Italic*Bold'
   bulletinBoardBase = Widget_Base(mcmcBase2)
   thisLabel = ' Phot+spec '
   label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)
   labelGeometry = Widget_Info(label, /GEOMETRY)
   labelYSize =  labelGeometry.ysize
   agnBase = Widget_Base(bulletinBoardBase, COLUMN=1, /FRAME, $
        YOFFSET=labelYSize/2, YPAD=10, XPAD=10, /EXCLUSIVE)
   label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)
; Widget jujutsu to create fancy label frames.

   stdButton = widget_button(agnBase, VALUE='Standard', UVALUE='PHOTSPEC_STD')
   mnestButton = widget_button(agnBase, VALUE='Bayesian avg.', UVALUE='PHOTSPEC_BAY')
   case (state.avg_method) of
   	0 : widget_control, stdButton, /SET_BUTTON
   	1 : widget_control, mnestButton, /SET_BUTTON
   endcase
   
mcmcBase2 = widget_base(mcmcBase, /ROW)
	
; Widget jujutsu to create fancy label frames.
   IF !VERSION.OS_FAMILY EQ 'WINDOWS' THEN fancyFont = 'Times*16*Italic*Bold'
   bulletinBoardBase = Widget_Base(mcmcBase2)
   thisLabel = ' Interpolation method '
   label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)
   labelGeometry = Widget_Info(label, /GEOMETRY)
   labelYSize =  labelGeometry.ysize
   agnBase = Widget_Base(bulletinBoardBase, COLUMN=1, /FRAME, $
        YOFFSET=labelYSize/2, YPAD=10, XPAD=10, /EXCLUSIVE)
   label = Widget_Label(bulletinBoardBase, VALUE=thisLabel, XOFFSET=5, FONT=fancyFont)
; Widget jujutsu to create fancy label frames.

   yesButton = widget_button(agnBase, VALUE='Neural network', UVALUE='NEURAL_INTERPOLATION')
   noButton = widget_button(agnBase, VALUE='Linear with DB', UVALUE='LINEAR_INTERPOLATION')
   case (state.neural_interpolation) of
   	0 : widget_control, yesButton, /SET_BUTTON
   	1 : widget_control, noButton, /SET_BUTTON
   endcase
	
	
; Observations
	horizobsBase = widget_base(bayesclumpyBase, /ROW, /ALIGN_LEFT)
	fileLabel = widget_label(horizobsBase, VALUE='Observations : ')
   state.fileText = widget_text(horizobsBase, VALUE=state.obsfile, /EDITABLE, XSIZE=70,YSIZE=1, UVALUE='OBS_FILE')
   fileButton = widget_button(horizobsBase, VALUE='Select file', UVALUE='OBS_FILE_DIALOG')
   plotSEDButton = widget_button(horizobsBase, VALUE='Plot SED', UVALUE='PLOT_SED')
   plotMAPButton = widget_button(horizobsBase, VALUE='Plot MAP', UVALUE='PLOT_MAP')
   plotinferenceButton = widget_button(horizobsBase, VALUE='Plot last inference', UVALUE='PLOT_INFERENCE')
   plotinferenceButton = widget_button(horizobsBase, VALUE='App. covering factor', UVALUE='CALCULATE_COVFACTOR')
   suggestobsButton = widget_button(horizobsBase, VALUE='Suggest observation', UVALUE='SUGGEST_OBSERVATION')
	inferenceButton = widget_button(horizobsBase, VALUE='Inference', UVALUE='DO_INFERENCE')
   
	if (state.reddening_law eq 0) then begin
		state.prior[7] = 'D'
		state.gaussian_prior_info[0,7] = 0.0
		state.gaussian_prior_info[1,7] = 0.0
	endif

; Try to read observations if they exist
	read_observations, state.obsfile, state
	
	widget_control, state.baseWidget, SET_UVALUE=state
	 
	return, state
end
