pro abre_ps, nombre, TODO=todo, ENCAPSULATED=encapsulated, COLOR=color, LANDSCAPE=landscape,$
	FONT=font, MITAD=mitad, SMALLFONT=smallfont, CENTER=center, XOFFSET=XOFFSET
	
	if (not keyword_set(XOFFSET)) then begin
		xoffset=0.5
	endif
	if (keyword_set(FONT)) then begin
	 	get_lun,u
	 	a = !p.charsize
	   openw,u,'~/idl/fuente.dat'
	 	printf,u,a
	 	close,u	
		free_lun,u
	 	!p.charsize=1.4
	endif
	if (keyword_set(SMALLFONT)) then begin
	 	get_lun,u
	 	a = !p.charsize
	   openw,u,'~/idl/fuente.dat'
	 	printf,u,a
	 	close,u	
		free_lun,u
	 	!p.charsize=1.0
	endif
	a = pswindow(/cm)
	set_plot,'ps',/copy
	encap = 0
	col = 0
	landsca = 0
	if (keyword_set(ENCAPSULATED)) then encap=1
	if (keyword_set(COLOR)) then col = 1
	if (keyword_set(LANDSCAPE)) then landsca = 1
	if (not keyword_set(TODO)) then begin
	 	  if (keyword_set(MITAD)) then begin
		   	device, filename=nombre, ENCAPSULATED=encap, COLOR=col, $
	 	  			LANDSCAPE=landsca, bits_per_pixel=8,xoffset=XOFFSET,yoffset=14.5,$
	 	  			xsize=20.5-XOFFSET,ysize=14	 	  			
		  endif else if (keyword_set(CENTER)) then begin		  
		   	device, filename=nombre, ENCAPSULATED=encap, COLOR=col, $
	 	  			LANDSCAPE=landsca, bits_per_pixel=8, $
					XSIZE=a.xsize, YSIZE=a.ysize, XOFFSET=a.xoffset, $
    	   		 YOFFSET=a.yoffset
		  endif else begin
		   	device, filename=nombre, ENCAPSULATED=encap, COLOR=col, $
	 	  			LANDSCAPE=landsca, bits_per_pixel=8
		  endelse
	 endif else begin
		device, filename=nombre,xoffset=XOFFSET,yoffset=1.5,xsize=19.5,ysize=26, ENCAPSULATED=encap,$
				COLOR=col, LANDSCAPE=landsca, bits_per_pixel=8
	 endelse
end


pro cierra_ps, FONT=font

; Reset the postcript window	

	device,/inches,xoffset=3./4,yoffset=5,xsize=7,ysize=5
	device,/close
	set_plot,'x'
	if (keyword_set(FONT)) then begin
	 	get_lun,u
	    openr,u,'~/idl/fuente.dat'
	 	readf,u,a
	 	close,u
		free_lun,u
	 	!p.charsize=a
	endif
end


pro extinction, file

; Read the Markov chains
	nlines = file_lines(strtrim(string(file),2)+'post_equal_weights.dat')
	chain = fltarr(10,nlines)
	
; Read the Markov chains
	openr,2,strtrim(string(file),2)+'post_equal_weights.dat'
   readf,2,chain
   close,2

   one_sigma = (1.d0-2.d0*(1.d0-gauss_pdf(1.)))*100.
	two_sigma = (1.d0-2.d0*(1.d0-gauss_pdf(2.)))*100.

	params = ['!7r!6','Y','N','q','!7s!6!dv!n','i','shift','ext']
	params_names = ['sigma','Y','N','q','tauv','i','shift','ext']

	N0=chain[2,*]
	tau=chain[4,*]
	sigma=chain[0,*]
	i=chain[5,*]
	
	Av=1.086*(N0*tau*exp(-(90-i)^2/sigma^2))
		
	
	hist=histog(Av,nbins=80,/plot)
	x=reform(hist[*,0])
	ystep=reform(hist[*,1])

 	abre_ps,'FIGURES/'+file+'_marginal_posterior_Av.eps',/encaps

	
	plot, x, smooth(ystep,3)/2.5e+4, psym=10, ycharsize=1.6,xcharsize=1.9, xran=[min(Av),max(Av)], xsty=1,$
	xtit='Av!ELOS!N',ytit='Normalized marginal posterior', thick=3
	
	
	conf_limits_ccf,x,yStep,est,errup,errdown,one_sigma
	print, 'Av=', est, ' +/-(1s) ', errup, errdown
	conf_limits_ccf,x,yStep,est,errup,errdown,two_sigma
	print, 'Av=', est, ' +/-(2s) ', errup, errdown
	print, 'Mode : ', x[cerca(yStep,max(yStep))]
	print, 'Mean : ', int_tabulated(x,x*yStep)/int_tabulated(x,yStep)

			ver, est, linestyle=2
			ver, x[cerca(yStep,max(yStep))]
			ver, est+errup, linestyle=1
			ver, est-errdown, linestyle=1

 	cierra_ps

	conf_limits_ccf_oneside,x,yStep,errup,one_sigma
	print, 'Upper limit for Av (1s)=', errup
	conf_limits_ccf_oneside,x,yStep,errup,two_sigma
	print, 'Upper limit for Av (2s)=', errup
	
	conf_limits_ccf_oneside,x,yStep,errup,100-one_sigma
	print, 'Lower limit for Av (1s)=', errup
	conf_limits_ccf_oneside,x,yStep,errup,100-two_sigma
	print, 'Lower limit for Av (2s)=', errup
	print, ' ' 


;	stop
	
END	



