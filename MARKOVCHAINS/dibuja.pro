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


pro bayesme_get_2d_prob,chain,x,y,z,z2,xmaxval=xmaxval,xminval=xminval,$
                        ymaxval=ymaxval,yminval=yminval

	if not keyword_set(xmaxval) then xmaxval = 0.0
	if not keyword_set(ymaxval) then ymaxval = 0.0
	if not keyword_set(xminval) then xminval = 0.0
	if not keyword_set(yminval) then yminval = 0.0

; Bins-X
	stepx = optbin(chain[0,*])
	xmin  = min(chain[0,*])
	if xminval ne 0.0 then xmin = xminval
	xmax  = max(chain[0,*])
	if xmaxval ne 0.0 then xmax = xmaxval
	nx    = round ( (xmax-xmin)/stepx )
	stepx = (xmax-xmin)/float(nx)

; Bins-Y
	stepy = optbin(chain[1,*])
	ymin  = min(chain[1,*])
	if yminval ne 0.0 then ymin = yminval
	ymax  = max(chain[1,*])
	if ymaxval ne 0.0 then ymax = ymaxval
	ny    = round ( (ymax-ymin)/stepy )
	stepy = (ymax-ymin)/float(ny)


; Variables
	x  = findgen(nx)*stepx + xmin
	y  = findgen(ny)*stepy + ymin
	z  = fltarr(nx,ny)
	z2 = fltarr(nx,ny)

; variables for the normalization
	nn = 50
	xx = findgen(nn)/float(nn-1) * (xmax-xmin) + xmin
	yy = findgen(nn)/float(nn-1) * (ymax-ymin) + ymin


; Doing histogram.
	sig = 1.2
	for i=0,nx-1 do begin
		print, i
   	for j=0,ny-1 do begin

; To save some time, I preselect those points within 6 sigmas.
        presel = where( abs(chain[0,*]-x[i]) le 6.0*(sig*stepx) or $
                        abs(chain[1,*]-y[j]) le 6.0*(sig*stepy), npresel)
        chosenx = chain[0,presel]
        choseny = chain[1,presel]

        cut = where(chosenx ge x[i]-stepx/2.0 and $
                    chosenx lt x[i]+stepx/2.0 and $
                    choseny ge y[j]-stepy/2.0 and $
                    choseny lt y[j]+stepy/2.0, ncut)

;         ; smoothing kernel: gaussian
        wei  = exp( -0.5 * ( (chain[0,presel]-x[i])^2./ (sig*stepx)^2. ) ) *  $
               exp( -0.5 * ( (chain[1,presel]-y[j])^2./ (sig*stepy)^2. ) )
        norm1 = int_tabulated(xx,exp( -0.5 * (xx-x[i])^2./(sig*stepx)^2.) )
        norm2 = int_tabulated(yy,exp( -0.5 * (yy-y[j])^2./(sig*stepy)^2.) )
        z[i,j] = total(wei) / (norm1*norm2)

        ; Clasico.
        z2[i,j] = ncut

    	endfor
	endfor

; Outputs
	z  = z/max(z)   ; using the gaussian kernel
	z2 = z2/max(z2) ; using the classical histogram

end


pro figura_paper_1d, file

	openr,2,file+'.hist1D'
	readf,2,nvar

   one_sigma = (1.d0-2.d0*(1.d0-gauss_pdf(1.)))*100.
	two_sigma = (1.d0-2.d0*(1.d0-gauss_pdf(2.)))*100.

	params = ['!7r!6','Y','N!i0!n','q','!7s!6!dV!n','i','shift','ext']
	params_names = ['sigma','Y','N','q','tauv','i','shift','ext']
	
	for i = 0, 6 do begin
 		abre_ps,'FIGURES/'+file+'_marginal_posterior_'+strtrim(string(i),2)+'.eps',/encaps
		readf,2,ind,n
		temp = fltarr(3,n)
		readf,2,temp
		x = reform(temp[0,*])
		yGauss = reform(temp[1,*])
		yStep = reform(temp[2,*])
		
		plot, x, smooth(yStep,3), psym=10, ycharsize=1.1,xcharsize=1.1, xran=[min(x),max(x)], xsty=1,$
			xtit=params[i],ytit='Normalized marginal posterior', thick=3

;		if (i ne 4 and i ne 3) then begin
			conf_limits_ccf,x,yStep,est,errup,errdown,one_sigma
			print, params_names[i]+'=', est, ' +/-(1s) ', errup, errdown
			ver, est, linestyle=2
			ver, x[cerca(yStep,max(yStep))]
			ver, est+errup, linestyle=1
			ver, est-errdown, linestyle=1
			conf_limits_ccf,x,yStep,est,errup,errdown,two_sigma
			print, params_names[i]+'=', est, ' +/-(2s) ', errup, errdown
			print, 'Mode : ', x[cerca(yStep,max(yStep))]
			print, 'Mean : ', int_tabulated(x,x*yStep)/int_tabulated(x,yStep)
;		endif
		
;		if (i eq 4 or i eq 3) then begin
			conf_limits_ccf_oneside,x,yStep,errup,one_sigma
			print, 'Upper limit for '+params_names[i]+'(1s)=', errup
			conf_limits_ccf_oneside,x,yStep,errup,two_sigma
			print, 'Upper limit for '+params_names[i]+'(2s)=', errup
			
			conf_limits_ccf_oneside,x,yStep,errup,100-one_sigma
			print, 'Lower limit for '+params_names[i]+'(1s)=', errup
			conf_limits_ccf_oneside,x,yStep,errup,100-two_sigma
			print, 'Lower limit for '+params_names[i]+'(2s)=', errup
			print, ' ' 
;		endif

 		cierra_ps

   endfor
   
   close,2

   !p.multi=0

   stop
end

pro figura_paper_2d, file

	
; Read the Markov chains
	openr,2,file+'.chain',/f77
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

   one_sigma = (1.d0-2.d0*(1.d0-gauss_pdf(1.)))*100.
	two_sigma = (1.d0-2.d0*(1.d0-gauss_pdf(2.)))*100.
	
	params = ['!7r!3','Y','N','q','!7s!3!dv!n','i','shift']
	params_names = ['sigma','Y','N','q','tauv','i','shift']

	left = [4, 4, 5, 1, 2, 2]
	right = [2, 5, 0, 3, 3, 1]

   for i = 0, n_elements(left)-1 do begin

   	abre_ps,'FIGURES/'+file+'_marginal_posterior2d_'+strtrim(string(params_names[left[i]]),2)+'_'+$
   		strtrim(string(params_names[right[i]]),2)+'.eps',/encaps
    	
   	chainx = chain[left[i],*]
   	chainy = chain[right[i],*]

   	deltax = max(chainx)-min(chainx)
   	deltay = max(chainy)-min(chainy)
   			
   	h = hist_2d(chainx,chainy,bin1=deltax/25,bin2=deltay/25,$
   		min1=min(chainx),max1=max(chainx),min2=min(chainy),max2=max(chainy))

		nx = n_elements(h[*,0])
		ny = n_elements(h[0,*])
				
   	x = findgen(nx)/(nx-1.d0)*deltax + min(chainx)
   	y = findgen(ny)/(ny-1.d0)*deltay + min(chainy)

   	h_smooth = filter_image(h,/iterate,smooth=5)

		tvframe, h_smooth/max(h_smooth)*100., xran=[min(x),max(x)], xsty=1, yran=[min(y),max(y)], ysty=1,$
			xtit=params[left[i]],ytit=params[right[i]], charsize=1.1, /bar, $
			btit='p('+strtrim(string(params[left[i]]),2)+','+$
   		strtrim(string(params[right[i]]),2)+')/p!dmax!n [%]'
   	contour, h_smooth/max(h_smooth)*100., x, y, levels=[100.d0-two_sigma,100.d0-one_sigma],/overpl,$
   		c_col=[255,255],c_thick=[3,3],c_line=[2,0],c_annotation=['95%','68%'],c_labels=[1,1]

   	cierra_ps
   			
   endfor

   close,2

   !p.multi=0

   stop
end
