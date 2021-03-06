;   16/01/2004   J.A.Rubi�
;   This program obtains the confidence limits for a likelihood function.
;   It reads the data from file "filein". If fromfile=1, then you can
;   directly supply the likelihood function through DATA.
;   It is based in the old conf_limits.pro, but now
;   the confidence limits are derived from the cumulative distribution
;   function.
;
;-
;

PRO conf_limits_ccf_oneside,xdata,ydata,errup,conf


; Keywords:
;if (NOT keyword_set(CONF)) then conf=-1


; Confidence level:
;if (conf lt 0) then begin
;    read,conf,prompt=' > Input confidence level required (%): '
;endif

xpdf_up  = conf/100.d0

if (xpdf_up  gt 1.0) then stop,'wrong value for CONF'

; Data:
a      = fltarr(2,n_elements(xdata))
a[0,*] = xdata[*]
a[1,*] = ydata[*]


; Variables:
xx = a[0,*]
yy = a[1,*] / max(a[1,*])


; Sorting:
ord = sort(xx)
xx = xx[ord]
yy = yy[ord]


; Interpolating by splines:
n    = n_elements(xx)
npto = 300
xmin = xx[0]
xmax = xx[n-1]

if (npto gt n) then begin
    
    x    = findgen(npto)/float(npto-1) * (xmax - xmin) + xmin
    ;y    = spline(xx,yy,x,10)
    y    = spline(xx,yy,x,1.0)
    n    = npto	

endif else begin

    x = xx
    y = yy

endelse

; Normalizing to the peak:
y = y/max(y)

; Computing the cumulative distribution function:
xF = x[1:n-1]
F  = fltarr(n-1)
norm = int_tabulated(x,y)
for i=1,n-1 do begin
	F[i-1] = int_tabulated(x[0:i],y[0:i]) / norm
endfor

; Computing the limits:
dist3 = abs( F - xpdf_up )
cut3  = where( dist3 eq min(dist3) )


upper = xF[cut3[0]]

errup   = upper


END

;   16/01/2004   J.A.Rubi�
;   This program obtains the confidence limits for a likelihood function.
;   It reads the data from file "filein". If fromfile=1, then you can
;   directly supply the likelihood function through DATA.
;   It is based in the old conf_limits.pro, but now
;   the confidence limits are derived from the cumulative distribution
;   function.
;
;-
;

PRO conf_limits_ccf,xdata,ydata, est,errup,errdown,conf


; Keywords:
;if (NOT keyword_set(CONF)) then conf=-1


; Confidence level:
;if (conf lt 0) then begin
;    read,conf,prompt=' > Input confidence level required (%): '
;endif

xpdf_low = 0.5 - conf/200.
xpdf_mid = 0.5
xpdf_up  = 0.5 + conf/200.

if (xpdf_low lt 0.0) then stop,'wrong value for CONF'
if (xpdf_up  gt 1.0) then stop,'wrong value for CONF'



; Data:
a      = fltarr(2,n_elements(xdata))
a[0,*] = xdata[*]
a[1,*] = ydata[*]


; Variables:
xx = a[0,*]
yy = a[1,*] / max(a[1,*])


; Sorting:
ord = sort(xx)
xx = xx[ord]
yy = yy[ord]


; Interpolating by splines:
n    = n_elements(xx)
npto = 300
xmin = xx[0]
xmax = xx[n-1]

if (npto gt n) then begin
    
    x    = findgen(npto)/float(npto-1) * (xmax - xmin) + xmin
    ;y    = spline(xx,yy,x,10)
    y    = spline(xx,yy,x,1.0)
    n    = npto	

endif else begin

    x = xx
    y = yy

endelse

; Normalizing to the peak:
y = y/max(y)



; Computing the cumulative distribution function:
xF = x[1:n-1]
F  = fltarr(n-1)
norm = int_tabulated(x,y)
for i=1,n-1 do begin
	F[i-1] = int_tabulated(x[0:i],y[0:i]) / norm
endfor



; Computing the limits:
dist1 = abs( F - xpdf_low )
cut1  = where( dist1 eq min(dist1) )

dist2 = abs( F - xpdf_mid )
cut2  = where( dist2 eq min(dist2) )

dist3 = abs( F - xpdf_up )
cut3  = where( dist3 eq min(dist3) )


lower = xF[cut1[0]]
est   = xF[cut2[0]]
upper = xF[cut3[0]]

errup   = upper - est
errdown = est - lower



;print,' (*) The result is: '
;print,est,' + ',errup,' - ',errdown



END


function near, x, value
 temp = x - value
 minim = min(abs(temp),ind)
 return, ind
end

pro bayes_get_2d_prob,chain,x,y,z,z2,xmaxval=xmaxval,xminval=xminval,$
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

; Read the Markov chains obtained with BayesCLUMPY and plot them, obtaining
; some information from them
; file -> name of the output file without extension
; e.g. readchain, 'test'
pro readchain, file

   one_sigma = (1.d0-2.d0*(1.d0-gauss_pdf(1.)))*100.
	two_sigma = (1.d0-2.d0*(1.d0-gauss_pdf(2.)))*100.

	params = ['!7r!6','Y','N','q','!7s!6!dv!n','i','shift','A!dV!n','z']
	params_names = ['sigma','Y','N','q','tauv','i','shift','extinction','redshift']

	openr,2,file+'.hist1D'	
   readf,2,nparam

	!p.multi = [0,3,3]
	for i = 0, nparam-1 do begin
		readf,2,loop,n
		x = fltarr(n)
		yStep = fltarr(n)
		yGauss = fltarr(n)
		temp = fltarr(3,n)
		readf,2,temp
		x = reform(temp[0,*])
		yGauss = reform(temp[1,*])
		yStep = reform(temp[2,*])
		
		plot, x, smooth(yStep,3), psym=10, charsize=1.6, xran=[min(x),max(x)], xsty=1,$
			xtit=params[i],ytit='Normalized marginal posterior', thick=3		

; GT		if (i ne 4 and i ne 3 and i ne 5) then begin

; Locate parameters with a Dirac prior
			step = x - shift(x,1)
			if (max(step) ne 0) then begin
				conf_limits_ccf,x,yStep,est,errup,errdown,one_sigma
				print, params_names[i]+'=', est, ' +/-(1s) ', errup, errdown
				conf_limits_ccf,x,yStep,est,errup,errdown,two_sigma
				print, params_names[i]+'=', est, ' +/-(2s) ', errup, errdown
				print, 'Mode : ', x[near(yStep,max(yStep))]
				print, 'Mean : ', int_tabulated(x,x*yStep)/int_tabulated(x,yStep)
; 			endif
;GT		endif
		
; GT		if (i eq 4 or i eq 3 or i eq 5) then begin
				conf_limits_ccf_oneside,x,yStep,errup,one_sigma
				print, 'Upper limit for '+params_names[i]+'(1s)=', errup
				conf_limits_ccf_oneside,x,yStep,errup,two_sigma
				print, 'Upper limit for '+params_names[i]+'(2s)=', errup
			endif
;GT		endif
   endfor
   
   close,2

   !p.multi=0

   stop
end
