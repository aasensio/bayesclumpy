;+
; NAME:
;    PSWINDOW
;
; PURPOSE:
;
;    This function is used to calculate the size of a PostScript
;    window that has the same aspect ratio (ratio of height to
;    width) as the current display graphics window. It creates
;    the largest possible PostScript output window with the
;    desired aspect ratio. This assures that graphics output
;    looks similar, if not identical, to PostScript output.
;
; AUTHOR:
;
;   FANNING SOFTWARE CONSULTING
;   David Fanning, Ph.D.
;   2642 Bradbury Court
;   Fort Collins, CO 80521 USA
;   Phone: 970-221-0438
;   E-mail: davidf@dfanning.com
;   Coyote's Guide to IDL Programming: http://www.dfanning.com/
;
; CATEGORY:
;
;    Graphics.
;
; CALLING SEQUENCE:
;
;    pageInfo = PSWINDOW()
;
; INPUTS:
;
;    None.
;
; KEYWORD PARAMETERS:
;
;    CM: Normally the structure value that is returned from this
;    function reports its values in inches. Setting this keyword
;    causes the return values to be in units of centimeters.
;
;    LANDSCAPE: Normally this function assumes a PostScript window
;    in Portrait mode. Setting this keyword assumes you want
;    the graphic in Landscape mode.
;
;    MARGIN:  The margin around the edges of the plot. The value must be
;    a floating point value between 0.0 and 0.5. It is expressed in
;    normalized coordinate units. The default margin is 0.15.
;
;    PAGESIZE: Set this keyword to a string indicating the type
;    of PostScript page size you want. Current values are "LETTER",
;    "LEGAL", and "A4". Default is "LETTER".
;
;    PRINTER: Set this keyword if the output will be used to
;    configure the PRINTER device, rather than the PS device.
;    (In the PRINTER device, offsets are always calculated from
;    the lower-left corner of the page and do not rotate in
;    Landscape mode, as they do with the PS device.) Note that
;    the PRINTER device is only able to accept these keywords
;    in IDL 5.1 and higher.
;
; OUTPUTS:
;
;    pageInfo: The output value is a named structure defined like
;    this:
;
;      pageInfo = {PSWINDOW_STRUCT, XSIZE:0.0, YSIZE:0.0, $
;         XOFSET:0.0, YOFFSET:0.0, INCHES:0, PORTRAIT:0, LANDSCAPE:0}
;
;    The units of the four size fields are inches unless the CM
;    keyword is set.
;
;    The output can be used to immediately configure the PostScript
;    or Printer device, like this:
;
;       Set_Plot, 'PS' ; or 'PRINTER'
;       Device, _Extra=pageInfo
;
; RESTRICTIONS:
;
;    The aspect ratio of the current graphics window is calculated
;    like this:
;
;       aspectRatio = FLOAT(!D.Y_VSIZE) / !D.X_VSIZE
;
; EXAMPLE:
;
;    To create a PostScript output window with the same aspect
;    ratio as the curently active display window, type:
;
;     pageInfo = PSWINDOW()
;     SET_PLOT, 'PS'
;     DEVICE, _Extra=pageInfo
;
;     To configure the PRINTER device:
;
;     pageInfo = PSWINDOW(/Printer)
;     SET_PLOT, 'PRINTER'
;     DEVICE, _Extra=pageInfo
;
; MODIFICATION HISTORY:
;
;    Written by: David Fanning, November 1996.
;       Fixed a bug in which the YOFFSET was calculated incorrectly
;          in Landscape mode. 12 Feb 97.
;       Took out a line of code that wasn't being used. 14 Mar 97.
;       Added correct units keyword to return structure. 29 JUN 98. DWF
;       Fixed a bug in how landscape offsets were calculated. 19 JUL 99. DWF.
;       Fixed a bug in the way margins were used to conform to my
;          original conception of the program. 19 JUL 99. DWF.
;       Added Landscape and Portrait fields to the return structure. 19 JUL 99. DWF.
;       Added PageSize keyword, changed MARGIN keyword, and completely
;          rewrote most of the intenal code. 9 FEB 2000. DWF.
;       Fixed a bug in how I calculated the aspect ratio. 1 MAR 2000. DWF.
;       Added PRINTER keyword to return proper offset values for the
;          PRINTER device, where the offset location is not rotated. 1 MAR 2000. DWF.
;-
;
;###########################################################################
;
; LICENSE
;
; This software is OSI Certified Open Source Software.
; OSI Certified is a certification mark of the Open Source Initiative.
;
; Copyright � 2000 Fanning Software Consulting
;
; This software is provided "as-is", without any express or
; implied warranty. In no event will the authors be held liable
; for any damages arising from the use of this software.
;
; Permission is granted to anyone to use this software for any
; purpose, including commercial applications, and to alter it and
; redistribute it freely, subject to the following restrictions:
;
; 1. The origin of this software must not be misrepresented; you must
;    not claim you wrote the original software. If you use this software
;    in a product, an acknowledgment in the product documentation
;    would be appreciated, but is not required.
;
; 2. Altered source versions must be plainly marked as such, and must
;    not be misrepresented as being the original software.
;
; 3. This notice may not be removed or altered from any source distribution.
;
; For more information on Open Source Software, visit the Open Source
; web site: http://www.opensource.org.
;
;###########################################################################



FUNCTION PSWINDOW_ASPECT, aspectRatio, MARGIN=margin, WindowAspect=wAspectRatio

ON_ERROR, 1

   ; Check for aspect ratio parameter and possibilities.

IF N_PARAMS() EQ 0 THEN aspectRatio = 1.0

IF aspectRatio EQ 0 THEN BEGIN
   MESSAGE, 'Aspect Ratio of 0. Changing to 1...', /Informational
   aspectRatio = 1.0
ENDIF

s = SIZE(aspectRatio)
IF s(s(0)+1) NE 4 THEN $
   MESSAGE, 'Aspect Ratio is not a FLOAT. Take care...', /Informational

   ; Check for margins.

IF N_ELEMENTS(margin) EQ 0 THEN margin = 0.15

   ; Error checking.

IF margin LT 0 OR margin GE 0.5 THEN $
   MESSAGE, 'The MARGIN keyword value must be between 0.0 and 0.5.'

   ; Calculate the aspect ratio of the current window.

IF N_Elements(wAspectRatio) EQ 0 THEN wAspectRatio = FLOAT(!D.Y_VSIZE) / !D.X_VSIZE

   ; Calculate normalized positions in window.

IF (aspectRatio LE wAspectRatio) THEN BEGIN
   xstart = margin
   ystart = 0.5 - (0.5 - margin) * (aspectRatio / wAspectRatio)
   xend = 1.0 - margin
   yend = 0.5 + (0.5 - margin) * (aspectRatio / wAspectRatio)
ENDIF ELSE BEGIN
   xstart = 0.5 - (0.5 - margin) * (wAspectRatio / aspectRatio)
   ystart = margin
   xend = 0.5 + (0.5 - margin) * (wAspectRatio / aspectRatio)
   yend = 1.0 - margin
ENDELSE

position = [xstart, ystart, xend, yend]

RETURN, position
END ; ----------------------------------------------------------------------------------



FUNCTION PSWINDOW, LANDSCAPE=landscape, CM=cm, MARGIN=margin, $
   PageSize=pagesize, Printer=printer

   ; Set up default values.

landscape = Keyword_Set(landscape)
cm = Keyword_Set(cm)
printer = Keyword_Set(printer)
inches = 1

   ; Get the page size.

IF N_Elements(pagesize) EQ 0 THEN pagesize = 'LETTER' $
   ELSE pagesize = StrUpCase(pagesize)
CASE pagesize OF
   'LETTER': BEGIN
      shortside = 8.5
      longside = 11.0
      ENDCASE
   'LEGAL': BEGIN
      shortside = 8.5
      longside = 14.0
      ENDCASE
    'A4': BEGIN
      shortside = 8.27
      longside = 11.7
      ENDCASE
    ELSE: BEGIN
      Message, 'Unknown page size. Using LETTER...', /Informational
      shortside = 8.5
      longside = 11.0
      ENDCASE
ENDCASE

   ; Need measurements in centimeters?

IF KEYWORD_SET(cm) THEN BEGIN
      shortside = shortside * 2.54
      longside = longside * 2.54
      inches = 0
ENDIF

   ; Determine the margin of the window on the page.

IF N_ELEMENTS(margin) EQ 0 THEN margin=0.15

   ; Get the aspect ratio of the current display window. Aspect ratio
   ; is ratio of xsize/ysize.

aspectRatio = FLOAT(!D.Y_VSIZE) / !D.X_VSIZE

   ; Get the aspect ratio of the page.

IF Keyword_Set(landscape) THEN pAspectRatio = shortside / longside $
   ELSE pAspectRatio = longside / shortside

   ; Get the position on the page for this window.

pos = PSWindow_Aspect(aspectRatio, Margin=margin, WindowAspect=pAspectRatio)

   ; Convert normalized position coordinates to size units.

IF KEYWORD_SET(landscape) THEN BEGIN
   IF printer THEN BEGIN
      xsize = (pos[2] - pos[0]) * longside
      ysize = (pos[3] - pos[1]) * shortside
      yoffset = pos[1] * shortside
      xoffset = pos[0] * longside
      landscape = 1
      portrait = 0
   ENDIF ELSE BEGIN
      xsize = (pos[2] - pos[0]) * longside
      ysize = (pos[3] - pos[1]) * shortside
      xoffset = pos[1] * shortside
      yoffset = longside - (pos[0] * longside)
      landscape = 1
      portrait = 0
   ENDELSE
ENDIF ELSE BEGIN
   xsize = (pos[2] - pos[0]) * shortside
   ysize = (pos[3] - pos[1]) * longside
   xoffset = pos[0] * shortside
   yoffset = pos[1] * longside
   landscape = 0
   portrait = 1
ENDELSE

   ; Return the proper DEVICE data structure.

RETURN, {PSWINDOW_STRUCT, XSIZE:xsize, YSIZE:ysize, $
   XOFFSET:xoffset, YOFFSET:yoffset, INCHES:inches, $
   PORTRAIT:portrait, LANDSCAPE:landscape}

END

;**********************************
; Close a postcript file
;**********************************
pro close_ps, FONT=font

; Reset the postcript window	

	device,/inches,xoffset=3./4,yoffset=5,xsize=7,ysize=5
	device,/close
	set_plot,'x'
end

;**********************************
; Open a postcript file
;**********************************
pro open_ps, nombre, TODO=todo, ENCAPSULATED=encapsulated, COLOR=color, LANDSCAPE=landscape,$
	FONT=font, MITAD=mitad, SMALLFONT=smallfont, CENTER=center, XOFFSET=XOFFSET
	
	if (not keyword_set(XOFFSET)) then begin
		xoffset=0.5
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


pro verx, x, _extra=e
	tam = !y.crange
 	ymin = tam(0)
 	ymax = tam(1)
 	if (!y.type eq 1) then begin
   	ymin = 10.d0^tam(0)
      ymax = 10.d0^tam(1)
 	endif
 	plots,x,ymin
 	plots,x,ymax,/continue, _extra=e
end

pro cwpal
;
	r=[0,255,234,0  ,0  ,255,0  ,255,255,0  ,175,255]
	g=[0,255,0  ,140,216,235,235,0  ,148,126,0  ,67 ]
	b=[0,255,0  ,234,0  ,0  ,228,200,0  ,0  ,201,67 ]
	tvlct, r, g, b
;
end

FUNCTION TSUM,X,Y,IMIN,IMAX               ;Trapezoidal summation
;+
; NAME:
;       TSUM
; PURPOSE:
;       Trapezoidal summation of the area under a curve. 
; EXPLANATION:
;       Adapted from the procedure INTEG in the IUE procedure library.  
;
; CALLING SEQUENCE:
;       Result = TSUM(y)
;              or
;       Result = TSUM( x, y, [ imin, imax ] )  
; INPUTS:
;       x = array containing monotonic independent variable.  If omitted, then
;               x is assumed to contain the index of the y variable.
;               x = lindgen( N_elements(y) ).
;       y = array containing dependent variable y = f(x)
;
; OPTIONAL INPUTS:
;       imin = scalar index of x array at which to begin the integration
;               If omitted, then summation starts at x[0].
;       imax = scalar index of x value at which to end the integration 
;               If omitted then the integration ends at x[npts-1].
;
; OUTPUTS:
;       result = area under the curve y=f(x) between x[imin] and x[imax].
;
; EXAMPLE:
;       IDL> x = [0.0,0.1,0.14,0.3] 
;       IDL> y = sin(x)
;       IDL> print,tsum(x,y)    ===>  0.0445843
;       
;       In this example, the exact curve can be computed analytically as 
;       1.0 - cos(0.3) = 0.0446635     
; PROCEDURE:
;       The area is determined of individual trapezoids defined by x[i],
;       x[i+1], y[i] and y[i+1].
;
;       If the data is known to be at all smooth, then a more accurate
;       integration can be found by interpolation prior to the trapezoidal
;       sums, for example, by the standard IDL User Library int_tabulated.pro.
; MODIFICATION HISTORY:
;       Written, W.B. Landsman, STI Corp. May 1986
;       Modified so X is not altered in a one parameter call Jan 1990
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Allow non-integer values of imin and imax  W. Landsman April 2001
;       Fix problem if only 1 parameter supplied W. Landsman June 2002
;-
; Set default parameters
 On_error,2
 npar = N_params()
   
 if npar EQ 1 then begin
    npts = N_elements(x)
    yy = x
    xx = lindgen(npts)
    ilo = 0   & imin = ilo
    ihi = npts-1 & imax = ihi
 endif else begin

   if ( npar LT 3 ) then imin = 0
   npts = min( [N_elements(x), N_elements(y)] )
   if ( npar LT 4 ) then imax = npts-1
   ilo = long(imin)
   ihi = long(imax)
   xx = x[ilo:ihi]
   yy = y[ilo:ihi]
   npts = ihi - ilo + 1
 endelse         
;
; Compute areas of trapezoids and sum result
;
  xdif = xx[1:*] - xx
  yavg =  ( yy[0:npts-2] + yy[1:npts-1] ) / 2.  
  sum = total( xdif*yavg ) 

; Now account for edge effects if IMIN or IMAX parameter are not integers

  hi = imax - ihi
  lo = imin - ilo
  if (ihi LT imax) then sum = sum + (x[ihi+1]-x[ihi])*hi* $
              (y[ihi] + (hi/2.) *(y[ihi+1] - y[ihi]) )
  if (ilo LT imin) then sum = sum - (x[ilo+1]-x[ilo])*lo* $
              (y[ilo] + (lo/2.) *(y[ilo+1] - y[ilo]) )
  return, sum

  end

;**********************************
; Verify consistency between executed Markov chain and the one to plot
;**********************************
function verify_consistency, type
	tmp = ''
	openr,2,'end.info'
	readf,2,tmp
	close,2

	verify_consistency = 1
	if (type eq 0 and strpos(tmp,'Multinest') ne -1) then verify_consistency = 0
	if (type eq 1 and strpos(tmp,'MCMC') ne -1) then verify_consistency = 0
	return, verify_consistency
end

;**********************************
; Verify if sampling has finished
;**********************************
function verify_finished
	tmp = ''
	openr,2,'end.info'
	readf,2,tmp
	close,2

	verify_finished = 0
	if (strpos(tmp,'finished') ne -1) then verify_finished = 1
	return, verify_finished
end

;**********************************
; Return the perc percetile of an array
;**********************************
function percentile, data, perc
	sidx = sort(data)
	ndata = n_elements(data)
	return, data[sidx[perc*ndata / 100]]
end

;**********************************
; Calculate the histogram of an array
;**********************************
function histog, a, nbins=nbins, binsize=binsize, plot=plot, optimbin=optimbin, min=min, max=max, _extra=_extra
	if (not keyword_set(min)) then begin
		mina = min(a)
	endif else begin
		mina = min
	endelse
	if (not keyword_set(max)) then begin
		maxa = max(a)
	endif else begin
		maxa = max
	endelse
	if (keyword_set(optimbin)) then begin
		binopt = optbin(a)
		h = histogram(a,binsize=binopt,min=mina,max=maxa)
		n = n_elements(h)
		x = findgen(n) * binopt + mina
		if (keyword_set(plot)) then begin
	 		plot,x,h, psym=10, _extra=_extra
	 	endif
		return,[[x],[h]]
	endif
		
	if (keyword_set(nbins)) then begin
		h = histogram(a,nbins=nbins,min=mina,max=maxa)	 
	 	x = findgen(nbins)/(nbins-1.d0) * (maxa-mina) + mina
	 	if (keyword_set(plot)) then begin
	 		plot,x,h, psym=10, _extra=_extra
	 	endif
	 	return,[[x],[h]]
	endif
	 
	if (keyword_set(binsize)) then begin
		h = histogram(a,binsize=binsize,min=mina,max=maxa)
		n = n_elements(h)
		x = findgen(n) * binsize + mina
		if (keyword_set(plot)) then begin
	 		plot,x,h, psym=10, _extra=_extra
	 	endif
		return,[[x],[h]]
	endif
	 	 
end

; Clean '' from a string
function clean, str
	temp = strpos(str,"'")
	if (temp ne -1) then begin
		strput,str,' ',temp
	endif
	temp = strpos(str,"'")
	if (temp ne -1) then begin
		strput,str,' ',temp
	endif
	return, strcompress(str, /remove_all)
end