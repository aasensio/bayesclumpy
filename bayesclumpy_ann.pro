;-----------------------------------------
; Return the AGN spectrum. Using Rowan-Robinson (1995)
;-----------------------------------------
function lambda_flambda_agn, lambda, cte
	lambdah = 0.01d0
   lambdau = 0.1d0
   lambdaRJ = 1.d0
   p = 0.5d0

   res = fltarr(n_elements(lambda))
        
   ind = where(lambda le lambdah)
	if (ind[0] ne -1) then begin
		res[ind] = cte*lambda[ind]^1.2 / lambdah^1.2
	endif

	ind = where(lambda gt lambdah and lambda le lambdau)
	if (ind[0] ne -1) then begin
		res[ind] = cte
	endif

	ind = where(lambda gt lambdau and lambda le lambdaRJ)
	if (ind[0] ne -1) then begin
		res[ind] = cte * lambda[ind]^(-p) / lambdau^(-p)
	endif

	ind = where(lambda gt lambdaRJ)
	if (ind[0] ne -1) then begin
		res[ind] = cte * lambdaRJ^(-p) / lambdau^(-p) * lambda[ind]^(-3) / lambdaRJ^(-3)
	endif

	return, res
end

; *****************************************
; Trapezoidal quadrature
; *****************************************
function tsum,x,y,imin,imax
	on_error,2
 	npar = n_params()
   
 	if npar eq 1 then begin
    	npts = n_elements(x)
    	yy = x
    	xx = lindgen(npts)
    	ilo = 0   & imin = ilo
    	ihi = npts-1 & imax = ihi
 	endif else begin

   	if ( npar lt 3 ) then imin = 0
   	npts = min( [n_elements(x), n_elements(y)] )
   	if ( npar lt 4 ) then imax = npts-1
   	ilo = long(imin)
   	ihi = long(imax)
   	xx = x[ilo:ihi]
   	yy = y[ilo:ihi]
   	npts = ihi - ilo + 1
 	endelse         
;
; compute areas of trapezoids and sum result
;
  	xdif = xx[1:*] - xx
  	yavg =  ( yy[0:npts-2] + yy[1:npts-1] ) / 2.  
  	sum = total( xdif*yavg ) 

; now account for edge effects if imin or imax parameter are not integers

  	hi = imax - ihi
  	lo = imin - ilo
  	if (ihi lt imax) then sum = sum + (x[ihi+1]-x[ihi])*hi* $
   	(y[ihi] + (hi/2.) *(y[ihi+1] - y[ihi]) )
  	if (ilo lt imin) then sum = sum - (x[ilo+1]-x[ilo])*lo* $
      (y[ilo] + (lo/2.) *(y[ilo+1] - y[ilo]) )
	return, sum
end
  
; *****************************************
; Activation function for the neural network
; *****************************************
function sigma, x	
	return, tanh(x)
end

; *****************************************
; Derivative of the activation function for the neural network
; *****************************************
function sigma_deriv, x	
	return, 1.d0-tanh(x)^2
end

; *****************************************
; Return the index of vector x closer to the given value
; *****************************************
function near, x, value
 temp = x - value
 minim = min(abs(temp),ind)
 return, ind
end

; *****************************************
; Evaluate a neural network given by the structure neuralnet and applied to a
; vector of input parameters
; *****************************************
function eval_neuralnetwork, neuralnet, input_data

	ndata = n_elements(input_data[0,*])
	ndim = n_elements(input_data[*,0])

	ninput = neuralnet.ninput
	nhidden = (*neuralnet.nhidden)
	noutput = neuralnet.noutput

	output = dblarr(noutput,ndata)
				
	for loop_input = 0L, ndata-1 do begin
		hidden = dblarr(nhidden)
		for i = 0, nhidden-1 do begin
			for j = 0, ninput-1 do begin
				hidden[i] = hidden[i] + (*neuralnet.input_hidden)[j,i] * input_data[j,loop_input]
			endfor
			hidden[i] = sigma(hidden[i] + (*neuralnet.bias)[i])
		endfor
		
		
		for i = 0, noutput-1 do begin
			for j = 0, nhidden-1 do begin
				output[i,loop_input] = output[i,loop_input] + (*neuralnet.hidden_output)[j,i] * hidden[j]
			endfor
		endfor
	endfor
	
	return, output
end

; *****************************************
; Evaluate a neural network given by the structure neuralnet and applied to a
; vector of input parameters
; *****************************************
function eval_neuralnetwork_derivative, neuralnet, input_data, which_par

	ndata = n_elements(input_data[0,*])
	ndim = n_elements(input_data[*,0])

	ninput = neuralnet.ninput
	nhidden = (*neuralnet.nhidden)
	noutput = neuralnet.noutput

	output = dblarr(noutput,ndata)
				
	for loop_input = 0L, ndata-1 do begin
		hidden = dblarr(nhidden)
		for i = 0, nhidden-1 do begin
			for j = 0, ninput-1 do begin
				hidden[i] = hidden[i] + (*neuralnet.input_hidden)[j,i] * input_data[j,loop_input]
			endfor
			hidden[i] = (*neuralnet.input_hidden)[which_par,i] * sigma_deriv(hidden[i] + (*neuralnet.bias)[i])			
		endfor
		
		
		for i = 0, noutput-1 do begin
			for j = 0, nhidden-1 do begin
				output[i,loop_input] = output[i,loop_input] + (*neuralnet.hidden_output)[j,i] * hidden[j]
			endfor
		endfor
	endfor
	
	return, output
end

; *****************************************
; Initialize some things of the DB. Value of the parameters and PCA eigenvectors
; *****************************************
function initialize_database
	file = 'DATABASE/compressed_database.bin'

	openr, 2, file, /f77

	nY = 0L
	nn0 = 0L
	nq = 0L
	ntauv = 0L
	nsig = 0L
	ni = 0L
	npca = 0L
	nlam = 0L
	
	readu, 2, nY, nN0, nq, ntauv, nsig, ni, npca, nlam

	temp = fltarr(nY)
	readu,2,temp
	temp = fltarr(nsig)
	readu,2,temp
	temp = fltarr(nn0)
	readu,2,temp
	temp = fltarr(nq)
	readu,2,temp
	temp = fltarr(ntauv)
	readu,2,temp
	temp = fltarr(ni)
	readu,2,temp

	database = {npca: npca, nlam: nlam, lambda: fltarr(nlam)}

	lambda = fltarr(nlam)
	readu,2,lambda

	database.lambda = lambda
	
	close, 2

	return, database
end


; *****************************************
; Initialize the neural network structure for the interpolation
; of the SED database
; It returns a structure with all the weights of the neural networks, the PCA
; coefficients and the normalization factors
; *****************************************
function initialize_neural_networks

; Read the number of networks to be used
	openr,2,'NETWORKS/neural_topologies.dat'
   readf,2,    Nnets
   print, 'Using ', fix(Nnets), ' neural networks'
   close,2

; Read PCA eigenvectors
   print, 'Restoring PCA eigenvectors...'
   res = file_test('NETWORKS/PCA_eigenvectors.idl')

   if (res eq 1) then begin
       restore,'NETWORKS/PCA_eigenvectors.idl'
   endif else begin
       res = dialog_message('NETWORKS/PCA_eigenvectors.idl file does not exist.'+string(10B)+$
           'You do not have a full distribution of BayesCLUMPY',/error)
       stop
   endelse

   nlambda = n_elements(evec[0,*])
   print, 'PCA vectors of length ', fix(nlambda)
	 
	net = {input_norm: fltarr(6), output_norm: fltarr(1), input_mean: fltarr(6),$
		output_mean: fltarr(1), PCAvector: fltarr(nlambda),$
		input_hidden: ptr_new(), hidden_output: ptr_new(), bias: ptr_new(), ninput: 6, $
		noutput: 1, nhidden: ptr_new(), nweights: ptr_new()}

	neural = {Nnets: Nnets, net: replicate(net, Nnets), meanSED: fltarr(nlambda), lambda: fltarr(nlambda)}

	neural.Nnets = Nnets
		
	if (res eq 1) then begin
		restore,'NETWORKS/PCA_eigenvectors.idl'
	endif else begin
		res = dialog_message('NETWORKS/PCA_eigenvectors.idl file does not exist.'+string(10B)+$
			'You do not have a full distribution of BayesCLUMPY',/error)
		stop
	endelse
	
	for i = 0, Nnets-1 do begin
		neural.net[i].PCAvector = evec[i,*]
	endfor

; Read normalization factors for all the networks
	print, 'Reading normalizations...'
	for i = 0, Nnets-1 do begin
		file = 'NETWORKS/normalization'+strtrim(string(i),2)+'.net'
		openr,2,file, error=err
		if (err eq 0) then begin
			temp = fltarr(2,6)
			readf,2,temp
			neural.net[i].input_mean = temp[0,*]
			neural.net[i].input_norm = temp[1,*]
			temp = fltarr(2,1)
			readf,2,temp
			neural.net[i].output_mean = temp[0,*]
			neural.net[i].output_norm = temp[1,*]
			close,2
		endif else begin
			res = dialog_message(file+' file does not exist.'+string(10B)+$
				'You do not have a full distribution of BayesCLUMPY',/error)
			stop
		endelse
	endfor

; Read mean SED
	print, 'Reading mean SED...'
	res = file_test('NETWORKS/meanSED.idl')
	
	if (res eq 1) then begin
		restore,'NETWORKS/meanSED.idl'
	endif else begin
		res = dialog_message('NETWORKS/meanSED.idl file does not exist.'+string(10B)+$
			'You do not have a full distribution of BayesCLUMPY',/error)
		stop
	endelse
	
	neural.meanSED = meanSED
	neural.lambda = lambda

; Read neural networks' weights
	print, 'Reading neural network weights...'
	for i = 0, Nnets-1 do begin
		file = 'NETWORKS/PCA'+strtrim(string(i),2)+'.net'
		openr,2,file, error=err
		if (err eq 0) then begin
			readf,2,ninput,nhidden,noutput,iteration
			
			neural.net[i].nhidden = ptr_new(nhidden)
			
			n = ninput*nhidden+nhidden*noutput+nhidden
			neural.net[i].nweights = ptr_new(n)

			x = dblarr(n)
			readf,2,x

			input_hidden = dblarr(ninput,nhidden)
			hidden_output = dblarr(nhidden,noutput)
			bias = dblarr(nhidden)
		
			loop = 0
			for k = 0, nhidden-1 do begin
				for l = 0, ninput-1 do begin
					input_hidden[l,k] = x[loop]
					loop = loop + 1
				endfor
			endfor
			
			for k = 0, nhidden-1 do begin
				for l = 0, noutput-1 do begin
					hidden_output[k,l] = x[loop]
					loop = loop + 1
				endfor
			endfor
			
			for k = 0, nhidden-1 do begin
				bias[k] = x[loop]
				loop = loop + 1
			endfor
				
			neural.net[i].input_hidden = ptr_new(input_hidden)
			neural.net[i].hidden_output = ptr_new(hidden_output)
			neural.net[i].bias = ptr_new(bias)
			close,2
		endif else begin
			res = dialog_message(file+' file does not exist.'+string(10B)+$
				'You do not have a full distribution of BayesCLUMPY',/error)
			stop
		endelse
	endfor

	return, neural
end

; *****************************************
; Return the SED for a given combination of parameters by using the
; interpolation neural network
; *****************************************
function neural_SED, neural, sigma, Y, N, q, tauv, angle, include_agn=include_agn, jansky=jansky,$
	out_agn = out_agn
	pars = 1.d0*[sigma, Y, N, q, tauv, angle]

	pars_normalized = pars

	PCA_coeffs = fltarr(neural.Nnets)
	output = neural.lambda * 0.d0

 	for loopnet = 0, neural.Nnets-1 do begin
 	
; Normalize input for this network
		for j = 0, neural.net[loopnet].ninput-1 do begin
			pars_normalized[j,*] = (pars[j,*]-neural.net[loopnet].input_mean[j]) / neural.net[loopnet].input_norm[j]
		endfor

; Use neural network to obtain PCA coefficients
		PCA_coeff = eval_neuralnetwork(neural.net[loopnet], pars_normalized)

; Apply inverse normalization
		PCA_coeff[0] = PCA_coeff[0] * neural.net[loopnet].output_norm + $
			neural.net[loopnet].output_mean

; Add coefficient times the eigenvector
		output = output + PCA_coeff[0] * neural.net[loopnet].PCAvector

	endfor

	factor = 1.d0
	if (keyword_set(jansky)) then begin
		c = 2.99792458d10
		factor = (neural.lambda*1.d-4 / c) / 1.d-26
	endif

; Add the mean SED and reverse the initial transformation to log10
	res = 10.d0^(output + neural.meanSED)
	
; If the AGN emission is included...
	if (keyword_set(include_agn)) then begin
		out_agn = lambda_flambda_agn(neural.lambda, 0.2784)
		res = res + out_agn
	endif
	
; Finally, transform to mJy
	return, factor*res
	
end

; *****************************************
; Return the SED for a given combination of parameters by using the
; interpolation neural network
; *****************************************
function neural_SED_derivative, neural, sigma, Y, N, q, tauv, angle, which_par, $
	include_agn=include_agn, jansky=jansky,out_agn = out_agn
	
	pars = 1.d0*[sigma, Y, N, q, tauv, angle]

	pars_normalized = pars

	PCA_coeffs = fltarr(neural.Nnets)
	output = neural.lambda * 0.d0
	output_deriv = neural.lambda * 0.d0

 	for loopnet = 0, neural.Nnets-1 do begin

; Normalize input for this network
		for j = 0, neural.net[loopnet].ninput-1 do begin
			pars_normalized[j,*] = (pars[j,*]-neural.net[loopnet].input_mean[j]) / neural.net[loopnet].input_norm[j]
		endfor

; Output
; Use neural network to obtain PCA coefficients
		PCA_coeff = eval_neuralnetwork(neural.net[loopnet], pars_normalized)

; Apply inverse normalization
		PCA_coeff[0] = PCA_coeff[0] * neural.net[loopnet].output_norm + $
			neural.net[loopnet].output_mean

; Add coefficient times the eigenvector
		output = output + PCA_coeff[0] * neural.net[loopnet].PCAvector

; Derivative
; Use neural network to obtain PCA coefficients
		PCA_coeff = eval_neuralnetwork_derivative(neural.net[loopnet], pars_normalized, which_par)

; Add coefficient times the eigenvector
		output_deriv = output_deriv + PCA_coeff[0] * neural.net[loopnet].PCAvector * $
			neural.net[loopnet].output_norm[0] / neural.net[loopnet].input_norm[which_par]

	endfor

	factor = 1.d0
	if (keyword_set(jansky)) then begin
		c = 2.99792458d10
		factor = (neural.lambda*1.d-4 / c) / 1.d-26
	endif

; Add the mean SED and reverse the initial transformation to log10
	res = 10.d0^(output + neural.meanSED) * alog(10.d0) * output_deriv
	
; If the AGN emission is included...
	if (keyword_set(include_agn)) then begin
		out_agn = lambda_flambda_agn(neural.lambda, 0.2784)
		res = res + out_agn
	endif
	
; Finally, transform to mJy
	return, factor*res
	
end

; *****************************************
; Read a set of filters
; *****************************************
function read_filters, which_filters
	info = {nlambdas: 0L, lambda: ptr_new(), transmission: ptr_new(), name: '', transmissionSED: dblarr(123),$
		central: 0.d0}
	nfilters = n_elements(which_filters)
	filters = {Nfilters: nfilters, info: replicate(info, nfilters)}

; Read the normalization constants for every filter
	
	openr,2,'FILTERS/normalizations.dat'
	readf,2,Ntot_filters
; 	print, 'Reading ', fix(Ntot_filters), ' filters'
	filter_names = strarr(Ntot_filters)
	normalizations = fltarr(Ntot_filters)
	temp = ''
	for i = 0, Ntot_filters-1 do begin
		readf,2,temp		
		res = strsplit(temp,' ',/extract)
		filter_names[i] = strsplit(res[0],"'",/extract)
		normalizations[i] = float(res[1])
	endfor
	close,2

	for i = 0, nfilters-1 do begin
		which_filters[i] = strsplit(which_filters[i],"'",/extract)
	endfor

	for i = 0, nfilters-1 do begin
		temp = 0.0
		file = 'FILTERS/'+which_filters[i]+'.res'
		filters.info[i].name = which_filters[i]
		openr,2,file
		readf,2,nlines,temp

		filters.info[i].nlambdas = nlines
		lambda = dblarr(nlines)
		trans = dblarr(nlines)		
		temp = dblarr(2,nlines)

		readf,2,temp
		close,2

; Locate normalization constant
		ind = where(filter_names eq which_filters[i])

		lambda = reform(temp[0,*])
		trans = reform(temp[1,*]) / normalizations[ind[0]]
		
; Avoid negative values for the filter transmission
		ind = where(trans lt 1.d-4)
		if (ind[0] ne -1) then begin
			trans[ind] = 0.d0
		endif
		
		filters.info[i].lambda = ptr_new(lambda)
		filters.info[i].transmission = ptr_new(trans)		

; Locate the central wavalength of the filter as a weighted mean
		filters.info[i].central = tsum((*filters.info[i].lambda),$
			(*filters.info[i].lambda)*(*filters.info[i].transmission)) /$
			tsum((*filters.info[i].lambda),(*filters.info[i].transmission))
	endfor

	return, filters
end

; *****************************************
; Return the integral of the SED times the filters for all input filters
; *****************************************
function filteredSED, neural, filters, sigma, Y, N, q, tauv, angle, z, shif, extinction_curve, include_agn=include_agn, jansky=jansky

; First return the SED
	SED = neural_SED(neural, sigma, Y, N, q, tauv, angle, include_agn=include_agn, jansky=jansky)
	SED = SED/1.d10 * 10.d0^(shif) * extinction_curve

	central_lambda = dblarr(filters.Nfilters)
	SED_filtered = dblarr(filters.Nfilters)
	
; Calculate the integral of the SED times the filter
	for i = 0, filters.Nfilters-1 do begin
		central_lambda[i] = filters.info[i].central

; Reinterpolate the SED to the wavelength axis of the filter
; This gives a better sampling since the SED is usually smooth and the filter extent
; in wavelength is much reduced than the SED
		SEDfilter = interpol(SED, neural.lambda*(1.d0+z), (*filters.info[i].lambda))

		SED_filtered[i] = tsum((*filters.info[i].lambda),SEDfilter*(*filters.info[i].transmission)) /$
			tsum((*filters.info[i].lambda),(*filters.info[i].transmission))
	endfor

	return, transpose([[central_lambda],[SED_filtered]])

end
