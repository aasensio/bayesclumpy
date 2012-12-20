; This example file shows how to read the samples from the posterior from the output of Bayesclumpy
; Just pass the root of the file used in Bayesclumpy

pro read_samples, file

; Define the parameters
	params = ['sigma','Y','N','q','tauv','angle','shift','extinction','redshift']
	nparam = n_elements(params)

; Open the file with the posterior samples (it is an ASCII file) and read the chains
	openr,2,file+'post_equal_weights.dat',error=err
	nlength = file_lines(file+'post_equal_weights.dat')
	chain = fltarr(nparam+1,nlength)
	readf,2,chain
	logposterior = reform(chain[nparam,*])
	chain = chain[0:nparam-1,*]
	close,2

; Do some plots
	!p.multi = [0,3,3]
	for i = 0, 8 do begin
		plot, chain[i,*], psym=3, ytit=params[i],charsize=1.4
	endfor
	!p.multi=0
	
end

pro test_read_samples
	read_samples, 'MARKOVCHAINS/circinus'
end