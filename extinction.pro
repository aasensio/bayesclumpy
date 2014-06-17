@bayesclumpy_routines
@read_samples
@read_seds
; This example shows how to estimate the posterior for the covering factor from the samples of the posterior
; Just pass the root of the file used in Bayesclumpy and the distance to the source in Mpc
; It also computes the total covering factor f2
pro extinction, file

	restore,'database.idl'
	lambda = database.lambda
	nlam = database.nlam
	
; Define the parameters
	params = ['sigma','Y','N','q','tauv','angle','shift','extinction','redshift']
	nparam = n_elements(params)

; Read the Markov chains
	read_samples, file, chain
	
	nlength = n_elements(chain[0,*])
			
; Calculate AGN luminosity
	N0 = chain[2,*]
	tau = chain[4,*]
	sigma = chain[0,*]
	i = chain[5,*]
	Av = 1.086*(N0*tau*exp(-(90-i)^2/sigma^2))

	h = histog(Av, nbin=40)
	plot, h[*,0], h[*,1] / max(h[*,1]), xtit='Extinction', ytit='Normalized posterior', psym=10
	perc = percentile(Av, [50.d0 - 50.d0*erf(1.d0/sqrt(2.d0)), 50.d0, 50.d0 + 50.d0*erf(1.d0/sqrt(2.d0))])
	verx, perc[0], line=1
	verx, perc[1], line=0
	verx, perc[2], line=1
	print, 'Apparent covering factor'
	print, 'Median = ', perc[1]
	print, '-1sigma = ', perc[1]-perc[0]
	print, '+1sigma = ', perc[2]-perc[1]
	
	stop

end


pro test_extinction
	extinction, 'MARKOVCHAINS/circinus'
end