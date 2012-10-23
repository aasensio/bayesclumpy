pro ks_test, file1, file2
; Read the Markov chains
	nparam = 9
	openr,2,file1+'post_equal_weights.dat',error=err
	if (err ne 0) then begin
		res = dialog_message('Error opening posterior samples.'+string(10B)+$
			'You should run again the inference',/error)
		return
	endif else begin
		nlength = file_lines(file1+'post_equal_weights.dat')
		chain1 = fltarr(nparam+1,nlength)		
		readf,2,chain1
		logposterior = reform(chain1[nparam,*])
		chain1 = chain1[0:nparam-1,*]
		close,2
	endelse

	openr,2,file2+'post_equal_weights.dat',error=err
	if (err ne 0) then begin
		res = dialog_message('Error opening posterior samples.'+string(10B)+$
			'You should run again the inference',/error)
		return
	endif else begin
		nlength = file_lines(file2+'post_equal_weights.dat')
		chain2 = fltarr(nparam+1,nlength)
		readf,2,chain2
		logposterior = reform(chain2[nparam,*])
		chain2 = chain2[0:nparam-1,*]
		close,2
	endelse

	pars = ['sigma','Y','N0','q','tau','i']

	for i = 0, n_elements(pars)-1 do begin
		sample1 = reform(chain1[i,*])
		sample2 = reform(chain2[i,*])
		kstwo, sample1, sample2, D, prob		
		print, pars[i], D, prob
		
	endfor

	stop

end
