pro plot_expected, file, window_left, window_right

	loadct,0,/silent
	t = ''
	nfilt = file_lines(file+'_information_gain.dat')
	filt_name = strarr(nfilt)
	lambda = dblarr(nfilt)
	information = dblarr(nfilt)
	openr,2,file+'_information_gain.dat'
	
	for i = 0, nfilt-1 do begin
		readf,2,t
		res = strsplit(t,' ',/extract)
		filt_name[i] = res[0]
		lambda[i] = double(res[1])
		information[i] = double(res[2])
	endfor

	wset, window_left
	plot, lambda, information, psym=8, xtit='Wavelength [!7l!6m]', ytit='Expected information', /xlog

	wset, window_right
	plot, findgen(nfilt+1), information, /nodata, xtit='Filter', ytit='Expected information', xsty=1, XTICKFORMAT='(A1)'
	colors = findgen(nfilt) / (nfilt-1.d0) * (255-70) + 70
	loadct,6,/silent
	for i = 0, nfilt-1 do begin	
		polyfill,[i,i,i+1,i+1],[0,information[i],information[i],0],col=colors[i]
		cwpal
		xyouts, i+1-0.25, 4, filt_name[i], orientation=90, charsize=1.1, col=5
		loadct,6,/silent
	endfor

	loadct,0,/silent
	
	close,2

end