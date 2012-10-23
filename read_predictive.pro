pro read_predictive, file
	openr,2,file
	readf,2,n
	nlines = file_lines(file)
	nfilter = fix(nlines / (n+1.d0))
	loop = 0
	c = ''
	pred = fltarr(nfilter,2,n)
	while (loop lt nfilter) do begin
		temp = fltarr(2,n)
		readf,2,temp
		pred[loop,*,*] = temp
		if (loop ne nfilter-1) then begin
			readf,2,c
		endif
		loop = loop + 1
	endwhile
	close,2
	!p.multi=[0,8,8]
	for i=0,nfilter-1 do plot,pred[i,0,*],pred[i,1,*],/xlog
	!p.multi=0
	stop
end
