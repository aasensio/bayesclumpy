; Run all the inferences for the nextfilter paper
pro nextfilter, doinference=doinference

	str = strarr(31)
	openr,2,'chain_original.cfg'
	readf,2,str
	close,2

	seds = ['ngc1566','ngc4151','sy1','ngc1068','circinus','sy2']
	sy_type = [1,1,1,2,2,2]
	prefix = ['','_below4','_above4']
	redshift = [0.005,0.0033,0.0,0.0038,0.0014,0.0]
	
	nseds = n_elements(seds)

	for i = 0, nseds-1 do begin
		for j = 0, 2 do begin
			name = seds[i]+prefix[j]

			str[12] = "'OBSERVATIONS/NEXTFILTER/"+name+".cat'"

			str[14] = "'MARKOVCHAINS/NEXTFILTER/"+name+"'"

			if (keyword_set(doinference)) then begin
; First do inference
				str[2] = '-2'

; If Sy1, then use AGN and extinction
				if (sy_type[i] eq 1) then begin
					str[4] = '9'
					str[16] = '1'
					str[18] = '6'
					str[29] = "'extinction'  0.00000   5.00000  'U'  0.00000  0.00000"
				endif else begin
					str[4] = '9'
					str[16] = '0'
					str[18] = '0'
					str[29] = "'extinction'  0.00000   5.00000  'D'  0.00000  0.00000"
				endelse

				str[30] = "'redshift'  0.00000   6.00000  'D'  "+strtrim(string(redshift[i]),2)+"  0.00000"

				openw,2,'chain.cfg'
				for k = 0, 30 do printf,2,str[k]
				close,2

				spawn,'./clumpy_mcmc'
			endif


; If Sy1, then use AGN and extinction
			if (sy_type[i] eq 1) then begin
				str[4] = '9'
				str[16] = '1'
				str[18] = '6'
				str[29] = "'extinction'  0.00000   5.00000  'U'  0.00000  0.00000"
			endif else begin
				str[4] = '9'
				str[16] = '0'
				str[18] = '0'
				str[29] = "'extinction'  0.00000   5.00000  'D'  0.00000  0.00000"
			endelse

			str[30] = "'redshift'  0.00000   6.00000  'D'  "+strtrim(string(redshift[i]),2)+"  0.00000"

; Then suggest new filter
			str[2] = '-3'
			
			openw,2,'chain.cfg'
			for k = 0, 30 do printf,2,str[k]
			close,2			
			spawn,'./clumpy_mcmc'
			
		endfor
	endfor
	stop
	
end