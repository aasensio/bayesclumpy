; Generate ALMA filters for channels 7 and 9
; Channel 7 : 800-1090 micron
; Channel 9 : 420-500 micron
pro gen_alma

; Channel 7
	lambda = findgen(100) / 99.d0 * (1090.d0-800.d0) + 800.d0	
	smx = long(100/24)
	xa = findgen(smx)
	spx = 0.5*(1-cos(2.0*!pi*xa/(2.0*smx-1.0)))
	transmission = replicate(0.8d0,100)
	transmission[0:smx-1] = spx
	transmission[100-smx:*] = reverse(spx)

	openw,2,'alma_ch7.res'
	printf,2,100,1.d0
	for i = 0, 99 do printf,2,lambda[i],transmission[i]
	close,2

; Channel 9
	lambda = findgen(100) / 99.d0 * (500.d0-420.d0) + 420.d0
	smx = long(100/24)
	xa = findgen(smx)
	spx = 0.5*(1-cos(2.0*!pi*xa/(2.0*smx-1.0)))
	transmission = replicate(0.8d0,100)
	transmission[0:smx-1] = spx
	transmission[100-smx:*] = reverse(spx)

	openw,2,'alma_ch9.res'
	printf,2,100,1.d0
	for i = 0, 99 do printf,2,lambda[i],transmission[i]
	close,2

	stop
end