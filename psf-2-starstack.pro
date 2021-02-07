path_to_data = '../../data/'
img_data_file = 'DEEP23-I1-sexbgsub-MJysr.fits'

bigfitsname = path_to_data+img_data_file
head_img = headfits(bigfitsname)
size_x = sxpar( head_img ,'NAXIS1')
size_y = sxpar( head_img ,'NAXIS2')
image_I1 = fltarr(size_x,size_y)

for i_rows = 0LL, (size_y - 1LL - 1LL)/2LL do begin
	i_rows_image = 2LL*i_rows
	if i_rows mod 2000 eq 1 then print, (size_y - 1LL - 1LL)/2LL - 1 - i_rows, ' <===== loading fits ', systime()
	image_i_rows = mrdfits(bigfitsname,0, ROWS = [i_rows_image, i_rows_image + 1], /sil)
	image_I1[*,i_rows_image:i_rows_image+1] = image_i_rows
endfor

;IMAGE_COV_I1 = image_I1-image_I1+1

extast, head_img, ast_img

size_x_I1 = ast_img.naxis[0]
size_y_I1 = ast_img.naxis[1]

pixel_scale_I1 = pixel_scale(path_to_data+img_data_file)

readcol, 'star_list.txt', x_star, y_star, flag

index = where(flag gt 0)

forprint, x_star[index], y_star[index], textout='psf_list.txt'

readcol, 'psf_list.txt', x_psf_list, y_psf_list

half_size_stamp_cut = 35

spawn, "rm -rf stamp_img_I1"
spawn, "mkdir stamp_img_I1"

for i_psf = 0, n_elements(x_psf_list) - 1 do begin

x_psf = x_psf_list[i_psf]
y_psf = y_psf_list[i_psf]

psf_stamp = image_I1[long(x_psf) - half_size_stamp_cut:long(x_psf) + half_size_stamp_cut, $
						 long(y_psf) - half_size_stamp_cut:long(y_psf) + half_size_stamp_cut]

writefits, 'stamp_img_I1/'+strtrim(i_psf,1)+'.fits', psf_stamp;, head_stamp

head_stamp = headfits('stamp_img_I1/'+strtrim(i_psf,1)+'.fits')

FXADDPAR, head_stamp, 'x_psf', x_psf, ' x in ip_subaru image'
FXADDPAR, head_stamp, 'y_psf', y_psf, ' y in ip_subaru image'
;FXADDPAR, head_stamp, 'cover number',IMAGE_COV_I1[long(x_psf),long(y_psf)], ' cover image number'

spawn, 'rm '+'stamp_img_I1/'+strtrim(i_psf,1)+'.fits'

writefits, 'stamp_img_I1/'+strtrim(i_psf,1)+'.fits', psf_stamp, head_stamp
print,'stamp_img_I1/'+strtrim(i_psf,1)+'.fits'

endfor



for i_number = 0, n_elements(x_psf_list) - 1 do begin
	image_name_string = strtrim(i_number,1)


print, image_name_string

galaxy_stamp_image = mrdfits('stamp_img_I1/'+image_name_string+'.fits',0)

galaxy_size = size(galaxy_stamp_image)

x_galaxy_size = galaxy_size[1]
y_galaxy_size = galaxy_size[2]


psf_size = 51

galaxy_recenter = findgen(psf_size,psf_size)

scale_size = 11

galaxy_temp_2 = findgen(x_galaxy_size,y_galaxy_size*scale_size)
galaxy_temp_interpol_2 = findgen(x_galaxy_size*scale_size,y_galaxy_size*scale_size)

for i_x = 0, x_galaxy_size-1 do begin	
	galaxy_temp_2[i_x,*]=transpose(interpol(transpose(galaxy_stamp_image[i_x,*]),findgen(y_galaxy_size),findgen(y_galaxy_size*scale_size)/(1.*scale_size)))
endfor
for i_y = 0, y_galaxy_size*scale_size -1 do begin	
	galaxy_temp_interpol_2[*,i_y] = interpol(galaxy_temp_2[*,i_y],findgen(x_galaxy_size),findgen(x_galaxy_size*scale_size)/(1.*scale_size))
endfor



galaxy_temp_1 = findgen(x_galaxy_size*scale_size,y_galaxy_size)
galaxy_temp_interpol_1 = findgen(x_galaxy_size*scale_size,y_galaxy_size*scale_size)

for i_y = 0, y_galaxy_size - 1 do begin
	galaxy_temp_1[*,i_y] = interpol(galaxy_stamp_image[*,i_y],findgen(x_galaxy_size),findgen(x_galaxy_size*scale_size)/(1.*scale_size))
endfor
	
for i_x = 0, x_galaxy_size*scale_size -1 do begin
	galaxy_temp_interpol_1[i_x,*]=transpose(interpol(transpose(galaxy_temp_1[i_x,*]),findgen(x_galaxy_size),findgen(x_galaxy_size*scale_size)/(1.*scale_size)))
endfor


galaxy_temp_interpol = (galaxy_temp_interpol_2 + galaxy_temp_interpol_1)/2.

	
;print, max(galaxy_temp_interpol_1)
;print, max(galaxy_temp_interpol_2)
;print, max(galaxy_temp_interpol)

size_galaxy_interp = size(galaxy_temp_interpol)

max_galaxy_interp = max(galaxy_temp_interpol, i_max)
i_max_x = i_max mod size_galaxy_interp[1]
i_max_y = i_max / size_galaxy_interp[1]

;print, galaxy_temp_interpol[i_max_x, i_max_y]




;help, galaxy_temp_interpol[i_max_x- long(((psf_size - 1)/2. )* scale_size) - (scale_size -1 )/2:i_max_x + long(((psf_size - 1)/2. )* scale_size) + (scale_size -1 )/2,i_max_y-long(((psf_size - 1)/2. )* scale_size) -(scale_size -1 )/2:i_max_y+long(((psf_size - 1)/2. )* scale_size)+(scale_size -1 )/2]

galaxy_recenter_interp = galaxy_temp_interpol[i_max_x- long(((psf_size - 1)/2. )* scale_size)-(scale_size -1 )/2:i_max_x + long(((psf_size - 1)/2. )* scale_size)+ (scale_size -1 )/2,i_max_y-long(((psf_size - 1)/2. )* scale_size)-(scale_size -1 )/2:i_max_y+long(((psf_size - 1)/2. )* scale_size)+(scale_size -1 )/2]


for i_psf = 0, psf_size - 1 do begin
	for j_psf = 0, psf_size - 1 do begin
		galaxy_recenter[i_psf, j_psf] = total(galaxy_recenter_interp[i_psf*scale_size:(i_psf+1)*scale_size-1,j_psf*scale_size:(j_psf+1)*scale_size-1])/(scale_size^2.)
	endfor
endfor

;------------------------------------------------
;substract the background by the mod value:
binsize = 0.0005
min_hist = 0
max_hist = max(galaxy_recenter)

hist_counts = histogram( galaxy_recenter, binsize = binsize, min = min_hist, max = max_hist, locations = l_hist)
plot, l_hist, hist_counts, ps = 10, xran = [0, 0.1]

print, l_hist(where(hist_counts eq max(hist_counts))), median(galaxy_recenter)

galaxy_recenter_mod = l_hist(where(hist_counts eq max(hist_counts)))
galaxy_recenter = galaxy_recenter - galaxy_recenter_mod[0]
;------------------------------------------------


psf_image = galaxy_recenter/total(galaxy_recenter)

surface, psf_image

writefits, 'stamp_img_I1/psf_'+image_name_string+'.fits', psf_image


wait, 0.1
;stop

endfor

readcol, 'psf_list.txt', x_psf_list, y_psf_list


;==============================================================================
;build an IDL code called PSF_generate.pro based on the above infor.
;In my mac, I can run this code by 	spawn, "echo  'PSF_generate'|idl"
;In my linux, I cannot. So I run this 'PSF_generate' individually
openw, lun_code, 'PSF_generate.pro', /get_lun
printf, lun_code, "pro PSF_generate"
for i_code = 0, n_elements(x_psf_list) - 1 do begin
printf, lun_code, "psf_stamp_image_"+strtrim(i_code,1)+" = mrdfits('stamp_img_I1/psf_"+strtrim(i_code,1)+".fits',0)
print, "psf_stamp_image_"+strtrim(i_code,1)+" = mrdfits('stamp_img_I1/psf_"+strtrim(i_code,1)+".fits',0)
;stop
endfor
printf, lun_code, "psf_stack_a = ($
for i_code = 0, n_elements(x_psf_list) - 2 do begin
printf, lun_code, "psf_stamp_image_"+strtrim(i_code,1)+" + $
print, "psf_stamp_image_"+strtrim(i_code,1)+" + $
endfor
printf, lun_code, "psf_stamp_image_"+strtrim(i_code, 1)+" )"
printf, lun_code, "psf_stack_a = psf_stack_a / total(psf_stack_a)
printf, lun_code, "psf_stack_m = findgen("+strtrim(psf_size, 2)+","+strtrim(psf_size, 2)+")
printf, lun_code, "for i_x = 0, "+strtrim(psf_size, 2)+'-1 do begin
printf, lun_code, "	for i_y = 0, "+strtrim(psf_size, 2)+'-1 do begin
printf, lun_code, "		psf_stack_m[i_x,i_y] = median($
printf, lun_code, "				[$
for i_code = 0, n_elements(x_psf_list) - 2 do begin
	printf, lun_code, "psf_stamp_image_"+strtrim(i_code,1)+"[i_x,i_y],$"
endfor
printf, lun_code, "psf_stamp_image_"+strtrim(i_code, 1)+"[i_x,i_y]], /EVEN)
printf, lun_code, "	endfor"
printf, lun_code, "endfor"
printf, lun_code, "writefits, 'psf_stamp_shift_coadd_m.fits',psf_stack_m/total(psf_stack_m)
printf, lun_code, "writefits, 'psf_stamp_shift_coadd_a.fits',psf_stack_a/total(psf_stack_a)
printf, lun_code, "end
free_lun, lun_code
;==============================================================================
spawn, "echo  'PSF_generate'|idl"
;==============================================================================

;then have a look at the psf:
psf_I1 = mrdfits('psf_stamp_shift_coadd_m.fits')
gauss_psf_I1 = gauss2dfit(psf_I1,para_I1)
psf_fwhm_I1 = 0.5*(para_I1[2]+para_I1[3])*pixel_scale_I1*2.*sqrt(2.*alog(2.0))

print, psf_fwhm_I1


end
