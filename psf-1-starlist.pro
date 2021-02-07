;load the image:

path_to_data = '../../data/'
img_data_file = 'DEEP23-I1-sexbgsub-MJysr.fits'

bigfitsname = path_to_data+img_data_file
img_head = headfits(bigfitsname)
size_x = sxpar( img_head ,'NAXIS1')
size_y = sxpar( img_head ,'NAXIS2')
image_I1 = fltarr(size_x,size_y)

for i_rows = 0LL, (size_y - 1LL - 1LL)/2LL do begin
	i_rows_image = 2LL*i_rows
	if i_rows mod 2000 eq 1 then print, (size_y - 1LL - 1LL)/2LL - 1 - i_rows, ' <===== loading fits ', systime()
	image_i_rows = mrdfits(bigfitsname,0, ROWS = [i_rows_image, i_rows_image + 1], /sil)
	image_I1[*,i_rows_image:i_rows_image+1] = image_i_rows
endfor

extast, img_head, ast_I1

size_x_I1 = ast_I1.naxis[0]
size_y_I1 = ast_I1.naxis[1]

; Pixel scale for this irac image is 0.6 arcsec
pixel_scale_I1 = pixel_scale(path_to_data+img_data_file)

; sex ../../data/DEEP23-I1-sexbgsub-MJysr.fits -c cat.sex -CATALOG_NAME DEEP23-ch1-raw.fits -MAG_ZEROPOINT 21.5816 -PHOT_APERTURES 4
; here is a rough catalog from SExtractor.


tab = mrdfits('DEEP23-I1-raw.fits',1)
;----------------------------------------------------------------
; Then we choose the point source:
; usually the point source have a different branch in mag_aper v.s. mag_auto plot
; and mag_auto v.s. FWHM_IMAGE plot.
; good point source candidates should have 
;	a sharp pixel value in the center, and few (or no) bright source nearby (e.g. within 5 or 10 FWHM) 
; 	no bad pixels
;	good coverage or exposure time
;	bright, but not satuated.
;	better have the center in one pixel to have a sharp peak in the center pixel,
;	or other requests
plot,tab.MAG_AUTO, tab.MAG_APER - tab.MAG_AUTO, ps = 3, yran = [0.3,1.5], xran = [14, 18]
stop

plot,tab.MAG_AUTO, tab.FWHM_IMAGE, ps = 3, yran = [3,5], xran = [14, 18]
ind_star = where(tab.FWHM_IMAGE gt 3.5 and tab.FWHM_IMAGE lt 4. $
and tab.MAG_AUTO gt 15.7 and tab.MAG_AUTO lt 16.5 $
and tab.MAG_APER - tab.MAG_AUTO gt 0.8 and tab.MAG_APER - tab.MAG_AUTO lt 1.0 $
and tab.CLASS_STAR gt 0.99)
help, ind_star

oplot,tab[ind_star].MAG_AUTO, tab[ind_star].FWHM_IMAGE, ps = 1, color = cgcolor('red')

stop
;----------------------------------------------------------------
;	now there is a star list.


; in case I would like to see the targets in ds9

; ra = tab[ind_star].alpha_J2000
;dec = tab[ind_star].delta_J2000
;color='cyan'
;openw, lun_reg, 'radec_I1_star.reg', /get_lun
;printf,lun_reg,'global color='+color+' font="helvetica 10 normal" select=1 edit=1 move=1 delete=1 include=1 fixed=0 source'
;for i_n = 0LL, N_ELEMENTS(ra)-1LL do begin
;	printf,lun_reg, strcompress('fk5;circle('+string(ra[i_n]) $
;			+','+string(dec[i_n])+',2.")',/remove_all)
;endfor
;free_lun, lun_reg


ad2xy, tab[ind_star].alpha_J2000, tab[ind_star].delta_J2000, ast_I1, x_star, y_star
flag = dblarr(n_elements(x_star))

IMAGE_HALF_SCALE = 45


;============================================================
;	I save the flag in time. The visual check process might be
;	break by worng input or something else. So before I continue
;	check the psf targets, I load the flag file from last check 
;	and start to visual check the rest of them.
;============================================================
delvar,x_star_load
delvar,y_star_load
delvar,flag_load
if file_search('star_list.txt') ne '' then readcol, 'star_list.txt', x_star_load, y_star_load, flag_load, format = 'f,f,f'
openw, lun, 'star_list.txt', /get_lun
;printf, lun, 'ra dec R_CMODEL mag_r mag_z mag_Y mag_NB921 flag', format = '(a)'
for i_load = 0, n_elements(x_star_load) - 1 do printf, lun, x_star_load[i_load], y_star_load[i_load], flag_load[i_load], format = '(f,f,f)'
flush, lun
;============================================================


print, 'You will have to check '+ strtrim(n_elements(x_star)-n_elements(x_star_load),2) + ' star by EYE'

window, 0, xs = 2000, ys = 500

for i_psf = 0LL+n_elements(x_star_load), n_elements(x_star) - 1LL do begin

cgimage, -alog10(image_I1[x_star[i_psf] - image_half_scale : x_star[i_psf] + image_half_scale, $
y_star[i_psf] - image_half_scale : y_star[i_psf] + image_half_scale]),/k, position = [0.,0,0.3,1];, $

surface, image_I1[x_star[i_psf] - image_half_scale : x_star[i_psf] + image_half_scale, $
y_star[i_psf] - image_half_scale : y_star[i_psf] + image_half_scale], position = [0.25,0,0.65,1], /noe, color = cgcolor('black')

surface, image_I1[x_star[i_psf] - image_half_scale /4. : x_star[i_psf] + image_half_scale/4., $
y_star[i_psf] - image_half_scale /4.: y_star[i_psf] + image_half_scale/4.], position = [0.6,0,1.,1], /noe, color = cgcolor('black')
					
wait, 0.1


ON_ERROR, 2
prompt_save = !prompt
rp = ''
read,rp,prompt=strtrim(n_elements(x_star) - i_psf, 2)+' left, press 1 to choose, 0 to delete          '
rp = strupcase(rp)
!prompt=prompt_save
if rp eq '1' then begin
	flag[i_psf] = rp
endif
if rp eq '0' then begin
	flag[i_psf] = rp
endif

if rp ne '1' and rp ne '0' then begin
	print, 'wrong answer'
	stop
endif

printf, lun, x_star[i_psf], y_star[i_psf], flag[i_psf], format = '(f,f,i)'
flush, lun

endfor

free_lun, lun




end
