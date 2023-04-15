pro lampda_re
; sample=mrdfits('./work5/picked_moph.fits',1)
; name=strtrim(sample.name,2)
; moph=strtrim(sample.moph,2)
; n=n_elements(name)  

; file1='drpall-v2_4_3_tar.fits'
data=mrdfits('massive_elliptical.fits',1)
name=strtrim(data.plateifu,2)
n=n_elements(name)

epsilon=fltarr(n)
lambda=fltarr(n)

; set_plot,'ps'
; device,filename='lampda_ep.ps',/color
; device, decomposed=0
; sauron_colormap
; zz=findgen(101)/100
; cgplot,zz,0.31*sqrt(zz),yrange=[0,1],xrange=[0,1],xtitle=cgSymbol('epsilon'),ytitle=cgSymbol('lambda'),charsize=1.5,charthick=3.5,thick=3
; cgoplot,[0,0.4],[0.08,0.18],color='purple',thick=3
; cgoplot,[0.4,0.4],[0.18,0],color='purple',thick=3

for i=0,n-1 do begin
;index=where(strtrim(massage.plateifu,2) eq name[i])
split=strsplit(name[i],'-',/extract)

file='./MPL-7/MPL-7/'+split[0]+'/'+split[1]+'/manga-'+name[i]+'-MAPS-HYB10-GAU-MILESHC.fits.gz'

if file_test(file) eq 0 then continue

  
pa = (data[i].NSA_SERSIC_PHI)/180.*!pi + !pi/2.;相对于x轴的旋转角（逆时针）
re = (data[i].NSA_SERSIC_TH50)
ba = (data[i].NSA_SERSIC_BA)

epsilon[i]=1-ba

cor = mrdfits(file,1)
mflux =mrdfits(file,3)
xcor = -cor[*,*,0]
ycor = cor[*,*,1]
snr = mrdfits(file,5)
vstar = mrdfits(file,15)
vstarmask = mrdfits(file,17)
sigma = mrdfits(file,18)
sigmamask = mrdfits(file,20)
sigma_cor=mrdfits(file,21)

sigma_fin2=sigma^2-sigma_cor^2
exam=where(sigma_fin2 lt 0)
sigma_fin2[exam]=0

igd = where(snr ge 3 and sigmamask eq 0 and vstarmask eq 0)

sysv=median(vstar[igd])
covel=vstar-sysv

rot1 = [[cos(pa), -sin(pa)], [sin(pa), cos(pa)]];逆时针
rot2 = [[cos(pa), sin(pa)], [-sin(pa), cos(pa)]];顺时针旋转
cor2= [[xcor[igd]], [ycor[igd]]] # rot2
 
  
yelip = re*ba*sqrt(1-(cor2[*,0]/re)^2) * (2*(cor2[*,1] gt 0)-1);画出re这条线
;lielip = igd[where(abs(cor2[*,1] - yelip) lt .7)];选择线上的点做计算
ielip = igd[where(abs(cor2[*,1]) le abs(yelip))];选择圈内的点做计算

l_up_rin = total(mflux[ielip]*sqrt(xcor[ielip]^2+ycor[ielip]^2)*abs(covel[ielip]))
l_down_rin = total(mflux[ielip]*sqrt(xcor[ielip]^2+ycor[ielip]^2)*sqrt(covel[ielip]^2 + sigma_fin2[ielip]))

lam_rin= l_up_rin/l_down_rin  

lambda[i]=lam_rin


; if moph[i] eq 'uncer' then cgoplot,[1-ba],[lam_rin],psym=14,color='green'
; if moph[i] eq 'spira' then cgoplot,[1-ba],[lam_rin],psym=46,color='blue'
; if moph[i] eq 'ellip' then cgoplot,[1-ba],[lam_rin],psym=16,color='red'
; if moph[i] eq '--' then cgoplot,[1-ba],[lam_rin],psym=14,color='black'

  
; xx=(findgen(101)-50.)/50.
; ell1=[[re*xx],[re*ba*sqrt(1-xx^2)]]
; ell2=[[re*xx],[-re*ba*sqrt(1-xx^2)]]
; rell1=ell1#rot1
; rell2=ell2#rot1

; display_pixels, xcor[igd], ycor[igd], covel[igd], range = [-100,100], title = 'V stars',xtitle='arcsec',$
;   ytitle='arcsec',charsize = 1.5,charthick=2
; cgoplot,rell1[*,0],rell1[*,1],color='red',thick=7
; cgoplot,rell2[*,0],rell2[*,1],color='red',thick=7
;cgcolorbar,range=[-50,50],/vertical,charthick=4,position=[0.76,0.17,0.79,0.92],/right,$
;  title='km/s',tlocation='right'

; oo=''
; read,oo
endfor

mwrfits,{epsi:epsilon,lamb:lambda},'elli_lamb.fits'
; device,/close_file
; set_plot,'x'

  end