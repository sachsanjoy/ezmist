import ezmist
age_min = 9.8
age_max = 10.0
delta_age = 0.01

feh_min = 0.25
feh_max = 0.46
delta_feh = 0.01

av_min = 0
av_max = 1.05
delta_av = 0.05

#grid_age_scale => ["linear", "log10"]
#grid_output_option => ["theory", "photometry"]
#grid_output => eg: "UBVRIplus", "PanSTARRS" , "SDSSugriz" ...
#grid_vvcrit => ['vvcrit0.0','vvcrit0.4']

ezmist.get_grid_isochrones(age_min, age_max, delta_age, feh_min, feh_max, delta_feh, av_min, av_max, delta_av,
                            grid_age_scale='log10',grid_output_option='photometry', grid_output='UBVRIplus', grid_vvcrit='vvcrit0.4',nprounds=3)

