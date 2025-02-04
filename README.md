EZMIST -- A python package that allows you to download MIST/MESA isochrones directly from their website
=======================================================================================================


This small package provides a direct interface to the MIST/MESA isochrone
webpage (http://waps.cfa.harvard.edu/MIST).
It compiles the URL needed to query the website and retrives the data into a
python variable.

This package has been tested on python 2.7 and python 3.

:version: 1
:author: MF

(this package is similar to EzPadova:  https://github.com/mfouesneau/ezpadova)

Installation
------------
Install with pip

```
pip install git+https://github.com/mfouesneau/ezmist
```
(`--user` if you want to install it in your user profile)

Manual installation

download the repository and run the setup

```python setup.py install```



EXAMPLE USAGE
-------------
* Example for downloading isochrone grid  in Age, FeH and Av ranges.
```python3
>>> import ezmist
>>> age_min = 9.8
>>> age_max = 10.30
>>> delta_age = 0.01

>>> feh_min = 0.25
>>> feh_max = 0.46
>>> delta_feh = 0.01

>>> av_min = 0
>>> av_max = 1.05
>>> delta_av = 0.05

>>> #grid_age_scale => ["linear", "log10"]
>>> #grid_output_option => ["theory", "photometry"]
>>> #grid_output => eg: "UBVRIplus", "PanSTARRS" , "SDSSugriz" ...
>>> #grid_vvcrit => ['vvcrit0.0','vvcrit0.4']

>>> ezmist.get_grid_isochrones(age_min, age_max, delta_age, feh_min, feh_max, delta_feh, av_min, av_max, delta_av,
                            grid_age_scale='log10',grid_output_option='photometry', grid_output='UBVRIplus', grid_vvcrit='vvcrit0.0',nprounds=3)


```

* Basic example of downloading a sequence of isochrones, plotting, saving
```python
>>> r = ezmist.get_t_isochrones(6.0, 7.0, 0.05, FeH_value=0.0, theory_output='full')
>>> import pylab as plt
>>> plt.scatter(r['logT'], r['logL'], c=r['logA'], edgecolor='None')
>>> plt.show()
>>> r.write('myiso.fits')
```

Note: MIST isochrone metallicities are defined in terms of [Fe/H] (not Z)

* getting only one isochrone
```python
>>> r = ezmist.get_one_isochrones(age=1e7, FeH=0.0,v_div_vcrit=0.0, age_scale='linear')

* getting synthetic isochrones as pandas dataframe

```python

>>> data=ezmist.get_one_isochrone(age=0.6e9,FeH=0.19,v_div_vcrit=0.0,age_scale='linear',output_option='photometry',output='UBVRIplus').to_pandas()
