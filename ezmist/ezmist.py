"""
EZMIST -- A python package that allows you to download MESA isochrones directly from the MIST
directly website

based on EZPADOVA

:version: 1.0
:author: MF
"""
from __future__ import print_function, unicode_literals, division

import sys
import os
import inspect
import time

if sys.version_info[0] > 2:
    py3k = True
    from urllib.parse import urlencode
    from urllib import request
    from urllib.request import urlopen
else:
    py3k = False
    from urllib import urlencode
    from urllib2 import urlopen

from io import StringIO, BytesIO
import zlib
import re
import json
from .simpletable import SimpleTable as Table


localpath = '/'.join(os.path.abspath(inspect.getfile(inspect.currentframe())).split('/')[:-1])

with open(localpath + '/mist.json') as f:
    _cfg = json.load(f)


# Help messages
# -------------

def file_type(filename, stream=False):
    """ Detect potential compressed file
    Returns the gz, bz2 or zip if a compression is detected, else None.
    """
    magic_dict = {"\x1f\x8b\x08": "gz",
                  "\x42\x5a\x68": "bz2",
                  "\x50\x4b\x03\x04": "zip",
                  b"\x50\x4b\x03\x04": "zip",
                  "PK\x03\x04": "zip",
                  b"PK\x03\x04": "zip",
                  }

    max_len = max(len(x) for x in magic_dict)
    if not stream:
        with open(filename) as f:
            file_start = f.read(max_len)
        for magic, filetype in magic_dict.items():
            if file_start.startswith(magic):
                return filetype
    else:
        for magic, filetype in magic_dict.items():
            if filename[:len(magic)] == magic:
                return filetype

    return None


# Build up URL request
# --------------------

def _get_url_args(**opts):
    """ Generates the query arguments given the selected options

    Parameters
    ----------
    opts: dict
        any field value

    Returns
    -------
    q: str
        string of arguments (joined by `&` char)
    """
    _opts = _cfg['defaults']
    _opts.update(**opts)

    # check None = ""
    q = []
    keys = _cfg["query_options"]

    for k in keys:
        val = _opts.get(k, "")
        if val is None:
            val = ""
        q.append("{key}={val}".format(key=k, val=val))

    return '&'.join(q)


def _extract_zip(zip_bytes):
    """ Extract the content of a zip file

    Parameters
    ----------
    zip_bytes: bytes
        string that contains the binary code

    Returns
    -------
    content:str
        ascii string contained in the zip code.
    """
    import io
    import zipfile
    fp = zipfile.ZipFile(io.BytesIO(zip_bytes))
    data = {name: fp.read(name) for name in fp.namelist()}
    if len(data) > 1:
        return data
    else:
        return data[list(data.keys())[0]]


def _query_website(q):
    """ Run the query on the website

    Parameters
    ----------
    q: str
        string of arguments (joined by `&` char)

    Returns
    -------
    r: str or bytes
        unzipped content of the query
    """

    url = _cfg["request_url"]
    print('Interrogating {0}...'.format(url))

    print('Request...', end='')
    if py3k:
        req = request.Request(url, q.encode('utf8'))
        print('done.')
        print("Reading content...", end='')
        c = urlopen(req).read().decode('utf8')
    else:
        c = urlopen(url, q).read()
    print('done.')

    try:
        fname = re.compile('<a href=".*">').findall(c)[0][9:-2]
    except Exception as e:
        print(e)
        raise RuntimeError("Something went wrong")

    furl = _cfg['download_url'] + fname

    print('Downloading data...{0}...'.format(furl), end='')
    if py3k:
        req = request.Request(furl)
        bf = urlopen(req)
    else:
        bf = urlopen(furl)
    r = bf.read()
    print("done.")


    typ = file_type(r, stream=True)
    # force format
    if (typ is None) & ('zip' in fname):
        typ = 'zip'
    if typ is not None:
        print(r[:100], type(r), bytes(r[:10]))
        print("decompressing archive (type={0})...".format(typ), end='')
        if 'zip' in typ:
            r = _extract_zip(bytes(r))
        else:
            r = zlib.decompress(bytes(r), 15 + 32)
        print("done.")

    return r


def _read_mist_iso_filecontent(data):

    """
    Reads in the isochrone file.

    Parameters
    ----------
    data: str or bytes
        content from the unzipped website query

    Returns
    -------
    t: Table
        table of the isochrones
    """
    import numpy as np

    try:
        try:
            f = data.decode('utf8').split('\n')
        except:
            f = data.split('\n')

        content = [line.split() for line in f]
        hdr = {'MIST': content[0][-1], 'MESA': content[1][-1]}
        abun = {content[3][i]:float(content[4][i]) for i in range(1,5)}
        hdr.update(**abun)
        hdr['ROT'] = float(content[4][-1])
        num_ages = int(content[6][-1])
        hdr['num_ages'] = num_ages

        #read one block for each isochrone
        iso_set = []
        counter = 0
        data = content[8:]

        # isochrone format
        for i_age in range(num_ages):

            #grab info for each isochrone
            _d = data[counter]
            num_eeps = int(_d[-2])
            num_cols = int(_d[-1])
            hdr_list = data[counter + 2][1:]
            if not py3k:
                # correcting for recfunctions not up to date for unicode dtypes
                hdr_list = [str(k) for k in hdr_list]
            formats = tuple([np.int32] + [np.float64 for i in range(num_cols - 1)])
            iso = np.zeros((num_eeps), {'names':tuple(hdr_list),'formats':tuple(formats)})

            #read through EEPs for each isochrone
            for eep in range(num_eeps):
                iso_chunk = data[3+counter+eep]
                iso[eep] = tuple(iso_chunk)

            iso_set.append(iso)

            counter += 3 + num_eeps + 2

        _data = np.lib.recfunctions.stack_arrays(iso_set, usemask=False)

        t = Table(_data, header=hdr)

        # make some aliases
        aliases = (('logL', 'log_L'),
                   ('logT', 'log_Teff'),
                   ('mass', 'star_mass'),
                   ('logg', 'log_g'))

        if 'log10_isochrone_age_yr' in t:
            aliases += (('logA', 'log10_isochrone_age_yr'),)
        else:
            aliases += (('age', 'isochrone_age_yr'),)

        for a, b in aliases:
            t.set_alias(a, b)
    except ValueError:
        buf = StringIO(data.decode('utf8'))
        t = Table(buf, dtype='dat')
        
    t.header['NAME'] = 'MIST/MESA isochrones'

    return t


# Convenient Functions
# --------------------

def simple_options(**kwargs):
    opts = _cfg['defaults']
    opts.update(kwargs)
    return opts


def get_standard_isochrone(ret_table=True, **kwargs):
    """ get the default isochrone set at a given time and [Fe/H]

    MIST standard age grid (107 ages for 5 < logAge < 10.3 in 0.05 dex steps)

    Parameters
    ----------

    ret_table: bool
        if set, return a eztable.Table object of the data

    **kwargs: other options

    Returns
    -------
    r: Table or str
        if ret_table is set, return a eztable.Table object of the data
        else return the string content of the data
    """
    opts = simple_options(
        age_type='standard',
        **kwargs)

    d = _get_url_args(**opts)

    r = _query_website(d)
    if ret_table is True:
        return _read_mist_iso_filecontent(r)
    else:
        return r

def get_one_isochrone(age, FeH,v_div_vcrit=0.0, age_scale='linear', ret_table=True, **kwargs):
    """ get one isochrone at a given time and [Fe/H]

    Parameters
    ----------

    age: float
        age of the isochrone (in yr)

    metal: float
        metalicity of the isochrone

    age_scale: str
        linear or log10 for units of age

    ret_table: bool
        if set, return a eztable.Table object of the data

    **kwargs: other options

    Returns
    -------
    r: Table or str
        if ret_table is set, return a eztable.Table object of the data
        else return the string content of the data
    """
    opts = simple_options(
        age_type='single',
        age_value=age,
        age_scale=age_scale,
        FeH_value=FeH,
        v_div_vcrit='vvcrit'+str(v_div_vcrit),
        **kwargs)

    d = _get_url_args(**opts)
#    print(d)
    r = _query_website(d)
    if ret_table is True:
        return _read_mist_iso_filecontent(r)
    else:
        return r


def get_t_isochrones(logt0, logt1, dlogt, age_scale='log10', ret_table=True, **kwargs):
    """ get a sequence of isochrones at constant Z

    Parameters
    ----------
    logt0: float
        minimal value of log(t/yr)

    logt1: float
        maximal value of log(t/yr)

    dlogt: float
        step in log(t/yr) for the sequence

    ret_table: bool
        if set, return a eztable.Table object of the data

    Returns
    -------
    r: Table or str
        if ret_table is set, return a eztable.Table object of the data
        else return the string content of the data
    """
    opts = simple_options(
        age_type='range',
        age_range_low=logt0,
        age_range_high=logt1,
        age_range_delta=dlogt,
        age_scale=age_scale,
        **kwargs)

    d = _get_url_args(**opts)

    r = _query_website(d)
    if ret_table is True:
        return _read_mist_iso_filecontent(r)
    else:
        return r

#        Download grid MIST
#--------------------------------
# only for python3
#--------------------------------
def grid_query_website(q):
    """ Run the query on the website

    Parameters
    ----------
    q: str
        string of arguments (joined by `&` char)

    Returns
    -------
    r: str or bytes
       passes the downloaded Filename
    """

    url = _cfg["request_url"]
    print('Interrogating {0}...'.format(url))

    print('Request...', end='')
    if py3k:
        req = request.Request(url, q.encode('utf8'))
        print('done.')
        print("Reading content...", end='')
        c = urlopen(req).read().decode('utf8')
    else:
        c = urlopen(url, q).read()

    try:
        fname = re.compile('<a href=".*">').findall(c)[0][9:-2]
    except Exception as e:
        print(e)
        raise RuntimeError("Something went wrong")

    furl = _cfg['download_url'] + fname

    print('Downloading data from MIST...{0}...'.format(furl), end='')
    if py3k:
        import requests
        rurl = requests.get(furl, allow_redirects=True)
        open(fname[4:], 'wb').write(rurl.content)
        print(fname[4:], ' File downloaded!')
        dwfname = fname[4:]
        return dwfname

def get_age_isochrones(logt0, logt1, dlogt, age_scale='log10', ret_table=True, **kwargs):
    """ get a sequence of isochrones at constant Z

    Parameters
    ----------
    logt0: float
        minimal value of log(t/yr)

    logt1: float
        maximal value of log(t/yr)

    dlogt: float
        step in log(t/yr) for the sequence

    ret_table: bool
        if set, return a eztable.Table object of the data

    Returns
    -------
    r: Table or str
        return the downloaded filename
    """
    opts = simple_options(
        age_type='range',
        age_range_low=logt0,
        age_range_high=logt1,
        age_range_delta=dlogt,
        age_scale=age_scale,
        **kwargs)

    d = _get_url_args(**opts)
    dwfname = grid_query_website(d)
    return dwfname

def extractallfn(zip_name,feh,av,vc):
    """ Extract the zipfile to make the grid and rename the files based on the parameters

    Parameters
    ----------
    zip_name: str
        name of the zipfile

    feh: float
        metalicity Fe/H

    av: float
        Av

    vc: str
        vvcritic
    """

    from zipfile import ZipFile
    import numpy as np
    print(zip_name)
    # Create a ZipFile Object and load sample.zip in it
    with ZipFile(zip_name, 'r') as zipObj:
        # Extract all the contents of zip file in different directory
        zipObj.extractall('.')
        print('mv '+zip_name[:-3]+'cmd '+'MIST_iso_feh'+str('%.2f'%feh)+'_av'+str('%.2f'%av)+'_'+vc+'.iso.cmd')
    os.system('mv '+zip_name[:-3]+'cmd '+'MIST_iso_feh'+str('%.2f'%feh)+'_av'+str('%.2f'%av)+'_'+vc+'.iso.cmd')
    os.system('rm -f '+zip_name)

def get_grid_isochrone(age_min,age_max,delta_age,feh_min,feh_max,delta_feh,av_min,av_max,delta_av,age_scale = 'log10',vc = 'vvcrit0.0',output_option = 'photometry',output = 'UBVRIplus',nprounds=3):
    """ get a grid of isochrones at varying age, FeH and Av

    Parameters
    ----------
    age_min: float
        minimal value of log(t/yr)

    age_max: float
        maximal value of log(t/yr)

    delta_age: float
        step in log(t/yr) for the sequence

    feh_min: float
        minimal value of metalicity

    feh_max: float
        maximal value of metalicity

    delta_feh: float
        step in feh for the sequence

    av_min: float
        minimal value of Av

    av_max: float
        maximal value of Av

    delta_Av: float
        step in Av for the sequence

    age_scale: str
        scaling fn of age
    
    vc : str
        critical velocity
    
    output_option: str
        default is 'photometry'
    
    output: str
       default is 'UBVRIplus'    

    nprounds: int  
        step precision default is 3
    """

    import numpy as np
    grid_age = np.round(np.arange(age_min,age_max,delta_age),nprounds) #age grid
    grid_FeH = np.round(np.arange(feh_min,feh_max,delta_feh),nprounds)  #metalicity grid
    grid_Av = np.round(np.arange(av_min,av_max,delta_av),nprounds)   #extintion grid
    
    print('Total number of grids to download : ', len(grid_Av)*len(grid_FeH))
    print('Total number of models to download : ', len(grid_Av)*len(grid_FeH)*len(grid_age))
    print(grid_Av.tolist())
    print(grid_FeH.tolist())
    print(grid_age.tolist())
    # grid loop
    for i in grid_Av: #fe
        for j in grid_FeH: #av
            dfname = get_age_isochrones(age_min, age_max, delta_age, age_scale = age_scale, ret_table=False,output_option=output_option,output=output,FeH_value=j, Av_value=i,v_div_vcrit=vc)
            extractallfn(dfname,j,i,vc)
            print('Current parameters FeH, Av : ', j , i)