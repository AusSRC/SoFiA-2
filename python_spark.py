'''
This is a small test script to test out Swig SoFiA-2

Created on 17Dec 2020
@author: ger063
'''
#!/usr/bin/env python

import os
import sys
import logging
import time
import threading
import sofia

import signal
PID = os.getpid()


from astropy.io import fits

DEBUG = True
FITS_HEADER_BLOCK_SIZE = 2880
FITS_HEADER_LINE_SIZE  =   80
FITS_HEADER_LINES      =   36
FITS_HEADER_KEYWORD_SIZE =  8
FITS_HEADER_KEY_SIZE     = 10
FITS_HEADER_VALUE_SIZE   = 70
FITS_HEADER_FIXED_WIDTH  = 20

def extractFromFits(fitsfile):
    ''' Return the FITS hdr as a dict and the FITS data
        as a flattened array of floats.
    '''
    hdr = {}
    data3D =[[[]]]
    with fits.open(fitsfile) as hdul:
        for key in hdul[0].header:
            hdr[key] = hdul[0].header[key]
            data3D = hdul[0].data
        hdul.info()
    return hdr,data3D.ravel()

def dict2FITSstr(hdr_dict):
    ''' Convert the header dict to a single string matching 
        the FITS format.
    '''
    str = ""
    key = ""
    hdrsize = 0
    for key in hdr_dict.keys():
        str += key[0:8] + " " * ((FITS_HEADER_KEYWORD_SIZE - len(key)) if len(key) < FITS_HEADER_KEYWORD_SIZE else 0)
        str += "= "
        if isinstance(hdr_dict[key],bool):
            if hdr_dict[key]:
                hdr_dict[key] = "TRUE"
            else:
                hdr_dict[key] = "FALSE"
        val = "%s" % hdr_dict[key]
        str += val[0:70] + " " * ((FITS_HEADER_VALUE_SIZE - len(val)) if len(val) < FITS_HEADER_VALUE_SIZE else 0)
        hdrsize += 80
    if not key == "END":
        str +=  "END" + " " * (FITS_HEADER_LINE_SIZE - 3)
        hdrsize += 80
    
    pad = (FITS_HEADER_BLOCK_SIZE - hdrsize) if (hdrsize < FITS_HEADER_BLOCK_SIZE) else (FITS_HEADER_BLOCK_SIZE - (hdrsize % FITS_HEADER_BLOCK_SIZE))
    for i in range(pad):
        str += " "
    hdrsize += pad
    return str,hdrsize

def run_SoFiA(par,hdr,data):
    ''' Run SoFiA2 shared lib '''
    if DEBUG:
        logging.info("Thread SOFIA: starting")
    sofia.mainline(par,hdr,data)
    if DEBUG:
        logging.info("Thread SOFIA: finishing")
    

if __name__ == "__main__":
    
    fitsfile = sys.argv[1]
    parameter_file = sys.argv[2]
    # Load some data into memory
    hdr,dataPtr = extractFromFits(fitsfile)
    hdrstr,hdrsize = dict2FITSstr(hdr)
    datalen = dataPtr.size
    path_to_par = "sofia.par"
    parsize = len(path_to_par)
    # pass data to sofia C library - we need to swap the byte order
    sofia.sofia_mainline(dataPtr,hdrstr,hdrsize,path_to_par,parsize)
    
