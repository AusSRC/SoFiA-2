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


from astropy.io import fits

DEBUG = True

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
    hdr,dataPtr = extractFromFits(fitsfile)
    datalen = dataPtr.size
#   hdr = None
#   data = None
#    sofia.mainline(parameter_file,hdr,data)
    sofia.mainline(dataPtr.newbyteorder())
    
#   if DEBUG:
#       logging.info("Main    : before creating thread")
#   sofia = threading.Thread(target=run_SoFiA, args=(parameter_file, hdr, data))
#   if DEBUG:
#       logging.info("Main    : before running thread")
#   sofia.start()
#   if DEBUG:
#       logging.info("SOFIA thread started")
#   sofia.join()
