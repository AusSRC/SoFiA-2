'''
This is a small test script to test out Swig SoFiA-2 as a shared library.

The sofia C library expects a pointer to the data, FITS header info and a parameter file.
The python script therefore needs to pass:
    1. a numpy 1D array of float32 (little-endian),
    2. a constructed string of header data identical to that stored in a FITS file,
    3. the path to the user parameter file - the parameter "input.source" should be set to MEM.

The function 'dict2FITSstr()' below, shows how to construct a valid FITS header string from 
a dictionary of header information. 

The function 'extractFromFits()' shows how to prepare a multi-dimensional numpy data array:
    data3D.ravel().astype('<f4')
    
The last line shows how to call the sofia library.

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
from builtins import TypeError
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
        # delete the END key, if there
        hdr.pop("END",None)
    return hdr,data3D.ravel().astype('<f4')

def dict2FITSstr(hdr_dict):
    ''' Convert the header dict to a single string matching 
        the FITS format.
    '''
    my_str = ""
    key = ""
    val = ""
    hdrsize = 0
    for key in hdr_dict.keys():
        my_str += key[0:FITS_HEADER_KEYWORD_SIZE] + " " * ((FITS_HEADER_KEYWORD_SIZE - len(key)) if len(key) < FITS_HEADER_KEYWORD_SIZE else 0)
        my_str += "= "
        if isinstance(hdr_dict[key],bool):
            if hdr_dict[key]:
                hdr_dict[key] = "T"
            else:
                hdr_dict[key] = "F"
        if isinstance(hdr_dict[key],str) and not (key in ["SIMPLEFITS_HEADER_VALUE_SIZE","EXTEND","COMMENT","HISTORY"]):
                val = "'%s'" % hdr_dict[key]
        else:
            val = "%s" % hdr_dict[key]
        my_str += val[0:FITS_HEADER_VALUE_SIZE] + " " * ((FITS_HEADER_VALUE_SIZE - len(val)) if len(val) < FITS_HEADER_VALUE_SIZE else 0)
        hdrsize += FITS_HEADER_LINE_SIZE
    if not key == "END":
        my_str +=  "END" + " " * (FITS_HEADER_LINE_SIZE - 3)
        hdrsize += FITS_HEADER_LINE_SIZE
    
    pad = (FITS_HEADER_BLOCK_SIZE - hdrsize) if (hdrsize < FITS_HEADER_BLOCK_SIZE) else (FITS_HEADER_BLOCK_SIZE - (hdrsize % FITS_HEADER_BLOCK_SIZE))
    for i in range(pad):
        my_str += " "
    hdrsize += pad
    return my_str,hdrsize

def usage():
    print("Usage: python3 python_spark.py <path-to-FITS-file> <path-to-parameter-file>")


if __name__ == "__main__":
   
#   if len(sys.argv)<3:
#       usage()
#       sys.exit(0)
         
    # Load some test data into memory
#   fitsfile = sys.argv[1]
    fitsfile = "sofia_test_datacube.fits"
    hdr,dataPtr = extractFromFits(fitsfile)
    # Format the header info appropriately
    hdrstr,hdrsize = dict2FITSstr(hdr)
    
#   path_to_par = sys.argv[2]
    path_to_par = "sofia.par"
    parsize = len(path_to_par)
    # pass off to sofia C library
    try:
        ret = sofia.sofia_mainline(dataPtr,hdrstr,hdrsize,path_to_par,parsize)
    except RuntimeError:
        print("Caught general exception\n")
        sys.exit()
    except StopIteration:
        print("Caught Null pointer\n")
        sys.exit()
    except MemoryError:
        print("Caught ALLOC error\n")
        sys.exit()
    except IndexError:
        print("Caught index range error\n")
        sys.exit()
    except IOError:
        print("Caught file error\n")
        sys.exit()
    except OverflowError:
        print("Caught integer overflow error\n")
        sys.exit()
    except TypeError:
        print("Caught user input error\n")
        sys.exit()
    except SystemExit:
        print("Caught no sources error\n")
        sys.exit()
    
    print("\nReturned to Python caller \n")
    print("With code %d\n" % ret[0])
    for i in range(len(ret)):
        print("\n",ret[i]);
    print("\n",ret[-1].tobytes())
    print(dataPtr)
