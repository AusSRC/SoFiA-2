#!/usr/bin/env python
# coding: utf-8


import numpy as np # linear algebra
import os
import math
import copy
import random
import sys

import time

import py4j.protocol  
from py4j.protocol import Py4JJavaError  
from py4j.java_gateway import JavaObject  
from py4j.java_collections import JavaArray, JavaList

from pyspark import RDD, SparkContext  
from pyspark.serializers import PickleSerializer, AutoBatchedSerializer
from pyspark.ml.feature import StringIndexer, StandardScaler,VectorAssembler
from pyspark.ml import Pipeline

from pyspark.sql import functions as F
import pyspark.sql.functions as f

from pyspark.sql.functions import udf, col
from pyspark.ml.linalg import Vectors, VectorUDT, DenseVector

from pyspark.sql.functions import rand
from pyspark.sql.types import ArrayType, StringType, FloatType,IntegerType, DataType, DoubleType, MapType, Row
from pyspark.sql.window import Window

from pyspark.mllib.evaluation import MulticlassMetrics

from astropy.io import fits
from astropy.table import Table
from astropy.visualization import astropy_mpl_style
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord


iType=IntegerType()
dType=DoubleType()
fType=FloatType()

aType=ArrayType(fType)


from time import time, clock
class Timer:
    """
    a simple class for printing time (s) since last call
    """
    def __init__(self):
        self.t0=time()
        self.t1=clock()
        self.elapsed = 1
        self.elapsedCPU = 1
    
    def start(self):
        self.t0=time()
        self.t1=clock()
        
    def stop(self):
        t1=time()
        t2=clock()
        print("Elapsed {:2.1f}s, CPU {:2.1f}s".format(t1-self.t0, t2-self.t1))
        self.elapsed = t1-self.t0
        self.elapsedCPU = t2-self.t1

timer=Timer()



def GetDetailArrays(sqlContext, fitsFilename, raType, decType, spectraType):
    sqlStmt=("""
        select sda_detail_array, sda_detail_type
        from Sparkfits_detail_arrays
        where sda_filename='{}'
    """).format(fitsFilename)

    detailsDF=sqlContext.sql(sqlStmt).persist()
    
    raArray=np.array(detailsDF.filter(f.col("sda_detail_type") == raType).select(f.col("sda_detail_array")).collect() )
    raArray=raArray.reshape(raArray.shape[2])

    decArray=np.array(detailsDF.filter(f.col("sda_detail_type") == decType).select(f.col("sda_detail_array")).collect() )
    decArray=decArray.reshape(decArray.shape[2])

    spectraArray=np.array(detailsDF.filter(f.col("sda_detail_type") == spectraType).select(f.col("sda_detail_array")).collect() )
    spectraArray=spectraArray.reshape(spectraArray.shape[2])
    
    detailsDF.unpersist()

    return raArray, decArray, spectraArray





def GetSubCube2(sqlContext, fitsFilename, decType,spectraType, raArray, raType, decArray, spectraArray, CubeSize, cacheTempTables=False):
    
    timer.start()
    if cacheTempTables:
        print("Caching intermediate dataframes...")
        
    ## Extract our parameters
    hiRa=CubeSize[0]
    loRa=CubeSize[1]
    hiDec=CubeSize[2]
    loDec=CubeSize[3]
    loFreq=CubeSize[4]
    hiFreq=CubeSize[5]
    
    ## Because the images table is partitioned on the band numbers, we extract the band number so we can leverage 
    ## the partitining of the images table
    bands = np.where((spectraArray >= loFreq) & (spectraArray <= hiFreq )) 
    loBand = min(bands[0])
    hiBand = max(bands[0])
    
    print("Parameters extracted")
    
    def CreateNewHeader():
        return True
    
    global raSelectRange
    raSelectRange=np.where(np.logical_and(raArray >= loRa, raArray <= hiRa ))
    raSelectList=raSelectRange[0].tolist()
    
    raHeaderIndex=raArray[raSelectList[0]]
    naxis1=len(raSelectList)
    
    print("Ra select range determined")
    #
    # Now we need to create the array or actual RA values to pass into SoFiA
    #
    
    raDF=sqlContext.sql("""
    with data as 
    (
        select posexplode(sda_detail_array)
        from sparkfits_detail_arrays
        where sda_filename='{}'
        and sda_detail_type='{}'
    )
    select pos as sda_index, float(col) as sda_ra from data

    """.format(fitsFilename,raType)      )
    
    raDF=raDF.filter((raDF.sda_ra <= hiRa ) & (raDF.sda_ra >= loRa )).persist()
    raDF.createOrReplaceTempView("RA")
    if cacheTempTables:
        raDF.cache().count()
        #sqlContext.catalog.cacheTable("RA") 
    
    print("Raw RA data extracted")
    #
    # This query pivots the independent RA values we want into a single list element we
    # can include as a Spark Column
    #
    
    filteredRaDF=sqlContext.sql("""
        with rawData as
        (
            select 1 grp,
            map(

                'sda_ra', sda_ra
            ) as kv
            from RA
        )
        select 
            grp,
            collect_list(float(a.kv['sda_ra'])) as ra
        from rawData a
        group by grp
    """).persist()
    print("Filtered RA data extracted")
    
    decDF=sqlContext.sql("""
    with data as 
    (
        select posexplode(sda_detail_array)
        from sparkfits_detail_arrays
        where sda_filename='{}'
        and sda_detail_type='{}'
    )
    select pos as sda_index, float(col) as sda_declination from data

    """.format(fitsFilename, decType)      ).persist()
    print("Raw Declination data extracted")
    
    freqDF=sqlContext.sql("""
    with data as 
    (
        select posexplode(sda_detail_array)
        from sparkfits_detail_arrays
        where sda_filename='{}'
        and sda_detail_type='{}'
    )
    select pos as sda_index, float(col) as sda_Frequency_hz from data

    """.format(fitsFilename, spectraType)      ).persist()
    
    print("Raw Frequency data extracted")
    

    decFilterDF   = decDF.filter((decDF.sda_declination <= hiDec ) & (decDF.sda_declination >= loDec )).persist()
    ## old calculation ###freqFilterDF  = freqDF.filter((freqDF.sda_Frequency_hz <= hiFreq ) & (freqDF.sda_Frequency_hz >= loFreq )).persist()
    freqFilterDF  = freqDF.filter((freqDF.sda_Frequency_hz <= hiFreq ) \
                              & (freqDF.sda_Frequency_hz >= loFreq ))\
    .select('sda_index', 'sda_frequency_hz', \
              F.ntile(200).over(Window.partitionBy().orderBy(freqDF['sda_index']) ).alias("bins")
             ).persist()

    # decFilterDF.registerTempTable("DECLINATIONS")
    PositionDF=decFilterDF.crossJoin(filteredRaDF).persist()
    
    PositionDF.createOrReplaceTempView("POSITIONS")
    if cacheTempTables:
        PositionDF.cache().count()
        #spark.catalog.cacheTable("POSITIONS") 
        
    freqFilterDF.createOrReplaceTempView("FREQUENCIES")
    if cacheTempTables:
        freqFilterDF.cache().count()
        #spark.catalog.cacheTable("FREQUENCIES") 
        
    print("Filtered Dec and Frequency data extracted")

    
    decHeaderIndex=decFilterDF.groupby().max('sda_declination').collect()[0].asDict()['max(sda_declination)']
    freqHeaderIndex=freqFilterDF.groupby().min('sda_Frequency_hz').collect()[0].asDict()['min(sda_Frequency_hz)']
    naxis2=decFilterDF.count()
    naxis4=freqFilterDF.count()
    
    myImageDF=sqlContext.sql("""
        select * 
        from sparkfits_images 
        where spi_filename='{}'
        and spi_band between {} and {}
        distribute by spi_band
        sort by spi_band, spi_index
    """.format(fitsFilename, loBand, hiBand)).persist()
    myImageDF.createOrReplaceTempView("IMAGES")
    if cacheTempTables:
        print("Caching Raw image data")
        myImageDF.cache().count()
        #spark.catalog.cacheTable("IMAGES") 

    print("Raw Images data extracted")
    
    subCubeSQLDF=sqlContext.sql("""
    select /*+ BROADCAST(POSITIONS), BROADCAST(FREQUENCIES) */
    spi_index, POSITIONS.ra, POSITIONS.sda_declination, IMAGES.spi_image, 
    IMAGES.spi_filename, IMAGES.spi_band, FREQUENCIES.sda_Frequency_hz,
    FREQUENCIES.bins
    from IMAGES
        INNER JOIN FREQUENCIES
            on IMAGES.spi_band == FREQUENCIES.sda_index
        INNER JOIN POSITIONS
            on IMAGES.spi_index == POSITIONS.sda_index
    """ ).persist() #.explain() .format(hiFreq, loFreq, hiDec, loDec)
    
    # We still have to apply the Ra selectin criteria
    
    subCubeSQLDF=subCubeSQLDF\
    .withColumn("raSelectRange", f.array(  [f.col("spi_image")[i]  for i in raSelectList  ] ))\
    .select("bins", "spi_index","ra", "sda_declination","spi_image","raSelectRange","spi_filename","spi_band","sda_Frequency_hz").persist()
    
    if cacheTempTables:
        print("Caching return subcube...")
        subCubeSQLDF.cache().count()
    
        print("Final base image dataframe created...removing intermediate dataframes from cache...")

        raDF.unpersist()
        filteredRaDF.unpersist()
        decDF.unpersist()
        freqDF.unpersist()
        decFilterDF.unpersist()
        freqFilterDF.unpersist()
        PositionDF.unpersist()
        myImageDF.unpersist()

        timer.stop()
    
    return subCubeSQLDF,raHeaderIndex,decHeaderIndex,freqHeaderIndex, naxis1,naxis2, naxis4

def FlattenDataFrame(sc, subCubeDF, fitsFilename):
    subCubeDF.createOrReplaceTempView("collatedImages")
    compress1=sc.sql("""
    with rawData as
    (
        select bins, sda_Frequency_hz, spi_index, ra as rightAscension,
        map(

            'dec', sda_declination
        ) as kva,
        map(
            'pixels', raSelectRange
        ) as kvi
        from collatedImages
        distribute by sda_Frequency_hz
        sort by sda_Frequency_hz, spi_index
    )
    select 
        sda_Frequency_hz as frequency,rightAscension,
        bins,
        collect_list(float(a.kva['dec']))as declination
        ,collect_list(array(a.kvi['pixels']))as pixs
    from rawData a
    group by frequency, rightAscension, bins
    """)
    
    compress1.createOrReplaceTempView("compress_one")
    
    compress2=sc.sql("""
    select
        bins,
        rightAscension, declination,
        collect_list(array(float(a.kvi['frequencies']))) as frequencies,
        collect_list(array(a.kva['pixs'])) as pixels
    from (
        select
            bins,
            rightAscension, declination,
            map('frequencies', frequency) as kvi,
            map('pixs', pixs) as kva
        from compress_one
        distribute by rightAscension
        sort by frequency
    ) a
    group by bins,
    rightAscension, declination
    """)
    
    originalHeader=sc.sql("""
    with rawData as
    (
        select grp, sfh_index, 
        map(

            'key', sfh_key
        ) as kva,
        map(
            'value', sfh_value
        ) as kvi
        from (
        --headers
            select 1 as grp, sfh_index, sfh_key, sfh_value 
            from sparkfits_fits_headers
            where sfh_fits_file='{}' 
            order by sfh_index
        ) a
        distribute by grp
        sort by grp, sfh_index
    )
    select 
        grp,
        collect_list(string(a.kva['key']))as keys
        ,collect_list(string(a.kvi['value']))as values
    from rawData a
    group by grp
    """.format(fitsFilename)).select("keys","values")
    
    compress3 = compress2.crossJoin(originalHeader)
    
    return compress3

def GetDataframeSize(sc, df):

    # Helper function to convert python object to Java objects
    def _to_java_object_rdd(rdd):  
        """ Return a JavaRDD of Object by unpickling
        It will convert each Python object into Java object by Pyrolite, whenever the
        RDD is serialized in batch or not.
        """
        rdd = rdd._reserialize(AutoBatchedSerializer(PickleSerializer()))
        return rdd.ctx._jvm.org.apache.spark.mllib.api.python.SerDe.pythonToJava(rdd._jrdd, True)

    # First you have to convert it to an RDD 
    JavaObj = _to_java_object_rdd(df.rdd)

    # Now we can run the estimator
    dfSize=sc._jvm.org.apache.spark.util.SizeEstimator.estimate(JavaObj)* 32 / 8e6
    return dfSize


# ====================================================
# Write results to Parquet table
# ====================================================

def writeResults(resultDF, vMode, vFormat, vTable):
    running=True
    i=0
    writeTries=10
    while running and i <= writeTries:
        try:
            resultDF.write.mode(vMode).format(vFormat).saveAsTable(vTable)
            running=False
        except Exception as e:
            ## caused by multi[ple inserts as singlerows
            logger.info("WARNING HIVE INSERT FAILURE...retrying...")

            if i == writeTries:
                logger.info("ERROR - HIVE INSERT FAILURE...{}".format(str(e)))
                running=False
        finally:
            pass
        
        i+=1

# ====================================================
#  The basic UDF that creates the arrays and/or the FITS object
# ====================================================

def CreateFITSSubCubeUDF(ra, decl, freq, pixels, keys, values, fitsFilename):    
    #
    # Functions we need
    #
    def RetrieveBaselineHeader(headerDict):
        # Create header file

        hdu=fits.PrimaryHDU()
        hduHeader=hdu.header

        def isFloat(string):
            try:
                float(string)
                return True
            except ValueError:
                return False
        #for row in hduData.rdd.toLocalIterator():
        for key, value in headerDict.items():

            if value=='True':
                v=True
            elif isFloat(value):
                v=float(value)
                pass
            elif value.isnumeric():
                v=int(value)
                pass
            else:
                v=value
            #print(row.key, v) #row.value, row.value.isnumeric(), isFloat(row.value), bool(row.value))
            hduHeader[key] = (v)
            pass

        return hduHeader
    
    def UpdateNewHeader(newheader, retrievedHeader, raIdx, decIdx, freqIdx):
        def isFloat(string):
            try:
                float(string)
                return True
            except ValueError:
                return False

        for i in np.arange(len(retrievedHeader)):
            try:
                newheader[list(retrievedHeader.keys())[i]]
                # print("Exists!", list(jheader.keys())[i], jheader[int(i)],  jheader.comments[int(i)])
            except Exception as e:
                # print(list(header.keys())[i], header[int(i)],  header.comments[int(i)])
                jkey = list(retrievedHeader.keys())[i]

                if isFloat(retrievedHeader[int(i)]):
                    jval=float(retrievedHeader[int(i)])
                    pass
                elif retrievedHeader[int(i)].isnumeric():
                    jval=int(row.value)
                    pass
                else:
                    jval = retrievedHeader[int(i)]


                #jcom = retrievedHeader.comments[int(i)]

                newheader[jkey] = (jval) #, jcom)
                
        wcs=WCS(retrievedHeader)
        crpix1,crpix2,crpix3,crpix4=wcs.wcs_world2pix(raIdx,decIdx ,1, freqIdx,1, ra_dec_order=True)

        newheader['CRPIX1']=float(crpix1)
        newheader['CRPIX2']=float(crpix2)
        newheader['CRPIX3']=float(1.0)
        newheader['CRPIX4']=float(crpix4)

        return newheader
    #
    # set up the dictionary to return the log messages
    # 
    #messages = {}
    #msgNum =0
    timer.start()
    try:
        #
        # Collect the arrays from the dataframe
        #

        msg = "Commencing data extraction from dataframe..."
        #newMsg = {msgNum:msg}
        #messages.update(newMsg)
        
        # For Sofia, numpy data arrays need to be Float32 little-endians, not Float64 'float64')

        ra=np.array(ra, dtype='<f4')
        declination=np.array(decl, dtype='<f4') #.reshape(r,sequence_len)
        frequencies=np.array(freq, dtype='<f4') #.reshape(r,sequence_len)
        pixels=np.array(pixels, dtype='<f4')
        keys=np.array(keys, dtype='str')
        values=np.array(values, dtype='str')

        #
        # Reshape the pixels array to the 3 dimensional numpy array
        #
        pixels=pixels.reshape(pixels.shape[0], pixels.shape[1], pixels.shape[2], pixels.shape[4])

        msg = "Data extracted from dataframe"

        #
        # Create the corner pixel indexes for the new header
        #

        raHeaderIndex=ra[0]
        decHeaderIndex=declination[0]
        freqHeaderIndex=frequencies[0]

        msg = "New position indexes created"

        #
        # Create the new FITS object
        #
         
        createFitsObject=False
        
        if createFitsObject:

            newhdu=fits.PrimaryHDU(data=pixels)
            newheader=newhdu.header

            msg = " new base fits object created "

            #
            # Get the baseline header
            #


            c=dict(zip(keys,values))

            baselineHeader=RetrieveBaselineHeader(c)

            msg="Baseline FITS header created."

            newHeader=UpdateNewHeader(newheader, baselineHeader, raHeaderIndex, decHeaderIndex, freqHeaderIndex)

            msg = "New base FITS header created "

            #
            # Create the new header
            #

            timer.start()
            newHeader=UpdateNewHeader(newheader, baselineHeader, raHeaderIndex, decHeaderIndex, freqHeaderIndex)

            msg = "New FITS header updated - getting size of FITS object"

            sizeOfFITS=sys.getsizeof(newhdu)
        
            msg = ("All complete, no errors for {} frequencies starting at {}. FITS size {} "\
                   .format(len(frequencies), str(freqHeaderIndex), str(sizeOfFITS) ) )
            pass
        else:
            msg = ("All complete, no errors for {} frequencies starting at {}. No FITS creation "\
                   .format(len(frequencies), str(freqHeaderIndex) ) )
            
    except Exception as e:
        msg = " ERROR! "
        if hasattr(e, 'message'):
            msg += str(e.message)[0:128]
        else:
            print(e)
            msg += e

    finally:
        timer.stop()
        return msg 