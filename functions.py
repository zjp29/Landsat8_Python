#############
# Zachery Porteous, Seamus Begley
# GSP 318 Lab # 11
# python functions 

import os
import rasterio
from rasterio import plot
from rasterio.merge import merge
import numpy
import matplotlib.pyplot as plt
import math

###### INSTRUCTIONS FOR RUNNING THESE FUNCTIONS #####

# use one of the test landsat folders that I submitted in the prior lab 

# call LS_file_parse(folder with landsat files) to get list of files

# then call MTL_Parse() to get the sun elevation and sun azimuth

# then call corrections() to radiometrically correct the DN's 

# most of the functions that we have made for calculating indices work but 
# the ones that use the sq root or multiplication are coming up with odd results. 
# I believe it has something to do with the data type

# if you want to try and run one of the functions just call it like this:
    
#    NDVI()
#    TDVI() etc...

# they do not write to files yet but that is supposedly an easy step in Rasterio

#####################

##### LIST OF GLOBAL VARIABLES #####

## bx_file = file name as string, where x is the band number
## MTL_name = metadata file name as string
## B_list = list of file names
## sun_azimuth
## sun_elevation
## bandx = radiometrically calibrated band, where x is the band number
## rad_cor_bands = list of radiometrically calibrated bands
## band_stats = dictionary of quickstats for each band

##### PARSING FUNCTIONS HERE #####
global folder_pa

def LS_file_parse(folder_path):
    """
    takes a folder path to the downloaded landsat file
    returns nothing and as a side effect assigns the first seven
    bands of the landsat image and the metadata
    txt file to separate variables

    NOTE: may add the panchromatic band later just in case user wants
    it for pan sharpening
    """
    global folder_pa
    folder_pa=folder_path
    os.chdir(folder_path)
    file_list=os.listdir(folder_path)
    # make all the variables global so they can
    # be used in other functions in the script
    global B1_file
    global B2_file
    global B3_file
    global B4_file
    global B5_file
    global B6_file
    global B7_file
    global MTL_name
    global B_list
    B_list=[]
    # assign all the bands and the MTL file to variables
    for x in file_list:
        if ".TIF" in x:
            if "B11" in x:
                file_list.remove(x) # remove bands 10 and 11 because they were
            elif "B10" in x:        # getting assigned to B1_file
                file_list.remove(x)
            elif "B1" in x:
                B1_file=x
                B_list.append(x)
            elif "B2" in x:
                B2_file=x
                B_list.append(x)
            elif "B3" in x:
                B3_file=x
                B_list.append(x)
            elif "B4" in x:
                B4_file=x
                B_list.append(x)
            elif "B5" in x:
                B5_file=x
                B_list.append(x)
            elif "B6" in x:
                B6_file=x
                B_list.append(x)
            elif "B7" in x:
                B7_file=x
                B_list.append(x)
        elif "MTL" in x:
            MTL_name=x
    return(B1_file,B2_file,B3_file,B4_file, # return all the global variables
          B5_file,B6_file,B7_file,MTL_name,B_list)

def MTL_Parse():
    """
    has no input but retrieves anything that we deem necessary from
    the MTL file and assigns it to a recognizable global variable
    as a side effect
    """
    # assign global variables for each output
    global sun_elevation
    global sun_azimuth
    needed_lines=[] # an array for the return values
    with open(MTL_name,"r") as MTL_file: # use with to open the metadata (closes upon unindent)
        needed_lines=[]
        MTL_lines=MTL_file.readlines()
        for x in MTL_lines: # iterate through the lines to find needed stuff
            if "SUN_AZIMUTH" in x:
                new=x.strip()
                needed_lines.append(new)
            elif "SUN_ELEVATION" in x:
                new=x.strip()
                needed_lines.append(new)
            elif "FILE_DATE" in x:
                new=x.strip()
                needed_lines.append(new)
        sun_azimuth=float(needed_lines[1][14:]) # get the stuff with string slices
        sun_elevation=float(needed_lines[2][16:])
        file_date=needed_lines[0][11:]
        return(sun_azimuth, sun_elevation,file_date)

##### END PARSING FUNCTIONS #####

##### RADIOMETRIC/ATMOSPHERIC CORRECTIONS #####
    
def corrections():
    """
    has no inputs but converts all seven of the stored bands to
    Top of atmosphere reflectance ==> Equation:
    
    TOA(without solar correction)=
    (.00002000(multiplicative rescaling factor)*DN's)+
    -.100000(additive rescaling factor)
    
    then TOA reflectance: ==> Equation:

    TOA = TOA/sin(sun_elevation)

    Also DOS will be in this function

    image_read = src.read(1)
    image_read_masked = numpy.ma.masked_array(image_read, mask=(image_read == 0))
    """
    global path
    path = folder_pa
    #os.mkdir(path)
    #os.chmod(path,777)
    global crs
    global width
    global height
    global affine
    global band1
    global band2 # these globals will be the numpy arrays
    global band3
    global band4
    global band5
    global band6
    global band7
    global rad_cor_bands # list of all the bands
    rad_cor_bands=[]
    band_param = rasterio.open(B_list[0],'r')
    height = band_param.height
    width = band_param.width
    affine = band_param.affine
    crs = band_param.crs
    band_param.close()
    for x in B_list: # open each band and perform the following calculation
        with rasterio.open(x,"r") as b:
            band=b.read(1)
            band_masked = numpy.ma.masked_array(band, mask=(band == 0))
            TOA_x=(.00002000*band_masked)+(-.100000)
            TOA=TOA_x/math.sin(sun_elevation)
            scaled=TOA*10000 # scale the numbers by 10000
            rad_cor_bands.append(scaled) # append the new arrays to a list
    band1=rad_cor_bands[0] # assign each band to a new variable
    band2=rad_cor_bands[1]
    band3=rad_cor_bands[2]
    band4=rad_cor_bands[3]
    band5=rad_cor_bands[4]
    band6=rad_cor_bands[5]
    band7=rad_cor_bands[6]
    #all_bands=merge(rad_cor_bands)




####################################

def QuickStats():
    """
    This function takes no arguments but outputs statistics like
    one would find in ENVI: min,max,average,std dev to a python dictionary
    
    will automatically be called when rad_calibration is called,
    but will not display until a button is pressed. **need min for DOS
    """
    global band_stats # dictionary of band values
    band_stats=[]
    for x in rad_cor_bands:
        band_min = numpy.min(x) # standard stat functions from numpy
        band_max = numpy.max(x)
        band_mean = numpy.mean(x)
        band_std_dev = numpy.std(x)
        band_stats.append({"min":band_min,"max":band_max,
                        "mean":band_mean,"std":band_std_dev}) # appending to Dict
    for x in band_stats: # print the dict
        print(x, end = "\n\n")


##### END CORRECTION FUNCTIONS #####

##### CALCULATING INDICES #####

def NDVI():
    """ 
    Bands 3 and 4 from LS_File_Parse enter the function to use in
    NDVI caulculation (NDVI = (Band3 - Band4) / (Band3 + Band4)
    code sourced from: http://www.loicdutrieux.net/pyLandsat/NDVI_calc.html -
    Credit to the creators 
    """    
    # Caulculate NDVI
    ndvi = (band3 - band4) / (band3 + band4) 
    # Plot and show graphic
    plt.imshow(ndvi)
    plt.show()
    new_NDVI = rasterio.open(path+'\\ndvi.tif','w', driver='Gtiff',
                          height=height, width=width,
                          count=1, dtype='float64', crs=crs,
                             transform = affine)
    new_NDVI.write(ndvi,1)
    new_NDVI.close()

    

def RVI(): # doesnt work with .plot but shows in ENVI
    """ 
    Bands 4 and 5 from LS_File_Parse enter the function to use in the
    Simple Ratio or Ratio Vegetation Index (RVI)
    ( Red / NIR )
    """ 
    # Calculate RVI
    rvi = (band4 / band5) 
    # Plot and show graphic
    plt.imshow(rvi)
    plt.show()
    new_rvi = rasterio.open(path+'\\rvi.tif','w', driver='Gtiff',
                          height=height, width=width,
                          count=1, dtype='float64', crs=crs,
                             transform = affine)
    new_rvi.write(rvi,1)
    new_rvi.close()

def DVI():
    """ 
    Bands 5 and 4 from LS_File_Parse enter the function
    to use in the Difference Vegetation Index (DVI)
    ( NIR - Red )
    """ 
    # Calculate DVI
    dvi = (band5 - band4) 
    # Plot and show graphic
    plt.imshow(dvi)
    plt.show()
    new_dvi = rasterio.open(path+'\\dvi.tif','w', driver='Gtiff',
                          height=height, width=width,
                          count=1, dtype='float64', crs=crs,
                             transform = affine)
    new_dvi.write(dvi,1)
    new_dvi.close()

def TDVI():# currently not working
    '''
    Transformed Difference Vegetation Index, good for urban areas where
    there is vegetation 
    '''
    tdvi = (1.5 * ((band5 - band4)/numpy.square((band5**2) + band4 + 0.5)))
    plot.show(tdvi)
    new_tdvi = rasterio.open(path+'\\tdvi.tif','w', driver='Gtiff',
                          height=height, width=width,
                          count=1, dtype='float64', crs=crs,
                             transform = affine)
    new_tdvi.write(tdvi,1)
    new_tdvi.close()

def ARVI():# currently not working
    '''
    Atmospherically Resistant Vegetation Index
    '''
    arvi = ((band5-((band4)*(band2)))/(band5 - (band4+band2)))
    plot.show(arvi)
    new_arvi = rasterio.open(path+'\\arvi.tif','w', driver='Gtiff',
                          height=height, width=width,
                          count=1, dtype='float64', crs=crs,
                             transform = affine)
    new_arvi.write(arvi,1)
    new_arvi.close()

def SAVI():
    '''
    Soil adjusted Vegetation Index --
    similar to NDVI but suppresses the effects of soil pixels,
    best used in areas with relatively
    sparse vegetation where soil is visible through canopy
    '''
    savi = (1.5 * (band5 - band4)/(band5 + band4 + 0.5))
    plot.show(savi)
    new_savi = rasterio.open(path+'\\savi.tif','w', driver='Gtiff',
                          height=height, width=width,
                          count=1, dtype='float64', crs=crs,
                             transform = affine)
    new_savi.write(tdvi,1)
    new_savi.close()
    
def OSAVI():
    '''
    Optimized soil adjusted vegetation index --
    based on SAVI, best used for sparsely covered vegetation
    where soil is visible through the
    canopy. 
    '''
    osavi = (band5 - band4)/(band5 + band4 + 0.16)
    plot.show(osavi)
    new_osavi = rasterio.open(path+'\\osavi.tif','w', driver='Gtiff',
                          height=height, width=width,
                          count=1, dtype='float64', crs=crs,
                             transform = affine)
    new_osavi.write(osavi,1)
    new_osavi.close()
    
def GOSAVI():
    """ 
    Bands 3 and 5 from LS_File_Parse enter the function to use
    in Green Optimized Soil Adjusted GOSAVI 
    ( = (NIR - Green) / (NIR + Green + 0.16) )
    """ 
    # Calculate GOSAVI
    gosavi = (band5 - band3) / (band5 + band3 + 0.16) 
    # Plot and show graphic
    plt.imshow(gosavi)
    plt.show()
    new_gosavi = rasterio.open(path+'\\gosavi.tif','w', driver='Gtiff',
                          height=height, width=width,
                          count=1, dtype='float64', crs=crs,
                             transform = affine)
    new_gosavi.write(gosavi,1)
    new_gosavi.close()

def MSAVI2(): # currently not working
    """ 
    Simpler version of MSAVI, reduces soil noise and increases
    dynamic range of the vegetation signal
    2 * NIR + 1 - math.sqrt(2 * NIR + 1)^2 - 8(NIR - Red) / 2
    """ 
    # Calculate MSAVI2
    msavi2 = (2 * band5 + 1 - numpy.square(2 * band5 + 1)**2 - (8*(band5 - band4) / 2)) 
    # Plot and show graphic
    plt.imshow(msavi2)
    plt.show()
    new_msavi2 = rasterio.open(path+'\\msavi2.tif','w', driver='Gtiff',
                          height=height, width=width,
                          count=1, dtype='float64', crs=crs,
                             transform = affine)
    new_msavi2.write(msavi2,1)
    new_msavi2.close()

def EVI():
    """ 
    Originally for MODIS improvement over NDVI by optimizing the
    vegetation signal in areas of high LAI. 
    2.5 * (NIR - Red) / (NIR + 6 * Red - 7.5 * Blue + 1)
    """ 
    # makes for use in LAI
    global evi 
    # Calculate DVI
    evi = 2.5 * ((band5 - band4) / band5 + 6 * band4 - 7.5 * band2 + 1)
    # Plot and show graphic
    plt.imshow(evi)
    plt.show()
    new_evi = rasterio.open(path+'\\evi.tif','w', driver='Gtiff',
                          height=height, width=width,
                          count=1, dtype='float64', crs=crs,
                             transform = affine)
    new_evi.write(evi,1)
    new_evi.close()

    
def LAI():
    """ 
    Used to estimate foliage cover and to forecast crop growth and yield, ideally
    bright pixel values should be masked out
    (3.618 * EVI - 0.118)
    """ 
    # Calculate LAI
    lai = (3.618 * evi - 0.118) 
    # Plot and show graphic
    plt.imshow(lai)
    plt.show()
    new_lai = rasterio.open(path+'\\lai.tif','w', driver='Gtiff',
                          height=height, width=width,
                          count=1, dtype='float64', crs=crs,
                             transform = affine)
    new_lai.write(lai,1)
    new_lai.close()
    
   
    
           
                
##### END INDICES #####



##### GRAPHICAL USER INTERFACE #####


##### END GRAPHICAL USER INTERFACE #####


