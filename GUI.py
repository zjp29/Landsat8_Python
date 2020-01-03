
######################
# Zachery Porteous - Seamus Begley
# GSP 318 Final Project

# Professor Graham
# 4/20/2019
######################
# rasterio and os imports
import sys
import os
import rasterio
from rasterio import plot
from rasterio.merge import merge
import numpy
import matplotlib.pyplot as plt
import math

# import the main classes and widgets from the Qt library
from PyQt5 import *
from PyQt5.QtWidgets import *
from PyQt5.QtGui import QPalette
from PyQt5.QtWidgets import QApplication, QPushButton

########################################################################################
# Main Dialog box class #2
########################################################################################
class Main(QWidget):
    """
    A GUI for rasterio Schtuff
    """
    def LS_file_parse(self):
        try: # Try except loop keeps the functions from crashing if there is an error
            """
            takes a folder path to the downloaded landsat file
            returns nothing and as a side effect assigns the first seven
            bands of the landsat image and the metadata
            txt file to separate variables
        
            NOTE: may add the panchromatic band later just in case user wants
            it for pan sharpening
            """
            global folder_path
            folder_path = self.TheLineEdit1.text()
            print(folder_path)
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
            self.Result1.setText(format(B1_file + '\n' + B2_file + '\n' +
                                        B3_file + '\n' + B4_file + '\n'
                                        + B5_file + '\n' + B6_file +'\n'
                                        + B7_file + '\n' + MTL_name))
            #return(B1_file,B2_file,B3_file,B4_file, # return all the global variables
              #    B5_file,B6_file,B7_file,MTL_name,B_list)
        except Exception as TheException: 
            self.ErrorReturn.setText("Error parsing, check input file and try again") # if an error occurs set ErrorReturn to this message
            
    def MTL_Parse(self):
        try:
            """
            has no input but retrieves anything that we deem necessary from
            the MTL file and assigns it to a recognizable global variable
            as a side effect
            """
            # assign global variables for each output
            global sun_elevation
            global sun_azimuth
            needed_lines=[] # an array for the return values
            with open(MTL_name,"r") as MTL_file: # use with to open the metadata
                needed_lines=[]
                MTL_lines=MTL_file.readlines()
                for x in MTL_lines: # iterate through the lines to find needed stuff
                    if "SUN_AZIMUTH" in x:
                        new=x.strip()
                        needed_lines.append(new)
                    elif "SUN_ELEVATION" in x:
                        new=x.strip()
                        needed_lines.append(new)
                sun_azimuth=float(needed_lines[0][14:]) # get the stuff with string slices
                sun_elevation=float(needed_lines[1][16:])
                self.Result2.setText(format('Sun Elevation: ' + str(sun_elevation) + '\n\n' + 'Sun Azimuth: ' 
                                            + str(sun_azimuth)+'\n'))
                #return(sun_azimuth, sun_elevation)   
        except Exception as TheException:
            self.ErrorReturn.setText("Error in the metadata parse operation - "+str(TheException)) # Same as above but add the error message with it          
            
    def corrections(self):
        try:
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
            path = folder_path  
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
            band_param = rasterio.open(B_list[0],'r') # opens a band to fetch the values needed to write to file
            height = band_param.height # height of the array
            width = band_param.width # width of the array
            affine = band_param.affine # geometric transformation used in writing
            crs = band_param.crs # gets the spatial reference of the band
            band_param.close()     # close the file   
            for x in B_list: # open each band and perform the following calculation
                with rasterio.open(x,"r") as b:
                    band=b.read(1)
                    band_masked = numpy.ma.masked_array(band, mask=(band == 0)) # masks the 0 values 
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
        except Exception as TheException:            
            self.ErrorReturn.setText("Issue with the corrections Input - "+str(TheException))
            
    def RVI(self):
        try:
            """ 
            Simple Ratio or Ratio Vegetation Index (RVI)
            Input: Bands 4 and 5 from LS_File_Parse
            Operation:( Red / NIR )
            Use: ratio of the wavelength with the highest reflectance for vegetation and the wavelength of the 
            deepest chlorophyll absorption
            """ 
            # Calculate RVI
            rvi = ( band5/band4 ) 
            # Plot and show graphic
            plt.imshow(rvi)
            plt.show()
            new_rvi = rasterio.open(path+'\\rvi.tif','w', driver='Gtiff',
                                      height=height, width=width,
                                      count=1, dtype='float64', crs=crs,
                                         transform = affine)
            new_rvi.write(rvi,1)
            new_rvi.close()        
        except Exception as TheException:            
            self.ErrorReturn.setText("Check operation and try again - "+str(TheException))
            
    def DVI(self):
        try:
            """ 
            Difference Vegetation Index (DVI)
            Input: Bands 5 and 4 from LS_File_Parse
            Operation: ( NIR - Red )
            Use: Distinguishes between soil and vegetation, does not account for the difference between 
            reflectance and radiance caused by atmospheric effects or shadows
            """ 
            # Calculate DVI
            dvi = (band5 - band4) 
            # Plot and show graphic
            plt.imshow(dvi)
            plt.show()
            # create the current array to .tif data using rasterio 
            new_dvi = rasterio.open(path+'\\dvi.tif','w', driver='Gtiff',
                                  height=height, width=width,
                                  count=1, dtype='float64', crs=crs,
                                     transform = affine)
            new_dvi.write(dvi,1) # write the tif data to a file
            new_dvi.close() # close the file            
        except Exception as TheException:
            self.ErrorReturn.setText("Check operation and try again - "+str(TheException))
            
    def TDVI(self):
        try:
            '''
            Transformed Difference Vegetation Index
            Input: Bands 4 and 5 from LS_File_Parse
            Operation: 1.5 * [(NIR-Red) / numpy.square((NIR**2) + Red +0.5)]
            Use: good for urban areas where there is vegetation 
            '''
            tdvi = (1.5 * ( (band5 - band4)/(numpy.sqrt((band5**2) + band4 + 0.5))))
            plot.show(tdvi)
            new_tdvi = rasterio.open(path+'\\tdvi.tif','w', driver='Gtiff',
                                  height=height, width=width,
                                  count=1, dtype='float64', crs=crs,
                                     transform = affine)
            new_tdvi.write(tdvi,1)
            new_tdvi.close()        
        except Exception as TheException:
            self.ErrorReturn.setText("Check operation and try again - "+str(TheException))
            
    def ARVI(self):
        try:
            '''
            Atmospherically Resistant Vegetation Index
            Input: Bands 4, 5, and 2 from LS_File_Parse
            Operation: (NIR-R*B)/(NIR+R*B)
            Use: good for urban areas where there is vegetation 
            '''
            arvi = (band5-(band4*band2))/(band5 + (band4*band2))
            plot.show(arvi)
            new_arvi = rasterio.open(path+'\\arvi.tif','w', driver='Gtiff',
                                  height=height, width=width,
                                  count=1, dtype='float64', crs=crs,
                                     transform = affine)
            new_arvi.write(arvi,1)
            new_arvi.close()        
        except Exception as TheException:
            self.ErrorReturn.setText("Check operation and try again - "+str(TheException))
            
    def SAVI(self):
        try:    
            '''
            Soil adjusted Vegetation Index
            Input: Band 4 and 5 from LS_File_Parse
            Operation: (1.5 * (NIR - Red)/(NIR + Red + 0.5))
            Use: similar to NDVI but suppresses the effects of soil pixels,
            best used in areas with relatively sparse vegetation where soil is visible through canopy
            '''
            savi = (1.5 * (band5 - band4)/(band5 + band4 + 0.5))
            plot.show(savi)
            new_savi = rasterio.open(path+'\\savi.tif','w', driver='Gtiff',
                                  height=height, width=width,
                                  count=1, dtype='float64', crs=crs,
                                     transform = affine)
            new_savi.write(savi,1)
            new_savi.close()        
        except Exception as TheException:
            self.ErrorReturn.setText("Check operation and try again - "+str(TheException))
            
    def OSAVI(self):
        try:
            '''
            Optimized soil adjusted vegetation index
            Input: Bands 4 and 5 from LS_File_Parse
            Operation: (NIR - Red)/(NIR + Red + 0.16)
            Use: based on SAVI, best used for sparsely covered vegetation
            where soil is visible through the canopy. 
            '''
            osavi = (band5 - band4)/(band5 + band4 + 0.16)
            plot.show(osavi)
            new_osavi = rasterio.open(path+'\\osavi.tif','w', driver='Gtiff',
                                  height=height, width=width,
                                  count=1, dtype='float64', crs=crs,
                                     transform = affine)
            new_osavi.write(osavi,1)
            new_osavi.close()        
        except Exception as TheException:
            self.ErrorReturn.setText("Check operation and try again - "+str(TheException))
            
    def GOSAVI(self):
        try:
            """ 
            Green Optimized Soil Adjusted GOSAVI 
            Input: Bands 3 and 5 from LS_File_Parse
            Operation: (NIR - Green) / (NIR + Green + 0.16)
            Use: Originally designed for color-infrared photography to predict nitrogen requirements for corn
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
        except Exception as TheException:
            self.ErrorReturn.setText("Check operation and try again - "+str(TheException))
            
    def MSAVI2(self):
        try:
            """ 
            Modified Secondary Soil-Adjusted Vegetation Index
            Input: Bands 4 and 5 from LS_File_Parse
            Operation: 2 * NIR + 1 - math.sqrt(2 * NIR + 1)^2 - 8(NIR - Red) / 2
            Use: Simpler version of MSAVI, reduces soil noise and increases
            dynamic range of the vegetation signal
            
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
        except Exception as TheException:
            self.ErrorReturn.setText("Check operation and try again - "+str(TheException))
            
    def EVI(self):
        try:
            """ 
            Enhanced Vegitation Index
            Input: Bands 2, 4, and 5 from LS_File_Parse
            Operation: 2.5 * (NIR - Red) / (NIR + 6 * Red - 7.5 * Blue + 1)
            Use: Originally for MODIS improvement over NDVI by optimizing the
            vegetation signal in areas of high LAI. 
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
        except Exception as TheException:
            self.ErrorReturn.setText("Check operation and try again - "+str(TheException))
            
    def LAI(self):
        try:        
            """ 
            Leaf Area Index
            Input: evi from EVI()
            Operation: (3.618 * EVI - 0.118)
            Use: to estimate foliage cover and to forecast crop growth and yield, ideally
            bright pixel values should be masked out
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
        except Exception as TheException:
            self.ErrorReturn.setText("Check operation and try again - "+str(TheException))
            
    def NDVI(self):
        try:
            """ 
            Normalized difference vegetation index
            Input: Bands 3 and 4 from LS_File_Parse
            Operation: (Band3 - Band4) / (Band3 + Band4)
            Use: Good for assessing living (green) vegetation, overall health of vegetation
            ####### code sourced from: http://www.loicdutrieux.net/pyLandsat/NDVI_calc.html -
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
            
        except Exception as TheException:
            self.ErrorReturn.setText("Check operation and try again - "+str(TheException))
        
    def CheckChecker(self):
        """
        Checks to see if each respective checkbox is equal to true
        Checkboxes are set to False initally and then set True when they are checked in the GUI
        If a checkbox is true, the function returns a call for the coresponding indice
        """
        
        if self.TheCheckBox_NDVI.isChecked() == True:
            self.NDVI()
        if self.TheCheckBox_TDVI.isChecked() == True:
            self.TDVI()
        if self.TheCheckBox_RVI.isChecked() == True:
            self.RVI()
        if self.TheCheckBox_DVI.isChecked() == True:
            self.DVI()
        if self.TheCheckBox_ARVI.isChecked() == True:
            self.ARVI()
        if self.TheCheckBox_SAVI.isChecked() == True:
            self.SAVI()
        if self.TheCheckBox_OSAVI.isChecked() == True:
            self.OSAVI()
        if self.TheCheckBox_GOSAVI.isChecked() == True:
            self.GOSAVI() 
        if self.TheCheckBox_MSAVI2.isChecked() == True:
            self.MSAVI2()
        if self.TheCheckBox_EVI.isChecked() == True:
            self.EVI()
            
        else:
            print("Nothing Selected")   
        
    def __init__(self):
        super().__init__()
            
        # Create a main layout to organize the rows and columns
        MainLayout = QVBoxLayout()
        self.setLayout(MainLayout)  # add the layout to the main dialog object
        
        ################################################### layouts
        
        #!!!!!!!!!!!!!! If you add Layouts EXPLICITLY say whats inside them !!!!!!!
        
        HorizontalLayoutstart=QHBoxLayout() # the input folder prompt, folder line edit box, and fetch button
        MainLayout.addLayout(HorizontalLayoutstart) 
        
        HorizontalLayoutstart2=QHBoxLayout() # Files Found button and LS_file_parse() results
        MainLayout.addLayout(HorizontalLayoutstart2) 
        
        HorizontalLayoutstart3=QHBoxLayout() # sun elevation label prompt
        MainLayout.addLayout(HorizontalLayoutstart3) 
        
        HorizontalLayoutstart4=QHBoxLayout() # get SE and SA button
        MainLayout.addLayout(HorizontalLayoutstart4) 
        
        HorizontalLayoutstart5=QHBoxLayout() # SE button results
        MainLayout.addLayout(HorizontalLayoutstart5)
        
        HorizontalLayoutstart5a=QHBoxLayout() # corrections() button
        MainLayout.addLayout(HorizontalLayoutstart5a)        
        
        HorizontalLayoutstart6=QHBoxLayout()# VI prompt
        MainLayout.addLayout(HorizontalLayoutstart6)        
        
        HorizontalLayout1=QHBoxLayout() # VI row #1
        MainLayout.addLayout(HorizontalLayout1)
        
        HorizontalLayout2=QHBoxLayout() # VI row # 2
        MainLayout.addLayout(HorizontalLayout2) 
        
        HorizontalLayout3=QHBoxLayout() # calculate indices button
        MainLayout.addLayout(HorizontalLayout3)         
                 
        HorrizontalErrorMessage=QHBoxLayout() # For displaying errors in INDICES (maybe for everything soon?)
        MainLayout.addLayout(HorrizontalErrorMessage)
        
        
        ############ Start Widgets ############################## folder edit box, LS_file_parse()
        
        # add a label for line edit   
        TheLabel = QLabel("Decompressed LS 8 Data Folder:")
        HorizontalLayoutstart.addWidget(TheLabel) # add the label
        
        # add a line edit box for folder
        self.TheLineEdit1 = QLineEdit(self)
        HorizontalLayoutstart.addWidget(self.TheLineEdit1)        
        
        # add a push button for the next window
        TheButton1=QPushButton("Fetch Files") # create the button and give it a text label
        HorizontalLayoutstart.addWidget(TheButton1) # add the button to the layout
        TheButton1.clicked.connect(self.LS_file_parse)
        
        # add a label for the result of LS_file_parse() 
        TheLabel = QLabel("Files Found:")
        HorizontalLayoutstart2.addWidget(TheLabel) # add the label 
        
        # LS_file_parse() results
        self.Result1=QLabel("")
        HorizontalLayoutstart2.addWidget(self.Result1) # will be the results of LS_file_parse() 
        
        ########################################## MTL_Parse() function button
        
        # add a label for line edit
        TheLabel_MTL = QLabel("The Sun Elevation is required for Radiometric Calibration...")
        HorizontalLayoutstart3.addWidget(TheLabel_MTL) # add the label
        
        # add a push button for the next window
        TheButton_MTL=QPushButton("Get Elevation and Azimuth") # create the button and give it a text label
        HorizontalLayoutstart4.addWidget(TheButton_MTL) # add the button to the layout
        TheButton_MTL.clicked.connect(self.MTL_Parse)
        
        self.Result2=QLabel("")
        HorizontalLayoutstart5.addWidget(self.Result2)
        
        # text before VI choices
        VI_label = QLabel("Please choose the Vegetation Indices you would like to calculate...")
        HorizontalLayoutstart6.addWidget(VI_label)
        
        ########################################## corrections() button
        
        # add a push button for the next window
        TheButton_corr=QPushButton("Radiometric Correction") # create the button and give it a text label
        HorizontalLayoutstart5a.addWidget(TheButton_corr) # add the button to the layout
        TheButton_corr.clicked.connect(self.corrections)        
        
        ########################################## Calculate Indices button
        
        TheButton_CI=QPushButton("Calculate Selected Indices") # create the button and give it a text label
        label45 = QLabel("tifs will be outputted to folder entered above...")
        HorizontalLayout3.addWidget(TheButton_CI) # add the button to the layout  
        HorizontalLayout3.addWidget(label45) 
        
        ############################################## Error message for indices (will pop-up if necessary)
        
        self.ErrorReturn = QLabel() # Returns from where it was called in an INDICIE
        HorrizontalErrorMessage.addWidget(self.ErrorReturn)
        
        #######################################################################add widgets here ROW #1
        
        NDVI_Label = QLabel("NDVI-->")
        HorizontalLayout1.addWidget(NDVI_Label)

        # add a check box
        self.TheCheckBox_NDVI = QCheckBox('') # creates a checkbox object 
        self.TheCheckBox_NDVI.setChecked(False) # sets the current value of setCheck() to false
        HorizontalLayout1.addWidget(self.TheCheckBox_NDVI) # add the checkbox to the layout
        
        TheButton_CI.clicked.connect(self.CheckChecker) # When the button is pushed, run checkchecker()
        
        TDVI_Label = QLabel("TDVI-->")
        HorizontalLayout1.addWidget(TDVI_Label)      
        
        # add a check box
        self.TheCheckBox_TDVI = QCheckBox('')
        self.TheCheckBox_TDVI.setChecked(False)
        HorizontalLayout1.addWidget(self.TheCheckBox_TDVI) # add the checkbox to the layout 
        
        RVI_Label = QLabel("RVI-->")
        HorizontalLayout1.addWidget(RVI_Label)              
        
        # add a check box
        self.TheCheckBox_RVI = QCheckBox('')
        self.TheCheckBox_RVI.setChecked(False)
        HorizontalLayout1.addWidget(self.TheCheckBox_RVI) # add the checkbox to the layout 
        
        DVI_Label = QLabel("DVI-->")
        HorizontalLayout1.addWidget(DVI_Label)              
        
        # add a check box
        self.TheCheckBox_DVI = QCheckBox('')
        self.TheCheckBox_DVI.setChecked(False)
        HorizontalLayout1.addWidget(self.TheCheckBox_DVI) # add the checkbox to the layout 
        
        ARVI_Label = QLabel("ARVI-->")
        HorizontalLayout1.addWidget(ARVI_Label)              
        
        # add a check box
        self.TheCheckBox_ARVI = QCheckBox('')
        self.TheCheckBox_ARVI.setChecked(False)
        HorizontalLayout1.addWidget(self.TheCheckBox_ARVI) # add the checkbox to the layout 
        
        ###################################################################################################### ROW #2
        
        SAVI_Label = QLabel("SAVI-->")
        HorizontalLayout2.addWidget(SAVI_Label)         
        
        # add a check box
        self.TheCheckBox_SAVI = QCheckBox('')
        self.TheCheckBox_SAVI.setChecked(False)
        HorizontalLayout2.addWidget(self.TheCheckBox_SAVI) # add the checkbox to the layout   
        
        OSAVI_Label = QLabel("OSAVI-->")
        HorizontalLayout2.addWidget(OSAVI_Label)         
        
        # add a check box
        self.TheCheckBox_OSAVI = QCheckBox('')
        self.TheCheckBox_OSAVI.setChecked(False)
        HorizontalLayout2.addWidget(self.TheCheckBox_OSAVI) # add the checkbox to the layout   
        
        GOSAVI_Label = QLabel("GOSAVI-->")
        HorizontalLayout2.addWidget(GOSAVI_Label)         
        
        # add a check box
        self.TheCheckBox_GOSAVI = QCheckBox('')
        self.TheCheckBox_GOSAVI.setChecked(False)
        HorizontalLayout2.addWidget(self.TheCheckBox_GOSAVI) # add the checkbox to the layout 
        
        MSAVI2_Label = QLabel("MSAVI2-->")
        HorizontalLayout2.addWidget(MSAVI2_Label)         
        
        # add a check box
        self.TheCheckBox_MSAVI2 = QCheckBox('')
        self.TheCheckBox_MSAVI2.setChecked(False)
        HorizontalLayout2.addWidget(self.TheCheckBox_MSAVI2) # add the checkbox to the layout
        
        EVI_Label = QLabel("EVI-->")
        HorizontalLayout2.addWidget(EVI_Label)         
        
        self.TheCheckBox_EVI = QCheckBox('')
        self.TheCheckBox_EVI.setChecked(False)
        HorizontalLayout2.addWidget(self.TheCheckBox_EVI) # add the checkbox to the layout         
                
        ######################################################################################################
        
        # setup the size of the window
        self.setGeometry(200, 300,500, 200) # x, y, width, height, in pixels
        self.setWindowTitle("LS8 Vegetation Indice Calculator") # title appears in the top bar of the dialog box
        
        self.show() # display the dialog to the user
        
        ########################################################################################
        # Main code
        
        ########################################################################################
        # Create the app and the dialog box and then run the app until the user closes the dialog
        
TheApplicationObject = QApplication(sys.argv)
MainDialogObject = Main()
sys.exit(TheApplicationObject.exec_())
