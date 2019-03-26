# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 13:31:50 2018
This script should be performed in ArcGIS Pro python console. 
To test algorithms, a python window can be opened from 
C:\Program Files\ArcGIS\Pro\bin\Python\envs\arcgispro-py3
To run a python script: 
"c:\Program Files\ArcGIS\Pro\bin\Python\Scripts\propy.bat" script.py

@author: zhenlinz
"""

import arcpy
arcpy.env.overwriteOutput = True
workspace = r'C:\Users\zhenlinz\Google Drive\1_Nutrient_Share\1_Projects_NUTRIENTS\Modeling\NOTES_ZZ\Projects\ModelAggregation\SuisunBay\AggProject'
arcpy.env.workspace = workspace

shapefilein = "AggProject/bathymetryPositivePlus"
Ncluster = 1000
shapefileout = shapefilein + str(Ncluster) + 'NH4'
arcpy.stats.SpatiallyConstrainedMultivariateClustering(shapefilein+'.shp', 
                                                       shapefileout, ["NH4"],  #["Depth","NO3","NH4"]
                                                       number_of_clusters=Ncluster,
                                                       spatial_constraints='CONTIGUITY_EDGES_ONLY')


#arcpy.stats.GenerateSpatialWeightsMatrix(shapefilein, "Index",
#                                         "example.swm", "CONTIGUITY_EDGES_ONLY ", "EUCLIDEAN")
#
#arcpy.ConvertSpatialWeightsMatrixtoTable_stats("example.swm","example.dbf")