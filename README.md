# delwaqAggregation
Delwaq model grid aggregation using ArcGIS-Pro

Model grid aggregation using ArcGIS-Pro:
Supply the scripts with a shapfile of full-resolution grid for Delwaq, add attributes, perform grid clustering based on the attributes (aggregation), and then correct for 
hollow polygons. 
The output is a shapefile with the aggregated grid. 

Workflow:
I’ve written 3 scripts to perform aggregated grid generation using GIS pro. 
1)	Aggregator_Preprocessing.py: adding attributes to the full resolution shapefile for attribute-dependent clustering. 
2)	AcrProAggregator.py: perform attribute-dependent clustering and should be run by
"c:\Program Files\ArcGIS\Pro\bin\Python\Scripts\propy.bat" ArcProAggregator.py in cmd
3)	Aggregator_Postprocessing.py: this script is mainly used to correct for the clustered polygons with holes. 
 
The instruction on the arcgis clustering method is here:
http://pro.arcgis.com/en/pro-app/tool-reference/spatial-statistics/spatially-constrained-multivariate-clustering.htm
And the detailed description of the method is here:
http://pro.arcgis.com/en/pro-app/tool-reference/spatial-statistics/how-spatially-constrained-multivariate-clustering-works.htm
The reference is: 
Assuncao, R. M., M. C. Neves, G. Camara, and C. Da Costa Freitas. 2006. "Efficient Regionalisation Techniques for Socio-economic Geographical Units using Minimum Spanning Trees" in International Journal of Geographical Information Science 20 (7): 797–811.
