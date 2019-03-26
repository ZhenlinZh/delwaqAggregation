# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 14:30:40 2018

@author: zhenlinz
"""

import geopandas as gpd
import xarray as xr 

gridshp = "./AggProject/bathymetryPositive.shp"
outgridshp = "./AggProject/bathymetryPositivePlus.shp"
datafile = r'Y:\zhenlin\dwaq\cascade\suisun_nutrient_cycling\ModelResults\dwaq_map6.nc'
varname = ['NO3','NH4']

def AddShpfileProperties(gridshp,outgridshp, datafile,varname,ts=1,pos='top'):
    """ 
    Not checked... 
    Add one or a list of variabiles from datafile to the 
    shapefile    
    gridshp: normally means the dwaq grid shape file
    outgridshp: output grid shape file
    datafile: normally means model results file
    varname: one or list of variables in the datafile;
    ts: time step
    pos: either 'top' or 'avg' over the depth. 
    """
        
    grid = gpd.read_file(gridshp)
    dset = xr.open_dataset(datafile)
    
    if pos == 'top':
        var =  dset[varname][ts,0,:].values
    elif pos == 'avg':
        var = dset[varname][ts,:,:].values.mean(axis='depth')
    else:
        raise ValueError(" 'pos' can either be 'top' or 'avg' ")
           
    grid[varname] = var    
    grid.crs = "+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs"
    grid.to_file(outgridshp)


AddShpfileProperties(gridshp,outgridshp, datafile,'NO3',ts=138,pos='top')
AddShpfileProperties(outgridshp,outgridshp, datafile,'NH4',ts=138,pos='top')