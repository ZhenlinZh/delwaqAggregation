# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 14:55:53 2018

@author: zhenlinz
"""

import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon,LineString, multipolygon
import copy
from shapely.ops import cascaded_union

shpfileIn = "./AggProject/bathymetryPositivePlus1000NH4.shp"
shpfileOut = "./AggProject/bathymetryPositivePlus1000NH4_out.shp"

gdf = gpd.read_file(shpfileIn)

aggGriddata = gdf.dissolve(by='CLUSTER_ID')
aggGriddata.drop( ['SOURCE_ID','NH4'],axis=1,inplace=True)


#%% Now check which clusters are polygons with a hole
gl = aggGriddata.geometry.length
gel = aggGriddata.geometry.exterior.length
ind = np.argwhere(gl!=gel)[:,0]
print("Total number of polygons need to be de-aggregated is: "+str(len(ind)))
clusteri = len(aggGriddata)
gdf_new = copy.deepcopy(gdf)
gdf_new['new_ID'] = 0
#%% find inner ring and outer ring of the polygon

def PolyRings(poly):
    """ This function finds the coordinates of inner and outer rings of the 
    supplied polygon obj
    To see the coordinates of inner or outer rings, do Outerring_coords[i]
    """
    Outerring_coords = poly.exterior.coords
    Innerring_coords = []
    Innerrings = poly.interiors
    for ring in Innerrings:
        Innerring_coords.append(ring.coords)
    return Outerring_coords, Innerring_coords

def PolyHollow(poly):
    """ quick algorithm to check if the polygon is hollow
    """
    if poly.exterior.length==poly.length:
        return False
    else:
        return True

def PolySplit(poly,cells):
    """ split a polygon vertically into 2 using a representative point 
    from the first hole
    """
    Outering_coords, Innerring_coords = PolyRings(poly) 
    # for each inner ring, create a polygon and find a representative point ( a 
    # point that's guaranteed to be within the polygon). We only need the x
    rpoint = []
    for ring in Innerring_coords:
        rpoint.append( Polygon(ring).representative_point() )
    indnewCluster1 = (cells.geometry.centroid.x>=rpoint[0].x) & (cells.within(poly))
    indnewCluster2 = (cells.geometry.centroid.x<rpoint[0].x) & (cells.within(poly))
    cl1 = cascaded_union(cells[indnewCluster1].geometry)
    cl2 = cascaded_union(cells[indnewCluster2].geometry)
    return cl1, cl2

def DissAggPoly(poly):
    """ dis-aggregate multiple polygons into multiple polygons
    """
    polyv = []
    if type(poly)==Polygon: # single polygon
        polyv.append(poly)
    else:
        for pi in poly:
            polyv.append(pi)
    return polyv
        
        
    
#%% Divid the polygon until no holes are in them.
for indi in ind:
    ID = indi+1
    poly =aggGriddata.geometry[ID]
    cells = gdf[gdf.CLUSTER_ID==ID]
    
    
    cl1, cl2 = PolySplit(poly,cells)
    cl1 = DissAggPoly(cl1)
    cl2 = DissAggPoly(cl2)
    polylist = np.concatenate((cl1,cl2))
       
    FurtherDivide = True
    polyv = []
    
    while FurtherDivide==True:
        polyinv = np.array([])
        for i, p in enumerate(polylist):          
            if PolyHollow(p)==True:
                # further division is required.
                cl3,cl4 = PolySplit(p,cells)
                cl3 = DissAggPoly(cl3)
                cl4 = DissAggPoly(cl4)
                polyinv = np.concatenate((polyinv,cl3,cl4))
            else:
                # this is one valid poly
                polyv.append(p)
        if len(polyinv)>0:
            polylist= copy.deepcopy(polyinv)
        else:
            FurtherDivide = False
    
    # create new_ID for each of the divided polygons. 
    for i in np.linspace(1,len(polyv)-1,len(polyv)-1).astype(int):
        gdf_new.loc[gdf_new.geometry.within(polyv[i]),'new_ID'] = i        


#%% output new shape file
aggGriddata_new = gdf_new.dissolve(by=['CLUSTER_ID','new_ID'])
aggGriddata_new.drop( ['SOURCE_ID','NH4'],axis=1,inplace=True)
nrows = aggGriddata_new.count().values[0]
aggGriddata_new['ID'] = np.linspace(1,nrows,nrows).astype(int)
aggGriddata_new.crs = "+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs"
aggGriddata_new.to_file(shpfileOut)

        
#%% create multipolygon object in matplotlib and plot       
def diaggPlot(poly,polyv):         
    from matplotlib.collections import PatchCollection
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    from matplotlib import cm
    pmap = []
    for p in polyv:
         pmap_i = patches.Polygon(np.transpose(p.exterior.xy))
         pmap.append(pmap_i)
    fig,ax = plt.subplots()
    colors = cm.jet( (np.linspace(1,len(polyv),len(polyv))/len(polyv)*255).astype(int))
    p1 = PatchCollection(pmap,alpha=0.8,edgecolor=[0.5]*3,facecolor = colors)
    plt.plot(poly.exterior.xy[0],poly.exterior.xy[1],'k')
    ax.add_collection(p1)
    for p in poly.interiors:
        plt.plot(p.xy[0],p.xy[1],'k')
    plt.axis('equal')



#%% Plot the results
import dwaq.PostProcessing as dpp
import matplotlib.pyplot as plt
import numpy as np
hgridFile = 'sal_temp_waqgeom.nc' 
labels = gdf['CLUSTER_ID'].values
n=labels.max()+1
perm=np.argsort(np.random.random(n))
labels=perm[labels]
fig,ax = plt.subplots()
dpp.CreateMapFromData(hgridFile,labels,fig=fig,cmap='jet',colorbar=None)
dpp.PlotGrid(hgridFile,fig=fig,plotz=False,facecolors='none',edgecolors='w',linewidth=0.05)
# South Bay limit
plt.xlim((5.5e5,5.8e5))
plt.ylim((4.15e6,4.18e6))

# San Pablo limit
plt.xlim((5.45e5,5.75e5))
plt.ylim((4.20e6,4.23e6))

# Suisun limit
plt.xlim((5.7e5,6.2e5))
plt.ylim((4.205e6,4.23e6))
plt.axis('off')
