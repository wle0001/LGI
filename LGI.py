#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 16:17:01 2024

@author: ellenbw
"""

import xarray as xr
import numpy as np
import rioxarray as rxr
import rasterio as rio
import geopandas as gpd
from pathlib import Path
from datetime import datetime, timedelta
#from tqdm.autonotebook import tqdm
#from tqdm import tqdm
from typing import Tuple
import glob

def get_shapefile()->Tuple[gpd.GeoDataFrame,gpd.GeoDataFrame]:
    """ Get shapefiles for state/county boundaries """
        
    # Assumes shapefiles in same directory as source code
    us_st = gpd.read_file(Path('s_08mr23/'))
    states = us_st[
        (us_st.STATE == 'AL')|\
        (us_st.STATE == 'FL')|\
        (us_st.STATE == 'GA')|\
        (us_st.STATE == 'NC')|\
        (us_st.STATE == 'TN')|\
        (us_st.STATE == 'MS')|\
        #(us_st.STATE == 'VA')|\
        (us_st.STATE == 'SC')]

    us_ct = gpd.read_file(Path('c_19se23/'))
    counties = us_ct[
        (us_ct.STATE == 'AL')|\
        (us_ct.STATE == 'FL')|\
        (us_ct.STATE == 'GA')|\
        (us_ct.STATE == 'NC')|\
        (us_ct.STATE == 'TN')|\
        (us_ct.STATE == 'MS')|\
        #(us_st.STATE == 'VA')|\
        (us_ct.STATE == 'SC')]
        
    return (states,counties)


def get_precip_al_21day(dates_in:list, 
                        st_shp:gpd.GeoDataFrame)-> xr.DataArray:
 #  data_dir = 'https://water.weather.gov/precip/downloads/'
    data_dir = '../precip/'
    #data_dir = 'nexrad/'
    
    ds_arr = []
    for date in dates_in:
        tif = data_dir+date.strftime('%Y%m%d')+'.tif'
        #tif = data_dir+'nws_precip_1day_{}_conus.tif'.format(date.strftime('%Y%m%d'))
    #for tif in glob.glob(data_dir+'*')[:21]:
        #print(f'{data_dir}{date:%Y}/{date:%m}/nws_precip_1day_{date:%Y%m%d}_conus.tif')
        #_ds = rxr.open_rasterio(f'{data_dir}{date:%Y}/{date:%m}/{date:%d}/'
                                #f'nws_precip_1day_{date:%Y%m%d}_conus.tif',
                                #masked=True).sel(band=1)
        _ds = rxr.open_rasterio(tif, masked = True).sel(band = 1)
                
        _ds = _ds.rio.clip(st_shp.geometry.values, st_shp.crs, drop=True).copy()
        #print(rio.open(tifList[0]).meta['nodata'])
        _ds = _ds.where(_ds != -9999.0)/25.4
        #_ds = _ds.where(_ds != -9999.0)
        ds_arr.append(_ds)
        #sys.exit()
    return xr.concat(ds_arr, dim='day')

def calc_weighted_precip(precip:xr.DataArray)->xr.DataArray:
    """ Calculate effective precip based on 21 day total """
    
    # weights for the past 21 days of precipitation
    wts = np.concatenate((np.ones(7),sorted((np.arange(14)+0.5)/14.0,reverse=True)))
    weights = xr.DataArray(wts, dims='day',coords={'day':precip.day})
    pcpn_weighted = precip.weighted(weights)
    
    # Calculate effective precipitation
    return pcpn_weighted.sum(("day"))


def calc_standard_pcpn(date:datetime)->float:
    """ Calculate standard amount of precipitation based on DOY """
    
    # Determine normal amount of effective precipitation based on day of year
    doy = int(f'{date:%j}')

    # Dictionary for adjustments for "standard" precip
    # - Based on code from Kevin Doty, Phillip Jones, Don Moss, & B Norris
    adj = {'a':(1.0,-0.75,240),'b':(250.0,0.75,235)}

    if 121<doy<250:
        return 2.0
    
    key = 'a' if doy <= 121 else 'b'

    std_val = 1.25 + adj[key][1]*np.cos((2*np.pi/adj[key][2])*(doy-adj[key][0]))
    return std_val



#################################################################################
#################################################################################
STATES, COUNTIES = get_shapefile()

DATES = gpd.pd.date_range(datetime(2000,2,1),datetime(2022,12,31))
for date in DATES:
    print(date)
    #sys.exit()
    
    dates = [date - timedelta(days=x) for x in range(21)]
        # DATES = [datetime.today() - timedelta(days=x) for x in range(21)]
    
        # Get shapefiles
    #STATES, COUNTIES = get_shapefile()
    
    # Get Precipitation Data
    TOTAL_PCPN_SE = get_precip_al_21day(dates, STATES)
    # TOTAL_PCPN_AL = TOTAL_PCPN_SE.rio.clip(STATES.iloc[0].geometry.values,
    #                                        STATES.iloc[0].crs, drop=True).copy()
    # Get mask for valid data within clip
    MASK = TOTAL_PCPN_SE.rio.reproject('EPSG:4326').sel(day=1).notnull()
    
    # Calculate effective precipitation
    PCPN_EFF_SE = calc_weighted_precip(TOTAL_PCPN_SE)
    # PCPN_EFF_AL = calc_weighted_precip(TOTAL_PCPN_AL)
    
    # Get "standard" precipitation amount based on day of year
    STD_VAL = calc_standard_pcpn(dates[0])
    
    #LGI calculation
    LGI_SE = PCPN_EFF_SE - STD_VAL
    # LGI_AL = PCPN_EFF_AL - STD_VAL
    
    # Write LGI to a CRS ad convert to EPSG:4326
    LGI_SE.rio.write_crs(TOTAL_PCPN_SE.rio.crs, inplace=True)
    LGI_WGS_SE = LGI_SE.rio.reproject('EPSG:4326').copy()
    LGI_WGS_SE = LGI_WGS_SE.where(MASK, np.nan).copy()
    
    LGI_WGS_SE.rio.write_nodata(LGI_WGS_SE.rio.nodata, encoded=True, inplace=True)
    LGI_WGS_SE.rio.to_raster(f'../LGI_data/LGI_{date:%Y%m%d}_AL.tif',mask=True)
    
    print(f'{date} : DONE!')

'''
sys.exit()



tifList = glob.glob('LGI_p/LGI_202206*.tif')
#sort so the dates are in order
tifList.sort()

#read the first file inot an array

tif = rio.open(tifList[0])
meta = tif.meta
tif1 = tif.read()


#split the name of the file to get the date
#i.e. here the file name is evap_20200901_01_00.tif
#i'm splitting bu "_" so I get a list that looks like this:
#    ['vic', 'al3/evap', '20080101', '01', '00.tif']
# then I take the 3 (index == 2) value...i.e. dts = 20080101
#dts = [tifList[0].split('_')[3]]
dts = [tifList[0].split('_')[1]]

# loop through the the list of tifs, turn each one into an array, then 
#add it to a stack

#Tif 1 is my starting array and has a shape of (1, 96, 71)
# the dimentions are (day, x, y)
# So I'm going to stack all the files on top of each other
# and end up with an array thay is (NumberOfDays, 96, 71)
tif_arr = []
for tif in tifList:
    
    #while I'm here I'm going to grab each date as well
    dts.append(tif.split('_')[1])
    print(tif.split('_')[2])
    
    #reading each file in the tifList[1:] (remember we already have the first one) 
    tif = rio.open(tif).read()
    
    #keep stacking them in Tif1
    tif1 = np.append(tif1,tif, axis = 0)
    #tif_arr.append(tif)

tif1 = np.where(tif1==meta['nodata'], np.nan, tif1)
tif1 = np.where(tif1<-10, np.nan, tif1)

#stack = xr.concat(tif_arr, dim='day')   

####
#reproject wierd prism files 
#need to put this in the script above to check first. 
ds1 = rio.open('prism/20000101.tif')

DATES = gpd.pd.date_range(datetime(2000,1,1),datetime(2022,12,31))
for d in DATES:
        dt = '{}'.format(d.strftime('%Y%m%d'))
        ds2 = rio.open('prism/{}.tif'.format(dt))
        print(dt)
        with rio.open('prism2/'+dt+'.tif', 'w',**ds1.meta) as dst:
                dst.write(ds2.read())
                
                
ds2 = rio.open('prism2/20211231.tif')        
with rio.open('prism/20211231.tif', 'w',**ds1.meta) as dst:
        dst.write(ds2.read())
        
        
        
'''      
    
