import geopandas as gpd
import pandas as pd
import re
import os
from urllib.request import urlopen
import urllib
import glob
import matplotlib.pyplot as plt
import numpy as np
import json
from shapely.geometry import Polygon

class NLDAS_Downloader():
    def __init__(self, start_date_utc, end_date_utc,start_hour_utc=0, end_hour_utc=23, parameter = 'APCPsfc'):
        '''
        NLDAS_Downloaßder class constructor
        :param start_date_utc: desired start date of time series
        :param start_hour_utc: desired start hour of time series, in UTC 
        :param end_date_utc: desired end date of time series
        :param end_hour_utc: desired end hour of time series, in UTC
        :param parameter: code for the NLDAS parameter you are interested in downloading
        '''

        self.start_date_utc = start_date_utc
        self.start_hour_utc = start_hour_utc
        self.end_date_utc = end_date_utc
        self.end_hour_utc = end_hour_utc
        self.parameter = parameter
        
    def watershed_from_gauge(self, usgs_path, gage_id_col, shp_out_path):
        '''
        Function to create a shapefile of all upstream watershed boundaries using USGS gage IDs
        :param usgs_path: path to a CSV file of numeric site numbers for USGS streamgages
        :param gage_id_col: a string specifying the name of which column in the file usgs_path contains the USGS gage IDs
        :param shp_out_path: path where the final shapefile of all watershed boundaries should be saved, ending with .shp
        '''
        usgs_info = pd.read_csv(usgs_path)
        usgs_info[gage_id_col] = usgs_info[gage_id_col].astype(str)
        usgs_info[gage_id_col] = usgs_info[gage_id_col].str.zfill(8)
        gages = usgs_info[gage_id_col]

        # Create empty geodataframe which will be populated in the for loop below
        polygons = gpd.GeoDataFrame(columns=['geometry', 'GAGE_ID'], geometry='geometry')

        for g in gages:
            # Pull boundary data from USGS API
            usgs_api = 'https://labs.waterdata.usgs.gov/api/nldi/linked-data/nwissite/USGS-'+g+'/basin?f=json'
            response = urlopen(usgs_api)
            boundary_data = json.load(response)

            # Pull geometry from the file and convert into a polygon
            coords = boundary_data['features'][0]['geometry']['coordinates'][0]

            longitude = []
            latitude = []

            for i in coords:
                longitude.append(i[0])
                latitude.append(i[1])


            polygon_geom = Polygon(zip(longitude, latitude))
            polygon = gpd.GeoDataFrame(index=[0], crs='EPSG:4326', geometry=[polygon_geom])
            polygon = polygon.assign(GAGE_ID=g)

            # Add all watersheds to the geodataframe we created before
            polygons = pd.concat([polygons, polygon])
        polygons.crs = 'epsg:4326'
        polygons.to_file(filename= shp_out_path, driver="ESRI Shapefile")
        
    def intersect_watershed(self, shp_path, grid_path):
        '''
        :param shp_path: filepath to a shapefile of one or more watershed boundaries
        :param grid_path: filepath to the shapefile of the NLDAS grid, which 1) is located in the sample data folder of the module and 
        2) can be downloaded here: https://ldas.gsfc.nasa.gov/faq/nldas
        '''
        self.shp_path = shp_path

        grid = gpd.read_file(grid_path)
        watershed = gpd.read_file(shp_path)

        # If there are multiple watersheds in the shapefile, keep track of which is which:
        watershed['save_index'] = watershed.index

        # Convert both shapefiles to the same PROJECTED coordinate system (Albers Equal Area)
        grid = grid.to_crs('EPSG:5070')
        watershed = watershed.to_crs('EPSG:5070')

        # Intersect grid with watershed(s)
        intsct = watershed.overlay(grid, how = 'intersection')

        # Format grid cell IDs as they will need to be formatted for the download URL
        new_id_col = []
        for i in range(len(intsct)):
            new_id = 'X'+str(intsct['NLDAS_X'][i]).zfill(3)+'-Y'+str(intsct['NLDAS_Y'][i]).zfill(3)
            new_id_col.append(new_id)

        intsct['nldas_id'] = new_id_col
        intsct['AREA']= intsct['geometry'].area

        self.intsct = intsct

    def url_builder(self):
        self.cells = list(self.intsct['nldas_id'])
        base_url = 'https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi?variable=NLDAS:NLDAS_FORA0125_H.002:'
        url_loc = '&location=NLDAS:'
        url_startDT = '&startDate='
        url_endDT = '&endDate='
        url_end = '&type=asc2'

        urls = []
        for i in self.cells:
            full_url = base_url+self.parameter+url_loc+i+url_startDT+self.start_date_utc+'T'+self.start_hour_utc+url_endDT+self.end_date_utc+'T'+self.end_hour_utc+url_end
            urls.append(full_url)

        self.urls = urls
        self.urls = list(set(self.urls))

    def download(self, out_dir):
        missing_cells = []
        self.out_dir = out_dir
        for i in list(set(self.urls)):
            try:
                with urlopen(i) as webpage:
                    content = webpage.read().decode()
                    content = content.split('\n', 39)[39]

                parameter = re.search('NLDAS_FORA0125_H.002:(.*)&location', i).group(1)
                cell = re.search('&location=NLDAS:(.*)&startDate', i).group(1)

                with open( '{}/{}_{}.csv'.format(out_dir, parameter, cell), 'w' ) as output:
                    output.write( content )
            except urllib.error.HTTPError as e:
                missing_cells.append(cell)           
                print('Failed to extract data at the following URL: {}'.format(i))
                continue

        self.missing_cells = missing_cells

    def area_weight_P(self, weight_dir):
        '''
        Function to calculate and output area weighted NLDAS precipitation timeseries
        :param weight_dir: path to directory where files with area weighted precipitation should be saved
        '''
        if not os.path.exists(weight_dir):
            os.mkdir(weight_dir)

        overlapping_area = self.intsct['AREA'] #m2
        overlapping_area = overlapping_area/1000000 #km2

        watershed = gpd.read_file(self.shp_path)
        watershed = watershed.to_crs('EPSG:5070')

        watershed['AREA'] = watershed['geometry'].area
        watershed_area = watershed['AREA']/1000000 #km2

        proportions = []
        watersheds = []
        
        for i in range(len(self.intsct)):
            proportion = overlapping_area[i]/watershed_area[self.intsct['save_index'][i]]
            proportions.append(proportion)

            if 'GAGE_ID' in self.intsct:
                watershed_id = self.intsct['GAGE_ID'][i]

            else:
                watershed_id = self.intsct['save_index'][i]
            watersheds.append(watershed_id)

        proportions_df = pd.DataFrame({'watershed_id': watersheds, 'nldas_id': self.cells, 'proportion_of_watershed': proportions})

        # Create empty time series of desired length to populate with area weighted precip values:
        start = self.start_date_utc+' '+self.start_hour_utc+':00:00'
        end = self.end_date_utc+' '+self.end_hour_utc+':00:00'
        full_timeseries = pd.date_range(start=start, end=end, freq='H')
        precip_values = [0]*len(full_timeseries)

        files = glob.glob(self.out_dir+'APCPsfc*.csv')

        # Split by watershed (if multiple)
        df_list = [d for _, d in proportions_df.groupby(['watershed_id'])]

        # For each watershed...
        for f in df_list:
            cells = list(f['nldas_id'])
            ws = str(max(f['watershed_id']))

            # Iterate through all data files...
            for i in files:
                cell = i[-13:-4]            
                # If the cell whose data is contained in file i intersects watershed f, read the data...
                if cell in cells:
                    prop = float(f['proportion_of_watershed'][f['nldas_id'] == cell].values[0])
                    d = open(i, mode='r')
                    data = d.readlines()[1:-1]

                    if len(precip_values) == len(data):

                        # Multiply the data by the proportional watershed area contained in the cell, and update precip timeseries
                        for u in range(len(precip_values)):
                            precip_values[u] += (float(data[u].split()[2]))*prop

                    else:
                        print('Error: length of time series does not match length of data; therefore data for NLDAS {} not included'.format(cell))
                        print('P: '+len(precip_values))
                        print('Data:'+len(data))
                        continue

            # Create/output timeseries dataframe:      
            merged_df = pd.DataFrame({'DateTime':full_timeseries, 'Precipitation (mm/hr)': precip_values})
            merged_df.to_csv('{}/area_weighted_watershed{}.csv'.format(weight_dir, ws))

            plt.plot(merged_df['DateTime'], merged_df['Precipitation (mm/hr)'])
            plt.xlabel('DateTime')
            plt.ylabel('Area Weighted Precipitation (mm/hr)')
            plt.savefig('{}/P_timeseries_ws{}.png'.format(weight_dir, ws))
            plt.clf()