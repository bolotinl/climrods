import geopandas as gpd
import pandas as pd
import re
import os
from urllib.request import urlopen
import urllib
import glob
import matplotlib.pyplot as plt
import numpy as np

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

        # Convert both shapefiles to the same coordinate system (Albers Equal Area)
        grid = grid.to_crs('+proj=aea +lat_1=29.5 +lat_2=42.5')
        watershed = watershed.to_crs('+proj=aea +lat_1=29.5 +lat_2=42.5')


        # Intersect grid with watershed(s)
        intsct = watershed.overlay(grid, how = 'intersection')

        # Format grid cell IDs as they will need to be formatted for the download URL
        new_id_col = []
        for i in range(len(intsct)):
            new_id = 'X'+str(intsct['NLDAS_X'][i]).zfill(3)+'-Y'+str(intsct['NLDAS_Y'][i]).zfill(3)
            new_id_col.append(new_id)

        intsct['nldas_id'] = new_id_col

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

    def download(self, out_dir):
        missing_cells = []
        for i in self.urls:
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

        overlapping_area = self.intsct.geometry.area #m2
        overlapping_area = overlapping_area/1000000 #km2

        watershed = gpd.read_file(self.shp_path)
        watershed_area = watershed['AREA']/1000000 #km2

        proportions = []
        watersheds = []
        
        for i in range(len(self.intsct)):
            proportion = overlapping_area[i]/watershed_area[self.intsct['save_index'][i]]
            proportions.append(proportion)
            watershed_id = self.intsct['save_index'][i]
            watersheds.append(watershed_id)

        proportions_df = pd.DataFrame({'watershed_id': watersheds, 'nldas_id': self.cells, 'proportion_of_watershed': proportions})

        # Create empty time series of desired length to populate with area weighted precip values:
        start = self.start_date_utc+' '+self.start_hour_utc+':00:00'
        end = self.end_date_utc+' '+self.end_hour_utc+':00:00'
        full_timeseries = pd.date_range(start=start, end=end, freq='H')
        precip_values = [0]*len(full_timeseries)

        files = glob.glob('./output/*.csv')

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


# TEST ------------------------------------------------------------------------------------------
# # Create object
# start_dt = '2001-01-01'
# end_dt = '2001-12-31'
# start_hr = '00'
# end_hr = '23'

# test = NLDAS_Downloader(start_dt, end_dt, start_hr, end_hr)

# os.chdir('/Volumes/GoogleDrive/My Drive/Overland Flow MS/Data/NLDAS')
# out_dir = "./output"
# grid_path = './NLDAS_Grid_Reference/NLDAS_Grid_Reference.shp'
# shp_path = './NLDAS_Grid_Reference/truckee_watersheds_0/truckee_HUC.shp'

# # Test function 1
# test.intersect_watershed(shp_path, grid_path)

# # Test function 2
# test.url_builder()

# # Test function 3
# out_dir = './output'
# test.download(out_dir)

# # Test function 4 
# weight_dir = './output/area_weighted/'
# test.area_weight_P(weight_dir)

