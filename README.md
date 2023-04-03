## NOTE: package is a work in progress

# climrods
 Python module for efficiently downloading hourly NLDAS forcing data and calculating area-weighted precipitation time series
 
 Currently, the tool requires that you can either provide 1) a shapefile for the watershed(s) you are interested in or 2) USGS streamgage iDs for watersheds you are interested in. Future development will provide functionality to delineate watersheds from a point within the tool.
 
 Currently, the tool only downloads and analyzes precipitation data, but other parameters will be included in the future.

The NLDAS grid shapefile is required for using this tool. It is available in the sample_data folder or can be downloaded here: https://ldas.gsfc.nasa.gov/faq/nldas

See main.py for an example importing and using the tool.
