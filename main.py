from climrods import *

## Example for if you ALREADY have a shapefile of your watersheds of interest
def main():
    start_dt = '2001-01-01'
    end_dt = '2001-12-31'
    start_hr = '00'
    end_hr = '23'

    test = NLDAS_Downloader(start_dt, end_dt, start_hr, end_hr)

    out_dir = "./output/existing_shapefile_example"
    grid_path = './sample_data/NLDAS_Grid_Reference.shp'
    shp_path = './sample_data/truckee_HUC.shp'

    # Test function 1
    test.intersect_watershed(shp_path, grid_path)

    # Test function 3
    test.url_builder()

    # Test function 4
    out_dir = './output'
    test.download(out_dir)

    # Test function 5 
    weight_dir = './output/existing_shapefile_example/area_weighted/'
    test.area_weight_P(weight_dir)

## Example for if you ONLY have the USGS Gage IDs of your watersheds of interest
def main():
    start_dt = '2001-01-01'
    end_dt = '2001-12-31'
    start_hr = '00'
    end_hr = '23'

    test = NLDAS_Downloader(start_dt, end_dt, start_hr, end_hr)

    grid_path = './sample_data/NLDAS_Grid_Reference.shp'
    usgs_path = './sample_data/Q_site_info.csv'
    shp_out_path = './output/usgs_gage_example/san_pedro_watersheds.shp'

    # Test function 1
    test.watershed_from_gauge(usgs_path, 'site_no', shp_out_path)

    # Test function 2
    test.intersect_watershed(shp_out_path, grid_path)

    # Test function 3
    test.url_builder()

    # Test function 4
    out_dir = './output/usgs_gage_example/'
    test.download(out_dir)

    # Test function 5 
    weight_dir = './output/usgs_gage_example/area_weighted/'
    test.area_weight_P(weight_dir)

if __name__ == '__main__': 
    main()
