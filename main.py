from climrods import *

def main():
    start_dt = '2001-01-01'
    end_dt = '2001-12-31'
    start_hr = '00'
    end_hr = '23'

    test = NLDAS_Downloader(start_dt, end_dt, start_hr, end_hr)

    out_dir = "./output"
    grid_path = './sample_data/NLDAS_Grid_Reference.shp'
    shp_path = './sample_data/truckee_HUC.shp'

    # Test function 1
    test.intersect_watershed(shp_path, grid_path)

    # Test function 2
    test.url_builder()

    # Test function 3
    out_dir = './output'
    test.download(out_dir)

    # Test function 4 
    weight_dir = './output/area_weighted/'
    test.area_weight_P(weight_dir)

if __name__ == '__main__': 
    main()
