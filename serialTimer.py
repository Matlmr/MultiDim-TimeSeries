import os
import csv
import rasterio
import numpy as np
import geopandas as gpd

from rasterstats import zonal_stats
from shapely.geometry import Point, Polygon, box


class SerialTimer():

    def __init__(self,
                 aoi_file=None,
                 aoi_box=[15.00, 48.00, 15.01, 48.01],
                 raster_dir='/home/volume2/Austria/Timeseries/',
                 polygon_file='/home/volume2/Austria/GroundTruth/invekos_schlaege_polygon.shp',
                 out_file='/home/mlamarre/Documents/FieldTimeSeries',
                 time_stamps=np.arange(1, 61),
                 dimensions=['BS.VH', 'BS.VV'],
                 quick_check=True):

        self.raster_dir = raster_dir
        self.polygon_file = polygon_file
        self.out_file = out_file
        if aoi_file:
            with open(aoi_file, mode='r') as file:
                csv_reader = csv.DictReader(file)
                for row in csv_reader:
                    self.aoi = Polygon(row['WKT'])
                    break
        else:
            self.aoi = box(minx=aoi_box[0], miny=aoi_box[1], maxx=aoi_box[2], maxy=aoi_box[3])
        self.time_stamps = time_stamps
        self.dimensions = dimensions
        self.quick_check = quick_check
        self.scanned = False

    def readPolygons(self):
        if os.path.isfile(self.polygon_file):
            print('Reading shape file...')
            self.polygons = gpd.read_file(self.polygon_file)
            self.poly_crs = self.polygons.crs
        else:
            print('ERROR: wrong shape file')

    def readRaster(self, time_stamp, dimension):
        filename = '.'.join([self.raster_dir + str(time_stamp), dimension, 'tif'])
        if os.path.isfile(filename):
            print('Reading raster file...')
            self.raster = rasterio.open(filename)
            self.rast_crs = self.raster.crs
            # Get the values of raster as numpy array and the transform of the raster
            self.array = self.raster.read(1)
            self.transform = self.raster.transform

        else:
            print('ERROR: wrong raster file')

    def selectPolygons(self):
        # Transform AoI to GeoDataFrame to change CRS and then back to Polygon to use within function
        print('Original coordinates: {}'.format(self.aoi.bounds))
        AOI = gpd.GeoDataFrame({'geometry': self.aoi}, index=[0], crs=self.rast_crs)
        AOI = AOI.to_crs(self.poly_crs)
        AOI = box(minx=AOI.bounds.loc[0]['minx'],
                  miny=AOI.bounds.loc[0]['miny'],
                  maxx=AOI.bounds.loc[0]['maxx'],
                  maxy=AOI.bounds.loc[0]['maxy'])
        print('Transformed coordinates: {}'.format(AOI.bounds))

        # Iterate over all the polygons and check if the 1st point is within AOI
        idxs = []
        for index, row in self.polygons.iterrows():
            if index % 2500 == 0:
                print('Scanned {} out of {} fields. ({}%)'.format(index, self.polygons.shape[0],
                                                                  index / self.polygons.shape[0] * 100), end='\r')
            if self.quick_check:
                if (Point(row.geometry.exterior.coords[0]).within(AOI)):
                    idxs.append(index)
            else:
                all_inside = True
                for point in row.geometry.exterior.coords:
                    if not point.within(AOI):
                        all_inside = False
                if all_inside:
                    idxs.append(index)
        print('{} fields were selected'.format(len(idxs)))

        # Keep the polygons of interest and change their CRS
        self.polys_OI = self.polygons.iloc[idxs]
        self.polys_OI = self.polys_OI.to_crs(self.rast_crs)

        # Drop irrelevant columns and reset indices
        self.polys_OI = self.polys_OI.drop(columns=['FS_KENNUNG', 'SL_FLAECHE'])
        self.polys_OI = self.polys_OI.reset_index(drop=True)

    def getZonalStats(self, time_stamp, dimension):
        # Selected polygons and get the zonal stats
        polys_stats = zonal_stats(self.polys_OI['geometry'], self.array, affine=self.transform, stats=['mean'])
        polys_stats = [v['mean'] for v in polys_stats]
        self.polys_OI[dimension + str(time_stamp)] = polys_stats

    def mergeColumns(self, dimension):
        # Merge columns of same time series
        columns_series = [dimension + str(t) for t in self.time_stamps]
        self.polys_OI[dimension] = self.polys_OI[columns_series].values.tolist()
        # Drop indivuals columns
        self.polys_OI = self.polys_OI.drop(columns=columns_series)

    def saveCSV(self):
        print('Writing DataFrame to CSV...')
        self.polys_OI.to_csv(self.out_file)

    def routine(self):

        # Read all the polygons
        self.readPolygons()

        # Iterate over dimensions and time
        for d in self.dimensions:
            for t in self.time_stamps:

                # Read each dimension- and time-specific raster
                self.readRaster(time_stamp=t, dimension=d)

                # Select the polygons only once
                if not self.scanned:
                    self.selectPolygons()
                    self.scanned = True

                # Get the information for all the polygons
                self.getZonalStats(time_stamp=t, dimension=d)

            # Merge the columns into a time series
            self.mergeColumns(dimension=d)

        # Save the dataframe
        self.saveCSV()