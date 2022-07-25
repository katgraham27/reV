# -*- coding: utf-8 -*-
"""
turbine packing module.
"""
from re import X
import numpy as np
from shapely.geometry import Polygon, MultiPolygon, Point, LineString, MultiLineString, MultiPoint
import shapely.geometry
from reV.bespoke.plotting_functions import get_xy
from reV.utilities.exceptions import WhileLoopPackingError

import matplotlib.pyplot as plt
from reV.bespoke.plotting_functions import plot_poly, plot_turbines,\
    plot_windrose

class PackTurbines():
    """Framework to maximize plant capacity in a provided wind plant area.
    """

    def __init__(self, min_spacing, safe_polygons, wd_bins, ws_bins, wind_dist, weight_x=0.003):
        """
        Parameters
        ----------
        min_spacing : float
            The minimum allowed spacing between wind turbines.
        safe_polygons : Polygon | MultiPolygon
            The "safe" area(s) where turbines can be placed without
            violating boundary, setback, exclusion, or other constraints.
        weight_x : float, optional
        """

        self.min_spacing = min_spacing
        self.safe_polygons = safe_polygons
        self.weight_x = weight_x
        # self.orig = orig
        # self.grid = grid
        # self.ang = np.radians(ang)
        #print("ang: ", self.ang)
        # KATQ ABOUT DEFAULT ARGS?
        self._wd_bins = wd_bins
        self._ws_bins = ws_bins
        self._wind_dist = wind_dist
        self.ang = None

        # turbine locations
        self.turbine_x = np.array([])
        self.turbine_y = np.array([])

        # history of self.leftover areas through iterations
        self.grid_history = []
        # history of turbine placements
        self.turbx_history = []
        self.turby_history = []
        # grid
        self.turbsgrid = 0

        # calculate the angle 
        wind_directions = np.arange(self._wd_bins[0], self._wd_bins[1],
                        self._wd_bins[2])
        wind_speeds = np.arange(self._ws_bins[0], self._ws_bins[1],
                        self._ws_bins[2])
        wind_frequencies = self._wind_dist
        sample_sizes = np.sum(wind_frequencies, axis = 0)
        #print("Wind freqs : ", wind_frequencies)
        #print("sample_sizes: ", sample_sizes)
        dominant_wind_direction = wind_directions[np.argmax(sample_sizes)]
        #mean_wind_direction = np.radians(np.sum(wind_directions*sample_sizes)/np.sum(sample_sizes))
        # Convert angle to unit circle (0 degrees is in the east direction and angles go counterclockwise
        #print("dominant wind direction (north is 0 degrees, clockwise): ", dominant_wind_direction)
        while (dominant_wind_direction > 360): 
            dominant_wind_direction = dominant_wind_direction - 360
            #print("theta: ", theta) 
        dominant_wind_direction = 360 - dominant_wind_direction + 90
        #print ("new theta: ", theta)      
        self.ang = np.radians(dominant_wind_direction)
        #print("dominant wind direction (unit circle): ", np.degrees(self.ang))

    # def dom_wind_dir(self):
    #     """Find the dominant wind direction
    #     Parameters
    #     ----------
    #     wind_directions : 1D array
    #         wind direction samples
    #     wind_speeds : 1D array
    #         wind speed samples
    #     wind_frequencies : 2D array
    #         frequency of wind direction and speed samples
    #         rows = wind speeds, cols = wind dirs
    #     """
    #     wind_directions = np.arange(self._wd_bins[0], self._wd_bins[1],
    #                     self._wd_bins[2])
    #     wind_speeds = np.arange(self._ws_bins[0], self._ws_bins[1],
    #                     self._ws_bins[2])
    #     wind_frequencies = self._wind_dist
    #     sample_sizes = np.sum(wind_frequencies, axis = 0)
    #     dominant_wind_direction = wind_directions[np.argmax(sample_sizes)]
    #     #mean_wind_direction = np.radians(np.sum(wind_directions*sample_sizes)/np.sum(sample_sizes))
    #     # Convert angle to unit circle (0 degrees is in the east direction and angles go counterclockwise
    #     while (dominant_wind_direction > 360): 
    #         dominant_wind_direction = dominant_wind_direction - 360
    #         #print("theta: ", theta) 
    #         dominant_wind_direction = 360 - dominant_wind_direction + 90
    #         #print ("new theta: ", theta)      
    #     return np.radians(dominant_wind_direction)


    # PJ's grid defining function 
    def discrete_grid(self, x_spacing,y_spacing,shear,rotation,center_x,center_y,boundary_setback,boundary_poly):
        """
        returns grid turbine layout. Assumes the turbines fill the entire plant area
        Args:
        x_spacing (Float): grid spacing in the unrotated x direction (m)
        y_spacing (Float): grid spacing in the unrotated y direction (m)
        shear (Float): grid shear (rad)
        rotation (Float): grid rotation (rad)
        center_x (Float): the x coordinate of the grid center (m)
        center_y (Float): the y coordinate of the grid center (m)
        boundary_poly (Polygon): a shapely Polygon of the wind plant boundary
        Returns
        return_x (Array(Float)): turbine x locations
        return_y (Array(Float)): turbine y locations
        """
        #shrunk_poly = boundary_poly.buffer(-boundary_setback)
        shrunk_poly = boundary_poly
        if shrunk_poly.area <= 0:
            return np.array([]), np.array([])
        # create grid
        minx, miny, maxx, maxy = shrunk_poly.bounds
        width = maxx-minx
        height = maxy-miny
        center_point = Point((center_x,center_y))
        poly_to_center = center_point.distance(shrunk_poly.centroid)
        width = np.max([width,poly_to_center])
        height = np.max([height,poly_to_center])
        nrows = int(np.max([width,height])/np.min([x_spacing,y_spacing]))*2 + 1
        ncols = nrows
        # xlocs = np.arange(-700,ncols+700)*x_spacing
        xlocs = np.arange(0,ncols)*x_spacing
        if (shear != 0):
            y_spacing = y_spacing*np.sin(shear)
        ylocs = np.arange(0,nrows)*y_spacing
        # ylocs = np.arange(-700,nrows+700)*y_spacing        
        row_number = np.arange(0,nrows)
        #row_number = np.arange(0,nrows+600)
        d = np.array([i for x in xlocs for i in row_number])
        layout_x = np.array([x for x in xlocs for y in ylocs]) + d*y_spacing*np.tan(shear)
        layout_y = np.array([y for x in xlocs for y in ylocs])
        # rotate
        rotate_x = np.cos(rotation)*layout_x - np.sin(rotation)*layout_y
        rotate_y = np.sin(rotation)*layout_x + np.cos(rotation)*layout_y
        # move center of grid
        rotate_x = (rotate_x - np.mean(rotate_x)) + center_x
        rotate_y = (rotate_y - np.mean(rotate_y)) + center_y
        # get rid of points outside of boundary polygon
        meets_constraints = np.zeros(len(rotate_x),dtype=bool)
        for i in range(len(rotate_x)):
            pt = Point(rotate_x[i],rotate_y[i])
            if shrunk_poly.contains(pt) or shrunk_poly.touches(pt):
                meets_constraints[i] = True
        # arrange final x,y points
        return_x = rotate_x[meets_constraints]
        return_y = rotate_y[meets_constraints]

        # return_x = xlocs
        # return_y = ylocs

        # return_x = rotate_x
        # return_y = rotate_y

        # return_x = layout_x
        # return_y = layout_y        
        return return_x, return_y

    def is_feas_point(self, point, poly):
        # SECOND VERSION BLAHBLAHBLAHBLAHBLAHBLAH
        bool = False
        if (poly == None):
            poly = self.leftover
        if (poly.geom_type == 'MultiPolygon' or poly.geom_type == 'MultiLineString'):
            for k in range (len(poly)):
                        if (poly[k].contains(point)):
                            bool = True
            
        else:
            if (poly.contains(point)):
                bool = True
            #print("leftover type: ", type(self.leftover))
        return bool


        # # FIRST VERSION BLAHBLAHBLAHBLAHBLAHBLAH
        # bool = False
        # if (poly == None):
        #     poly = self.leftover
        # if poly.geom_type == 'MultiPolygon':
        #     for k in range (len(poly)):
        #                 if (poly[k].contains(point)):
        #                     bool = True
            
        # elif poly.geom_type == 'Polygon':
        #     if (poly.contains(point)):
        #         bool = True
        # else:
        #     print("leftover type: ", type(self.leftover))
        #     # raise IOError('Shape is not a polygon.')
        # return bool

    def recurse_grid(self, pointx, pointy):
        theta = self.ang - np.pi/6
        radius = self.min_spacing
        for i in range(0,6):
            xposition = pointx + radius*np.cos(theta+i*np.pi/3)
            yposition = pointy + radius*np.sin(theta+i*np.pi/3)
            new_point = Point(xposition, yposition)
            if (self.is_feas_point(point = new_point, poly = self.leftover)):
                            # PLACE TURBINE. 
                            self.turbine_x = np.append(self.turbine_x,
                                                    xposition)
                            self.turbine_y = np.append(self.turbine_y,
                                                    yposition)
                            new_turbine = new_point.buffer(self.min_spacing, resolution = 100)
                            self.leftover = self.leftover.difference(new_turbine)
                            if isinstance(self.leftover, Polygon):
                                self.leftover = MultiPolygon([self.leftover])
                            #self.grid_history.append(self.leftover)
                            #self.turbx_history.append(self.turbine_x)
                            #self.turby_history.append(self.turbine_y) 
                            self.turbsgrid +=1
                            self.recurse_grid(xposition,yposition)

    def pack_turbines_poly(self):
        """Fast packing algorithm that maximizes plant capacity in a
        provided wind plant area. Sets the the optimal locations to
        self.turbine_x and self.turbine_y
        """

        if self.safe_polygons.area > 0.0:
            can_add_more = True
            self.leftover = MultiPolygon(self.safe_polygons)
            self.grid_history.append(self.leftover)
            iters = 0

            # -----------------------------------------------------------------
            # PJ'S GRID -------------------------------------------------------
            # -----------------------------------------------------------------
            if(False):
                #Get all polygons in the safe_polygon MultiPolygon, and save bounds
                allpolys = list(self.leftover)
                cell_width = self.min_spacing
                cell_height = self.min_spacing
                # for hex grid, shear should  be pi/3. For square grid, shear should be 0
                grid_shear = np.pi/3
                grid_rotation = -np.pi/6 + self.ang
                minx = 1000000
                miny = 1000000
                maxx = -1
                maxy = -1
                for i in range (len(allpolys)): 
                    tempminx, tempminy, tempmaxx, tempmaxy = allpolys[i].bounds
                    if (tempminx<minx):
                        minx = tempminx
                    if (tempminy<miny):
                        miny = tempminy
                    if (tempmaxx>maxx):
                        maxx = tempmaxx
                    if (tempmaxy>maxy):
                        maxy = tempmaxy
                bbox = Polygon([[minx,miny], [minx, maxy], [maxx, maxy], [maxx, miny]])
                midx = (minx + maxx)/2
                midy = (miny + maxy)/2
                xpoints, ypoints = self.discrete_grid(x_spacing = cell_width, 
                                            y_spacing = cell_height,
                                            shear = grid_shear,
                                            rotation = grid_rotation,
                                            center_x = midx, 
                                            center_y = midy, 
                                            boundary_setback = 0,
                                            boundary_poly = bbox)
                # Place turbine if grid point is feasible
                for i in range (len(xpoints)):
                    for k in range (len(allpolys)):
                        if (allpolys[k].contains(Point(xpoints[i], ypoints[i]))):
                            
                            self.turbine_x = np.append(self.turbine_x,
                                                    xpoints[i])
                            self.turbine_y = np.append(self.turbine_y,
                                                    ypoints[i])
                            new_turbine = Point(xpoints[i], ypoints[i]).buffer(self.min_spacing)                                 
                            # Update leftover 
                            self.leftover = self.leftover.difference(new_turbine)
                            if isinstance(self.leftover, Polygon):
                                self.leftover = MultiPolygon([self.leftover])
                            # self.grid_history.append(self.leftover)
                            # self.turbx_history.append(self.turbine_x)
                            # self.turby_history.append(self.turbine_y)
                # KAT KAT KAT KAT KAT KAT KAT. THE FOLLWOING PRINT STATEMENT IS IMPORTANT!!!
                # print("number of turbines grid-packed: ", len(self.turbine_x))
                # Fill in the rest of the space
                #can_add_more = False
                if self.leftover.area > 0:
                    can_add_more = True
                while can_add_more:
                    iters += 1
                    if iters > 10000:
                        msg = ('Too many points placed in packing algorithm')
                        raise WhileLoopPackingError(msg)     
                    if self.leftover.area > 0:
                        nareas = len(self.leftover.geoms)
                        areas = np.zeros(len(self.leftover.geoms))
                        points = 0
                        # for i in range(nareas):
                        #     areas[i] = self.leftover.geoms[i].area
                        #     if (areas[i] != 0):
                        #         points += len(self.leftover.geoms[i].exterior.coords[:])
                        #     else:
                        #         points += len(self.leftover.geoms[i].coords[:])
                        #         print("hi")

                        m = min(i for i in areas if i > 0)
                        ind = np.where(areas == m)[0][0]
                        # smallest_area = self.leftover.geoms[np.argmin(areas)]
                        smallest_area = self.leftover.geoms[ind]
                        exterior_coords = smallest_area.exterior.coords[:]
                        x, y = get_xy(exterior_coords)
                        metric = self.weight_x*x + y
                        index = np.argmax(metric)
                        # Place turbine
                        self.turbine_x = np.append(self.turbine_x,
                                            x[index])
                        self.turbine_y = np.append(self.turbine_y,
                                            y[index])
                        new_turbine = Point(x[index],
                                            y[index]
                                            ).buffer(self.min_spacing, resolution=100)

                    else:
                        break
                    
                    self.leftover = self.leftover.difference(new_turbine)
                    if isinstance(self.leftover, Polygon):
                        self.leftover = MultiPolygon([self.leftover])
                    # self.grid_history.append(self.leftover)
                    # self.turbx_history.append(self.turbine_x)
                    # self.turby_history.append(self.turbine_y)


            # -----------------------------------------------------------------
            # RECURSIVE GRID --------------------------------------------------
            # -----------------------------------------------------------------
            if (True):
                if self.leftover.area > 0:
                    can_add_more = True
                else:
                    can_add_more = False
                
                while can_add_more:
                    iters += 1
                    if iters > 10000:
                            msg = ('Too many points placed in packing algorithm')
                            raise WhileLoopPackingError(msg)
                    
                    if self.leftover.area > 0:
                        nareas = len(self.leftover.geoms)
                        areas = np.zeros(len(self.leftover.geoms))
                        points = 0
                        # Code to understand how much the first recursion places:
                        # if (iters == 2): 
                        #     print("self.turbsgrid after 1 iter: ", self.turbsgrid)
                        #     ax = plot_poly(self.leftover)
                        #     ax.set_title("asdf")
                        #     plt.show()
                        for i in range(nareas):
                            areas[i] = self.leftover.geoms[i].area
                        #     if self.leftover.geoms[i].geom_type == Polygon: 
                        #         points += len(self.leftover.geoms[i].exterior.coords[:])
                        #     else:
                        #         points += len(self.leftover.geoms[i].coords[:])
                        # print the total number of exterior points in the polygon
                        #print(points)
                        m = min(i for i in areas if i > 0)
                        ind = np.where(areas == m)[0][0]
                        # smallest_area = self.leftover.geoms[np.argmin(areas)]
                        smallest_area = self.leftover.geoms[ind]
                        exterior_coords = smallest_area.exterior.coords[:]
                        x, y = get_xy(exterior_coords)
                        metric =self.weight_x*x + y
                        index = np.argmax(metric)
                        self.turbine_x = np.append(self.turbine_x,
                                            x[index])
                        self.turbine_y = np.append(self.turbine_y,
                                            y[index])
                        new_turbine = Point(x[index],
                                            y[index]
                                            ).buffer(self.min_spacing, resolution=100)                      
                        
                        self.leftover = self.leftover.difference(new_turbine)
                        if isinstance(self.leftover, Polygon):
                            #print("type : ", type(self.leftover))
                            self.leftover = MultiPolygon([self.leftover])
                        # KAT KAT KAT KAT KAT KAT KAT. THE FOLLWOING PRINT STATEMENT WAS LEFT IN!!!
                        # if (not type(self.leftover) == MultiPolygon and not type(self.leftover) == Polygon):
                        #     print("leftova type: ", type(self.leftover))
                        # self.grid_history.append(self.leftover)
                        # self.turbx_history.append(self.turbine_x)
                        # self.turby_history.append(self.turbine_y) 
                        
                        self.recurse_grid(pointx = x[index], pointy = y[index])
                    # else: 
                    #     if isinstance(self.leftover, MultiLineString) or isinstance(self.leftover, MultiPoint):
                            
                    
                    # elif len(self.leftover.geoms)>0:
                    #     ngeoms = len(self.leftover.geoms)
                    #     for i in range(ngeoms):

                    else: 
                        # KAT KAT KAT KAT KAT KAT KAT. THE FOLLWOING PRINT STATEMENTS WERE LEFT IN!!!
                        # print("self.turbsgrid: ", self.turbsgrid)
                        # print("iters: ", iters)
                        # print("type of the final leftover area: ", type(self.leftover))
                        # print(self.leftover)
                        break


            # -----------------------------------------------------------------
            # ORIGINAL --------------------------------------------------------
            # -----------------------------------------------------------------
            if (False):
                while can_add_more:                                    
                    iters += 1 
                    if iters > 10000:
                        msg = ('Too many points placed in packing algorithm')
                        raise WhileLoopPackingError(msg)
                    
                    if self.leftover.area > 0:
                        nareas = len(self.leftover.geoms)
                        areas = np.zeros(len(self.leftover.geoms))
                        points = 0
                        for i in range(nareas):
                            areas[i] = self.leftover.geoms[i].area
                            points += len(self.leftover.geoms[i].exterior.coords[:])
                        # print the total number of exterior points in the polygon
                        #print(points)
                        m = min(i for i in areas if i > 0)
                        ind = np.where(areas == m)[0][0]
                        # smallest_area = self.leftover.geoms[np.argmin(areas)]
                        smallest_area = self.leftover.geoms[ind]
                        exterior_coords = smallest_area.exterior.coords[:]
                        x, y = get_xy(exterior_coords)
                        metric = self.weight_x * x + y
                        index = np.argmax(metric)
                        self.turbine_x = np.append(self.turbine_x,
                                            x[index])
                        self.turbine_y = np.append(self.turbine_y,
                                            y[index])
                        new_turbine = Point(x[index],
                                            y[index]
                                            ).buffer(self.min_spacing, resolution=200)                      

                    else:
                        break
                    self.leftover = self.leftover.difference(new_turbine)
                    if isinstance(self.leftover, Polygon):
                        self.leftover = MultiPolygon([self.leftover])
                    # self.grid_history.append(self.leftover)
                    # self.turbx_history.append(self.turbine_x)
                    # self.turby_history.append(self.turbine_y)




    def clear(self):
        """Reset the packing algorithm by clearing the x and y turbine arrays
        """
        self.turbine_x = np.array([])
        self.turbine_y = np.array([])


def smallest_area_with_tiebreakers(g):
    """_summary_

    This function helps break ties in the area of two different
    geometries using their exterior coordinate values.

    Parameters
    ----------
    g : _type_
        A geometry object with an `area` and an
        `exterior.coords` coords attribute.

    Returns
    -------
    tuple
        Tuple with the following elements:
            - area of the geometry
            - minimum exterior coordinate (southwest)
            - maximum exterior coordinate (northeast)
    """
    return g.area, min(g.exterior.coords), max(g.exterior.coords)
