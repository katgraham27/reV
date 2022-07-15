# -*- coding: utf-8 -*-
"""
An example single run to get bespoke wind plant layout
"""
import numpy as np
import matplotlib.pyplot as plt
from reV.bespoke.bespoke import BespokeSinglePlant
from reV.bespoke.plotting_functions import plot_poly, plot_turbines,\
    plot_windrose
from reV import TESTDATADIR
from reV.supply_curve.tech_mapping import TechMapping

import json
import os
import shutil
import tempfile

from shapely.geometry import Polygon, MultiPolygon, Point
import shapely.geometry

SAM = os.path.join(TESTDATADIR, 'SAM/i_windpower.json')
EXCL = os.path.join(TESTDATADIR, 'ri_exclusions/ri_exclusions.h5')
RES = os.path.join(TESTDATADIR, 'wtk/ri_100_wtk_{}.h5')
TM_DSET = 'techmap_wtk_ri_100'
AGG_DSET = ('cf_mean', 'cf_profile')

# note that this differs from the
EXCL_DICT = {'ri_srtm_slope': {'inclusion_range': (None, 5),
                               'exclude_nodata': False},
             'ri_padus': {'exclude_values': [1],
                          'exclude_nodata': False},
             'ri_reeds_regions': {'inclusion_range': (None, 400),
                                  'exclude_nodata': False}}

with open(SAM, 'r') as f:
    SAM_SYS_INPUTS = json.load(f)

SAM_SYS_INPUTS['wind_farm_wake_model'] = 2
SAM_SYS_INPUTS['wind_farm_losses_percent'] = 0
del SAM_SYS_INPUTS['wind_resource_filename']
TURB_RATING = np.max(SAM_SYS_INPUTS['wind_turbine_powercurve_powerout'])
SAM_CONFIGS = {'default': SAM_SYS_INPUTS}

def windrose(bsp):
    # plot_windrose(wind directions, wind speeds, wind frequencies,ax=None, colors=None )
    # wind directions = np.arange(bsp._wd_bins[0], bsp._wd_bins[1],bsp._wd_bins[2])
    # wd_dist is the join probability distribution between wind speeds and directions ??? and it adds up to 1. 
    # think kat think !! OHHHH i think i get it !!! 
    # wd_bins:
            # wd_bins : tuple
            # 3-entry tuple with (start, stop, step) for the winddirection
            # binning of the wind joint probability distribution. The stop value
            # is inclusive, so ws_bins=(0, 360, 90) would result in four bins
            # with bin edges (0, 90, 180, 270, 360).
    # so this makes sense ^. and they just copy that over as self._wd_bins 
    # similar with _ws_bins
    # np.arange(start, stop, step). 
    ax = plot_windrose(np.arange(bsp._wd_bins[0], bsp._wd_bins[1],
                       bsp._wd_bins[2]),
                       np.arange(bsp._ws_bins[0], bsp._ws_bins[1],
                       bsp._ws_bins[2]),
                       bsp._wind_dist)
    ax.set_title("wind rose")
    

def dom_wind_dir(bsp):
    """plot windrose
    Parameters
    ----------
    wind_directions : 1D array
        wind direction samples
    wind_speeds : 1D array
        wind speed samples
    wind_frequencies : 2D array
        frequency of wind direction and speed samples
        rows = wind speeds, cols = wind dirs
    """
    wind_directions = np.arange(bsp._wd_bins[0], bsp._wd_bins[1],
                       bsp._wd_bins[2])
    wind_speeds = np.arange(bsp._ws_bins[0], bsp._ws_bins[1],
                       bsp._ws_bins[2])
    wind_frequencies = bsp._wind_dist
    print("wind dir shape: ", np.shape(wind_directions))
    print("wind speeds shape: ", np.shape(wind_speeds))
    print("wind freqs shape: ", np.shape(wind_frequencies))
    print("wind dir : ", wind_directions)
    print("wind speeds : ", wind_speeds)
    print("wind freqs : ", wind_frequencies)
    # check to see if each wind direction has the same sample size -- THEY DONT 
    sample_sizes = np.sum(wind_frequencies, axis = 0)
    print("sample sizes: ", sample_sizes)
    # wind_weights = np.zeros((len(wind_speeds),len(wind_directions)))
    # for i in range(len(wind_speeds)):
    #     for j in range(len(wind_directions)):
    #         wind_weights[i,j] = wind_frequencies[i,j]*wind_speeds[i]
    # wind_weights = np.sum(wind_weights, axis = 0)
    # index = np.argmax(wind_weights)
    # dominant_wind_direction = wind_directions[index]

    dominant_wind_direction = wind_directions[np.argmax(sample_sizes)]
    print("dominant_wind_direction in degrees: ", dominant_wind_direction)
    #mean_wind_direction = np.sum(wind_weights*wind_directions)/np.sum(wind_directions)
    mean_wind_direction = np.sum(wind_directions*sample_sizes)/np.sum(sample_sizes)
    # broadcast_wind_directions = np.ones(4,8)*wind_directions 
    # print("broadcast_wind_directions: ", broadcast_wind_directions)
    print("mean_wind_direction: ", mean_wind_direction)
    # NOTE: wind_directions is in degrees**. 
    # NOTE: should i find the mean dominant wind direction, or just the most dominant? 

    


# NOTE THE RESULTS WERE NOT PASSED IN IF THIS WAS A FUNCITON, JENK CODE
def poly(bsp, title):
    #ax = plot_poly(results["packing_polygons"])
    packing_polygons = bsp.plant_optimizer.packing_polygons
    ax = plot_poly(packing_polygons)
    ax = plot_turbines(bsp.plant_optimizer.x_locations,
                       bsp.plant_optimizer.y_locations,
                       50, ax=ax)
    #ax = plot_poly(MultiPolygon(bsp.plant_optimizer.grid_cells))
    # ax.axis("equal")
    # ax.set_title("packed polys, packed points")
    ax.set_title(title)
    n= len(bsp.plant_optimizer.x_locations)
    print(title,": ", n)
    #lot_poly(MultiPolygon(bsp1.plant_optimizer.grid_cells))



def plotdelay(bsp):
        n= len(bsp.plant_optimizer.x_locations)
        #print("n, ", n)
        plt.figure() 
        plt.axis("square")
        for i in range(n):
            plt.gca().cla()
            plot_poly(bsp.plant_optimizer.grid_history[i], ax = plt.gca())
            plot_turbines(bsp.plant_optimizer.turbx_history[i], bsp.plant_optimizer.turby_history[i], 50, ax=plt.gca())
            plt.pause(0.7)
            #plt.show() 
        #plt.show()

def turb_compare(bsp1,bsp2):
        n= len(bsp1.plant_optimizer.x_locations)
        packing_polygons = bsp1.plant_optimizer.packing_polygons
        print("n, ", n)
        plt.figure() 
        plt.axis("square")
        #ax = plot_poly(results1["packing_polygons"])
        ax = plot_poly(packing_polygons)
        ax = plot_turbines(bsp1.plant_optimizer.x_locations,
                       bsp1.plant_optimizer.y_locations,
                       50, ax=ax, color = "r")
        ax = plot_turbines(bsp2.plant_optimizer.x_locations,
                       bsp2.plant_optimizer.y_locations,
                       50, ax=ax, color = "b")
        ax.axis("equal")
        ax.set_title("packed polys, packed points")            
        plt.show()

def convert_ang(theta): 
    while (theta > 360): 
        theta = theta - 360
    print("theta: ", theta) 
    theta = 360 - theta + 90
    print ("new theta: ", theta)
    # KATNOTE - i do not know a clear mapping between the clockwise and offset wind angles and a regular angle.
    # but i do know that there are only 8 options for wind angles so i could just hard code those mappings ?

    
# def cost_function(x):
#     """dummy cost function"""
#     R = 0.1
#     return 200 * x * np.exp(-x / 1E5 * R + (1 - R))


# def objective_function(aep, cost):
#     """dummy objective function"""
#     return cost / aep


if __name__ == "__main__":

    capital_cost_function = """200 * system_capacity * np.exp(-system_capacity /
        1E5 * 0.1 + (1 - 0.1))"""
    objective_function = "capital_cost / aep"
    variable_operating_cost_function = "0.0"
    fixed_operating_cost_function = "0.0"

    output_request = ('system_capacity', 'cf_mean', 'cf_profile')
    gid = 35
    with tempfile.TemporaryDirectory() as td:
        excl_fp = os.path.join(td, 'ri_exclusions.h5')
        res_fp = os.path.join(td, 'ri_100_wtk_{}.h5')
        shutil.copy(EXCL, excl_fp)
        shutil.copy(RES.format(2012), res_fp.format(2012))
        shutil.copy(RES.format(2013), res_fp.format(2013))
        res_fp = res_fp.format('*')

        TechMapping.run(excl_fp, RES.format(2012), dset=TM_DSET, max_workers=1)
        orig_bsp = BespokeSinglePlant(gid, excl_fp, res_fp, TM_DSET,
                                 SAM_SYS_INPUTS,
                                 objective_function, capital_cost_function,
                                 fixed_operating_cost_function,
                                 variable_operating_cost_function,
                                 ga_kwargs={'max_time': 20},
                                 excl_dict=EXCL_DICT,
                                 output_request=output_request
                                 )
        # pj_grid_bsp = BespokeSinglePlant(gid, excl_fp, res_fp, TM_DSET,
        #                          SAM_SYS_INPUTS,
        #                          objective_function, capital_cost_function,
        #                          fixed_operating_cost_function,
        #                          variable_operating_cost_function,
        #                          ga_kwargs={'max_time': 20},
        #                          excl_dict=EXCL_DICT,
        #                          output_request=output_request
        #                          )
        # recurse_bsp = BespokeSinglePlant(gid, excl_fp, res_fp, TM_DSET,
        #                          SAM_SYS_INPUTS,
        #                          objective_function, capital_cost_function,
        #                          fixed_operating_cost_function,
        #                          variable_operating_cost_function,
        #                          ga_kwargs={'max_time': 20},
        #                          excl_dict=EXCL_DICT,
        #                          output_request=output_request
        #                          )

        results1 = orig_bsp.run_plant_optimization()
        # results2 = pj_grid_bsp.run_plant_optimization()
        # results3 = recurse_bsp.run_plant_optimization()
        print("ORIGINAL:::::")
        orig_bsp.plant_optimizer.define_exclusions()
        #packing_polygons = bsp1.plant_optimizer.packing_polygons
        #full_polyons = bsp1.plant_optimizer.full_polygons
        orig_bsp.plant_optimizer.initialize_packing()
        # print("PJ GRID:::")
        # pj_grid_bsp.plant_optimizer.define_exclusions()
        # pj_grid_bsp.plant_optimizer.initialize_packing()
        # print("RECURSE GRID:::")
        # recurse_bsp.plant_optimizer.define_exclusions()
        # recurse_bsp.plant_optimizer.initialize_packing()
    #poly(bsp1)
    #plt.show()
    

    #rotor_diameter = bsp.sam_sys_inputs["wind_turbine_rotor_diameter"]
                    
    # windrose(bsp)
    #poly(bsp1, "recursive grid method")
    poly(orig_bsp, "Original")
    # poly(pj_grid_bsp, "PJ grid function")
    # poly(recurse_bsp, "Recurse grid function")
    #plotdelay(bsp1)
    #plotdelay(bsp1)
    #plt.plot(bsp1.plant_optimizer.grid_cells)
    #plt.show()
    #plt.figure(2)
    #plt.plot(bsp1.plant_optimizer.grid_cells)
    #plt.show()
    # for i in range(len(bsp1.plant_optimizer.grid_cells)):
    #       plot_poly(bsp1.plant_optimizer.grid_cells[i])
    #     print(bsp1.plant_optimizer.grid_cells[i].exterior.coords)
    #     plt.plot(bsp1.plant_optimizer.grid_cells[i].exterior.coords)
    # turb_compare(bsp1,bsp2)
    #plot_poly(MultiPolygon(bsp1.plant_optimizer.grid_cells))

    # dom_wind_dir(bsp1)
    # windrose(bsp1)
    plt.show()

    # convert_ang(45)
    # convert_ang(90)
    # convert_ang(180)
    # convert_ang(270)
    # convert_ang(315)


