# -*- coding: utf-8 -*-
"""
PyTest file for Wind generation in Rhode Island.

Created on Thu Nov 29 09:54:51 2018

@author: gbuster
"""

import os
import pytest
import pandas as pd
import numpy as np
import shutil
import tempfile

from reV.handlers.outputs import Outputs
from reV.generation.generation import Gen
from reV.config.project_points import ProjectPoints
from reV import TESTDATADIR
from reV.utilities.exceptions import SAMExecutionError

from rex import Resource
from rex.utilities.utilities import mean_irrad


def test_forecast():
    """Test several forecast features implemented in reV gen including
    site_data for timezone and elevation input and gid_map for forecast meta
    mapping"""

    res_files_source = TESTDATADIR + '/nsrdb/ri_100_nsrdb_2012.h5'
    sam_files = TESTDATADIR + '/SAM/i_pvwattsv7.json'
    with tempfile.TemporaryDirectory() as td:
        res_file = os.path.join(td, 'temp_2012.h5')
        shutil.copy(res_files_source, res_file)

        with Outputs(res_file, mode='a') as f:
            meta = f.meta
            meta = meta.drop(['timezone', 'elevation'], axis=1)
            del f._h5['meta']
            f._meta = None
            f.meta = meta

        with Outputs(res_file, mode='r') as f:
            assert 'timezone' not in f.meta
            assert 'elevation' not in f.meta

        with Resource(res_file) as res:
            ghi = res['ghi']

        points = ProjectPoints(slice(0, 5), sam_files, 'pvwattsv7',
                               res_file=res_file)
        output_request = ('cf_mean', 'ghi_mean')
        site_data = pd.DataFrame({'gid': np.arange(5), 'timezone': -5,
                                  'elevation': 0})
        gid_map = {0: 20, 1: 20, 2: 50, 3: 51, 4: 51}

        # test that this raises an error with missing timezone
        with pytest.raises(SAMExecutionError):
            _ = Gen.reV_run('pvwattsv7', points, sam_files, res_file,
                            max_workers=1, sites_per_worker=3,
                            out_fpath=None, output_request=output_request)

        gen1 = Gen.reV_run('pvwattsv7', points, sam_files, res_file,
                           max_workers=1, sites_per_worker=3,
                           site_data=site_data,
                           out_fpath=None, output_request=output_request)

        for i in range(5):
            assert np.allclose(gen1.out['ghi_mean'][i], mean_irrad(ghi[:, i]),
                               atol=0.0, rtol=0.001)

        gen2 = Gen.reV_run('pvwattsv7', points, sam_files, res_file,
                           max_workers=1, sites_per_worker=3,
                           site_data=site_data, gid_map=gid_map,
                           out_fpath=None, output_request=output_request)

        for i in range(5):
            j = gid_map[i]
            assert np.allclose(gen2.out['ghi_mean'][i], mean_irrad(ghi[:, j]),
                               atol=0.0, rtol=0.001)
