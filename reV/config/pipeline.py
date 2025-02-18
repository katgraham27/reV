# -*- coding: utf-8 -*-
"""
reV pipeline config

Created on May 28 2019

@author: gbuster
"""
import os

from reV.config.base_analysis_config import AnalysisConfig
from reV.utilities.exceptions import ConfigError, PipelineError


class PipelineConfig(AnalysisConfig):
    """SAM-based analysis config (generation, lcoe, etc...)."""

    def __init__(self, config):
        """
        Parameters
        ----------
        config : str | dict
            File path to config json (str), serialized json object (str),
            or dictionary with pre-extracted config.
        """

        super().__init__(config, run_preflight=False)
        self._check_pipeline()
        self._parse_dirout()
        self._check_dirout_status()

    def _check_pipeline(self):
        """Check pipeline steps input. ConfigError if bad input."""

        if 'pipeline' not in self:
            raise ConfigError('Could not find required key "pipeline" in the '
                              'pipeline config.')

        if not isinstance(self.pipeline, list):
            raise ConfigError('Config arg "pipeline" must be a list of '
                              '(command, f_config) pairs, but received "{}".'
                              .format(type(self.pipeline)))

        for di in self.pipeline:
            for f_config in di.values():
                if not os.path.exists(f_config):
                    raise ConfigError('Pipeline step depends on non-existent '
                                      'file: {}'.format(f_config))

    def _parse_dirout(self):
        """Parse pipeline steps for common dirout and unique job names."""

        dirouts = []
        names = []
        for di in self.pipeline:
            for f_config in di.values():
                config = AnalysisConfig(f_config, check_keys=False)
                dirouts.append(config.dirout)

                if 'name' in config:
                    names.append(config.name)

        # The following two exceptions are just extra protection -
        # the user can no longer change "dirout" or "name" as of 2/23/2022.
        # All dirouts should now automatically be set to "./" and job names
        # should be inferrred from the directory name
        if len(set(dirouts)) != 1:
            raise ConfigError('Pipeline steps must have a common output '
                              'directory but received {} different '
                              'directories.'.format(len(set(dirouts))))
        else:
            self.dirout = dirouts[0]

        if len(set(names)) != len(names):
            raise ConfigError('Pipeline steps must have a unique job names '
                              'directory but received {} duplicate names.'
                              .format(len(names) - len(set(names))))

    def _check_dirout_status(self):
        """Check unique status file in dirout."""

        if os.path.exists(self.dirout):
            for fname in os.listdir(self.dirout):
                if (fname.endswith('_status.json')
                        and fname != '{}_status.json'.format(self.name)):
                    msg = ('Cannot run pipeline "{}" in directory '
                           '{}. Another pipeline appears to have '
                           'been run here with status json: {}'
                           .format(self.name, self.dirout, fname))
                    raise PipelineError(msg)

    @property
    def pipeline(self):
        """Get the pipeline steps.

        Returns
        -------
        pipeline : list
            reV pipeline run steps. Should be a list of (command, config)
            pairs.
        """

        return self['pipeline']

    @property
    def logging(self):
        """Get logging kwargs for the pipeline.

        Returns
        -------
        dict
        """
        return self.get('logging', {"log_file": None, "log_level": "INFO"})

    @property
    def hardware(self):
        """Get argument specifying which hardware the pipeline is being run on.

        Defaults to "eagle" (most common use of the reV pipeline)

        Returns
        -------
        hardware : str
            Name of hardware that this pipeline is being run on.
            Defaults to "eagle".
        """
        return self.get('hardware', 'eagle')

    @property
    def status_file(self):
        """Get status file path.

        Returns
        -------
        _status_file : str
            reV status file path.
        """
        if self.dirout is None:
            raise ConfigError('Pipeline has not yet been initialized.')

        return os.path.join(self.dirout, '{}_status.json'.format(self.name))
