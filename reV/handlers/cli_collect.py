# -*- coding: utf-8 -*-
# pylint: disable=all
"""
File collection CLI entry points.
"""
import click
import logging
import os
import pprint
import time

from reV.config.collection import CollectionConfig
from reV.handlers.collection import Collector
from reV.pipeline.status import Status
from reV.utilities import ModuleName
from reV import __version__

from rex.utilities.cli_dtypes import STR, STRLIST, INT
from rex.utilities.hpc import SLURM
from rex.utilities.loggers import init_mult
from rex.utilities.utilities import get_class_properties

logger = logging.getLogger(__name__)


@click.group()
@click.version_option(version=__version__)
@click.option('--name', '-n', default=os.path.basename(os.getcwd()),
              type=STR, show_default=True, help='reV Collect job name.')
@click.option('-v', '--verbose', is_flag=True,
              help='Flag to turn on debug logging. Default is not verbose.')
@click.pass_context
def main(ctx, name, verbose):
    """reV Collect Command Line Interface"""
    ctx.ensure_object(dict)
    ctx.obj['NAME'] = name
    ctx.obj['VERBOSE'] = verbose


@main.command()
def valid_config_keys():
    """
    Echo the valid Collect config keys
    """
    click.echo(', '.join(get_class_properties(CollectionConfig)))


@main.command()
@click.option('--config_file', '-c', required=True,
              type=click.Path(exists=True),
              help='reV collection configuration json file.')
@click.option('-v', '--verbose', is_flag=True,
              help='Flag to turn on debug logging. Default is not verbose.')
@click.pass_context
def from_config(ctx, config_file, verbose):
    """Run reV gen from a config file."""

    # Instantiate the config object
    config = CollectionConfig(config_file)

    # take name from config
    name = ctx.obj['NAME'] = config.name

    # Enforce verbosity if logging level is specified in the config
    if config.log_level == logging.DEBUG:
        verbose = True

    # make output directory if does not exist
    if not os.path.exists(config.dirout):
        os.makedirs(config.dirout)

    # initialize loggers.
    init_mult(name, config.log_directory, modules=[__name__, 'reV', 'rex'],
              verbose=verbose)

    # Initial log statements
    logger.info('Running reV collection from config file: "{}"'
                .format(config_file))
    logger.info('Target output directory: "{}"'.format(config.dirout))
    logger.info('Target logging directory: "{}"'.format(config.log_directory))
    logger.info('Target collection patterns: "{}"'
                .format(config.collect_pattern))
    logger.info('The following project points were specified: "{}"'
                .format(config.get('project_points', None)))
    logger.debug('The full configuration input is as follows:\n{}'
                 .format(pprint.pformat(config, indent=4)))

    # set config objects to be passed through invoke to direct methods
    ctx.obj['LOG_DIR'] = config.log_directory
    ctx.obj['DSETS'] = config.dsets
    ctx.obj['PROJECT_POINTS'] = config.project_points
    ctx.obj['PURGE_CHUNKS'] = config.purge_chunks
    ctx.obj['VERBOSE'] = verbose

    zip_iter = zip(config.fn_out_names, config.collect_pattern)
    for fn_out, collect_pattern in zip_iter:
        ctx.obj['NAME'] = name + '_{}'.format(fn_out.replace('.h5', ''))
        ctx.obj['H5_FILE'] = os.path.join(config.dirout, fn_out)
        ctx.obj['COLLECT_PATTERN'] = collect_pattern

        if config.execution_control.option == 'local':
            status = Status.retrieve_job_status(
                config.dirout,
                module=ModuleName.COLLECT,
                job_name=ctx.obj['NAME']
            )
            if status != 'successful':
                Status.add_job(
                    config.dirout, module=ModuleName.COLLECT,
                    job_name=ctx.obj['NAME'], replace=True,
                    job_attrs={'hardware': 'local',
                               'fout': fn_out,
                               'dirout': config.dirout})
                ctx.invoke(collect)

        elif config.execution_control.option in ('eagle', 'slurm'):
            ctx.invoke(collect_slurm,
                       alloc=config.execution_control.allocation,
                       memory=config.execution_control.memory,
                       walltime=config.execution_control.walltime,
                       feature=config.execution_control.feature,
                       conda_env=config.execution_control.conda_env,
                       module=config.execution_control.module,
                       stdout_path=os.path.join(config.log_directory,
                                                'stdout'),
                       sh_script=config.execution_control.sh_script,
                       verbose=verbose)


@main.group()
@click.option('--h5_file', '-f', required=True, type=click.Path(),
              help='H5 file to be collected into.')
@click.option('--collect_pattern', '-cp', required=True, type=str,
              help='Unix-style /filepath/pattern*.h5 to get a list of '
              'files to collect into h5_file')
@click.option('--project_points', '-pp', type=STR, required=False,
              help='Project points file representing the full '
              'collection scope. This doesnt have to be provided '
              'if points list is to be ignored (collect all data '
              'in h5_files without checking that all gids are there)')
@click.option('--dsets', '-ds', required=True, type=STRLIST,
              help='Dataset names to be collected.')
@click.option('--log_dir', '-ld', type=STR, default='./logs',
              show_default=True,
              help='Directory to put log files.')
@click.option('-p', '--purge_chunks', is_flag=True,
              help='Flag to delete chunked files after collection.')
@click.option('-v', '--verbose', is_flag=True,
              help='Flag to turn on debug logging.')
@click.pass_context
def direct(ctx, h5_file, collect_pattern, project_points, dsets,
           log_dir, purge_chunks, verbose):
    """Main entry point for collection with context passing."""
    ctx.obj['H5_FILE'] = h5_file
    ctx.obj['COLLECT_PATTERN'] = collect_pattern
    ctx.obj['PROJECT_POINTS'] = project_points
    ctx.obj['DSETS'] = dsets
    ctx.obj['LOG_DIR'] = log_dir
    ctx.obj['PURGE_CHUNKS'] = purge_chunks
    ctx.obj['VERBOSE'] = verbose


@direct.command()
@click.option('-v', '--verbose', is_flag=True,
              help='Flag to turn on debug logging.')
@click.pass_context
def collect(ctx, verbose):
    """Run collection on local worker."""

    name = ctx.obj['NAME']
    h5_file = ctx.obj['H5_FILE']
    collect_pattern = ctx.obj['COLLECT_PATTERN']
    project_points = ctx.obj['PROJECT_POINTS']
    dsets = ctx.obj['DSETS']
    log_dir = ctx.obj['LOG_DIR']
    purge_chunks = ctx.obj['PURGE_CHUNKS']
    verbose = any([verbose, ctx.obj['VERBOSE']])

    # initialize loggers for multiple modules
    init_mult(name, log_dir, modules=[__name__, 'reV.handlers.collection'],
              verbose=verbose, node=True)

    for key, val in ctx.obj.items():
        logger.debug('ctx var passed to collection method: "{}" : "{}" '
                     'with type "{}"'.format(key, val, type(val)))

    logger.info('Collection is being run for "{}" with job name "{}" '
                'and collection pattern: {}. Target output path is: {}'
                .format(dsets, name, collect_pattern, h5_file))
    t0 = time.time()

    Collector.collect(h5_file, collect_pattern, project_points, dsets[0])

    if len(dsets) > 1:
        for dset_name in dsets[1:]:
            Collector.add_dataset(h5_file, collect_pattern, dset_name)

    if purge_chunks:
        Collector.purge_chunks(h5_file, collect_pattern, project_points)
    else:
        Collector.move_chunks(h5_file, collect_pattern, project_points)

    runtime = (time.time() - t0) / 60
    logger.info('Collection completed in: {:.2f} min.'.format(runtime))

    # add job to reV status file.
    status = {'dirout': os.path.dirname(h5_file),
              'fout': os.path.basename(h5_file), 'job_status': 'successful',
              'runtime': runtime,
              'finput': collect_pattern}
    Status.make_job_file(
        os.path.dirname(h5_file), ModuleName.COLLECT, name, status
    )


def get_node_cmd(name, h5_file, collect_pattern, project_points, dsets,
                 log_dir='./logs/', purge_chunks=False, verbose=False):
    """Make a reV collection local CLI call string.

    Parameters
    ----------
    name : str
        reV collection jobname.
    h5_file : str
        Path to .h5 file into which data will be collected
    collect_pattern : str
        Unix-style /filepath/pattern*.h5 to get a list of files to
        collect into h5_file
    project_points : str | slice | list | pandas.DataFrame
        Project points that correspond to the full collection of points
        contained in the .h5 files to be collected
    dsets : list
        List of datasets (strings) to be collected.
    log_dir : str
        Log directory.
    purge_chunks : bool
        Flag to delete the chunked files after collection.
    verbose : bool
        Flag to turn on DEBUG logging

    Returns
    -------
    cmd : str
        Single line command line argument to call the following CLI with
        appropriately formatted arguments based on input args:
            python -m reV.handlers.cli_collect [args] collect
    """
    # make a cli arg string for direct() in this module
    args = ['-f {}'.format(SLURM.s(h5_file)),
            '-cp {}'.format(SLURM.s(collect_pattern)),
            '-pp {}'.format(SLURM.s(project_points)),
            '-ds {}'.format(SLURM.s(dsets)),
            '-ld {}'.format(SLURM.s(log_dir)),
            ]

    if purge_chunks:
        args.append('-p')

    if verbose:
        args.append('-v')

    # Python command that will be executed on a node
    # command strings after cli v7.0 use dashes instead of underscores
    cmd = ('python -m reV.handlers.cli_collect -n {} direct {} collect'
           .format(SLURM.s(name), ' '.join(args)))
    logger.debug('Creating the following command line call:\n\t{}'.format(cmd))

    return cmd


@direct.command()
@click.option('--alloc', '-a', required=True, type=STR,
              help='SLURM allocation account name.')
@click.option('--memory', '-mem', default=None, type=INT,
              show_default=True,
              help='SLURM node memory request in GB. Default is None')
@click.option('--walltime', '-wt', default=1.0, type=float,
              show_default=True,
              help='SLURM walltime request in hours. Default is 1.0')
@click.option('--feature', '-l', default=None, type=STR,
              show_default=True,
              help=('Additional flags for SLURM job. Format is "--qos=high" '
                    'or "--depend=[state:job_id]". Default is None.'))
@click.option('--conda_env', '-env', default=None, type=STR,
              show_default=True,
              help='Conda env to activate')
@click.option('--module', '-mod', default=None, type=STR, show_default=True,
              help='Module to load')
@click.option('--stdout_path', '-sout', default='./out/stdout', type=str,
              show_default=True,
              help='Subprocess standard output path. Default is ./out/stdout')
@click.option('--sh_script', '-sh', default=None, type=STR,
              show_default=True,
              help='Extra shell script commands to run before the reV call.')
@click.option('-v', '--verbose', is_flag=True,
              help='Flag to turn on debug logging. Default is not verbose.')
@click.pass_context
def collect_slurm(ctx, alloc, memory, walltime, feature, conda_env, module,
                  stdout_path, sh_script, verbose):
    """Run collection on HPC via SLURM job submission."""

    name = ctx.obj['NAME']
    h5_file = ctx.obj['H5_FILE']
    collect_pattern = ctx.obj['COLLECT_PATTERN']
    log_dir = ctx.obj['LOG_DIR']
    project_points = ctx.obj['PROJECT_POINTS']
    dsets = ctx.obj['DSETS']
    purge_chunks = ctx.obj['PURGE_CHUNKS']
    verbose = any([verbose, ctx.obj['VERBOSE']])

    slurm_manager = ctx.obj.get('SLURM_MANAGER', None)
    if slurm_manager is None:
        slurm_manager = SLURM()
        ctx.obj['SLURM_MANAGER'] = slurm_manager

    cmd = get_node_cmd(name, h5_file, collect_pattern, project_points, dsets,
                       log_dir=log_dir, purge_chunks=purge_chunks,
                       verbose=verbose)
    if sh_script:
        cmd = sh_script + '\n' + cmd

    status = Status.retrieve_job_status(os.path.dirname(h5_file),
                                        module=ModuleName.COLLECT,
                                        job_name=name, hardware='slurm',
                                        subprocess_manager=slurm_manager)

    msg = 'Collect CLI failed to submit jobs!'
    if status == 'successful':
        msg = ('Job "{}" is successful in status json found in "{}", '
               'not re-running.'
               .format(name, os.path.dirname(h5_file)))
    elif 'fail' not in str(status).lower() and status is not None:
        msg = ('Job "{}" was found with status "{}", not resubmitting'
               .format(name, status))
    else:
        logger.info('Running reV collection on SLURM with node name "{}", '
                    'collecting data to "{}" from pattern "{}".'
                    .format(name, h5_file, collect_pattern))
        out = slurm_manager.sbatch(cmd,
                                   alloc=alloc,
                                   memory=memory,
                                   walltime=walltime,
                                   feature=feature,
                                   name=name,
                                   stdout_path=stdout_path,
                                   conda_env=conda_env,
                                   module=module)[0]
        if out:
            msg = ('Kicked off reV collection job "{}" (SLURM jobid #{}).'
                   .format(name, out))

        # add job to reV status file.
        Status.add_job(
            os.path.dirname(h5_file), module=ModuleName.COLLECT,
            job_name=name, replace=True,
            job_attrs={'job_id': out, 'hardware': 'slurm',
                       'fout': os.path.basename(h5_file),
                       'dirout': os.path.dirname(h5_file)})

    click.echo(msg)
    logger.info(msg)


if __name__ == '__main__':
    try:
        main(obj={})
    except Exception:
        logger.exception('Error running reV Collection CLI')
        raise
