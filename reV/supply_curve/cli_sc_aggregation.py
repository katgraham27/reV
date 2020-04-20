# -*- coding: utf-8 -*-
"""
reV Supply Curve Aggregation command line interface (cli).
"""
import os
import click
import logging
import pprint
import time
import h5py

from reV.config.supply_curve_configs import SupplyCurveAggregationConfig
from reV.pipeline.status import Status
from reV.supply_curve.tech_mapping import TechMapping
from reV.supply_curve.sc_aggregation import SupplyCurveAggregation

from rex.utilities.execution import SLURM
from rex.utilities.cli_dtypes import STR, INT, FLOAT, FLOATLIST, STRFLOAT
from rex.utilities.loggers import init_mult
from rex.utilities.utilities import dict_str_load

logger = logging.getLogger(__name__)


@click.command()
@click.option('--config_file', '-c', required=True,
              type=click.Path(exists=True),
              help='reV exclusions configuration json file.')
@click.option('-v', '--verbose', is_flag=True,
              help='Flag to turn on debug logging. Default is not verbose.')
@click.pass_context
def from_config(ctx, config_file, verbose):
    """Run reV SC aggregation from a config file."""
    name = ctx.obj['NAME']

    # Instantiate the config object
    config = SupplyCurveAggregationConfig(config_file)

    # take name from config if not default
    if config.name.lower() != 'rev':
        name = config.name

    # Enforce verbosity if logging level is specified in the config
    if config.log_level == logging.DEBUG:
        verbose = True

    # initialize loggers
    init_mult(name, config.logdir, modules=[__name__, 'reV.config',
                                            'reV.utilities'],
              verbose=verbose)

    # Initial log statements
    logger.info('Running reV supply curve aggregation from config '
                'file: "{}"'.format(config_file))
    logger.info('Target output directory: "{}"'.format(config.dirout))
    logger.info('Target logging directory: "{}"'.format(config.logdir))
    logger.debug('The full configuration input is as follows:\n{}'
                 .format(pprint.pformat(config, indent=4)))

    if config.execution_control.option == 'local':
        status = Status.retrieve_job_status(config.dirout,
                                            'supply-curve-aggregation',
                                            name)
        if status != 'successful':
            Status.add_job(
                config.dirout, 'supply-curve-aggregation', name, replace=True,
                job_attrs={'hardware': 'local',
                           'fout': '{}.csv'.format(name),
                           'dirout': config.dirout})
            ctx.invoke(main,
                       name=name,
                       excl_fpath=config.excl_fpath,
                       gen_fpath=config.gen_fpath,
                       res_fpath=config.res_fpath,
                       tm_dset=config.tm_dset,
                       excl_dict=config.excl_dict,
                       check_excl_layers=config.check_excl_layers,
                       res_class_dset=config.res_class_dset,
                       res_class_bins=config.res_class_bins,
                       cf_dset=config.cf_dset,
                       lcoe_dset=config.lcoe_dset,
                       data_layers=config.data_layers,
                       resolution=config.resolution,
                       power_density=config.power_density,
                       area_filter_kernel=config.area_filter_kernel,
                       min_area=config.min_area,
                       friction_fpath=config.friction_fpath,
                       friction_dset=config.friction_dset,
                       out_dir=config.dirout,
                       log_dir=config.logdir,
                       verbose=verbose)

    elif config.execution_control.option in ('eagle', 'slurm'):

        ctx.obj['NAME'] = name
        ctx.obj['EXCL_FPATH'] = config.excl_fpath
        ctx.obj['GEN_FPATH'] = config.gen_fpath
        ctx.obj['RES_FPATH'] = config.res_fpath
        ctx.obj['TM_DSET'] = config.tm_dset
        ctx.obj['EXCL_DICT'] = config.excl_dict
        ctx.obj['CHECK_LAYERS'] = config.check_excl_layers
        ctx.obj['RES_CLASS_DSET'] = config.res_class_dset
        ctx.obj['RES_CLASS_BINS'] = config.res_class_bins
        ctx.obj['CF_DSET'] = config.cf_dset
        ctx.obj['LCOE_DSET'] = config.lcoe_dset
        ctx.obj['DATA_LAYERS'] = config.data_layers
        ctx.obj['RESOLUTION'] = config.resolution
        ctx.obj['POWER_DENSITY'] = config.power_density
        ctx.obj['AREA_FILTER_KERNEL'] = config.area_filter_kernel
        ctx.obj['MIN_AREA'] = config.min_area
        ctx.obj['FRICTION_FPATH'] = config.friction_fpath
        ctx.obj['FRICTION_DSET'] = config.friction_dset
        ctx.obj['OUT_DIR'] = config.dirout
        ctx.obj['LOG_DIR'] = config.logdir
        ctx.obj['VERBOSE'] = verbose

        ctx.invoke(slurm,
                   alloc=config.execution_control.alloc,
                   memory=config.execution_control.node_mem,
                   feature=config.execution_control.feature,
                   walltime=config.execution_control.walltime,
                   conda_env=config.execution_control.conda_env,
                   module=config.execution_control.module)


@click.group(invoke_without_command=True)
@click.option('--name', '-n', default='agg', type=STR,
              help='Job name. Default is "agg".')
@click.option('--excl_fpath', '-ef', type=STR, required=True,
              help='Exclusions file (.h5).')
@click.option('--gen_fpath', '-gf', type=STR, required=True,
              help='reV generation/econ output file.')
@click.option('--res_fpath', '-rf', type=STR, default=None,
              help='Resource file, required if techmap dset is to be created.')
@click.option('--tm_dset', '-tm', type=STR, required=True,
              help='Dataset in the exclusions file that maps the exclusions '
              'to the resource being analyzed.')
@click.option('--excl_dict', '-exd', type=STR, default=None,
              help='String representation of a dictionary of exclusion '
              'LayerMask arguments {layer: {kwarg: value}} where layer is a '
              'dataset in excl_fpath and kwarg can be "inclusion_range", '
              '"exclude_values", "include_values", "use_as_weights", '
              '"exclude_nodata", and/or "weight".')
@click.option('--check_excl_layers', '-cl', is_flag=True,
              help=('run a pre-flight check on each exclusion layer to '
                    'ensure they contain un-excluded values'))
@click.option('--res_class_dset', '-cd', type=STR, default=None,
              help='Dataset to determine the resource class '
              '(must be in gen_fpath).')
@click.option('--res_class_bins', '-cb', type=FLOATLIST, default=None,
              help='List of resource class bin edges.')
@click.option('--cf_dset', '-cf', type=STR, default='cf_mean-means',
              help='Dataset containing capacity factor values to aggregate.')
@click.option('--lcoe_dset', '-lc', type=STR, default='lcoe_fcr-means',
              help='Dataset containing lcoe values to aggregate.')
@click.option('--data_layers', '-d', type=STR, default=None,
              help='String representation of a dictionary of additional data '
              'layers to include in the aggregation e.g. '
              '{"slope": {"dset": "srtm_slope", "method": "mean"}}')
@click.option('--resolution', '-r', type=INT, default=64,
              help='Number of exclusion points along a squares edge to '
              'include in an aggregated supply curve point.')
@click.option('--power_density', '-pd', type=STRFLOAT, default=None,
              help='Power density in MW/km2 or filepath to variable power '
              'density csv file. None will attempt to infer a constant '
              'power density from the generation meta data technology.')
@click.option('--area_filter_kernel', '-afk', type=STR, default='queen',
              help='Contiguous area filter kernel name ("queen", "rook").')
@click.option('--min_area', '-ma', type=FLOAT, default=None,
              help='Contiguous area filter minimum area, default is None '
              '(No minimum area filter).')
@click.option('--friction_fpath', '-ff', type=STR, default=None,
              help='Optional h5 filepath to friction surface data. '
              'Must match the exclusion shape/resolution and be '
              'paired with the --friction_dset input arg.')
@click.option('--friction_dset', '-fd', type=STR, default=None,
              help='Optional friction surface dataset in friction_fpath.')
@click.option('--out_dir', '-o', type=STR, default='./',
              help='Directory to save aggregation summary output.')
@click.option('--log_dir', '-ld', type=STR, default='./logs/',
              help='Directory to save aggregation logs.')
@click.option('-v', '--verbose', is_flag=True,
              help='Flag to turn on debug logging. Default is not verbose.')
@click.pass_context
def main(ctx, name, excl_fpath, gen_fpath, res_fpath, tm_dset, excl_dict,
         check_excl_layers, res_class_dset, res_class_bins, cf_dset, lcoe_dset,
         data_layers, resolution, power_density, area_filter_kernel, min_area,
         friction_fpath, friction_dset, out_dir, log_dir, verbose):
    """reV Supply Curve Aggregation Summary CLI."""

    ctx.ensure_object(dict)
    ctx.obj['NAME'] = name
    ctx.obj['EXCL_FPATH'] = excl_fpath
    ctx.obj['GEN_FPATH'] = gen_fpath
    ctx.obj['RES_FPATH'] = res_fpath
    ctx.obj['TM_DSET'] = tm_dset
    ctx.obj['EXCL_DICT'] = excl_dict
    ctx.obj['CHECK_LAYERS'] = check_excl_layers
    ctx.obj['RES_CLASS_DSET'] = res_class_dset
    ctx.obj['RES_CLASS_BINS'] = res_class_bins
    ctx.obj['CF_DSET'] = cf_dset
    ctx.obj['LCOE_DSET'] = lcoe_dset
    ctx.obj['DATA_LAYERS'] = data_layers
    ctx.obj['RESOLUTION'] = resolution
    ctx.obj['POWER_DENSITY'] = power_density
    ctx.obj['AREA_FILTER_KERNEL'] = area_filter_kernel
    ctx.obj['MIN_AREA'] = min_area
    ctx.obj['FRICTION_FPATH'] = friction_fpath
    ctx.obj['FRICTION_DSET'] = friction_dset
    ctx.obj['OUT_DIR'] = out_dir
    ctx.obj['LOG_DIR'] = log_dir
    ctx.obj['VERBOSE'] = verbose

    if ctx.invoked_subcommand is None:
        t0 = time.time()
        init_mult(name, log_dir, modules=[__name__, 'reV.supply_curve'],
                  verbose=verbose)

        with h5py.File(excl_fpath, mode='r') as f:
            dsets = list(f)
        if tm_dset not in dsets:
            try:
                TechMapping.run(excl_fpath, res_fpath, tm_dset)
            except Exception as e:
                logger.exception('TechMapping process failed. Received the '
                                 'following error:\n{}'.format(e))
                raise e

        if isinstance(excl_dict, str):
            excl_dict = dict_str_load(excl_dict)

        if isinstance(data_layers, str):
            data_layers = dict_str_load(data_layers)

        try:
            summary = SupplyCurveAggregation.summary(
                excl_fpath, gen_fpath, tm_dset,
                excl_dict=excl_dict,
                res_class_dset=res_class_dset,
                res_class_bins=res_class_bins,
                cf_dset=cf_dset,
                lcoe_dset=lcoe_dset,
                data_layers=data_layers,
                resolution=resolution,
                power_density=power_density,
                area_filter_kernel=area_filter_kernel,
                min_area=min_area,
                friction_fpath=friction_fpath,
                friction_dset=friction_dset,
                check_excl_layers=check_excl_layers)

        except Exception as e:
            logger.exception('Supply curve Aggregation failed. Received the '
                             'following error:\n{}'.format(e))
            raise e

        fn_out = '{}.csv'.format(name)
        fpath_out = os.path.join(out_dir, fn_out)
        summary.to_csv(fpath_out)

        runtime = (time.time() - t0) / 60
        logger.info('Supply curve aggregation complete. '
                    'Time elapsed: {:.2f} min. Target output dir: {}'
                    .format(runtime, out_dir))

        finput = [excl_fpath, gen_fpath]
        if res_fpath is not None:
            finput.append(res_fpath)

        # add job to reV status file.
        status = {'dirout': out_dir, 'fout': fn_out,
                  'job_status': 'successful',
                  'runtime': runtime,
                  'finput': finput}
        Status.make_job_file(out_dir, 'supply-curve-aggregation', name, status)


def get_node_cmd(name, excl_fpath, gen_fpath, res_fpath, tm_dset, excl_dict,
                 check_excl_layers, res_class_dset, res_class_bins, cf_dset,
                 lcoe_dset, data_layers, resolution, power_density,
                 area_filter_kernel, min_area, friction_fpath, friction_dset,
                 out_dir, log_dir, verbose):
    """Get a CLI call command for the SC aggregation cli."""

    args = ('-n {name} '
            '-ef {excl_fpath} '
            '-gf {gen_fpath} '
            '-rf {res_fpath} '
            '-tm {tm_dset} '
            '-exd {excl_dict} '
            '-cd {res_class_dset} '
            '-cb {res_class_bins} '
            '-cf {cf_dset} '
            '-lc {lcoe_dset} '
            '-d {data_layers} '
            '-r {resolution} '
            '-pd {power_density} '
            '-afk {area_filter_kernel} '
            '-ma {min_area} '
            '-ff {friction_fpath} '
            '-fd {friction_dset} '
            '-o {out_dir} '
            '-ld {log_dir} '
            )

    args = args.format(name=SLURM.s(name),
                       excl_fpath=SLURM.s(excl_fpath),
                       gen_fpath=SLURM.s(gen_fpath),
                       res_fpath=SLURM.s(res_fpath),
                       tm_dset=SLURM.s(tm_dset),
                       excl_dict=SLURM.s(excl_dict),
                       res_class_dset=SLURM.s(res_class_dset),
                       res_class_bins=SLURM.s(res_class_bins),
                       cf_dset=SLURM.s(cf_dset),
                       lcoe_dset=SLURM.s(lcoe_dset),
                       data_layers=SLURM.s(data_layers),
                       resolution=SLURM.s(resolution),
                       power_density=SLURM.s(power_density),
                       area_filter_kernel=SLURM.s(area_filter_kernel),
                       min_area=SLURM.s(min_area),
                       friction_fpath=SLURM.s(friction_fpath),
                       friction_dset=SLURM.s(friction_dset),
                       out_dir=SLURM.s(out_dir),
                       log_dir=SLURM.s(log_dir),
                       )

    if check_excl_layers:
        args += '-cl '

    if verbose:
        args += '-v '

    cmd = 'python -m reV.supply_curve.cli_sc_aggregation {}'.format(args)
    return cmd


@main.command()
@click.option('--alloc', '-a', required=True, type=STR,
              help='SLURM allocation account name.')
@click.option('--walltime', '-wt', default=1.0, type=float,
              help='SLURM walltime request in hours. Default is 1.0')
@click.option('--feature', '-l', default=None, type=STR,
              help=('Additional flags for SLURM job. Format is "--qos=high" '
                    'or "--depend=[state:job_id]". Default is None.'))
@click.option('--memory', '-mem', default=None, type=INT,
              help='SLURM node memory request in GB. Default is None')
@click.option('--module', '-mod', default=None, type=STR,
              help='Module to load')
@click.option('--conda_env', '-env', default=None, type=STR,
              help='Conda env to activate')
@click.option('--stdout_path', '-sout', default=None, type=STR,
              help='Subprocess standard output path. Default is in out_dir.')
@click.pass_context
def slurm(ctx, alloc, walltime, feature, memory, module, conda_env,
          stdout_path):
    """slurm (Eagle) submission tool for reV supply curve aggregation."""

    name = ctx.obj['NAME']
    excl_fpath = ctx.obj['EXCL_FPATH']
    gen_fpath = ctx.obj['GEN_FPATH']
    res_fpath = ctx.obj['RES_FPATH']
    tm_dset = ctx.obj['TM_DSET']
    excl_dict = ctx.obj['EXCL_DICT']
    check_excl_layers = ctx.obj['CHECK_LAYERS']
    res_class_dset = ctx.obj['RES_CLASS_DSET']
    res_class_bins = ctx.obj['RES_CLASS_BINS']
    cf_dset = ctx.obj['CF_DSET']
    lcoe_dset = ctx.obj['LCOE_DSET']
    data_layers = ctx.obj['DATA_LAYERS']
    resolution = ctx.obj['RESOLUTION']
    power_density = ctx.obj['POWER_DENSITY']
    area_filter_kernel = ctx.obj['AREA_FILTER_KERNEL']
    min_area = ctx.obj['MIN_AREA']
    friction_fpath = ctx.obj['FRICTION_FPATH']
    friction_dset = ctx.obj['FRICTION_DSET']
    out_dir = ctx.obj['OUT_DIR']
    log_dir = ctx.obj['LOG_DIR']
    verbose = ctx.obj['VERBOSE']

    if stdout_path is None:
        stdout_path = os.path.join(log_dir, 'stdout/')

    cmd = get_node_cmd(name, excl_fpath, gen_fpath, res_fpath,
                       tm_dset, excl_dict, check_excl_layers,
                       res_class_dset, res_class_bins,
                       cf_dset, lcoe_dset, data_layers, resolution,
                       power_density, area_filter_kernel, min_area,
                       friction_fpath, friction_dset,
                       out_dir, log_dir, verbose)

    status = Status.retrieve_job_status(out_dir, 'supply-curve-aggregation',
                                        name)
    if status == 'successful':
        msg = ('Job "{}" is successful in status json found in "{}", '
               'not re-running.'
               .format(name, out_dir))
    else:
        logger.info('Running reV SC aggregation on SLURM with '
                    'node name "{}"'.format(name))
        slurm = SLURM(cmd, alloc=alloc, memory=memory,
                      walltime=walltime, feature=feature,
                      name=name, stdout_path=stdout_path,
                      conda_env=conda_env, module=module)
        if slurm.id:
            msg = ('Kicked off reV SC aggregation job "{}" '
                   '(SLURM jobid #{}).'
                   .format(name, slurm.id))
            Status.add_job(
                out_dir, 'supply-curve-aggregation', name, replace=True,
                job_attrs={'job_id': slurm.id, 'hardware': 'eagle',
                           'fout': '{}.csv'.format(name), 'dirout': out_dir})
        else:
            msg = ('Was unable to kick off reV SC job "{}". '
                   'Please see the stdout error messages'
                   .format(name))
    click.echo(msg)
    logger.info(msg)


if __name__ == '__main__':
    try:
        main(obj={})
    except Exception:
        logger.exception('Error running reV SC aggregation CLI.')
        raise