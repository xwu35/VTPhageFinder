#!/usr/bin/env python3

import sys
import os
import subprocess
import click
from datetime import datetime

version = "1.0.0"
@click.version_option(version, "--version", "-v")

@click.command(
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=150),
    help='Usage:\n python VTPhageFinder --reads_dir <reads directory> --sample_info <sample information table> --output_dir <output directory> '
    '--reference_genome <reference genome sequences> --prophage_region <Prophage coordinates>'
)
@click.option(
    '--reads_dir',
    required=True,
    type=click.Path(dir_okay=True, exists=True, resolve_path=True),
    help='Reads directory'
)
@click.option(
    '--sample_info',
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
    help='Sample information table (tab separated).' 
    ' The table must contain three columns (sample, R1, R2)'
)
@click.option(
    '--output_dir',
    default="VTPhageFinder_OUTPUT",
    type=click.Path(dir_okay=True, resolve_path=True),
    show_default=True,
    help=('Output directory')
)
@click.option(
    '--reference_genome',
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
    help=('Reference genome sequences')
)
@click.option(
    '--prophage_region',
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
    help=('Prophage coordinates on the host genome')
)
@click.option(
    '--prophage_identity',
    default=100,
    type=int,
    show_default=True,
    help=('Percent of identity to remove reads mapped to the prophage region'),
)
@click.option(
    '--nonprophage_identity',
    default=98,
    type=int,
    show_default=True,
    help=('Percent of identity to remove reads mapped to the non-prophage region')
)
@click.option(
    '--mapper',
    default="minimap2",
    type=str,
    show_default=True,
    help=('Mapping software; available options are: minimap2, bowtie2')
)
@click.option(
    '--step',
    default='assemble',
    type=str,
    show_default=True,
    help=('Steps to run; available options are: fastqc, preprocess, assemble')
)
@click.option(
    '--dryrun',
    is_flag=True,
    default=False,
    show_default=True,
    help='Check rules to run and files to produce'
)
@click.option(
    '--conda_envs',
    default='',
    show_default=True,
    help='Directory to store conda environments.'
    ' By default, the ".snakemake" directory relative to the invocation directory is used'
)
@click.option(
    '--profile',
    default='slurm',
    show_default=True,
    help='Snakemake profile for cluster execution'
)

def run_vtphagefinder(reads_dir, sample_info, output_dir, reference_genome, 
          prophage_region, prophage_identity, nonprophage_identity,
          mapper, step, dryrun, conda_envs, profile):

          # write run log if it is not a dry run
          if not dryrun:
            os.makedirs(output_dir, exist_ok=True)
            logfile = os.path.join(output_dir, f"{os.path.basename(output_dir)}_run.log")
            with open(logfile, "w") as log:
                log.write("================VTPhageFinder run log==============\n")
                log.write(f"Start time: {datetime.now()}\n")
                log.write(f"VTPhageFinder version: {version}\n")
                log.write(f"Raw reads direcotry: {reads_dir}\n")
                log.write(f"Sample table: {sample_info}\n")
                log.write(f"Results directory: {output_dir}\n")
                log.write(f"Host genome: {reference_genome}\n")
                log.write(f"Host prophage coordinates: {prophage_region}\n")
                log.write(f"Host prophage region identity cutoff: {prophage_identity}\n")
                log.write(f"Host non-prophage region identity cutoff: {nonprophage_identity}\n")
                log.write(f"Mappine software: {mapper}")
          
          cmd = (
            'snakemake --snakefile workflow/Snakefile '
                '--use-conda --conda-frontend mamba '
                '{conda_envs} '
                '--profile {profile} --rerun-incomplete ' 
                '--printshellcmds --nolock --show-failed-logs '
                '{dryrun} '
                '--config reads_dir={reads} sample_info={meta} '
                    'results_dir={results} host_genome={host} '
                    'host_prophage_region={region} prophage_identity={propid} '
                    'nonprophage_identity={nonpropid} '
                    'mapper={mapper} step={step}'
          ).format(
            conda_envs='' if conda_envs=='' else '--conda-prefix {}'.format(conda_envs),
            profile=profile,
            dryrun='--dryrun' if dryrun else '',
            reads=reads_dir,
            meta=sample_info,
            results=output_dir,
            host=reference_genome,
            region=prophage_region,
            propid=prophage_identity,
            nonpropid=nonprophage_identity,
            mapper=mapper,
            step=step
          )

          # run snakemake with command-line config
          try:
            subprocess.run(cmd, check=True, shell=True)
          except subprocess.CalledProcessError:
            print("Snakemake failed. see log for details.", file=sys.stderr)
            sys.exit(1)

if __name__ == "__main__":
  run_vtphagefinder()
