#!/usr/bin/env python3
import click
from lollipop import __version__

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

from .generate_mutlist import generate_mutlist
from .deconvolute import deconvolute
from .getmutations_from_basecount import from_basecount


@click.group()
def getmutations():
    """
    Get mutations from a single sample
    """
    pass


getmutations.add_command(from_basecount)


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(__version__)
def cli():
    pass


cli.add_command(generate_mutlist)
cli.add_command(getmutations)
cli.add_command(deconvolute)

if __name__ == "__main__":
    cli()
