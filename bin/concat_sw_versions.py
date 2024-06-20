#!/usr/bin/env python
"""Concatinate software versions."""

import click
import yaml
from yaml import Loader
import pandas as pd
from pathlib import Path


def get_versions(version_obj: dict[str: dict]) -> dict[str, str]:
    workflow_name = list(version_obj.keys())[0].split(":")[-1]
    raw_softwares = list(version_obj.values())[0]
    # add workflow name to the list of all softwares
    softwares = {}
    for sw, version_info in raw_softwares.items():
        version_info["workflow"] = workflow_name
        softwares[sw] = version_info
        # get container
        if "http" not in version_info["container"]:
            version_info["container"] = None
    return softwares


@click.command()
@click.option("-o", "--output", type=click.File("w"), help="Path to write output file to.")
@click.argument("version_files", nargs=-1)
def cli(output, version_files):
    """Concatinate the versions of softwares."""

    all_versions = {}
    for file in version_files:
        with open(file) as vfile:
            sw_version = get_versions(yaml.load(vfile, Loader=Loader))
            # combine new sw versions with existing sw versions
            all_versions = {**all_versions, **sw_version}

    # convert version dict to csv tables
    df = (pd.DataFrame
        .from_dict(all_versions, orient="index")
        .drop("workflow", axis=1)
        .fillna("-")
    )
    df.index.name = "software"
    df.reset_index(inplace=True)
    df.sort_values("software", inplace=True)
    df.columns = [col.capitalize() for col in df.columns]
    # export to csv
    df.to_csv(output, index=False)
    click.secho(f"Wrote output file: {output.name}", fg="green")


if __name__ == "__main__":
    cli()