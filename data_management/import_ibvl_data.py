import os
import sys
import click
import json
import logging
import requests
import pandas as pd
import numpy as np

from pprint import pprint

from constants import (
    BASE_URL,
    LOGGING_FORMAT,
    TABLES,
    TABLE_POST_DICT
)


# Set up logging
logging.basicConfig(format=LOGGING_FORMAT, stream=sys.stderr, level=logging.INFO)


def make_request(request_url, request_type, data=None):
    """
    """
    if request_type == "GET":
        response = requests.get(
            request_url,
            auth=(
                os.environ.get("ORACLE_USERNAME"),
                os.environ.get("ORACLE_PASSWORD")
            )
        )

        if response.status_code == requests.codes.ok:
            logging.info(response.json())
        else:
            logging.error(response.reason)
        status_code = response.status_code
    elif request_type == "POST":
        response = requests.post(
            request_url,
            auth=(
                os.environ.get("ORACLE_USERNAME"),
                os.environ.get("ORACLE_PASSWORD")
            ),
            data=json.dumps(data),
            headers={
                'Accept': 'application/json',
                'Content-Type': 'application/json',
            }
        )
        if response.status_code == requests.codes.ok:
            logging.info(response.text)
        else:
            logging.error(response.reason)
        status_code = response.status_code
    
    return status_code


def prep_df_data(df, endpoint):
    """
    """
    data = {}
    for col in df.index:
        if df[col] != "NONE":
            data[col] = df[col]

    status = make_request(endpoint, "POST", data=data)


@click.command()
@click.option("--dir", default=os.getcwd(), help="Full path to the directory containing files to import")
def import_dir_data(dir):
    """
    """
    # Check that the directory exists
    try:
        assert os.path.isdir(dir)
    except AssertionError:
        logging.error("The directory {} does not exist -- please provide the correct directory".format(dir))
        return
    
    # Get all of the files in the directory
    dir_files = os.listdir(dir)

    # Search for files that match the template
    for table in TABLES:
        filename = ".".join([table, "tsv"])
        if filename in dir_files:
            filepath = os.path.join(dir, filename)
            df = pd.read_csv(filepath, sep="\t")

            df = df.fillna('NONE')

            # Post data
            df.apply(prep_df_data, axis=1, args=(TABLE_POST_DICT[table],))


if __name__ == '__main__':
    import_dir_data()
    