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
    LOGGING_FORMAT,
    TABLES,
    DEV_TABLE_POST_DICT,
    PROD_TABLE_POST_DICT,
    ID_MAPPER
)


# Set up logging
def get_logger(    
    LOG_FORMAT     = LOGGING_FORMAT,
    LOG_NAME       = '',
    LOG_FILE_INFO  = 'out.log',
    LOG_FILE_ERROR = 'out.err'):

    log           = logging.getLogger(LOG_NAME)
    log_formatter = logging.Formatter(LOG_FORMAT)

    # Comment this out to suppress console output
    #stream_handler = logging.StreamHandler()
    #stream_handler.setFormatter(log_formatter)
    #log.addHandler(stream_handler)

    # Setup info file
    file_handler_info = logging.FileHandler(LOG_FILE_INFO, mode='w')
    file_handler_info.setFormatter(log_formatter)
    file_handler_info.setLevel(logging.INFO)
    log.addHandler(file_handler_info)

    # Setup error file
    file_handler_error = logging.FileHandler(LOG_FILE_ERROR, mode='w')
    file_handler_error.setFormatter(log_formatter)
    file_handler_error.setLevel(logging.ERROR)
    log.addHandler(file_handler_error)

    log.setLevel(logging.INFO)

    return log

logger = get_logger()


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
            logger.info(response.json())
        else:
            logger.error(response.reason)
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
            logger.info(response.text)
        else:
            logger.error(response.reason)
        status_code = response.status_code
    
    return status_code


def prep_df_data(df, endpoint):
    """
    """
    data = {}
    for col in df.index:
        if df[col] != "NONE":
            data[col] = df[col]

    df["status"] = make_request(endpoint, "POST", data=data)
    return df


def set_unique(df, unique_fields):
    """
    """
    unique_ids = []
    for unique_field in unique_fields:
        unique_ids.append(str(df[unique_field]))

    unique_id = "_".join(unique_ids)
    df["unique_id"] = unique_id
    return df


@click.command()
@click.option("--dir", default=os.getcwd(), help="Full path to the directory containing files to import")
@click.option("--prod", is_flag=True)
def import_dir_data(dir, prod):
    """
    """
    # Get the base path
    if prod:
        TABLE_POST_DICT = PROD_TABLE_POST_DICT
    else:
        TABLE_POST_DICT = DEV_TABLE_POST_DICT

    # Check that the directory exists
    try:
        assert os.path.isdir(dir)
    except AssertionError:
        logger.error("The directory {} does not exist -- please provide the correct directory".format(dir))
        return
    
    # Get all of the files in the directory
    dir_files = os.listdir(dir)

    err_df = pd.DataFrame()

    # Search for files that match the template
    for table in TABLES:
        filename = ".".join([table, "tsv"])
        if filename in dir_files:
            logger.info("IMPORTING {}".format(filename))
            filepath = os.path.join(dir, filename)
            df = pd.read_csv(filepath, sep="\t")

            df = df.fillna('NONE')

            # Post data
            df = df.apply(prep_df_data, axis=1, args=(TABLE_POST_DICT[table],))

            # Get the summary of failed imports
            subset = df.loc[df['status'] != requests.codes.ok]
            if subset.empty:
                continue
            subset = subset.apply(set_unique, axis=1, args=(ID_MAPPER[table],))
            tmp_df = subset[["unique_id", "status"]]
            tmp_df["table"] = table
            
            # Concat to output dataframe
            err_df = pd.concat([err_df, tmp_df], ignore_index=True)

    err_df.to_csv("import_summary.csv", index=False)


if __name__ == '__main__':
    import_dir_data()
    