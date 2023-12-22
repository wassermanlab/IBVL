import json
from sqlalchemy import (
    create_engine,
    text,
    MetaData,
    Table,
    Column,
    Integer,
    String,
    delete,
    Float,
    UniqueConstraint,
    ForeignKeyConstraint,
)
from sqlalchemy.exc import DataError, IntegrityError, ProgrammingError
import pandas as pd
import numpy as np
import logging
import warnings
import sys
import os
from dotenv import load_dotenv
from datetime import datetime
from sqlalchemy import text

load_dotenv()

# get command line arguments
rootDir = os.environ.get("ORACLE_TABLE_PATH")
chunk_size = int(os.environ.get("CHUNK_SIZE"))
verbose = os.environ.get("VERBOSE") == "true"
dbConnectionString = os.environ.get("DB")

if rootDir == None:
    print("No root directory specified")
    exit()

data_issue_logger = logging.getLogger("data_issues")
data_issue_logger.setLevel(logging.WARNING)
data_issue_handler = logging.FileHandler(f"data_issues-{str(datetime.now())}.log")
data_issue_handler.setLevel(logging.WARNING)
data_issue_logger.addHandler(data_issue_handler)


print(rootDir)

if False and verbose:
    engine = create_engine(
        dbConnectionString,
        echo=True,
        pool_pre_ping=True,
        pool_recycle=3600,
        connect_args={"autocommit": True},
    )
else:
    engine = create_engine(
        dbConnectionString,
        pool_pre_ping=True,
        pool_recycle=3600,
        connect_args={"autocommit": True},
    )

mydb = engine.connect()

metadata = MetaData()

use_n_chromosomes = 0


def transform_header(header):
    # Example: Convert all column names to uppercase
    return [col.upper() for col in header]


def readTSV(file, info, dtype={}):
#    df = pd.read_csv(file, sep=info["separator"], dtype=dtype, na_values=["NA"], keep_default_na=False)
    df = pd.read_csv(file, sep=info["separator"], dtype=dtype)
    df.rename(columns={"All_info$variant": "variant"}, inplace=True)
    df.columns = [col.upper() for col in df.columns]
    return df


def inspectTSV(file):
    total_rows = 0
    separator = "\t"
    small_read = pd.read_csv(file, sep="\t", nrows=4, header=None)

    num_columns = small_read.shape[1]
    if num_columns == 1 and "gene" not in file:
        separator = " "
        small_read = pd.read_csv(file, sep=" ", nrows=4, header=None)
        if small_read.shape[1] == 1:
            print("Could not determine separator for " + file)
            quit()

    columns = [col.upper() for col in small_read.values.tolist()[0]]

    for chunk in pd.read_csv(file, sep="\t", chunksize=chunk_size):
        total_rows += len(chunk)

    return {
        "total_rows": total_rows,
        "num_columns": num_columns,
        "columns": columns,
        "types": {},
        "separator": separator,
    }


pk_maps = {}
next_id_maps = {}

try:
    # Load from file
    with open("pk_maps.json", "r") as f:
        pk_maps = json.load(f)

    with open("next_id_maps.json", "r") as f:
        next_id_maps = json.load(f)
except FileNotFoundError:
    print("No pk_maps.json file found. Starting from scratch.")
    pass

tables = {}


# TODO: chunkify the insert operation: catch exception on bulk, then break into smaller chunks and iterate
# add persisting maps
# add boolean to skip emptying table before importing
def persist_maps():
    # Save to file
    with open("pk_maps.json", "w") as f:
        json.dump(pk_maps, f)

    with open("next_id_maps.json", "w") as f:
        json.dump(next_id_maps, f)


def chunks(l, n):
    """Yield successive n-sized chunks from list l."""
    for i in range(0, len(l), n):
        yield l[i : i + n]


def resolve_PK(referencedModel, name):
    try:
        result = pk_maps[referencedModel][name]
        return result
    except KeyError:
        #        print("Could not find PK in "+referencedModel+": "+str(name))
        # query db?
        return None


def import_chr_into_table(file, file_info, action_info):
    name = action_info.get("name")
    fk_map = action_info.get("fk_map")
    pk_lookup_col = action_info.get("pk_lookup_col")
    filters = action_info.get("filters") or {}

    missingRefCount = 0
    table = tables.get(name)
    if table is None:
        columns = file_info["columns"]
        if "DO_COMPOUND_FK" in fk_map:
            columns.remove("VARIANT")
            columns.remove("TRANSCRIPT")
            columns.append("VARIANT_TRANSCRIPT")
        table = Table(name.upper(), metadata, autoload_with=engine)
        tables[name] = table

    # for column in table.columns:
    #    print(f"{name} Column Name: {column.name}, Column Type: {column.type}")

    #    warnings.filterwarnings("ignore")
    types_dict = {}
    for column in table.columns:
        # convert sql types to pandas types
        if isinstance(column.type, Integer):
            types_dict[column.name] = "Int64"
        elif isinstance(column.type, String):
            types_dict[column.name] = "str"
        elif isinstance(column.type, Float):
            types_dict[column.name] = "float64"
        else:
            pass

    df = readTSV(file, file_info, dtype=types_dict)
    df.replace(np.nan, None, inplace=True)
#    df = df.where(pd.notnull(df), None)
#    df.fillna(None)
    #   warnings.resetwarnings()
    data_list = []
    for index, row in df.iterrows():
        data = row.to_dict()
        pk = next_id_maps[name] + index
        data["ID"] = pk

        skip = False
        for col, filter in filters.items():
            data[col] = filter(data[col])
        for fk_col, fk_model in fk_map.items():
            map_key = None
            resolved_pk = None
            debug_row = None
            if fk_col == "DO_COMPOUND_FK":
                debug_row = data.copy()
                v_id = resolve_PK("variants", data["VARIANT"])
                t_id = resolve_PK("transcripts", data["TRANSCRIPT"])
                map_key = "-".join([str(t_id), str(v_id)])
                del data["VARIANT"]
                del data["TRANSCRIPT"]
                resolved_pk = resolve_PK("variants_transcripts", map_key)
                fk_col = "VARIANT_TRANSCRIPT"
            else:
                map_key = data[fk_col]
                if map_key == "NA":
                    data[fk_col] = None
                else:
                    resolved_pk = resolve_PK(fk_model, map_key)
            if map_key is not None and verbose:
                print(
                    "resolved "
                    + fk_model
                    + "."
                    + data[fk_col]
                    + " to "
                    + str(resolved_pk)
                )
            if resolved_pk is not None:
                data[fk_col] = resolved_pk
            else:
                # MOVE TO error.log
                data_issue_logger.warning(
                    "Could not resolve using map key "
                    + map_key
                    + " for "
                    + fk_model
                    + "."
                    + fk_col
                    + " in "
                )
                if (debug_row is not None):
                    data_issue_logger.warning(debug_row)
                else:
                    data_issue_logger.warning(data)
                missingRefCount += 1
                skip = True
        if skip:
            continue
        data_list.append(data)

    with engine.connect() as connection:
        successCount = 0
        failCount = 0
        duplicateCount = 0
        successful_chunks = 0
        fail_chunks = 0

        for chunk in chunks(data_list, chunk_size):
            try:
                connection.execute(table.insert(), chunk)
                # chunk worked
                successful_chunks += 1
                successCount += len(chunk)
            except Exception as e:
                #                print(e)
                fail_chunks += 1
                for row in chunk:
                    try:
                        connection.execute(table.insert(), row)
                        successCount += 1

                    except DataError as e:
                        data_issue_logger.warn(e)
                        failCount += 1
                        quit()
                    except IntegrityError as e:
                        msg = str(e)
                        if "Duplicate" in msg:
                            duplicateCount += 1
                            successCount += 1
                        else:
                            failCount += 1
                            data_issue_logger.warn(e)
#                            quit()
                    except Exception as e:
                        data_issue_logger.warn(e)
                        failCount += 1

            if pk_lookup_col is not None:
                pk_map = {}
                for data in chunk:
                    # record the PKS for each row that was added
                    if isinstance(pk_lookup_col, list):
#                        print(pk_lookup_col)
#                        print(data)
                        map_key = "-".join([str(data[col]) for col in pk_lookup_col])
                    elif isinstance(pk_lookup_col, str):
                        map_key = data[pk_lookup_col]
                    pk_map[map_key] = data["ID"]
                    if verbose:
                        print("added " + name + "." + map_key + " to pk map")
                pk_maps[name].update(pk_map)

    next_id_maps[name] += file_info["total_rows"]

    return {
        "success": successCount,
        "fail": failCount,
        "missingRef": missingRefCount,
        "duplicate": duplicateCount,
        "successful_chunks": successful_chunks,
        "fail_chunks": fail_chunks,
    }


# map of import functions
model_import_actions = {
    "genes": {
        "name": "genes",
        "pk_lookup_col": "SHORT_NAME",
        "fk_map": {},
        "empty_first": True
        #        "skip":True
    },
    "transcripts": {
        "name": "transcripts",
        "pk_lookup_col": "TRANSCRIPT_ID",
        "fk_map": {"GENE": "genes"},
        "empty_first": True,
        "filters": {
            "TRANSCRIPT_TYPE": lambda x: x.replace("RefSeq", "R"),
        },
    },
    "variants": {
        "name": "variants",
        "pk_lookup_col": "VARIANT_ID",
        "fk_map": {},
        "empty_first": True
        #        "skip":True
    },
    "variants_transcripts": {
        "name": "variants_transcripts",
        "pk_lookup_col": ["TRANSCRIPT", "VARIANT"],
        "fk_map": {"TRANSCRIPT": "transcripts", "VARIANT": "variants"},
        "empty_first": True
        #        "skip":True
    },
    "variants_annotations": {
        "name": "variants_annotations",
        "pk_lookup_col": None,
        "fk_map": {"DO_COMPOUND_FK": "for variants_transcripts"},
        "empty_first": True,
        "filters":{
            "HGVSP": lambda x: x.replace("%3D","=")
        }
    },
    "variants_consequences": {
        "name": "variants_consequences",
        "pk_lookup_col": None,
        "fk_map": {"DO_COMPOUND_FK": "for variants_transcripts"},
        "empty_first": True,
    },
    "sv_consequences": {
        "name": "sv_consequences",
        "pk_lookup_col": None,
        "fk_map": {"GENE": "genes", "VARIANT": "variants"},
        "empty_first": True,
    },
    "snvs": {
        "name": "snvs",
        "pk_lookup_col": None,
        "fk_map": {"VARIANT": "variants"},
        "empty_first": True,
    },
    "svs": {
        "name": "svs",
        "pk_lookup_col": None,
        "fk_map": {"VARIANT": "variants"},
        "empty_first": True,
    },
    "svs_ctx": {
        "name": "svs_ctx",
        "pk_lookup_col": None,
        "fk_map": {"VARIANT": "variants"},
        "empty_first": True,
    },
    "str": {
        "name": "str",
        "pk_lookup_col": None,
        "fk_map": {"VARIANT": "variants"},
        "empty_first": True,
    },
    "mts": {
        "name": "mts",
        "pk_lookup_col": None,
        "fk_map": {"VARIANT": "variants"},
        "empty_first": True,
    },
    "genomic_ibvl_frequencies": {
        "name": "genomic_ibvl_frequencies",
        "pk_lookup_col": None,
        "fk_map": {"VARIANT": "variants"},
        "empty_first": True,
    },
    "genomic_gnomad_frequencies": {
        "name": "genomic_gnomad_frequencies",
        "pk_lookup_col": None,
        "fk_map": {"VARIANT": "variants"},
        "empty_first": True,
    },
    "mt_ibvl_frequencies": {
        "name": "mt_ibvl_frequencies",
        "pk_lookup_col": None,
        "fk_map": {"VARIANT": "variants"},
        "empty_first": True,
    },
    "mt_gnomad_frequencies": {
        "name": "mt_gnomad_frequencies",
        "pk_lookup_col": None,
        "fk_map": {"VARIANT": "variants"},
        "empty_first": True,
    },
}


def report_counts(counts):
    percent_success = "N/A"
    if counts["successCountTotal"] + counts["failCountTotal"] > 0:
        percent_success = (
            100
            * counts["successCountTotal"]
            / (counts["successCountTotal"] + counts["failCountTotal"])
        )
    print(
        str(percent_success)
        + "% success. Count: "
        + str(counts["successCountTotal"])
        + ", Total fail: "
        + str(counts["failCountTotal"])
        + ". Total missing refs: "
        + str(counts["missingRefCountTotal"])
        + ". Total duplicates: "
        + str(counts["duplicatesCountTotal"])
        + ". Total successful chunks: "
        + str(counts["successful_chunks"])
        + ". Total fail chunks: "
        + str(counts["fail_chunks"])
    )


def main():
    now = datetime.now()
    counts = {}
    counts["successCountTotal"] = 0
    counts["failCountTotal"] = 0
    counts["missingRefCountTotal"] = 0
    counts["duplicatesCountTotal"] = 0
    counts["successful_chunks"] = 0
    counts["fail_chunks"] = 0

    for modelName, action_info in model_import_actions.items():
        model_counts = {}
        model_counts["successCountTotal"] = 0
        model_counts["failCountTotal"] = 0
        model_counts["missingRefCountTotal"] = 0
        model_counts["duplicatesCountTotal"] = 0
        model_counts["successful_chunks"] = 0
        model_counts["fail_chunks"] = 0

        if action_info.get("skip"):
            continue
        modelNow = datetime.now()
        # TODO: consider running multiple times / retrying
        if modelName not in pk_maps:
            pk_maps[modelName] = {}
        if modelName not in next_id_maps:
            next_id_maps[modelName] = 1
        if action_info.get("empty_first"):
            print("Emptying table " + modelName)
            pk_maps[modelName] = {}
            next_id_maps[modelName] = 1
            with engine.connect() as connection:
                connection.execute(text("DELETE FROM " + modelName.upper()))
                # reset table id counter
                connection.execute(
                    text("ALTER TABLE " + modelName.upper() + " AUTO_INCREMENT = 1")
                )

        # get files from folder
        # sort alphabetically but with numbers in order
        sorted_files = sorted(
            os.listdir(rootDir + "/" + modelName),
            # key=lambda x: x.reverse(),
            #            reverse=True
        )
        #        sorted_files = sorted_files[:use_n_chromosomes]
        for file in sorted_files:
            if file.endswith(".tsv"):
                targetFile = rootDir + "/" + modelName + "/" + file
                file_info = inspectTSV(targetFile)
                print(
                    "\nimporting "
                    + modelName
                    + " ("
                    + targetFile.split("/")[-1]
                    + "). Expecting "
                    + str(file_info["total_rows"])
                    + " rows..."
                )
                # print(targetFile)

                results = import_chr_into_table(
                    targetFile,
                    file_info,
                    action_info,
                )
                successCount = results.get("success") or 0
                failCount = results.get("fail") or 0
                missingRefCount = results.get("missingRef") or 0
                duplicateCount = results.get("duplicate") or 0
                successful_chunks = results.get("successful_chunks") or 0
                fail_chunks = results.get("fail_chunks") or 0

                percent_success = "N/A"
                if successCount + failCount > 0:
                    percent_success = 100 * successCount / (successCount + failCount)
                print(
                    "Success: "
                    + str(successCount)
                    + ", Fail: "
                    + str(failCount)
                    + ". total: "
                    + str(file_info["total_rows"])
                    + ". missing refs: "
                    + str(missingRefCount)
                    + ". Success rate: "
                    + str(percent_success)
                    + "%. Duplicates: "
                    + str(duplicateCount)
                    + ". Successful chunks: "
                    + str(successful_chunks)
                    + ". Fail chunks: "
                    + str(fail_chunks)
                )
                persist_maps()
                if successCount == 0:
                    print("No rows were imported. Exiting.")
                    quit()

                model_counts["successCountTotal"] += successCount
                model_counts["failCountTotal"] += failCount
                model_counts["missingRefCountTotal"] += missingRefCount
                model_counts["duplicatesCountTotal"] += duplicateCount
                model_counts["successful_chunks"] += successful_chunks
                model_counts["fail_chunks"] += fail_chunks

                counts["successCountTotal"] += successCount
                counts["failCountTotal"] += failCount
                counts["missingRefCountTotal"] += missingRefCount
                counts["duplicatesCountTotal"] += duplicateCount
                counts["successful_chunks"] += successful_chunks
                counts["fail_chunks"] += fail_chunks

        print(
            "Finished importing "
            + modelName
            + ". Took this much time: "
            + str(datetime.now() - modelNow)
        )
        report_counts(model_counts)
    print("finished importing IBVL. Time Taken: " + str(datetime.now() - now))
    report_counts(counts)


main()
