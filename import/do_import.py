import json
from natsort import natsorted
from sqlalchemy import (
    create_engine,
    text,
    MetaData,
    Table,
    Column,
    Integer,
    String,
    select,
    func,
    delete,
    Float,
    UniqueConstraint,
    ForeignKeyConstraint,
)
from sqlalchemy.exc import DataError, IntegrityError, ProgrammingError
import pandas as pd
import numpy as np
import signal
import sys
import os
from dotenv import load_dotenv
from datetime import datetime
from sqlalchemy import text

from import_utils import *
from sqlalchemy.orm import sessionmaker
import traceback


load_dotenv()

# get command line arguments
rootDir = os.environ.get("ORACLE_TABLE_PATH")
chunk_size = int(os.environ.get("CHUNK_SIZE"))
#verbose = os.environ.get("VERBOSE") == "true"
dbConnectionString = os.environ.get("DB")
copy_maps_from_job = os.environ.get("COPY_MAPS_FROM_JOB")
isDevelopment = os.environ.get("ENVIRONMENT") != "production"
schema = os.environ.get("SCHEMA_NAME")

if rootDir == None:
    print("No root directory specified")
    exit()

data_issue_logger = None
output_logger = None

print(rootDir)

engine = None

metadata = MetaData()

job_dir = ""
maps_load_dir = ""
# map of import functions
model_import_actions = {
    "genes": {
        "name": "genes",
        "pk_lookup_col": "SHORT_NAME",
        "fk_map": {},
        "empty_first": True,
        "filters": {
            "SHORT_NAME": lambda x: x.upper()
        },
        #"skip":True
    },
    "transcripts": {
        "name": "transcripts",
        "pk_lookup_col": "TRANSCRIPT_ID",
        "fk_map": {"GENE": "genes"},
        "empty_first": True,
        "filters": {
            "TRANSCRIPT_TYPE": lambda x: x.replace("RefSeq", "R"),
        },
        #"skip":True
    },
    "variants": {
        "name": "variants",
        "pk_lookup_col": "VARIANT_ID",
        "fk_map": {},
        "empty_first": True,
        #"skip":True
    },
    "variants_transcripts": {
        "name": "variants_transcripts",
        "pk_lookup_col": ["TRANSCRIPT", "VARIANT"],
        "fk_map": {"TRANSCRIPT": "transcripts", "VARIANT": "variants"},
        "empty_first": True,
        #"skip":True
    },
    "variants_annotations": {
        "name": "variants_annotations",
        "pk_lookup_col": None,
        "fk_map": {"DO_COMPOUND_FK": "for variants_transcripts"},
        "empty_first": True,
        "filters":{
            "HGVSP": lambda x: x.replace("%3D","=")
        },
        #"skip":True
    },
    "variants_consequences": {
        "name": "variants_consequences",
        "pk_lookup_col": None,
        "fk_map": {"DO_COMPOUND_FK": "for variants_transcripts"},
        "empty_first": True,
        #"skip":True,
    },
    "sv_consequences": {
        "name": "sv_consequences",
        "pk_lookup_col": None,
        "fk_map": {"GENE": "genes", "VARIANT": "variants"},
        "empty_first": True,
        #"skip":True
    },
    "snvs": {
        "name": "snvs",
        "pk_lookup_col": None,
        "fk_map": {"VARIANT": "variants"},
        "empty_first": True,
        "filters":{
            "DBSNP_ID": lambda x: x.split('&')[0] if x is not None else None
        },
        #"skip":True
    },
    "svs": {
        "name": "svs",
        "pk_lookup_col": None,
        "fk_map": {"VARIANT": "variants"},
        "empty_first": True,
        #"skip":True
    },
    "svs_ctx": {
        "name": "svs_ctx",
        "pk_lookup_col": None,
        "fk_map": {"VARIANT": "variants"},
        "empty_first": True,
        #"skip":True
    },
    "str": {
        "name": "str",
        "pk_lookup_col": None,
        "fk_map": {"VARIANT": "variants"},
        "empty_first": True,
        #"skip":True
    },
    "mts": {
        "name": "mts",
        "pk_lookup_col": None,
        "fk_map": {"VARIANT": "variants"},
        "empty_first": True,
        #"skip":True
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


pk_maps = {}
next_id_maps = {}
tables = {}

def load_maps(models=[]):
    try:
        # Load from files
        for modelName in models:
            with open(maps_load_dir + "/" + modelName + "_pk_map.json", "r") as f:
                pk_maps[modelName] = json.load(f)
            with open(maps_load_dir + "/" + modelName + "_next_id.json", "r") as f:
                next_id_maps[modelName] = json.load(f)
            log_output("loaded map for " + modelName +". number of records: " + str(len(pk_maps[modelName])))

    except FileNotFoundError:
        pk_maps[modelName] = {}
        pass

def append_to_map(modelName, key, value):
    if modelName not in pk_maps:
        log_output("")
        load_maps(models=[modelName])
    if key not in pk_maps[modelName]:
        pk_maps[modelName][key] = value

def persist_and_unload_maps():
    try:
        for modelName, pk_map in pk_maps.items():
            log_output("saving pk map for " + modelName)
            with open(os.path.join(job_dir,modelName+"_pk_map.json"), "w") as f:
                json.dump(pk_map, f)
        for modelName, next_id in next_id_maps.items():
            with open(os.path.join(job_dir, modelName+"_next_id.json"), "w") as f:
                json.dump(next_id, f)
    except Exception as e:
        log_data_issue("Error saving maps")
        quit()
    pk_maps.clear()
    log_output("cleared the pk maps")

def resolve_PK(referencedModel, name):
    try:
        result = pk_maps[referencedModel][name.upper()]
        return result
    except KeyError:
        return None

def get_table(model):
    global tables
    if model in tables:
        return tables[model]
    else:
        if isinstance(schema, str) and len(schema) > 0:
            table = Table(model, metadata, schema=schema)
        else:
            table = Table(model.upper(), metadata, autoload_with=engine)
        tables[model] = table
        return table
    
def inject(model, data, map_key):
# need to dynamically inject the single obj that was missing from original data
    pk = None
    table = get_table(model)
    with engine.connect() as connection:
        try:
            result = connection.execute(table.insert(), data)
            connection.commit()
            pk = result.inserted_primary_key[0]
            append_to_map(model, map_key, pk)
            next_id_maps[model]  = pk + 1
            log_data_issue(f"dynamically added to {model}: {data}")
        except IntegrityError as e:
            log_data_issue("a dynamically injected obj had an integrity error.")
            log_data_issue(e)
#                                quit() # LATER: comment this out
        except Exception as e:
            log_data_issue("a dynamically injected obj had an error.")
            log_data_issue(e)
#                                quit() # LATER: comment this out?
    return pk

def import_file(file, file_info, action_info):
    name = action_info.get("name")
    fk_map = action_info.get("fk_map")
    pk_lookup_col = action_info.get("pk_lookup_col")
    filters = action_info.get("filters") or {}

    missingRefCount = 0
    table = get_table(name)
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
                map_key = data[fk_col].upper()
                if map_key == "NA":
                    data[fk_col] = None
                else:
                    resolved_pk = resolve_PK(fk_model, map_key)
                    ## resolved PK was not found from maps, so.. if it's a gene, we could dynamically inject
                    if (resolved_pk == None and fk_col == "GENE" and name == "transcripts"):
                        resolved_pk = inject("genes",{"SHORT_NAME":map_key}, map_key)
                    elif (resolved_pk == None and fk_col == "VARIANT" and name in ["sv_consequences", "svs", "snvs", "mts"]):
                        
                        if (name == "sv_consequences" or name == "svs"):
                            var_type = "SV"
                        elif (name == "snvs"):
                            var_type = "SNV"
                        elif (name == "mts"):
                            var_type = "MT"
                        resolved_pk = inject("variants",{"VARIANT_ID":map_key, "VAR_TYPE": var_type}, map_key)
            if map_key is not None and False:
                log_output(
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
                log_data_issue(
                    "Missing "
                    + fk_col
                    + " "
                    + map_key
                    + " referenced from "
                    + name
                )
                if (debug_row is not None):
                    log_data_issue(debug_row)
                else:
                    log_data_issue(data)
                missingRefCount += 1
                skip = True
        if skip:
            continue
        data_list.append(data)

    # dispose of df to save ram
    del df
    with engine.connect() as connection:
        successCount = 0
        failCount = 0
        duplicateCount = 0
        successful_chunks = 0
        fail_chunks = 0

        for chunk in chunks(data_list, chunk_size):
            try:
                connection.execute(table.insert(), chunk)
                #commit
                connection.commit()
                # chunk worked
                successful_chunks += 1
                successCount += len(chunk)
            except Exception as e:
                #                print(e)
                fail_chunks += 1
                for row in chunk:
                    try:
                        connection.execute(table.insert(), row)
                        connection.commit()
                        successCount += 1

                    except DataError as e:
                        log_data_issue(e)
                        failCount += 1
#                        quit()
                    except IntegrityError as e:
                        msg = str(e)
                        if "Duplicate" in msg:
                            duplicateCount += 1
                            successCount += 1
                        else:
                            failCount += 1
                            log_data_issue(e)
#                            quit()
                    except Exception as e:
                        log_data_issue(e)
                        failCount += 1

            if pk_lookup_col is not None:
                pk_map = {}
                for data in chunk:
                    # record the PKS for each row that was added
                    if isinstance(pk_lookup_col, list):
#                        log_output(pk_lookup_col)
#                        log_output(data)
                        map_key = "-".join([str(data[col]) for col in pk_lookup_col])
                    elif isinstance(pk_lookup_col, str):
                        map_key = data[pk_lookup_col].upper()
                    
                    if map_key not in pk_map:
                        pk_map[map_key] = data["ID"]
                    if False:
                        log_output("added " + name + "." + map_key + " to pk map")
                for key in pk_map:
                    append_to_map(name, key.upper(), pk_map[key])

    next_id_maps[name] += file_info["total_rows"]

    return {
        "success": successCount,
        "fail": failCount,
        "missingRef": missingRefCount,
        "duplicate": duplicateCount,
        "successful_chunks": successful_chunks,
        "fail_chunks": fail_chunks,
    }

def cleanup(sig, frame):
    global engine, pk_maps, next_id_maps, tables, metadata, data_issue_logger, output_logger
    print('cleaning up ...')
    persist_and_unload_maps()
    engine.dispose()
    #garbage collect
    del pk_maps
    del next_id_maps
    del tables
    del metadata
    del data_issue_logger
    del output_logger
    print('done')
    sys.exit(0)

signal.signal(signal.SIGINT, cleanup)

def start(db_engine):


    # Assuming 'engine' is your Engine object
    global job_dir, maps_load_dir, engine, schema
    engine = db_engine

    if isinstance(schema,str) and len(schema) > 0:
        metadata.reflect(bind=engine, schema=schema)
    else:
        metadata.reflect(bind=engine)
    Session = sessionmaker(bind=engine)
    jobs_dir = os.path.abspath(os.path.join("", "jobs"))
    os.makedirs(jobs_dir, exist_ok=True)
    os.makedirs(os.path.join(jobs_dir, "1"), exist_ok=True)
    
    without_hidden = [f for f in os.listdir(jobs_dir) if not f.startswith('.')]
    last_job = int(natsorted(without_hidden)[-1])
    if (os.listdir(os.path.join(jobs_dir, str(last_job))) == []):
        job_dir = os.path.join(jobs_dir, str(last_job))
    else:
        job_dir = os.path.join(jobs_dir, str(last_job + 1))
    os.makedirs(job_dir, exist_ok=True)
    os.chmod(job_dir, 0o777)  # Set read and write permissions for the directory
    print("using job dir " + job_dir)
    setup_loggers(job_dir)

    if copy_maps_from_job is not None and copy_maps_from_job != "":
        maps_load_dir = os.path.join(jobs_dir, copy_maps_from_job)
    else:
        maps_load_dir = job_dir

    now = datetime.now()
    counts = {}
    counts["success"] = 0
    counts["fail"] = 0
    counts["missingRef"] = 0
    counts["duplicate"] = 0
    counts["successful_chunks"] = 0
    counts["fail_chunks"] = 0

    for modelName, action_info in model_import_actions.items():
        model_counts = {}
        model_counts["success"] = 0
        model_counts["fail"] = 0
        model_counts["missingRef"] = 0
        model_counts["duplicate"] = 0
        model_counts["successful_chunks"] = 0
        model_counts["fail_chunks"] = 0

        if action_info.get("skip"):
            continue

        referenced_models = action_info.get("fk_map").values()
        if "DO_COMPOUND_FK" in action_info.get("fk_map"):
            referenced_models = ["variants_transcripts", "variants", "transcripts"]
        load_maps(models=referenced_models)
        modelNow = datetime.now()
        
        # if modelName not in pk_maps:
        #     pk_maps[modelName] = {}
        if modelName not in next_id_maps:
             next_id_maps[modelName] = 1
        if action_info.get("empty_first"):
            log_output("Emptying table " + modelName)
            table = get_table(modelName)


            # Assuming 'table' is your Table object
            # Replace 'ID' with your actual column name
            with Session() as session:
                max_id = session.query(func.max(table.columns['ID'])).scalar()

            print("max id found to be "+str(max_id))
            with engine.connect() as connection:

                # Split the deletion into smaller chunks
                chunk_size = 1000
                offset = 0
                while True:
                    try:
                        delete_stmt = table.delete().where(table.c.ID <= max_id - offset).where(table.c.ID > max_id - chunk_size - offset)
                        connection.execute(delete_stmt)
                        offset += chunk_size
                        connection.commit()
#                        print("did delete a chunk")
                    except ProgrammingError as e:
                        print(e)
                        break
                    except Exception as e:
                        print("error emptying table " + modelName)
                        print(e)
                        break
        sorted_files = natsorted(
            [f for f in os.listdir(rootDir + "/" + modelName) if not f.startswith('.')],
        )
        for file in sorted_files:
            if file.endswith(".tsv"):
                targetFile = rootDir + "/" + modelName + "/" + file
                file_info = inspectTSV(targetFile)
                log_output(
                    "\nimporting "
                    + modelName
                    + " ("
                    + targetFile.split("/")[-1]
                    + "). Expecting "
                    + str(file_info["total_rows"])
                    + " rows..."
                )
                # log_output(targetFile)
                if (file_info["total_rows"] == 0):
                    log_output("Skipping empty file")
                    continue
                results = import_file(
                    targetFile,
                    file_info,
                    action_info,
                )
                if results["success"] == 0:
                    log_output("No rows were imported.")
                    
                for key in ["success", "fail", "missingRef", "duplicate", "successful_chunks", "fail_chunks"]:
                    model_counts[key] += results[key]
                    counts[key] += results[key]

                report_counts(results)

        log_output(
            "\nFinished importing "
            + modelName
            + ". Took this much time: "
            + str(datetime.now() - modelNow)
        )
        report_counts(model_counts)
        this_model_index = list(model_import_actions.keys()).index(modelName)
        if this_model_index + 1 < len(model_import_actions.keys()):
            leftover_models = list(model_import_actions.keys())[this_model_index+1:]
            log_output("\nmodels left still: " + str(leftover_models) + "\n")
        
        persist_and_unload_maps()
    log_output("finished importing IBVL. Time Taken: " + str(datetime.now() - now))
    report_counts(counts)
    cleanup(None, None)

