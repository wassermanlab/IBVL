from sqlalchemy import (
    create_engine,
    Sequence,
    text,
    event,
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
import sys
import os
from dotenv import load_dotenv
from datetime import datetime
import traceback
import oracledb

load_dotenv()

# get command line arguments
rootDir = os.environ.get("PIPELINE_OUTPUT_PATH")
chunk_size = int(os.environ.get("CHUNK_SIZE"))
instant_client = os.environ.get("ORACLE_INSTANT_CLIENT_PATH")
oracle_wallet = os.environ.get("ORACLE_WALLET_PATH")
schema = os.environ.get("SCHEMA_NAME")

if isinstance(instant_client, str) and isinstance(oracle_wallet, str) and len(instant_client) > 0 and len(oracle_wallet) > 0:
    print("initializing oracle client...")
    oracledb.init_oracle_client(lib_dir=instant_client, config_dir=oracle_wallet)


if rootDir == None:
    print("No root directory specified")
    exit()

dbConnectionString = os.environ.get("DB")
container = os.environ.get("DB_CONTAINER")


print("connecting...")

engine = create_engine(dbConnectionString, echo=True, pool_pre_ping=True, pool_recycle=3600)

def check_for_bail(question):
    doContinue = input(question)
    if (doContinue != "y"):
        print("exiting...")
        engine.dispose()
        quit()

with engine.connect() as connection:

    print("looks like connecting worked. will try to reflect ")

    metadata = MetaData()

    metadata.reflect(bind=engine, schema=schema)

    print("reflected tables:" + str(metadata.sorted_tables))

    check_for_bail("continue to table test? (y/n): ")

    table_name = 'genes'

    print("initializing Table for "+table_name+"...")

    if isinstance(schema, str) and len(schema) > 0:
#        print("setting schema = "+schema+"...")
        test_table = Table(
            table_name,
            metadata,
#            autoload_with=engine,
            schema=schema
        )
    else:
        print("no schema provided...")
        test_table = Table(
            table_name,
            metadata,
            Column("ID", Integer, primary_key=True),
            Column("SHORT_NAME", String(30), nullable=False),
        )

    try:
        result = connection.execute(test_table.select())
        print("table exists")
        print(result.fetchone())
    except Exception:
        traceback.print_exc()


    check_for_bail("continue to write table test? (y/n): ")

    n = input("name?")


    try:
        result = connection.execute(test_table.insert(), {"short_name": n})
        #commit
        connection.commit()
        print("inserted ")
    except Exception:
        traceback.print_exc()

engine.dispose()