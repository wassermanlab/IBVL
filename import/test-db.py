from sqlalchemy import (
    create_engine,
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
rootDir = os.environ.get("ORACLE_TABLE_PATH")
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
if isinstance(container, str) and len(container) > 0:
    print("hooking into connect event to set container = "+container)
    def set_container(dbapi_connection, connection_record):
        cursor = dbapi_connection.cursor()
        cursor.execute("ALTER SESSION SET CONTAINER = "+container)
        cursor.close()
        print("container was set.")
    engine = create_engine(dbConnectionString, echo=True, pool_pre_ping=True, pool_recycle=3600 )
    event.listen(engine, 'connect', set_container)
else:
    engine = create_engine(dbConnectionString)

def check_for_bail(question):
    doContinue = input(question)
    if (doContinue != "y"):
        print("exiting...")
        engine.dispose()
        quit()

with engine.connect() as connection:

    print("looks like connecting worked")

    check_for_bail("continue to read test? (y/n): ")

    input_name = input("enter table name (default: GENES): ")
    table_name = input_name or "GENES"


    try:
        df = pd.read_sql_query("SELECT * FROM "+table_name, engine)
        print(df)
    except Exception:
        pass
    check_for_bail("continue to table test? (y/n): ")

    metadata = MetaData()

    print("initializing Table for "+table_name+"...")

    if isinstance(schema, str) and len(schema) > 0:
        print("setting schema = "+schema+"...")
        test_table = Table(
            table_name,
            metadata,
            Column("ID", Integer, primary_key=True),
            Column("SHORT_NAME", String(30), nullable=False),
            schema=schema,
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


    try:
        result = connection.execute(test_table.insert(), [{"SHORT_NAME": "TEST"}])
        print("inserted ")
    except Exception:
        traceback.print_exc()

engine.dispose()