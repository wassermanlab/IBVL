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

load_dotenv()

# get command line arguments
rootDir = os.environ.get("ORACLE_TABLE_PATH")
chunk_size = int(os.environ.get("CHUNK_SIZE"))
verbose = os.environ.get("VERBOSE") == "true"

if rootDir == None:
    print("No root directory specified")
    exit()

dbConnectionString = os.environ.get("DB")
container = os.environ.get("DB_CONTAINER")






def set_container(dbapi_connection, connection_record):
    cursor = dbapi_connection.cursor()
    cursor.execute("ALTER SESSION SET CONTAINER = "+container)
    cursor.close()

print("connecting...")
if (container != None):
    engine = create_engine(dbConnectionString, echo=True, pool_pre_ping=True, pool_recycle=3600 )
    event.listen(engine, 'connect', set_container)
else:
    engine = create_engine(dbConnectionString, pool_pre_ping=True, pool_recycle=3600, connect_args={'autocommit': True})

mydb = engine.connect()

print("looks like connecting worked")

doContinue = input("continue? (y/n): ")

if doContinue != "y":
    print("exiting...")
    engine.dispose()
    quit()

print("testing read...")

df = pd.read_sql_query("SELECT * FROM GENES", engine)
print(df)



doContinue = input("continue? (y/n): ")

if doContinue != "y":
    print("exiting...")
    engine.dispose()
    quit()

metadata = MetaData()

test_table = Table(
    "TEST",
    metadata,
    Column("ID", Integer, primary_key=True),
    Column("SHORT_NAME", String(30), nullable=False)
)

print("creating test table...")
metadata.create_all(engine)

print("test complete")

engine.dispose()