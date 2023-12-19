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

print("connecting...")
if (verbose):
    engine = create_engine(dbConnectionString, echo=True, pool_pre_ping=True, pool_recycle=3600, connect_args={'autocommit': True})
else:
    engine = create_engine(dbConnectionString, pool_pre_ping=True, pool_recycle=3600, connect_args={'autocommit': True})

mydb = engine.connect()

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