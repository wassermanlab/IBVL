from do_import import start
from dotenv import load_dotenv
import os
from sqlalchemy import create_engine
from sqlalchemy_utils import database_exists, create_database
from sqlalchemy import (
    create_engine
)

load_dotenv()

dbConnectionString = os.environ.get("DB")
isDevelopment = os.environ.get("ENVIRONMENT") == "development"

db_name = dbConnectionString.split("/")[-1]

if isDevelopment:
    # create the database if it doesn't exist
    engine = create_engine(dbConnectionString, echo=True)
    if database_exists(engine.url):
        #assume already has structure
        pass
    else:
        create_database(dbConnectionString)
        import tables
    engine.dispose()

engine = create_engine(
    dbConnectionString,
    pool_pre_ping=True,
    pool_recycle=3600,
    connect_args={"autocommit": True},
)
start(engine)   

exit()

