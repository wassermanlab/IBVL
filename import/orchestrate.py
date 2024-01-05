from do_import import start
from dotenv import load_dotenv
import os
from sqlalchemy import create_engine
from sqlalchemy_utils import database_exists, create_database
from sqlalchemy import (
    create_engine,
    event
)

load_dotenv()

dbConnectionString = os.environ.get("DB")
container = os.environ.get("DB_CONTAINER")
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

def set_container(dbapi_connection, connection_record):
    cursor = dbapi_connection.cursor()
    cursor.execute("ALTER SESSION SET CONTAINER = "+container)
    cursor.close()

print("connecting...")
if isinstance(container, str) and len(container) > 0:
    engine = create_engine(dbConnectionString, echo=True, pool_pre_ping=True, pool_recycle=3600 )
    event.listen(engine, 'connect', set_container)
else:
    engine = create_engine(
        dbConnectionString,
        pool_pre_ping=True,
        pool_recycle=3600,
        connect_args={"autocommit": True},
    )
start(engine)   

exit()

