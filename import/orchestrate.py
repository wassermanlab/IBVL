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
verbose = os.environ.get("VERBOSE") == "true" or os.environ.get("VERBOSE") == "True"
container = os.environ.get("DB_CONTAINER")
isDevelopment = os.environ.get("ENVIRONMENT") != "production"

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


print("connecting...")
if isinstance(container, str) and len(container) > 0:
    # ORACLE (prod)

    print("hooking into connect event to set container = "+container)
    def set_container(dbapi_connection, connection_record):
        cursor = dbapi_connection.cursor()
        cursor.execute("ALTER SESSION SET CONTAINER = "+container)
        cursor.close()
        print("container was set.")
    engine = create_engine(dbConnectionString, echo=verbose, pool_pre_ping=True, pool_recycle=3600)
    event.listen(engine, 'connect', set_container)
else:
    # MYSQL (dev)
    engine = create_engine(
        dbConnectionString,
        echo=verbose,
        pool_pre_ping=True,
        pool_recycle=3600
    )
start(engine)   

exit()

