# IBVL

## Repo organization

This repo is intended for Wasserman lab members working on data processing for IBVL.

The repo is organized in sub-folders depending on the different aspects of the data processing.

### Metadata tracking
Concerning metadata tracking (tracking ofinformation associated with each sample).

One of the considered tool to track metadata is OpenCGA, refer to [the openCGA folder](https://github.com/scorreard/IBVL/tree/main/opencga) for more information.

### Nextflow Scripts
Concerning the scripts used to generate the IBVL.

The [Nextflow wrapper](https://www.nextflow.io) is used to allow treacability and reproducibility, to review / comment the scripts, refer to [the script folder](https://github.com/scorreard/IBVL/tree/main/Nextflow_script)

## Import directory

How to run an import:
  1) copy the `import/.env-sample` file to `import/.env` and set values appropriately
  2) (optional) if you need to, run `python tables.py` to create the tables (database should be empty before this)
  3) `python orchestrate.py` will kick off the migration

The script creates a directory called "jobs", and a directory inside that called "1" the first time, "2" the second time, eg. 

Each of these job folders has working data for the migration and two output logs (one for errors, one for progress). The working data is just (for each model) a file with the latest primary key, and a reverse lookup map for entity id (eg gene or variant or transcript id) to primary key.

### Import environment vars
  - `ORACLE_TABLE_PATH` - (not necessarily just for oracle destinations) the full path to the directory containing pipeline output files
  - `COPY_MAPS_FROM_JOB` - The script maintains maps in order to resolve primary keys, they are persisted as json to the job directory, named using an incrementing number. If a job fails and you want to use the maps from a previous run, enter the run's job folder number as the value of this environment variable
  - `SCHEMA_NAME` - for an Oracle destination db, the schema name goes here.
  - `START_AT_MODEL` - to pick up after a previous migration run left off, you can enter the model name here, and the script will skip to that model (it runs in the order of keys as defined in the `model_import_actions` map)
  - (`START_AT_FILE`) - for convenience, you can also skip to a particular file in the first model dir imported, using natural sorting. Be very careful if using this in production as it will lead to false duplicates unless the primary key for new row insertions is corrected.
