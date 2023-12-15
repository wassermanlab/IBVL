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

print(rootDir)

if (verbose):
    engine = create_engine(dbConnectionString, echo=True, pool_pre_ping=True, pool_recycle=3600, connect_args={'autocommit': True})
else:
    engine = create_engine(dbConnectionString, pool_pre_ping=True, pool_recycle=3600, connect_args={'autocommit': True})

mydb = engine.connect()

metadata = MetaData()

dirs = [
    "genes",
    "transcripts"
    ]
#chromosomes = ['chr1','chr2','chr3']

chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'MT', 'str', 'sv_MEI_chr1', 'sv_MEI_chr2', 'sv_MEI_chr3', 'sv_MEI_chr4', 'sv_MEI_chr5', 'sv_MEI_chr6', 'sv_MEI_chr7', 'sv_MEI_chr8', 'sv_MEI_chr9', 'sv_MEI_chr10', 'sv_MEI_chr11', 'sv_MEI_chr12', 'sv_MEI_chr13', 'sv_MEI_chr14', 'sv_MEI_chr15', 'sv_MEI_chr16', 'sv_MEI_chr17', 'sv_MEI_chr18', 'sv_MEI_chr19', 'sv_MEI_chr20', 'sv_MEI_chr21', 'sv_MEI_chr22', 'sv_MEI_chrX', 'sv_MEI_chrY']
# Define the "GENES" table
genes_table = Table(
    "GENES",
    metadata,
    Column("ID", Integer, primary_key=True),
    Column("SHORT_NAME", String(30), nullable=False)
#    UniqueConstraint("SHORT_NAME", name="UNIQUE"),
)


transcripts_table = Table(
    "TRANSCRIPTS",
    metadata,
    Column("TRANSCRIPT_ID", String(255)),
    Column("ID", Integer, primary_key=True),
    Column("GENE", Integer),
    Column("TRANSCRIPT_TYPE", String(1)),
    Column("TSL", String(255)),
#    UniqueConstraint("TRANSCRIPT_ID", name="TRANSCRIPTS_UNIQUE"),
 #   ForeignKeyConstraint("GENE", name="TRANSCRIPTS_GENES_FK", ondelete="CASCADE", table="GENES"),
)
def transform_header(header):
    # Example: Convert all column names to uppercase
    return [col.upper() for col in header]

def readTSV(file):
    df = pd.read_csv(file,sep="\t")#,header=0,names=transform_header)
    df.columns = [col.upper() for col in df.columns]
    return df

def inspectTSV(file):
    total_rows = 0
    num_columns = pd.read_csv(file, sep='\t', nrows=4, header=None).shape[1]

    for chunk in pd.read_csv(file, sep='\t', chunksize=chunk_size):
        total_rows += len(chunk)

    return {
        "total_rows": total_rows,
        "num_columns": num_columns
    }

pk_maps = {};
next_id_maps = {};

for model in dirs:
    pk_maps[model] = {}
    next_id_maps[model] = 1

def import_genes(file, file_info):
    df = readTSV(file)

    with engine.connect() as connection:
#        connection.execute(delete(genes_table))
#        print(file_info)
        successCount = 0
        failCount = 0
        for index,row in df.iterrows():
            data = row.to_dict()
            pk = next_id_maps["genes"] + index
            data["ID"] = pk
#            print(data)

            cmd = genes_table.insert().values(data)
            try:
                connection.execute(cmd)
                successCount +=1
                pk_maps["genes"][data["SHORT_NAME"]] = pk
            except DataError as e:
                print(e)
                failCount +=1
            except IntegrityError as e:
                print(e)
                failCount +=1
            except Exception as e:
                print(e)
                failCount +=1
        next_id_maps["genes"] += file_info["total_rows"]
    return {
        "success": successCount,
        "fail": failCount
    }

def resolve_PK(referencedModel, name):
    try:
        return pk_maps[referencedModel][name]
    except KeyError:
#        print("Could not find PK in "+referencedModel+": "+str(name))
        # query db?
        return None

def import_transcripts(file, file_info):

    df = readTSV(file)
    df.fillna('NA', inplace=True)
    with engine.connect() as connection:
        successCount = 0
        failCount = 0

        for index,row in df.iterrows():
            data = row.to_dict()
            pk = next_id_maps["transcripts"] + index
            data["ID"] = pk
            
            if(data["GENE"] == "NA"):
                data["GENE"] = None
            else:
                resolved_pk = resolve_PK("genes", data["GENE"])
                if resolved_pk != None:
                    data["GENE"] = resolved_pk
                else:
                    failCount +=1
                    continue
            cmd = transcripts_table.insert().values(data)
            try:
                connection.execute(cmd)
                successCount +=1
                pk_maps["transcripts"][data["TRANSCRIPT_ID"]] = pk
            except DataError as e:
                print(e)
                failCount +=1
            except IntegrityError as e:
                print(e)
                failCount +=1
            except ProgrammingError as e:
                print(e)
                exit()
            except Exception as e:
                print(e)
                failCount +=1
        next_id_maps["transcripts"] += file_info["total_rows"]
    return {
        "success": successCount,
        "fail": failCount
    }
# map of import functions
importFunctions = {
    "genes": import_genes,
    "transcripts": import_transcripts
    }

def emptyTables():
    with engine.connect() as connection:
            
        for modelName in dirs:
            connection.execute(text("DELETE FROM "+modelName.upper()))
    #        mydb.execute("DELETE FROM "+modelName.upper())


emptyTables()
now = datetime.now()
for chromosome in chromosomes:
    successCount = 0
    failCount = 0
    for modelName in dirs:
        # get all files in folder
        targetFile = rootDir + "/" + modelName + "/" + modelName + "_" + chromosome + ".tsv"
        if importFunctions[modelName] == None:
            print("No import function for " + modelName)
            exit()
        else:
    #        mydb.drop(modelName, if_exists="cascade")
    #       mydb.execute(createStatements[modelName])
            
            if os.path.isfile(targetFile):
                print("importing "+chromosome+" to Table: "+modelName+"...")
                #print(targetFile)

                file_info = inspectTSV(targetFile)
                results = importFunctions[modelName](targetFile, file_info)
                print("finished importing "+chromosome+" into "+modelName+". Success: "+str(results["success"])+", Fail: "+str(results["fail"])+". total: "+str(file_info["total_rows"]))
            else:
                print("File not found: "+targetFile)
    print("Chromosome imported: "+chromosome)

print("Finished importing all chromosomes. Took this much time: "+str(datetime.now() - now))
# Read TSV
# df = pd.read_csv('genes/genes_chr1.tsv', sep='\t')

# Create cursor
# cursor = mydb.cursor()

# Create table
# mydb.execute("CREATE TABLE genes (gene_id INT AUTO_INCREMENT PRIMARY KEY, gene_name VARCHAR(255), gene_start INT, gene_end INT, gene_strand VARCHAR(255), gene_biotype VARCHAR(255), gene_chromosome VARCHAR(255))")


# Insert data

# df.to_sql('genes', con=mydb, if_exists='append', index=False)
