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
dbConnectionString = os.environ.get("DB")

if rootDir == None:
    print("No root directory specified")
    exit()


print(rootDir)

if False and verbose:
    engine = create_engine(
        dbConnectionString,
        echo=True,
        pool_pre_ping=True,
        pool_recycle=3600,
        connect_args={"autocommit": True},
    )
else:
    engine = create_engine(
        dbConnectionString,
        pool_pre_ping=True,
        pool_recycle=3600,
        connect_args={"autocommit": True},
    )

mydb = engine.connect()

metadata = MetaData()

dirs = ["genes", "transcripts", "severities", "variants"]
chromosomes = ["chrY", "chrX", "chr10"]
# chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'MT', 'str', 'sv_MEI_chr1', 'sv_MEI_chr2', 'sv_MEI_chr3', 'sv_MEI_chr4', 'sv_MEI_chr5', 'sv_MEI_chr6', 'sv_MEI_chr7', 'sv_MEI_chr8', 'sv_MEI_chr9', 'sv_MEI_chr10', 'sv_MEI_chr11', 'sv_MEI_chr12', 'sv_MEI_chr13', 'sv_MEI_chr14', 'sv_MEI_chr15', 'sv_MEI_chr16', 'sv_MEI_chr17', 'sv_MEI_chr18', 'sv_MEI_chr19', 'sv_MEI_chr20', 'sv_MEI_chr21', 'sv_MEI_chr22', 'sv_MEI_chrX', 'sv_MEI_chrY']
# Define the "GENES" table

# max_rows=5


def transform_header(header):
    # Example: Convert all column names to uppercase
    return [col.upper() for col in header]


def readTSV(file):
    df = pd.read_csv(file, sep="\t")  # ,header=0,names=transform_header)
    df.columns = [col.upper() for col in df.columns]
    return df


def inspectTSV(file):
    total_rows = 0
    small_read = pd.read_csv(file, sep="\t", nrows=4, header=None)
    num_columns = small_read.shape[1]
    columns = [col.upper() for col in small_read.values.tolist()[0]]

    for chunk in pd.read_csv(file, sep="\t", chunksize=chunk_size):
        total_rows += len(chunk)

    return {"total_rows": total_rows, "num_columns": num_columns, "columns": columns}


pk_maps = {}
next_id_maps = {}
tables = {}

for model in dirs:
    pk_maps[model] = {}
    next_id_maps[model] = 1


def import_chr_into_table(file, file_info, name, pk_lookup_col=None, fk_map={}):
    df = readTSV(file)
    df.fillna("NA", inplace=True)

    print(file_info)

    table = tables.get(name)
    if table == None:
        table = Table(
            name.upper(),
            metadata,
            Column("ID"),
            *(Column(col) for col in file_info["columns"])
        )
        tables[name] = table

    with engine.connect() as connection:
        #        connection.execute(delete(genes_table))
        #        print(file_info)
        successCount = 0
        failCount = 0
        missingRefCount = 0
        duplicateCount = 0
        for index, row in df.iterrows():
            data = row.to_dict()
            pk = next_id_maps[name] + index
            data["ID"] = pk
            if verbose:
                print(data)

            skip = False
            for fk_col, fk_model in fk_map.items():
                # print("found fk map: "+fk_col+" -> "+fk_model+" in "+name)
                if data[fk_col] == "NA":
                    data[fk_col] = None
                else:
                    resolved_pk = resolve_PK(fk_model, data[fk_col])
                    if verbose:
                        print(
                            "resolved "
                            + fk_model
                            + "."
                            + data[fk_col]
                            + " to "
                            + str(resolved_pk)
                        )
                    if resolved_pk != None:
                        data[fk_col] = resolved_pk
                    else:
                        missingRefCount += 1
                        skip = True
            if skip:
                continue
            cmd = table.insert().values(data)
            try:
                connection.execute(cmd)
                successCount += 1
                if pk_lookup_col != None:
                    pk_maps[name][data[pk_lookup_col]] = pk
                    if verbose:
                        print(
                            "added " + name + "." + data[pk_lookup_col] + " to pk map"
                        )
            except DataError as e:
                print(e)
                failCount += 1
                quit()
            except IntegrityError as e:
                msg = str(e)
                if "Duplicate" in msg:
                    duplicateCount += 1
                    successCount += 1
                else:
                    failCount += 1
                    print(e)
                    quit()  # temp
            except Exception as e:
                print(e)
                failCount += 1
        next_id_maps[name] += file_info["total_rows"]
    return {
        "success": successCount,
        "fail": failCount,
        "missingRef": missingRefCount,
        "duplicate": duplicateCount,
    }


def resolve_PK(referencedModel, name):
    try:
        return pk_maps[referencedModel][name]
    except KeyError:
        #        print("Could not find PK in "+referencedModel+": "+str(name))
        # query db?
        return None


# map of import functions
model_operation_definitions = {
    "genes": {
        "f": lambda file, file_info: import_chr_into_table(
            file, file_info, "genes", "SHORT_NAME"
        ),
        "file_suffixes": chromosomes,
    },
    "transcripts": {
        "f": lambda file, file_info: import_chr_into_table(
            file, file_info, "transcripts", "TRANSCRIPT_ID", {"GENE": "genes"}
        ),
        "file_suffixes": chromosomes,
    },
    "variants": {
        "f": lambda file, file_info: import_chr_into_table(
            file, file_info, "variants", "VARIANT_ID"
        ),
        "file_suffixes": ["MT", "str"],
    },
    "mt_gmonad_frequencies": {
        "f": lambda file, file_info: import_chr_into_table(
            file,
            file_info,
            "mt_gmonad_frequencies",
            "VARIANT_ID",
            {"VARIANT": "variants"},
        ),
        "files": ["mt_gnomad_frequencies"],
    },
    "mt_ibvl_frequencies": {
        "f": lambda file, file_info: import_chr_into_table(
            file,
            file_info,
            "mt_ibvl_frequencies",
            "VARIANT_ID",
            {"VARIANT": "variants"},
        ),
        "files": ["mt_ibvl_frequencies"],
    },
    "mts": {
        "f": lambda file, file_info: import_chr_into_table(
            file, file_info, "mts", "VARIANT_ID", {"VARIANT": "variants"}
        ),
        "files": ["mts"],
    },
    "str": {
        "f": lambda file, file_info: import_chr_into_table(
            file, file_info, "str", "VARIANT_ID", {"VARIANT": "variants"}
        ),
        "files": ["str"],
    },
    "variants_transcripts": {
        "f": lambda file, file_info: import_chr_into_table(
            file,
            file_info,
            "variants_transcripts",
            "VARIANT_ID",
            {"TRANSCRIPT": "transcripts"},
        ),
        "file_suffixes": chromosomes,
    },
    "variants_annotations": {
        "f": lambda file, file_info: import_chr_into_table(
            file,
            file_info,
            "variants_annotations",
            "VARIANT_ID",
            {"TRANSCRIPT": "transcripts"},
        ),
        "file_suffixes": chromosomes,
    },
    "variants_consequences": {
        "f": lambda file, file_info: import_chr_into_table(
            file,
            file_info,
            "variants_consequences",
            "VARIANT_ID",
            {"TRANSCRIPT": "transcripts"},
        ),
        "file_suffixes": chromosomes,
    },
    "sv_consequences": {
        "f": lambda file, file_info: import_chr_into_table(
            file,
            file_info,
            "sv_consequences",
            "VARIANT_ID",
            {"TRANSCRIPT": "transcripts"},
        ),
        "file_suffixes": chromosomes,
    },
    "snvs": {
        "f": lambda file, file_info: import_chr_into_table(
            file, file_info, "snvs", "VARIANT_ID", {"TRANSCRIPT": "transcripts"}
        ),
        "file_suffixes": chromosomes,
    },
    "svs": {
        "f": lambda file, file_info: import_chr_into_table(
            file, file_info, "svs", "VARIANT_ID", {"TRANSCRIPT": "transcripts"}
        ),
        "file_suffixes": chromosomes,
    },
    "svs_ctx": {
        "f": lambda file, file_info: import_chr_into_table(
            file, file_info, "svs_ctx", "VARIANT_ID", {"TRANSCRIPT": "transcripts"}
        ),
        "file_suffixes": chromosomes,
    },
    "str": {
        "f": lambda file, file_info: import_chr_into_table(
            file, file_info, "str", "VARIANT_ID", {"TRANSCRIPT": "transcripts"}
        ),
        "file_suffixes": chromosomes,
    },
    "mts": {
        "f": lambda file, file_info: import_chr_into_table(
            file, file_info, "mts", "VARIANT_ID", {"TRANSCRIPT": "transcripts"}
        ),
        "file_suffixes": chromosomes,
    },
    "genomic_ibvl_frequencies": {
        "f": lambda file, file_info: import_chr_into_table(
            file,
            file_info,
            "genomic_ibvl_frequencies",
            "VARIANT_ID",
            {"VARIANT": "variants"},
        ),
        "file_suffixes": chromosomes,
    },
    "genomic_gnomad_frequencies": {
        "f": lambda file, file_info: import_chr_into_table(
            file,
            file_info,
            "genomic_gnomad_frequencies",
            "VARIANT_ID",
            {"VARIANT": "variants"},
        ),
        "file_suffixes": chromosomes,
    },
    "mt_ibvl_frequencies": {
        "f": lambda file, file_info: import_chr_into_table(
            file,
            file_info,
            "mt_ibvl_frequencies",
            "VARIANT_ID",
            {"VARIANT": "variants"},
        ),
        "file_suffixes": ["MT"],
    },
    "mt_gnomad_frequencies": {
        "f": lambda file, file_info: import_chr_into_table(
            file,
            file_info,
            "mt_gnomad_frequencies",
            "VARIANT_ID",
            {"VARIANT": "variants"},
        ),
        "file_suffixes": ["MT"],
    },
}


def emptyTables():
    with engine.connect() as connection:
        for modelName in dirs:
            connection.execute(text("DELETE FROM " + modelName.upper()))
    #        mydb.execute("DELETE FROM "+modelName.upper())


def try_import_directory(modelName, targetFile, import_function, counts):
    successCount = 0
    failCount = 0
    missingRefCount = 0
    if import_function == None:
        print("No import function for " + modelName)
        exit()
    #        mydb.drop(modelName, if_exists="cascade")
    #       mydb.execute(createStatements[modelName])

    if os.path.isfile(targetFile):
        print("importing into " + modelName + " (" + targetFile + ")...")
        # print(targetFile)

        file_info = inspectTSV(targetFile)
        results = import_function(targetFile, file_info)
        successCount += results["success"]
        failCount += results["fail"]
        missingRefCount += results["missingRef"]
        duplicateCount = results["duplicate"]
        percent_success = 100 * successCount / (successCount + failCount)
        print(
            "finished importing "
            + targetFile
            + "\nSuccess: "
            + str(results["success"])
            + ", Fail: "
            + str(results["fail"])
            + ". total: "
            + str(file_info["total_rows"])
            + ". missing refs: "
            + str(results["missingRef"])
            + ". Success rate: "
            + str(percent_success)
            + "%. Duplicates: "
            + str(duplicateCount)
        )
    else:
        print("File not found: " + targetFile)
    counts["successCountTotal"] += successCount
    counts["failCountTotal"] += failCount
    counts["missingRefCountTotal"] += missingRefCount
    counts["duplicatesCountTotal"] += duplicateCount


def main():
    emptyTables()
    now = datetime.now()
    counts = {}
    counts["successCountTotal"] = 0
    counts["failCountTotal"] = 0
    counts["missingRefCountTotal"] = 0
    counts["duplicatesCountTotal"] = 0
    for op in model_operation_definitions:
        modelName = op
        import_function = model_operation_definitions[op]["f"]
        print(modelName)
        print(import_function)
        file_suffixes = model_operation_definitions[op]["file_suffixes"]
        files = model_operation_definitions[op].get("files") or []
        for single_file in [
            rootDir + "/" + modelName + "/" + file_name + ".tsv" for file_name in files
        ]:
            try_import_directory(modelName, single_file, import_function, counts)
        for targetFile in [
            rootDir + "/" + modelName + "/" + modelName + "_" + suffix + ".tsv"
            for suffix in file_suffixes
        ]:
            try_import_directory(modelName, targetFile, import_function, counts)

    print(
        "Finished importing all chromosomes. Took this much time: "
        + str(datetime.now() - now)
    )
    percent_success = (
        100
        * counts["successCountTotal"]
        / (counts["successCountTotal"] + counts["failCountTotal"])
    )
    print(
        str(percent_success)
        + "% of objects were imported. Total success: "
        + str(counts["successCountTotal"])
        + ", Total fail: "
        + str(counts["failCountTotal"])
        + ". Total missing refs: "
        + str(counts["missingRefCountTotal"])
        + ". Total duplicates: "
        + str(counts["duplicatesCountTotal"])
    )
    # Read TSV
    # df = pd.read_csv('genes/genes_chr1.tsv', sep='\t')

    # Create cursor
    # cursor = mydb.cursor()

    # Create table
    # mydb.execute("CREATE TABLE genes (gene_id INT AUTO_INCREMENT PRIMARY KEY, gene_name VARCHAR(255), gene_start INT, gene_end INT, gene_strand VARCHAR(255), gene_biotype VARCHAR(255), gene_chromosome VARCHAR(255))")

    # Insert data

    # df.to_sql('genes', con=mydb, if_exists='append', index=False)


main()
