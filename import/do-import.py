import array
import string
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
import warnings
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

chromosomes = ["chrY", "chrX", "chr10"]
# chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'MT', 'str', 'sv_MEI_chr1', 'sv_MEI_chr2', 'sv_MEI_chr3', 'sv_MEI_chr4', 'sv_MEI_chr5', 'sv_MEI_chr6', 'sv_MEI_chr7', 'sv_MEI_chr8', 'sv_MEI_chr9', 'sv_MEI_chr10', 'sv_MEI_chr11', 'sv_MEI_chr12', 'sv_MEI_chr13', 'sv_MEI_chr14', 'sv_MEI_chr15', 'sv_MEI_chr16', 'sv_MEI_chr17', 'sv_MEI_chr18', 'sv_MEI_chr19', 'sv_MEI_chr20', 'sv_MEI_chr21', 'sv_MEI_chr22', 'sv_MEI_chrX', 'sv_MEI_chrY']
# Define the "GENES" table

use_n_chromosomes=0


def transform_header(header):
    # Example: Convert all column names to uppercase
    return [col.upper() for col in header]


def readTSV(file, info):
    df = pd.read_csv(file, sep=info["separator"])  # ,header=0,names=transform_header)
    df.rename(columns={"All_info$variant": "variant"}, inplace=True)
    df.columns = [col.upper() for col in df.columns]
    #  set column "a" to instead be "b"
    return df


def inspectTSV(file):
    total_rows = 0
    separator = "\t"
    small_read = pd.read_csv(file, sep="\t", nrows=4, header=None)

    num_columns = small_read.shape[1]
    if (num_columns == 1 and "gene" not in file):
        separator = " "
        small_read = pd.read_csv(file, sep=" ", nrows=4, header=None)
        if (small_read.shape[1] ==1):
            print("Could not determine separator for "+file)
            quit()
        
    columns = [col.upper() for col in small_read.values.tolist()[0]]

    for chunk in pd.read_csv(file, sep="\t", chunksize=chunk_size):
        total_rows += len(chunk)

    return {"total_rows": total_rows, "num_columns": num_columns, "columns": columns, "types": {}, "separator": separator}


pk_maps = {}
next_id_maps = {}
tables = {}

def import_chr_into_table(file, file_info, name, pk_lookup_col=None, fk_map={}):
    warnings.filterwarnings("ignore")
    df = readTSV(file, file_info)
    df.fillna("NA", inplace=True)
    warnings.resetwarnings()

#    print(file_info)

    table = tables.get(name)
    if table == None:
        columns = file_info["columns"]
        if "DO_COMPOUND_FK" in fk_map:
            columns.remove("VARIANT")
            columns.remove("TRANSCRIPT")
            columns.append("VARIANT_TRANSCRIPT")
        table = Table(
            name.upper(),
            metadata,
            Column("ID"),
            *(Column(col) for col in columns)
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
#                print("found fk map: "+fk_col+" -> "+fk_model+" in "+name)
                map_key = None
                resolved_pk = None
                if fk_col == "DO_COMPOUND_FK":
                    map_key = "|".join([ data["VARIANT"], data["TRANSCRIPT"]])
                    del data["VARIANT"]
                    del data["TRANSCRIPT"]
                    resolved_pk = resolve_PK("variants_transcripts", map_key)
                    fk_col = "VARIANT_TRANSCRIPT"
                else:
                    map_key = data[fk_col]
                if map_key == "NA":
                    data[fk_col] = None # FIXME causes missing refs
                else:
                    resolved_pk = resolve_PK(fk_model, map_key)
                if map_key != None and verbose:
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
                map_key = None
                if isinstance(pk_lookup_col,list):
                    map_key = "|".join([data[col] for col in pk_lookup_col])
                elif isinstance(pk_lookup_col, str):
                    map_key = data[pk_lookup_col]
                if pk_lookup_col != None:
                    pk_maps[name][map_key] = pk
                    if verbose:
                        print(
                            "added " + name + "." + map_key + " to pk map"
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
#        "skip":True
    },
    "transcripts": {
        "f": lambda file, file_info: import_chr_into_table(
            file, file_info, "transcripts", "TRANSCRIPT_ID", {"GENE": "genes"}
        )
    },
    "variants": {
        "f": lambda file, file_info: import_chr_into_table(
            file, file_info, "variants", "VARIANT_ID"
        ),
#        "skip":True
    },
    "variants_transcripts": {
        "f": lambda file, file_info: import_chr_into_table(
            file,
            file_info,
            "variants_transcripts",
            ["TRANSCRIPT", "VARIANT"],
            {"TRANSCRIPT": "transcripts", "VARIANT": "variants"},
        ),
#        "skip":True
    },
    "variants_annotations": {
          # RESOLVE variant AND transcript value to PK in variants_transcripts
        "f": lambda file, file_info: import_chr_into_table(
            file,
            file_info,
            "variants_annotations",
            None,
            {"DO_COMPOUND_FK":"for variants and annotations"},
        )
    },
    "variants_consequences": {
          # RESOLVE variant AND transcript value to PK in variants_transcripts
        "f": lambda file, file_info: import_chr_into_table(
            file,
            file_info,
            "variants_consequences",
            None,
            {"DO_COMPOUND_FK":"for variants and annotations"}
        )
    },
    "sv_consequences": {
        "f": lambda file, file_info: import_chr_into_table(
            file,
            file_info,
            "sv_consequences",
            None,
            {"GENE": "genes", "VARIANT": "variants"},
        ),
#        "preprocess": lambda file, file_info: convert(file, file_info),
    },
    "snvs": {
        "f": lambda file, file_info: import_chr_into_table(
            file, file_info, "snvs", None, {"VARIANT": "variants"}
        )
    },
    "svs": {
        "f": lambda file, file_info: import_chr_into_table(
            file, file_info, "svs", None, {"VARIANT": "variants"}
        )
    },
    "svs_ctx": {
        "f": lambda file, file_info: import_chr_into_table(
            file, file_info, "svs_ctx", None, {"VARIANT": "variants"}
        )
    },
    "str": {
        "f": lambda file, file_info: import_chr_into_table(
            file, file_info, "str", None, {"VARIANT": "variants"}
        )
    },
    "mts": {
        "f": lambda file, file_info: import_chr_into_table(
            file, file_info, "mts", None, {"VARIANT": "variants"}
        )
    },
    "genomic_ibvl_frequencies": {
        "f": lambda file, file_info: import_chr_into_table(
            file,
            file_info,
            "genomic_ibvl_frequencies",
            None,
            {"VARIANT": "variants"},
        )
    },
    "genomic_gnomad_frequencies": {
        "f": lambda file, file_info: import_chr_into_table(
            file,
            file_info,
            "genomic_gnomad_frequencies",
            None,
            {"VARIANT": "variants"},
        )
    },
    "mt_ibvl_frequencies": {
        "f": lambda file, file_info: import_chr_into_table(
            file,
            file_info,
            "mt_ibvl_frequencies",
            None,
            {"VARIANT": "variants"},
        )
    },
    "mt_gnomad_frequencies": {
        "f": lambda file, file_info: import_chr_into_table(
            file,
            file_info,
            "mt_gnomad_frequencies",
            None,
            {"VARIANT": "variants"},
        )
    },
}


def try_import_file(modelName, targetFile, import_function, counts):
    successCount = 0
    failCount = 0
    missingRefCount = 0
    duplicateCount = 0
    if import_function == None:
        print("No import function for " + modelName)
        exit()
    #        mydb.drop(modelName, if_exists="cascade")
    #       mydb.execute(createStatements[modelName])

    if os.path.isfile(targetFile):
        print("\nimporting into " + modelName + " (" + targetFile.split("/")[-1] + ")...")
        # print(targetFile)

        file_info = inspectTSV(targetFile)
        results = import_function(targetFile, file_info)
        successCount += results["success"]
        failCount += results["fail"]
        missingRefCount += results["missingRef"]
        duplicateCount = results["duplicate"]

        percent_success = "N/A"
        if successCount + failCount > 0:
            percent_success = 100 * successCount / (successCount + failCount)
        print(
            "finished importing "
            + targetFile.split("/")[-1]
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

def report_counts(counts):
    percent_success = "N/A"
    if counts["successCountTotal"] + counts["failCountTotal"] > 0:
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
def main():
    now = datetime.now()
    counts = {}
    counts["successCountTotal"] = 0
    counts["failCountTotal"] = 0
    counts["missingRefCountTotal"] = 0
    counts["duplicatesCountTotal"] = 0
    for modelName in model_operation_definitions:
        if(model_operation_definitions[modelName].get("skip")):
            continue
        modelNow = datetime.now()
        # TODO: consider running multiple times / retrying
        pk_maps[modelName] = {}
        next_id_maps[modelName] = 1
        with engine.connect() as connection:
            connection.execute(text("DELETE FROM " + modelName.upper()))
        modelName = modelName
        import_function = model_operation_definitions[modelName]["f"]
        # get files from folder
        #sort alphabetically but with numbers in order
        sorted_files = sorted(
            os.listdir(rootDir + "/" + modelName),
            #key=lambda x: x.reverse(),
            reverse=True
        )
#        sorted_files = sorted_files[:use_n_chromosomes]
        for file in sorted_files:
            if file.endswith(".tsv"):
                file_path = rootDir + "/" + modelName + "/" + file
       #     rootDir + "/" + modelName + "/" + file_name + ".tsv" for file_name in files
                try_import_file(modelName, file_path, import_function, counts)

        print(
            "Finished importing "+modelName+". Took this much time: "
            + str(datetime.now() - modelNow)
        )
        report_counts(counts)
    print("finished importing IBVL. Time Taken: " + str(datetime.now() - now))
    report_counts(counts)

main()
