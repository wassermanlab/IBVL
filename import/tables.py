from sqlalchemy import (
    create_engine,
    MetaData,
    Table,
    Column,
    Integer,
    String,
    Float,
    UniqueConstraint,
    ForeignKey,
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
# Replace 'sqlite:///:memory:' with your actual database connection string

if (False and verbose):
    engine = create_engine(dbConnectionString, echo=True, pool_pre_ping=True, pool_recycle=3600, connect_args={'autocommit': True})
else:
    engine = create_engine(dbConnectionString, pool_pre_ping=True, pool_recycle=3600, connect_args={'autocommit': True})

metadata = MetaData()

# Define the "GENES" table
genes_table = Table(
    "GENES",
    metadata,
    Column("ID", Integer, primary_key=True),
    Column("SHORT_NAME", String(30), nullable=False),
    UniqueConstraint("SHORT_NAME", name="UNIQUE"),
)


transcripts_table = Table(
    "TRANSCRIPTS",
    metadata,
    Column("TRANSCRIPT_ID", String(255)),
    Column("ID", Integer, primary_key=True),
    Column("GENE", Integer, ForeignKey("GENES.ID", ondelete="CASCADE")),
    Column("TRANSCRIPT_TYPE", String(1)),
    Column("TSL", String(255)),
    UniqueConstraint("TRANSCRIPT_ID", name="TRANSCRIPTS_UNIQUE"),
#    ForeignKey("GENE", name="TRANSCRIPTS_GENES_FK", ondelete="CASCADE", table="GENES"),
)

# Define the "GENOMIC_GNOMAD_FREQUENCIES" table
genomic_gnomad_frequencies_table = Table(
    "GENOMIC_GNOMAD_FREQUENCIES",
    metadata,
    Column("ID", Integer, primary_key=True),
    Column("VARIANT", Integer, ForeignKey("VARIANTS.ID", ondelete="CASCADE")),
    Column("AF_TOT", Float),
    Column("AC_TOT", Integer),
    Column("AN_TOT", Integer),
    Column("HOM_TOT", Integer),
    UniqueConstraint("VARIANT", name="GENO_GNOMAD_FREQ_UNIQUE"),
)

# Define the "GENOMIC_IBVL_FREQUENCIES" table
genomic_ibvl_frequencies_table = Table(
    "GENOMIC_IBVL_FREQUENCIES",
    metadata,
    Column("ID", Integer, primary_key=True),
    Column("VARIANT", Integer, ForeignKey("VARIANTS.ID", ondelete="CASCADE")),
    Column("AF_TOT", Float),
    Column("AF_XX", Float),
    Column("AF_XY", Float),
    Column("AC_TOT", Integer),
    Column("AC_XX", Integer),
    Column("AC_XY", Integer),
    Column("AN_TOT", Integer),
    Column("AN_XX", Integer),
    Column("AN_XY", Integer),
    Column("HOM_TOT", Integer),
    Column("HOM_XX", Integer),
    Column("HOM_XY", Integer),
    Column("QUALITY", Integer),
    UniqueConstraint("VARIANT", name="GENO_IBVL_FREQ_UNIQUE")
)

mt_gnomad_frequencies_table = Table(
    "MT_GNOMAD_FREQUENCIES",
    metadata,
    Column("ID", Integer, primary_key=True),
    Column("VARIANT", Integer, ForeignKey("VARIANTS.ID", ondelete="CASCADE")),
    Column("AN", Integer),
    Column("AC_HOM", Integer),
    Column("AC_HET", Integer),
    Column("AF_HOM", Float),
    Column("AF_HET", Float),
    Column("MAX_HL", Integer),
    UniqueConstraint("VARIANT", name="MT_GNOMAD_FREQ_UNIQUE")
)

mt_ibvl_frequencies_table = Table(
    "MT_IBVL_FREQUENCIES",
    metadata,
    Column("ID", Integer, primary_key=True),
    Column("VARIANT", Integer, ForeignKey("VARIANTS.ID", ondelete="CASCADE")),
    Column("AN", Integer),
    Column("AC_HOM", Integer),
    Column("AC_HET", Integer),
    Column("AF_HOM", Float),
    Column("AF_HET", Float),
    Column("HL_HIST", String(30)),
    Column("MAX_HL", Integer),
    UniqueConstraint("VARIANT", name="MT_IBVL_FREQ_UNIQUE"),
)


mts_table = Table(
    "MTS",
    metadata,
    Column("VARIANT", Integer, ForeignKey("VARIANTS.ID", ondelete="CASCADE")),
    Column("ID", Integer, primary_key=True),
    Column("POS", Integer),
    Column("REF", String(1)),
    Column("ALT", String(30)),
    Column("UCSC_URL", String(511)),
    Column("MITOMAP_URL", String(511)),
    Column("GNOMAD_URL", String(511)),
    Column("DBSNP_ID", String(30)),
    Column("DBSNP_URL", String(511)),
    Column("CLINVAR_URL", String(511)),
    Column("CLINVAR_VCV", Integer),
    UniqueConstraint("VARIANT", name="MTS_UNIQUE"),
)

snvs_table = Table(
    "SNVS",
    metadata,
    Column("VARIANT", Integer, ForeignKey("VARIANTS.ID", ondelete="CASCADE")),
    Column("ID", Integer, primary_key=True),
    Column("TYPE", String(30)),
    Column("LENGTH", Integer),
    Column("CHR", Integer),
    Column("POS", Integer),
    Column("REF", String(255)),
    Column("ALT", String(255)),
    Column("CADD_INTR", String(255)),
    Column("CADD_SCORE", Integer),
    Column("DBSNP_ID", String(30)),
    Column("DBSNP_URL", String(511)),
    Column("UCSC_URL", String(511)),
    Column("ENSEMBL_URL", String(511)),
    Column("CLINVAR_VCV", Integer),
    Column("CLINVAR_URL", String(511)),
    Column("GNOMAD_URL", String(511)),
    Column("SPLICE_AI", Integer),
    UniqueConstraint("VARIANT", name="SNVS_UNIQUE")
)

# Define the "STR" table
str_table = Table(
    "STR",
    metadata,
    Column("VARIANT", Integer, ForeignKey("VARIANTS.ID", ondelete="CASCADE")),
    Column("ID", Integer, primary_key=True),
    Column("REPEAT_UNIT", String(20)),
    Column("MIN_N_REPEAT", Integer),
    Column("MAX_N_REPEAT", Integer),
    Column("ALLELE_DISTRIB", String(255)),
    Column("REFERENCE_REGION", String(40)),
    UniqueConstraint("VARIANT", name="STR_UNIQUE")
)

# Define the "SV_CONSEQUENCES" table
sv_consequences_table = Table(
    "SV_CONSEQUENCES",
    metadata,
    Column("ID", Integer, primary_key=True),
    Column("GENE", Integer, ForeignKey("GENES.ID", ondelete="CASCADE")),
    Column("VARIANT", Integer, ForeignKey("VARIANTS.ID", ondelete="CASCADE")),
    Column("CONSEQUENCE", String(255)),
    UniqueConstraint("VARIANT", name="SV_CONSEQUENCES_UNIQUE")
)

svs_table = Table(
    "SVS",
    metadata,
    Column("VARIANT", Integer, ForeignKey("VARIANTS.ID", ondelete="CASCADE")),
    Column("ID", Integer, primary_key=True),
    Column("CHR1", Integer),
    Column("CHR1_POS1", Integer),
    Column("CHR1_POS2", Integer),
    Column("SV_TYPE", String(30)),
    Column("SV_LENGTH", Integer),
    Column("ALGORITHM", String(30)),
    Column("UCSC_URL", String(511)),
    Column("GNOMAD_ID", String(30)),
    Column("GNOMAD_URL", String(511)),
    UniqueConstraint("VARIANT", name="SVS_UNIQUE"),
#    ForeignKey("VARIANT", name="SVS_FK", ondelete="CASCADE", table="VARIANTS"),
)

svs_ctx_table = Table(
    "SVS_CTX",
    metadata,
    Column("VARIANT", Integer, ForeignKey("VARIANTS.ID", ondelete="CASCADE")),
    Column("ID", Integer, primary_key=True),
    Column("CHR2", Integer),
    Column("CHR2_POS1", Integer),
    Column("UCSC_URL2", String(511)),
    Column("GNOMAD_ID2", String(30)),
    Column("GNOMAD_URL2", String(511)),
    UniqueConstraint("VARIANT", name="SVS_CTX_UNIQUE"),
#    ForeignKey("VARIANT", name="SVS_CTX_FK", ondelete="CASCADE", table="VARIANTS"),
)


variants_table = Table(
    "VARIANTS",
    metadata,
    Column("VARIANT_ID", String(255)),
    Column("ID", Integer, primary_key=True),
    Column("VAR_TYPE", String(30)),
    UniqueConstraint("VARIANT_ID", name="VARIANTS_UNIQUE"),
)

variants_annotations_table = Table(
    "VARIANTS_ANNOTATIONS",
    metadata,
    Column("ID", Integer, primary_key=True),
    Column("HGVSP", String(255)),
    Column("POLYPHEN", String(255)),
    Column("SIFT", String(255)),
    Column("VARIANT_TRANSCRIPT", Integer, ForeignKey("VARIANTS_TRANSCRIPTS.ID", ondelete="CASCADE")),
    UniqueConstraint("VARIANT_TRANSCRIPT", name="VARIANTS_ANNOTATIONS_UNIQUE"),
)
# ... (Previous code)

# Define the "VARIANTS_CONSEQUENCES" table
variants_consequences_table = Table(
    "VARIANTS_CONSEQUENCES",
    metadata,
    Column("ID", Integer, primary_key=True),
    Column("SEVERITY", Integer, ForeignKey("SEVERITIES.ID", ondelete="CASCADE")),
    Column("VARIANT_TRANSCRIPT", Integer, ForeignKey("VARIANTS_TRANSCRIPTS.ID", ondelete="CASCADE")),
    UniqueConstraint("VARIANT_TRANSCRIPT", name="VARIANTS_CONSEQUENCES_UNIQUE"),
)


variants_transcripts_table = Table(
    "VARIANTS_TRANSCRIPTS",
    metadata,
    Column("TRANSCRIPT", Integer, ForeignKey("TRANSCRIPTS.ID", ondelete="CASCADE")),
    Column("ID", Integer, primary_key=True),
    Column("VARIANT", Integer, ForeignKey("VARIANTS.ID", ondelete="CASCADE")),
    Column("HGVSC", String(255)),
    UniqueConstraint("TRANSCRIPT", "VARIANT", name="VARIANTS_TRANSCRIPTS_UNIQUE"),
)

severities_table = Table(
    "SEVERITIES",
    metadata,
    Column("SEVERITY_NUMBER", Integer),
    Column("ID", Integer, primary_key=True),
    Column("CONSEQUENCE", String(255)),
    UniqueConstraint("SEVERITY_NUMBER", name="SEVERITIES_UNIQUE"),
)

metadata.create_all(engine)