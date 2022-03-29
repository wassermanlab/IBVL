import os

# TODO: Load dotenv here

# URLs
BASE_URL = os.environ.get("DEV_BASE_API")

# Logging
LOGGING_FORMAT = "%(asctime)s - %(levelname)s - %(lineno)d - %(message)s"

# Template files
TABLES = [
    "genes",
    "variants",
    "snvs",
    "svs",
    "svs_ctx",
    "mts",
    "genomic_gnomad_frequencies",
    "genomic_ibvl_frequencies",
    "mt_gnomad_frequencies",
    "mt_ibvl_frequencies",
    "transcripts",
    "variants_transcripts",
    "variants_annotations",
    "severities",
    "variants_consequences",
    "sv_consequences"
]

TABLE_POST_DICT = {
    "genes": os.path.join(BASE_URL, "genes"),
    "variants": os.path.join(BASE_URL, "variants"),
    "snvs": os.path.join(BASE_URL, "snvs"),
    "svs": os.path.join(BASE_URL, "svs"),
    "svs_ctx": os.path.join(BASE_URL, "svs-ctx"),
    "mts": os.path.join(BASE_URL, "mts"),
    "genomic_gnomad_frequencies": os.path.join(BASE_URL, "gen-gnomad-freq"),
    "genomic_ibvl_frequencies": os.path.join(BASE_URL, "gen-ibvl-freq"),
    "mt_gnomad_frequencies": os.path.join(BASE_URL, "mt-gnomad-freq"),
    "mt_ibvl_frequencies": os.path.join(BASE_URL, "mt-ibvl-freq"),
    "transcripts": os.path.join(BASE_URL, "transcripts"),
    "variants_transcripts": os.path.join(BASE_URL, "variant-transcripts"),
    "variants_annotations": os.path.join(BASE_URL, "variant-annotations"),
    "severities": os.path.join(BASE_URL, "severity"),
    "variants_consequences": os.path.join(BASE_URL, "variant-con"),
    "sv_consequences": os.path.join(BASE_URL, "sv-con")
}