import os

# TODO: Load dotenv here

# URLs
DEV_BASE_URL = os.environ.get("DEV_BASE_API")
PROD_BASE_URL = os.environ.get("PROD_BASE_API")

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

DEV_TABLE_POST_DICT = {
    "genes": os.path.join(DEV_BASE_URL, "genes"),
    "variants": os.path.join(DEV_BASE_URL, "variants"),
    "snvs": os.path.join(DEV_BASE_URL, "snvs"),
    "svs": os.path.join(DEV_BASE_URL, "svs"),
    "svs_ctx": os.path.join(DEV_BASE_URL, "svs-ctx"),
    "mts": os.path.join(DEV_BASE_URL, "mts"),
    "genomic_gnomad_frequencies": os.path.join(DEV_BASE_URL, "gen-gnomad-freq"),
    "genomic_ibvl_frequencies": os.path.join(DEV_BASE_URL, "gen-ibvl-freq"),
    "mt_gnomad_frequencies": os.path.join(DEV_BASE_URL, "mt-gnomad-freq"),
    "mt_ibvl_frequencies": os.path.join(DEV_BASE_URL, "mt-ibvl-freq"),
    "transcripts": os.path.join(DEV_BASE_URL, "transcripts"),
    "variants_transcripts": os.path.join(DEV_BASE_URL, "variant-transcripts"),
    "variants_annotations": os.path.join(DEV_BASE_URL, "variant-annotations"),
    "severities": os.path.join(DEV_BASE_URL, "severity"),
    "variants_consequences": os.path.join(DEV_BASE_URL, "variant-con"),
    "sv_consequences": os.path.join(DEV_BASE_URL, "sv-con")
}

PROD_TABLE_POST_DICT = {
    "genes": os.path.join(PROD_BASE_URL, "genes"),
    "variants": os.path.join(PROD_BASE_URL, "variants"),
    "snvs": os.path.join(PROD_BASE_URL, "snvs"),
    "svs": os.path.join(PROD_BASE_URL, "svs"),
    "svs_ctx": os.path.join(PROD_BASE_URL, "svs-ctx"),
    "mts": os.path.join(PROD_BASE_URL, "mts"),
    "genomic_gnomad_frequencies": os.path.join(PROD_BASE_URL, "gen-gnomad-freq"),
    "genomic_ibvl_frequencies": os.path.join(PROD_BASE_URL, "gen-ibvl-freq"),
    "mt_gnomad_frequencies": os.path.join(PROD_BASE_URL, "mt-gnomad-freq"),
    "mt_ibvl_frequencies": os.path.join(PROD_BASE_URL, "mt-ibvl-freq"),
    "transcripts": os.path.join(PROD_BASE_URL, "transcripts"),
    "variants_transcripts": os.path.join(PROD_BASE_URL, "variant-transcripts"),
    "variants_annotations": os.path.join(PROD_BASE_URL, "variant-annotations"),
    "severities": os.path.join(PROD_BASE_URL, "severity"),
    "variants_consequences": os.path.join(PROD_BASE_URL, "variant-con"),
    "sv_consequences": os.path.join(PROD_BASE_URL, "sv-con")
}

ID_MAPPER = {
    "genes": ["short_name"],
    "variants": ["variant_id"],
    "snvs": ["variant"],
    "svs": ["variant"],
    "svs_ctx": ["variant"],
    "mts": ["variant"],
    "genomic_gnomad_frequencies": ["variant"],
    "genomic_ibvl_frequencies": ["variant"],
    "mt_gnomad_frequencies": ["variant"],
    "mt_ibvl_frequencies": ["variant"],
    "transcripts": ["transcript_id"],
    "variants_transcripts": ["transcript", "variant"],
    "variants_annotations": ["transcript", "variant"],
    "severities": ["severity_number"],
    "variants_consequences": ["transcript", "variant"],
    "sv_consequences": ["variant"],
}