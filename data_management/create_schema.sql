ALTER TABLE SNVS DROP CONSTRAINT SNVS_FK;
ALTER TABLE SVS DROP CONSTRAINT SVS_FK;
ALTER TABLE SVS_CTX DROP CONSTRAINT SVS_CTX_FK;
ALTER TABLE STR DROP CONSTRAINT STR_FK;
ALTER TABLE MTS DROP CONSTRAINT MTS_FK;
ALTER TABLE GENOMIC_GNOMAD_FREQUENCIES DROP CONSTRAINT GENOMIC_GNOMAD_FREQUENCIES_FK;
ALTER TABLE GENOMIC_IBVL_FREQUENCIES DROP CONSTRAINT GENOMIC_IBVL_FREQUENCIES_FK;
ALTER TABLE MT_GNOMAD_FREQUENCIES DROP CONSTRAINT MT_GNOMAD_FREQUENCIES_FK;
ALTER TABLE MT_IBVL_FREQUENCIES DROP CONSTRAINT MT_IBVL_FREQUENCIES_FK;
ALTER TABLE TRANSCRIPTS DROP CONSTRAINT TRANSCRIPTS_GENES_FK;
ALTER TABLE VARIANTS_TRANSCRIPTS DROP CONSTRAINT VARIANTS_TRANSCRIPTS_FK;
ALTER TABLE VARIANTS_TRANSCRIPTS DROP CONSTRAINT VARIANTS_TRANSCRIPTS_FK2;
ALTER TABLE VARIANTS_TRANSCRIPTS DROP CONSTRAINT VARIANTS_TRANSCRIPTS_UNIQUE;
ALTER TABLE VARIANTS_ANNOTATIONS DROP CONSTRAINT VARIANTS_ANNOTATIONS_FK;
ALTER TABLE VARIANTS_CONSEQUENCES DROP CONSTRAINT VARIANTS_CONSEQUENCES_FK;
ALTER TABLE VARIANTS_CONSEQUENCES DROP CONSTRAINT VARIANTS_CONSEQUENCES_FK2;
ALTER TABLE SV_CONSEQUENCES DROP CONSTRAINT SV_CONSEQUENCES_VAR_FK;
ALTER TABLE SV_CONSEQUENCES DROP CONSTRAINT SV_CONSEQUENCES_GENE_FK;
ALTER TABLE GENES DROP CONSTRAINT GENES_PK;
ALTER TABLE VARIANTS DROP CONSTRAINT VARIANTS_PK;
ALTER TABLE SNVS DROP CONSTRAINT SNVS_PK;
ALTER TABLE SVS DROP CONSTRAINT SVS_PK;
ALTER TABLE SVS_CTX DROP CONSTRAINT SVS_CTX_PK;
ALTER TABLE STR DROP CONSTRAINT STR_PK;
ALTER TABLE MTS DROP CONSTRAINT MTS_PK;
ALTER TABLE GENOMIC_GNOMAD_FREQUENCIES DROP CONSTRAINT GENOMIC_GNOMAD_FREQUENCIES_PK;
ALTER TABLE GENOMIC_IBVL_FREQUENCIES DROP CONSTRAINT GENOMIC_IBVL_FREQUENCIES_PK;
ALTER TABLE MT_GNOMAD_FREQUENCIES DROP CONSTRAINT MT_GNOMAD_FREQUENCIES_PK;
ALTER TABLE MT_IBVL_FREQUENCIES DROP CONSTRAINT MT_IBVL_FREQUENCIES_PK;
ALTER TABLE TRANSCRIPTS DROP CONSTRAINT TRANSCRIPTS_PK;
ALTER TABLE VARIANTS_TRANSCRIPTS DROP CONSTRAINT VARIANTS_TRANSCRIPTS_PK;
ALTER TABLE VARIANTS_ANNOTATIONS DROP CONSTRAINT VARIANTS_ANNOTATIONS_PK;
ALTER TABLE SEVERITIES DROP CONSTRAINT SEVERITIES_PK;
ALTER TABLE VARIANTS_CONSEQUENCES DROP CONSTRAINT VARIANTS_CONSEQUENCES_PK;
ALTER TABLE SV_CONSEQUENCES DROP CONSTRAINT SV_CONSEQUENCES_PK;

DROP SEQUENCE GENES_SEQ;
DROP SEQUENCE VARIANTS_SEQ;
DROP SEQUENCE SNVS_SEQ;
DROP SEQUENCE SVS_SEQ;
DROP SEQUENCE SVS_CTX_SEQ;
DROP SEQUENCE STR_SEQ;
DROP SEQUENCE MTS_SEQ;
DROP SEQUENCE GENOMIC_GNOMAD_FREQUENCIES_SEQ;
DROP SEQUENCE GENOMIC_IBVL_FREQUENCIES_SEQ;
DROP SEQUENCE MT_GNOMAD_FREQUENCIES_SEQ;
DROP SEQUENCE MT_IBVL_FREQUENCIES_SEQ;
DROP SEQUENCE TRANSCRIPTS_SEQ;
DROP SEQUENCE VARIANTS_TRANSCRIPTS_SEQ;
DROP SEQUENCE VARIANTS_ANNOTATIONS_SEQ;
DROP SEQUENCE SEVERITIES_SEQ;
DROP SEQUENCE VARIANTS_CONSEQUENCES_SEQ;
DROP SEQUENCE SV_CONSEQUENCES_SEQ;



DROP TABLE GENES;
DROP TABLE VARIANTS;
DROP TABLE SNVS;
DROP TABLE SVS;
DROP TABLE SVS_CTX;
DROP TABLE STR;
DROP TABLE MTS;
DROP TABLE GENOMIC_GNOMAD_FREQUENCIES;
DROP TABLE GENOMIC_IBVL_FREQUENCIES;
DROP TABLE MT_GNOMAD_FREQUENCIES;
DROP TABLE MT_IBVL_FREQUENCIES;
DROP TABLE TRANSCRIPTS;
DROP TABLE VARIANTS_TRANSCRIPTS;
DROP TABLE VARIANTS_ANNOTATIONS;
DROP TABLE SEVERITIES;
DROP TABLE VARIANTS_CONSEQUENCES;
DROP TABLE SV_CONSEQUENCES;

-- Create Genes
CREATE SEQUENCE GENES_SEQ;

CREATE TABLE GENES 
   (	"ID" NUMBER NOT NULL ENABLE,  
	"SHORT_NAME" VARCHAR2(30) NOT NULL ENABLE UNIQUE, 
	 CONSTRAINT "GENES_PK" PRIMARY KEY ("ID")
  USING INDEX  ENABLE
   );

-- Create Variants
CREATE SEQUENCE VARIANTS_SEQ;

CREATE TABLE VARIANTS 
   (	"VARIANT_ID" VARCHAR2(100) UNIQUE, 
        "ID" NUMBER NOT NULL ENABLE,
	    "VAR_TYPE" VARCHAR2(30), 
	    CONSTRAINT "VARIANTS_PK" PRIMARY KEY ("ID") USING INDEX ENABLE
   );

-- Create SNVs
CREATE SEQUENCE SNVS_SEQ;


CREATE TABLE SNVS
   (	"VARIANT" NUMBER UNIQUE, 
    "ID" NUMBER NOT NULL ENABLE,
	"TYPE" VARCHAR2(30), 
	"LENGTH" NUMBER, 
	"CHR" VARCHAR2(2), 
	"POS" NUMBER, 
	"REF" VARCHAR2(60), 
	"ALT" VARCHAR2(60), 
    "CADD_INTR" VARCHAR2(255),
	"CADD_SCORE" NUMBER, 
	"DBSNP_ID" VARCHAR2(30), 
	"DBSNP_URL" VARCHAR2(511), 
	"UCSC_URL" VARCHAR2(511), 
	"ENSEMBL_URL" VARCHAR2(511), 
    "CLINVAR_VCV" NUMBER,
	"CLINVAR_URL" VARCHAR2(511), 
	"GNOMAD_URL" VARCHAR2(511), 
    "SPLICE_AI" NUMBER,
	 CONSTRAINT "SNVS_PK" PRIMARY KEY ("ID")
  USING INDEX  ENABLE
   );

ALTER TABLE  "SNVS" ADD CONSTRAINT "SNVS_FK" FOREIGN KEY ("VARIANT")
	  REFERENCES  "VARIANTS" ("ID") ON DELETE CASCADE ENABLE;


-- Create SVs
CREATE SEQUENCE SVS_SEQ;

CREATE TABLE SVS
   (	"VARIANT" NUMBER UNIQUE,
    "ID" NUMBER NOT NULL ENABLE, 
	"CHR1" VARCHAR2(2), 
	"CHR1_POS1" NUMBER, 
	"CHR1_POS2" NUMBER, 
	"SV_TYPE" VARCHAR2(30), 
	"SV_LENGTH" NUMBER, 
	"ALGORITHM" VARCHAR2(100), 
	"UCSC_URL" VARCHAR2(511), 
	"GNOMAD_ID" VARCHAR2(100), 
	"GNOMAD_URL" VARCHAR2(511), 
	 CONSTRAINT "SVS_PK" PRIMARY KEY ("ID")
  USING INDEX  ENABLE
   );

ALTER TABLE  "SVS" ADD CONSTRAINT "SVS_FK" FOREIGN KEY ("VARIANT")
	  REFERENCES  "VARIANTS" ("ID") ON DELETE CASCADE ENABLE;

-- Create SVs CTX
CREATE SEQUENCE SVS_CTX_SEQ;

CREATE TABLE SVS_CTX
   (	"VARIANT" NUMBER UNIQUE,
    "ID" NUMBER NOT NULL ENABLE, 
	"CHR2" VARCHAR2(2), 
	"CHR2_POS1" NUMBER, 
	"UCSC_URL2" VARCHAR2(511), 
    "GNOMAD_ID2" VARCHAR2(100),
    "GNOMAD_URL2" VARCHAR2(511),
	 CONSTRAINT "SVS_CTX_PK" PRIMARY KEY ("ID")
  USING INDEX  ENABLE
   );

ALTER TABLE  "SVS_CTX" ADD CONSTRAINT "SVS_CTX_FK" FOREIGN KEY ("VARIANT")
	  REFERENCES  "VARIANTS" ("ID") ON DELETE CASCADE ENABLE;

-- Create STR
CREATE SEQUENCE STR_SEQ;

CREATE TABLE STR 
   (	"VARIANT" NUMBER UNIQUE,
    "ID" NUMBER NOT NULL ENABLE, 
	"REPEAT_UNIT" VARCHAR2(20), 
	"MIN_N_REPEAT" NUMBER, 
	"MAX_N_REPEAT" NUMBER, 
    "ALLELE_DISTRIB" VARCHAR2(255),
    "REFERENCE_REGION" VARCHAR2(40),
	 CONSTRAINT "STR_PK" PRIMARY KEY ("ID")
  USING INDEX  ENABLE
   );

ALTER TABLE  "STR" ADD CONSTRAINT "STR_FK" FOREIGN KEY ("VARIANT")
	  REFERENCES  "VARIANTS" ("ID") ON DELETE CASCADE ENABLE;

-- Create Mts
CREATE SEQUENCE MTS_SEQ;

CREATE TABLE MTS
   (	"VARIANT" NUMBER UNIQUE, 
    "ID" NUMBER NOT NULL ENABLE,
	"POS" NUMBER, 
	"REF" VARCHAR2(60), 
	"ALT" VARCHAR2(60), 
	"UCSC_URL" VARCHAR2(511), 
	"MITOMAP_URL" VARCHAR2(511), 
	"GNOMAD_URL" VARCHAR2(511), 
	"DBSNP_ID" VARCHAR2(30), 
	"DBSNP_URL" VARCHAR2(511), 
	"CLINVAR_URL" VARCHAR2(511), 
    "CLINVAR_VCV" NUMBER,
	 CONSTRAINT "MTS_PK" PRIMARY KEY ("ID")
  USING INDEX  ENABLE
   );

ALTER TABLE  "MTS" ADD CONSTRAINT "MTS_FK" FOREIGN KEY ("VARIANT")
	  REFERENCES  "VARIANTS" ("ID") ON DELETE CASCADE ENABLE;

-- Create Genomic gnomAD Frequencies
CREATE SEQUENCE GENOMIC_GNOMAD_FREQUENCIES_SEQ;


CREATE TABLE GENOMIC_GNOMAD_FREQUENCIES
   (	"ID" NUMBER NOT NULL ENABLE, 
	"VARIANT" NUMBER UNIQUE, 
	"AF_TOT" NUMBER, 
	"AC_TOT" NUMBER, 
	"AN_TOT" NUMBER, 
	"HOM_TOT" NUMBER, 
	 CONSTRAINT "GENOMIC_GNOMAD_FREQUENCIES_PK" PRIMARY KEY ("ID")
  USING INDEX  ENABLE
   );

ALTER TABLE  "GENOMIC_GNOMAD_FREQUENCIES" ADD CONSTRAINT "GENOMIC_GNOMAD_FREQUENCIES_FK" FOREIGN KEY ("VARIANT")
	  REFERENCES  "VARIANTS" ("ID") ON DELETE CASCADE ENABLE;


-- Create Genomic IBVL Frequencies
CREATE SEQUENCE GENOMIC_IBVL_FREQUENCIES_SEQ;


CREATE TABLE GENOMIC_IBVL_FREQUENCIES 
   (	"ID" NUMBER NOT NULL ENABLE, 
	"VARIANT" NUMBER UNIQUE, 
	"AF_TOT" NUMBER, 
	"AF_XX" NUMBER, 
	"AF_XY" NUMBER, 
	"AC_TOT" NUMBER, 
	"AC_XX" NUMBER, 
	"AC_XY" NUMBER, 
	"AN_TOT" NUMBER, 
	"AN_XX" NUMBER, 
	"AN_XY" NUMBER, 
	"HOM_TOT" NUMBER, 
	"HOM_XX" NUMBER, 
	"HOM_XY" NUMBER, 
	"QUALITY" NUMBER, 
	 CONSTRAINT "GENOMIC_IBVL_FREQUENCIES_PK" PRIMARY KEY ("ID")
  USING INDEX  ENABLE
   );

ALTER TABLE  "GENOMIC_IBVL_FREQUENCIES" ADD CONSTRAINT "GENOMIC_IBVL_FREQUENCIES_FK" FOREIGN KEY ("VARIANT")
	  REFERENCES  "VARIANTS" ("ID") ON DELETE CASCADE ENABLE;


-- Create Mt gnomAD Frequencies
CREATE SEQUENCE MT_GNOMAD_FREQUENCIES_SEQ;


CREATE TABLE MT_GNOMAD_FREQUENCIES 
   (	"ID" NUMBER NOT NULL ENABLE, 
	"VARIANT" NUMBER UNIQUE, 
	"AN" NUMBER, 
	"AC_HOM" NUMBER, 
	"AC_HET" NUMBER, 
	"AF_HOM" NUMBER, 
	"AF_HET" NUMBER, 
	"MAX_HL" NUMBER, 
	 CONSTRAINT "MT_GNOMAD_FREQUENCIES_PK" PRIMARY KEY ("ID")
  USING INDEX  ENABLE
   );

ALTER TABLE  "MT_GNOMAD_FREQUENCIES" ADD CONSTRAINT "MT_GNOMAD_FREQUENCIES_FK" FOREIGN KEY ("VARIANT")
	  REFERENCES  "VARIANTS" ("ID") ON DELETE CASCADE ENABLE;

-- Create Mt IBVL Frequencies
CREATE SEQUENCE MT_IBVL_FREQUENCIES_SEQ;


CREATE TABLE MT_IBVL_FREQUENCIES 
   (	"ID" NUMBER NOT NULL ENABLE,
	"VARIANT" NUMBER UNIQUE, 
	"AN" NUMBER, 
	"AC_HOM" NUMBER, 
	"AC_HET" NUMBER, 
	"AF_HOM" NUMBER, 
	"AF_HET" NUMBER, 
	"HL_HIST" VARCHAR2(30), 
	"MAX_HL" NUMBER, 
	 CONSTRAINT "MT_IBVL_FREQUENCIES_PK" PRIMARY KEY ("ID")
  USING INDEX  ENABLE
   );

ALTER TABLE  "MT_IBVL_FREQUENCIES" ADD CONSTRAINT "MT_IBVL_FREQUENCIES_FK" FOREIGN KEY ("VARIANT")
	  REFERENCES  "VARIANTS" ("ID") ON DELETE CASCADE ENABLE;

-- Create Transcripts
CREATE SEQUENCE TRANSCRIPTS_SEQ;

CREATE TABLE  TRANSCRIPTS
   (	"TRANSCRIPT_ID" VARCHAR2(100) UNIQUE,
    "ID" NUMBER NOT NULL ENABLE, 
	"GENE" NUMBER, 
	"TRANSCRIPT_TYPE" VARCHAR2(1), 
    "TSL" VARCHAR2(12),
	 CONSTRAINT "TRANSCRIPTS_PK" PRIMARY KEY ("ID")
  USING INDEX  ENABLE
   );

ALTER TABLE  "TRANSCRIPTS" ADD CONSTRAINT "TRANSCRIPTS_GENES_FK" FOREIGN KEY ("GENE")
	  REFERENCES  "GENES" ("ID") ON DELETE CASCADE ENABLE;

-- Create Variants Transcripts
CREATE SEQUENCE VARIANTS_TRANSCRIPTS_SEQ;


CREATE TABLE  VARIANTS_TRANSCRIPTS
   (	"TRANSCRIPT" NUMBER, 
    "ID" NUMBER NOT NULL ENABLE,
	"VARIANT" NUMBER, 
	"HGVSC" VARCHAR2(100), 
	 CONSTRAINT "VARIANTS_TRANSCRIPTS_PK" PRIMARY KEY ("ID")
  USING INDEX  ENABLE
   );

ALTER TABLE  "VARIANTS_TRANSCRIPTS" ADD CONSTRAINT "VARIANTS_TRANSCRIPTS_FK" FOREIGN KEY ("TRANSCRIPT")
	  REFERENCES  "TRANSCRIPTS" ("ID") ON DELETE CASCADE ENABLE;

ALTER TABLE  "VARIANTS_TRANSCRIPTS" ADD CONSTRAINT "VARIANTS_TRANSCRIPTS_FK2" FOREIGN KEY ("VARIANT")
	  REFERENCES  "VARIANTS" ("ID") ON DELETE CASCADE ENABLE;

ALTER TABLE "VARIANTS_TRANSCRIPTS" ADD CONSTRAINT "VARIANTS_TRANSCRIPTS_UNIQUE" UNIQUE("TRANSCRIPT", "VARIANT");


-- Create Variants Annotations
CREATE SEQUENCE VARIANTS_ANNOTATIONS_SEQ;


CREATE TABLE VARIANTS_ANNOTATIONS 
   (	"ID" NUMBER NOT NULL ENABLE, 
	"HGVSP" VARCHAR2(100), 
	"POLYPHEN" VARCHAR2(100), 
	"SIFT" VARCHAR2(100), 
	"VARIANT_TRANSCRIPT" NUMBER UNIQUE, 
	 CONSTRAINT "VARIANTS_ANNOTATIONS_PK" PRIMARY KEY ("ID")
  USING INDEX  ENABLE
   );

ALTER TABLE  "VARIANTS_ANNOTATIONS" ADD CONSTRAINT "VARIANTS_ANNOTATIONS_FK" FOREIGN KEY ("VARIANT_TRANSCRIPT")
	  REFERENCES  "VARIANTS_TRANSCRIPTS" ("ID") ON DELETE CASCADE ENABLE;


-- Create Severities
CREATE SEQUENCE SEVERITIES_SEQ;


CREATE TABLE SEVERITIES
   (	"SEVERITY_NUMBER" NUMBER UNIQUE, 
    "ID" NUMBER NOT NULL ENABLE,
	"CONSEQUENCE" VARCHAR2(40), 
	 CONSTRAINT "SEVERITIES_PK" PRIMARY KEY ("ID")
  USING INDEX  ENABLE
   );


-- Create Variant Consequences
CREATE SEQUENCE VARIANTS_CONSEQUENCES_SEQ;


CREATE TABLE  VARIANTS_CONSEQUENCES
   (	"ID" NUMBER NOT NULL ENABLE, 
	"SEVERITY" NUMBER, 
	"VARIANT_TRANSCRIPT" NUMBER UNIQUE, 
	 CONSTRAINT "VARIANTS_CONSEQUENCES_PK" PRIMARY KEY ("ID")
  USING INDEX  ENABLE
   );

ALTER TABLE  "VARIANTS_CONSEQUENCES" ADD CONSTRAINT "VARIANTS_CONSEQUENCES_FK" FOREIGN KEY ("VARIANT_TRANSCRIPT")
	  REFERENCES  "VARIANTS_TRANSCRIPTS" ("ID") ON DELETE CASCADE ENABLE;

ALTER TABLE  "VARIANTS_CONSEQUENCES" ADD CONSTRAINT "VARIANTS_CONSEQUENCES_FK2" FOREIGN KEY ("SEVERITY")
	  REFERENCES  "SEVERITIES" ("ID") ON DELETE CASCADE ENABLE;


-- Create SV Consequences
CREATE SEQUENCE SV_CONSEQUENCES_SEQ;

CREATE TABLE  SV_CONSEQUENCES 
   (	"ID" NUMBER NOT NULL ENABLE, 
	"GENE" NUMBER, 
	"VARIANT" NUMBER UNIQUE NOT NULL ENABLE, 
	"CONSEQUENCE" VARCHAR2(100), 
	 CONSTRAINT "SV_CONSEQUENCES_PK" PRIMARY KEY ("ID")
  USING INDEX  ENABLE
   );

ALTER TABLE  "SV_CONSEQUENCES" ADD CONSTRAINT "SV_CONSEQUENCES_VAR_FK" FOREIGN KEY ("VARIANT")
	  REFERENCES  "VARIANTS" ("ID") ENABLE;

ALTER TABLE  "SV_CONSEQUENCES" ADD CONSTRAINT "SV_CONSEQUENCES_GENE_FK" FOREIGN KEY ("GENE")
	  REFERENCES  "GENES" ("ID") ON DELETE CASCADE ENABLE;






CREATE OR REPLACE EDITIONABLE TRIGGER "bi_GENES" 
  BEFORE INSERT ON "GENES"              
  FOR EACH ROW 
BEGIN  
  IF :new."ID" IS NULL THEN
    SELECT "GENES_SEQ".nextval INTO :new."ID" FROM sys.dual;
  END IF;
END;
/

ALTER TRIGGER "bi_GENES" ENABLE;


CREATE OR REPLACE EDITIONABLE TRIGGER "bi_VARIANTS" 
  BEFORE INSERT ON VARIANTS              
  FOR EACH ROW 
BEGIN  
  IF :new.ID IS NULL THEN
    SELECT VARIANTS_SEQ.nextval INTO :new.ID FROM sys.dual;
  END IF;
END;
/

ALTER TRIGGER "bi_VARIANTS" ENABLE;


CREATE OR REPLACE EDITIONABLE TRIGGER "bi_SNVS" 
  BEFORE INSERT ON "SNVS"              
  FOR EACH ROW 
BEGIN  
  IF :new."ID" IS NULL THEN
    SELECT "SNVS_SEQ".nextval INTO :new."ID" FROM sys.dual;
  END IF;
END;
/

ALTER TRIGGER "bi_SNVS" ENABLE;


CREATE OR REPLACE EDITIONABLE TRIGGER "bi_SVS" 
  BEFORE INSERT ON "SVS"              
  FOR EACH ROW 
BEGIN  
  IF :new."ID" IS NULL THEN
    SELECT "SVS_SEQ".nextval INTO :new."ID" FROM sys.dual;
  END IF;
END;
/

ALTER TRIGGER "bi_SVS" ENABLE;


CREATE OR REPLACE EDITIONABLE TRIGGER "bi_SVS_CTX" 
  BEFORE INSERT ON "SVS_CTX"              
  FOR EACH ROW 
BEGIN  
  IF :new."ID" IS NULL THEN
    SELECT "SVS_CTX_SEQ".nextval INTO :new."ID" FROM sys.dual;
  END IF;
END;
/

ALTER TRIGGER "bi_SVS_CTX" ENABLE;


CREATE OR REPLACE EDITIONABLE TRIGGER "bi_STR" 
  BEFORE INSERT ON "STR"              
  FOR EACH ROW 
BEGIN  
  IF :new."ID" IS NULL THEN
    SELECT "STR_SEQ".nextval INTO :new."ID" FROM sys.dual;
  END IF;
END;
/

ALTER TRIGGER "bi_STR" ENABLE;


CREATE OR REPLACE EDITIONABLE TRIGGER "bi_MTS" 
  BEFORE INSERT ON "MTS"              
  FOR EACH ROW 
BEGIN  
  IF :new."ID" IS NULL THEN
    SELECT "MTS_SEQ".nextval INTO :new."ID" FROM sys.dual;
  END IF;
END;
/

ALTER TRIGGER "bi_MTS" ENABLE;


CREATE OR REPLACE EDITIONABLE TRIGGER "bi_GENOMIC_GNOMAD_FREQUENCIES" 
  BEFORE INSERT ON "GENOMIC_GNOMAD_FREQUENCIES"              
  FOR EACH ROW 
BEGIN  
  IF :new."ID" IS NULL THEN
    SELECT "GENOMIC_GNOMAD_FREQUENCIES_SEQ".nextval INTO :new."ID" FROM sys.dual;
  END IF;
END;
/

ALTER TRIGGER "bi_GENOMIC_GNOMAD_FREQUENCIES" ENABLE;


CREATE OR REPLACE EDITIONABLE TRIGGER "bi_GENOMIC_IBVL_FREQUENCIES" 
  BEFORE INSERT ON "GENOMIC_IBVL_FREQUENCIES"              
  FOR EACH ROW 
BEGIN  
  IF :new."ID" IS NULL THEN
    SELECT "GENOMIC_IBVL_FREQUENCIES_SEQ".nextval INTO :new."ID" FROM sys.dual;
  END IF;
END;
/

ALTER TRIGGER "bi_GENOMIC_IBVL_FREQUENCIES" ENABLE;


CREATE OR REPLACE EDITIONABLE TRIGGER "bi_MT_GNOMAD_FREQUENCIES" 
  BEFORE INSERT ON "MT_GNOMAD_FREQUENCIES"              
  FOR EACH ROW 
BEGIN  
  IF :new."ID" IS NULL THEN
    SELECT "MT_GNOMAD_FREQUENCIES_SEQ".nextval INTO :new."ID" FROM sys.dual;
  END IF;
END;
/

ALTER TRIGGER "bi_MT_GNOMAD_FREQUENCIES" ENABLE;


CREATE OR REPLACE EDITIONABLE TRIGGER "bi_MT_IBVL_FREQUENCIES" 
  BEFORE INSERT ON "MT_IBVL_FREQUENCIES"              
  FOR EACH ROW 
BEGIN  
  IF :new."ID" IS NULL THEN
    SELECT "MT_IBVL_FREQUENCIES_SEQ".nextval INTO :new."ID" FROM sys.dual;
  END IF;
END;
/

ALTER TRIGGER "bi_MT_IBVL_FREQUENCIES" ENABLE;


CREATE OR REPLACE EDITIONABLE TRIGGER "bi_TRANSCRIPTS" 
  BEFORE INSERT ON "TRANSCRIPTS"              
  FOR EACH ROW 
BEGIN  
  IF :new."ID" IS NULL THEN
    SELECT "TRANSCRIPTS_SEQ".nextval INTO :new."ID" FROM sys.dual;
  END IF;
END;
/

ALTER TRIGGER "bi_TRANSCRIPTS" ENABLE;


CREATE OR REPLACE EDITIONABLE TRIGGER "bi_VARIANTS_TRANSCRIPTS" 
  BEFORE INSERT ON "VARIANTS_TRANSCRIPTS"              
  FOR EACH ROW 
BEGIN  
  IF :new."ID" IS NULL THEN
    SELECT "VARIANTS_TRANSCRIPTS_SEQ".nextval INTO :new."ID" FROM sys.dual;
  END IF;
END;
/

ALTER TRIGGER "bi_VARIANTS_TRANSCRIPTS" ENABLE;


CREATE OR REPLACE EDITIONABLE TRIGGER "bi_VARIANTS_ANNOTATIONS" 
  BEFORE INSERT ON "VARIANTS_ANNOTATIONS"              
  FOR EACH ROW 
BEGIN  
  IF :new."ID" IS NULL THEN
    SELECT "VARIANTS_ANNOTATIONS_SEQ".nextval INTO :new."ID" FROM sys.dual;
  END IF;
END;
/

ALTER TRIGGER "bi_VARIANTS_ANNOTATIONS" ENABLE;


CREATE OR REPLACE EDITIONABLE TRIGGER "bi_SEVERITIES" 
  BEFORE INSERT ON "SEVERITIES"              
  FOR EACH ROW 
BEGIN  
  IF :new."ID" IS NULL THEN
    SELECT "SEVERITIES_SEQ".nextval INTO :new."ID" FROM sys.dual;
  END IF;
END;
/

ALTER TRIGGER "bi_SEVERITIES" ENABLE;


CREATE OR REPLACE EDITIONABLE TRIGGER "bi_VARIANTS_CONSEQUENCES" 
  BEFORE INSERT ON "VARIANTS_CONSEQUENCES"              
  FOR EACH ROW 
BEGIN  
  IF :new."ID" IS NULL THEN
    SELECT "VARIANTS_CONSEQUENCES_SEQ".nextval INTO :new."ID" FROM sys.dual;
  END IF;
END;
/

ALTER TRIGGER "bi_VARIANTS_CONSEQUENCES" ENABLE;


CREATE OR REPLACE EDITIONABLE TRIGGER "bi_SV_CONSEQUENCES" 
  BEFORE INSERT ON "SV_CONSEQUENCES"              
  FOR EACH ROW 
BEGIN  
  IF :new."ID" IS NULL THEN
    SELECT "SV_CONSEQUENCES_SEQ".nextval INTO :new."ID" FROM sys.dual;
  END IF;
END;
/

ALTER TRIGGER "bi_SV_CONSEQUENCES" ENABLE;

