-- Create Genes
CREATE SEQUENCE GENES_SEQ  
    MINVALUE 1 
    MAXVALUE 9999999999999999999999999999 
    INCREMENT BY 1 
    START WITH 1 
    CACHE 20 
    NOORDER  
    NOCYCLE  
    NOKEEP  
    GLOBAL

CREATE TABLE  "GENES" 
   (	"ID" NUMBER NOT NULL ENABLE,  
	"SHORT_NAME" VARCHAR2(30) NOT NULL ENABLE UNIQUE, 
	 CONSTRAINT "GENES_PK" PRIMARY KEY ("ID")
  USING INDEX  ENABLE
   )
/

CREATE OR REPLACE EDITIONABLE TRIGGER  "bi_GENES" 
  before insert on "GENES"              
  for each row 
begin  
  if :new."ID" is null then
    select "GENES_SEQ".nextval into :new."ID" from sys.dual;
  end if;
end;

/
ALTER TRIGGER  "bi_GENES" ENABLE
/

-- Create Variants
CREATE SEQUENCE VARIANTS_SEQ  
    MINVALUE 1 
    MAXVALUE 9999999999999999999999999999 
    INCREMENT BY 1 
    START WITH 1 
    CACHE 20 
    NOORDER  
    NOCYCLE  
    NOKEEP  
    GLOBAL;
/

CREATE TABLE VARIANTS 
   (	"VARIANT_ID" VARCHAR2(255) UNIQUE, 
        "ID" NUMBER NOT NULL ENABLE,
	    "VAR_TYPE" VARCHAR2(30), 
	    CONSTRAINT "VARIANTS_PK" PRIMARY KEY ("ID") USING INDEX ENABLE
   )
/

CREATE OR REPLACE EDITIONABLE TRIGGER  "bi_VARIANTS" 
    before insert on VARIANTS              
        for each row 
            begin  
                if :new.ID is null then
                    select VARIANTS_SEQ.nextval into :new.ID from sys.dual;
                end if;
            end;

/
ALTER TRIGGER  "bi_VARIANTS" ENABLE
/

-- Create SNVs
CREATE SEQUENCE SNVS_SEQ  
    MINVALUE 1 
    MAXVALUE 9999999999999999999999999999 
    INCREMENT BY 1 
    START WITH 1 
    CACHE 20 
    NOORDER  
    NOCYCLE  
    NOKEEP  
    GLOBAL
/

CREATE TABLE  "SNVS" 
   (	"VARIANT" NUMBER UNIQUE, 
    "ID" NUMBER NOT NULL ENABLE,
	"TYPE" VARCHAR2(30), 
	"LENGTH" NUMBER, 
	"CHR" NUMBER, 
	"POS" NUMBER, 
	"REF" VARCHAR2(255), 
	"ALT" VARCHAR2(255), 
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
   )
/
ALTER TABLE  "SNVS" ADD CONSTRAINT "SNVS_FK" FOREIGN KEY ("VARIANT")
	  REFERENCES  "VARIANTS" ("ID") ON DELETE CASCADE ENABLE
/

CREATE OR REPLACE EDITIONABLE TRIGGER  "bi_SNVS" 
  before insert on "SNVS"              
  for each row 
begin  
  if :new."ID" is null then
    select "SNVS_SEQ".nextval into :new."ID" from sys.dual;
  end if;
end;

/
ALTER TRIGGER  "bi_SNVS" ENABLE
/

-- Create SVs
CREATE SEQUENCE SVS_SEQ  
    MINVALUE 1 
    MAXVALUE 9999999999999999999999999999 
    INCREMENT BY 1 
    START WITH 1 
    CACHE 20 
    NOORDER  
    NOCYCLE  
    NOKEEP  
    GLOBAL
/

CREATE TABLE  "SVS" 
   (	"VARIANT" NUMBER UNIQUE,
    "ID" NUMBER NOT NULL ENABLE, 
	"CHR1" NUMBER, 
	"CHR1_POS1" NUMBER, 
	"CHR1_POS2" NUMBER, 
	"TYPE" VARCHAR2(30), 
	"SV_LENGTH" NUMBER, 
	"ALGORITHM" VARCHAR2(30), 
	"UCSC_URL" VARCHAR2(511), 
	"GNOMAD_ID" VARCHAR2(30), 
	"GNOMAD_URL" VARCHAR2(511), 
	 CONSTRAINT "SVS_PK" PRIMARY KEY ("ID")
  USING INDEX  ENABLE
   )
/
ALTER TABLE  "SVS" ADD CONSTRAINT "SVS_FK" FOREIGN KEY ("VARIANT")
	  REFERENCES  "VARIANTS" ("ID") ON DELETE CASCADE ENABLE
/

CREATE OR REPLACE EDITIONABLE TRIGGER  "bi_SVS" 
  before insert on "SVS"              
  for each row 
begin  
  if :new."ID" is null then
    select "SVS_SEQ".nextval into :new."ID" from sys.dual;
  end if;
end;

/
ALTER TRIGGER  "bi_SVS" ENABLE
/

-- Create SVs CTX
CREATE SEQUENCE SVS_CTX_SEQ  
    MINVALUE 1 
    MAXVALUE 9999999999999999999999999999 
    INCREMENT BY 1 
    START WITH 1 
    CACHE 20 
    NOORDER  
    NOCYCLE  
    NOKEEP  
    GLOBAL
/

CREATE TABLE  "SVS_CTX" 
   (	"VARIANT" NUMBER UNIQUE,
    "ID" NUMBER NOT NULL ENABLE, 
	"CHR2" NUMBER, 
	"CHR2_POS1" NUMBER, 
	"UCSC_URL2" VARCHAR2(511), 
    "GNOMAD_ID2" VARCHAR2(30),
    "GNOMAD_URL2" VARCHAR2(511),
	 CONSTRAINT "SVS_CTX_PK" PRIMARY KEY ("ID")
  USING INDEX  ENABLE
   )
/
ALTER TABLE  "SVS_CTX" ADD CONSTRAINT "SVS_CTX_FK" FOREIGN KEY ("VARIANT")
	  REFERENCES  "VARIANTS" ("ID") ON DELETE CASCADE ENABLE
/

CREATE OR REPLACE EDITIONABLE TRIGGER  "bi_SVS_CTX" 
  before insert on "SVS_CTX"              
  for each row 
begin  
  if :new."ID" is null then
    select "SVS_CTX_SEQ".nextval into :new."ID" from sys.dual;
  end if;
end;

/
ALTER TRIGGER  "bi_SVS_CTX" ENABLE
/

-- Create STR
CREATE SEQUENCE STR_SEQ  
    MINVALUE 1 
    MAXVALUE 9999999999999999999999999999 
    INCREMENT BY 1 
    START WITH 1 
    CACHE 20 
    NOORDER  
    NOCYCLE  
    NOKEEP  
    GLOBAL
/

CREATE TABLE  "STR" 
   (	"VARIANT" NUMBER UNIQUE,
    "ID" NUMBER NOT NULL ENABLE, 
	"REPEAT_UNIT" NUMBER, 
	"MIN_N_REPEAT" NUMBER, 
	"MAX_N_REPEAT" NUMBER, 
    "ALLELE_DISTRIB" VARCHAR2(255),
	 CONSTRAINT "STR_PK" PRIMARY KEY ("ID")
  USING INDEX  ENABLE
   )
/
ALTER TABLE  "STR" ADD CONSTRAINT "STR_FK" FOREIGN KEY ("VARIANT")
	  REFERENCES  "VARIANTS" ("ID") ON DELETE CASCADE ENABLE
/

CREATE OR REPLACE EDITIONABLE TRIGGER  "bi_STR" 
  before insert on "STR"              
  for each row 
begin  
  if :new."ID" is null then
    select "STR_SEQ".nextval into :new."ID" from sys.dual;
  end if;
end;

/
ALTER TRIGGER  "bi_STR" ENABLE
/

-- Create Mts
CREATE SEQUENCE MTS_SEQ  
    MINVALUE 1 
    MAXVALUE 9999999999999999999999999999 
    INCREMENT BY 1 
    START WITH 1 
    CACHE 20 
    NOORDER  
    NOCYCLE  
    NOKEEP  
    GLOBAL
/

CREATE TABLE  "MTS" 
   (	"VARIANT" NUMBER UNIQUE, 
    "ID" NUMBER NOT NULL ENABLE,
	"POS" NUMBER, 
	"REF" VARCHAR2(1), 
	"ALT" VARCHAR2(30), 
	"UCSC_URL" VARCHAR2(511), 
	"MITOMAP_URL" VARCHAR2(511), 
	"GNOMAD_URL" VARCHAR2(511), 
	"DBSNP_ID" VARCHAR2(30), 
	"DBSNP_URL" VARCHAR2(511), 
	"CLINVAR_URL" VARCHAR2(511), 
    "CLINVAR_VCV" NUMBER,
	 CONSTRAINT "MTS_PK" PRIMARY KEY ("ID")
  USING INDEX  ENABLE
   )
/
ALTER TABLE  "MTS" ADD CONSTRAINT "MTS_FK" FOREIGN KEY ("VARIANT")
	  REFERENCES  "VARIANTS" ("ID") ON DELETE CASCADE ENABLE
/

CREATE OR REPLACE EDITIONABLE TRIGGER  "bi_MTS" 
  before insert on "MTS"              
  for each row 
begin  
  if :new."ID" is null then
    select "MTS_SEQ".nextval into :new."ID" from sys.dual;
  end if;
end;

/
ALTER TRIGGER  "bi_MTS" ENABLE
/

-- Create Genomic gnomAD Frequencies
CREATE SEQUENCE GENOMIC_GNOMAD_FREQUENCIES_SEQ  
    MINVALUE 1 
    MAXVALUE 9999999999999999999999999999 
    INCREMENT BY 1 
    START WITH 1 
    CACHE 20 
    NOORDER  
    NOCYCLE  
    NOKEEP  
    GLOBAL
/

CREATE TABLE  "GENOMIC_GNOMAD_FREQUENCIES" 
   (	"ID" NUMBER NOT NULL ENABLE, 
	"VARIANT" NUMBER UNIQUE, 
	"AF_TOTAL" NUMBER, 
	"AC" NUMBER, 
	"AN" NUMBER, 
	"HOM_ALT" NUMBER, 
	 CONSTRAINT "GENOMIC_GNOMAD_FREQUENCIES_PK" PRIMARY KEY ("ID")
  USING INDEX  ENABLE
   )
/
ALTER TABLE  "GENOMIC_GNOMAD_FREQUENCIES" ADD CONSTRAINT "GENOMIC_GNOMAD_FREQUENCIES_FK" FOREIGN KEY ("VARIANT")
	  REFERENCES  "VARIANTS" ("ID") ON DELETE CASCADE ENABLE
/

CREATE OR REPLACE EDITIONABLE TRIGGER  "bi_GENOMIC_GNOMAD_FREQUENCIES" 
  before insert on "GENOMIC_GNOMAD_FREQUENCIES"              
  for each row 
begin  
  if :new."ID" is null then
    select "GENOMIC_GNOMAD_FREQUENCIES_SEQ".nextval into :new."ID" from sys.dual;
  end if;
end;

/
ALTER TRIGGER  "bi_GENOMIC_GNOMAD_FREQUENCIES" ENABLE
/

-- Create Genomic IBVL Frequencies
CREATE SEQUENCE GENOMIC_IBVL_FREQUENCIES_SEQ  
    MINVALUE 1 
    MAXVALUE 9999999999999999999999999999 
    INCREMENT BY 1 
    START WITH 1 
    CACHE 20 
    NOORDER  
    NOCYCLE  
    NOKEEP  
    GLOBAL
/

CREATE TABLE  "GENOMIC_IBVL_FREQUENCIES" 
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
   )
/
ALTER TABLE  "GENOMIC_IBVL_FREQUENCIES" ADD CONSTRAINT "GENOMIC_IBVL_FREQUENCIES_FK" FOREIGN KEY ("VARIANT")
	  REFERENCES  "VARIANTS" ("ID") ON DELETE CASCADE ENABLE
/

CREATE OR REPLACE EDITIONABLE TRIGGER  "bi_GENOMIC_IBVL_FREQUENCIES" 
  before insert on "GENOMIC_IBVL_FREQUENCIES"              
  for each row 
begin  
  if :new."ID" is null then
    select "GENOMIC_IBVL_FREQUENCIES_SEQ".nextval into :new."ID" from sys.dual;
  end if;
end;

/
ALTER TRIGGER  "bi_GENOMIC_IBVL_FREQUENCIES" ENABLE
/

-- Create Mt gnomAD Frequencies
CREATE SEQUENCE MT_GNOMAD_FREQUENCIES_SEQ  
    MINVALUE 1 
    MAXVALUE 9999999999999999999999999999 
    INCREMENT BY 1 
    START WITH 1 
    CACHE 20 
    NOORDER  
    NOCYCLE  
    NOKEEP  
    GLOBAL
/

CREATE TABLE  "MT_GNOMAD_FREQUENCIES" 
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
   )
/
ALTER TABLE  "MT_GNOMAD_FREQUENCIES" ADD CONSTRAINT "MT_GNOMAD_FREQUENCIES_FK" FOREIGN KEY ("VARIANT")
	  REFERENCES  "VARIANTS" ("ID") ON DELETE CASCADE ENABLE
/

CREATE OR REPLACE EDITIONABLE TRIGGER  "bi_MT_GNOMAD_FREQUENCIES" 
  before insert on "MT_GNOMAD_FREQUENCIES"              
  for each row 
begin  
  if :new."ID" is null then
    select "MT_GNOMAD_FREQUENCIES_SEQ".nextval into :new."ID" from sys.dual;
  end if;
end;

/
ALTER TRIGGER  "bi_MT_GNOMAD_FREQUENCIES" ENABLE
/

-- Create Mt IBVL Frequencies
CREATE SEQUENCE MT_IBVL_FREQUENCIES_SEQ  
    MINVALUE 1 
    MAXVALUE 9999999999999999999999999999 
    INCREMENT BY 1 
    START WITH 1 
    CACHE 20 
    NOORDER  
    NOCYCLE  
    NOKEEP  
    GLOBAL
/

CREATE TABLE  "MT_IBVL_FREQUENCIES" 
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
   )
/
ALTER TABLE  "MT_IBVL_FREQUENCIES" ADD CONSTRAINT "MT_IBVL_FREQUENCIES_FK" FOREIGN KEY ("VARIANT")
	  REFERENCES  "VARIANTS" ("ID") ON DELETE CASCADE ENABLE
/

CREATE OR REPLACE EDITIONABLE TRIGGER  "bi_MT_IBVL_FREQUENCIES" 
  before insert on "MT_IBVL_FREQUENCIES"              
  for each row 
begin  
  if :new."ID" is null then
    select "MT_IBVL_FREQUENCIES_SEQ".nextval into :new."ID" from sys.dual;
  end if;
end;

/
ALTER TRIGGER  "bi_MT_IBVL_FREQUENCIES" ENABLE
/

-- Create Transcripts
CREATE SEQUENCE TRANSCRIPTS_SEQ  
    MINVALUE 1 
    MAXVALUE 9999999999999999999999999999 
    INCREMENT BY 1 
    START WITH 1 
    CACHE 20 
    NOORDER  
    NOCYCLE  
    NOKEEP  
    GLOBAL
/

CREATE TABLE  "TRANSCRIPTS" 
   (	"TRANSCRIPT_ID" VARCHAR2(255) UNIQUE,
    "ID" NUMBER NOT NULL ENABLE, 
	"GENE" NUMBER, 
	"TRANSCRIPT_TYPE" VARCHAR2(1), 
    "TSL" VARCHAR2(255),
	 CONSTRAINT "TRANSCRIPTS_PK" PRIMARY KEY ("ID")
  USING INDEX  ENABLE
   )
/
ALTER TABLE  "TRANSCRIPTS" ADD CONSTRAINT "TRANSCRIPTS_GENES_FK" FOREIGN KEY ("GENE")
	  REFERENCES  "GENES" ("ID") ON DELETE CASCADE ENABLE
/

CREATE OR REPLACE EDITIONABLE TRIGGER  "bi_TRANSCRIPTS" 
  before insert on "TRANSCRIPTS"              
  for each row 
begin  
  if :new."ID" is null then
    select "TRANSCRIPTS_SEQ".nextval into :new."ID" from sys.dual;
  end if;
end;

/
ALTER TRIGGER  "bi_TRANSCRIPTS" ENABLE
/

-- Create Variants Transcripts
CREATE SEQUENCE VARIANTS_TRANSCRIPTS_SEQ  
    MINVALUE 1 
    MAXVALUE 9999999999999999999999999999 
    INCREMENT BY 1 
    START WITH 1 
    CACHE 20 
    NOORDER  
    NOCYCLE  
    NOKEEP  
    GLOBAL
/

CREATE TABLE  "VARIANTS_TRANSCRIPTS" 
   (	"TRANSCRIPT" NUMBER, 
    "ID" NUMBER NOT NULL ENABLE,
	"VARIANT" NUMBER, 
	"HGVSC" VARCHAR2(255), 
	 CONSTRAINT "VARIANTS_TRANSCRIPTS_PK" PRIMARY KEY ("ID")
  USING INDEX  ENABLE
   )
/
ALTER TABLE  "VARIANTS_TRANSCRIPTS" ADD CONSTRAINT "VARIANTS_TRANSCRIPTS_FK" FOREIGN KEY ("TRANSCRIPT")
	  REFERENCES  "TRANSCRIPTS" ("ID") ON DELETE CASCADE ENABLE
/
ALTER TABLE  "VARIANTS_TRANSCRIPTS" ADD CONSTRAINT "VARIANTS_TRANSCRIPTS_FK2" FOREIGN KEY ("VARIANT")
	  REFERENCES  "VARIANTS" ("ID") ON DELETE CASCADE ENABLE
/
ALTER TABLE "VARIANTS_TRANSCRIPTS" ADD CONSTRAINT "VARIANTS_TRANSCRIPTS_UNIQUE" UNIQUE("TRANSCRIPT", "VARIANT")
/
CREATE OR REPLACE EDITIONABLE TRIGGER  "bi_VARIANTS_TRANSCRIPTS" 
  before insert on "VARIANTS_TRANSCRIPTS"              
  for each row 
begin  
  if :new."ID" is null then
    select "VARIANTS_TRANSCRIPTS_SEQ".nextval into :new."ID" from sys.dual;
  end if;
end;

/
ALTER TRIGGER  "bi_VARIANTS_TRANSCRIPTS" ENABLE
/

-- Create Variants Annotations
CREATE SEQUENCE VARIANTS_ANNOTATIONS_SEQ  
    MINVALUE 1 
    MAXVALUE 9999999999999999999999999999 
    INCREMENT BY 1 
    START WITH 1 
    CACHE 20 
    NOORDER  
    NOCYCLE  
    NOKEEP  
    GLOBAL
/

CREATE TABLE  "VARIANTS_ANNOTATIONS" 
   (	"ID" NUMBER NOT NULL ENABLE, 
	"HGVSP" VARCHAR2(255), 
	"POLYPHEN" VARCHAR2(255), 
	"SIFT" VARCHAR2(255), 
	"VARIANT_TRANSCRIPT" NUMBER UNIQUE, 
	 CONSTRAINT "VARIANTS_ANNOTATIONS_PK" PRIMARY KEY ("ID")
  USING INDEX  ENABLE
   )
/
ALTER TABLE  "VARIANTS_ANNOTATIONS" ADD CONSTRAINT "VARIANTS_ANNOTATIONS_FK" FOREIGN KEY ("VARIANT_TRANSCRIPT")
	  REFERENCES  "VARIANTS_TRANSCRIPTS" ("ID") ON DELETE CASCADE ENABLE
/

CREATE OR REPLACE EDITIONABLE TRIGGER  "bi_VARIANTS_ANNOTATIONS" 
  before insert on "VARIANTS_ANNOTATIONS"              
  for each row 
begin  
  if :new."ID" is null then
    select "VARIANTS_ANNOTATIONS_SEQ".nextval into :new."ID" from sys.dual;
  end if;
end;

/
ALTER TRIGGER  "bi_VARIANTS_ANNOTATIONS" ENABLE
/

-- Create Severities
CREATE SEQUENCE SEVERITIES_SEQ  
    MINVALUE 1 
    MAXVALUE 9999999999999999999999999999 
    INCREMENT BY 1 
    START WITH 1 
    CACHE 20 
    NOORDER  
    NOCYCLE  
    NOKEEP  
    GLOBAL
/

CREATE TABLE  "SEVERITIES" 
   (	"SEVERITY_NUMBER" NUMBER UNIQUE, 
    "ID" NUMBER NOT NULL ENABLE,
	"CONSEQUENCE" VARCHAR2(255), 
	 CONSTRAINT "SEVERITIES_PK" PRIMARY KEY ("ID")
  USING INDEX  ENABLE
   )
/

CREATE OR REPLACE EDITIONABLE TRIGGER  "bi_SEVERITIES" 
  before insert on "SEVERITIES"              
  for each row 
begin  
  if :new."ID" is null then
    select "SEVERITIES_SEQ".nextval into :new."ID" from sys.dual;
  end if;
end;

/
ALTER TRIGGER  "bi_SEVERITIES" ENABLE
/

-- Create Variant Consequences
CREATE SEQUENCE VARIANTS_CONSEQUENCES_SEQ  
    MINVALUE 1 
    MAXVALUE 9999999999999999999999999999 
    INCREMENT BY 1 
    START WITH 1 
    CACHE 20 
    NOORDER  
    NOCYCLE  
    NOKEEP  
    GLOBAL
/

CREATE TABLE  "VARIANTS_CONSEQUENCES" 
   (	"ID" NUMBER NOT NULL ENABLE, 
	"SEVERITY" NUMBER, 
	"VARIANT_TRANSCRIPT" NUMBER UNIQUE, 
	 CONSTRAINT "VARIANTS_CONSEQUENCES_PK" PRIMARY KEY ("ID")
  USING INDEX  ENABLE
   )
/
ALTER TABLE  "VARIANTS_CONSEQUENCES" ADD CONSTRAINT "VARIANTS_CONSEQUENCES_FK" FOREIGN KEY ("VARIANT_TRANSCRIPT")
	  REFERENCES  "VARIANTS_TRANSCRIPTS" ("ID") ON DELETE CASCADE ENABLE
/
ALTER TABLE  "VARIANTS_CONSEQUENCES" ADD CONSTRAINT "VARIANTS_CONSEQUENCES_FK2" FOREIGN KEY ("SEVERITY")
	  REFERENCES  "SEVERITIES" ("ID") ON DELETE CASCADE ENABLE
/

CREATE OR REPLACE EDITIONABLE TRIGGER  "bi_VARIANTS_CONSEQUENCES" 
  before insert on "VARIANTS_CONSEQUENCES"              
  for each row 
begin  
  if :new."ID" is null then
    select "VARIANTS_CONSEQUENCES_SEQ".nextval into :new."ID" from sys.dual;
  end if;
end;

/
ALTER TRIGGER  "bi_VARIANTS_CONSEQUENCES" ENABLE
/

-- Create SV Consequences
CREATE SEQUENCE SV_CONSEQUENCES_SEQ  
    MINVALUE 1 
    MAXVALUE 9999999999999999999999999999 
    INCREMENT BY 1 
    START WITH 1 
    CACHE 20 
    NOORDER  
    NOCYCLE  
    NOKEEP  
    GLOBAL
/

CREATE TABLE  "SV_CONSEQUENCES" 
   (	"ID" NUMBER NOT NULL ENABLE, 
	"GENE" NUMBER, 
	"VARIANT" NUMBER UNIQUE NOT NULL ENABLE, 
	"CONSEQUENCE" VARCHAR2(255), 
	 CONSTRAINT "SV_CONSEQUENCES_PK" PRIMARY KEY ("ID")
  USING INDEX  ENABLE
   )
/
ALTER TABLE  "SV_CONSEQUENCES" ADD CONSTRAINT "SV_CONSEQUENCES_VAR_FK" FOREIGN KEY ("VARIANT")
	  REFERENCES  "VARIANTS" ("ID") ENABLE
/
ALTER TABLE  "SV_CONSEQUENCES" ADD CONSTRAINT "SV_CONSEQUENCES_GENE_FK" FOREIGN KEY ("GENE")
	  REFERENCES  "GENES" ("ID") ON DELETE CASCADE ENABLE
/

CREATE OR REPLACE EDITIONABLE TRIGGER  "bi_SV_CONSEQUENCES" 
  before insert on "SV_CONSEQUENCES"              
  for each row 
begin  
  if :new."ID" is null then
    select "SV_CONSEQUENCES_SEQ".nextval into :new."ID" from sys.dual;
  end if;
end;

/
ALTER TRIGGER  "bi_SV_CONSEQUENCES" ENABLE
/

-- TODO: Update creation of Users here
CREATE SEQUENCE   TEST_USERS_SEQ  
    MINVALUE 1 
    MAXVALUE 999999999999999999999999999 
    INCREMENT BY 1 
    START WITH 1 
    NOCACHE  
    NOORDER  
    NOCYCLE  
    NOKEEP  
    GLOBAL

CREATE TABLE  "TEST_USERS" 
   (	"USER_ID" NUMBER, 
	"USER_NAME" VARCHAR2(255) NOT NULL ENABLE, 
	"PASSWORD" VARCHAR2(255) NOT NULL ENABLE, 
	 PRIMARY KEY ("USER_ID")
  USING INDEX  ENABLE, 
	 CONSTRAINT "USERS_U1" UNIQUE ("USER_NAME")
  USING INDEX  ENABLE
   )
/

CREATE OR REPLACE EDITIONABLE TRIGGER  "TEST_BI_USERS" 
  before insert on test_users 
      for each row 
          begin 
            -- Get a unique sequence value to use as the primary key
            SELECT "TEST_USERS_SEQ".nextval INTO :new.user_id from sys.dual; 
            -- Make sure to save the username in upper case
            :new.user_name := upper(:new.user_name); 
            -- Hash the password so we are not saving clear text
            :new.password := hash_password(upper(:new.user_name), :new.password); 
            --DBMS_OUTPUT.PUT_LINE('HI');
          end;

/
ALTER TRIGGER  "TEST_BI_USERS" ENABLE
/

CREATE OR REPLACE EDITIONABLE TRIGGER  "TEST_BU_USERS" 
  before update on test_users 
  for each row 
begin 
  -- Make sure to save the user name in upper case
  :new.user_name := upper(:new.user_name); 
  -- If the new password is not null 
  if :new.password is not null then 
    -- Make sure to hash the password so it is not stored in clear text
    :new.password := hash_password(upper(:new.user_name), :new.password); 
  -- If the password is empty
  else 
    -- Keep the old hashed password. We don't want a blank password.
    :new.password := :old.password; 
  end if; 
end;

/
ALTER TRIGGER  "TEST_BU_USERS" ENABLE
/
