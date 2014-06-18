/*
 Navicat SQLite Data Transfer

 Source Server         : RawDB
 Source Server Version : 3008004
 Source Database       : main

 Target Server Version : 3008004
 File Encoding         : utf-8

 Date: 06/10/2014 14:45:40 PM
*/

PRAGMA foreign_keys = false;

-- ----------------------------
--  Table structure for HMM_Data
-- ----------------------------
DROP TABLE IF EXISTS "HMM_Data";
CREATE TABLE "HMM_Data" (
	 "HMM_Model" text(35,0) NOT NULL,
	 "HMM_Family" TEXT(35,0),
	PRIMARY KEY("HMM_Model")
);

-- ----------------------------
--  Table structure for HMM_Hits
-- ----------------------------
DROP TABLE IF EXISTS "HMM_Hits";
CREATE TABLE "HMM_Hits" (
	 "Hit_HASH" text(40,0) NOT NULL,
	 "Protein_Accession" text(25,0) NOT NULL,
	 "HMM_Model" text(35,0) NOT NULL,
	 "HMM_Score" integer(4,0) NOT NULL,
	 "HMM_E_Value" real(4,10) NOT NULL,
	 "Ali_From" integer(4,0) NOT NULL,
	 "Ali_To" integer(4,0) NOT NULL,
	 "HMM_From" integer(4,0) NOT NULL,
	 "HMM_To" integer(4,0) NOT NULL,
	 "HMM_Coverage" real(2,10) NOT NULL,
	PRIMARY KEY("Hit_HASH"),
	CONSTRAINT "fk_Proteins" FOREIGN KEY ("Protein_Accession") REFERENCES "Proteins" ("Protein_Accession") ON DELETE CASCADE ON UPDATE CASCADE,
	CONSTRAINT "fk_HMM_Data" FOREIGN KEY ("HMM_Model") REFERENCES "HMM_Data" ("HMM_Model") ON DELETE CASCADE ON UPDATE CASCADE
);

-- ----------------------------
--  Table structure for Organisms
-- ----------------------------
DROP TABLE IF EXISTS "Organisms";
CREATE TABLE "Organisms" (
	 "Organism_Accession" text(25,0) NOT NULL,
	 "Accession_Type" text(150,0) NOT NULL,
	 "Organism_Description" text(300,0) NOT NULL,
	 "Source" text(300,0) NOT NULL,
	 "Organism_Phylogeny" text(300,0) NOT NULL,
	 "Sequence_Length" integer(13631488,0),
	PRIMARY KEY("Organism_Accession")
);

-- ----------------------------
--  Table structure for Proteins
-- ----------------------------
DROP TABLE IF EXISTS "Proteins";
CREATE TABLE "Proteins" (
	 "Protein_Accession" text(25,0) NOT NULL,
	 "Organism_Accession" text(25,0) NOT NULL,
	 "Locus" text(25,0) NOT NULL,
	 "Start" integer(300,0) NOT NULL,
	 "End" integer(300,0) NOT NULL,
	 "Strand" integer(1,0) NOT NULL,
	 "FASTA_Sequence" text(5500,0) NOT NULL,
	PRIMARY KEY("Protein_Accession"),
	CONSTRAINT "fk_Organisms" FOREIGN KEY ("Organism_Accession") REFERENCES "Organisms" ("Organism_Accession") ON DELETE CASCADE ON UPDATE CASCADE
);

PRAGMA foreign_keys = true;
