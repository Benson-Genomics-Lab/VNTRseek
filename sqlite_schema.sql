-- Pragmas used when creating tables in the database only
PRAGMA synchronous = OFF;
PRAGMA journal_mode = MEMORY;
PRAGMA user_version = 1;
BEGIN TRANSACTION;
CREATE TABLE IF NOT EXISTS `clusterlnk` (
  `clusterid` integer NOT NULL,
  `repeatid` integer NOT NULL,
  `direction` char(1) NOT NULL,
  `reserved` integer NOT NULL,
  `reserved2` integer NOT NULL,
  PRIMARY KEY (`clusterid`,`repeatid`),
  UNIQUE (`repeatid`)
);
CREATE TABLE IF NOT EXISTS `clusters` (
  `cid` integer NOT NULL,
  `minpat` integer NOT NULL,
  `maxpat` integer NOT NULL,
  `refs_flank_undist` integer DEFAULT NULL,
  `refs_flank_dist` integer DEFAULT NULL,
  `refs_flank_undist_l` integer DEFAULT NULL,
  `refs_flank_undist_r` integer DEFAULT NULL,
  `refs_flank_undist_lr` integer DEFAULT NULL,
  `repeatcount` integer NOT NULL,
  `refcount` integer NOT NULL,
  `variability` integer NOT NULL DEFAULT '0',
  `assemblyreq` integer DEFAULT NULL,
  `profdensity` float DEFAULT NULL,
  `flankdensity` float DEFAULT NULL,
  `mcpattern` varchar(5000) DEFAULT NULL,
  `aveentropy` float DEFAULT NULL,
  PRIMARY KEY (`cid`)
);
CREATE TABLE IF NOT EXISTS `fasta_reads` (
  `sid` integer  NOT NULL,
  `head` varchar(500) NOT NULL,
  `dna` varchar(8000) DEFAULT NULL,
  `qual` varchar(8000) DEFAULT NULL,
  PRIMARY KEY (`sid`),
  UNIQUE (`head`)
);
-- sequence data moved to external db
CREATE TABLE IF NOT EXISTS `fasta_ref_reps` (
  `rid` integer NOT NULL,
  `comment` varchar(500) DEFAULT NULL,
  `flank_disting` integer DEFAULT NULL,
  `entropy` float NOT NULL DEFAULT '0',
  `has_support` integer DEFAULT NULL,
  `span1` integer DEFAULT NULL,
  `spanN` integer DEFAULT NULL,
  `homez_same` integer DEFAULT NULL,
  `homez_diff` integer DEFAULT NULL,
  `hetez_same` integer DEFAULT NULL,
  `hetez_diff` integer DEFAULT NULL,
  `hetez_multi` integer DEFAULT NULL,
  `support_vntr` integer DEFAULT NULL,
  `support_vntr_span1` integer DEFAULT NULL,
  `alleles_sup` integer DEFAULT NULL,
  `allele_sup_same_as_ref` integer DEFAULT NULL,
  PRIMARY KEY (`rid`),
  UNIQUE (`rid`,`comment`)
);
CREATE TABLE IF NOT EXISTS `flank_connection` (
  `refid` integer NOT NULL,
  `clusterid` integer NOT NULL,
  `fcomponentsize` integer NOT NULL,
  `fcomponentid` integer NOT NULL,
  `paramsetid` integer NOT NULL,
  PRIMARY KEY (`refid`,`paramsetid`),
  CONSTRAINT `flank_connection_ibfk_1` FOREIGN KEY (`paramsetid`) REFERENCES `flank_params` (`paramsetid`)
);
CREATE TABLE IF NOT EXISTS `flank_params` (
  `paramsetid` integer NOT NULL PRIMARY KEY,
  `flength` integer NOT NULL,
  `ferrors` integer NOT NULL
);
CREATE TABLE IF NOT EXISTS `map` (
  `refid` integer NOT NULL,
  `readid` integer NOT NULL,
  `reserved` integer NOT NULL,
  `reserved2` integer NOT NULL,
  `bbb` integer NOT NULL DEFAULT '0',
  PRIMARY KEY (`refid`,`readid`)
);
CREATE TABLE IF NOT EXISTS `rank` (
  `refid` integer NOT NULL,
  `readid` integer NOT NULL,
  `score` float DEFAULT NULL,
  `ties` integer NOT NULL DEFAULT '0',
  `refdir` char(1) NOT NULL,
  PRIMARY KEY (`refid`,`readid`)
);
CREATE TABLE IF NOT EXISTS `rankflank` (
  `refid` integer NOT NULL,
  `readid` integer NOT NULL,
  `score` float DEFAULT NULL,
  `ties` integer NOT NULL DEFAULT '0',
  PRIMARY KEY (`refid`,`readid`)
);
CREATE TABLE IF NOT EXISTS `replnk` (
  `rid` integer NOT NULL,
  `sid` integer  NOT NULL,
  `first` integer NOT NULL,
  `last` integer NOT NULL,
  `patsize` integer NOT NULL,
  `copynum` float NOT NULL,
  `pattern` text NOT NULL,
  `profile` text COLLATE BINARY,
  `profilerc` text COLLATE BINARY,
  `profsize` integer DEFAULT NULL,
  PRIMARY KEY (`rid`)
);
CREATE TABLE IF NOT EXISTS `stats` (
  `id` integer NOT NULL,
  `MAP_ROOT` varchar(500) DEFAULT NULL,
  `PARAM_TRF` varchar(500) DEFAULT NULL,
  `PARAM_PROCLU` varchar(500) DEFAULT NULL,
  `FOLDER_FASTA` varchar(500) DEFAULT NULL,
  `FOLDER_PROFILES` varchar(500) DEFAULT NULL,
  `FOLDER_PROFILES_CLEAN` varchar(500) DEFAULT NULL,
  `FOLDER_REFERENCE` varchar(500) DEFAULT NULL,
  `FILE_REFERENCE_LEB` varchar(500) DEFAULT NULL,
  `FILE_REFERENCE_SEQ` varchar(500) DEFAULT NULL,
  `NUMBER_READS` integer DEFAULT NULL,
  `NUMBER_TRS_IN_READS` integer DEFAULT NULL,
  `NUMBER_TRS_IN_READS_GE7` integer DEFAULT NULL,
  `NUMBER_READS_WITHTRS` integer DEFAULT NULL,
  `NUMBER_READS_WITHTRS_GE7` integer DEFAULT NULL,
  `NUMBER_READS_WITHTRS_GE7_AFTER_REDUND` integer DEFAULT NULL,
  `NUMBER_TRS_IN_READS_AFTER_REDUND` integer DEFAULT NULL,
  `NUMBER_REF_TRS` integer DEFAULT NULL,
  `NUMBER_REFS_TRS_AFTER_REDUND` integer DEFAULT NULL,
  `CLUST_NUMBER_OF_PROCLU_CLUSTERS` integer DEFAULT NULL,
  `CLUST_NUMBER_OF_PROCLU_CLUSTERS_BEFORE_REJOIN` integer DEFAULT NULL,
  `CLUST_NUMBER_OF_EXACTPAT_CLUSTERS` integer DEFAULT NULL,
  `CLUST_NUMBER_OF_REF_REPS_IN_CLUSTERS` integer DEFAULT NULL,
  `CLUST_NUMBER_OF_READ_REPS_IN_CLUSTERS` integer DEFAULT NULL,
  `CLUST_LARGEST_NUMBER_OF_TRS_IN_PROCLU_CLUSTER` integer DEFAULT NULL,
  `CLUST_LARGEST_NUMBER_OF_REFS_IN_PROCLU_CLUSTER` integer DEFAULT NULL,
  `CLUST_LARGEST_PATRANGE_IN_PROCLU_CLUSTER` integer DEFAULT NULL,
  `CLUST_LARGEST_NUMBER_OF_TRS_IN_EXACTPAT_CLUSTER` integer DEFAULT NULL,
  `CLUST_LARGEST_NUMBER_OF_REFS_IN_EXACTPAT_CLUSTER` integer DEFAULT NULL,
  `CLUST_NUMBER_OF_REFS_WITH_PREDICTED_VNTR` integer DEFAULT NULL,
  `CLUST_NUMBER_OF_CLUSTERS_WITH_PREDICTED_VNTR` integer DEFAULT NULL,
  `NUMBER_REFS_VNTR_SPAN_N` integer DEFAULT NULL,
  `NUMBER_REFS_SINGLE_REF_CLUSTER` integer DEFAULT NULL,
  `NUMBER_REFS_SINGLE_REF_CLUSTER_WITH_READS_MAPPED` integer DEFAULT NULL,
  `NUMBER_REFS_SINGLE_REF_CLUSTER_WITH_NO_READS_MAPPED` integer DEFAULT NULL,
  `NUMBER_MAPPED` integer DEFAULT NULL,
  `NUMBER_RANK` integer DEFAULT NULL,
  `NUMBER_RANKFLANK` integer DEFAULT NULL,
  `INTERSECT_RANK_AND_RANKFLANK` integer DEFAULT NULL,
  `INTERSECT_RANK_AND_RANKFLANK_BEFORE_PCR` integer DEFAULT NULL,
  `BBB_WITH_MAP_DUPS` integer DEFAULT NULL,
  `BBB` integer DEFAULT NULL,
  `RANK_EDGES_OVERCUTOFF` integer DEFAULT NULL,
  `RANK_REMOVED_SAMEREF` integer DEFAULT NULL,
  `RANK_REMOVED_SAMESEQ` integer DEFAULT NULL,
  `RANK_REMOVED_PCRDUP` integer DEFAULT NULL,
  `RANKFLANK_EDGES_INSERTED` integer DEFAULT NULL,
  `RANKFLANK_REMOVED_SAMEREF` integer DEFAULT NULL,
  `RANKFLANK_REMOVED_SAMESEQ` integer DEFAULT NULL,
  `RANKFLANK_REMOVED_PCRDUP` integer DEFAULT NULL,
  `TIME_MYSQLCREATE` integer DEFAULT NULL,
  `TIME_TRF` integer DEFAULT NULL,
  `TIME_RENUMB` integer DEFAULT NULL,
  `TIME_REDUND` integer DEFAULT NULL,
  `TIME_PROCLU` integer DEFAULT NULL,
  `TIME_JOINCLUST` integer DEFAULT NULL,
  `TIME_DB_INSERT_REFS` integer DEFAULT NULL,
  `TIME_DB_INSERT_READS` integer DEFAULT NULL,
  `TIME_WRITE_FLANKS` integer DEFAULT NULL,
  `TIME_MAP_FLANKS` integer DEFAULT NULL,
  `TIME_MAP_REFFLANKS` integer DEFAULT NULL,
  `TIME_MAP_INSERT` integer DEFAULT NULL,
  `TIME_EDGES` integer DEFAULT NULL,
  `TIME_INDEX_PCR` integer DEFAULT NULL,
  `TIME_PCR_DUP` integer DEFAULT NULL,
  `TIME_MAP_DUP` integer DEFAULT NULL,
  `TIME_VNTR_PREDICT` integer DEFAULT NULL,
  `TIME_ASSEMBLYREQ` integer DEFAULT NULL,
  `TIME_REPORTS` integer DEFAULT NULL,
  `TIME_DB_CONVERSION_AND_READ_COMPRESSION` integer DEFAULT NULL,
  `DATE_MYSQLCREATE` text DEFAULT NULL,
  `DATE_TRF` text DEFAULT NULL,
  `DATE_RENUMB` text DEFAULT NULL,
  `DATE_REDUND` text DEFAULT NULL,
  `DATE_PROCLU` text DEFAULT NULL,
  `DATE_JOINCLUST` text DEFAULT NULL,
  `DATE_DB_INSERT_REFS` text DEFAULT NULL,
  `DATE_DB_INSERT_READS` text DEFAULT NULL,
  `DATE_WRITE_FLANKS` text DEFAULT NULL,
  `DATE_MAP_FLANKS` text DEFAULT NULL,
  `DATE_MAP_REFFLANKS` text DEFAULT NULL,
  `DATE_MAP_INSERT` text DEFAULT NULL,
  `DATE_EDGES` text DEFAULT NULL,
  `DATE_INDEX_PCR` text DEFAULT NULL,
  `DATE_PCR_DUP` text DEFAULT NULL,
  `DATE_MAP_DUP` text DEFAULT NULL,
  `DATE_VNTR_PREDICT` text DEFAULT NULL,
  `DATE_ASSEMBLYREQ` text DEFAULT NULL,
  `DATE_REPORTS` text DEFAULT NULL,
  `DATE_DB_CONVERSION_AND_READ_COMPRESSION` text DEFAULT NULL,
  `ERROR_STEP` integer NOT NULL DEFAULT '0',
  `ERROR_DESC` varchar(500) NOT NULL DEFAULT '',
  `ERROR_CODE` integer NOT NULL DEFAULT '0',
  `N_MIN_SUPPORT` integer NOT NULL DEFAULT '0',
  `MIN_FLANK_REQUIRED` integer NOT NULL DEFAULT '0',
  `MAX_FLANK_CONSIDERED` integer NOT NULL DEFAULT '0',
  PRIMARY KEY (`id`)
);
CREATE TABLE IF NOT EXISTS `vntr_support` (
  `refid` integer NOT NULL,
  `copies` integer NOT NULL,
  `sameasref` integer NOT NULL,
  `support` integer NOT NULL DEFAULT '0',
  `copiesfloat` float NOT NULL,
  `representative` integer DEFAULT NULL,
  PRIMARY KEY (`refid`,`copies`)
);
CREATE INDEX IF NOT EXISTS "idx_map_read_index" ON "map" (`readid`);
CREATE INDEX IF NOT EXISTS "idx_replnk_sid" ON "replnk" (`sid`);
CREATE INDEX IF NOT EXISTS "idx_vntr_support_read_index" ON "vntr_support" (`representative`);
CREATE INDEX IF NOT EXISTS "idx_vntr_support_copies_support" ON "vntr_support" (`copies`,`support`);
--CREATE INDEX IF NOT EXISTS "idx_fasta_ref_reps_head" ON "fasta_ref_reps" (`head`);
CREATE INDEX IF NOT EXISTS "idx_flank_connection_paramsetid" ON "flank_connection" (`paramsetid`);
DELETE FROM stats;
INSERT INTO stats DEFAULT VALUES;
END TRANSACTION;
