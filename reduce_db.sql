-- This is an sqlite3 sql script to convert a VNTRseek sqlite3 database
-- to a smaller and more efficient version to use with VNTRview
-- 
-- USAGE:
-- $ sqlite3 source_database < conversion_script (this script)
-- eg. 
-- $ sqlite3 hg00096.db < newdb.sql

-- Note that the name of the new database is hardcoded below and this has
-- to be made into a variable


-- ATTACH DATABASE 'temp_reduced.db' as 'newdb';
-- ATTACH DATABASE '{modified_database}' as 'newdb';

BEGIN TRANSACTION;

--*******************************************************
-- rank table
-- fields ties, refdir not used in any queries
-- multi-field primary key, use WITHOUT ROWID
-- ***field score is only used in mapped_reads queries for concatenated fields named Info (in some queries), or unnamed (in other queries); the Info field seems to be used for file naming in the fasta section of server.py

DROP TABLE IF EXISTS newdb.rank;

--CREATE TABLE newdb.rank (
--  refid integer NOT NULL,
--  readid integer NOT NULL,
--  score float DEFAULT NULL,
  -- ties integer NOT NULL DEFAULT 0,
  -- refdir char(1) NOT NULL,
--  PRIMARY KEY (refid, readid)
--) WITHOUT ROWID;

--INSERT INTO  newdb.rank (refid, readid, score) SELECT refid, readid, score FROM rank;

--*******************************************************
-- rankflank table
-- field ties is not used in any query
-- multi-field primary key, use WITHOUT ROWID
-- ***field score is only used in mapped_reads queries for concatenated fields named Info (in some queries), or unnamed (in other queries); the Info field seems to be used for file naming in the fasta section of server.py

DROP TABLE IF EXISTS newdb.rankflank;

--CREATE TABLE newdb.rankflank (
--  refid integer NOT NULL,
--  readid integer NOT NULL,
--  score float DEFAULT NULL,
  -- ties integer NOT NULL DEFAULT 0,
--  PRIMARY KEY (refid, readid)
--) WITHOUT ROWID;

-- INSERT INTO  newdb.rankflank (refid, readid, score) SELECT refid, readid, score FROM rankflank;

--*******************************************************
-- rankandrankflank table
-- new table, to avoid join of rank and rankflank in reference.sql query
-- may drop this table

DROP TABLE IF EXISTS newdb.rankandrankflank;

-- CREATE TABLE newdb.rankandrankflank (
--   refid integer NOT NULL,
--   readid integer NOT NULL,
--   -- rankscore float DEFAULT NULL,  -- do not need, not used in reference.sql
--   -- rankflankscore float DEFAULT NULL, -- do not need, not used in reference.sql
--   PRIMARY KEY (refid,readid)
-- ) WITHOUT ROWID;
--
-- INSERT INTO  newdb.rankandrankflank (refid, readid) SELECT refid, readid FROM rank JOIN rankflank USING (refid, readid);

--*******************************************************
-- map table
-- field bbb needed to determine which reads are mapped bbb,
-- -- KA: 8/22
--     "but is never used in the new database once refidcounts table is created, so remove"
--    I don't think this is true
-- field reserved is 1 if the mapping is for a vntr
-- field reserved2 not used in any queries
-- multi-field primary key, use WITHOUT ROWID
-- XXX might be XXX is combined with other tables with same multi-field key
-- Analysis of queries shows that map has lookup of refid and a join using readid, but only after finding the refid, so can drop readid index

DROP TABLE IF EXISTS newdb.map;

--CREATE TABLE newdb.map (
--  refid integer NOT NULL,
--  readid integer NOT NULL,
--  reserved integer NOT NULL,
  -- reserved2 integer NOT NULL,
--  bbb integer NOT NULL DEFAULT 0,
--  PRIMARY KEY (refid, readid)
--) WITHOUT ROWID;

--INSERT INTO  newdb.map (refid, readid, reserved, bbb) SELECT refid, readid, reserved, bbb FROM map;
--INSERT INTO  newdb.map (refid, readid, reserved) SELECT refid, readid, reserved FROM map;

--CREATE INDEX "idx_map_read_index" ON "map" (`readid`);

--*******************************************************
-- scores new table
-- combines map, rank, and rankflank tables
-- field bbb needed to determine which reads are mapped bbb
-- field reserved is 1 if the mapping is for a vntr

-- multi-field primary key, use WITHOUT ROWID

DROP TABLE IF EXISTS newdb.scores;

CREATE TABLE newdb.scores (
  refid integer NOT NULL,
  readid integer NOT NULL,
  reserved integer NOT NULL,
  bbb integer NOT NULL DEFAULT 0,
  rank_score float DEFAULT NULL,
  rankflank_score float DEFAULT NULL,
  PRIMARY KEY (refid, readid)
) WITHOUT ROWID;

INSERT INTO newdb.scores (refid, readid, reserved, bbb, rank_score, rankflank_score)
SELECT refid, readid, reserved, bbb, rank.score, rankflank.score
FROM map LEFT JOIN rank using(refid, readid) LEFT JOIN rankflank USING(refid, readid);


--*******************************************************
-- refidcounts new table
-- used to speed up references query
-- and to add BBB mappings without ties

DROP TABLE IF EXISTS newdb.refidcounts;

CREATE TABLE newdb.refidcounts (
  refid integer NOT NULL,
  mapcount integer,
  rankcount integer,
  rankflankcount integer,
  bbbcount integer,
  support text,
  PRIMARY KEY (refid)
);


INSERT INTO newdb.refidcounts (refid, mapcount, rankcount, rankflankcount, bbbcount, support) 
SELECT rid, 
(SELECT COUNT(*) FROM map WHERE map.refid=rid) as mapcount,
(SELECT COUNT(*) FROM rank WHERE rank.refid=rid) as rankcount,
(SELECT COUNT(*) FROM rankflank WHERE rankflank.refid=rid) as rankflankcount,
(SELECT COUNT(*) FROM map WHERE map.refid=rid and bbb=1) as bbbcount,
(SELECT GROUP_CONCAT(vntr_support.copies || "(" || support || ")") FROM vntr_support WHERE vntr_support.refid=-rid ORDER BY vntr_support.copies) as support
FROM fasta_ref_reps; 

--*******************************************************
-- clusterlnk table
-- field reserved is VNTR_copy_diff in reference.sql
-- field direction is used in vcf_row...sql and others
-- ***field reserved2 is always zero and is used in mapped_align_reads_format...sql files, but is unnamed
-- clusterid is not used for joins
-- repeatid is unique and thus can be primary key
-- remove multi-field primary key

DROP TABLE IF EXISTS newdb.clusterlnk;

CREATE TABLE newdb.clusterlnk (
  clusterid integer NOT NULL,
  repeatid integer NOT NULL,
  direction char(1) NOT NULL,
  reserved integer NOT NULL,
  -- reserved2 integer NOT NULL,
  --PRIMARY KEY (clusterid,repeatid),
  --UNIQUE (repeatid)
  PRIMARY KEY (repeatid)
);

INSERT INTO  newdb.clusterlnk (clusterid, repeatid, direction, reserved) SELECT clusterid, repeatid, direction, reserved FROM clusterlnk;

-- doesn't work for some reason, syntax error at .
-- CREATE INDEX newdb.idx_clusterlnk_clusterid_index ON clusterlnk (clusterid);

--*******************************************************
-- replnk table
-- all joins use rid
-- no joins use sid, so index is unneeded
-- fields profile, profilerc, profsize, patsize are not used in any queries

DROP TABLE IF EXISTS newdb.replnk;

CREATE TABLE newdb.replnk (
  rid integer NOT NULL,
  sid integer  NOT NULL,
  first integer NOT NULL,
  last integer NOT NULL,
  -- patsize integer NOT NULL,
  copynum float NOT NULL,
  pattern text NOT NULL,
  -- profile text COLLATE BINARY,
  -- profilerc text COLLATE BINARY,
  -- profsize integer DEFAULT NULL,
  PRIMARY KEY (rid)
);
--CREATE INDEX "idx_replnk_sid" ON "replnk" (sid);

INSERT INTO  newdb.replnk (rid, sid, first, last, copynum, pattern) SELECT rid, sid, first, last, copynum, pattern FROM replnk;

--*******************************************************
-- fasta_reads table
-- field head is not used in a join, unique index can be removed
-- field qual is not used in any query

DROP TABLE IF EXISTS newdb.fasta_reads;

CREATE TABLE newdb.fasta_reads (
  sid integer  NOT NULL,
  head varchar(500) NOT NULL,
  dna varchar(8000) DEFAULT NULL,
  -- qual varchar(8000) DEFAULT NULL,
  PRIMARY KEY (sid)
  -- UNIQUE (head)
);

INSERT INTO  newdb.fasta_reads (sid, head, dna) SELECT sid, head, dna FROM fasta_reads;

--*******************************************************
-- clusters table
-- the following fields are empty: refs_flank_undist, refs_flank_dist, refs_flank_undist_l, refs_flank_undist_lr, refs_flank_undist_r
-- the following fields are empty, but used in cluster.sql and cluster_count.sql: assemblyreq, flankdensity, mcpattern, aveentropy
-- they are not used in the cluster.html page, so remove and let it fail
-- they are used in the filtering and ordering by parts of the html pages, so remove

DROP TABLE IF EXISTS newdb.clusters;

CREATE TABLE newdb.clusters (
  cid integer NOT NULL,
  minpat integer NOT NULL,
  maxpat integer NOT NULL,
  -- refs_flank_undist integer DEFAULT NULL,
  -- refs_flank_dist integer DEFAULT NULL,
   -- refs_flank_undist_l integer DEFAULT NULL,
   -- refs_flank_undist_r integer DEFAULT NULL,
   -- refs_flank_undist_lr integer DEFAULT NULL,
  repeatcount integer NOT NULL,
  refcount integer NOT NULL,
  variability integer NOT NULL DEFAULT '0',
  -- assemblyreq integer DEFAULT NULL,
  profdensity float DEFAULT NULL,
  -- flankdensity float DEFAULT NULL,
  -- mcpattern varchar(5000) DEFAULT NULL,
  -- aveentropy float DEFAULT NULL,
  PRIMARY KEY (cid)
);

INSERT INTO  newdb.clusters (cid, minpat, maxpat, repeatcount, refcount, variability, profdensity) SELECT cid, minpat, maxpat, repeatcount, refcount, variability, profdensity FROM clusters;

--*******************************************************
-- vntr_support table
-- all joins use refid, but it is not unique
-- multi-field primary key, use WITHOUT ROWID
-- no joins use representative, copies, or support so remove indexes
-- ***field copiesfloat is not used in any query

DROP TABLE IF EXISTS newdb.vntr_support;

CREATE TABLE newdb.vntr_support (
  refid integer NOT NULL,
  copies integer NOT NULL,
  sameasref integer NOT NULL,
  support integer NOT NULL DEFAULT '0',
  copiesfloat float NOT NULL,
  representative integer DEFAULT NULL,
  PRIMARY KEY (refid,copies)
) WITHOUT ROWID;

-- CREATE INDEX "idx_vntr_support_read_index" ON "vntr_support" (`representative`);
-- CREATE INDEX "idx_vntr_support_copies_support" ON "vntr_support" (`copies`,`support`);

INSERT INTO  newdb.vntr_support (refid, copies, sameasref, support, copiesfloat, representative) SELECT refid, copies, sameasref, support, copiesfloat, representative FROM vntr_support;

--*******************************************************
-- fasta_ref_reps table (sample)
-- all joins use rid
-- unique index not required
-- fields not used: flank_disting, has_support, span1, spanN, support_vntr, support_vntr_span1
-- fields maybe will be used in future: homez..., hetez...  
-- field not required: comment

DROP TABLE IF EXISTS newdb.fasta_ref_reps;

CREATE TABLE newdb.fasta_ref_reps (
  rid integer NOT NULL,
  -- comment varchar(500) DEFAULT NULL,
  -- flank_disting integer DEFAULT NULL,
  entropy float NOT NULL DEFAULT '0',
  -- has_support integer DEFAULT NULL,
  -- span1 integer DEFAULT NULL,
  -- spanN integer DEFAULT NULL,
  homez_same integer DEFAULT NULL,
  homez_diff integer DEFAULT NULL,
  hetez_same integer DEFAULT NULL,
  hetez_diff integer DEFAULT NULL,
  hetez_multi integer DEFAULT NULL,
  -- support_vntr integer DEFAULT NULL,
  -- support_vntr_span1 integer DEFAULT NULL,
  alleles_sup integer DEFAULT NULL,
  allele_sup_same_as_ref integer DEFAULT NULL,
  PRIMARY KEY (rid)
  --UNIQUE (rid,comment)
);

INSERT INTO  newdb.fasta_ref_reps (rid, entropy, homez_same, homez_diff, hetez_same, hetez_diff, hetez_multi, alleles_sup, allele_sup_same_as_ref) SELECT rid, entropy, homez_same, homez_diff, hetez_same, hetez_diff, hetez_multi, alleles_sup, allele_sup_same_as_ref FROM fasta_ref_reps;

--*******************************************************
-- stats table

DROP TABLE IF EXISTS newdb.stats;

CREATE TABLE newdb.stats (
  id integer NOT NULL,
  MAP_ROOT varchar(500) DEFAULT NULL,
  PARAM_TRF varchar(500) DEFAULT NULL,
  PARAM_PROCLU varchar(500) DEFAULT NULL,
  FOLDER_FASTA varchar(500) DEFAULT NULL,
  FOLDER_PROFILES varchar(500) DEFAULT NULL,
  FOLDER_PROFILES_CLEAN varchar(500) DEFAULT NULL,
  FOLDER_REFERENCE varchar(500) DEFAULT NULL,
  FILE_REFERENCE_LEB varchar(500) DEFAULT NULL,
  FILE_REFERENCE_SEQ varchar(500) DEFAULT NULL,
  NUMBER_READS integer DEFAULT NULL,
  NUMBER_TRS_IN_READS integer DEFAULT NULL,
  NUMBER_TRS_IN_READS_GE7 integer DEFAULT NULL,
  NUMBER_READS_WITHTRS integer DEFAULT NULL,
  NUMBER_READS_WITHTRS_GE7 integer DEFAULT NULL,
  NUMBER_READS_WITHTRS_GE7_AFTER_REDUND integer DEFAULT NULL,
  NUMBER_TRS_IN_READS_AFTER_REDUND integer DEFAULT NULL,
  NUMBER_REF_TRS integer DEFAULT NULL,
  NUMBER_REFS_TRS_AFTER_REDUND integer DEFAULT NULL,
  CLUST_NUMBER_OF_PROCLU_CLUSTERS integer DEFAULT NULL,
  CLUST_NUMBER_OF_PROCLU_CLUSTERS_BEFORE_REJOIN integer DEFAULT NULL,
  CLUST_NUMBER_OF_EXACTPAT_CLUSTERS integer DEFAULT NULL,
  CLUST_NUMBER_OF_REF_REPS_IN_CLUSTERS integer DEFAULT NULL,
  CLUST_NUMBER_OF_READ_REPS_IN_CLUSTERS integer DEFAULT NULL,
  CLUST_LARGEST_NUMBER_OF_TRS_IN_PROCLU_CLUSTER integer DEFAULT NULL,
  CLUST_LARGEST_NUMBER_OF_REFS_IN_PROCLU_CLUSTER integer DEFAULT NULL,
  CLUST_LARGEST_PATRANGE_IN_PROCLU_CLUSTER integer DEFAULT NULL,
  CLUST_LARGEST_NUMBER_OF_TRS_IN_EXACTPAT_CLUSTER integer DEFAULT NULL,
  CLUST_LARGEST_NUMBER_OF_REFS_IN_EXACTPAT_CLUSTER integer DEFAULT NULL,
  CLUST_NUMBER_OF_REFS_WITH_PREDICTED_VNTR integer DEFAULT NULL,
  CLUST_NUMBER_OF_CLUSTERS_WITH_PREDICTED_VNTR integer DEFAULT NULL,
  NUMBER_REFS_VNTR_SPAN_N integer DEFAULT NULL,
  NUMBER_REFS_SINGLE_REF_CLUSTER integer DEFAULT NULL,
  NUMBER_REFS_SINGLE_REF_CLUSTER_WITH_READS_MAPPED integer DEFAULT NULL,
  NUMBER_REFS_SINGLE_REF_CLUSTER_WITH_NO_READS_MAPPED integer DEFAULT NULL,
  NUMBER_MAPPED integer DEFAULT NULL,
  NUMBER_RANK integer DEFAULT NULL,
  NUMBER_RANKFLANK integer DEFAULT NULL,
  INTERSECT_RANK_AND_RANKFLANK integer DEFAULT NULL,
  INTERSECT_RANK_AND_RANKFLANK_BEFORE_PCR integer DEFAULT NULL,
  BBB_WITH_MAP_DUPS integer DEFAULT NULL,
  BBB integer DEFAULT NULL,
  RANK_EDGES_OVERCUTOFF integer DEFAULT NULL,
  RANK_REMOVED_SAMEREF integer DEFAULT NULL,
  RANK_REMOVED_SAMESEQ integer DEFAULT NULL,
  RANK_REMOVED_PCRDUP integer DEFAULT NULL,
  RANKFLANK_EDGES_INSERTED integer DEFAULT NULL,
  RANKFLANK_REMOVED_SAMEREF integer DEFAULT NULL,
  RANKFLANK_REMOVED_SAMESEQ integer DEFAULT NULL,
  RANKFLANK_REMOVED_PCRDUP integer DEFAULT NULL,
  TIME_MYSQLCREATE integer DEFAULT NULL,
  TIME_TRF integer DEFAULT NULL,
  TIME_RENUMB integer DEFAULT NULL,
  TIME_REDUND integer DEFAULT NULL,
  TIME_PROCLU integer DEFAULT NULL,
  TIME_JOINCLUST integer DEFAULT NULL,
  TIME_DB_INSERT_REFS integer DEFAULT NULL,
  TIME_DB_INSERT_READS integer DEFAULT NULL,
  TIME_WRITE_FLANKS integer DEFAULT NULL,
  TIME_MAP_FLANKS integer DEFAULT NULL,
  TIME_MAP_REFFLANKS integer DEFAULT NULL,
  TIME_MAP_INSERT integer DEFAULT NULL,
  TIME_EDGES integer DEFAULT NULL,
  TIME_INDEX_PCR integer DEFAULT NULL,
  TIME_PCR_DUP integer DEFAULT NULL,
  TIME_MAP_DUP integer DEFAULT NULL,
  TIME_VNTR_PREDICT integer DEFAULT NULL,
  TIME_ASSEMBLYREQ integer DEFAULT NULL,
  TIME_REPORTS integer DEFAULT NULL,
  TIME_DB_CONVERSION_AND_READ_COMPRESSION integer DEFAULT NULL,
  DATE_MYSQLCREATE text DEFAULT NULL,
  DATE_TRF text DEFAULT NULL,
  DATE_RENUMB text DEFAULT NULL,
  DATE_REDUND text DEFAULT NULL,
  DATE_PROCLU text DEFAULT NULL,
  DATE_JOINCLUST text DEFAULT NULL,
  DATE_DB_INSERT_REFS text DEFAULT NULL,
  DATE_DB_INSERT_READS text DEFAULT NULL,
  DATE_WRITE_FLANKS text DEFAULT NULL,
  DATE_MAP_FLANKS text DEFAULT NULL,
  DATE_MAP_REFFLANKS text DEFAULT NULL,
  DATE_MAP_INSERT text DEFAULT NULL,
  DATE_EDGES text DEFAULT NULL,
  DATE_INDEX_PCR text DEFAULT NULL,
  DATE_PCR_DUP text DEFAULT NULL,
  DATE_MAP_DUP text DEFAULT NULL,
  DATE_VNTR_PREDICT text DEFAULT NULL,
  DATE_ASSEMBLYREQ text DEFAULT NULL,
  DATE_REPORTS text DEFAULT NULL,
  DATE_DB_CONVERSION_AND_READ_COMPRESSION text DEFAULT NULL,
  ERROR_STEP integer NOT NULL DEFAULT '0',
  ERROR_DESC varchar(500) NOT NULL DEFAULT '',
  ERROR_CODE integer NOT NULL DEFAULT '0',
  N_MIN_SUPPORT integer NOT NULL DEFAULT '0',
  MIN_FLANK_REQUIRED integer NOT NULL DEFAULT '0',
  MAX_FLANK_CONSIDERED integer NOT NULL DEFAULT '0',
  PRIMARY KEY (id)
);

-- copy table without naming fields
INSERT INTO newdb.stats SELECT * FROM stats;

--*******************************************************

COMMIT;

ANALYZE newdb;

-- can't vacuum from within a transaction
-- for vacuum, need to have enough space in the /var/tmp directory
-- During a VACUUM, SQLite creates a temporary file that is approximately the same size as the original database. It does this in order to maintain the database ACID properties. SQLite uses a directory to hold the temporary files it needs while doing the VACUUM. In order to determine what directory to use, it descends a hierarchy looking for the first directory that has the proper access permissions.For Unix (and Linux) the hierarchy is:

--    Whatever is specified by the SQLITE_TMPDIR environment variable,
--    Whatever is specified by the TMPDIR environment variable,
--    /var/tmp,
--    /usr/tmp,
--    /tmp, and finally
--    The current working directory.

-- possible fix (should put these at top)
-- may not work because temp_store_directory is deprecated
-- pragma temp_store = 1;
-- pragma temp_store_directory = '/var/www/VNTRview_1.09/db';
VACUUM newdb;

DETACH newdb;