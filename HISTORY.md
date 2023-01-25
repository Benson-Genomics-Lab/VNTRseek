# Version 2.0.3 - Jan 25, 2023

- added proper default in `defaults.vs.cnf` for `REFERENCE_INDIST_PRODUCE`
  - fixed blank entry in complete `$OUTPUT_DIR/vntr_${RUN_NAME}` config
- fixed sqlite version printing
- now removes incomplete `encoded_reads` table
- removed ALGNURL field from vcf output
- added stats time and date entries for step 20 `CONVERSION_AND_READ_COMPRESSION`
- fixed GEN_CONFIG option
- Fixed trf2proclu
- Revised updaterefs.pl
  - READ_LENGTH wired in
  - print_latex, print_distr, commify, nowhitespace removed
    - latex, graph production code extracted
  - removed `VNTRPIPE_` from database in vcf output
  - removed SERVER argument
- Isolated non-pipeline data summary scripts into perl/lib/Summarize
- Re-configured input file parsing
  - support .fa, .fq, .sam, .cram extensions
  - allow mixing .fa and .fasta, .fq and .fastq
- Added cram handling to documentation
- Dropped bedtools requirement, it's not used
- Changed default `OUTPUT_DIR` to `.`
  - simplify output structure for simple use cases
- Changed reduced_db.sql script calling to be from DBI
  - extracted `ATTACH DB` line to remove renaming step
    - previously created relatively (meaning in install dir, probably not writable!)
    - using final db name will be absolute (meaning in output_dir)
- Added last last last thing to migrate date and time stamps for step 20 to reduced database
- Merged `extra_index.pl` into `pcr_dup.pl`
  - Merged step 15 into 14 in `vntrseek.pl`
  - pending `stats` table adjustment for defunct `INDEX_PCR` stats
- Removed file greps from `pcr_dup.pl`
- Removed some non-capturing grouping from greps in insert_reads.pl
- Wrote ouput structure description
- Internalized edlib

# Dec 19, 2022
- fixed imports in `compress_reads.pl`

# Version 2.0.2 - Nov 3, 2022

- added `CLEANUP` option to remove intermediates on success
- changed `DBSUFFIX` to `RUN_NAME`
- changed `OUTPUT_ROOT` to `OUTPUT_DIR`
- redirected `run_edges.pl` and `run_flnkcomp.pl` output
  - removed now unused `TMPDIR` option
- removed unused `CLEAN` option
- removed obselete `stats.pl`
- consolidated error checks after subscript calls into subprocedure
   - replaced `SetError` sub with `FlagError`
- removed unused tmp paramater from `run_edges.pl`
- hid (commented out) `map_dup.pl` output file `DBSUFFIX.map_dup.txt`
- stopped `all.clusters` duplication in `insert_reads.pl`
  - removed now unused rotatedfolder parameter
- removed second run of `set_db_stats.pl`
- fixed time stamps around all steps
- raised `SeqReader` out of `VNTRseek` folder
- revised all SQLite statements for consistency, readability, and efficiency
- created step 20 for database reduction steps
  - combination `scores` table from `map`, `rank`, and `rankflank`
  - implemented read compression in perl, potential new `compressed_reads` table
  - header reduction for fasta reads matching standard Illumina format
- added [edlib](https://github.com/KylerAKA-BU/edlib) to build
  - modified lots of files for this
- changed how `SeqReader.pm` locates `seqtk`
- created `FinishStep` function to reduce code repetition
- removed statistics and blocks for empty steps (6, 7, 11, 18)

---
# Version 2.0.1 - Jul 26, 2022

- fixed `NPROCESSES` internally so config is accurate
- made support for relative input, output, and temp paths
- added default relative paths for output and temp paths
- separated `DBSUFFIX` and config paths
- made blank `DBSUFFIX` field in `defaults.vs.cnf`
- added `DBSUFFIX` to validation
- removed `BACKEND` from `defaults.vs.cnf`
- removed unnecssary entries in `%opts`/`%VSCNF_FILE`
  - `REFERENCE` extensions
  - `BACKEND`
- restructured parameter backwards compatibility
- unified parameter default value handling
- revised config validation messages
- sorted parameter handling sections
- cleaned `print_config` and usage
- added `GEN_CONFIG` option
- revised start and end step argument parsing
- fixed subscript config loading and unified preambles
- reduced stderr spam in steps 1, 4, 8, 15, 17
- reduced stdout spam in steps 3, 5, 10, 14, 15, 16
  - and in c scripts: `redund2.c`, `joinc.c`, `flankalign.c`
- rewired `map_dup.pl` and `run_flankcomp.pl` to free up `STDOUT`
- restructured `run_trf_ng.pl` to match other subscripts
  - argument passing, output
- changed `STDERR`s to `STDOUT` where appropriate
  - warns, timestamps
- fixed VNTR vcf output skipping first entry
- (Dr. Benson) fixed redund bug mishandling profiles
- added output database minimization to cleanup
- removed unused arguments in subscript calls
  - notably insert_reads.pl
- added output preclearing to all steps
- moved final output to top level
- removed many redundant output checks after `vs_db_insert` calls
- removed `SERVER` symlinking
- removed obselete scripts
- cleaned out obselescent commented code
- updated help files

---
# Version 2.0.0 - Jul 11, 2022

- cleaned up of `run_edges.pl` to reduce file io
- instated new `redund2.c` to reduce time
- removed overwrite of config file every run
- merged disaparte changes in stale dev branches

---
# Version 1.10.0-rc.4 - 2019

- See previous releases for older changes
