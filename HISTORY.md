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
- resturctured parameter backwards compatibility
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
