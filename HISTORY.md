---
# ver 2.0.1 - Jul 21, 2022

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
- fixed subscript config loading and unified preambles
- reduced stderr spam in steps 1, 4, 8, 15, 17
- reduced stdout spam in steps 3, 5, 10, 14, 15, 16
  - c scripts: redund.c, joinc.c, flankalign.c
- rewired `map_dup.pl` and `run_flankcomp.pl` to free up `STDOUT`
- changed `STDERR`s to `STDOUT`
- fixed VNTR vcf output skipping first entry
- determined new redund has affected fasta_reads table
- added output database minimization to cleanup
- removed obselete scripts
- cleaned out obselescent commented code
- updated help files

---
# ver 2.0.0 - Jul 11, 2022

- Cleanup of run_edges.pl to reduce file io
- New redund.c to reduce time
- No longer overwrites config file every run
- merge of disaparte changes in stale dev branches

---
# ver 1.10.0-rc.4 - 2019

- See previous releases for changes