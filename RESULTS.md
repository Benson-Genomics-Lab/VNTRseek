# Output Description

```
{OUTPUT_DIR}/
    vntr_{RUN_NAME}/
        {RUN_NAME}.db
        {RUN_NAME}_rl{READ_LENGTH}.db

        {RUN_NAME}.vs.cnf

        {RUN_NAME}.span2.vcf
        {RUN_NAME}.allwithsupport.span2.vcf

        data_out/
            [TRF, trf2proclu results]
        data_out_clean/
            [data processing]
```

VNTRseek nests all results for a particular data set
in a folder called `vntr_{RUN_NAME}`. For large scale processing
over multiple data sets, an optional parameter `OUTPUT_DIR` is exposed
for relocating the run folder. This way, the `OUTPUT_DIR` can be
set in a config used for many data sets to redirect them all.
By default, `vntr_{RUN_NAME}` will be created wherever
VNTRseek is called.

The `{RUN_NAME}.db` is the database constructed during the main processing,
containing the complete fasta_reads and various other processing data.

The `{RUN_NAME}_rl{READ_LENGTH}.db` is the final output db.
It has fewer tables, compressed fasta reads,
and if the input used standard illumina headers,
shortened headers.

The `{RUN_NAME}.vs.cnf` contains a config containing the 
last settings used for this data set. It can be provided with
the `--CONFIG` option to reproduce the last run.

The vcf files contain the called VNTRs and any ref TRs with support,
respectively.

The `data_out*` folders contains intermediate files
for data used across multiple steps.
