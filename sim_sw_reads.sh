#!/bin/bash

## ./sim_sw_reads.sh ref trs read flank

REFERENCE="$1"
# REFERENCE=/project/vntrseek/share/GRCh38/GRCh38.fa
SEQ_FILE="$2"
READ_LENGTH="$3"
FLANK_LENGTH="$4"
BUFFER=5 # buffer (in bases) around start/end positions

if [ "$READ_LENGTH" -lt $((FLANK_LENGTH + BUFFER)) ]
then
	echo "Read length must be longer than the flank length + $BUFFER";
	exit;
fi

if [ ! -e "$REFERENCE".fai ]
then
	echo "Creating FASTA index for reference..."
	samtools faidx "$REFERENCE"
fi
# OFFSET=0

TMPFILE=$(mktemp)
# For each TR, produce two windows:
# - one where the read ends at the TR start and goes until the read ends at TR+flank
# - one where the read starts at the TR start-flank and goes until the read starts at TR end
# Start position is just $2-(R+B)-1 since bedtools start pos is 0-based
awk -F "," -v R="$READ_LENGTH" -v F="$FLANK_LENGTH" -v B=BUFFER '(NR>1){
	if (($3-$2)+1+2*(F+B) < R ) {
		print $5 "\t" ($2-(R+B)-1) "\t" ($2+B) "\t" $1;
		# Break up window into two for smaller TRs
		print $5 "\t" ($2-(F+B)-1) "\t" ($3+B+R) "\t" $1;
	}
	else {
		print $5 "\t" ($2-(R+B)-1) "\t" ($3+B+R) "\t" $1;
	}
}' "$SEQ_FILE" > "$TMPFILE"

bedtools getfasta -fi "$REFERENCE" -bed "$TMPFILE" | awk 'NR%2==0{ print}' | paste "$TMPFILE" - | \
awk -v R="$READ_LENGTH" -v F="$FLANK_LENGTH" '{
	if(length($5)>=R){
		for (i=0; i<=length($5)-R; i++) {
			print ">" $1 ":" ($2+i) "-" ($2+i+R) "_" $4 "\n" substr($5,i,R)
		}
	}
}'

rm "$TMPFILE"
