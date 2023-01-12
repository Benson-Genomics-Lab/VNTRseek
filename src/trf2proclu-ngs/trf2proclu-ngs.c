/****************************************************************
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Please direct all questions related to this program or bug reports or updates
 *to Yevgeniy Gelfand (ygelfand@bu.edu)
 *
 ****************************************************************/

#include <assert.h>
#include <ctype.h>
#include <getopt.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

/* converts TRF .dat file -> PROCLU .leb32 file */
/* version 1.00 */
/* started by Yevgeniy Gelfand on July 09, 2009 */

//#define _WIN_32_YES

#include "../libs/easylife/easylife.h"
#include "patupdt.h"
#include "profile.h"

// if this is set to 1, profiles are not rotated during sorting (must correspond
// to -i flag in redund.c)
#define IDENTICAL_ONLY ( 1 )

// Warning:
// When modifying LBI_MAX_HEADER_SIZE, must also modify TRF_HEADER_FORMAT.
#define LBI_MAX_HEADER_SIZE 200
#define TRF_HEADER_FORMAT "@%200[^\n\r]"
#define LBI_MAX_FLANK_SIZE 1000

// Warning:
// When modifying LBI_MAX_PATTERN_SIZE or LBI_MAX_SEQUENCE_SIZE, must also
// modify TRF_RECORD_FORMAT.
#define LBI_MAX_PATTERN_SIZE 50000
#define LBI_MAX_SEQUENCE_SIZE 520000
#define TRF_RECORD_FORMAT \
    "%d %d %d %f %d %d %d %d %d %d %d %d %f %50000s %520000s %1000s %1000s"

typedef struct {

    char *         header;
    unsigned char *pattern;
    unsigned char *sequence;
    unsigned char *left;
    unsigned char *right;

    int *    minRepresentation;
    PROFILE *prof;
    PROFILE *profrc;
    int      minrlen;

    int   id;
    int   matchperc;
    int   indelperc;
    int   score;
    int   acount;
    int   ccount;
    int   gcount;
    int   tcount;
    int   acgtCount;
    float copynum;
    float conserved;
    float entropy;

    int firstindex;
    int lastindex;
    int period;
    int patsize;
    int spanning;

} REP_STRUCT;

static char *usage =
  "Usage: trf2proclu-ngs.exe -f <num> -m <num> -s <num> -i <num> -o <string> "
  "-p <num> "
  "-l <num> [input.dat]\nWhere:\n\t-f specifies the id assigned to the "
  "first record in the file (must be greater than or equal to "
  "1),\n\t-m must be equal to the matching weight parameter of the "
  "corresponding TRF run,\n\t-s must be equal to the mismatch penalty "
  "parameter of the corresponding TRF run,\n\t-i must be equal to the "
  "indel penalty parameter of the corresponding TRF run,\n\t-o "
  "specifies the prefix of the output files,\n\t-p specifies "
  "minimum patsize to keep TR,\n\t-l specifies minimum flanksize "
  "(either side) to keep TR\ntrf2proclu-ngs.exe reads a DAT file output by the "
  "TRF and "
  "produces an LEB36 file and an index file that contains records from "
  "the DAT file each preceded by a unique id.\n";

/********************************	doCriticalErrorAndQuit
 * ****************************/
void doCriticalErrorAndQuit( const char *format, ... ) {
    /* needed for EasyLife error handling */
    va_list argp;

    if ( format == NULL )
        exit( -1 );

    fputs( "ERROR: ", stderr );
    va_start( argp, format );
    fprintf( stderr, format, argp );
    va_end( argp );
    exit( -1 );
}

/********************************* pattern_sort
 * ************************************/
int pattern_sort( const void *item1, const void *item2 ) {

    int p1, p2;

    p1 = ( (REP_STRUCT *) item1 )->patsize;
    p2 = ( (REP_STRUCT *) item2 )->patsize;

    if ( p1 > p2 )
        return 1;
    else if ( p1 < p2 )
        return -1;
    else {
        // tie resolution by id
        p1 = ( (REP_STRUCT *) item1 )->id;
        p2 = ( (REP_STRUCT *) item2 )->id;

        if ( p1 > p2 )
            return 1;
        else if ( p1 < p2 )
            return -1;
    }

    return 0;
}

/*******************************************************************************************/
int pindcmp( int *p1, int *p2, int len1, int len2 ) {

    int i, len;

    if ( len1 > len2 )
        return 1;
    else if ( len1 < len2 )
        return -1;

    len = min( len1, len2 );

    for ( i = 0; i < len; i++ ) {
        if ( p1[i] > p2[i] )
            return 1;
        else if ( p1[i] < p2[i] )
            return -1;
    }

    return 0;
}

/*******************************************************************************************/
int arsize_and_min_rep_sort( const void *item1, const void *item2 ) {

    REP_STRUCT *p1, *p2;

    p1 = ( (REP_STRUCT *) item1 );
    p2 = ( (REP_STRUCT *) item2 );

    /*
        if (p1->acgtCount > p2->acgtCount)
            return 1;
        else if (p1->acgtCount < p2->acgtCount)
            return -1;
    */

    return pindcmp(
      p1->minRepresentation, p2->minRepresentation, p1->minrlen, p2->minrlen );
}

/*******************************************************************************************/
int *pintdup( int *src, int size ) {
    int i, *res;

    res = smalloc( sizeof( int ) * size );

    for ( i = 0; i < size; i++ ) {
        res[i] = src[i];
    }

    return res;
}

int RC = 0;

/*******************************************************************************************/
int *MinimumRepresentation( int *indptr, int len, int *indptrrc, int lenrc,
  int *minrlen, int identical_only ) {

    int  i, j;
    int *src, *srcrc, *tar, tmp;

    RC = 0;

    if ( NULL == indptr || NULL == indptrrc )
        return NULL;

    src   = pintdup( indptr, len );
    srcrc = pintdup( indptrrc, lenrc );

    tar = pintdup( indptr, len );

    /* added on feb 6, 2013 */
    if ( identical_only ) {

        // minimum of the forward and reverse profiles
        if ( pindcmp( srcrc, tar, lenrc, len ) < 0 ) {
            free( tar );
            tar = pintdup( srcrc, lenrc );
            RC  = 1;
        }

    } else {

        // FORWARD - test each rotation and find the minimum one
        for ( i = 1; i < len; i++ ) {

            // printf("\n\nstart: %s",src);

            // rotate by 1
            tmp = src[0];

            for ( j = 1; j < len; j++ ) {
                src[j - 1] = src[j];
            }

            src[len - 1] = tmp;

            // printf("\nend: %s comp: %d\n\n",src,strcmp(src,tar));
            // exit(1);

            if ( pindcmp( src, tar, len, len ) < 0 ) {
                free( tar );
                tar = pintdup( src, len );
            }
        }

        // REVERSE - test each rotation and find the minimum one
        for ( i = 1; i < lenrc; i++ ) {

            // printf("\n\nstart: %s",src);

            // rotate by 1
            tmp = srcrc[0];

            for ( j = 1; j < len; j++ ) {
                srcrc[j - 1] = srcrc[j];
            }

            srcrc[lenrc - 1] = tmp;

            // printf("\nend: %s comp: %d\n\n",src,strcmp(src,tar));
            // exit(1);

            if ( pindcmp( srcrc, tar, lenrc, len ) < 0 ) {
                free( tar );
                tar = pintdup( srcrc, lenrc );
                RC  = 1;
            }
        }
    }

    // save rotated
    free( src );
    free( srcrc );

    *minrlen = ( 0 == RC ) ? len : lenrc;

    return tar;
}

/********************************	main
 * ********************************************/

// Flag set by --ngs command line option
static int verbose = 0;

int main( int argc, char **argv ) {

    short *pc;
    int    lentemp, theid = 0, comps[5];
    int    i, c, j, error, firstindex, lastindex, period, patsize, matchperc,
      indelperc, score, acount, ccount, gcount, tcount, cRes, lenleft, lenright,
      headerlen, rcyes;
    int   thecount;
    float copynum = 0, entropy = 0;
    FILE *fp, *fph;
    // FILE *         *logfp; // reenable, perhaps with some verbose mode
    char *         header, *flankleft, *flankright;
    unsigned char *sequence, *pattern;
    int            startid = -1000000, match = -1000000, mismatch = -1000000,
        indel = -1000000, minperiod = -1000000, minflank = -1000000, *sm;
    REPEAT      rep;
    REP_STRUCT *repPtr;
    EASY_LIST * repList;
    EASY_NODE * tnode;
    PROFILE *   profptr;
    char *      outfile_prefix;
    char *      leb36file;
    char *      indexfileh;

    header     = smalloc( sizeof( *header ) * ( LBI_MAX_HEADER_SIZE + 1 ) );
    pattern    = smalloc( sizeof( *pattern ) * ( LBI_MAX_PATTERN_SIZE + 1 ) );
    sequence   = smalloc( sizeof( *sequence ) * ( LBI_MAX_SEQUENCE_SIZE + 1 ) );
    flankleft  = smalloc( sizeof( *flankleft ) * ( LBI_MAX_FLANK_SIZE + 1 ) );
    flankright = smalloc( sizeof( *flankright ) * ( LBI_MAX_FLANK_SIZE + 1 ) );

    while ( 1 ) {
        static struct option long_options[] = {// These options set a flag.
          {"verbose", no_argument, &verbose, 1},
          // These options don't set a flag.
          // We distinguish them by their indices.
          {"help", no_argument, NULL, 'h'},
          {"first", required_argument, NULL, 'f'},
          {"match", required_argument, NULL, 'm'},
          {"mismatch", required_argument, NULL, 's'},
          {"indel", required_argument, NULL, 'i'},
          {"minperiod", required_argument, NULL, 'p'},
          {"minflank", required_argument, NULL, 'l'},
          {"output", required_argument, NULL, 'o'}, {0, 0, NULL, 0}};
        int option_index = 0; // getopt_long() stores the option index here
        c                = getopt_long(
          argc, argv, "hf:m:s:i:p:l:o:", long_options, &option_index );

        if ( c == -1 )
            break; // detect the end of the options

        switch ( c ) {
        case 0:
            // If this option set a flag, do nothing else now.
            // if (long_options[option_index].flag != 0) break;
            // fprintf(stderr, "option %s", long_options[option_index].name);
            // if (optarg) fprintf (stderr, " with arg %s", optarg);
            // fprintf(stderr, "\n");
            break;

        case 'h':
            fprintf( stderr, "%s", usage );
            return ( 3 );

        case 'f':
            startid = atoi( optarg );

            if ( startid < 1 )
                fputs(
                  "Starting id must be greater than or equal to 1. Aborting.\n",
                  stderr );

            break;

        case 'm':
            match = atoi( optarg );
            break;

        case 'p':
            minperiod = atoi( optarg );
            break;

        case 'l':
            minflank = atoi( optarg );
            break;

        case 's':
            mismatch = atoi( optarg );

            if ( mismatch < 0 )
                mismatch = -mismatch;

            break;

        case 'i':
            indel = atoi( optarg );

            if ( indel < 0 )
                indel = -indel;

            break;

        case 'o':
            outfile_prefix = optarg;
            // if (indexfile == NULL) {
            //     perror("Could not allocate memory to store indexfile name");
            //     return (4);
            // }

            // strncpy(indexfile, optarg, strlen(optarg));
            break;

        case '?': // getopt_long already printed an error message
            break;

        default: // should never happen
            assert( 0 );
        }
    }

    if ( outfile_prefix == NULL || strlen( outfile_prefix ) == 0 ) {
        fputs( "No name given for output files. Aborting.\n", stderr );
        fprintf( stderr, "%s", usage );
        return ( 5 );
    }

    if ( minperiod <= -1000000 ) {
        fputs( "No minperiod provided. Aborting.\n", stderr );
        fprintf( stderr, "%s", usage );
        return ( 6 );
    }

    if ( minflank <= -1000000 ) {
        fputs( "No minflank provided. Aborting.\n", stderr );
        fprintf( stderr, "%s", usage );
        return ( 7 );
    }

    if ( match <= -1000000 ) {
        fputs( "No match weight provided. Aborting.\n", stderr );
        fprintf( stderr, "%s", usage );
        return ( 8 );
    }

    if ( mismatch <= -1000000 ) {
        fputs( "No mismatch weight provided. Aborting.\n", stderr );
        fprintf( stderr, "%s", usage );
        return ( 9 );
    }

    if ( indel <= -1000000 ) {
        fputs( "No indel weight provided. Aborting.\n", stderr );
        fprintf( stderr, "%s", usage );
        return ( 10 );
    }

    if ( startid <= -1000000 ) {
        fputs( "No startng id provided. Aborting.\n", stderr );
        fprintf( stderr, "%s", usage );
        return ( 11 );
    }

    int nonopt =
      argc - optind; // determine the number of non-option ARGV elements

    if ( nonopt == 0 ) {
        // read from standard input if no non-option arguments specified
        fp = stdin;
    } else if ( nonopt == 1 ) {
        // read from file given by the non-option argument
        fp = fopen( argv[argc - 1], "r" );

        if ( fp == NULL )
            fprintf( stderr,
              "Unable to open file '%s' for reading. Aborting.\n",
              argv[argc - 1] );

        fprintf( stderr, "Reading input file '%s'... ", argv[argc - 1] );
    } else if ( nonopt > 1 ) {
        fputs( "Too many non-option arguments specified. Aborting.\n", stderr );
        return ( 12 );
    } else {         // optind < 0
        assert( 0 ); // should never happen
    }

    // // Open log file
    // char *logfile = calloc( strlen( indexfile ) + 1, sizeof( char ) );
    // strncpy( logfile, indexfile, strlen( indexfile ) - 6 );
    // // fprintf(stderr, "Log file name: %s\n", logfile);
    // strncat( logfile, ".log", 5 );
    // // fprintf(stderr, "Log file name: %s\n", logfile);
    // logfp = fopen( logfile, "w" );

    // if ( logfp == NULL )
    //     fprintf( stderr,
    //       "Unable to open log file %s for writing. TRF progress will not be "
    //       "written out.\n",
    //       logfile );

    // // Flush log file output immediately by disabling buffer
    // setbuf( logfp, NULL );

    // if (optind < argc) {
    //	fprintf (stderr, "non-option ARGV-elements: ");
    //	while (optind < argc)
    //		fprintf (stderr, "%s ", argv[optind++]);
    //	putchar ('\n');
    //}
    // if (nonopt == 1) {
    //	fprintf(stderr, "Closing input file '%s'\n", argv[argc - 1]);
    //	fclose(fp);
    //}
    // abort();

    // skip data file header portion
    // for(i=0;i<14;i++) {
    //	fgets(pattern,LBI_MAX_PATTERN_SIZE,fp);
    //	if(pattern[0]=='P') // get the parameters
    //	{
    //		sscanf(pattern,"%s %d %d %d %d %d
    //%d",sequence,&match,&mismatch,&indel,&PM,&PI,&minscore);
    //	}
    //}

    // read records from file
    repList  = EasyListCreate( NULL, free );
    error    = 0;
    thecount = 0;
    theid = startid; // TODO: this must be incremented in steps of one for every

    // allocate the similarity matrix
    sm = CreateSubstitutionMatrix( match, -mismatch );

    while ( !feof( fp ) ) {
        if ( ferror( fp ) ) {
            fclose( fp );
            fputs( "Error reading from file. Aborting.\n", stderr );
            return ( 13 );
        }

        // scan TRF header
        i = fscanf( fp, TRF_HEADER_FORMAT, header );

        if ( i != 1 ) {
            fclose( fp );
            fputs( "Could not read TRF header. No TRs in input?\n", stderr );
            return ( 14 );
        }

        // scan TRF record(s)
        while ( 17 == fscanf( fp, TRF_RECORD_FORMAT, &firstindex, &lastindex,
                        &period, &copynum, &patsize, &matchperc, &indelperc,
                        &score, &acount, &ccount, &gcount, &tcount, &entropy,
                        pattern, sequence, flankleft, flankright ) ) {

            lenleft  = ( flankleft[0] == '.' ) ? 0 : strlen( flankleft );
            lenright = ( flankright[0] == '.' ) ? 0 : strlen( flankright );

            // printf("%s\n%s\n%s\n\n", header, pattern, sequence);
            /* add to list */
            thecount++;
            repPtr         = (REP_STRUCT *) scalloc( 1, sizeof( REP_STRUCT ) );
            repPtr->header = strdup( header );
            repPtr->id     = theid;
            repPtr->firstindex = firstindex;
            repPtr->lastindex  = lastindex;
            repPtr->period     = period;
            repPtr->copynum    = copynum;
            repPtr->patsize    = patsize;
            repPtr->matchperc  = matchperc;
            repPtr->indelperc  = indelperc;
            repPtr->score      = score;
            repPtr->acount     = acount;
            repPtr->ccount     = ccount;
            repPtr->gcount     = gcount;
            repPtr->tcount     = tcount;
            repPtr->entropy    = entropy;
            repPtr->pattern    = (unsigned char *) strdup( (char *) pattern );
            repPtr->sequence   = (unsigned char *) strdup( (char *) sequence );
            repPtr->spanning   = 1;
            repPtr->prof       = NULL;
            repPtr->profrc     = NULL;

            if ( flankleft[0] == '.' )
                repPtr->left = (unsigned char *) strdup( "" );
            else
                repPtr->left = (unsigned char *) strdup( flankleft );

            if ( flankright[0] == '.' )
                repPtr->right = (unsigned char *) strdup( "" );
            else
                repPtr->right = (unsigned char *) strdup( flankright );

            /* if header ends with _RCYES, we need to rccomp the results first
             * and eliminate duplicates */
            headerlen = strlen( repPtr->header );
            rcyes     = 0;

            if ( headerlen >= 6 &&
                 0 == strcmp( repPtr->header + headerlen - 6, "_RCYES" ) ) {
                int            seqlen = 0;
                char *         src;
                unsigned char *strtemp, *strtemp2;

                rcyes = 1;

                for ( src = repPtr->header + headerlen - 7;
                      src != repPtr->header; src-- ) {
                    if ( *src == '_' ) {
                        *src   = '\0';
                        seqlen = strtol( src + 1, NULL, 10 );
                        break;
                    }
                }

                if ( 0 == seqlen ) {
                    fputs( "No sequence len passed though header of the RC "
                           "read. Aborting.\n",
                      stderr );
                    return ( 15 );
                }

                repPtr->firstindex = seqlen - lastindex + 1;
                repPtr->lastindex  = seqlen - firstindex + 1;
                firstindex         = repPtr->firstindex;
                lastindex          = repPtr->lastindex;

                repPtr->acount = tcount;
                repPtr->tcount = acount;
                repPtr->ccount = gcount;
                repPtr->gcount = ccount;
                acount         = repPtr->acount;
                tcount         = repPtr->tcount;
                ccount         = repPtr->ccount;
                gcount         = repPtr->gcount;
                free( repPtr->pattern );
                repPtr->pattern = GetReverseComplement( pattern );
                free( repPtr->sequence );
                repPtr->sequence = GetReverseComplement( sequence );

                strtemp  = GetReverseComplement( repPtr->left );
                strtemp2 = GetReverseComplement( repPtr->right );
                free( repPtr->left );
                repPtr->left = strtemp2;
                free( repPtr->right );
                repPtr->right = strtemp;
            }

            /* FOR SORTING BY MINREPRESENTATION FOR PIPELINE */
            // set member variables
            rep.updates       = 0;
            rep.concensussize = repPtr->patsize;
            rep.firstindex    = repPtr->firstindex;
            rep.lastindex     = repPtr->lastindex;
            rep.pattern = (unsigned char *) strdup( (char *) repPtr->pattern );
            rep.subsequence =
              (unsigned char *) strdup( (char *) repPtr->sequence );
            rep.profile = NULL;
            rep.proflen = 0;

            // get profile
            cRes = ComputeBetterPattern( &rep, 0, -indel, sm );

            if ( -1 == cRes ) {
                fputs( "Memory allocation failed on ComputeBetterPattern(). "
                       "Aborting.\n",
                  stderr );
                return ( 16 );
            }

            // fix
            repPtr->patsize = rep.concensussize;
            repPtr->copynum =
              (double) ( rep.gappedlength ) / (double) ( rep.proflen );
            repPtr->conserved =
              rep.countmatch /
              (double) ( rep.countmatch + rep.countmismatch + rep.countindel );
            repPtr->matchperc = rep.percentmatch;
            repPtr->indelperc = rep.percentindels;
            repPtr->score     = rep.score;
            free( repPtr->pattern );
            repPtr->pattern = (unsigned char *) strdup( (char *) rep.pattern );

            int         found = 0, complen;
            REP_STRUCT *tempPtr;

            for ( tnode = repList->tail; tnode != NULL; tnode = tnode->prev ) {
                tempPtr = (REP_STRUCT *) ( tnode->item );
                complen = strlen( tempPtr->header );
                complen = min( complen, headerlen );

                if ( 0 != memcmp( repPtr->header, tempPtr->header, complen ) ) {
                    break;
                }

                if ( repPtr->firstindex == tempPtr->firstindex &&
                     repPtr->lastindex == tempPtr->lastindex &&
                     repPtr->patsize == tempPtr->patsize ) {
                    found = 1;
                    break;
                }
            }

            if ( found ) {
                free( rep.subsequence );
                free( rep.pattern );
                free( rep.profile );
                continue;
            }

            /* skip periods less than minperiod or with flanks less than
             * minflank */
            if ( repPtr->patsize < minperiod || lenleft < minflank ||
                 lenright < minflank ) {
                repPtr->spanning = 0;
            }

            /* insert into list */
            EasyListInsertTail( repList, repPtr );

            // minrep
            if ( 1 /*repPtr->spanning*/ ) {

                repPtr->prof = profptr =
                  (PROFILE *) smalloc( sizeof( PROFILE ) );

                profptr->key     = repPtr->id;
                profptr->copynum = repPtr->copynum;
                profptr->patlen  = repPtr->patsize;
                profptr->a = profptr->c = profptr->g = profptr->t = profptr->n =
                  0;
                profptr->proflen = rep.proflen;
                profptr->indices = smalloc( rep.proflen * sizeof( int ) );

                // records counts from sequence
                lentemp = repPtr->lastindex - repPtr->firstindex + 1;

                for ( i = 0; i < lentemp; i++ ) {
                    if ( toupper( repPtr->sequence[i] ) == 'A' )
                        profptr->a++;
                    else if ( toupper( repPtr->sequence[i] ) == 'C' )
                        profptr->c++;
                    else if ( toupper( repPtr->sequence[i] ) == 'G' )
                        profptr->g++;
                    else if ( toupper( repPtr->sequence[i] ) == 'T' )
                        profptr->t++;
                    else
                        profptr->n++;
                }

                profptr->acgtCount =
                  profptr->a + profptr->c + profptr->g + profptr->t;
            }

            pc = (short *) rep.profile;

            for ( i = 0; i < rep.proflen; i++ ) {
                for ( j = 0; j < 5; j++ ) {
                    comps[j] = pc[j];
                }

                // get index;
                profptr->indices[i] = GetCompositionId( comps );
                pc += 6;
            }

            free( rep.subsequence );
            free( rep.pattern );
            free( rep.profile );

            /* <<<- RC PROFILE */
            if ( 1 /*repPtr->spanning*/ ) {

                repPtr->profrc = profptr =
                  (PROFILE *) smalloc( sizeof( PROFILE ) );

                profptr->key = repPtr->id;

                // set member variables
                rep.updates       = 0;
                rep.concensussize = repPtr->patsize;
                rep.firstindex    = repPtr->firstindex;
                rep.lastindex     = repPtr->lastindex;
                rep.pattern       = GetReverseComplement( repPtr->pattern );
                rep.subsequence   = GetReverseComplement( repPtr->sequence );

                if ( NULL == rep.pattern || NULL == rep.subsequence ) {
                    fputs( "Memory allocation failed on allocatin a revese "
                           "complement. Aborting.\n",
                      stderr );
                    return ( 20 );
                }

                rep.profile = NULL;
                rep.proflen = 0;

                // get RC profile
                cRes = ComputeBetterPattern( &rep, 0, -indel, sm );

                if ( -1 == cRes ) {
                    fputs( "Memory allocation failed on "
                           "ComputeBetterPattern(). Aborting.\n",
                      stderr );
                    return ( 21 );
                }

                // fill out RC indices and set other vars
                pc = (short *) rep.profile;

                // fix
                profptr->patlen = rep.concensussize;
                profptr->copynum =
                  (double) ( rep.gappedlength ) / (double) ( rep.proflen );
                profptr->proflen   = rep.proflen;
                profptr->indices   = smalloc( rep.proflen * sizeof( int ) );
                profptr->acgtCount = repPtr->prof->acgtCount;

                for ( i = 0; i < rep.proflen; i++ ) {
                    for ( j = 0; j < 5; j++ ) {
                        comps[j] = pc[j];
                    }

                    // get index;
                    profptr->indices[i] = GetCompositionId( comps );
                    pc += 6;
                }

                // free memory
                free( rep.subsequence );
                free( rep.pattern );
                free( rep.profile );

                repPtr->minRepresentation = MinimumRepresentation(
                  repPtr->prof->indices, repPtr->prof->proflen,
                  repPtr->profrc->indices, repPtr->profrc->proflen,
                  &( repPtr->minrlen ), IDENTICAL_ONLY );

                if ( NULL == repPtr->minRepresentation ) {
                    fputs(
                      "Minimum representation is NULL. Aborting.\n", stderr );
                    return ( 100 );
                }

                repPtr->acgtCount = repPtr->prof->acgtCount;

                /* to limit memory usage */
                if ( !repPtr->spanning ) {
                    free( repPtr->prof->indices );
                    free( repPtr->profrc->indices );
                    free( repPtr->prof );
                    repPtr->prof = NULL;
                    free( repPtr->profrc );
                    repPtr->profrc = NULL;
                }
            }

            /* END FOR SORTING BY MINREPRESENTATION FOR PIPELINE */

            theid++;

            // if ( ( thecount % 10000 ) == 0 )
            //     fprintf( logfp, "TRF output lines processed: %d\n", thecount
            //     );
        }
    }

    fclose( fp );
    // fclose( logfp );

    /* sort becaue PROCLU expects a file sorted on pattern size ASC */
    /* PIPELINE: arsize_and_min_rep_sort is now used to sort so that redundancy
     * program can cleanup later */
    if ( verbose )
        fputs( "done\nsorting... ", stderr );

    EasyListQuickSort( repList, arsize_and_min_rep_sort );

    /* output indexed .dat, leb36 files */
    if ( verbose )
        fprintf( stderr,
          "done (%zu items)\nwriting indexed .dat and leb36 files... ",
          repList->size );

    leb36file = calloc( strlen( outfile_prefix ) + 7, sizeof( *leb36file ) );
    snprintf(
      leb36file, strlen( outfile_prefix ) + 7, "%s.leb36", outfile_prefix );
    fp = fopen( leb36file, "w" );

    if ( fp == NULL ) {
        fputs( "Unable to open leb36 file for writing. Aborting.\n", stderr );
        return ( 17 );
    }

    theid = startid;
    EASY_STRING_HASH *headerHash =
      EasyStringHashCreate( repList->size, NULL, free );
    size_t readCount = 0, ge7ReadCount = 0, TRCount = 0, ge7TRCount = 0;
    for ( tnode = repList->head; tnode != NULL; tnode = tnode->next ) {
        repPtr = (REP_STRUCT *) ( tnode->item );
        /*fprintf(fp,"%d %d %d %d %.1f %d %d %d %d %d %d %d %d %.2f %s
           %s\n",repPtr->id, repPtr->firstindex, repPtr->lastindex,
           repPtr->period, repPtr->copynum, repPtr->patsize, repPtr->matchperc,
                        repPtr->indelperc, repPtr->score, repPtr->acount,
                        repPtr->ccount, repPtr->gcount, repPtr->tcount,
                        repPtr->entropy, repPtr->pattern, repPtr->sequence );*/

        if ( repPtr->spanning ) {
            repPtr->id        = theid;
            repPtr->prof->key = theid;
            theid++;
            /* output index to standard output */
            fprintf( stdout, "%d\t%s\t%d\t%d\t%.1f\t%d\t%s\n", repPtr->id,
              repPtr->header, repPtr->firstindex, repPtr->lastindex,
              repPtr->copynum, repPtr->patsize, repPtr->pattern );
            /* output leb36 to file */
            WriteProfileWithRC( fp, repPtr->prof, repPtr->profrc,
              (char *) repPtr->left, (char *) repPtr->right );
        }

        // We've not seen this header before. Set and increment the
        // ge7ReadCount.
        int *hval;
        if ( ( hval = EasyStringHashGet( headerHash, repPtr->header ) ) ==
             NULL ) {
            readCount++;
            ge7ReadCount +=
              ( ( repPtr->patsize >= 7 ) && ( repPtr->spanning ) );
            int *flag = malloc( sizeof( int ) );
            *flag     = ( ( repPtr->patsize >= 7 ) && ( repPtr->spanning ) );
            EasyStringHashSet( headerHash, repPtr->header, flag );
        }
        // We've seen this header before. Only increment ge7ReadCount
        // if value in hash was not already set to 1.
        else {
        // We’ve seen this header before. Only increment ge7ReadCount
        // if value in hash was not already set to 1.
        // Fixed by G. Benson 11.12.22 to set hval correctly, i.e.,
        //       repeat has patternsize ≥ 7
        //       and is spanning
        // Then use hval to increment ge7ReadCount
            if ( *hval == 0 ) {
                *hval = ( ( repPtr->patsize >= 7 ) && ( repPtr->spanning ) );
                ge7ReadCount += *hval;
            }

            //ge7ReadCount += ( *hval == 0 );
            //*hval = 1;
        }

        TRCount++;
        // Increment the ge7TRCount if TR pattern >= 7 and spanning
        ge7TRCount += ( ( repPtr->patsize >= 7 ) && ( repPtr->spanning ) );

        // fprintf( fph, "%d\t%s\t%d\t%d\t%.1f\t%d\t%s\n", repPtr->id,
        //   repPtr->header, repPtr->firstindex, repPtr->lastindex,
        //   repPtr->copynum, repPtr->patsize, repPtr->pattern );
    }

    fclose( fp );

    /* this is a file summarizing all results, including nonspanning and
     * all patsizes, needed for statistics at the end */
    indexfileh = calloc( strlen( outfile_prefix ) + 11, sizeof( *indexfileh ) );
    snprintf( indexfileh, strlen( outfile_prefix ) + 11, "%s.indexhist",
      outfile_prefix );
    fph = fopen( indexfileh, "w" );

    if ( fph == NULL ) {
        fputs(
          "Unable to open indexhist file for writing. Aborting.\n", stderr );
        return ( 18 );
    }
    fprintf( fph, "%zu\t%zu\t%zu\t%zu\n", ge7TRCount, TRCount, ge7ReadCount,
      readCount );
    fclose( fph );

    // /* output leb36 to standard output */
    // if ( verbose )
    //     fputs( "done\ngenerating profiles and writing leb36 file... ", stderr
    //     );

    // for ( tnode = repList->head; tnode != NULL; tnode = tnode->next ) {
    //     repPtr = (REP_STRUCT *) ( tnode->item );

    //     /* only output spanning ones */
    //     if ( 0 == repPtr->spanning )
    //         continue;

    //     // write profile to standard output
    //     WriteProfileWithRC( stdout, repPtr->prof, repPtr->profrc,
    //       (char *) repPtr->left, (char *) repPtr->right );
    // }

    // free data
    free( sm );
    EasyListDestroy( repList );

    if ( verbose )
        fputs( "done\n", stderr );

    /* return last used id */
    // return (thecount >= 1) ? (theid - 1) : (-2);
    return ( thecount >= 1 ) ? 0 : 1;
}
