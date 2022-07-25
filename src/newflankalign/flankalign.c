/****************************************************************************************************

   flankalign.c - creates .MAP files (using narrowband alignment)

   Usage: ./flankalign outdir xmldir inputfile MAXERRORS MAXTHREADS
PATLEN_SIZE_ERR

****************************************************************************************************/

#include <stdio.h>

#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#define max( a, b ) ( ( ( a ) >= ( b ) ) ? ( a ) : ( b ) )
#define min( a, b ) ( ( ( a ) <= ( b ) ) ? ( a ) : ( b ) )

#include "../libs/easylife/easylife.h"

//#include "narrowbandDistanceAlignment.h"
#include "bitwise edit distance alignment multiple word no end penalty.h"

//#include "aln.h"
//#include "aln_nb.h"

//#define PER_READ_STATS
//#define PRINT_ALIGNMENTS
int LARGESERRORALLOWED = 0;
int MAXEL1             = 0;
int MAXER1             = 0;
int MAXEL2             = 0;
int MAXER2             = 0;

#define MATCH_SCORE ( 1 )
#define MISM_PEN ( -3 )
#define GAP_OPEN_PEN ( -11 )
#define GAP_EXT_PEN ( -4 )

int MAXERRORS  = ( 0 );
int MAXTHREADS = ( 1 );

#define MAXLINESIZE ( 10000000 )

int PATLEN_SIZE_ERR      = 0;
int MAX_FLANK_CONSIDERED = 1000;

char *Complementascii =
  NULL; /* precomputed array used for the reverse complement function */
int ThreadCounter = 0;

/********************************       init_complement_ascii
 * ****************************/
char *init_complement_ascii( void )

{
    /* complement has 256 entries so that finding the entries for A, C, G and T
     * which are complement ascii values */
    /* require no calculation */
    int   i;
    char *ca = (char *) calloc( 256, sizeof( int ) );

    if ( ca == NULL ) {
        return NULL;
    }

    for ( i = 0; i < 256; i++ )
        ca[i] = i;

    ca['A'] = 'T';
    ca['C'] = 'G';
    ca['G'] = 'C';
    ca['T'] = 'A';
    ca['a'] = 't';
    ca['c'] = 'g';
    ca['g'] = 'c';
    ca['t'] = 'a';

    return ca;
}

/*******************************************************************************************/
char *GetReverse( char *original ) {
    char *sourceptr, *destinptr;
    int   length;
    char *buffer;

    /* find out how long the string is */
    length = strlen( original );

    /* allocate the new string */
    buffer = (char *) malloc( sizeof( char ) * ( length + 1 ) );
    if ( buffer == NULL ) {
        return NULL;
    }

    /* reverse pattern */
    sourceptr  = original;
    destinptr  = buffer + length; /* position at termination */
    *destinptr = '\0';            /* terminate string */
    destinptr--;
    while ( *sourceptr != '\0' ) {
        ( *destinptr ) = *sourceptr;
        destinptr--;
        sourceptr++;
    }

    /* return buffer containg reverse complement */
    return buffer;
}

/***********************************************************************/
char *GetComplement( char *original ) {
    char *sourceptr, *destinptr;
    int   length;
    char *buffer;

    /* safety check */
    if ( Complementascii == NULL ) {
        Complementascii = init_complement_ascii();
        if ( NULL == Complementascii )
            return NULL;
    }

    /* find out how long the string is */
    length = strlen( original );

    /* allocate the new string */
    buffer = (char *) malloc( sizeof( char ) * ( length + 1 ) );
    if ( buffer == NULL ) {
        return NULL;
    }

    /* complement pattern */
    sourceptr = original;
    destinptr = buffer;
    while ( *sourceptr != '\0' ) {
        ( *destinptr ) = Complementascii[( *sourceptr )];
        destinptr++;
        sourceptr++;
    }
    *destinptr = '\0'; /* terminate string */

    /* return buffer containg reverse complement */
    return buffer;
}

/*******************************************************************************************/
void doCriticalErrorAndQuit( const char *format, ... ) {

    va_list argp;

    printf( "\n\nERROR!!!!:\n" );
    if ( format == NULL )
        exit( 0 );
    ;

    va_start( argp, format );
    vprintf( format, argp );
    fflush( stdout );
    va_end( argp );

    exit( 1 );
}

typedef struct {
    long long int id;
    int           patsize, leftlen, rightlen;
    char *        left;
    char *        right;
    char *        leftcp;
    char *        rightcp;
} FLANK;

/*******************************************************************************************/
void flankDestroy( void *in ) {

    FLANK *g = (FLANK *) in;

    free( g->left );
    free( g->right );
    free( g->leftcp );
    free( g->rightcp );
    free( g );
}

/*******************************************************************************************/
int __patsizeCmp( const void *item1, const void *item2 ) {

    int d1, d2;

    d1 = ( (FLANK *) item1 )->patsize;
    d2 = ( (FLANK *) item2 )->patsize;

    if ( d1 > d2 )
        return 1;
    else if ( d1 < d2 )
        return -1;

    return 0;
}

/*******************************************************************************************/
int main( int argc, char *argv[] ) {

    char *outdir, *xmldir, *filename, *buffer, *src1, *src2, *src3, *src4,
      *src5, *dst, tempbuf[2000];
    FILE *        rd_fp;
    long long int rid;
    long long int processed = 0;
    char *        hasdata   = NULL;
    int           a1, a2, status, i, *PARR;

    /* cmd arguments */
    if ( argc < 8 )
        doCriticalErrorAndQuit(
          "\n\nFlankAlign - Please use ./flankalign  outdir xmldir inputfile "
          "MAXERRORS MAX_FLANK_CONSIDERED MAXTHREADS PATLEN_SIZE_ERR\n\n" );

    outdir   = argv[1];
    xmldir   = argv[2]; // never used
    filename = argv[3];

    if ( argc >= 5 ) {
        MAXERRORS = atoi( argv[4] );
    }

    if ( argc >= 6 ) {
        MAX_FLANK_CONSIDERED = atoi( argv[5] );
    }

    if ( MAX_FLANK_CONSIDERED < 1 ) {
        doCriticalErrorAndQuit(
          "\n\nFlankAlign - ? should be more than or equal to one!\n\n" );
    }

    if ( argc >= 7 ) {
        MAXTHREADS = atoi( argv[6] );
    }

    if ( MAXTHREADS < 1 ) {
        doCriticalErrorAndQuit( "\n\nFlankAlign - maxthreads should be more "
                                "than or equal to one!\n\n" );
    }

    if ( argc >= 8 ) {
        PATLEN_SIZE_ERR = atoi( argv[7] );
    }

    if ( PATLEN_SIZE_ERR < 1 ) {
        doCriticalErrorAndQuit(
          "\n\nFlankAlign - PATLEN_SIZE_ERR must be at least 1%%!\n\n" );
    }

    PARR = (int *) smalloc( sizeof( int ) * MAXTHREADS );

    /* alignment options */
    // if (NULL==(sm=init_sm(MATCH_SCORE,MISM_PEN)))
    //   doCriticalErrorAndQuit("Memory error. Aborting!\n\n");

    if ( Complementascii == NULL ) {
        Complementascii = init_complement_ascii();
        if ( NULL == Complementascii )
            return 1;
    }

    /* enter out dir */
    if ( chdir( outdir ) ) { // should this be negated?
        mkdir( outdir, 0755 );

        if ( chdir( outdir ) ) {

            doCriticalErrorAndQuit(
              "\n\nFlankAlign - could not enter or create output dir!\n\n" );
        }
    }

    /* parse file one cluster at a time */
    buffer = smalloc( sizeof( char ) * ( MAXLINESIZE + 1 ) );

    if ( filename[0] == '-' )
        rd_fp = stdin;
    else
        rd_fp = fopen( filename, "r" );
    if ( !rd_fp )
        doCriticalErrorAndQuit(
          "\n\nFlankAlign - Could not open input file. Aborting!\n\n" );

    while ( ( hasdata = fgets( buffer, MAXLINESIZE, rd_fp ) ) != NULL ) {

        if ( buffer[0] == '@' ) {

            /* read cluster info */
            a1 = a2 = 0;
            sscanf( buffer, "@(%d_%d)", &a1, &a2 );
            if ( a1 <= 0 || a2 <= 0 ) {
                printf( "Can't read cluster number." );
                exit( 1 );
            }
            break;
        }
    };

    if ( hasdata )
        while ( 1 ) {

            EASY_LIST *ref_list, *read_list;
            EASY_NODE *nd1, *nd2, *windowstart;
            int        start, end, length, i;
            FLANK *    templ_flank_ptr;

            ref_list  = EasyListCreate( NULL, flankDestroy );
            read_list = EasyListCreate( NULL, flankDestroy );

            /* read cluster contents one line at a time */
            while ( 1 ) {

                rid = strtol( buffer, NULL, 10 );

                if ( rid < 0 ) {

                    // printf("rid: %ld\n",rid);

                    src1 = buffer;
                    while ( *src1 != ',' && *src1 != '\0' ) {
                        src1++;
                    }
                    if ( *src1 == ',' ) {
                        *src1 = '\0';
                        src1++;
                    }

                    src2 = src1;
                    while ( *src2 != ',' && *src2 != '\0' ) {
                        *src2 = toupper( *src2 );
                        if ( NULL == strchr( "ACGT", *src2 ) ) {
                            *src2 = 'N';
                        }
                        src2++;
                    }
                    if ( *src2 == ',' ) {
                        *src2 = '\0';
                        src2++;
                    }

                    src3 = src2;
                    while ( *src3 != ',' && *src3 != '\0' ) {
                        src3++;
                    }
                    if ( *src3 == ',' ) {
                        *src3 = '\0';
                        src3++;
                    }

                    src4 = src3;
                    while ( *src4 != ',' && *src4 != '\0' ) {
                        *src4 = toupper( *src4 );
                        if ( NULL == strchr( "ACGT", *src4 ) ) {
                            *src4 = 'N';
                        }
                        src4++;
                    }
                    if ( *src4 == ',' ) {
                        *src4 = '\0';
                        src4++;
                    }

                    src5 = src4 + strlen( src4 ) - 1;
                    while ( *src5 == '\n' || *src5 == '\r' ) {
                        *src5 = '\0';
                        src5--;
                    }

                    templ_flank_ptr = (FLANK *) smalloc( sizeof( FLANK ) );

                    templ_flank_ptr->left     = GetReverse( src1 );
                    templ_flank_ptr->right    = strdup( src3 );
                    templ_flank_ptr->leftlen  = strlen( src1 );
                    templ_flank_ptr->rightlen = strlen( src3 );
                    templ_flank_ptr->leftcp   = NULL;
                    templ_flank_ptr->rightcp  = NULL;

                    if ( !templ_flank_ptr->left || !templ_flank_ptr->right )
                        doCriticalErrorAndQuit(
                          "\n\nFlankAlign - memory error 1. Aborting!\n\n" );

                    templ_flank_ptr->patsize = strlen( src4 );
                    templ_flank_ptr->id      = rid;

                    EasyListInsertTail( ref_list, templ_flank_ptr );

                } else if ( rid > 0 ) {

                    // printf("rid: %ld\n",rid);

                    src1 = buffer;
                    while ( *src1 != ',' && *src1 != '\0' )
                        src1++;
                    if ( *src1 == ',' ) {
                        *src1 = '\0';
                        src1++;
                    }

                    src2 = src1;
                    while ( *src2 != ',' && *src2 != '\0' )
                        src2++;
                    if ( *src2 == ',' ) {
                        *src2 = '\0';
                        src2++;
                    }

                    src3 = src2;
                    while ( *src3 != ',' && *src3 != '\0' ) {
                        src3++;
                    }
                    if ( *src3 == ',' ) {
                        *src3 = '\0';
                        src3++;
                    }

                    src4 = src3;
                    while ( *src4 != ',' && *src4 != '\0' ) {
                        *src4 = toupper( *src4 );
                        if ( NULL == strchr( "ACGT", *src4 ) ) {
                            *src4 = 'N';
                        }
                        src4++;
                    }
                    if ( *src4 == ',' ) {
                        *src4 = '\0';
                        src4++;
                    }

                    src5 = src4 + strlen( src4 ) - 1;
                    while ( *src5 == '\n' || *src5 == '\r' ) {
                        *src5 = '\0';
                        src5--;
                    }

                    templ_flank_ptr = (FLANK *) smalloc( sizeof( FLANK ) );

                    start = strtol( src1, NULL, 10 );
                    if ( ( errno == ERANGE ) || ( errno != 0 && start == 0 ) ) {
                        perror( "strtol" );
                        exit( EXIT_FAILURE );
                    }

                    end = strtol( src2, NULL, 10 );
                    if ( ( errno == ERANGE ) || ( errno != 0 && end == 0 ) ) {
                        perror( "strtol" );
                        exit( EXIT_FAILURE );
                    }

                    length = strlen( src3 );

                    // fprintf(stderr, "start: %d end: %d length: %d\n\n",start,
                    // end, length);
                    // exit(1);

                    templ_flank_ptr->left =
                      calloc( start + 2, sizeof( *templ_flank_ptr->left ) );
                    templ_flank_ptr->right = calloc(
                      ( length - end ) + 2, sizeof( *templ_flank_ptr->right ) );

                    if ( !templ_flank_ptr->left )
                        doCriticalErrorAndQuit(
                          "\n\nFlankAlign - memory error 2 in left flank. rid: "
                          "%d, start: %d end: %d length: %d, start+2: %d, "
                          "length-end+2: %d. Aborting!\n\n",
                          rid, start, end, length, start + 2,
                          length - end + 2 );

                    if ( !templ_flank_ptr->right )
                        doCriticalErrorAndQuit(
                          "\n\nFlankAlign - memory error 2 in right flank. "
                          "rid: %d, start: %d end: %d length: %d, start+2: %d, "
                          "length-end+2: %d. Aborting!\n\n",
                          rid, start, end, length, start + 2,
                          length - end + 2 );

                    dst = templ_flank_ptr->left;
                    for ( i = start - 2; i >= 0; i-- ) {
                        *dst = src3[i];
                        dst++;
                    }
                    *dst = 0;

                    dst = templ_flank_ptr->right;
                    for ( i = end; i < length; i++ ) {
                        *dst = src3[i];
                        dst++;
                    }
                    *dst = 0;

                    templ_flank_ptr->patsize = strlen( src4 );
                    templ_flank_ptr->id      = rid;

                    /* FOR FLANKSHORT */
                    templ_flank_ptr->leftlen = strlen( templ_flank_ptr->left );
                    templ_flank_ptr->rightlen =
                      strlen( templ_flank_ptr->right );

                    if ( templ_flank_ptr->leftlen > MAX_FLANK_CONSIDERED )
                        templ_flank_ptr->left[MAX_FLANK_CONSIDERED] = '\0';
                    if ( templ_flank_ptr->rightlen > MAX_FLANK_CONSIDERED )
                        templ_flank_ptr->right[MAX_FLANK_CONSIDERED] = '\0';
                    /* ENDO OF FLANKSHORT */

                    templ_flank_ptr->leftlen = strlen( templ_flank_ptr->left );
                    templ_flank_ptr->rightlen =
                      strlen( templ_flank_ptr->right );

                    templ_flank_ptr->leftcp =
                      GetComplement( templ_flank_ptr->left );
                    templ_flank_ptr->rightcp =
                      GetComplement( templ_flank_ptr->right );

                    EasyListInsertTail( read_list, templ_flank_ptr );
                }

                if ( ( hasdata = fgets( buffer, MAXLINESIZE, rd_fp ) ) ==
                       NULL ||
                     buffer[0] == '@' )
                    break;
            }

            /* process cluster*/
            processed++;
            ThreadCounter++;

            // if (processed>=10) exit(0);

            pid_t pID = fork();

            if ( pID == 0 ) { // child
                // Code only executed by child process

                // BASICALIGNPAIR     *ba;
                FILE * fp;
                FLANK *read_ptr, *ref_ptr;
                char   tempname[100];
                int    comma;

                /* Redirect standard files to /dev/null */
                // freopen( "/dev/null", "r", stdin);
                // freopen( "/dev/null", "w", stdout);
                // freopen( "/dev/null", "w", stderr);

                EasyListQuickSort( read_list, __patsizeCmp );
                EasyListQuickSort( ref_list, __patsizeCmp );

                sprintf( tempname, "%d_%d.map", a1, a2 );
                fp = fopen( tempname, "w" );

                if ( !fp )
                    doCriticalErrorAndQuit( "\n\nFlankAlign - can't open "
                                            "output file. Aborting!\n\n" );

                fprintf(
                  fp, "@%d_%d	1	1	1	1	1", a1, a2 );

                windowstart = ref_list->head;
                for ( nd1 = read_list->head; nd1 != NULL; nd1 = nd1->next ) {

                    read_ptr = (FLANK *) EasyListItem( nd1 );
                    fprintf( fp,
                      "\n%lld	'	%d	%d	1	1	"
                      "1	",
                      read_ptr->id, read_ptr->leftlen, read_ptr->rightlen );

                    comma = 0;
                    for ( nd2 = windowstart; nd2 != NULL; nd2 = nd2->next ) {

                        int lerr, rerr;
                        int lerr1, rerr1;
                        int lerr2, rerr2;

                        int errors, errors1, errors2;
                        int sumerr1, sumerr2, refleft, refright;

                        double sizeerror;

                        ref_ptr = (FLANK *) EasyListItem( nd2 );

                        /* if sizes are too different, stop */
                        if ( ref_ptr->patsize > read_ptr->patsize ) {
                            sizeerror =
                              ( ref_ptr->patsize / (double) read_ptr->patsize -
                                1.0 ) *
                              100;
                        } else {
                            sizeerror =
                              ( read_ptr->patsize / (double) ref_ptr->patsize -
                                1.0 ) *
                              100;
                        }

                        if ( sizeerror > PATLEN_SIZE_ERR ) {

                            if ( ref_ptr->patsize < read_ptr->patsize ) {
                                windowstart = nd2;
                                continue;
                            } else {
                                break;
                            }
                        }

                        /*****************************************************************************************************************/

                        /* if 0 passed for maxerror, calculate based on flank
                         * length */
                        if ( 0 != MAXERRORS ) {
                            MAXEL1             = MAXERRORS;
                            MAXER1             = MAXERRORS;
                            MAXEL2             = MAXERRORS;
                            MAXER2             = MAXERRORS;
                            LARGESERRORALLOWED = MAXERRORS;
                        }

                        /* shorten reference from 1000 to slightly more than
                         * readlen */
                        if ( 0 == MAXERRORS ) {
                            MAXEL1 =
                              min( 8, (int) ( 0.4 * read_ptr->leftlen + .01 ) );
                            MAXER1 = min(
                              8, (int) ( 0.4 * read_ptr->rightlen + .01 ) );
                            LARGESERRORALLOWED = max( MAXEL1, MAXER1 );
                        }
                        refleft  = min( ref_ptr->leftlen,
                          read_ptr->leftlen + LARGESERRORALLOWED + 2 );
                        refright = min( ref_ptr->rightlen,
                          read_ptr->rightlen + LARGESERRORALLOWED + 2 );

                        /* If ref is from ends of chromosome it could be
                         * shorter. Thats why readlen will be shortened to
                         * reflen */
                        lerr1 = Edit_Distance_multiple_word_NoEndPenaltySeq1(
                          ref_ptr->left, read_ptr->left, refleft,
                          ( refleft > read_ptr->leftlen ) ? read_ptr->leftlen
                                                          : refleft );
                        rerr1 = Edit_Distance_multiple_word_NoEndPenaltySeq1(
                          ref_ptr->right, read_ptr->right, refright,
                          ( refright > read_ptr->rightlen ) ? read_ptr->rightlen
                                                            : refright );

                        if ( lerr1 < 0 || rerr1 < 0 )
                            doCriticalErrorAndQuit(
                              "\n\nFlankAlign - wrong sizes for narrowband "
                              "alignment, seq1 must be larger. Aborting!\n\n" );

                        /* try aligning to the compliment of the other flank
                         * instead */
                        /*****************************************************************************************************************/

                        /* shorten reference from 1000 to slightly more than
                         * readlen */
                        if ( 0 == MAXERRORS ) {
                            MAXEL2 = min(
                              8, (int) ( 0.4 * read_ptr->rightlen + .01 ) );
                            MAXER2 =
                              min( 8, (int) ( 0.4 * read_ptr->leftlen + .01 ) );
                            LARGESERRORALLOWED = max( MAXEL2, MAXER2 );
                        }
                        refleft  = min( ref_ptr->leftlen,
                          read_ptr->rightlen + LARGESERRORALLOWED + 2 );
                        refright = min( ref_ptr->rightlen,
                          read_ptr->leftlen + LARGESERRORALLOWED + 2 );

                        /* If ref is from ends of chromosome it could be
                         * shorter. Thats why readlen will be shortened to
                         * reflen */
                        lerr2 = Edit_Distance_multiple_word_NoEndPenaltySeq1(
                          ref_ptr->left, read_ptr->rightcp, refleft,
                          ( refleft > read_ptr->rightlen ) ? read_ptr->rightlen
                                                           : refleft );
                        rerr2 = Edit_Distance_multiple_word_NoEndPenaltySeq1(
                          ref_ptr->right, read_ptr->leftcp, refright,
                          ( refright > read_ptr->leftlen ) ? read_ptr->leftlen
                                                           : refright );

                        if ( lerr2 < 0 || rerr2 < 0 )
                            doCriticalErrorAndQuit(
                              "\n\nFlankAlign - wrong sizes for narrowband "
                              "alignment, seq1 must be larger. Aborting!\n\n" );

                        // printf("%d - %d - %s - %s - %d\n",ref_ptr->leftlen,
                        // read_ptr->leftlen,ref_ptr->left, read_ptr->left,
                        // lerr); printf("%d - %d - %s - %s -
                        // %d\n",ref_ptr->rightlen,
                        // read_ptr->rightlen,ref_ptr->right, read_ptr->right,
                        // rerr);

                        /*****************************************************************************************************************/
                        /*****************************************************************************************************************/
                        /*****************************************************************************************************************/

                        /* if passes criteria in EITHER direction */
                        if ( ( lerr1 <= MAXEL1 && rerr1 <= MAXER1 ) ||
                             ( lerr2 <= MAXEL2 && rerr2 <= MAXER2 ) ) {

                            /* if passes criteria in BOTH direction */
                            if ( ( lerr1 <= MAXEL1 && rerr1 <= MAXER1 ) &&
                                 ( lerr2 <= MAXEL2 && rerr2 <= MAXER2 ) ) {

                                /* print #errors with smallest sum */
                                sumerr1 = lerr1 + rerr1;
                                sumerr2 = lerr2 + rerr2;

                                if ( sumerr1 <= sumerr2 ) {
                                    lerr = lerr1;
                                    rerr = rerr1;
                                } else {
                                    lerr = lerr2;
                                    rerr = rerr2;
                                }

                                /* if passes criteria in FORWARD direction */
                            } else if ( lerr1 <= MAXEL1 && rerr1 <= MAXER1 ) {

                                lerr = lerr1;
                                rerr = rerr1;

                                /* if passes criteria in REVERSE direction */
                            } else {

                                lerr = lerr2;
                                rerr = rerr2;
                            }

                            /* print comma separated refs with errors */
                            if ( 1 == comma ) {
                                fprintf( fp, "," );
                            }
                            fprintf(
                              fp, "%lld:%d:%d", ref_ptr->id, lerr, rerr );
                            comma = 1;
                        }
                    }
                }

                fclose( fp );

                _exit( 0 );

            } else if ( pID < 0 ) { // failed to fork
                doCriticalErrorAndQuit(
                  "\n\nFlankAlign - failed to fork code=%d (%s). Aborting!\n\n",
                  errno, strerror( errno ) );
                exit( 1 );
            }

            PARR[ThreadCounter - 1] = pID;

            // Code below only executed by parent process
            //printf( "processing: %d_%d (pid: %d) threadcounter: %d\n", a1, a2,
            //  pID, ThreadCounter );

            /* child has its own copy now, so it's ok to destroy for parent */
            EasyListDestroy( ref_list );
            EasyListDestroy( read_list );

            /* if full batch scheduled, wait for one to finish and schedule
             * another */
            if ( ThreadCounter >= MAXTHREADS ) {

                pid_t retid = waitpid( PARR[0], &status, 0 );

                if ( -1 == retid )
                    doCriticalErrorAndQuit(
                      "\n\nFlankAlign - wait function returned -1. Aborting!" );
                //else
                //    printf( "\tterminated: %d\n", retid );

                for ( i = 1; i < ThreadCounter; i++ ) {
                    PARR[i - 1] = PARR[i];
                }

                ThreadCounter--;
            }

            /* end of file? */
            if ( !hasdata )
                break;

            /* read cluster info */
            a1 = a2 = 0;
            sscanf( buffer, "@(%d_%d)", &a1, &a2 );
            if ( a1 <= 0 || a2 <= 0 ) {
                printf( "Can't read cluster number. \n%s\n", buffer );
                exit( 1 );
            }
        }

    /* wait for all child processes to finish */
    for ( i = 0; i < ThreadCounter; i++ ) {

        pid_t retid = waitpid( PARR[i], &status, 0 );

        if ( -1 == retid )
            doCriticalErrorAndQuit(
              "\n\nFlankAlign - wait function returned -1. Aborting!" );
        //else
        //    printf( "\tterminated: %d\n", retid );
    }
    ThreadCounter = 0;

    /* done */
    fclose( rd_fp );
    printf( "Processed: %lld\n", processed);

    return 0;
}
