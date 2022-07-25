#include <dirent.h>
#include <malloc.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/types.h>
#include <time.h>

#include "../libs/easylife/easylife.h"
#include "ghash.h"

#define max3( a, b, c )                                           \
    ( ( ( a ) >= ( b ) ) ? ( ( ( a ) >= ( c ) ) ? ( a ) : ( c ) ) \
                         : ( ( ( b ) >= ( c ) ) ? ( b ) : ( c ) ) )
#define min3( a, b, c )                                           \
    ( ( ( a ) <= ( b ) ) ? ( ( ( a ) <= ( c ) ) ? ( a ) : ( c ) ) \
                         : ( ( ( b ) <= ( c ) ) ? ( b ) : ( c ) ) )
#define max( a, b ) ( ( ( a ) >= ( b ) ) ? ( a ) : ( b ) )
#define min( a, b ) ( ( ( a ) <= ( b ) ) ? ( a ) : ( b ) )

typedef struct _rlink {
    int cid, refs, reads;

    int * idlist; // reads
    char *dirlist;

    int * negidlist; // refs
    char *negdirlist;

    struct _rlink *next;
} RLINK;

int LoadClustersFromFile( char *filename, RLINK *citem );

int     clustercount = 0, totclusters = 0, finclusters = 0;
GSHASH *HASHSEEN = NULL;

int TotalRefs         = 0;
int TotalRefsWritten  = 0;
int TotalReads        = 0;
int TotalReadsWritten = 0;

/*******************************************************************************************/
int name_cmp( const void *item1, const void *item2 ) {

    return strcmp( (char *) item1, (char *) item2 );
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

/*******************************************************************************************/
int *intdup( int a ) {

    int *b = (int *) scalloc( 1, sizeof( int ) );

    if ( b == NULL ) {
        printf( "\nERROR: Insuficient memory to create an integer\n\n" );
        exit( 1 );
    }

    *b = a;

    return b;
}

/*******************************************************************************************/

typedef struct _rd {
    int  id;
    char dir;
} RD;

RD *RDCreate( int a, int dup ) {

    RD *rdptr = (RD *) scalloc( 1, sizeof( RD ) );

    if ( rdptr == NULL ) {
        printf( "\nERROR: Insuficient memory to create RD structure.\n\n" );
        exit( 1 );
    }

    rdptr->id  = a;
    rdptr->dir = dup ? '\"' : '\'';

    return rdptr;
}

/*******************************************************************************************/
int intcmp( const void *item1, const void *item2 ) {

    int a, b;

    a = *(int *) item1;
    b = *(int *) item2;

    if ( a > b )
        return 1;
    else if ( a < b )
        return -1;

    return 0;
}

/*******************************************************************************************/
int main( int argc, char **argv ) {
    FILE *fpo;
    char  bigtempbuf[2000];
    int   k, i, j, filescreated = 1, ccount, comma, flip, flipset, depth = 0,
                 maxdepth = 0;
    time_t         startTime;
    struct dirent *de = NULL;
    DIR *          d  = NULL;
    long long int  nwritten, nread;
    RLINK *        CLIST = NULL, *clist2 = NULL, *temp = NULL, *tmlookup = NULL,
          *temphead = NULL, *LASTMERGED = NULL;
    EASY_NODE *fl1;
    EASY_LIST *FILE_LIST = EasyListCreate( NULL, NULL );

    // verify parameter count
    if ( argc < 3 ) {
        printf(
          "\njoin_clusters v1.00 - joins all .CLU files in a directory based "
          "on negative indices found in a cluster, assigns unique ids." );
        printf( "\nauthors: Yevgeniy Gefland" );
        printf( "\n\nUsage:\n\n%s INPUTFOLDER outputfile\n", argv[0] );
        exit( 1 );
    }

    // open out file
    fpo = fopen( argv[2], "w" );

    if ( NULL == fpo ) {
        printf( "\n\nError, cannot open file '%s' for writing!\n", argv[2] );
        exit( 1 );
    }

    // read file(s)
    d = opendir( argv[1] );

    if ( d != NULL ) {

        while ( ( de = readdir( d ) ) != NULL ) {

            if ( strlen( de->d_name ) > 6 &&
                 0 == strcasecmp(
                        ".clu", de->d_name + strlen( de->d_name ) - 4 ) ) {

                bigtempbuf[0] = 0;
                strcat( bigtempbuf, argv[1] );

                if ( '/' != bigtempbuf[strlen( bigtempbuf ) - 1] )
                    strcat( bigtempbuf, "/" );

                strcat( bigtempbuf, de->d_name );

                //printf( "\nreading %s...", bigtempbuf );
                ccount = LoadClustersFromFile( bigtempbuf, NULL );
                //printf( "%d clusters\n", ccount );
                totclusters += ccount;
            }
        }

        closedir( d );
    }

    // allocate cluster stucture
    if ( 0 == totclusters )
        return 0;

    CLIST = (RLINK *) calloc( totclusters, sizeof( RLINK ) );

    if ( NULL == CLIST ) {
        printf( "\n\nCould not allocate CLIST. Aborting!\n" );
        exit( 1 );
    }

    // hashseen
    HASHSEEN = CreateSingleHash( 2 * totclusters );

    if ( NULL == HASHSEEN ) {
        printf( "\n\nCould not allocate CLIST. Aborting!\n" );
        exit( 1 );
    }

    // load clusters(s)
    d            = opendir( argv[1] );
    clustercount = 0;

    if ( d != NULL ) {

        while ( ( de = readdir( d ) ) != NULL ) {

            if ( strlen( de->d_name ) > 6 &&
                 0 == strcasecmp(
                        ".clu", de->d_name + strlen( de->d_name ) - 4 ) ) {

                bigtempbuf[0] = 0;
                strcat( bigtempbuf, argv[1] );

                if ( '/' != bigtempbuf[strlen( bigtempbuf ) - 1] )
                    strcat( bigtempbuf, "/" );

                strcat( bigtempbuf, de->d_name );

                EasyListInsertHead( FILE_LIST, strdup( bigtempbuf ) );
            }
        }

        closedir( d );
    }

    /* sort to make deterministic */
    EasyListQuickSort( FILE_LIST, name_cmp );

    /* load */
    for ( fl1 = FILE_LIST->head; fl1 != NULL; fl1 = fl1->next ) {
        char *strtmp = (char *) EasyListItem( fl1 );
        //printf( "\nprocessing %s...", strtmp );
        ccount = LoadClustersFromFile( strtmp, CLIST );
        //printf( "%d clusters\n", ccount );
    }

    /* join clusters CHANGED Nov 10, 2010 because of directionality problem
     * detected in pipeline */
    finclusters = totclusters;

    for ( i = 0; i < totclusters; i++ ) {

        //printf( "\n%d %.2lf%%", i, i * 100 / (double) totclusters );

        LASTMERGED = NULL;

        for ( j = 0; j < CLIST[i].refs; j++ ) {

            // printf(" %d",CLIST[i].negidlist[j]);

            flip = flipset = 0;

            if ( NULL != ( clist2 = GetSingleHashItem(
                             HASHSEEN, CLIST[i].negidlist[j] ) ) ) {

                if ( 0 != clist2->cid &&
                     clist2 !=
                       ( CLIST + i ) ) { // copy unmarked only and not self

                    if ( NULL == LASTMERGED )
                        CLIST[i].next = clist2;
                    else
                        LASTMERGED->next = clist2;

                    clist2->cid = 0; // mark
                    finclusters--;

                    // change hashseen pointers of all linked lists as well
                    for ( temp = clist2; temp != NULL; temp = temp->next ) {

                        LASTMERGED = temp;

                        for ( k = 0; k < temp->refs; k++ ) {
                            SetSingleHashItem(
                              HASHSEEN, temp->negidlist[k], CLIST + i );

                            // flip?
                            if ( !flipset &&
                                 temp->negidlist[k] == CLIST[i].negidlist[j] ) {
                                flipset = 1;
                                flip    = ( temp->negdirlist[k] !=
                                         CLIST[i].negdirlist[j] );

                                //if ( flip ) {
                                //    printf( "\nflipped c:%d (%d)!", i,
                                //      CLIST[i].negidlist[j] );
                                //}
                            }
                        }
                    }

                    // flip if needed
                    if ( flip ) {
                        for ( temp = clist2; temp != NULL; temp = temp->next ) {
                            for ( k = 0; k < temp->refs; k++ ) {
                                temp->negdirlist[k] =
                                  ( temp->negdirlist[k] == '\'' ) ? '\"' : '\'';
                            }

                            for ( k = 0; k < temp->reads; k++ ) {
                                temp->dirlist[k] =
                                  ( temp->dirlist[k] == '\'' ) ? '\"' : '\'';
                            }
                        }
                    }
                }
            }

            SetSingleHashItem( HASHSEEN, CLIST[i].negidlist[j], CLIST + i );
        }
    }

    // printf("\nGot here!");
    // exit(1);

    /* write output and double CHECK for dups */
    printf( "Writing final cluster file...\n" );
    DestroySingleHash( HASHSEEN );
    HASHSEEN = CreateSingleHash( 2 * totclusters );

    for ( i = 0; i < totclusters; i++ ) {
        if ( 0 != CLIST[i].cid ) {
            // fprintf(fpo,"@");

            temphead = CLIST + i;
            comma    = 0;

            // write refs, while skipping dups in negidlists
            for ( temp = CLIST + i; temp != NULL; temp = temp->next ) {

                for ( k = 0; k < temp->refs; k++ ) {

                    if ( NULL != ( tmlookup = GetSingleHashItem(
                                     HASHSEEN, temp->negidlist[k] ) ) ) {

                        // we got a duplicate but it's already in the same
                        // cluster so it's ok
                        if ( temphead == tmlookup ) {

                            // duplicate in diff clustrs? that shouldn't happen
                            // at this point
                        } else {
                            printf( "\n\nError: duplicate ids found in "
                                    "clusters! (%d)\n",
                              temp->negidlist[k] );
                            // exit(1);
                        }
                    }

                    if ( temphead != tmlookup ) {
                        if ( comma ) {
                            fprintf( fpo, "," );
                        }

                        fprintf( fpo, "%d%c", temp->negidlist[k],
                          temp->negdirlist[k] );
                        TotalRefsWritten++;
                        comma = 1;
                    }

                    SetSingleHashItem(
                      HASHSEEN, temp->negidlist[k], CLIST + i );
                }
            }

            // write reads
            depth = 0;

            for ( temp = CLIST + i; temp != NULL; temp = temp->next ) {
                depth++;

                for ( k = 0; k < temp->reads; k++ ) {
                    fprintf( fpo, ",%d%c", temp->idlist[k], temp->dirlist[k] );
                    TotalReadsWritten++;
                }
            }

            if ( depth > maxdepth ) {
                maxdepth = depth;
            }

            fprintf( fpo, "\n" );
        }
    }

    fclose( fpo );

    printf( "Total clusters read: %d, after join: %d, total refs(written) %d "
            "(%d),  total reads(written) %d (%d) maxdepth=%d\n",
      totclusters, finclusters, TotalRefs, TotalRefsWritten, TotalReads,
      TotalReadsWritten, maxdepth );
    return 0;
}

/*******************************************************************************************/
/* this function loads clusters from file produced by PROCLU */
int LoadClustersFromFile( char *filename, RLINK *cstart ) {
    int   i, sourcesetid, cluserid, direction = -1, clustersread = 0, pos, neg;
    RD *  a;
    int   repeatkey, repeatcount;
    FILE *fp;
    void *vptr;
    char  symbol, clustername[100], comment[510], pattern[2010] = "";
    int   maxsize, minsize, representative;
    EASY_LIST *elist;
    EASY_NODE *nof1;
    RLINK *    citem = NULL;

    /* open cluster file for reading */
    fp = fopen( filename, "r" );

    if ( fp == NULL ) {
        printf( "\n\nUnable to open intermediate cluster file!" );
        exit( 1 );
        return 0;
    }

    /* loop for every cluster */
    while ( 1 ) {
        /* read the next symbol (ignore white space) */
        i = fscanf( fp, " %c", &symbol );

        if ( i != 1 )
            break;

        if ( symbol != '@' ) {
            printf( "Invalid format detected in cluster file 1!" );
            exit( 1 );
        }

        /* increase cluster count */
        clustercount++;
        clustersread++;

        if ( NULL != cstart && clustercount > totclusters ) {
            printf( "\n\nError: Too many clusters? clustercount=%d "
                    "totclusters=%d\n\n",
              clustercount, totclusters );
            exit( 1 );
        }

        if ( NULL != cstart ) {
            citem      = cstart + clustercount - 1;
            citem->cid = clustercount;

            /* create temp list */
            elist = EasyListCreate( NULL, free );
        }

        /* loop for every member */
        representative = 0;
        repeatcount    = 0;

        while ( 1 ) {
            /* read number */
            i = fscanf( fp, "%d", &repeatkey );

            if ( i != 1 ) {
                printf( "Invalid format detected in cluster file 2!" );
                exit( 1 );
            }

            repeatcount++;

            if ( repeatcount == 1 )
                representative = repeatkey;

            /* read symbol */
            i = fscanf( fp, " %c", &symbol );

            if ( i != 1 || ( symbol != '\'' && symbol != '"' ) ) {
                printf(
                  "\n\nUnexpected value qualifier detected in cluster file!" );
                exit( 1 );
            }

            /* assign direction */
            if ( symbol == '\'' )
                direction = 0;
            else
                direction = 1;

            /* do something with it */
            // printf(" %d%c",repeatkey,symbol);
            if ( NULL != cstart ) {
                EasyListInsertHead( elist, RDCreate( repeatkey, direction ) );
            }

            /* read a character */
            symbol = 0;
            i      = fscanf( fp, " %c", &symbol );

            /* if no more characters then break */
            if ( i != 1 )
                break;

            /* if comment, read it and  then break */
            if ( symbol == '#' ) {
                fscanf( fp, "%s", comment );
                // fscanf(fp," %c",&symbol); // scan extra character
                break;
            }

            /* put character back */
            ungetc( symbol, fp );

            /* if character was a digit continue otherwise break */
            if ( ( symbol >= '0' && symbol <= '9' ) || symbol == '-' )
                continue;
            else
                break;
        }

        // printf("\n");
        if ( NULL != citem ) {

            /* count pos and neg */
            pos = 0;
            neg = 0;

            for ( nof1 = elist->head; nof1 != NULL; nof1 = nof1->next ) {
                a = (RD *) EasyListItem( nof1 );

                if ( a->id >= 0 )
                    pos++;
                else
                    neg++;
            }

            TotalReads += pos;
            TotalRefs += neg;

            citem->idlist     = smalloc( pos * sizeof( int ) );
            citem->dirlist    = smalloc( pos * sizeof( char ) );
            citem->negidlist  = smalloc( neg * sizeof( int ) );
            citem->negdirlist = smalloc( neg * sizeof( char ) );
            citem->next       = NULL;

            /* insert into RLINK structure and free temp list */
            pos = 0;
            neg = 0;

            for ( nof1 = elist->head; nof1 != NULL; nof1 = nof1->next ) {
                a = (RD *) EasyListItem( nof1 );

                if ( a->id >= 0 ) {
                    citem->idlist[pos]  = a->id;
                    citem->dirlist[pos] = a->dir;
                    pos++;
                } else {
                    citem->negidlist[neg]  = a->id;
                    citem->negdirlist[neg] = a->dir;
                    neg++;
                }
            }

            citem->refs  = neg;
            citem->reads = pos;
            EasyListDestroy( elist );
        }
    }

    /* close file */
    fclose( fp );

    return clustersread;
}
