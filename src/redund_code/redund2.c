/*
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Lesser General Public
 *   License as published by the Free Software Foundation; either
 *   version 2.1 of the License, or (at your option) any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with this library; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *   ************************************************************
 *
 *   Please direct all questions related to this program or bug reports or
 * updates to Yevgeniy Gelfand (ygelfand@bu.edu)
 *
 *   Authors: Yevgeniy Gelfand
 */

/* CHANGES

  Apr 4 - changed minrepresentation to use reverse compliment from the file
  instead of RC of the forward concensus

        - reading profile now exits correctly on error

*/


//compile with:
//gcc redund2.c easylife.c -lsqlite3 -o redund2
//run with:
//'./redund2 /Users/gary/mac2008/programs/github/VNTRseek/src/redund2_code/data/ all.leb36 -i'

//#define RECORDS_PER_FILE (5000)
#define RECORDS_PER_FILE ( 100000 )
#define OUTPREFIX "reads"
#define MAX_BUFFER_SIZE 1000

/***************************************************************
    redund.c    :   Program that takes a file with a list of
                    profiles and remove redundancy based on alphabetic rotation
   of profile

                Usage:

            redund.exe INPUTFILE OUTPUTFILE [-s]
            OR
            redund.exe INPUTFOLDER OUTPUTFILE [-n]


            This also produces an index file called "outfile".removed that has
   an index of the preserved repeat followed by ids of removed repeats on each
   line. Preserved repeats that have no removed redundant repeats ARE NOT in
   that file.

            If -s switch is provided, this simply sorts the file. You will need
   to call this again to actually remove redundancy.

            Use -n switch to make the program output a single file (not broken
   up in multiples.)

    VERSION     :   1.00

*/

#include <dirent.h>
#include <errno.h>
#include <libgen.h>
//#include <malloc.h>
#include <sqlite3.h>
#include <strings.h>
#include <sys/types.h>

//#define _WIN_32_YES

#ifndef _WIN_32_YES
#define __int64 long long
#define LARGE_INTEGER time_t
#define I64d ld
#endif

/***************************************************************/
//#define EASY_LIFE_DEBUG

//#define PROCLU_DEBUG

#include <stdio.h>

#ifdef _WIN_32_YES
#include <windows.h>
#endif

#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

char *Complementascii = NULL;

#include "../libs/easylife/easylife.h"
//#include "easylife.h"
#include "profile.h"

#include <sys/resource.h>

/*******************************************************************************************/
typedef struct{
    PROFILE *prof, *profrc;
    int dir;
} RECORD_STRUCT;

typedef struct {
    char *   inputfile;
    char *   outputfile;
    char *   outputfile2;
    char *   d_name;
    FILE *   in;
    PROFILE *prof, *profrc;
    int *    minRepresentation;
    int      minrlen;
    int      dir;
    RECORD_STRUCT **buffer;
    int buffercount;
    int bufferindex;
    int inputclosed;
} FITEM_STRUCT;

typedef struct {
    PROFILE *prof, *profrc;
    int *    minRepresentation;
    int      minrlen;
} PITEM_STRUCT;

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
int pindcmp( int *p1, int *p2, int len1, int len2 ) {

//compares profile lengths and profile integer values (translated previously from leb36 format)
//-1 = first smaller
// +1 = second smaller
//0 = equal

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
int name_cmp( const void *item1, const void *item2 ) {

    FITEM_STRUCT *fiptr1, *fiptr2;

    fiptr1 = ( (FITEM_STRUCT *) item1 );
    fiptr2 = ( (FITEM_STRUCT *) item2 );

    return strcmp( fiptr1->d_name, fiptr2->d_name );
}

/*******************************************************************************************/
int arsize_and_min_rep_cmp( const void *item1, const void *item2 ) {
//compares min representation lengths and the representations themselves
//-1 = first smaller or only second null
// +1 = second smaller or only first null
//0 = equal or both null

    FITEM_STRUCT *fiptr1, *fiptr2;

    fiptr1 = ( (FITEM_STRUCT *) item1 );
    fiptr2 = ( (FITEM_STRUCT *) item2 );

    // 0th - NULL?
    if ( NULL == fiptr1->minRepresentation &&
         NULL == fiptr2->minRepresentation )
        return 0;

    if ( NULL == fiptr1->minRepresentation )
        return 1;

    if ( NULL == fiptr2->minRepresentation )
        return -1;

    /*
        // 1st - arlen
        if (fiptr1->prof->acgtCount > fiptr2->prof->acgtCount)
            return 1;
        else if (fiptr1->prof->acgtCount < fiptr2->prof->acgtCount)
            return -1;
    */
    // 2nd
    return pindcmp( fiptr1->minRepresentation, fiptr2->minRepresentation,
      fiptr1->minrlen, fiptr2->minrlen );
}

/*******************************************************************************************/
int arsize_and_min_rep_cmp_pitem( const void *item1, const void *item2 ) {

    PITEM_STRUCT *fiptr1, *fiptr2;

    fiptr1 = ( (PITEM_STRUCT *) item1 );
    fiptr2 = ( (PITEM_STRUCT *) item2 );

    // 0th - NULL?
    if ( NULL == fiptr1->minRepresentation &&
         NULL == fiptr2->minRepresentation )
        return 0;

    if ( NULL == fiptr1->minRepresentation )
        return 1;

    if ( NULL == fiptr2->minRepresentation )
        return -1;

    /*
        // 1st - arlen
        if (fiptr1->prof->acgtCount > fiptr2->prof->acgtCount)
            return 1;
        else if (fiptr1->prof->acgtCount < fiptr2->prof->acgtCount)
            return -1;
    */

    // 2nd
    return pindcmp( fiptr1->minRepresentation, fiptr2->minRepresentation,
      fiptr1->minrlen, fiptr2->minrlen );
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

/*******************************************************************************************/
int ftouch( char *in ) {

    FILE *fp;

    if ( ( fp = fopen( in, "r" ) ) == NULL )
        return 0;

    fclose( fp );
    return 1;
}

int RC = 0;

/*******************************************************************************************/
int *MinimumRepresentation( int *indptr, int len, int *indptrrc, int lenrc,
  int *minrlen, int identical_only ) {

//finds the minimum profile and profile length from the profile and reverse complement profile
//when identical_only == 1, no profile rotations are performed
//this is the case with the current VNTRseek

    int  i, j;
    int *src, *srcrc, *tar,
      // *sourceptr,
      // *destinptr,
      tmp //,
      // *indices
      ;

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

/*******************************************************************************************/
//restores the heap after adding a new element to the top
void siftDown(EASY_ARRAY *array, int start, int n) {
    FITEM_STRUCT *temp, *fiptr1, *fiptr2;
    int swap, i, j, test;

    swap = start;
    i=2*start+1; //left child
    while (i < n) {
        //put smallest element on top
        //test left child
        fiptr1 = (FITEM_STRUCT *)EasyArrayItem(array,swap);
        fiptr2 = (FITEM_STRUCT *)EasyArrayItem(array,i);

        test = arsize_and_min_rep_cmp(fiptr2, fiptr1);
        //printf("\nGot here in siftDown");

        if (test < 0){
            swap=i;
            fiptr1 = fiptr2;
        }
        //test right child
        j=i+1;
        if (j < n){
            fiptr2 = (FITEM_STRUCT *)EasyArrayItem(array,j);
            test = arsize_and_min_rep_cmp(fiptr2, fiptr1);
            if(test < 0) {
                swap = j;
                fiptr1 = fiptr2;
            }
        }
        //break if no change
        if (swap == start) break;
        //else swap
        temp = EasyArrayItem(array,start);
        //have to set these directly because there's no EasyArraySet function;
        array->array[start] = fiptr1;
        array->array[swap] = temp;
        //printf("\nGot here in siftDown");
        //reset and continue
        start = swap;
        i = 2*start+1;
    }
    //printf("\nGot here end of siftDown");
}

/*******************************************************************************************/
void heapify(EASY_ARRAY *array, int n) {
    //n is length of array[0..n-1]
    int start;

    for(start=(n-2)/2; start>=0; start--) {
        siftDown(array,start,n);
    }
    //printf("\nGot here end of heapify");
}

    /*******************************************************************************************/
void readRecordsFromFileToBuffer(FITEM_STRUCT *fiptr) {
    int  bn, buffercounter=0;
    RECORD_STRUCT * bufferptr, *tempbufferptr;

    //allocate buffer pointer memory on reading first records
    if (fiptr->buffer == NULL)
        fiptr->buffer = (RECORD_STRUCT **)calloc(MAX_BUFFER_SIZE, sizeof(RECORD_STRUCT *));

    //free records memory on subsequent rounds before allocating new records
    if(fiptr->buffer[0] != NULL)
        free(fiptr->buffer[0]);

    //allocate new records memory
    tempbufferptr = (RECORD_STRUCT *)calloc(MAX_BUFFER_SIZE, sizeof(RECORD_STRUCT)); 
    for(bn = 0; bn < MAX_BUFFER_SIZE; bn++) {
        fiptr->buffer[bn] = tempbufferptr;
        tempbufferptr++;
    }

    //printf("\n\nBefore read records for file %s, bufferindex = %d, buffercount = %d",fiptr->inputfile, fiptr->bufferindex, fiptr->buffercount);

    //read records
    //printf("\nReading records from %s",fiptr->inputfile);
    while (!feof(fiptr->in) && buffercounter<MAX_BUFFER_SIZE) {
        bufferptr = fiptr->buffer[buffercounter];
        bufferptr->prof = ReadProfileWithRC( fiptr->in, &( bufferptr->profrc ));
        //test, because could be empty line
        if (bufferptr->prof != NULL) {
            buffercounter++;
        }
    }

    //set the counts
    fiptr->buffercount = buffercounter;
    fiptr->bufferindex = 0;
    //printf("\nRead %d records",fiptr->buffercount);

    //printf("\nAfter read records for file %s, bufferindex = %d, buffercount = %d",fiptr->inputfile, fiptr->bufferindex, fiptr->buffercount);

    //close the file if finished
    if(fiptr->buffercount == 0) {
        //printf("\nNo more records, closing file %s.",fiptr->inputfile);
        fclose(fiptr->in);
        fiptr->inputclosed = 1;
    }

}

/*******************************************************************************************/
void loadRecordFromBufferAndGetMinRepresetation (FITEM_STRUCT * fiptr, int identical_only) {

    int i = fiptr->bufferindex;

    //load record
    fiptr->prof = fiptr->buffer[i]->prof;
    fiptr->profrc = fiptr->buffer[i]->profrc;

    /*** this set of dir may be incorrect ***/
    //	fiptr->dir = fiptr->buffer[i]->dir;
    //probably don't need this
    //fiptr->buffer[i] = NULL;

    fiptr->bufferindex++;

    //compute minRepresentation
    if ( NULL != fiptr->prof ) {
        fiptr->minRepresentation =  MinimumRepresentation( fiptr->prof->indices, fiptr->prof->proflen,
                fiptr->profrc->indices, fiptr->profrc->proflen,
                &( fiptr->minrlen ), identical_only );

    /*** added RC here ***/
        fiptr->dir = RC;
    }
    else {
        fiptr->minRepresentation = NULL;
    }
}
/*******************************************************************************************/
int main( int argc, char **argv ) {
    int   SINGLE_OUTFILE, SORT_ONLY, IDENTICAL_ONLY;
    FILE *fpto, *fpto2;
    char *bigtempbuf, *inputfile, *outputdname, *outputbname, *outputfile,
      *outputfile2, *outdb, *err_msg = 0;
    int            i, filescreated = 1, rc;
    time_t         startTime;
    FITEM_STRUCT * fiptr, *fiptr2, *lastwrite;
    PITEM_STRUCT * piptr;
    EASY_ARRAY *   FARRAY = NULL;
    struct dirent *de     = NULL;
    DIR *          d      = NULL;
    long long int  nwritten, nread;
    sqlite3 *      db;
    sqlite3_stmt * res, *pStmt;
    struct rlimit old_lim, lim, new_lim; 
    int filecounter, buffercounter, bn, softlimit;
    RECORD_STRUCT *tempbufferptr;
    struct rusage *usage;

    usage = (struct rusage *)calloc(1, sizeof(struct rusage));

    // verify parameter count
    if ( argc < 3 ) {
        printf(
          "\nREDUND v2.00 - removes redundancy based on rotated profiles." );
        printf( "\nauthors: Yevgeniy Gefland, Gary Benson" );
        printf( "\n\nUsage:\n\n%s INPUTFILE OUTPUTFILE [-s]\n\tOR\n%s "
                "INPUTFOLDER OUTPUFFILE\n",
          argv[0], argv[0] );
        printf( "   -s options will sort the file (only for single file input) "
                "\n\n\n" );
        printf( "   -n options will output single file (will break in multiple "
                "otherwise) \n\n\n" );
        printf( "   -i options will only remove identical profiles, without "
                "rotating \n\n\n" );

        exit( 1 );
    }

    init_complement_ascii();

    // parse parameters
    size_t indirlen = strlen( argv[1] );
    if ( '/' == argv[1][strlen( argv[1] ) - 1] )
        indirlen--;
    inputfile = strndup( argv[1], indirlen );
    if ( !inputfile ) {
        perror( "Error copying input directory" );
        exit( 1 );
    }

    IDENTICAL_ONLY = 0;

    if ( argc >= 4 &&
         ( 0 == strcmp( "-I", argv[3] ) || 0 == strcmp( "-i", argv[3] ) ) ) {
        IDENTICAL_ONLY = 1;
    }

    if ( argc >= 5 &&
         ( 0 == strcmp( "-I", argv[4] ) || 0 == strcmp( "-i", argv[4] ) ) ) {
        IDENTICAL_ONLY = 1;
    }

    if ( argc >= 6 &&
         ( 0 == strcmp( "-I", argv[5] ) || 0 == strcmp( "-i", argv[5] ) ) ) {
        IDENTICAL_ONLY = 1;
    }

    SINGLE_OUTFILE = 0;

    if ( argc >= 4 &&
         ( 0 == strcmp( "-N", argv[3] ) || 0 == strcmp( "-n", argv[3] ) ) ) {
        SINGLE_OUTFILE = 1;
    }

    if ( argc >= 5 &&
         ( 0 == strcmp( "-N", argv[4] ) || 0 == strcmp( "-n", argv[4] ) ) ) {
        SINGLE_OUTFILE = 1;
    }

    if ( argc >= 6 &&
         ( 0 == strcmp( "-N", argv[5] ) || 0 == strcmp( "-n", argv[5] ) ) ) {
        SINGLE_OUTFILE = 1;
    }

    SORT_ONLY = 0;

    if ( argc >= 4 &&
         ( 0 == strcmp( "-S", argv[3] ) || 0 == strcmp( "-s", argv[3] ) ) ) {
        SORT_ONLY      = 1;
        SINGLE_OUTFILE = 1;
    }

    if ( argc >= 5 &&
         ( 0 == strcmp( "-S", argv[4] ) || 0 == strcmp( "-s", argv[4] ) ) ) {
        SORT_ONLY      = 1;
        SINGLE_OUTFILE = 1;
    }

    if ( argc >= 6 &&
         ( 0 == strcmp( "-S", argv[5] ) || 0 == strcmp( "-s", argv[5] ) ) ) {
        SORT_ONLY      = 1;
        SINGLE_OUTFILE = 1;
    }

    char *tmp   = strdup( argv[2] );
    outputbname = strdup( basename( tmp ) );
    outputdname = strdup( dirname( tmp ) );
    free( tmp );
    outputfile  = calloc( strlen( outputdname ) + strlen( outputbname ) + 10,
      sizeof( *outputfile ) );
    outputfile2 = calloc( strlen( outputdname ) + strlen( outputbname ) + 19,
      sizeof( *outputfile2 ) );
    outdb       = calloc(
      strlen( outputdname ) + strlen( outputbname ) + 4, sizeof( *outdb ) );

    if ( SINGLE_OUTFILE ) {
        sprintf( outputfile, "%s", argv[2] );
        sprintf( outputfile2, "%s.rotindex", argv[2] );
        sprintf( outdb, "%s.db", argv[2] );
    } else {
        sprintf( outputfile, "%s/1.%s", outputdname, outputbname );
        sprintf( outputfile2, "%s/1.%s.rotindex", outputdname, outputbname );
        sprintf( outdb, "%s/%s.db", outputdname, outputbname );

        if ( getenv( "DEBUG" ) && strcmp( "1", getenv( "DEBUG" ) ) == 0 ) {
            fprintf(
              stderr, "Dirname: %s, basename: %s\n", outputdname, outputbname );
            fprintf( stderr, "outputfile %s\n", outputfile );
            fprintf( stderr, "outputfile2 %s\n", outputfile2 );
        }
    }


    //Started 8/25/20 Gary Benson
    //changes to convert to a merge sort of the existing leb36 files which are already sorted
    //Basic outline
    //  Set open file limit high enough for the number of files plus some extra
    //  Open all the files
    //  read first MAX_BUFFER_SIZE records from each file into a buffer for each file
    //  move the first record from each buffer into a heap
    //  heapify to find the smallest record and write to the output file
    //  reload a record from the buffer which had the smallest record
    //  if a buffer is empty, read the next 10,000 records from that buffer's file
    //  if the file is empty,
    //    close the file
    //    do something in the heap to get the next smallest record
    //    keep repeating until the heap is empty




    // create a file array
    //expands when necessary
    FARRAY = EasyArrayCreate( 1000, NULL, NULL );

     startTime = time( NULL );
    //open directory and get all files
    // list file(s)
    d = opendir( inputfile );


    if ( d != NULL ) {
        // fprintf( stderr, "inputfile last char: %c\n",
        //   *( inputfile + strlen( inputfile ) - 1 ) );
        //bigtempbuf holds name of input directory
        bigtempbuf = calloc( indirlen + 30, sizeof( *bigtempbuf ) );  
        // fprintf( stderr, "inputfile: %s\n", inputfile );

        filecounter=0;
        while ( ( de = readdir( d ) ) != NULL ) {

            if ( strlen( de->d_name ) > 17 &&
                 0 == strcasecmp( ".leb36.renumbered",
                        de->d_name + strlen( de->d_name ) - 17 ) ) {
                //fiptr is a FITEM_STRUCT which contains too much information
                //only needs inputfile path and name, file name, file pointer
                filecounter++;
                fiptr = smalloc( sizeof( FITEM_STRUCT ) );

                // strcpy( bigtempbuf, inputfile );
                sprintf( bigtempbuf, "%s/%s", inputfile, de->d_name );

                fiptr->inputfile = strdup( bigtempbuf );
                fiptr->d_name    = strdup( de->d_name );
                EasyArrayInsert( FARRAY, fiptr );
            }
        }

        closedir( d );

        // Set open file limit high enough for the number of files plus some extra

        //get existing limits
        getrlimit(RLIMIT_NOFILE, &old_lim);

        softlimit = filecounter+1000;
        if (softlimit > old_lim.rlim_max) {
            fprintf( stderr, "Open file limit will be exceeded for .leb36.renumbered files; %d files to be read and %d + 1000 > hard limit = %llu\n", softlimit, softlimit,  old_lim.rlim_max );
            exit (1);
        }

        //adjust soft limit
        lim.rlim_cur = softlimit; 
        lim.rlim_max = old_lim.rlim_max;
        setrlimit(RLIMIT_NOFILE, &lim);

        //debug
        getrlimit(RLIMIT_NOFILE, &new_lim);
        ////print to confirm
        printf("Open file limits:\nOriginal: soft: %llu, hard: %llu\nNew: soft: %llu, hard: %llu\n",old_lim.rlim_cur, old_lim.rlim_max, new_lim.rlim_cur, new_lim.rlim_max);

    }
    else {  //input 1 not a directory, means single file, not used in VNTRseek anymore

        if ( ftouch( inputfile ) ) {

            // -s options to sort single file
            if ( SORT_ONLY ) {

                EASY_LIST *profileList;
                EASY_NODE *nof1;
                FILE *     fpi;

                fpi = fopen( inputfile, "r" );

                if ( NULL == fpi ) {
                    printf( "\nERROR: Unable to open input file '%s'\n\n",
                      inputfile );
                    exit( 1 );
                }

                profileList = EasyListCreate( NULL, NULL );

                while ( 1 ) {

                    piptr = smalloc( sizeof( PITEM_STRUCT ) );

                    piptr->prof = ReadProfileWithRC( fpi, &( piptr->profrc ) );
                    piptr->minRepresentation = NULL;

                    if ( NULL == piptr->prof || NULL == piptr->profrc ) {
                        free( piptr );
                        break;
                    }

                    // find minimum representation and put in list
                    piptr->minRepresentation = MinimumRepresentation(
                      piptr->prof->indices, piptr->prof->proflen,
                      piptr->profrc->indices, piptr->profrc->proflen,
                      &( piptr->minrlen ), IDENTICAL_ONLY );

                    if ( NULL == piptr->minRepresentation ) {
                        printf( "\nERROR: minrepresentation is NULL!\n\n" );
                        exit( 1 );
                    }

                    EasyListInsertTail( profileList, piptr );
                }

                // open the output file for writing
                fpto = fopen( outputfile, "w" );

                if ( fpto == NULL ) {
                    printf( "\nERROR: Unable to open output file '%s'\n\n",
                      outputfile );
                    exit( 1 );
                }

                rc = sqlite3_open( outdb, &db );

                if ( rc != SQLITE_OK ) {

                    fprintf( stderr, "Cannot open database: %s\n",
                      sqlite3_errmsg( db ) );
                    sqlite3_close( db );

                    return 1;
                }

                rc = sqlite3_exec(
                  db, "PRAGMA synchronous = OFF", 0, 0, &err_msg );

                if ( rc != SQLITE_OK ) {
                    fprintf( stderr, "SQL error: %s\n", err_msg );

                    sqlite3_free( err_msg );
                    sqlite3_close( db );

                    return 1;
                }

                rc = sqlite3_exec( db,
                  "CREATE TABLE minreporder (`rid` integer PRIMARY KEY,"
                  "`idx` integer)",
                  0, 0, &err_msg );

                if ( rc != SQLITE_OK ) {
                    fprintf( stderr, "SQL error: %s\n", err_msg );

                    sqlite3_free( err_msg );
                    sqlite3_close( db );

                    return 1;
                }

                char *sql = "INSERT INTO minreporder VALUES(?, ?)";

                rc = sqlite3_prepare_v2( db, sql, -1, &pStmt, 0 );

                if ( rc != SQLITE_OK ) {
                    fprintf( stderr, "Cannot prepare statement: %s\n",
                      sqlite3_errmsg( db ) );

                    return 1;
                }

                // sort and output
                fclose( fpi );
                EasyListQuickSort( profileList, arsize_and_min_rep_cmp_pitem );
                i = 0;

                for ( nof1 = profileList->head; nof1 != NULL;
                      nof1 = nof1->next ) {
                    piptr = (PITEM_STRUCT *) EasyListItem( nof1 );
                    WriteProfileWithRC( fpto, piptr->prof, piptr->profrc );
                    sqlite3_bind_int( pStmt, 1, piptr->prof->key );
                    sqlite3_bind_int( pStmt, 2, i++ );
                    sqlite3_step( pStmt );
                    sqlite3_reset( pStmt );
                }

                sqlite3_finalize( pStmt );
                sqlite3_close( db );

                fclose( fpto );
                printf( "\n\nredund.exe: Input File sorted by minimum "
                        "representation. Now call this program again without "
                        "-s switch to remove redundancy. \n\n" );
                fflush( stdout );

                exit( 0 );
            }

            fiptr            = smalloc( sizeof( FITEM_STRUCT ) );
            fiptr->inputfile = strdup( inputfile );
            fiptr->d_name =
              strdup( inputfile ); // this won't be used in this case
            EasyArrayInsert( FARRAY, fiptr );

        } else {
            printf(
              "\nERROR: Unable to read directory or file '%s'\n\n", inputfile );
            exit( 1 );
        }
    }

    // first order by file name, because readdir is not sorted by default
    EasyArrayQuickSort( FARRAY, name_cmp );

    // open input file(s) for reading
    printf( "Opening files.\n" );

    if ( 0 == FARRAY->size ) {
        printf( "\nERROR: No input files found in the directory.\n\n" );
        exit( 1 );
    }

    for ( i = 0; i < filecounter; i++ ) {
        fiptr     = (FITEM_STRUCT *) EasyArrayItem( FARRAY, i );
        fiptr->in = fopen( fiptr->inputfile, "r" );
        fiptr->inputclosed = 0;

        if ( NULL == fiptr->in ) {
            printf( "\nERROR opening input file '%s': %s\n", fiptr->inputfile,
              strerror( errno ) );
            exit( 1 );
        }

        //printf( "\t%s\n", fiptr->inputfile );
    }

    printf("There were %d files opened.\n",filecounter);

    // open the output file for writing
    fpto = fopen( outputfile, "w" );

    if ( fpto == NULL ) {
        printf( "\nERROR: Unable to open output file '%s'.\n\n", outputfile );
        exit( 1 );
    }

    // open the index file for writing
    fpto2 = fopen( outputfile2, "w" );

    if ( fpto2 == NULL ) {
        printf( "\nERROR: Unable to open index file '%s'\n\n", outputfile2 );
        exit( 1 );
    }

    //read first MAX_BUFFER_SIZE records from each file into a buffer for each file
    for ( i = 0; i < filecounter; i++ ) {
        fiptr = (FITEM_STRUCT *) EasyArrayItem( FARRAY, i );
        //read records
        readRecordsFromFileToBuffer(fiptr);
        if (fiptr->buffercount == 0){
            printf( "\nERROR: Input file %s was empty\n\n", fiptr->inputfile);
            exit( 1 );
        }
    }

    //getrusage(RUSAGE_SELF, usage);
    //printf("\nMemory usage: %ld",usage->ru_maxrss);

    //put smallest record in each buffer directly into the top level of fiptr for each FARRAY element and compute the minRepresentation
    for ( i = 0; i < filecounter; i++ ) {
        fiptr = (FITEM_STRUCT *) EasyArrayItem( FARRAY, i);
        loadRecordFromBufferAndGetMinRepresetation (fiptr, IDENTICAL_ONLY); 
    }

    //test that filecounter and FARRAY->size are the same
    if (filecounter!=FARRAY->size){
        printf("\nError: filecounter: %d not equal to FARRAY->size: %zu",filecounter, FARRAY->size);
        exit(1);
    }

    //printf("\nGot here");

    //treat FARRAY as the heap
    //heapify FARRAY
    heapify(FARRAY,filecounter);


    //set up for loop
    nread     = 0;   //***this has to be moved up to incorporate the records read before
    nwritten  = 0;
    lastwrite = NULL;

    while (1) { //break on empty heap

        //get smallest record for writing
        //printf("\nGot here before first fiptr");
        fiptr = (FITEM_STRUCT *) EasyArrayItem( FARRAY, 0 );

        //test if done
        //***alternately, test for filecounter == 0
        if (fiptr->prof == NULL)
            break;

        nread++;

        // decide how to write and write profiles
        // first if identifies duplicates of the last written record
        // test if all the below occur
        //   1: there is a last written record
        //   2: for this record and the last written one
        //     the minimum profiles match
        //   3: for this record and the last written one
        //     either the two forward profiles match,
        //     the two reverse profiles match,
        //     or there is a match between one profile and one reverse profile
        // The last condition seems to be redundant, but probably assures that the two profiles are not NULLwhich would return a match for condition 2. 

        if ((lastwrite != NULL) &&
            (0 == arsize_and_min_rep_cmp( lastwrite, fiptr ) ) &&

            // (ADDED apr 18, 2013)
            ( ( 0 == pindcmp( fiptr->prof->indices, lastwrite->prof->indices,
                        fiptr->prof->proflen, lastwrite->prof->proflen ) &&
                0 == pindcmp( fiptr->profrc->indices,
                        lastwrite->profrc->indices, fiptr->profrc->proflen,
                        lastwrite->profrc->proflen ) ) ||
             ( 0 == pindcmp( fiptr->prof->indices, lastwrite->profrc->indices,
                        fiptr->prof->proflen, lastwrite->profrc->proflen ) &&
                0 == pindcmp( fiptr->profrc->indices, lastwrite->prof->indices,
                        fiptr->profrc->proflen,
                        lastwrite->prof->proflen ) )
            ) ) 
            {

            fprintf( fpto2, " %d%c", fiptr->prof->key,
              ( 0 == fiptr->dir ) ? '\'' : '\"' ); //add duplicate index entry

            // in proclu 1.87 we need all profiles written

            WriteProfileWithRC( fpto, fiptr->prof, fiptr->profrc );

            //don't free lastwrite here, because still in use

        } else {

            //else condition if for non-duplicates
            nwritten++;

            // write out

            WriteProfileWithRC( fpto, fiptr->prof, fiptr->profrc );
            //printf("\nGot here after first write");
            //fflush(stdout);


            if ( NULL != lastwrite ) {
                fprintf( fpto2, "\n" ); //start new preserved index entry
            }

            fprintf( fpto2, "%d%c", fiptr->prof->key,
              ( 0 == fiptr->dir ) ? '\'' : '\"' ); // add preserved index entry

                    // make sure new minrep is equal or larger than old minrep


            //free lastwrite here, because no longer used


            if ( NULL != lastwrite ) {          //*** may have to free other things
                FreeProfile( lastwrite->prof );
                FreeProfile( lastwrite->profrc );
                free( lastwrite->minRepresentation );
                free( lastwrite );
            }


            //set lastwrite for comparison with next smallest record
            //*** may not need all this information
            lastwrite = smalloc( sizeof( FITEM_STRUCT ) );
            lastwrite->inputfile   = fiptr->inputfile;
            lastwrite->outputfile  = fiptr->outputfile;
            lastwrite->outputfile2 = fiptr->outputfile2;
            lastwrite->in = fiptr->in;
            lastwrite->prof = CopyProfile( fiptr->prof );
            lastwrite->profrc = CopyProfile( fiptr->profrc );
            lastwrite->minRepresentation = pintdup( fiptr->minRepresentation, fiptr->minrlen );
            lastwrite->minrlen = fiptr->minrlen;

            //printf("\nGot here after lastwrite assignment");
            //fflush(stdout);

            // open new file?
            //shouldn't this test be before the record is written?
            if ( !SINGLE_OUTFILE && ( nwritten % ( RECORDS_PER_FILE ) ) == 0 ) {
                //printf("\nGot here just after test for multiple output files");
                //fflush(stdout);
                fclose( fpto );
                fclose( fpto2 );

                filescreated++;

                // open the output file for writing
                sprintf( outputfile, "%s/%d.%s", outputdname, filescreated,
                  outputbname );
                fpto = fopen( outputfile, "w" );

                if ( fpto == NULL ) {
                    printf( "\nERROR: Unable to open output file '%s'\n\n",
                      outputfile );
                    exit( 1 );
                }

                // open the index file for writing
                sprintf( outputfile2, "%s/%d.%s.rotindex", outputdname,
                  filescreated, outputbname );
                fpto2 = fopen( outputfile2, "w" );

                if ( fpto2 == NULL ) {
                    printf( "\nERROR: Unable to open output file '%s'\n\n",
                      outputfile2 );
                    exit( 1 );
                }
            }
        }


        //at this point, free the profiles and minRepresentation
        FreeProfile( fiptr->prof );
        FreeProfile( fiptr->profrc );
        free( fiptr->minRepresentation );

        //debug test
        //printf("\nGot here before test of null fiptr");
        //fflush(stdout);
        //if (fiptr == NULL) {
        //    printf("\nfiptr is null");
        //    fflush(stdout);
        //}
        //printf("\nGot here after test of null fiptr");
        //fflush(stdout);

         //load new profiles or put last item in heap at top and shorten heap
         //check buffer not empty
         if (fiptr->bufferindex != fiptr->buffercount){
                loadRecordFromBufferAndGetMinRepresetation (fiptr, IDENTICAL_ONLY);
         }
         else
         {
            //buffer is empty
            //read records from file into the buffer
            //debug print
            readRecordsFromFileToBuffer(fiptr);

            //if some records were read, load new profile
            if(fiptr->buffercount != 0){
                loadRecordFromBufferAndGetMinRepresetation (fiptr, IDENTICAL_ONLY);

            }
            else
            {
                //no more records for this fiptr
                //no profiles or minRepresentations remain
                free(fiptr->buffer[0]);
                free(fiptr->buffer);
                free(fiptr);

                //get the last FARRAY item so it can be moved to the top
                //it already has the profile and minRepresentation information 
                fiptr = EasyArrayItem(FARRAY, filecounter-1);
                FARRAY->array[0] = fiptr;

                //reduce filecounter by 1
                filecounter--;
                //printf("\nFiles left = %d",filecounter);

                //test if heap empty and if so, break out of while
                if(filecounter == 0)
                    break;
            }

         }

        //fix heap
        siftDown(FARRAY, 0, filecounter);

    }


    free( outputbname );
    free( outputdname );
    free( outputfile );
    free( outputfile2 );
    free( outdb );

    printf( "%llu profiles read, %llu profiles marked nonredundant. (time: "
            "%ld seconds)\n",
    nread, nwritten, time( NULL ) - startTime );
    fflush( stdout );

    getrusage(RUSAGE_SELF, usage);
    printf("Memory usage: %ld\n",usage->ru_maxrss);
    //free(usage);

    return 0;
}

/*******************************************************************************************/
void doCriticalErrorAndQuit( const char *format, ... ) {

    va_list argp;

    printf( "\n\nERROR!!!!:\n" );

    if ( format == NULL )
        exit( 0 );

    va_start( argp, format );
    vprintf( format, argp );
    fflush( stdout );
    va_end( argp );

    exit( 1 );
}
