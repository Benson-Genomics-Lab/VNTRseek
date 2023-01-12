#include <stdio.h>
#include <string.h>
#include "edlib.h"

int main( int argc, char **argv) {
    // parse args
    char query[1000] = "", target[1000] = "";

    strcpy(query, argv[1]);
    strcpy(target, argv[2]);

    // perform alignmnet
    EdlibAlignConfig cnf = edlibNewAlignConfig(
        -1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0);

    EdlibAlignResult result = edlibAlign(
        query, strlen(query), target, strlen(target), cnf);

    int alignStart  = result.startLocations[0];
    int alignEnd    = result.endLocations[0];
    unsigned char* alignment = result.alignment;
    int alignLength = result.alignmentLength;

    // print start location
    printf("%d\n", alignStart);


    // print NICE alignment
    int tIdx = -1;
    int qIdx = -1;

    tIdx = alignEnd;
    int i;
    for (i = 0; i < alignLength; i++) {
        if (alignment[i] != EDLIB_EDOP_INSERT)
            tIdx--;
    }

    // target
    for (i = 0; i < alignLength; i++) {
        if (alignment[i] == EDLIB_EDOP_INSERT)
            printf("-");
        else
            printf("%c", target[++tIdx]);
    }
    printf("\n");

    // query
    for (i = 0; i < alignLength; i++) {
        if (alignment[i] == EDLIB_EDOP_DELETE)
            printf("-");
        else
            printf("%c", query[++qIdx]);
    }
    printf("\n");

    edlibFreeAlignResult(result);
}
