#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "FileReader.h"

int main(int argc, char ** argv)
{
    printf("raw2vec raw-file xdim ydim zdim\n");
    if (argc < 5)
    {
        exit(1);
    }

    int dim[3];
    dim[0] = atoi(argv[2]);
    dim[1] = atoi(argv[3]);
    dim[2] = atoi(argv[4]);

    float *data = ReadStaticDataRawNoHeader(argv[1], dim);

    char str[1024];
    strcpy(str, argv[1]);
    strcat(str, ".vec");

    FILE *fp = fopen(str, "wb");
    fwrite(dim, 3, 4, fp);
    fwrite(data, dim[0]*dim[1]*dim[2]*3, 4, fp);
    fclose(fp);

    printf("File saved: %s\n", str);
    return 0;
}
