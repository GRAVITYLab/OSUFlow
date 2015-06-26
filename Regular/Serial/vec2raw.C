#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "FileReader.h"

int main(int argc, char ** argv)
{
    printf("vec2raw vec-file\n");

    int dim[3];

    float *data = ReadStaticDataRaw(argv[1], dim);

    printf("Dimension: %d %d %d\n", dim[0], dim[1], dim[2]);

    char str[1024];
    strcpy(str, argv[1]);
    strcat(str, ".raw");

    FILE *fp = fopen(str, "wb");
    fwrite(data, dim[0]*dim[1]*dim[2]*3, 4, fp);
    fclose(fp);

    printf("File saved: %s\n", str);
    return 0;
}
