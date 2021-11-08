#include <errno.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#define STL_CHECK_LENGTH 4096

// Possible return values:
//  -1 if the STL is ASCII, the number of triangles cannot be determined
//  value > 0 means binary, the value is the number of triangles
int stl_is_binary(const char* filename)
{
    // Get the file size
    struct stat st;
    if (stat(filename, &st))
    {
        fprintf(stderr, "Stat error: %s\n", filename);
        fprintf(stderr, "errno: %d\n", errno);
        abort();
    }
    size_t size = st.st_size;

    // The minimum file size of a binary STL is the 80 byte header + 4 bytes
    // in sum 84 bytes. If the file size is smaller, then it cannot be a
    // valid binary STL file
    if (size < 84)
    {
        return -1;
    }

    FILE *fp = fopen(filename, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Could not open: %s\n", filename);
        fprintf(stderr, "errno: %d\n", errno);
        abort();
    }

    // Check for non-printable ASCII characters in the first bytes of the file
    unsigned char tmpbuf[STL_CHECK_LENGTH];
    size_t nbytes = fread(&tmpbuf, 1, sizeof(tmpbuf), fp);
    size_t nascii = 0;
    for(size_t i = 0; i < nbytes; ++i)
    {
        // 32 is the first printable ASCII character
        // 127 is the DEL character, and not considered printable
        // 126 is thus the last printable character
        // Additional control characters \n \r and \t considered "printable"
        if ((tmpbuf[i] >= 32 && tmpbuf[i] <= 126) || tmpbuf[i] == '\n' || \
            tmpbuf[i] == '\r' || tmpbuf[i] == '\t')
        {
            ++nascii;
        }
    }

    // Read number of triangles - 80 bytes beyond file start
    fseek(fp, 80, SEEK_SET);
    unsigned int ntri;
    fread(&ntri, sizeof(ntri), 1, fp);
    fclose(fp);

    // If *only* encountered ASCII characters, this is most probably an ASCII
    // file
    if (nascii == nbytes)
    {
        return -1;
    }

    // If ntri cannot be safely casted to an int MGLET cannot use this STL
    if (ntri > 0x7FFFFFFF)
    {
        fprintf(stderr, "Too many triangles in STL: %s\n", filename);
        fprintf(stderr, "  Number of triangles: %d\n", ntri);
        fprintf(stderr, "  Exceed MGLET maximum: %d\n", 0x7FFFFFFF);
        return 0;
    }

    // Additionally, if the file size correspond to the expected file size
    // based on the number of triangles, we are absolutely positively
    // sure it is a binary file
    size_t expected_size = 84 + 50*ntri;
    if (size == expected_size)
    {
        return ntri;
    }

    // If the file size is not matching that of a binary STL, then we should
    // print a warning, but nevertheless flag it as binary, since it did
    // contain non-printable characters
    fprintf(stderr, "Error in detecting ASCII/binary STL: %s\n", filename);
    fprintf(stderr, "  File length: %ld\n", size);
    fprintf(stderr, "  Read %ld, found %ld valid chars\n", nbytes, nascii);
    fprintf(stderr, "  If this is binary, it has %d triangles\n", ntri);
    fprintf(stderr, "  Expected filesize is then %ld\n", expected_size);
    fprintf(stderr, "Now assuming binary mode and hope for the best...\n");
    return ntri;
}
