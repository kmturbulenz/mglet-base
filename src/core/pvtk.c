#include <dirent.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>  // PATH_MAX

bool endswith(const char* base, const char* suffix);

void create_pvtk(const char *path, const char *prefix) {

    size_t n_pre = strlen(prefix);

    // Open directory
    DIR *fd;
    if ((fd = opendir(path)) == NULL) {
        printf("pvtk_dir: can't open %s\n", path);
        return;
    }

    // Count the relevant files
    struct dirent *dp;
    int nfiles = 0;
    while ((dp = readdir(fd)) != NULL) {
        if (strncmp(dp->d_name, prefix, n_pre) == 0 &&
                endswith(dp->d_name, ".vtk")) {
            nfiles = nfiles + 1;
        }
    }

    // If no files are found - return
    if (nfiles <= 0) {
        return;
    }

    // Make filename
    char fname[PATH_MAX];
    int nchar = snprintf(fname, PATH_MAX, "%s/%s.pvtk", path, prefix);
    if (nchar < 0 || nchar >= PATH_MAX) {
        printf("pvtk_dir: could not make filename\n");
        return;
    }

    // Write file
    FILE *fh;
    if ((fh = fopen(fname, "w")) == NULL) {
        printf("pvtk_dir: can't open file %s", fname);
        return;
    }

    fprintf(fh, "<File version=\"pvtk-1.0\" dataType=\"vtkUnstructuredGrid\" ");
    fprintf(fh, "numberOfPieces=\"%i\">\n", nfiles);

    // Iterating over files and adding them in file
    rewinddir(fd);
    while ((dp = readdir(fd)) != NULL) {
        if (strncmp(dp->d_name, prefix, n_pre) == 0 &&
                endswith(dp->d_name, ".vtk")) {
            fprintf(fh, "  <Piece fileName=\"%s\" />\n", dp->d_name);
        }
    }

    fprintf(fh, "</File>\n");

    fclose(fh);
    closedir(fd);
}

bool endswith(const char* base, const char* suffix) {
    size_t blen = strlen(base);
    size_t slen = strlen(suffix);
    return (blen >= slen) && (0 == strcmp(base + blen - slen, suffix));
}
