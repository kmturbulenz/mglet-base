#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

void create_directory_f(const char *path) {
    struct stat st = {0};

    if (stat(path, &st) == -1) {
        mkdir(path, (mode_t)0777);
    }
}
