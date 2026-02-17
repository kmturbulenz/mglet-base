#include <stddef.h>
#include <stdint.h>
#include <sys/resource.h>

void getrusage_c(long long *maxmem, int *ierr)
{
    if (maxmem == NULL || ierr == NULL) {
        if (ierr) *ierr = 1;
        return;
    }

    struct rusage ru;
    if (getrusage(RUSAGE_SELF, &ru) != 0) {
        *ierr = -1;
        *maxmem = -1;
        return;
    }

    *maxmem = (long long)ru.ru_maxrss;

    *ierr = 0;
}
