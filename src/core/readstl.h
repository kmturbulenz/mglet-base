#ifndef __READSTL_H__
#define __READSTL_H__

#include <string>       // std::string
#include <vector>       // std::vector

extern "C" {
#include <ISO_Fortran_binding.h>
}

extern "C" {
    void readstl(CFI_cdesc_t* data, const char* filename, int* ierr);
}

int stl_is_binary(const std::string& filename);
int stl_read_ascii(std::vector<float>& vertices, const std::string& filename);
int stl_read_binary(std::vector<float>& vertices, const std::string& filename);

inline char* skip_whitespace(const char* p, const char* endp)
{
    char* q = const_cast<char*>(p);
    while (q < endp)
    {
        if (*q == ' ' || *q == '\t')
        {
            ++q;
        }
        else
        {
            break;
        }
    }
    return q;
}

#endif /* __READSTL_H__ */
