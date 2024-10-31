#include "readstl.h"

#include <algorithm>    // std::min
#include <cassert>      // assert
#include <cerrno>       // errno
#include <charconv>     // std::from_chars (C++17)
#include <cstdio>       // std::sscanf
#include <fstream>      // std::ifstream
#include <iostream>     // std::cerr, std::cout
#include <string>       // std::string
#include <vector>       // std::vector

#include <sys/stat.h>   // stat
extern "C" {
#include <ISO_Fortran_binding.h>
}

#define STL_CHECK_LENGTH 4096


extern "C" {
    void readstl(CFI_cdesc_t* data, const char* filename, int* ierr)
    {
        // Sanity checks of passed Fortran structure
        assert (data->rank == 3);
        assert (data->type == CFI_type_float || data->type == CFI_type_double);
        assert (data->attribute == CFI_attribute_allocatable);

        auto fname = std::string(filename);
        int is_binary = stl_is_binary(fname);

        if (is_binary < 0)
        {
            *ierr = is_binary;
            return;
        }

        // For holding the vertices that was read from the file
        std::vector<float> vertices;

        // Read the file
        if (is_binary == 0)
        {
            *ierr = stl_read_ascii(vertices, fname);
        }
        else
        {
            *ierr = stl_read_binary(vertices, fname);
        }

        // When an error occurred, we return immediately
        if (*ierr)
        {
            return;
        }

        // Number of vertices that was read and number of triangles from this
        const size_t nverts = vertices.size();
        const size_t ntri = nverts/9;

        const std::string type = (is_binary == 0) ? "ASCII" : "binary";
        std::cout << " STL: " << fname << " is " << type
                  << ", containing " << ntri << " triangles" << std::endl;

        // If already allocated - free
        if (data->base_addr) CFI_deallocate(data);

        // Allocate the Fortran array
        const CFI_index_t lower_bounds[3] = {1, 1, 1};
        const CFI_index_t upper_bounds[3] =
            {3, 3, static_cast<CFI_index_t>(ntri)};
        *ierr = CFI_allocate(data, lower_bounds, upper_bounds, 0);
        if (*ierr) {
            std::cerr << "Error allocating ntri: " << ntri << std::endl;
            std::cerr << "Got error: " << *ierr << std::endl;
            return;
        }

        // Copy data from C++ vector to Fortran array
        if (data->type == CFI_type_float)
        {
            float* ptr = reinterpret_cast<float*>(data->base_addr);
            for (size_t i = 0; i < nverts; ++i)
            {
                ptr[i] = vertices[i];
            }
        }
        else if (data->type == CFI_type_double)
        {
            double* ptr = reinterpret_cast<double*>(data->base_addr);
            for (size_t i = 0; i < nverts; ++i)
            {
                ptr[i] = vertices[i];
            }
        }

        *ierr = 0;
        return;
    }
}


// Return number of triangles when the STL is binary, 0 when ASCII, and
// negative in case of error
int stl_is_binary(const std::string& filename)
{
    // Get the file size
    // Better alternative: std::filesystem::file_size from C++17
    struct stat st;
    if (stat(filename.c_str(), &st))
    {
        std::cerr << "Stat error: " << filename << std::endl;
        std::cerr << "errno: " << errno << std::endl;
        std::cerr << std::flush;
        return -1;
    }
    size_t size = st.st_size;

    // The minimum file size of a binary STL is the 80 byte header + 4 bytes
    // in sum 84 bytes. If the file size is smaller, then it cannot be a
    // valid binary STL file
    if (size < 84)
    {
        return 0;
    }

    std::ifstream file(filename, std::ios::binary);
    if (!file)
    {
        std::cerr << "Could not open: " << filename << std::endl;
        std::cerr << "errno: " << errno << std::endl;
        std::cerr << std::flush;
        return -2;
    }

    // Read STL_CHECK_LENGTH bytes from file
    unsigned char tmpbuf[STL_CHECK_LENGTH + 1] = {0};

    // GCC does not like that we read beyond end of file, so we limit the
    // number of bytes read
    const size_t maxlen = std::min(static_cast<size_t>(STL_CHECK_LENGTH), size);
    file.read(reinterpret_cast<char*>(tmpbuf), maxlen);
    size_t nbytes = file.gcount();

    size_t nascii = 0;
    for (size_t i = 0; i < nbytes; ++i)
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

    // If *only* encountered ASCII characters, this is most probably an ASCII
    // file
    if (nascii == nbytes)
    {
        return 0;
    }

    // Read number of triangles - 80 bytes beyond file start
    unsigned int ntri;
    file.seekg(80, std::ios_base::beg);
    file.read(reinterpret_cast<char*>(&ntri), sizeof(ntri));

    // If ntri cannot be safely casted to an int MGLET cannot use this STL
    if (ntri > 0x7FFFFFFF)
    {
        std::cerr << "Too many triangles in STL: " << filename << std::endl;
        std::cerr << "  Number of triangles: " << ntri << std::endl;
        std::cerr << "  Exceed MGLET maximum: " << 0x7FFFFFFF << std::endl;
        std::cerr << std::flush;
        return -3;
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
    std::cerr << "Error in detecting ASCII/binary STL: " << filename
              << std::endl;
    std::cerr << "  File length: " << size << std::endl;
    std::cerr << "  Read " << nbytes << ", found " << nascii << " valid chars"
              << std::endl;
    std::cerr << "  If this is binary, it has " << ntri << " triangles"
              << std::endl;
    std::cerr << "  Expected filesize is then " << expected_size
              << std::endl;
    std::cerr << "Now assuming binary mode and hope for the best..."
              << std::endl;
    std::cerr << std::flush;
    return ntri;
}


// Read ASCII STL file
int stl_read_ascii(std::vector<float>& vertices, const std::string& filename)
{
    std::ifstream file(filename);
    if (!file)
    {
        std::cerr << "Could not open: " << filename << std::endl;
        std::cerr << "errno: " << errno << std::endl;
        return -1;
    }

    int lineno = 0;
    std::string line;
    while (std::getline(file, line))
    {
        // Increment line number counter
        ++lineno;

        // Skip leading whitespace
        size_t pos = line.find_first_not_of(" \t");
        if (pos == std::string::npos)
        {
            continue;
        }
        line = line.substr(pos);

        // Skip "facet" and "endfacet" lines
        if (line.substr(0, 5) == "facet" || line.substr(0, 8) == "endfacet")
        {
            continue;
        }

        // Skip "outer loop" and "endloop" lines
        if (line.substr(0, 10) == "outer loop" || line.substr(0, 7) == "endloop")
        {
            continue;
        }

        // We read the information on the "vertex" lines
        if (line.substr(0, 6) == "vertex")
        {
            float x, y, z;

#if __cpp_lib_to_chars >= 201611L
            // C++17 std::from_chars is the fastest way to convert string to
            // float - but we manually has to skip whitespace
            const char* p = line.c_str() + 7;
            const char* endp = line.c_str() + line.size();

            std::from_chars_result ans;

            char* q = skip_whitespace(p, endp);
            ans = std::from_chars(q, endp, x);

            q = skip_whitespace(ans.ptr, endp);
            ans = std::from_chars(q, endp, y);

            q = skip_whitespace(ans.ptr, endp);
            std::from_chars(q, endp, z);
#else
            // Read line using std::sscanf - performance-wise this is in
            // between using std::istringstream and std::from_chars
            if (std::sscanf(line.c_str() + 7, "%f %f %f", &x, &y, &z) != 3)
            {
                std::cerr << "Error reading vertex line in ASCII STL: "
                          << filename << std::endl;
                std::cerr << "  Line number: " << lineno << std::endl;
                std::cerr << "  Line: " << line << std::endl;
                return -2;
            }
#endif
            vertices.push_back(x);
            vertices.push_back(y);
            vertices.push_back(z);

            continue;
        }

        // Empty lines and solid/endolid are more rare - therefore check these
        // in the end

        // Skip empty lines
        if (line.empty())
        {
            continue;
        }

        // Skip "solid" and "endsolid" lines
        if (line.substr(0, 5) == "solid" || line.substr(0, 8) == "endsolid")
        {
            continue;
        }

        // If we reach this point, we have an unexpected line in the file
        std::cerr << "Unexpected line in ASCII STL: " << filename << std::endl;
        std::cerr << "  Line number: " << lineno << std::endl;
        std::cerr << "  Line: " << line << std::endl;
    }

    return 0;
}


// Read binary STL file
int stl_read_binary(std::vector<float>& vertices, const std::string& filename)
{
    std::ifstream file(filename, std::ios::binary);
    if (!file)
    {
        std::cerr << "Could not open: " << filename << std::endl;
        std::cerr << "errno: " << errno << std::endl;
        std::cerr << std::flush;
        return -1;
    }

    // Skip the 80 byte header
    file.seekg(80);

    // Read number of triangles
    unsigned int ntri;
    file.read(reinterpret_cast<char*>(&ntri), sizeof(ntri));

    // Reserve space in vector
    vertices.reserve(9*ntri);

    // Read all triangles
    for (unsigned int i = 0; i < ntri; ++i)
    {
        struct __attribute__((packed)) triangle_t {
            float normal[3];
            float v1[3];
            float v2[3];
            float v3[3];
            unsigned short attribute;
        } triangle;

        file.read(reinterpret_cast<char*>(&triangle), sizeof(triangle));

        size_t nbytes = file.gcount();
        if (nbytes != sizeof(triangle))
        {
            std::cerr << "Error reading binary STL: " << filename << std::endl;
            std::cerr << "  Expected " << sizeof(triangle) << " bytes, got "
                      << nbytes << " bytes" << std::endl;
            return -2;
        }

        vertices.push_back(triangle.v1[0]);
        vertices.push_back(triangle.v1[1]);
        vertices.push_back(triangle.v1[2]);
        vertices.push_back(triangle.v2[0]);
        vertices.push_back(triangle.v2[1]);
        vertices.push_back(triangle.v2[2]);
        vertices.push_back(triangle.v3[0]);
        vertices.push_back(triangle.v3[1]);
        vertices.push_back(triangle.v3[2]);
    }

    return 0;
}
