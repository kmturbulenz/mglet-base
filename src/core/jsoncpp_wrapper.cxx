#include "jsoncpp_wrapper.h"
#define JSON_USE_IMPLICIT_CONVERSIONS 0
#include <nlohmann/json.hpp>

#include <string>
#include <fstream>
#include <iostream>
#include <type_traits>

#include <ISO_Fortran_binding.h>

using json = nlohmann::json;

struct jsoncppc {
    void *obj;
    int refcounter;
};

// Dirty hack for Intel compiler on Centos 7
//
// When using the devtoolset's with GNU compilers, we have access to the
// last header files with C++17 and the corresponding libraries. When using
// the Intel compilers on Centos 7 it uses the system headers, which is
// from GNU 4.8... These does not support C++17. Therefore this hack is needed.
//
// The warning "inline variables are a C++17 extension" can be silenced
// with the flag -Wno-c++17-extensions
#if __cplusplus < 201703
//template< class T, class U >
//inline constexpr bool is_same_v = std::is_same<T, U>::value;
#else
//using std :: is_same_v;
#endif

// https://nachtimwald.com/2017/08/18/wrapping-c-objects-in-c/
extern "C" {
    jsoncppc_t* json_from_file(const char* filename, int* ierr) {
        *ierr = 0;

        // Allocate and initialize new json object
        jsoncppc_t* jsonc = (typeof(jsonc))malloc(sizeof(*jsonc));
        json* obj = new json();

        // Parse file
        std::ifstream fh(filename);
        try
        {
            //fh >> *obj;
            *obj = json::parse(fh, nullptr, true, true);
        }
        catch (json::exception& e)
        {
            std::cerr << "JSON error:\n";
            std::cerr << e.what() << '\n';
            *ierr = e.id;
        }

        jsonc->obj = obj;
        jsonc->refcounter = 1;
        return jsonc;
    }

    jsoncppc_t* json_from_json(jsoncppc_t* jsonc, const char* key, int* ierr) {
        *ierr = 0;

        // Parent json object
        json* obj = static_cast<json *>(jsonc->obj);

        // Allocate new json object
        jsoncppc_t* jsonc_new = (typeof(jsonc))malloc(sizeof(*jsonc));

        // Fetch value and set to pointer
        try
        {
            json::json_pointer json_ptr(key);
            // The below line that is commented out is to make a deep copy
            // json* value = new json(obj->at(json_ptr));
            json* value = &obj->at(json_ptr);
            jsonc_new->obj = value;
            jsonc_new->refcounter = jsonc->refcounter + 1;
        }
        catch (json::exception& e)
        {
            std::cerr << "JSON error:\n";
            std::cerr << e.what() << '\n';
            *ierr = e.id;
        }

        return jsonc_new;
    }

    void json_free(jsoncppc_t* jsonc, int* ierr) {
        *ierr = 0;
        if (jsonc == NULL) {
            *ierr = 1;
            return;
        }

        jsonc->refcounter = jsonc->refcounter - 1;
        if (jsonc->refcounter < 1) {
            delete static_cast<json *>(jsonc->obj);
        }

        free(jsonc);
    }

    void json_dump(jsoncppc_t* jsonc, CFI_cdesc_t* res, int* ierr) {
        *ierr = 0;
        if (jsonc == NULL || res == NULL) {
            std::cerr << "FATAL JSON error - NULL given\n";
            *ierr = 1;
            return;
        }
        json* obj = static_cast<json *>(jsonc->obj);

        // Dump JSON to C++ string, check length
        std::string jsondump = obj->dump(4, ' ', true);
        const CFI_index_t strlen = jsondump.length();  // This is w.o NULL

        // Sanity checks of passed Fortran structure
        assert (res->rank == 1);
        assert (res->type == CFI_type_char);
        assert (res->elem_len == 1);
        assert (res->attribute == CFI_attribute_allocatable);

        // If already allocated - free
        if (res->base_addr) CFI_deallocate(res);

        // Allocate to length of string
        const CFI_index_t lower_bounds[1] = {1};  // Remember this is Fortran
        const CFI_index_t upper_bounds[1] = {strlen};
        *ierr = CFI_allocate(res, lower_bounds, upper_bounds, 1);
        if (*ierr) {
            std::cerr << "Error allocating string: " << strlen << std::endl;
            return;
        }

        // Copy string into allocated memory
        strncpy((char*)res->base_addr, jsondump.c_str(), strlen);
    }

    void json_get_int(jsoncppc_t* jsonc, const char* key,
            int* val, int* ierr) {
        json_get_number(jsonc, key, val, ierr);
    }

    void json_set_int(jsoncppc_t* jsonc, const char* key,
            const int* val, int* ierr) {
        json_set_number(jsonc, key, val, ierr);
    }

    void json_get_int64(jsoncppc_t* jsonc, const char* key,
            int64_t* val, int* ierr) {
        json_get_number(jsonc, key, val, ierr);
    }

    void json_get_float(jsoncppc_t* jsonc, const char* key,
            float* val, int* ierr) {
        json_get_number(jsonc, key, val, ierr);
    }

    void json_get_double(jsoncppc_t* jsonc, const char* key,
            double* val, int* ierr) {
        json_get_number(jsonc, key, val, ierr);
    }

    void json_get_bool(jsoncppc_t* jsonc, const char* key,
            _Bool* val, int* ierr) {
        json_get_number(jsonc, key, val, ierr);
    }

    void json_get_int_arr(jsoncppc_t* jsonc, const char* key,
            int* arr, const size_t length, int* ierr) {
        json_get_arr(jsonc, key, arr, length, ierr);
    }

    void json_get_int64_arr(jsoncppc_t* jsonc, const char* key,
            int64_t* arr, const size_t length, int* ierr) {
        json_get_arr(jsonc, key, arr, length, ierr);
    }

    void json_get_float_arr(jsoncppc_t* jsonc, const char* key,
            float* arr, const size_t length, int* ierr) {
        json_get_arr(jsonc, key, arr, length, ierr);
    }

    void json_get_double_arr(jsoncppc_t* jsonc, const char* key,
            double* arr, const size_t length, int* ierr) {
        json_get_arr(jsonc, key, arr, length, ierr);
    }

    void json_get_char(jsoncppc_t* jsonc, const char* key,
            char* cval, const size_t maxlen, size_t* length, int* ierr) {
        std::string value;
        json_get_number(jsonc, key, &value, ierr);
        if (*ierr == 0) {
            *length = value.length();
            if (*length > maxlen) {
                std::cerr << "Reading string from key: " << key << "\n";
                std::cerr << "String too long, got: " << value << "\n";
                std::cerr << "Length: " << *length
                          << " maxlen: " << maxlen << "\n";
                *ierr = 1;
                return;
            }
            std::char_traits<char>::copy(cval, value.c_str(), *length);
        }
    }

    void json_get_size(jsoncppc_t* jsonc, const char* key,
            size_t* size, int* ierr) {
        *ierr = 0;
        if (jsonc == NULL || key == NULL || size == NULL || ierr == NULL) {
            std::cerr << "FATAL JSON error - NULL given\n";
            *ierr = 1;
            return;
        }

        json* obj = static_cast<json *>(jsonc->obj);
        try
        {
            json::json_pointer json_ptr(key);
            json value = obj->at(json_ptr);
            *size = value.size();
        }
        catch (json::exception& e)
        {
            std::cerr << "JSON error:\n";
            std::cerr << e.what() << '\n';
            *ierr = e.id;
        }
    }


    void json_exists(jsoncppc_t* jsonc, const char* key,
            _Bool* exists, int* type, int* ierr) {

        *ierr = 0;
        *exists = false;
        *type = -1;
        if (jsonc == NULL || key == NULL || exists == NULL || ierr == NULL) {
            std::cerr << "FATAL JSON error - NULL given\n";
            *ierr = 1;
            return;
        }

        json* obj = static_cast<json *>(jsonc->obj);
        try
        {
            json::json_pointer json_ptr(key);
            json value = obj->at(json_ptr);
            *exists = true;
            *type = (int)value.type();
        }
        catch (json::out_of_range& e)
        {
            // Do nothing - false is default value
        }
        catch (json::exception& e)
        {
            std::cerr << "JSON error:\n";
            std::cerr << e.what() << '\n';
            *ierr = e.id;
        }
    }
}


template<typename T>
void json_get_number(jsoncppc_t* jsonc, const char* key, T* val, int* ierr) {
    *ierr = 0;
    if (jsonc == NULL || key == NULL || val == NULL || ierr == NULL) {
        std::cerr << "FATAL JSON error - NULL given" << std::endl;
        *ierr = 1;
        return;
    }

    json* obj = static_cast<json *>(jsonc->obj);

    try
    {
        json::json_pointer json_ptr(key);
        // *val = (*obj)[json_ptr];
        // *val = obj->operator[](json_ptr);
        // *val = obj->at(json_ptr).get<int>();
        json value = obj->at(json_ptr);

        if ((std::is_same<T, int64_t>::value || std::is_same<T, int>::value)
                && !value.is_number_integer()) {
            std::cerr << "Invalid datatype for key: " << key << "\n";
            std::cerr << "Expected integer, got: " << value << "\n";
            *ierr = 1;
            return;
        }
        else if ((std::is_same<T, float>::value ||
                std::is_same<T, double>::value)
                && !value.is_number_float()) {
            std::cerr << "Invalid datatype for key: " << key << "\n";
            std::cerr << "Expected real, got: " << value << "\n";
            *ierr = 1;
            return;
        }
        else if ((std::is_same<T, bool>::value ||
                std::is_same<T, _Bool>::value)
                && !value.is_boolean()) {
            std::cerr << "Invalid datatype for key: " << key << "\n";
            std::cerr << "Expected logical, got: " << value << "\n";
            *ierr = 1;
            return;
        }
        else if (std::is_same<T, std::string>::value && !value.is_string()) {
            std::cerr << "Invalid datatype for key: " << key << "\n";
            std::cerr << "Expected string, got: " << value << "\n";
            *ierr = 1;
            return;
        }
        *val = value.get<T>();
    }
    catch (json::out_of_range& e)
    {
        if (e.id == 403) {
            // Key does not exist
            *ierr = -1;
        }
        else {
            std::cerr << "JSON error:\n";
            std::cerr << e.what() << '\n';
            *ierr = e.id;
        }
    }
    catch (json::exception& e)
    {
        std::cerr << "JSON error:\n";
        std::cerr << e.what() << '\n';
        *ierr = e.id;
    }
}


template<typename T>
void json_set_number(jsoncppc_t* jsonc, const char* key,
        const T* val, int* ierr) {
    *ierr = 0;
    if (jsonc == NULL || key == NULL || val == NULL || ierr == NULL) {
        std::cerr << "FATAL JSON error - NULL given\n";
        *ierr = 1;
        return;
    }

    json* obj = static_cast<json *>(jsonc->obj);

    // Check is entry exists
    _Bool exists;
    int type;
    json_exists(jsonc, key, &exists, &type, ierr);
    if (*ierr != 0) {
        return;
    }

    json::json_pointer json_ptr(key);

    // If entry exists - update it - but only if type match
    // (i.e. cannot set an integer to a float)
    if (exists) {
        json value = obj->at(json_ptr);

        if ((std::is_same<T, int64_t>::value || std::is_same<T, int>::value)
                && !value.is_number_integer()) {
            std::cerr << "Invalid datatype for key: " << key << "\n";
            std::cerr << "Expected integer, got: " << value << "\n";
            *ierr = 1;
            return;
        }
        else if ((std::is_same<T, float>::value ||
                std::is_same<T, double>::value)
                && !value.is_number_float()) {
            std::cerr << "Invalid datatype for key: " << key << "\n";
            std::cerr << "Expected real, got: " << value << "\n";
            *ierr = 1;
            return;
        }
        else if ((std::is_same<T, bool>::value ||
                std::is_same<T, _Bool>::value)
                && !value.is_boolean()) {
            std::cerr << "Invalid datatype for key: " << key << "\n";
            std::cerr << "Expected logical, got: " << value << "\n";
            *ierr = 1;
            return;
        }
        else if (std::is_same<T, std::string>::value && !value.is_string()) {
            std::cerr << "Invalid datatype for key: " << key << "\n";
            std::cerr << "Expected string, got: " << value << "\n";
            *ierr = 1;
            return;
        }

        try
        {
            value = *val;
        }
        catch (json::exception& e)
        {
            std::cerr << "JSON error:\n";
            std::cerr << e.what() << '\n';
            *ierr = e.id;
        }
    }
    // Create it
    else {
        try
        {
            obj->operator[](json_ptr) = *val;
        }
        catch (json::exception& e)
        {
            std::cerr << "JSON error:\n";
            std::cerr << e.what() << '\n';
            *ierr = e.id;
        }
    }
}


template<typename T>
void json_get_arr(jsoncppc_t* jsonc, const char* key,
        T* arr, const size_t length, int* ierr) {
    *ierr = 0;
    if (jsonc == NULL || key == NULL || arr == NULL || ierr == NULL) {
        std::cerr << "FATAL JSON error - NULL given\n";
        *ierr = 1;
        return;
    }

    json* obj = static_cast<json *>(jsonc->obj);
    try
    {
        json::json_pointer json_ptr(key);
        json value = obj->at(json_ptr);
        if (!value.is_array()) {
            std::cerr << "Invalid datatype for key: " << key << "\n";
            std::cerr << "Expected array, got: " << value << "\n";
            *ierr = 1;
            return;
        }

        size_t count = 0;
        for (auto elem: value) {
            if (count >= length) {
                *ierr = -1;
                break;
            }
            arr[count++] = (T)elem;
        }
    }
    catch (json::exception& e)
    {
        std::cerr << "JSON error:\n";
        std::cerr << e.what() << '\n';
        *ierr = e.id;
    }
}
