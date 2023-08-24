#ifndef __JSONCPP_C_H__
#define __JSONCPP_C_H__

#include <cstddef>
#include <cstdbool>
#include <cstdint>

#include <ISO_Fortran_binding.h>

#ifdef __cplusplus
extern "C" {
#endif

struct jsoncppc;
typedef struct jsoncppc jsoncppc_t;

jsoncppc_t* json_from_file(const char*, int*);
jsoncppc_t* json_from_json(jsoncppc_t*, const char*, int*);
void json_free(jsoncppc_t*, int*);
void json_dump(jsoncppc_t*, CFI_cdesc_t*, int*);

void json_get_int(jsoncppc_t*, const char*, int*, int*);
void json_set_int(jsoncppc_t*, const char*, const int*, int*);
void json_get_int64(jsoncppc_t*, const char*, int64_t*, int*);
void json_set_int64(jsoncppc_t*, const char*, const int64_t*, int*);
void json_get_float(jsoncppc_t*, const char*, float*, int*);
void json_set_float(jsoncppc_t*, const char*, const float*, int*);
void json_get_double(jsoncppc_t*, const char*, double*, int*);
void json_set_double(jsoncppc_t*, const char*, const double*, int*);
void json_get_bool(jsoncppc_t*, const char*, _Bool*, int*);
void json_set_bool(jsoncppc_t*, const char*, const _Bool*, int*);
void json_get_char(jsoncppc_t*, const char*, char*, const size_t,
    size_t*, int*);
void json_set_char(jsoncppc_t*, const char*, const char*, int*);

void json_get_int_arr(jsoncppc_t*, const char*, int*, const size_t, int*);
void json_get_int64_arr(jsoncppc_t*, const char*, int64_t*, const size_t, int*);
void json_get_float_arr(jsoncppc_t*, const char*, float*, const size_t, int*);
void json_get_double_arr(jsoncppc_t*, const char*, double*, const size_t, int*);

void json_get_size(jsoncppc_t*, const char*, size_t* size, int*);
void json_exists(jsoncppc_t*, const char*, _Bool*, int*, int*);

#ifdef __cplusplus
}
#endif


template<typename T>
void json_get_number(jsoncppc_t*, const char*, T*, int*);

template<typename T>
void json_set_number(jsoncppc_t*, const char*, const T*, int*);

template<typename T>
void json_get_arr(jsoncppc_t*, const char*, T*, const size_t, int*);

#endif /* __JSONCPP_C_H__ */
