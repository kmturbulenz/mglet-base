#include <stdint.h>

#if defined(_MGLET_PROFILE_ROCTX_)
#include <rocprofiler-sdk-roctx/roctx.h>
#elif defined(_MGLET_PROFILE_NVTX_)
#include <nvtx3/nvToolsExt.h>

static const uint32_t colors[] = {
    0xff4e79a7,
    0xfff28e2b,
    0xffe15759,
    0xff76b7b2,
    0xff59a14f,
    0xffedc948,
    0xffb07aa1,
    0xffff9da7,
    0xff9c755f,
    0xffbab0ab,
    0xff86bcb6,
    0xffd37295,
    0xff8a8acb
};

static const uint32_t num_colors = sizeof(colors) / sizeof(colors[0]);

static void nvtx_range_push_colored(const char *name, uint32_t color_id) {
    color_id %= num_colors;

    nvtxEventAttributes_t event_attrib = {
        .version = NVTX_VERSION,
        .size = NVTX_EVENT_ATTRIB_STRUCT_SIZE,
        .colorType = NVTX_COLOR_ARGB,
        .color = colors[color_id],
        .messageType = NVTX_MESSAGE_TYPE_ASCII,
        .message.ascii = name
    };

    nvtxRangePushEx(&event_attrib);
}

// Uses the FNV-1a hash algorithm to generate a color ID from a string
// https://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function
static uint32_t color_id_from_name(const char *s) {
    const uint32_t fnv_offset_basis = 0x811c9dc5u;
    const uint32_t fnv_prime = 0x01000193u;
    uint32_t hash = fnv_offset_basis;

    while (*s != '\0') {
        hash ^= (unsigned char)*s++;
        hash *= fnv_prime;
    }

    return hash % num_colors;
}
#endif

void profile_range_push_impl(const char *msg) {
#if defined(_MGLET_PROFILE_NVTX_)
    nvtx_range_push_colored(msg, color_id_from_name(msg));
#elif defined(_MGLET_PROFILE_ROCTX_)
    roctxRangePushA(msg);
#else
    (void)msg;
#endif
}

void profile_range_pop_impl() {
#if defined(_MGLET_PROFILE_NVTX_)
    nvtxRangePop();
#elif defined(_MGLET_PROFILE_ROCTX_)
    roctxRangePop();
#endif
}

void profile_mark_impl(const char *msg) {
#if defined(_MGLET_PROFILE_NVTX_)
    nvtxMarkA(msg);
#elif defined(_MGLET_PROFILE_ROCTX_)
    roctxMarkA(msg);
#else
    (void)msg;
#endif
}
