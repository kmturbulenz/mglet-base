#pragma once

#include <cstdint>

// This is basically duplicating mglet_precision.h from the core library.
// The reason is that otherwise we would need a core library for the core
// library that can be linked to the backend library.

#ifdef _MGLET_INT64_
using mgletint = std::int64_t;
#else
using mgletint = std::int32_t;
#endif

#ifdef _MGLET_IFK64_
using mgletifk = std::int64_t;
#else
using mgletifk = std::int32_t;
#endif

#ifdef _MGLET_DOUBLE_PRECISION_
using mgletreal = double;
#else
using mgletreal = float;
#endif
