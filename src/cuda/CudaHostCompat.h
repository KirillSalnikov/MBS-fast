#pragma once

// CUDA 12.x declares sinpi/cospi/rsqrt differently from glibc 2.43. Keep
// GNU extensions disabled for nvcc, then restore the two pthread declarations
// expected by libstdc++ builds configured against glibc GNU extensions.
#if defined(__linux__) && defined(__cplusplus)
#include <pthread.h>
#include <time.h>

#if defined(__GLIBC__) && __GLIBC_PREREQ(2, 30) && !defined(__USE_GNU)
extern "C" int pthread_cond_clockwait(
    pthread_cond_t* cond,
    pthread_mutex_t* mutex,
    clockid_t clockId,
    const struct timespec* absoluteTime) noexcept;

extern "C" int pthread_mutex_clocklock(
    pthread_mutex_t* mutex,
    clockid_t clockId,
    const struct timespec* absoluteTime) noexcept;
#endif
#endif
