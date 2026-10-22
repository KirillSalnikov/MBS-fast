#pragma once

// Controller-only upgrades may retain an older native solver. A registered CLI
// option is not proof that main/RunConfig/the tracer implement its products.
// Enable this only for a build whose complete native prefix path was verified.
#ifndef MBS_NATIVE_ORIENTATION_PREFIX_OUTPUT
#define MBS_NATIVE_ORIENTATION_PREFIX_OUTPUT 0
#endif
#if MBS_NATIVE_ORIENTATION_PREFIX_OUTPUT != 0 && MBS_NATIVE_ORIENTATION_PREFIX_OUTPUT != 1
#error MBS_NATIVE_ORIENTATION_PREFIX_OUTPUT must be 0 or 1
#endif
namespace NativeWorkerCapabilities {
constexpr bool OrientationPrefixOutput(){return MBS_NATIVE_ORIENTATION_PREFIX_OUTPUT != 0;}
}
