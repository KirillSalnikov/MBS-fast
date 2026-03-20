#pragma once

#include <vector>
#include "compl.hpp"
#include "JonesMatrix.h"
#include "Handler.h"

/// A single cached beam: all geometry normalized by D_ref so that
/// multi-size recomputation only needs to rescale phases.
struct CachedBeam {
    // Precomputed 2D aperture vertices (in aperture coordinate system, already centered)
    // These are in absolute units at tracing time; we normalize by D_ref below.
    double vx_norm[32], vy_norm[32]; // 2D vertex coords / D_ref
    double slope_yx[32], slope_xy[32]; // edge slopes (dimensionless, scale-invariant)
    double intercept_y_norm[32], intercept_x_norm[32]; // intercepts / D_ref
    bool edge_valid_x[32], edge_valid_y[32];
    int nVertices;
    bool edgeDataValid;

    // Scale-invariant optical data
    Matrix2x2c J;              // Jones matrix (Fresnel coefficients, angle-dependent only)
    double opticalPath_norm;   // optical path / D_ref
    double projLength_norm;    // projected length / D_ref
    double area_norm;          // area / D_ref^2
    bool isExternal;           // lastFacetId == INT_MAX (shadow beam)

    // Precomputed polarization data (direction-independent)
    BeamPolData polData;

    // Beam direction (unit vector, scale-invariant)
    double dirx, diry, dirz;

    // Aperture coordinate axes (unit vectors, scale-invariant)
    double horAx, horAy, horAz;
    double verAx, verAy, verAz;
    double normx, normy, normz;
    double cenx_norm, ceny_norm, cenz_norm; // center / D_ref
};

/// All beams from one orientation
struct OrientationBeams {
    std::vector<CachedBeam> beams;
    double weight; // orientation weight (1/N or sin(beta)/sum)
    double incomingEnergy; // incident energy contribution for this orientation
};

/// The full beam cache: trace once, reuse for multiple sizes
class BeamCache {
public:
    std::vector<OrientationBeams> orientations;
    double D_ref; // reference diameter (= x_ref * wavelength / PI) used during tracing

    void clear() { orientations.clear(); }
    size_t totalBeams() const {
        size_t n = 0;
        for (const auto& o : orientations) n += o.beams.size();
        return n;
    }
};
