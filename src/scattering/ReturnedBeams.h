#pragma once

#include "../Beam.h"
#include <cstring>
#include <vector>

// A traced aperture can contain up to MAX_VERTEX_NUM vertices, but usually
// occupies only a fraction of that array. Keep terminal rays in trace order
// and store their actual vertices contiguously at full geometry precision.
class ReturnedBeams
{
    struct Ray
    {
        Light light;
        Track track;
        Matrix2x2c jones;
        double opticalPath, front;
        int reflections, facet, vertexCount;
        Location location;
        size_t firstVertex;
#ifdef _DEBUG
        std::vector<Polygon> pols;
        std::vector<Vector3f> dirs;
        std::vector<double> ops;
#endif
    };
    std::vector<Ray> rays;
    std::vector<Point3f> vertices;

public:
    void Clear()
    {
        rays.clear();
        vertices.clear();
    }

    size_t size() const { return rays.size(); }

    void Append(const Beam &beam)
    {
        if (beam.nVertices < 0 || beam.nVertices > MAX_VERTEX_NUM)
            throw std::runtime_error("invalid returned aperture vertex count");
        rays.emplace_back();
        Ray &ray = rays.back();
        ray.light = static_cast<const Light &>(beam);
        ray.track.locations = beam.locations;
#ifdef _TRACK_ALLOW
        ray.track.id = beam.id;
#endif
        ray.jones = beam.J;
        ray.opticalPath = beam.opticalPath;
        ray.front = beam.front;
        ray.reflections = beam.nActs;
        ray.facet = beam.lastFacetId;
        ray.location = beam.location;
        ray.vertexCount = beam.nVertices;
        ray.firstVertex = vertices.size();
        vertices.insert(vertices.end(), beam.arr, beam.arr + beam.nVertices);
#ifdef _DEBUG
        ray.pols = beam.pols;
        ray.dirs = beam.dirs;
        ray.ops = beam.ops;
#endif
    }

    // Each iterator reconstructs one mutable Beam for the existing PO
    // preprocessing pass. Changes to that beam stay local to the pass.
    class Iterator
    {
        const ReturnedBeams &storage;
        size_t index;
        Beam beam;
    public:
        Iterator(const ReturnedBeams &s, size_t i) : storage(s), index(i) {}
        Iterator(const Iterator &other) : storage(other.storage), index(other.index) {}
        bool operator!=(const Iterator &other) const { return index != other.index; }
        void operator++() { ++index; }
        Beam &operator*()
        {
            const Ray &ray = storage.rays[index];
            static_cast<Light &>(beam) = ray.light;
            beam.locations = ray.track.locations;
#ifdef _TRACK_ALLOW
            beam.id = ray.track.id;
#endif
            beam.J = ray.jones;
            beam.opticalPath = ray.opticalPath;
            beam.front = ray.front;
            beam.nActs = ray.reflections;
            beam.lastFacetId = ray.facet;
            beam.location = ray.location;
            beam.nVertices = ray.vertexCount;
            if (ray.vertexCount)
                std::memcpy(beam.arr, storage.vertices.data() + ray.firstVertex,
                            ray.vertexCount * sizeof(Point3f));
#ifdef _DEBUG
            beam.pols = ray.pols;
            beam.dirs = ray.dirs;
            beam.ops = ray.ops;
#endif
            return beam;
        }
    };
    Iterator begin() const { return Iterator(*this, 0); }
    Iterator end() const { return Iterator(*this, size()); }
};
