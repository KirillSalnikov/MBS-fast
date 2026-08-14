#pragma once

#include "geometry_lib.h"

#include <stdexcept>

/**
 * @brief The Polygon struct
 * Convex polygon
 */
class Polygon
{
public:
    Point3f arr[MAX_VERTEX_NUM];
    int nVertices = 0;

    Polygon() {}
    explicit Polygon(int vertexCount) : nVertices(vertexCount)
    {
        if (vertexCount < 0 || vertexCount > MAX_VERTEX_NUM)
            throw std::out_of_range("invalid polygon vertex count");
    }
    Polygon(const Polygon &other);
    Polygon(Polygon &&other);

    void AddVertex(const Point3f &v);

    void Concat(const Polygon &other)
    {
        for (int i = 0; i < other.nVertices; ++i)
        {
            AddVertex(other.arr[i]);
        }
    }

    void InsertVertex(int index, const Point3f &v)
    {
        if (nVertices >= MAX_VERTEX_NUM)
            throw std::runtime_error(
                "polygon exceeds the 64-vertex geometry limit");
        if (index < 0 || index > nVertices)
            throw std::out_of_range("polygon insertion index is out of range");
        ++nVertices;

        for (int i = nVertices-1; i > index; --i)
        {
            arr[i] = arr[i-1];
        }

        arr[index] = v;
    }

    void DeleteVertex(int index)
    {
        if (index < 0 || index >= nVertices)
            throw std::out_of_range("polygon deletion index is out of range");
        for (int i = index; i + 1 < nVertices; ++i)
        {
            arr[i] = arr[i+1];
        }

        --nVertices;
    }

    static void InverseVertexOrder(Polygon &polygon)
    {
        Point3f tmp;
        int size = polygon.nVertices-1;

        for (int vi = 0; vi <= size/2; ++vi)
        {
            tmp = polygon.arr[vi];
            polygon.arr[vi] = polygon.arr[size-vi];
            polygon.arr[size-vi] = tmp;
        }
    }

    Polygon & operator = (const Polygon &other);
    Polygon & operator = (Polygon &&other);
    friend std::ostream & operator << (std::ostream &os, const Polygon &beam);

    void Clear()
    {
        nVertices = 0;
    }

    double Area() const;
    Point3f Center() const;
    Point3f Normal() const;
};

class Polygon512
{
public:
    Point3f arr[512];
    size_t nVertices = 0;

    Polygon512() {}

    void AddVertex(const Point3f &v)
    {
        if (nVertices >= 512)
            throw std::runtime_error(
                "temporary polygon exceeds the 512-vertex geometry limit");
        arr[nVertices++] = v;
    }

    void Concat(const Polygon &other)
    {
        for (int i = 0; i < other.nVertices; ++i)
        {
            AddVertex(other.arr[i]);
        }
    }
};

class PolygonArray
{
public:
    Polygon arr[MAX_POLYGON_NUM];
    size_t size = 0;

    void Push(const Polygon &p)
    {
        if (size >= MAX_POLYGON_NUM)
        {
            throw std::runtime_error(
                "polygon clipping produced more than 512 disjoint pieces; "
                "increase MAX_POLYGON_NUM or simplify the particle geometry");
        }
        arr[size++] = p;
    }

    Polygon &Pop()
    {
        if (size == 0)
            throw std::out_of_range("cannot pop an empty polygon array");
        return arr[--size];
    }

    void Clear()
    {
        size = 0;
    }
};
