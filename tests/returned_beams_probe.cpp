#include "../src/scattering/ReturnedBeams.h"
#include <cassert>

int main()
{
    ReturnedBeams storage;
    for (int count : {0, 3, MAX_VERTEX_NUM}) {
        Beam ray;
        ray.nVertices = count;
        for (int i = 0; i < count; ++i)
            ray.arr[i] = Point3f(i + .25, -i, count + .75, -i - .5);
        ray.direction = Point3f(0, 0, 1);
        ray.polarizationBasis = Point3f(0, 1, 0);
        ray.J.Fill(complex(count, -count));
        ray.nActs = 7;
        ray.location = Location::Out;
        ray.locations = 65;
        ray.lastFacetId = count + 3;
        ray.id = count + 2;
        ray.opticalPath = 12.25;
        ray.front = -1.5;
        storage.Append(ray);
    }
    for (int pass = 0; pass < 2; ++pass) {
        auto it = storage.begin();
        for (int count : {0, 3, MAX_VERTEX_NUM}) {
            assert(it != storage.end());
            Beam &ray = *it;
            assert(ray.nVertices == count && ray.lastFacetId == count + 3);
            assert(ray.nActs == 7 && ray.location == Location::Out && ray.locations == 65);
#ifdef _TRACK_ALLOW
            assert(ray.id == count + 2);
#endif
            assert(ray.front == -1.5 && ray.opticalPath == 12.25);
            assert(ray.direction.cz == 1 && ray.polarizationBasis.cy == 1);
            assert(abs(ray.J.m12 - complex(count, -count)) == 0);
            for (int i = 0; i < count; ++i) {
                const auto &point = ray.arr[i];
                assert(point.cx == i + .25 && point.cy == -i && point.cz == count + .75);
                assert(point.d_param == -i - .5);
                ray.arr[i] = Point3f(0, 0, 0);
            }
            ray.J.Fill(complex(0, 0));
            ++it;
        }
        assert(!(it != storage.end()));
    }
    storage.Clear();
    assert(storage.size() == 0);
}
