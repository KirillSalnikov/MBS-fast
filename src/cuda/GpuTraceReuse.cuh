// Bounded device-resident raw trace cache. The exact input bytes include the
// rotated geometry, initial fields, paths, depth and trace cutoff. Observation
// theta/phi and Jones-to-Mueller accumulation are deliberately recomputed.
namespace {
struct TraceCacheEntry {
    std::string key;int device=0;size_t bytes=0,count=0;
    Input* inputs=nullptr;Output* output=nullptr;std::vector<Status> status;
    ~TraceCacheEntry(){if(inputs||output){int before=0;cudaGetDevice(&before);cudaSetDevice(device);cudaFree(inputs);cudaFree(output);cudaSetDevice(before);}}
};
struct TraceCache {
    std::mutex mutex;size_t limit=0,bytes=0;
    std::list<std::unique_ptr<TraceCacheEntry>> entries;
};
TraceCache& RawCache(){static TraceCache cache;return cache;}
}
void SetCacheLimit(size_t bytes){auto& cache=RawCache();std::lock_guard<std::mutex> lock(cache.mutex);cache.limit=bytes;
    while(cache.bytes>bytes&&!cache.entries.empty()){cache.bytes-=cache.entries.back()->bytes;cache.entries.pop_back();}}
