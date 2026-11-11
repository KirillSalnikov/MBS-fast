#include <sys/stat.h>
#include <time.h>
#include <cerrno>
// Deterministic time fixtures. Normal legacy cases see stable files older
// than two seconds without slowing this test with sleeps. Coarse-stamp cases
// independently hold inode timestamps fixed across actual byte mutations.
static bool fixtureStamp=false,fixtureInode=false,fixedClock=false,failedClock=false;
static struct timespec fixtureTime{},fixtureNow{};
static int TestStat(const char* path,struct stat* result){
    const int rc=::stat(path,result);
    if(rc==0&&fixtureStamp){result->st_mtim=fixtureTime;result->st_ctim=fixtureTime;if(fixtureInode)result->st_ino=42;}
    return rc;
}
static int TestClock(clockid_t kind,struct timespec* result){
    if(failedClock){errno=EIO;return -1;}
    if(fixedClock){*result=fixtureNow;return 0;}
    const int rc=::clock_gettime(kind,result);if(rc==0)result->tv_sec+=2;return rc;
}
#define stat(...) TestStat(__VA_ARGS__)
#define clock_gettime(...) TestClock(__VA_ARGS__)
#include "../src/CheckedFileCache.h"
#undef clock_gettime
#undef stat
#include <cassert>
#include <fcntl.h>
#include <unistd.h>
#include <cstdio>
#include <iostream>
static void Write(const std::string& p,const std::string& s){std::ofstream f(p,std::ios::binary);f<<s;}
int main(){
    char temp[]="/tmp/mbs-checked-files-XXXXXX";assert(mkdtemp(temp));const std::string root=temp,p=root+"/a",q=root+"/b";
    // A real Kika ext4/WSL failure: distinct newly created job paths can
    // immediately reuse inode+size+ctime+mtime, even at nanosecond resolution.
    // Full logical paths disambiguate mature aliases; young same-path writes
    // require the timestamp quarantine as well, not only a longer cache key.
    fixtureStamp=true;fixtureInode=true;fixedClock=true;
    fixtureTime={1000,250000000};fixtureNow={1002,250000000};
    {
        CheckedFileCache matureAliases(8192);Write(p,"alpha");Write(q,"bravo");
        assert(matureAliases.Text(p)=="alpha");assert(matureAliases.Text(q)=="bravo");
        const auto n=matureAliases.opens;
        assert(matureAliases.Text(p)=="alpha"&&matureAliases.Text(q)=="bravo"&&matureAliases.opens==n);
        assert(matureAliases.Bytes()<=matureAliases.Limit());
    }
    fixtureNow={1000,254000000};
    {
        CheckedFileCache young(8192);Write(p,"alpha");const auto prior=young.Hash(p);Write(p,"bravo");
        assert(young.Text(p)=="bravo"&&young.Hash(p)!=prior);
        bool rejected=false;try{young.VerifiedText(p,prior);}catch(const std::exception&){rejected=true;}assert(rejected);
        const auto n=young.opens;assert(young.Text(p)=="bravo"&&young.opens==n+1);
        fixtureNow={1001,249999999};young.Hash(p);const auto beforeBoundary=young.opens;young.Hash(p);assert(young.opens==beforeBoundary+1);
        fixtureNow={1001,250000000};young.Hash(p);const auto mature=young.opens;young.Hash(p);assert(young.opens==mature);
        failedClock=true;young.Hash(p);assert(young.opens==mature+1);failedClock=false;
        fixtureNow={999,0};const auto future=young.opens;young.Hash(p);assert(young.opens==future+1);
    }
    // A young value must not become reusable simply because time passes:
    // another write can occur in its coarse tick without any intervening Get.
    fixtureNow={1000,254000000};
    {
        CheckedFileCache delayed(8192);Write(p,"alpha");const auto prior=delayed.Hash(p);
        assert(delayed.Bytes()==0);Write(p,"bravo");fixtureNow={1002,254000000};
        assert(delayed.Text(p)=="bravo"&&delayed.Hash(p)!=prior);
        const auto opened=delayed.opens;assert(delayed.Text(p)=="bravo"&&delayed.opens==opened);
        // A clock move into the young interval must also evict a previously
        // retained value; a subsequent same-stamp write cannot revive it.
        fixtureNow={1000,254000000};assert(delayed.Text(p)=="bravo"&&delayed.Bytes()==0);
        Write(p,"delta");fixtureNow={1002,254000000};assert(delayed.Text(p)=="delta");
    }
    fixtureStamp=fixtureInode=fixedClock=false;unlink(q.c_str());
    CheckedFileCache cache(8192);Write(p,"hello");assert(cache.Hash(p)=="a430d84680aabd0b");
    assert(cache.Text(p)=="hello");const auto opens=cache.opens;
    assert(cache.Text(p)=="hello"&&cache.Hash(p)=="a430d84680aabd0b"&&cache.opens==opens);
    assert(cache.VerifiedText(p,"a430d84680aabd0b")=="hello"&&cache.opens==opens);
    bool wrongDigest=false;try{cache.VerifiedText(p,"wrong");}catch(const std::exception&){wrongDigest=true;}assert(wrongDigest);
    // A same-size edit with restored mtime must still invalidate via ctime.
    struct stat old;assert(stat(p.c_str(),&old)==0);usleep(2000);Write(p,"world");
    const timespec times[2]={old.st_atim,old.st_mtim};assert(utimensat(AT_FDCWD,p.c_str(),times,0)==0);
    assert(cache.Text(p)=="world"&&cache.Hash(p)!="a430d84680aabd0b");
    wrongDigest=false;try{cache.VerifiedText(p,"a430d84680aabd0b");}catch(const std::exception&){wrongDigest=true;}assert(wrongDigest);
    assert(cache.VerifiedText(p,cache.Hash(p))=="world");
    // Atomic replacement, hard links, missing paths, empty/binary/large data.
    Write(q,"hello");assert(rename(q.c_str(),p.c_str())==0);assert(cache.Text(p)=="hello");
    assert(link(p.c_str(),q.c_str())==0);assert(cache.Text(q)=="hello");Write(q,"");assert(cache.Text(p).empty());
    std::string binary(200000,'x');binary[100]='\0';Write(p,binary);assert(cache.Text(p)==binary&&CheckedFileCache::Read(p)==binary);
    const auto large=cache.opens;cache.Text(p);assert(cache.opens==large+1); // large text is not retained
    const auto largeHash=cache.Hash(p);assert(cache.VerifiedText(p,largeHash)==binary);
    const std::pair<size_t,const char*> boundaries[]={
        {0,"cbf29ce484222325"},{1,"af63ec4c860207bc"},{8191,"36f7cfa0bd404bb6"},
        {8192,"a76195219440c325"},{65535,"81b9ec1d3a206bb6"},
        {65536,"8f57fca9c5172325"},{65537,"a9a5a879e65107bc"}};
    for(const auto& sample:boundaries){
        const std::string boundary(sample.first,'q');Write(p,boundary);
        assert(cache.Text(p)==boundary&&cache.Hash(p)==sample.second);
        assert(cache.VerifiedText(p,sample.second)==boundary);
    }
    // procfs advertises a zero-length regular file but exposes bytes on read.
    assert(cache.Text("/proc/self/cmdline")==CheckedFileCache::Read("/proc/self/cmdline"));
    if(access("/sys/devices/system/cpu/online",R_OK)==0)
        assert(cache.Text("/sys/devices/system/cpu/online")==CheckedFileCache::Read("/sys/devices/system/cpu/online"));
    for(int i=0;i<100;++i){const auto f=root+"/"+std::to_string(i);Write(f,std::string(100,'a'+i%20));assert(cache.Text(f).size()==100);assert(cache.Bytes()<=cache.Limit());unlink(f.c_str());}
    const auto marker=root+"/complete.txt";Write(marker,"receipt");cache.Hash(marker);const auto markerOpens=cache.opens;
    assert(cache.Text(marker)=="receipt"&&cache.opens==markerOpens);unlink(marker.c_str());
    unlink(p.c_str());bool missing=false;try{cache.Hash(p);}catch(const std::exception&){missing=true;}assert(missing);
    unlink(q.c_str());rmdir(root.c_str());std::cout<<"checked file cache: freshness, binary data, eviction and memory bound passed\n";
}
