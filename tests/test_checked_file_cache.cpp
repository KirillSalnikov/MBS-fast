#include "../src/CheckedFileCache.h"
#include <cassert>
#include <fcntl.h>
#include <unistd.h>
#include <cstdio>
#include <iostream>
static void Write(const std::string& p,const std::string& s){std::ofstream f(p,std::ios::binary);f<<s;}
int main(){
    char temp[]="/tmp/mbs-checked-files-XXXXXX";assert(mkdtemp(temp));const std::string root=temp,p=root+"/a",q=root+"/b";
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
    for(int i=0;i<100;++i){const auto f=root+"/"+std::to_string(i);Write(f,std::string(100,'a'+i%20));assert(cache.Text(f).size()==100);assert(cache.Bytes()<=cache.Limit());unlink(f.c_str());}
    const auto marker=root+"/complete.txt";Write(marker,"receipt");cache.Hash(marker);const auto markerOpens=cache.opens;
    assert(cache.Text(marker)=="receipt"&&cache.opens==markerOpens);unlink(marker.c_str());
    unlink(p.c_str());bool missing=false;try{cache.Hash(p);}catch(const std::exception&){missing=true;}assert(missing);
    unlink(q.c_str());rmdir(root.c_str());std::cout<<"checked file cache: freshness, binary data, eviction and memory bound passed\n";
}
