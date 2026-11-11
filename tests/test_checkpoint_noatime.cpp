#include "../src/CheckedFileCache.h"
#include <cassert>
#include <iostream>

static void Put(const std::string& path,const std::string& bytes){
    std::ofstream out(path,std::ios::binary);out.write(bytes.data(),bytes.size());out.close();assert(out);
}
static struct stat Cold(const std::string& path){
    const timespec times[2]={{946684800,123},{978307200,456}};
    assert(utimensat(AT_FDCWD,path.c_str(),times,0)==0);
    struct stat value{};assert(stat(path.c_str(),&value)==0);return value;
}
static void Unchanged(const std::string& path,const struct stat& before){
    struct stat after{};assert(stat(path.c_str(),&after)==0);
    assert(before.st_atim.tv_sec==after.st_atim.tv_sec&&before.st_atim.tv_nsec==after.st_atim.tv_nsec);
    assert(PackedJobStore::Same(before,after));
}
template<class F>static bool Reject(F f){try{f();return false;}catch(const std::exception&){return true;}}
int main(){
    char temporary[]="/tmp/mbs-noatime-XXXXXX";assert(mkdtemp(temporary));const std::string root=temporary;
    std::string bytes(131073,'x');bytes[0]='\0';bytes[65536]='\xff';bytes.back()='\n';
    const auto path=root+"/data";Put(path,bytes);const auto before=Cold(path);
    CheckedFileCache cache(8192);const auto digest=PackedJobStore::Hex(PackedJobStore::Digest(bytes));
    assert(cache.Hash(path)==digest&&cache.Text(path)==bytes&&cache.VerifiedText(path,digest)==bytes);
    assert(CheckedFileCache::Read(path)==bytes);Unchanged(path,before);
    {CheckpointRead::File file(path.c_str());assert(file.fd>=0);
#ifdef O_NOATIME
        assert(fcntl(file.fd,F_GETFL)&O_NOATIME);
#endif
        assert(fcntl(file.fd,F_GETFD)&FD_CLOEXEC);
    }
    // A real foreign-owned readable file exercises EPERM, not an artificial
    // success path. Reading must still work without requesting extra privilege.
    struct stat foreign{};assert(stat("/etc/passwd",&foreign)==0);
    if(geteuid()!=foreign.st_uid){
#ifdef O_NOATIME
        int forbidden=open("/etc/passwd",O_RDONLY|O_NOATIME);assert(forbidden<0&&errno==EPERM);
#endif
        CheckpointRead::File file("/etc/passwd");assert(file.fd>=0);char byte;
        assert(CheckpointRead::Read(file.fd,&byte,1)==1);
#ifdef O_NOATIME
        assert(!(fcntl(file.fd,F_GETFL)&O_NOATIME));
#endif
        assert(!CheckedFileCache::Read("/etc/passwd").empty());
    }
    // Missing file, read-on-directory and O_NOFOLLOW errors remain failures.
    errno=0;assert(CheckpointRead::Open((root+"/missing").c_str())<0&&errno==ENOENT);
    assert(Reject([&]{CheckedFileCache::Read(root+"/missing");}));
    assert(Reject([&]{CheckedFileCache::Read(root);}));
    assert(symlink(path.c_str(),(root+"/link").c_str())==0);
    errno=0;assert(CheckpointRead::Open((root+"/link").c_str(),O_NOFOLLOW)<0&&errno==ELOOP);
    // Same-size modification with restored mtime must still fail old hashes.
    usleep(2000);bytes[1]='z';Put(path,bytes);const timespec restored[2]={before.st_atim,before.st_mtim};
    assert(utimensat(AT_FDCWD,path.c_str(),restored,0)==0);
    assert(Reject([&]{cache.VerifiedText(path,digest);}));
    assert(cache.Text(path)==bytes);
    // Archive index/member reads also preserve old atime and still authenticate
    // the exact original bytes through the ordinary raw-first logical paths.
    assert(mkdir((root+"/jobs").c_str(),0700)==0);const auto job=root+"/jobs/j123";
    assert(mkdir(job.c_str(),0700)==0);Put(job+"/complete.txt","receipt\n");Put(job+"/data",bytes);
    PackedJobStore::Pack(job);const auto archive=PackedJobStore::Shard(root,"j123");const auto packBefore=Cold(archive);
    CheckedFileCache packed(8192);assert(packed.Text(job+"/data")==bytes);
    assert(CheckedFileCache::Read(job+"/complete.txt")=="receipt\n");Unchanged(archive,packBefore);
    std::cout<<"PASS checkpoint noatime: raw/packed bytes, timestamps, foreign-owner EPERM fallback, errors, ctime corruption guard\n";
}
