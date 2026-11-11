#pragma once
#include <cerrno>
#include <fcntl.h>
#include <unistd.h>

// Reading immutable checkpoints must not dirty an inode merely to record its
// access time. In particular, relatime can trigger one Btrfs COW update for
// every old file during a dependency-graph restore. Foreign-owned readable
// files still work: retry without O_NOATIME only when that flag is forbidden.
// This changes neither file contents nor the mtime/ctime integrity checks.
namespace CheckpointRead {
inline int Open(const char* path,int extraFlags=0){
    const int flags=O_RDONLY|O_CLOEXEC|extraFlags;
    int fd;
#ifdef O_NOATIME
    do{fd=open(path,flags|O_NOATIME);}while(fd<0&&errno==EINTR);
    if(fd>=0||errno!=EPERM)return fd;
#endif
    do{fd=open(path,flags);}while(fd<0&&errno==EINTR);
    return fd;
}
inline ssize_t Read(int fd,void* data,size_t size){
    ssize_t result;do{result=read(fd,data,size);}while(result<0&&errno==EINTR);
    return result;
}
struct File {
    int fd;
    explicit File(const char* path,int extraFlags=0):fd(Open(path,extraFlags)){}
    ~File(){if(fd>=0)close(fd);}
    File(const File&)=delete;
    File& operator=(const File&)=delete;
};
}
