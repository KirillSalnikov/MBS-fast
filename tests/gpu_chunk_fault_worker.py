#!/usr/bin/env python3
"""Test-only status injection. Successful requests execute the real CPU worker."""
import json,os,socket,struct,sys
from pathlib import Path
real=os.environ['TASK_REAL_WORKER'];marker=Path(os.environ['TASK_GPU_RETRY_MARKER'])
sequence=[x for x in os.environ.get('TASK_GPU_RETRY_SEQUENCE','').split(',') if x]
used=int(marker.read_text()) if marker.exists() else 0
if len(sys.argv)!=3 or sys.argv[1]!='--fullauto-worker-server' or used>=len(sequence):
    os.execv(real,[real]+sys.argv[1:])
marker.write_text(str(used+1));fault=sequence[used];channel=socket.socket(fileno=int(sys.argv[2]))
def receive(n):
    data=b''
    while len(data)<n:
        block=channel.recv(n-len(data))
        if not block:raise SystemExit(2)
        data+=block
    return data
def number(fmt):return struct.unpack(fmt,receive(struct.calcsize(fmt)))
count,=number('<I');words=[receive(number('<I')[0]).decode() for _ in range(count)]
stream=words[0]=='MBS_WORKER_MEMORY_REQUEST_V2'
if stream:
    root=words[1];request=words[3:]
    magic,version,count,total=number('<IIIQ');assert magic==0x4d425353 and version==2
    for _ in range(count):
        sizepath,size,digest,descriptor=number('<IQQQ');receive(sizepath)
        while size:
            block=min(size,65536);receive(block);size-=block
else:request=words
with marker.with_suffix('.requests.jsonl').open('a') as record:record.write(json.dumps({'pid':os.getpid(),'status':fault,'words':request})+'\n')
log=request[0];folder=str(Path(log).parent)
files={log:f'Injected typed GPU chunk retry status {fault}\n'.encode(),folder+'/result/result.dat':b'POISON INCOMPLETE MATRIX MUST NOT BE ACCEPTED\n'}
if not stream:
    for path,data in files.items():p=Path(path);p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data)
if fault=='eof':channel.close();raise SystemExit(0)
if fault=='partial76':channel.sendall(struct.pack('<i',76));channel.close();raise SystemExit(0)
channel.sendall(struct.pack('<i',int(fault)))
if stream:
    def fnv(data,seed=14695981039346656037):
        for byte in data:seed=((seed^byte)*1099511628211)&0xffffffffffffffff
        return seed
    members=[(path[len(root)+1:].encode(),data) for path,data in files.items()]
    channel.sendall(struct.pack('<IIIQ',0x4d425353,2,len(members),sum(len(path)+len(data) for path,data in members)))
    for path,data in members:
        digest=fnv(data);descriptor=fnv(struct.pack('<QQ',len(data),digest),fnv(b'\0',fnv(path)))
        channel.sendall(struct.pack('<IQQQ',len(path),len(data),digest,descriptor)+path+data)
channel.close()
