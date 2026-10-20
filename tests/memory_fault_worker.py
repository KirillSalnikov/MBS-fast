#!/usr/bin/env python3
"""Protocol fault injection only: production binary performs successful work."""
import os,sys,socket,struct
from pathlib import Path
binary=os.environ['TASK_REAL_WORKER'];mode=os.environ.get('TASK_FAULT','none');marker=Path(os.environ['TASK_FAULT_MARKER'])
if len(sys.argv)==3 and sys.argv[1]=='--fullauto-worker-server' and mode!='none' and (mode=='always' or not marker.exists()):
 marker.write_text('injected\n');s=socket.socket(fileno=int(sys.argv[2]))
 def recv(n):
  b=b''
  while len(b)<n:
   a=s.recv(n-len(b))
   if not a:raise SystemExit(2)
   b+=a
  return b
 n=struct.unpack('I',recv(4))[0];words=[recv(struct.unpack('I',recv(4))[0]).decode() for _ in range(n)]
 Path(words[0]).write_text('Injected resource failure for recovery test\n');s.sendall(struct.pack('i',75 if mode!='physics' else 1));s.close();sys.exit(0)
os.execv(binary,[binary]+sys.argv[1:])
