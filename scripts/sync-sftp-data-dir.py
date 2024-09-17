#!/usr/bin/python3

# Exmemple usage
# ./scripts/sync-sftp-data-dir.py ${HOME}/data sftp://carrard@inti.ocre.cea.fr:/ccc/home/cont001/xstampdev/xstampdev/data
# ./scripts/hash-data-dir.sh /usr/local/xstampv2/data

import sys
import urllib.parse
import paramiko
import getpass
import tempfile
import os
import random
import time

def read_md5_file(filepath):
  db={}
  for l in open(filepath).readlines():
    (name,H,date) = l.split()
    db[name] = (H,date)
  return db

def write_md5_file(db,filepath):
  f=open(filepath,"w")
  for k in db:
    f.write( "%s %s %s\n" % (k,db[k][0],db[k][1]) )
  f.close()

# local dir and sftp URL
localdir = sys.argv[1]
remotedir = sys.argv[2]

print("Connect to remote directory %s ..."%remotedir)
u = urllib.parse.urlparse( remotedir )
p=u.port
if not p: p=22

t = paramiko.Transport( (u.hostname,p) )
t.connect(username=u.username,password=getpass.getpass())
sftp = paramiko.SFTPClient.from_transport(t)

print("remote chdir to %s"%u.path)
sftp.chdir(u.path)
os.chdir(localdir)

# read local MD5 data base
localmd5 = localdir+"/data_dir.md5"
print("read local MD5 DB",localmd5)
local_db = read_md5_file(localmd5)

# wait for distant DB to be ready and acquire MD5 data base
remotemd5 = tempfile.mktemp()
while not "data_dir.md5" in sftp.listdir():
  time.sleep(2)
print("lock remote MD5 data base")
sftp.get("data_dir.md5",remotemd5)
sftp.remove("data_dir.md5")
remote_db = read_md5_file(remotemd5)

def file_exists(obj,pfx,name):
  e=True
  if pfx: name=pfx+'/'+name
  try:
    obj.stat(name)
  except:
    e=False
  return e

def remote_file_exists(name):
  return file_exists(sftp,None,name)

def local_file_exists(name):
  return file_exists(os,localdir,name)

def copy_to_remote(name):
  #print("copy to remote ",name)
  remote_db[name] = local_db[name]
  if remote_file_exists(name):
    #print("remote rm",name)
    sftp.remove(name)
  if remote_file_exists(name+".xz"):
    #print("remote rm",name+".xz")
    sftp.remove(name+".xz")
  if not local_file_exists(name):
    #print("local file is compressed")
    name=name+".xz"
  #print("put",name)
  sftp.put(name,name)

def copy_from_remote(name):
  #print("copy from remote",name)
  local_db[name] = remote_db[name]
  if local_file_exists(name):
    #print("local rm",name)
    os.remove(name)
  if local_file_exists(name+".xz"):
    #print("local rm",name+".xz")
    os.remove(name+".xz")
  if remote_file_exists(name+".xz"):
    #print("remote file is compressed")
    name=name+".xz"
  #print("get",name)
  sftp.get(name,name)

for k in local_db:
  if k in remote_db.keys():
    if local_db[k][0] != remote_db[k][0]:
      if local_db[k][1] > remote_db[k][1]:
        print("%-40s UPDATE LOCAL -> REMOTE" % k)
        copy_to_remote(k)
      else:
        print("%-40s UPDATE REMOTE -> LOCAL" % k)
        copy_from_remote(k)
  else:
    print("%-40s COPY LOCAL -> REMOTE" % k)
    copy_to_remote(k)

for k in remote_db:
  if not k in local_db.keys():
    print("%-40s COPY REMOTE -> LOCAL" % k)
    copy_from_remote(k)

# update MD5 data bases (remote and local)
print("remote write MD5 DB")
os.remove(remotemd5)
write_md5_file(remote_db,remotemd5)
sftp.put(remotemd5,"data_dir.md5")

print("local write MD5 DB")
os.remove(localmd5)
write_md5_file(local_db,localmd5)

print("done")
