import os

def run_tests(runcmd,optlist,ctx,tag=""):
  if len(optlist)==0:
    for s in ctx["tagrepl"]: tag = tag.replace(s[0],s[1])
    i = ctx["testidx"] + 1
    N = ctx["nbtests"]
    print("echo \"%d / %d : %s\"" % (i,N,tag))
    print("echo \"%d / %d : %s\"" % (i,N,runcmd))
    logfile="log"+tag+".txt"
    resultdir="result"+tag
    print("mkdir -p %s " % resultdir )
    print("%s | tee %s/%s" % (runcmd,resultdir,logfile) )
    print("echo \"copy results to %s ...\"" % resultdir )
    print("mv -f %s %s/" % (ctx["tracefiles"],resultdir) )
    ctx["testidx"] += 1
  else:
    kv = optlist[0]
    rem = optlist[1:]
    opt = kv[0]
    #if tag != "": tag = tag+","
    for val in kv[1]:
      ncmd = "%s %s %s" % ( runcmd , opt , val )
      run_tests(ncmd,rem,ctx,tag+str(opt)+"="+str(val))

def number_of_combinations(l):
  N=1
  for kv in l: N=N*len(kv[1])
  return N


def run_variation_sets(ctx,tag=""):
  opsets = [ [ x for x in ol.items() ] for ol in ctx["opsets"] ]
  ctx["testidx"] = 0
  ctx["nbtests"] = 0
  for optlist in opsets: ctx["nbtests"] += number_of_combinations(optlist)
  print("echo \"total number of tests = %d\"" % ctx["nbtests"] )
  for optlist in opsets: run_tests(ctx["runcmd"],optlist,ctx,tag)


def run_variations_recurse(ctx, cmdformat, runvals, runopts, tag=""):
  if len(runopts)==0:
      if ctx.keys().__contains__("runoptfunc"):
          runvals=ctx["runoptfunc"]( * tuple(runvals) )
      ctx["runcmd"] = cmdformat % tuple(runvals)
      print("\n\n# base command = ",ctx["runcmd"])
      run_variation_sets(ctx,tag)
  else:
      optvals=runopts[0][1]
      optname=runopts[0][0]
      for optvalue in optvals:
          run_variations_recurse(ctx,cmdformat,runvals+[optvalue],runopts[1:],tag+"_"+optname+"="+str(optvalue))

def run_variations(ctx):
    run_variations_recurse(ctx,ctx["runcmd"],[],ctx["runopts"])

# here after is an exemple usage of this module to run a series of tests

'''
#!/usr/bin/python3

from run_tests_funcs import *

# ------------------ test case paremeters ------------------------

case = "microjet.msp"
endit = 600040
opts = "--profiling-summary true"
part = "milan-bxi"
nmpi = 1
ncores = 128
ht = 2


# ------------------ run configuration ------------------------

ctx = {}
ctx["opsets"] = [ \
        { "--set-update_particle_neighbors":["ref"] , "--set-chunk_neighbors-chunk_size":["4","8"] , "--set-rebuild_amr-sub_grid_density":["6"] , "--set-rebuild_amr-enforced_ordering":["1","2"] } , \
  { "--set-update_particle_neighbors":["ref"] , "--set-chunk_neighbors-chunk_size":["4","8"] , "--set-rebuild_amr-sub_grid_density":["12"] , "--set-rebuild_amr-enforced_ordering":["4"] } , \
  { "--set-update_particle_neighbors":["lw"] , "--set-chunk_neighbors_lightweight-chunk_size":["4","8"] , "--set-rebuild_amr-sub_grid_density":["16","24","32"] , "--set-rebuild_amr-enforced_ordering":["2","4"] } \
]
ctx["tagrepl"] = [("--set-","_"),("update_particle_neighbors","cn"),("chunk_neighbors-chunk_size","chk"),("chunk_neighbors_lightweight-chunk_size","chk"), ("rebuild_amr-sub_grid_density","dens"),("rebuild_amr-enforced_ordering","sr") ]
ctx["runcmd"] = "ccc_mprun -p%s -n%d -c%d /usr/bin/env OMP_NUM_THREADS=%d ./xstampv2 %s --set-global-simulation_end_iteration %d %s" % (part,nmpi,ncores,ncores*ht,case,endit,opts)
ctx["tracefiles"] = "trace.*"

# ----------------------------------------------------------------------------

run_variation_sets( ctx )
'''

