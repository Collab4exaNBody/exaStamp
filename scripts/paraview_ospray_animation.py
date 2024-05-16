### Get local MPI rank and set frames to render ###
###########################################


### High quality path tracing render view ###
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [800, 450]
renderView1.AxesGrid = 'Grid Axes 3D Actor'
renderView1.CenterOfRotation = [300.0, 150.0, 75.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [306.422975733372, 249.71760525318993, 817.8908525794668]
renderView1.CameraFocalPoint = [300.00000000000006, 150.00000000000006, 74.99999999999997]
renderView1.CameraViewUp = [-0.00831927521262538, 0.9910863877537779, -0.13296075236363908]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 343.693177121688
renderView1.LegendGrid = 'Legend Grid Actor'
renderView1.EnableRayTracing = 1
renderView1.BackEnd = 'OSPRay pathtracer'
renderView1.SamplesPerPixel = 64
renderView1.OSPRayMaterialLibrary = materialLibrary1
#################################################

### Save animation ###
import os
import sys

def CustomSaveAnimation(FRAME_START,FRAME_END,WORKER_ID,NB_WORKERS):
  print("ParaView animation worker %d / %d : renders frames [ %d ; %d [" % (WORKER_ID,NB_WORKERS,FRAME_START,FRAME_END) )
  for I in range(FRAME_START+WORKER_ID,FRAME_END,NB_WORKERS):
      filename = "images/img.%04d.png" % I
      if os.path.isfile(filename):
          print("Worker %d skips file %s" % (WORKER_ID,filename) )
      else:
          print("Worker #%d renders image %d" % (WORKER_ID,I) )
          animationScene1.AnimationTime = float(I)
          SaveScreenshot( filename , GetLayout() , ImageResolution=[1920,1080] )
      sys.stdout.flush()

CustomSaveAnimation( 0 , 100 , int(os.getenv("SLURM_PROCID")) , int(os.getenv("SLURM_NTASKS")) )
#####################

