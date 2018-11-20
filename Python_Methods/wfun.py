# Inputs are:
#t = time
#  w0 = angular velocity
#  torque = total amount of torque on CubeSat
#  sc = Struct describing CubeSat
# inertia method needs to be written, use numpy to find derivative in numpy
def wfun(t,w0, torque, sc):
  dwdt = wfun(t,w0, torque, sc)
  w0=np.reshape(w0, (1,-1))
  torque=np.reshape(torque, (1,-1))
  dwdt = sc.inertia/(torque.getH()-np.cross(w0.getH(),(sc.inertia*w0.getH())).getH())
  return(dwdt)
