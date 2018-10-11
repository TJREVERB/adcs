# Inputs are:
#t = time
#  w0 = angular velocity
#  torque = total amount of torque on CubeSat
#  sc = Struct describing CubeSat
# inertia method needs to be written, use numpy to find derivative in numpy
def wfun(t,w0, torque, sc):
dwdt = wfun(t,w0, torque, sc)



function [dwdt] = wfun(t,w0,torque,sc)

w0 = w0(:)';
torque = torque(:)';
dwdt = sc.inertia/(torque'-cross(w0',(sc.inertia*w0'))');
