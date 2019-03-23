def getthetaerr(q):
  q = q.getH()
  thetaerr = 2*(q[0:3]/q[4])
  return thetaerr
