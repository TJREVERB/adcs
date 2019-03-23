from .q_inv import q_inv
from .q_mult import q_mult
from numpy import linalg

def getqerr(q,qref)
  qrefinv = q_inv(qref)
  qerr = q_mult(q,qref)
  qerr = q_mult(qrefinv,q)
  qerr = qerr/linalg.norm(qerr)
  return qerr
