# Generated with SMOP  0.41
from libsmop import *
#

    ##mate01.m
## steel, elastic perfectly plastic E=29000Ksi, Sigma_y = 50ksi
## arguments(2)
##               isw: task
##               e:   strain


@function
def mate01(isw=None,e=None,*args,**kwargs):
    varargin = mate01.varargin
    nargin = mate01.nargin

    ee=abs(e)
    ss=sign(e)
    epsy=0.0017
    sigy=50
#----------------------------------------------------------------------
#      Tangent
#----------------------------------------------------------------------
    if isw == 1:
        if (ee <= epsy):
            r=sigy / epsy
        else:
            r=0
        #----------------------------------------------------------------------
#    stress
#----------------------------------------------------------------------
    else:
        if isw == 2:
            if (ee < epsy):
                r=dot(sigy,e) / epsy
            else:
                if (ss == 1):
                    r=copy(sigy)
                else:
                    r=- sigy
