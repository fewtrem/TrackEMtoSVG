'''
Created on 17 Nov 2016

@author: s1144899
'''
import ctypes, numpy as np
lib = ctypes.cdll.LoadLibrary('../cVersion2.so')
cdoCirc = lib.doCirc
caddSphere = lib.addSphere
def interpBez(bezPoints,t):
    tm = 1-t
    return tm*tm*tm*bezPoints[0]+3*tm*tm*t*bezPoints[1]+3*t*t*tm*bezPoints[2]+t*t*t*bezPoints[3]
def interpBezDer(bezPoints,t):
    tm = 1-t
    return 3*tm*tm*(bezPoints[1]-bezPoints[0])+6*tm*t*(bezPoints[2]-bezPoints[1])+3*t*t*(bezPoints[3]-bezPoints[2])
def doCrossCircC(i,t,outputImg,radius,bezCurvePoints):
    #print t
    centrePoint = interpBez(bezCurvePoints[i],t)
    perpDir = interpBezDer(bezCurvePoints[i],t)
    sumD = np.sum(np.power(perpDir,2))
    # check it has a length:
    if sumD !=0:
        rotAmount = 0.01
        incAmount = 0.5
        # Calculatable things:
        outputSize = outputImg.shape
        tDir = perpDir.copy()
        tMax = np.argmax(np.abs(perpDir))
        tDir[tMax] = -perpDir[tMax]
        initDir = np.cross(tDir,perpDir)
        initDirN = initDir/np.sqrt(np.sum(np.power(initDir,2)))
        perpDirN = perpDir/np.sqrt(np.sum(np.power(perpDir,2)))
        #print tMax
        #print tDir
        #print centrePoint
        #print initDirN
        #print perpDir
        #print np.dot(initDirN,perpDir)
        # Pointers:
        perpDirNP = (ctypes.c_double * len(perpDirN))(*perpDirN)
        initDirNP = (ctypes.c_double * len(initDirN))(*initDirN)
        centrePointP = (ctypes.c_double * len(centrePoint))(*centrePoint)
        outputSizeP = (ctypes.c_int * len(outputSize))(*outputSize)
        # Do the function:
        cdoCirc(perpDirNP,
                initDirNP,
                centrePointP,
                ctypes.c_void_p(outputImg.ctypes.data),
                outputSizeP,
                ctypes.c_double(rotAmount),
                ctypes.c_double(radius),
                ctypes.c_double(incAmount))
def Cfill(outputImg,bezCurvePoints,widths):
    outputImg = np.ascontiguousarray(outputImg)
    rotAmount = 0.01
    incAmount = 0.5
    print "  Doing the points:"
    for i in range(len(bezCurvePoints)):
        print " ",bezCurvePoints[i],
        for t in np.arange(0.0,1.0,0.001):
            radius = t*widths[i]+(1-t)*widths[i+1]
            doCrossCircC(i,t,outputImg,radius,bezCurvePoints)
        if i!=0:
            radius = widths[i]
            centrePointP = (ctypes.c_double * len(bezCurvePoints[i][0]))(*bezCurvePoints[i][0])
            outputSizeP = (ctypes.c_int * len(outputImg.shape))(*outputImg.shape)
            caddSphere(centrePointP,
            ctypes.c_void_p(outputImg.ctypes.data),
            outputSizeP,
            ctypes.c_double(rotAmount),
            ctypes.c_double(radius),
            ctypes.c_double(incAmount))
    print ""