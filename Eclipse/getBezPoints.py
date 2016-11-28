'''
Created on 17 Nov 2016

@author: s1144899
'''
import re, numpy as np
def getBezPoints(thisPipe):
    # points (x,y):
    pVals = thisPipe.attrib['d']
    # z value:
    zVals = [int(a) for a in thisPipe.attrib['layer_ids'].split(',')]
    # width:
    wVals = [float(a) for a in thisPipe.attrib['p_width'].split(',')]
    # get the x,y points:
    mSpl = re.split('[a-zA-Z ]',pVals)
    # get the pipe's tranform info for the bounding box around it:
    transformInfo = re.findall(r'[\(].*[\)]',thisPipe.attrib['transform'])[0].replace('(','').replace(')','').split(',')
    shifter = np.array([float(transformInfo[4]),float(transformInfo[5])])
    # get the x,y point values:
    # there are 2 values for the end points and 3 for the middle points:
    outList = []
    for tmSpl in mSpl:
        if bool(re.match(tmSpl,'.+,.+'))==False:
            nSpl = tmSpl.split(',')
            outList.append(np.array([float(nSpl[0])+shifter[0],float(nSpl[1])+shifter[1],0.0]))
    # add in the z values:
    for z in range(len(zVals)):
        midList = z*3
        for i in range(3):
            thisL = midList+i-1
            if thisL>=0 and thisL<len(outList):
                outList[thisL][2] = zVals[z]/3
    # now make a list of the BezCurve points:
    bezCurvePoints = []
    i = 0
    while i < len(outList):
        #sanity check:
        if i+3<len(outList):
            bezCurvePoints.append([outList[i],outList[i+1],outList[i+2],outList[i+3]])
        i+=3
    return bezCurvePoints,wVals