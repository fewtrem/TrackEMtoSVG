'''
Created on 17 Nov 2016

@author: s1144899
'''
theXMLPath = "../untitled.xml"
# import the XML:
import xml.etree.ElementTree as ET, numpy as np,nrrd
from getBezPoints import getBezPoints
from cFunc import Cfill
tree = ET.parse(theXMLPath)
root = tree.getroot()
pipes = root.findall('t2_layer_set')[0].findall('t2_pipe')
# Get the original dimensions:
origImDims = [0,0,0]
layerSet = root.findall('t2_layer_set')[0]
origImDims[0] = int(float(layerSet.attrib['layer_width']))
origImDims[1] = int(float(layerSet.attrib['layer_height']))
origImDims[2] = len(layerSet.findall('t2_layer'))
outputImg = np.zeros(np.array(origImDims)+np.array((20,20,20)),dtype=np.uint8)
for thisPipe in pipes:
    print "---pipe: ",thisPipe.attrib['title'],"---"
    bezCurvePoints,widths = getBezPoints(thisPipe)
    # Add the pipes:
    Cfill(outputImg,bezCurvePoints,widths)
nrrd.write("/media/s1144899/OS/Temp/trackEM2/c_ouput5.nrrd",outputImg*255)