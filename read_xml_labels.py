#!/usr/bin/python

import os
import sys

import xml.etree.ElementTree as ET

def read_xml_labels(labelfn):
    tree = ET.parse(labelfn)
    labs = tree.getroot()
    a = labs[1]
    b = a.getchildren()

    labs = {"ind": [], "name": []}
    for tb in b:
        labs["ind"].append(int(tb.attrib["index"]))
        labs["name"].append(tb.text)

    return labs
