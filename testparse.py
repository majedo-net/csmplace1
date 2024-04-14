import sys, os
import pyverilog
from pyverilog.vparser.parser import parse

netfile = ['csmplace1/test.vg']

ast, directives = parse(netfile)

for item in ast.description.definitions[0].items:
    if type(item)==pyverilog.vparser.ast.InstanceList:
        print(item.module)