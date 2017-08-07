#!/usr/bin/python
# runs ./findAlmost with randomBranching=true for a specified budget 50 times
# and saves the data in a file named based on the input file

import sys
import subprocess
from multiprocessing import Pool
import re

input_dir = '../test_graphs/'
output_dir = './random_results/'
replications = 50

if len(sys.argv) != 3:
    print "usage: 'testRandomBranching.py <inputgraph> <budget>'"
    sys.exit(2)

def runFindAlmost(counter):
    inputGraph = sys.argv[1]
    budget = sys.argv[2]
    findAlmost_output = subprocess.check_output( 
            ["../findAlmost", "--budget=" + budget, "--randomBranching=true", input_dir + inputGraph, "outfile"] )
    print findAlmost_output
    nodes_line = re.search('Created\s*\d\d*\s*100.0%', findAlmost_output).group(0)
    nodes = int(re.search('\d\d*', nodes_line).group(0))
    print nodes
    return nodes


pool = Pool()
nodes = pool.map(runFindAlmost, range(replications))
pool.close()
pool.join()

print "nodes = ", nodes

outFile = open(output_dir + sys.argv[1] + '.budget' + sys.argv[2] + '.random_tests', 'w')

for count in nodes:
    print >>outFile, count

outFile.close()
