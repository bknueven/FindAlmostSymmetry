#!/usr/bin/python

import subprocess
import re

test_instances = [ ('games120.col', 2, 5), \
                   ('miles500.col', 2, 2), \
                   ('miles750.col', 2, 2), \
                   ('miles1000.col', 2, 2), \
                   ('miles1500.col', 2, 2), \
                   ('le450_15b.col', 6, 7), \
                   ('le450_25b.col', 4, 5), \
                   ('ran10_100_a.bliss', 9, 9), \
                   ('ran10_100_b.bliss', 9, 9), \
                   ('ran10_100_c.bliss', 9, 9), \
                   ('ran10_100_d.bliss', 9, 9), \
                 ]


for f, low, high in test_instances:
    for B in range(low,high+1):
        print "running ./testRandomBranching.py " + str(f) + " " + str(B)
        subprocess.call(["./testRandomBranching.py", str(f), str(B)])


for f, low, high in test_instances:
    for B in range(low,high+1):
        call_list = ["../findAlmost", "--budget=" + str(B), "../test_graphs/" + str(f), "outfile"]
        print "running " + ' '.join(call_list)
        findAlmost_output = subprocess.check_output( call_list )
        print findAlmost_output
        nodes_line = re.search('Created\s*\d\d*\s*100.0%', findAlmost_output).group(0)
        nodes = int(re.search('\d\d*', nodes_line).group(0))
        print nodes
        with open( './node_counts/'+str(f)+'.budget'+str(B)+'.nodecount', 'w' ) as outFile:
            print >>outFile, nodes


for f, low, high in test_instances:
    for B in range(low,high+1):
        call_list = ["../findAlmost", "--budget=" + str(B), "--localBranching=true", "../test_graphs/" + str(f), "outfile"]
        print "running " + ' '.join(call_list)
        findAlmost_output = subprocess.check_output( call_list )
        print findAlmost_output
        nodes_line = re.search('Created\s*\d\d*\s*100.0%', findAlmost_output).group(0)
        nodes = int(re.search('\d\d*', nodes_line).group(0))
        print nodes
        with open( './node_counts/'+str(f)+'.budget'+str(B)+'.localBranching.nodecount', 'w' ) as outFile:
            print >>outFile, nodes
