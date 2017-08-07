#!/usr/bin/python3

from subprocess import STDOUT, call, TimeoutExpired
from time import sleep

runs =  [ ('games120.col',7),
          ('miles250.col',5),
          ('miles500.col',5),
          ('miles750.col',5),
          ('miles1000.col',5),
          ('miles1500.col',3),
          ('le450_5b.col',18),
          ('le450_15b.col',15),
          ('le450_25b.col',9),
          ('ran10_100_a.bliss',10),
          ('ran10_100_b.bliss',10),
          ('ran10_100_c.bliss',10),
          ('ran10_100_d.bliss',10),
          ]

## for each file, budget above

wait_seconds = 60

test_graph_dir = '../../test_graphs/'

out_dir = './solver_output/'

for f,budget in runs:
            run_list = ['unbuffer', '../../findAlmost', 
                    '--budget=%d'%budget, '--trackEdges=True',
                    test_graph_dir+f, 'outfile','2>&1','>',
                    out_dir+f+'.budget%d'%budget+'.trackedEdges.out' ]
            print('Running '+' '.join(run_list))
            #print('\tpass %d'%(i+1))
            call(' '.join(run_list), stderr=STDOUT, shell=True)

            print('Waiting for %d seconds before next job...\n'%wait_seconds) 
            sleep(wait_seconds)
