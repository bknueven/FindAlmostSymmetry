#!/usr/bin/python3

from subprocess import STDOUT, call, TimeoutExpired
from time import sleep

runs =  [ ('games120.col',6),
          ('miles250.col',5),
          ('miles500.col',5),
          ('miles750.col',5),
          ('miles1000.col',5),
          ('miles1500.col',3),
          ('ran10_100_a.bliss',10),
          ('ran10_100_b.bliss',10),
          ('ran10_100_c.bliss',10),
          ('ran10_100_d.bliss',10),
          ]

## for each file above

wait_seconds = 60

threads = [1, 2, 4, 8, 16, 32,]

test_graph_dir = '../test_graphs/'

out_dir = './scaling_tests/'

for f,budget in runs:
    for thread_num in threads:
            run_list = ['unbuffer', 'mpirun', '-np', str(thread_num), '../findAlmost', 
                    '--budget=%d'%budget,
                    test_graph_dir+f, 'outfile','2>&1','>',
                    out_dir+f+'.budget%d'%budget+'.threads%d'%thread_num ]
            print('Running '+' '.join(run_list))
            #print('\tpass %d'%(i+1))
            call(' '.join(run_list), stderr=STDOUT, shell=True)

            print('Waiting for %d seconds before next job...\n'%wait_seconds) 
            sleep(wait_seconds)
