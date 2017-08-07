#!/usr/bin/python3

from subprocess import STDOUT, call, TimeoutExpired
from time import sleep
import re
import gc

files = [ 'games120.col',
          'miles250.col',
          'miles500.col',
          'miles750.col',
          'miles1000.col',
          'miles1500.col',
          'le450_5b.col',  
          'le450_15b.col',
          'le450_25b.col',
          'ran10_100_a.bliss',
          'ran10_100_b.bliss',
          'ran10_100_c.bliss',
          'ran10_100_d.bliss',
          ]

## for each file above
max_budgets = [ 9, 7, 8, 8, 7, 5, 19, 17, 10,
                12, 12, 13,
               12 ]

repeat_times = 1
time_out_min = 60

wait_seconds = 10

test_graph_dir = '../test_graphs/'

out_dir = './auto_output/'
out_dir_tee = './auto_output_tee/'

for f,max_budget in zip(files, max_budgets):
    for budget in range(max_budget+1):
        for i in range(repeat_times):
            run_list = ['unbuffer', 'timeout', '%d'%((time_out_min+1)*60), 
                    'mpirun', '-np', '32', '../findAlmost', 
                    '--budget=%d'%budget, '--earlyOutputMinutes=%d'%(time_out_min-1),
                    '--maxWallMinutes=%d'%time_out_min, '--output='+out_dir+f+'.budget%d'%budget+'.run%d'%i+'.sol.txt',
                    test_graph_dir+f, 
                    '2>&1','>',
                    out_dir_tee+f+'.budget%d'%budget+'.run%d'%i]
            print('Running '+' '.join(run_list))
            #print('\tpass %d'%(i+1))
            call(' '.join(run_list), stderr=STDOUT, shell=True)

            print('Waiting for %d seconds before next job...\n'%wait_seconds) 
            sleep(wait_seconds)
