#!/usr/bin/python3
import re


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
max_budgets = [ 9, 7, 8, 8, 7, 5, 19, 17, 10, 12, 12, 13, 12 ]

repeat_times = 5
time_out_min = 60

wait_seconds = 60

test_graph_dir = '../test_graphs/'

out_dir = './table2output/'

results_time = {}
results_sol = {}
for f,max_budget in zip(files, max_budgets):
    for budget in range(max_budget+1):
        results_time[f, budget] = []
        results_sol[f, budget] = []
        for i in range(repeat_times):
            with open(out_dir+f+'.budget%d'%budget+'.run%d'%i, 'r') as out_file:
                out_file_str = out_file.read()                
                ## get objective value
                match = re.search('Best Solution:  Value = \d*', out_file_str)
                if match == None:
                    ## do something to get best obj vale
                    ## find all incumbent updates
                    match = re.findall('inc=\d\d*', out_file_str)
                    ## get the very last one
                    last_reported_inc = match.pop()
                    ## get the objective value
                    obj_val = re.search('\d\d*', last_reported_inc).group(0)
                else:
                    obj_val_line = match.group(0)
                    obj_val = re.search('\d\d*', obj_val_line).group(0)
                results_sol[f,budget].append(int(obj_val))

                ## get run time
                match = re.search('Average search time \(Wall clock\)\s*\d*\.\d seconds.', out_file_str)
                if match == None:
                    time = 'NaN'
                else:
                    time_line = match.group(0)
                    time = re.search('\d\d*\.\d', time_line).group(0)
                results_time[f,budget].append(float(time))


