#!/usr/bin/python3
import re
import pandas as pd

from table2_collect_results import results_sol

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

optimal_obj = {}
ground_orbits = {}
for f,max_budget in zip(files, max_budgets):
    ground_orbits[f] = min(results_sol[f, 0])
    for budget in range(1, max_budget+1):
        optimal_obj[f, budget] = min(results_sol[f, budget])


test_graph_dir = '../test_graphs/'

out_dir = './heuristicoutput/'

results_time = {}
results_sol = {}
for f,max_budget in zip(files, max_budgets):
    for budget in range(1, max_budget):
        for boolian in ['True', 'False']:
            with open(out_dir+f+'.budget%d'%budget+'.localBranching%s'%boolian, 'r') as out_file:
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
                if ground_orbits[f] == optimal_obj[f, budget]:
                    results_sol[f, budget, eval(boolian)] = '-'
                else:
                    results_sol[f,budget,eval(boolian)] =  (ground_orbits[f] -  int(obj_val)) / \
                                                    ( ground_orbits[f] - optimal_obj[f, budget] )

                ## get run time
                match = re.search('Wall clock total time\s*=\s*\d*\.\d seconds', out_file_str)
                if match == None:
                    time = 'NaN'
                else:
                    time_line = match.group(0)
                    time = re.search('\d\d*\.\d', time_line).group(0)
                results_time[f,budget,eval(boolian)] = (float(time))

#invert dictionaries for writing

sol_results = {}
time_results = {}
for boolian in [True, False]:
    sol_results[boolian] = {}
    time_results[boolian] = {}
    for f,max_budget in zip(files, max_budgets):
        sol_results[boolian][f] = {}
        time_results[boolian][f] = {}
        for budget in range(1, max_budget):
            sol_results[boolian][f][budget] = results_sol[f,budget,boolian]
            time_results[boolian][f][budget] = results_time[f,budget,boolian]


pd.DataFrame(sol_results[False]).T.to_excel('heuristic_results_global_branching.xlsx')
pd.DataFrame(sol_results[True]).T.to_excel('heuristic_results_local_branching.xlsx')
pd.DataFrame(time_results[False]).T.to_excel('heuristic_results_global_branching_time.xlsx')
pd.DataFrame(time_results[True]).T.to_excel('heuristic_results_local_branching_time.xlsx')


