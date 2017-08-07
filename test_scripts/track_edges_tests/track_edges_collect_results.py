#!/usr/bin/python3

import pandas as pd
import csv

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

default_top_num = 10

for f,b in runs:

    file_name = f+'.trackedEdges.budget'+str(b)+'.txt'
    counts = pd.read_csv(file_name, header=None)[0]
    percentages = counts.value_counts()/counts.count()

    ## sort so slicing works right
    percentages.sort_index(inplace=True)

    top_num = min(default_top_num, percentages.count())

    with open(file_name+'.summary.csv', 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter='\t')
        for i in range(top_num):
            csvwriter.writerow([str(i+1), str(percentages[i])])
        ## combine all others
        csvwriter.writerow([str(top_num+1)+'+', str(percentages[top_num:].sum())])

    del counts, percentages
