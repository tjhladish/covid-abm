#!/usr/bin/python
from glob import glob
from sys import argv

merged_file = "ppb-v2-merged.csv"
fo = open(merged_file, 'w')
fo.write('date,' + ','.join(['ppb' + str(i) for i in range(1000)]) + '\n')

all_data = []
dates = []

for filename in glob(argv[1] + '/*.csv'):
    data = []
    header = True
    for line in open(filename):
        if header:
            header = False
            continue
        p = line.strip().split(',')
        data.append(p[1])
        if len(dates) < len(data):
            dates.append(p[0])
    all_data.append(data)

for i in range(len(dates)):
    fo.write(dates[i] + ',' + ','.join([all_data[j][i] for j in range(len(all_data))]) + '\n')

fo.close()
