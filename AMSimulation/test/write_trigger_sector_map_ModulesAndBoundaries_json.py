#!/usr/bin/env python

import json

mymap = {}
with open("../data/trigger_sector_map.csv", "r") as f:
    for line in f:
        if not line[0].isdigit():
            continue
        values = line.split(",")

        # Convert to int
        values = [int(x) for x in values]

        key = (values[0]-1)*8 + (values[1]-1)
        values = sorted(values[2:])
        mymap[key] = values

assert(len(mymap) == 6*8)
json.dump(mymap, open("../data/trigger_sector_map.json", "w"), sort_keys=True)


mymap = {}
with open("../data/trigger_sector_boundaries.csv", "r") as f:
    for line in f:
        if not line[0].isdigit():
            continue
        values = line.split(",")
        assert(len(values) == 6)

        # Convert to int or float
        values = [float(x) if "." in x else int(x) for x in values]

        key = values[0]*100 + values[1]
        values = values[2:]
        mymap[key] = values

json.dump(mymap, open("../data/trigger_sector_boundaries.json", "w"), sort_keys=True)
