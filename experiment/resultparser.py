#! pyhthon3 resultparser.py
from enum import Enum
import os
import csv
import shutil
import numpy as np
import math 
import matplotlib.pyplot as plt

d=os.path.dirname(os.getcwd())
opt_file = open(d+"/benchmark/maxcut.opts", "r")
ls = opt_file.readlines()

opt_dict = {}

mclasses = ["pm1", "w01", "t2g", "t3g"]

instances = []


for line in ls:
   line = line.split(" ")
   entry = {}
   opt_dict[line[0]] = int(line[1])
   instances.append(line[0])



log_path = d + "/result"
logs = os.listdir(log_path)

entries = []


def extract_result(entry_, file_path):
   ls = open(file_path).readlines()
   entry={}
   for (i, l) in enumerate(ls):
      l = l.strip()
      l = l.split(": ")
      val = l[1]
      if i == 1:
         val = float(val)
      if i == 2:
         val = float(val)
      elif i == 4:
         val = int(val)
      entry[l[0]] = val
   return entry


for log in logs:
   entry = extract_result(entry, log_path + "/"+log)
   entry["mclass"] = entry["instance"][0:3]
   entries.append(entry)


def Stat(instance, algorithm, level, mclass):
   return {"instance": instance, "setting": algorithm, "level": level, "mclass": mclass,  "total": 0, "time": 0.0, "time_lst": [], "gap": 0.0,  "gap_lst": []} 

display_keys = [ "gap", "time"]


SheraliAdams1 = ("SheraliAdams", 1)
Lasserre1 = ("Lasserre", 1)
Submodular1 = ("Submodular", 1)
Submodular2 = ("Submodular", 2)
Submodular3 = ("Submodular", 3)
Submodular4 = ("Submodular", 4)
settings = (SheraliAdams1, Lasserre1, Submodular1, Submodular2, Submodular3, Submodular4)

for entry in entries:
   opt = opt_dict[entry["instance"]]
   bd = entry["obj"]
   entry["gap"] = abs(bd - opt) / bd 
   try:
      assert bd > opt
   except AssertionError as e:
      print(e, " ", bd, " ", opt)


def add(stat, entry):
   stat["total"] += 1
   stat["time_lst"].append(entry["time"])
   stat["gap_lst"].append(entry["gap"])


def SGM(lst, total, bias):
   return np.exp(np.sum([np.log(ele + bias) for ele in lst]) / total) - bias


def avgStat(stat):
   stat["time"] = SGM(stat["time_lst"], stat["total"], 1)
   stat["gap"] = SGM(stat["gap_lst"] , stat["total"], 0.1)


def printStat(stat):
   #print(setting, pclass, stat)
   s = [ str(round(stat[display_key], 1) if display_key is "gap" else round(stat[display_key], 3)) + " & " for display_key in display_keys ]
   print("".join(s), end="")



stats = {}
allstat = {}
for (algorithm, level) in settings:
   print(algorithm, " ", level, " & ", end="")
   allstat[(algorithm, level)] = Stat("all", algorithm, level, "all")
   for mclass in mclasses:
      stats[(algorithm, level, mclass)] = Stat("all", algorithm, level, mclass) 
      for entry in entries: 
         if entry["mclass"] == mclass and entry["algo"] == algorithm and entry["level"] == level:
            add(stats[(algorithm, level, mclass)], entry)   
            add(allstat[(algorithm, level)], entry)             
      avgStat(stats[(algorithm, level, mclass)])
      printStat(stats[(algorithm, level, mclass)])
   avgStat(allstat[(algorithm, level)] )     
   printStat(allstat[(algorithm, level)])
   print(" \\\\ ")


def toName(setting):
   return str(setting[0]) + " " + str(setting[1])

pairs = ((SheraliAdams1, Submodular1), (Lasserre1, Submodular1), (Submodular2, Submodular1),  (Submodular3, Submodular1),  (Submodular4, Submodular1), (Lasserre1, Submodular4),)


data= {}
k = 0
for i in range(0,2):
    for j in range(0,3):
      pair = pairs[k]
      data[(k, pair[0])] = []
      data[(k, pair[1])] = []
      k += 1

print(data.keys(), pairs)

for (k, (setting0, setting1)) in enumerate(pairs):
   print(k, setting0, setting1)
   for instance in instances:
         for entry in entries:
               if entry["instance"] == instance and entry["algo"] == setting0[0] and entry["level"] == setting0[1]:
                  data[(k, setting0)].append(entry["gap"])
               if entry["instance"] == instance and entry["algo"] == setting1[0] and entry["level"] == setting1[1]:
                  data[(k, setting1)].append(entry["gap"])


fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(10, 6))  # define the figure and subplots



num_instance = len(instances)
for i in range(0,2):
    for j in range(0,3):
      k = 3 * i + j
      pair = pairs[k]
      setting0 = pair[0]
      setting1 = pair[1]
      wins = [0,0]
      maxd = 100 if k == 0 else 0.5
      for l in range(num_instance):
         if data[(k, setting0)][l] < data[(k, setting1)][l]:
            wins[0] += 1
         else:
            wins[1] += 1
      axes[i,j].scatter(data[(k, setting0)], data[(k, setting1)], color = 'blue', marker = '+')
      axes[i,j].plot([0, maxd], [0, maxd], color = 'green')
      axes[i,j].set_xlabel(toName(setting0) +" wins " + str(wins[0]))
      axes[i,j].set_ylabel(toName(setting1) + " wins " +  str(wins[1]))

fig.tight_layout()
plt.savefig('scatter_bpo.pdf') 