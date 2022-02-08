#!/usr/bin/python

import sys,os,glob,re

myfilename = sys.argv[1]
lines = open(myfilename,"r").readlines()
start_pattern = re.compile("\(occupancy\).*bond orbital.*/.*coefficients.*/.*hybrids")
end_pattern = re.compile("nho.*directionality.*and.*bond.*bending")
for i in range(len(lines)):
    if re.search(start_pattern,lines[i].lower()):
        start = i
    elif re.search(end_pattern,lines[i].lower()):
        end = i

lp_pattern = re.compile("[0-9]+\\. \\(([+\\-]{0,1}[\\d]+(?:\\.[\\d]+)*)\\).*LP \( [0-9]+\).*")
bd_pattern = re.compile("[0-9]+\\. \\(([+\\-]{0,1}[\\d]+(?:\\.[\\d]+)*)\\).*BD \( [0-9]+\).*- .*")
c3line_pattern = re.compile("[0-9]+\\. \\(([+\\-]{0,1}[\\d]+(?:\\.[\\d]+)*)\\).*3C \( [0-9]+\).*- .*- .*")
num_pattern = re.compile("([0-9]+)")
bond_pattern = re.compile(" ([0-9]+) .*-.* ([0-9]+) ")
c3_pattern = re.compile(" ([0-9]+) .*-.* ([0-9]+) .*-.* ([0-9]+) ")
lp_list = []
bd_list = []
c3_list = []
for i in range(start,end):
    if re.search(lp_pattern,lines[i]):
        myline = lines[i].split("s(")[0].split(")")
        myline = myline[len(myline)-1]
        lp_list.append(int(num_pattern.search(myline).group(1)))
    elif re.search(bd_pattern,lines[i]):
        myline = lines[i].replace("\r","").replace("\n","").replace("-"," - ").split(")")
        myline = myline[len(myline)-1] + " "
        matches = bond_pattern.search(myline)
        bd_list.append([int(matches.group(1)),int(matches.group(2))])
    elif re.search(c3line_pattern,lines[i]):
        myline = lines[i].replace("\r","").replace("\n","").replace("-"," - ").split(")")
        myline = myline[len(myline)-1] + " "
        matches = c3_pattern.search(myline)
        c3_list.append([int(matches.group(1)),int(matches.group(2)),int(matches.group(3))])
chooselines = ["$CHOOSE","ALPHA"]
# PROCESS LONE PAIRS
if len(lp_list)>0:
    lp_groups = []
    for i in range(len(lp_list)):
        ngroup = -1
        for j in range(len(lp_groups)):
            if lp_groups[j][0] == lp_list[i]:
                ngroup = j
        if ngroup == -1:
            lp_groups.append([lp_list[i], 1])
        else:
            lp_groups[ngroup][1] += 1
    first = True
    chooselines.append("LONE")
    for item in lp_groups:
        if first:
            chooselines[len(chooselines)-1] += " %d %d" % (item[0], item[1])
            first = False
        else:
            chooselines.append("      %d %d" % (item[0], item[1]))
    chooselines.append("END")

# PROCESS 2C BONDS
c2_groups = []
if len(bd_list)>0:
    for i in range(len(bd_list)):
        ngroup = -1
        for j in range(len(c2_groups)):
            if c2_groups[j][0] == bd_list[i][0] and c2_groups[j][1] == bd_list[i][1]:
                ngroup = j
        if ngroup == -1:
            c2_groups.append([bd_list[i][0], bd_list[i][1], 1])
        else:
            c2_groups[ngroup][2] += 1
    first = True
    chooselines.append("BOND")
    for item in c2_groups:
        if item[2] == 1:
            newline = "S %d %d" % (item[0], item[1])
        elif item[2] == 2:
            newline = "D %d %d" % (item[0], item[1])
        elif item[2] == 3:
            newline = "T %d %d" % (item[0], item[1])
        
        if first:
            chooselines[len(chooselines)-1] += " " + newline
            first = False
        else:
            chooselines.append("     " + newline)
    chooselines.append("END")

# PROCESS 3C BONDS
if len(c3_list)>0:
    c3_groups = []
    for i in range(len(c3_list)):
        ngroup = -1
        for j in range(len(c3_groups)):
            if c3_groups[j][0] == c3_list[i][0] and c3_groups[j][1] == c3_list[i][1] and c3_groups[j][2] == c3_list[i][2]:
                ngroup = j
        if ngroup == -1:
            c3_groups.append([c3_list[i][0], c3_list[i][1], c3_list[i][2], 1])
        else:
            c3_groups[ngroup][3] += 1
    chooselines.append("3CBOND")
    first = True
    for item in c3_groups:
        if item[3] == 1:
            newline = "S %d %d %d" % (item[0], item[1], item[2])
        elif item[3] == 2:
            newline = "D %d %d %d" % (item[0], item[1], item[2])
        elif item[3] == 3:
            newline = "T %d %d %d" % (item[0], item[1], item[2])
            
        if first:
            chooselines[len(chooselines)-1] += " " + newline
            first = False
        else:
            chooselines.append("       " + newline)
    chooselines.append("END")
chooselines.append("$END")
print("\n".join(chooselines))
