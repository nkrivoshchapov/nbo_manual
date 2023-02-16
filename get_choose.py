#!/usr/bin/python

import sys
import os
import glob
import re

NBO_HEADER = re.compile("Gaussian NBO Version 3\.1")
START_PATTERN = re.compile("\(occupancy\).*bond orbital.*/.*coefficients.*/.*hybrids")
END_PATTERN = re.compile("nho.*directionality.*and.*bond.*bending")

ALPHA_PATTERN = re.compile("\*         Alpha spin orbitals         \*")
BETA_PATTERN = re.compile("\*         Beta  spin orbitals         \*")

LP_PATTERN = re.compile("[0-9]+\\. \\(([+\\-]{0,1}[\\d]+(?:\\.[\\d]+)*)\\).*LP \( [0-9]+\).*")
BD_PATTERN = re.compile("[0-9]+\\. \\(([+\\-]{0,1}[\\d]+(?:\\.[\\d]+)*)\\).*BD \( [0-9]+\).*- .*")
C3LINE_PATTERN = re.compile("[0-9]+\\. \\(([+\\-]{0,1}[\\d]+(?:\\.[\\d]+)*)\\).*3C \( [0-9]+\).*- .*- .*")
NUM_PATTERN = re.compile("([0-9]+)")
BOND_PATTERN = re.compile(" ([0-9]+) .*-.* ([0-9]+) ")
C3_PATTERN = re.compile(" ([0-9]+) .*-.* ([0-9]+) .*-.* ([0-9]+) ")

ORDER_MATCH = {1: 'S', 2: 'D', 3: 'T', 4: 'Q'}
for i in range(5, 10):
    ORDER_MATCH[i] = str(i)


def get_main_start(lines):
    start = None
    for i, line in enumerate(lines):
        if re.search(NBO_HEADER, line) and re.search(NBO_HEADER, lines[i - 3]):
            start = i
    assert start is not None
    return start


def is_unrestricted(lines):
    alpha_found = False
    beta_found = False
    for i, line in enumerate(lines):
        if re.search(ALPHA_PATTERN, line):
            alpha_found = True
        elif re.search(BETA_PATTERN, line):
            beta_found = True
    assert (alpha_found and beta_found) or (not alpha_found and not beta_found), \
           f"I'm confused. alpha_found {alpha_found}, beta_found {beta_found}."
    return alpha_found


def get_start_end(lines):
    start = None
    end = None
    alpha_read = False
    beta_read = False

    def add_idx(target, value):
        if alpha_read:
            assert target is None
            return {'alpha': value}
        elif beta_read:
            assert 'alpha' in target
            target['beta'] = value
            return target
        else:
            assert target is None
            return value

    for i, line in enumerate(lines):
        if re.search(START_PATTERN, line.lower()):
            start = add_idx(start, i)
        elif re.search(END_PATTERN, line.lower()):
            end = add_idx(end, i)
        elif re.search(ALPHA_PATTERN, line):
            alpha_read = True
            beta_read = False # Redundant
        elif re.search(BETA_PATTERN, line):
            beta_read = True
            alpha_read = False
    assert start is not None
    assert end is not None
    return start, end


def get_nbo_lists(lines, start, end):
    nbo_lists = {
        'LP': [],
        'BD': [],
        'C3': []
    }
    for i in range(start, end):
        if re.search(LP_PATTERN, lines[i]):
            myline = lines[i].split("s(")[0].split(")")
            myline = myline[len(myline)-1]
            nbo_lists['LP'].append(int(NUM_PATTERN.search(myline).group(1)))
        elif re.search(BD_PATTERN, lines[i]):
            myline = lines[i].replace("\r", "").replace("\n", "").replace("-", " - ").split(")")
            myline = myline[len(myline)-1] + " "
            matches = BOND_PATTERN.search(myline)
            nbo_lists['BD'].append([int(matches.group(1)), int(matches.group(2))])
        elif re.search(C3LINE_PATTERN, lines[i]):
            myline = lines[i].replace("\r", "").replace("\n", "").replace("-", " - ").split(")")
            myline = myline[len(myline)-1] + " "
            matches = C3_PATTERN.search(myline)
            nbo_lists['C3'].append([int(matches.group(1)), int(matches.group(2)), int(matches.group(3))])
    assert len(nbo_lists['LP']) > 0 or len(nbo_lists['BD']) > 0 or len(nbo_lists['C3']) > 0, "NBO lists are empty"
    return nbo_lists


def gen_lp_section(nbo_lists):
    if len(nbo_lists['LP']) == 0:
        return []
    
    lp_counts = {}
    for i, atom_idx in enumerate(nbo_lists['LP']):
        if atom_idx not in lp_counts:
            lp_counts[atom_idx] = 1
        else:
            lp_counts[atom_idx] += 1
    
    reslines = []
    reslines.append("LONE")
    first = True
    for atom_idx, lp_count in lp_counts.items():
        if first:
            reslines[len(reslines) - 1] += f" {atom_idx} {lp_count}"
            first = False
        else:
            reslines.append(f"     {atom_idx} {lp_count}")
    reslines.append("END")
    return reslines


def gen_bd_section(nbo_lists):
    if len(nbo_lists['BD']) == 0:
        return []
    
    c2_counts = {}
    itemlist = []
    for atomA, atomB in nbo_lists['BD']:
        if (atomA, atomB) not in itemlist:
            c2_counts[(atomA, atomB)] = 1
            itemlist.append((atomA, atomB))
        else:
            c2_counts[(atomA, atomB)] += 1
    
    reslines = []
    reslines.append("BOND")
    first = True
    for bond, count in c2_counts.items():
        newline = "{} {} {}".format(ORDER_MATCH[count], *bond)
        if first:
            reslines[len(reslines) - 1] += " " + newline
            first = False
        else:
            reslines.append("     " + newline)
    reslines.append("END")
    return reslines


def gen_c3_section(nbo_lists):
    if len(nbo_lists['C3']) == 0:
        return []

    c3_counts = {}
    itemlist = []
    for atomA, atomB, atomC in nbo_lists['C3']:
        if (atomA, atomB, atomC) not in itemlist:
            c3_counts[(atomA, atomB, atomC)] = 1
            itemlist.append((atomA, atomB, atomC))
        else:
            c3_counts[(atomA, atomB, atomC)] += 1
    
    reslines = []
    reslines.append("3CBOND")
    first = True
    for bond, count in c3_counts.items():
        newline = "{} {} {} {}".format(ORDER_MATCH[count], *bond)

        if first:
            reslines[len(reslines) - 1] += " " + newline
            first = False
        else:
            reslines.append("       " + newline)
    reslines.append("END")
    return reslines


if __name__ == "__main__":
    myfilename = sys.argv[1]
    lines = open(myfilename, 'r').readlines()

    main_start = get_main_start(lines) # This is the line where NBO analysis starts
    lines = lines[main_start:]

    chooselines = ["$CHOOSE", "ALPHA"]
    if is_unrestricted(lines):
        start, end = get_start_end(lines)
        sections = {}
        for spin in ['alpha', 'beta']:
            nbo_lists = get_nbo_lists(lines, start[spin], end[spin])

            sections[spin] = []
            sections[spin] += gen_lp_section(nbo_lists)
            sections[spin] += gen_bd_section(nbo_lists)
            sections[spin] += gen_c3_section(nbo_lists)
        chooselines += sections['alpha']
        chooselines.append("END")
        chooselines.append("BETA")
        chooselines += sections['beta']
        chooselines.append("END")
    else:
        start, end = get_start_end(lines)
        nbo_lists = get_nbo_lists(lines, start, end)

        chooselines += gen_lp_section(nbo_lists)
        chooselines += gen_bd_section(nbo_lists)
        chooselines += gen_c3_section(nbo_lists)
        
    chooselines.append("$END")
    print("\n".join(chooselines))
