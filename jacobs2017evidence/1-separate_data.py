#!/usr/bin/env python

import re
import os
import errno

filepath = "../pnas.1705772114.sd01.txt"


with open(filepath) as f:
    datalines = f.readlines()

try:
    os.mkdir("data")
except OSError as exc:
    if exc.errno != errno.EEXIST:
        raise
    pass

begin_end = []
filename_re = re.compile("==> ([a-zA-Z0-9_./]+) <==")

begin_line = 1
end_line = None
filename = filename_re.findall(datalines[0])[0]
for line_index, line in enumerate(datalines):
    match = filename_re.match(line)
    if match:
        if line_index == 0:
            continue
        new_filename = match.groups()[0]
        end_line = line_index
        with open(filename, 'w') as f:
            f.writelines(datalines[begin_line:line_index])
        begin_line = line_index + 1
        filename = new_filename
