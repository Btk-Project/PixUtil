#!/usr/bin/env python3
import os

s  = "/* Generated by cat_all.py */"

s += "#ifndef _BTK_PROJECT_PIXUTIL_ALL_HPP_\n"
s += "#define _BTK_PROJECT_PIXUTIL_ALL_HPP_\n"
s += "\n"

with open("pixutil.hpp","r") as f:
    s += f.read()

for file in os.listdir():
    if not file.startswith("pixutil_"):
        continue
    with open(file, "r") as f:
        s += f.read()
        s += "\n"
#Remove #include "pixutil.hpp"
s = s.replace("#include \"pixutil.hpp\"", "")

s += "#endif //_BTK_PROJECT_PIXUTIL_ALL_HPP_\n"

with open("pixutils_all.hpp","w") as f:
    f.write(s)