# This script generates a basic header file from the sirius_api.cpp file
# A lot of assumptions are made: all functions are void, we only care about the extern "C",
# etc. In the end, the idea is to provide an input for Clang.jl wrapping

import re
import sys
 
#Only read the file from the extern "C" section beginning
is_extern = False
content = ""
with open("sirius_api.cpp", "r") as myfile:
    for line in myfile:
        if not is_extern:
            if line.startswith("extern"):
                is_extern = True

        if is_extern:
            content += line

#We want C-style complex types, not C++
content = content.replace("std::complex<double>", "double complex")

#A check that the handler naming convention is respected:
#handler for ctx, gs_handler, kp_handler and H0_handerl for the rest
handlers = {"sim_ctx": "handler__", "ks": "ks_handler__",
            "gs": "gs_handler__", "H0": "H0_handler__"}
for prefix, vname in handlers.items():
    pattern = r'get_{}\([^)]*handler__[^)]*\);*'.format(prefix)
    matches = re.findall(pattern, content)
    for match in matches:
        if not "({});".format(vname) in match:
            sys.exit("sirius_api.cpp: some handler name does not follow conventions ({}).".format(prefix))
            
#A regex that matches all SIRIUS API calls: returns a void,
#starts with sirius_, and allows for line breaks
pattern = r'\bvoid\s+(?:\s*\n\s*)?sirius_\w+\s*\([^{]*\{'
matches = re.findall(pattern, content, re.DOTALL)

signatures = []
for match in matches:
    signatures.append(match.strip().replace("\n{", ";\n\n"))

with open("sirius_c_headers.h", "w") as myfile:
    myfile.write("#include <stdbool.h>\n")
    myfile.write("#include <complex.h>\n\n")
    for signature in signatures:
        myfile.write(signature)