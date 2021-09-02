# from funcs import *
import os
KEEPSAKES = os.listdir()
path = os.path.join(os.getcwd(), "examples/NIMROD/NIMRODsilica/NIMRODmesoporoussilica.txt")
with open(path, "r") as input:
    lines = input.read()
lines = lines.replace("/Users/dtb73/Gudrun2014/Gudrun/run/NIMROD/NIMRODsilica/", os.path.join(os.getcwd(), "examples/NIMROD/NIMRODsilica/"))
lines = lines.replace("/Users/dtb73/Gudrun2014/Gudrun/run/NIMROD/NIMRODsilica/NIMRODsilicadata/", os.path.join(os.getcwd(), "examples/NIMROD/NIMRODsilica/NIMRODsilicadata/"))
lines = lines.replace("/Users/dtb73/Gudrun2014/Gudrun", os.path.abspath("."))
lines = lines.replace("/Users/dtb73/Gudrun2014/Gudrun/StartupFiles/SLS", os.path.abspath('bin/StartupFiles/SLS'))
with open(path, "w") as output:
    output.write(lines)
