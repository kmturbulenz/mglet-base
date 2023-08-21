import json

# Read entire expression in a string
with open('expression.txt', mode='r') as fh:
    expression = fh.read()

# Read JSON parameters.json
with open('parameters.json', mode='r') as fh:
    prmdict = json.load(fh)

# Set expression for u and w
prmdict['flow']['uinf'][0] = expression + "uu;\n"
prmdict['flow']['uinf'][2] = expression + "ww;\n"

# Write JSON to file
with open('parameters.json', mode='w') as fh:
    json.dump(prmdict, fh, indent=4)
