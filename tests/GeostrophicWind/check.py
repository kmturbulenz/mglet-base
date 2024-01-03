#
# Checking if the geostrophic equlibrium between Coriolis "force"
# and pressure gradient has been maintained throughout the run.
#
logfile = open('LOGS/uvwbulk.log', 'r')
lines = logfile.readlines()

# ubulk is the third entry in the log file
nubulk = 2
lubulk = []

# reading lines and extracting ubulk
for line in lines[1:]:
    lubulk.append( float(line.split()[nubulk]) )

# throwing error if ubulk has changed
if ( max(lubulk) - min( lubulk ) ) > 0.00001:
  raise ValueError("GeostrophicWind: ubulk should not change during the test. Error.")
