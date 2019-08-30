#!/usr/bin/env python3

# Author: Anna Schroder Lassen
# Date: August 29, 2019

#####################################################
# Takes a file with pathways and names of xlsx-files
# and combines them in a new excel sheet as output.
# The output is also saved as a csv file.
# Finally, the data is printed to the screen.
#####################################################

# Import modules
import sys
import pandas as pd

# Get file with excel sheets
if len(sys.argv) == 1:
    filename = input('Enter a file containing pathways and names of excel files: ')
elif len(sys.argv) == 2:
    filename = sys.argv[1]
    
# Open file
try:
    file_list = open(filename, 'r')
except IOError as error:
    print('Could not open the chosen file, reason: ', str(error))
    sys.exit(1)
    
# Create data frame
df = pd.DataFrame()

# Add sheets to data frame
print('Adding excel sheets ...')
for file in file_list:
    try:
        data = pd.read_excel(file[:-1])
        df = df.append(data)
        print('Added sheet', file[:-1])
    except FileNotFoundError as error:
        print('Could not find sheet:', file[:-1])
        
# Close infile
file_list.close()

# Convert data frame to excel and csv
df.to_excel('/home/projects/HT2_leukngs/people/s173461/output/all_excel_sheets.xlsx')
df.to_csv('/home/projects/HT2_leukngs/people/s173461/output/all_excel_sheets.csv')

# Show excel sheet
print(df)