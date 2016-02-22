"""
Author: Eric J. Ma  

Purpose: To merge PDFs together in an automated fashion.
"""


import os

from PyPDF2 import PdfFileReader, PdfFileMerger

files_dir = os.getcwd()

all_files = list()
# Add in main text file.
main_text = [f for f in os.listdir(files_dir) if 'Draft Text' in f and 'pdf' in f]
all_files.extend(main_text)

# Add in Figures.
figures = sorted([f for f in os.listdir(files_dir) if 'Figure ' in f and 'pdf' in f])
all_files.extend(figures)


# Merge the files
merger = PdfFileMerger()
for f in all_files:
	merger.append(PdfFileReader(f), 'rb')

merger.write('EJM_Reassortment.pdf')