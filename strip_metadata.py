import os

"""
Author: Eric J. Ma

Purpose: To strip metadata from a PDF file, simply open it up using PyPDF2, and save it again.
"""

from PyPDF2 import PdfFileReader, PdfFileMerger

files_dir = os.getcwd()

if 'stripped' not in os.listdir(files_dir):
	os.mkdir('stripped')

for f in os.listdir(files_dir):
	if f.split('.')[-1] == 'pdf':
		merger = PdfFileMerger()
		merger.append(PdfFileReader(f, 'rb'))
		merger.write('stripped/{0}'.format(f))