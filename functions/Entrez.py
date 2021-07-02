#!/usr/bin/python3

"""Entrez

  Version 1.0 for Py3+ 10/05/2021
  
  Paulo Cseri Ricardo
  Contact Info: cseri.bio@gmail.com
"""

def Entrez(tax):
	from Bio import Entrez
	taxid = tax
	#set email
	Entrez.email = "youremail@gmail.com"
	handle = Entrez.efetch(db="taxonomy", id=r.taxid, mode="text", rettype="xml")
	records = Entrez.read(handle)
	lineage = records[0]['Lineage']
	lineage = lineage.split("; ")
	return lineage
