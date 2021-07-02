#!/usr/bin/python3

"""NCBITaxa

  Version 1.0 for Py3+ 10/05/2021
  
  Paulo Cseri Ricardo
  Contact Info: cseri.bio@gmail.com
"""

def NCBITaxa(tax):
	from ete3 import NCBITaxa
	ncbi = NCBITaxa()
	taxid = tax
	lineage = ncbi.get_lineage(taxid)
	names = ncbi.get_taxid_translator(lineage)
	u_taxa = []
	for t in lineage:
		u_taxa.append(names[t])
	t_taxa = []
	for u in u_taxa:
		t_taxa.append(str(u))
	return t_taxa
