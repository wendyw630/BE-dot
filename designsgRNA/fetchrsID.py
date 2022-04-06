############################################################
__author__ = "Shiran Abadi"
################### fetch genomic_seq  ######################
from Bio.Seq import Seq
from xml.dom.minidom import parseString
import urllib.request
import os, re

# from Bio.Alphabet import generic_dna 

def download_file(url, filename):
	# print >> sys.stderr, "filename:" + filename
	if not os.path.exists(filename):
		f = urllib.request.urlretrieve(url, filename)
	return filename


def get_xml_filename(genome, chrom, startpos, endpos, cache_dir):
	current_filename = "{0}/{1}_{2}_{3}_{4}.xml".format(cache_dir, genome, chrom, str(startpos), str(endpos))
	return current_filename


def get_dna_coordinates_xmlfile(genome, chromosome, startpos, endpos, cache_dir):
	url = "http://genome.ucsc.edu/cgi-bin/das/{0}/dna?segment={1}:{2},{3}" 

	# chr15:65637530,65637553"
	current_filename = get_xml_filename(genome, chromosome, startpos, endpos, cache_dir)
	# current_filename='user_files/results.csv'
	current_url = url.format(genome, chromosome, startpos, endpos)
	print("current_url:\n"+current_url)
	print("current_filename:\n" + current_filename)
	'''
	current_url
	http://genome.ucsc.edu/cgi-bin/das/hg38/dna?segment=16:77367589,77367639 
	current_filename
	/hg38_16_77367589_77367639.xml

	'''
	download_file(current_url, current_filename)
	return current_filename

### the main func():
## input: genome, chrom, startpos, endpos
## output: DNA-seq
def fetch_dna_coordinates(genome, chrom, startpos, endpos, cache_dir):
	current_filename = get_dna_coordinates_xmlfile(genome, chrom, startpos, endpos, cache_dir)
	with open(current_filename) as fp:
		xmldata = parseString(fp.read())
	seq = re.sub("\s", "", xmldata.childNodes[1].childNodes[1].childNodes[1].childNodes[0].data.upper())
	return seq

'''
def get_seq_by_orientation(seq, strand):
	"""
		returns the original sequence if it's the plus strand, o/w returns the reverse-complement
	"""
	if strand == "-":
		# seqobject = Seq(seq, generic_dna) ### wzl_del
		seqobject = Seq(seq)  ### wzl_add
		seq = str(Seq.reverse_complement(seqobject))
	return seq
'''

#############################################################
######################  single_rsid  ########################
import copy
from Bio import Entrez, SeqIO

#######   Entrez.efetch()
def get_rsIDsnp(snp_id, temp_dir):
	"""
	This function takes as input a snp identifiers and a temporary directory to save temp files to
	and returns a parsed dictionary of their data from Entrez.
	"""
	Entrez.email = 'A.N.Other@example.com'

	try:
		handle = Entrez.efetch(db='SNP', id=snp_id)  
		snp_info = handle.read()
		handle.close()
	except:
		print("Failed to fetch SNP.")
		return None

	r = {} # Return dictionary variable
	field_names = ["SNP_ID", "CLINICAL_SIGNIFICANCE", "GENE_ID", "CHRPOS", "FXN_CLASS", "DOCSUM"]
	for field_name in field_names:
		try:
			r[field_name] = re.search("<" + field_name + ">(.*?)</" + field_name + ">", snp_info).group(1)
		except:
			r[field_name] = None

	res = {}
	res["formal_id"] = r["SNP_ID"]
	res["ClinicalSignificance"] = r["CLINICAL_SIGNIFICANCE"]
	res["GENE_ID"] = r["GENE_ID"]
	try:
		loc_attributes = re.search("(.*):(.*)", r["CHRPOS"])
	except:
		return []
	res["chromosome"] = loc_attributes.group(1)
	res["coordinate"] = int(loc_attributes.group(2))
	res["fxnClass"] = r["FXN_CLASS"]

	res["genomic_flanking"] = \
			fetch_dna_coordinates("hg38", res["chromosome"],
			                              res["coordinate"] - 50, res["coordinate"] + 50, temp_dir)

	## DOCSUM may have multiple mutations due to alternative splicing, and might look like:
	## <DOCSUM>HGVS=NC_000011.10:g.13492506G&gt;A,NC_000011.10:g.13492506G&gt;T,NC_000011.9:g.13514053G&gt;
	# A,NC_000011.9:g.13514053G&gt;T,NG_008962.1:g.8515C&gt;T,NG_008962.1:g.8515C&gt;A,NM_000315.4:c.247C&gt;
	# T,NM_000315.4:c.247C&gt;A,NM_000315.3:c.247C&gt;T,NM_000315.3:c.247C&gt;A,NM_000315.2:c.247C&gt;
	# T,NM_000315.2:c.247C&gt;A,NM_001316352.1:c.343C&gt;T,NM_001316352.1:c.343C&gt;
	# A,NP_000306.1:p.Arg83Ter,NP_001303281.1:p.Arg115Ter|SEQ=[G/A/T]|GENE=PTH:5741</DOCSUM>
	## So I take every record that begins with a nucleotide and then NM (mRNA record) and compute the reading frame according to the position that follows.
	## For example, A,NM_001316352.1:c.343C&gt means that there was an SNV to A at position 343 of the mRNA

	m = re.finditer("([ACGT]),(NM_.*?):c\.([0-9\-]+)", r["DOCSUM"])
	orig_nuc = re.search("SEQ\=\[([ACGT])", r["DOCSUM"]).group(1)  
	mutations = {}
	for match in m:
		try:
			mut = (match.group(1), int(match.group(3)))
			mutations[mut] = mutations.get(mut, []) + [match.group(2)]
			### wzl_note: eg. of mutations:
			### {('C', 605): ['NM_199355.4', 'NM_199355.3', 'NM_199355.2', 'NM_139054.2', 'NM_139054.1'],
			# ('C', 85): ['NM_001326358.2', 'NM_001326358.1']}
		except:
			pass # if we reached here, the coordinate is in the UTRs (noncoding part that might be in the NM query)

	# retrieve a valid codon from one of the mrna sequences
	#wzl_note: for v in ['NM_199355.4', 'NM_199355.3', 'NM_199355.2', 'NM_139054.2', 'NM_139054.1']:
	all_snps = []
	for k, v in mutations.items():
		try:
			handle = Entrez.efetch(db='nucleotide', id=",".join(v), rettype="gb", retmode="text")  
			# "gb": a seq format, like fasta
			response_item = SeqIO.parse(handle, "gb").__next__() #only need one
			for feature in response_item.features:
				if feature.type == "CDS":
					cds_start_index = feature.location.start
					cds_end_index = feature.location.end
					break

			mrna = response_item.seq[cds_start_index:cds_end_index]
			snp_idx = k[1]-1

			all_snps.append(copy.deepcopy(res))
			all_snps[-1]["5UTR"] = k[1] < 0
			all_snps[-1]["3UTR"] = k[1] > (cds_end_index - cds_start_index)

			if not (all_snps[-1]["5UTR"] or all_snps[-1]["3UTR"]): 
				orig_codon = mrna[(snp_idx // 3) * 3:(snp_idx // 3 + 1) * 3]
				all_snps[-1]["original_codon"] = str(orig_codon)
				all_snps[-1]["reading_frame"] = snp_idx % 3 + 1
				all_snps[-1]["orientation"] = "+" if mrna[snp_idx]==orig_nuc else "-"
				all_snps[-1]["SNP"] = mrna[snp_idx] + ">" + k[0]
				all_snps[-1]["CDS_flanking"] = mrna[max(0,snp_idx-50):min(snp_idx+52, len(mrna))]
			else:
				all_snps[-1]["SNP"] = response_item.seq[cds_start_index+snp_idx+1] + ">" + k[0]

			handle.close()
		except:
			pass
	return all_snps

"""
if __name__ == '__main__':
	snp_id = "rs118192143"
	temp_dir = ""
	snp_info = get_rsIDsnp(snp_id, temp_dir)
	print (snp_info)
"""


# eg. of snp_info: rs28937596
"""
[{'formal_id': '28937596', 'ClinicalSignificance': 'pathogenic', 'GENE_ID': '8893', 'chromosome': '3',
  'coordinate': 184144111, 'fxnClass': 'missense_variant,coding_sequence_variant,non_coding_transcript_variant',
  'genomic_flanking': 'CTTTCTTCCATAGCTGCTAAAGGCCTGGAGCCCTGTTTTTAGGAACTACAT', '5UTR': False, '3UTR': False,
  'original_codon': 'TGG', 'reading_frame': 1, 'orientation': '+', 'SNP': 'T>C',
  'CDS_flanking': Seq('CCTGCTGCTTCCTCTGCTAAAGGCCTGGAGCCCTGTTTTTAGGAACTACATA')}]
"""
# eg. of snp_info: rs118192143
"""
[{'formal_id': '118192143', 'ClinicalSignificance': 'not-provided,pathogenic', 'GENE_ID': '6261', 'chromosome': '19', 'coordinate': 38580395, 'fxnClass': 'coding_sequ
ence_variant,genic_downstream_transcript_variant,missense_variant', 'genomic_flanking': 'CTGGTGATGACCGTGGGCCTTCTGGCGGTGGTCGTCTACCTGTACACCGTG', '5UTR': False, '3UTR':
False, 'original_codon': 'GCG', 'reading_frame': 2, 'orientation': '+', 'SNP': 'C>T', 'CDS_flanking': Seq('CTGGTGATGACCGTGGGCCTTCTGGCGGTGGTCGTCTACCTGTACACCGTGG')}, 
{'formal_id': '118192143', 'ClinicalSignificance': 'not-provided,pathogenic', 'GENE_ID': '6261', 'chromosome': '19', 'coordinate': 38580395, 'fxnClass': 'coding_sequenc
e_variant,genic_downstream_transcript_variant,missense_variant', 'genomic_flanking': 'CTGGTGATGACCGTGGGCCTTCTGGCGGTGGTCGTCTACCTGTACACCGTG', '5UTR': False, '3UTR': Fal
se, 'original_codon': 'GCG', 'reading_frame': 2, 'orientation': '+', 'SNP': 'C>G', 'CDS_flanking': Seq('CTGGTGATGACCGTGGGCCTTCTGGCGGTGGTCGTCTACCTGTACACCGTGG')}, 
{'formal_id': '118192143', 'ClinicalSignificance': 'not-provided,pathogenic', 'GENE_ID': '6261', 'chromosome': '19', 'coordinate': 38580395, 'fxnClass': 'coding_sequence_v
ariant,genic_downstream_transcript_variant,missense_variant', 'genomic_flanking': 'CTGGTGATGACCGTGGGCCTTCTGGCGGTGGTCGTCTACCTGTACACCGTG', '5UTR': False, '3UTR': False,
 'original_codon': 'GCG', 'reading_frame': 2, 'orientation': '+', 'SNP': 'C>T', 'CDS_flanking': Seq('CTGGTGATGACCGTGGGCCTTCTGGCGGTGGTCGTCTACCTGTACACCGTGG')}, 
 {'formal_id': '118192143', 'ClinicalSignificance': 'not-provided,pathogenic', 'GENE_ID': '6261', 'chromosome': '19', 'coordinate': 38580395, 'fxnClass': 'coding_sequence_vari
ant,genic_downstream_transcript_variant,missense_variant', 'genomic_flanking': 'CTGGTGATGACCGTGGGCCTTCTGGCGGTGGTCGTCTACCTGTACACCGTG', '5UTR': False, '3UTR': False, 'o
riginal_codon': 'GCG', 'reading_frame': 2, 'orientation': '+', 'SNP': 'C>G', 'CDS_flanking': Seq('CTGGTGATGACCGTGGGCCTTCTGGCGGTGGTCGTCTACCTGTACACCGTGG')}]
"""
