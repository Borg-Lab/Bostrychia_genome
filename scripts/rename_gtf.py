#!/usr/bin/env python3

import sys

def select_second_character(word):
	spt=word[0].split('_')
	return (int(spt[0]),int(spt[1]))


def sortChromosomes(data):
	with open(data,"r") as file, open('sorted_renamed_chrom.gtf', 'w') as output:
		dictLines={}
		line_num=0
		for line in file:
			fields=line.strip().split('\t')
			chrom_name=fields[0]
			chrom_name_split=chrom_name.split("_")
			chrom_number=chrom_name_split[2]
			chrom_number_digit=("%02d" % (int(chrom_number),))
			fields[0] = f'Chr{chrom_number_digit}'

			line_num += 1
			name=fields[0].split("r")
			Chr_number=name[1]
			feature_type=fields[2]

			dictLines.update({f'{Chr_number}_{line_num}':('\t'.join(fields) + '\n')})
		od = dict(sorted(dictLines.items(), key=select_second_character))
		for k, v in od.items(): output.write(v)

def get_id(spt):
	chs=spt[0]
	if "Vm" in chs:
		chr=42
	elif "Uf" in chs:
		chr=43
	else:
		chr=int(spt[0].replace('Chr',''))
	if spt[2]=='gene':
		return (chr,spt[-1].strip())
	if spt[2]=='transcript':
		return (chr,spt[-1].split('.')[0].strip())
	sp2=spt[-1].split(';')
	for stx in sp2:
		if 'gene_id' in stx:
			return (chr,stx.split('"')[1].strip())

def preprocess_gtf(input_gtf):
	with open(input_gtf, 'r') as input_file, open("adjusted_gene_order.gtf", 'w') as output_file:
		total = {}
		header = ""
		for line in input_file:
			if line.startswith('#'):
				output_file.write(line)
				continue

			fields = line.strip().split('\t')
			spt = line.split('\t')
			chrom_name=fields[0]
			chrom_name_split=chrom_name.split("r")
			chrom_number=chrom_name_split[1]
			feature_type = fields[2]
			attributes = fields[8]

			name = spt[2]
			endc = spt[-1].split(';')
			id = get_id(spt)
			if id not in total:
				total[id] = []
			total[id].append(spt)

		chrcount={}
		for genename in total:
			chromosome = genename[0]
			gene_name_old = genename[1]
			if chromosome not in chrcount:
				chrcount[chromosome] = 0
			chrcount[chromosome] += 10
			genenumber = chrcount[chromosome]
			newname = f'g_{genenumber}'
			for i in range(len(total[genename])):
				entry = total[genename]
				for spt in range(len(entry)):
					total[genename][spt][-1] = total[genename][spt][-1].replace(gene_name_old,newname)

		for id in total:
			for line in total[id]:
				output_file.write('\t'.join(line))

def changeGeneIDs(file,output_gtf):
	with open(file, 'r') as input_file, open(output_gtf, 'w') as output_file:
		for line in input_file:

			fields = line.strip().split('\t')
			spt = line.split('\t')
			chrom_name=fields[0]
			chrom_name_split=chrom_name.split("r")
			chrom_number=chrom_name_split[1]
			feature_type = fields[2]
			attributes = fields[8]

			if feature_type == 'gene':
				gene_split=attributes.split("_")
				gene_id_short=gene_split[1]
				gene_id=("%06d" % ((int(gene_id_short)),))
				attributes_mod = f'gene_id "Bm{chrom_number}g{gene_id}";'
				fields[8] = attributes_mod
				output_file.write('\t'.join(fields) + '\n')

			elif feature_type == 'transcript':
				transcript_split=attributes.split(".")
				transcript_id=transcript_split[1]
				attributes_mod = f'transcript_id "Bm{chrom_number}g{gene_id}.{transcript_id}"; gene_id "Bm{chrom_number}g{gene_id}";'
				fields[8] = attributes_mod
				output_file.write('\t'.join(fields) + '\n')

			else:
				transcript_split=attributes.split(".")
				transcript_id_long=transcript_split[1].split("\"")
				transcript_id=transcript_id_long[0]
				attributes_mod = f'transcript_id "Bm{chrom_number}g{gene_id}.{transcript_id}"; gene_id "Bm{chrom_number}g{gene_id}";'
				fields[8]=attributes_mod
				output_file.write('\t'.join(fields) + '\n')

preprocess_gtf(sys.argv[1])
changeGeneIDs('adjusted_gene_order.gtf', 'preprocessed.gtf')
