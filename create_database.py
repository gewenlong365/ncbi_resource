#!/usr/bin/env python
'''
      Copyright (C) Hangzhou
      作    者: 葛文龙
      通讯邮件: gwl9505@163.com
      脚本名称: create_database.py.pl
      版    本: 1.0
      创建日期: 2023年04月20日
'''
import csv
import sqlite3
import zlib
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
import glob
import os
import re

# 创建数据库
conn = sqlite3.connect('ncbi_resource.sqlite')
c = conn.cursor()

# 将assembly_summary_refseq.txt存入assembly_summary_refseq表
# c.execute('''CREATE TABLE assembly_summary_refseq
#              (assembly_accession text PRIMARY KEY   NOT NULL,
#              bioproject text,
#              biosample text,
#              wgs_master text,
#              refseq_category text,
#              taxid int,
#              species_taxid int,
#              organism_name text,
#              infraspecific_name text,
#              isolate text,
#              version_status text,
#              assembly_level text,
#              release_type text,
#              genome_rep text,
#              seq_rel_date text,
#              asm_name text,
#              submitter text,
#              gbrs_paired_asm text,
#              paired_asm_comp text,
#              ftp_path text,
#              excluded_from_refseq text,
#              relation_to_type_material text,
#              asm_not_live_date text)''')


# with open('assembly_summary_refseq.txt', 'r') as f:
#     reader = csv.reader(filter(lambda row: row[0] != '#', f), delimiter='\t')
#     for row in reader:
        # c.execute("INSERT INTO assembly_summary_refseq VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,?)", row)
# 查询句式：SELECT taxid FROM assembly_summary_refseq where assembly_accession = 'GCF_000002865.3';

# 将assembly_summary_genbank.txt存入assembly_summary_genbank表
# c.execute('''CREATE TABLE assembly_summary_genbank
#              (assembly_accession text PRIMARY KEY   NOT NULL,
#              bioproject text,
#              biosample text,
#              wgs_master text,
#              refseq_category text,
#              taxid int,
#              species_taxid int,
#              organism_name text,
#              infraspecific_name text,
#              isolate text,
#              version_status text,
#              assembly_level text,
#              release_type text,
#              genome_rep text,
#              seq_rel_date text,
#              asm_name text,
#              submitter text,
#              gbrs_paired_asm text,
#              paired_asm_comp text,
#              ftp_path text,
#              excluded_from_refseq text,
#              relation_to_type_material text,
#              asm_not_live_date text)''')
#
# with open('assembly_summary_genbank.txt', 'r') as f:
#     reader = csv.reader(filter(lambda row: row[0] != '#', f), delimiter='\t')
#     for row in reader:
#         c.execute("INSERT INTO assembly_summary_genbank VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,?)", row)
# 查询句式：SELECT taxid FROM assembly_summary_genbank where assembly_accession = 'GCA_000001895.4';


# 将nodes.dmp存入nodes表
# c.execute('''CREATE TABLE nodes
#              (tax_id INTEGER PRIMARY KEY,
#               parent_tax_id INTEGER,
#               rank TEXT,
#               embl_code TEXT,
#               division_id INTEGER,
#               inherited_div_flag INTEGER,
#               genetic_code_id INTEGER,
#               inherited_GC_flag INTEGER,
#               mitochondrial_genetic_code_id INTEGER,
#               inherited_MGC_flag INTEGER,
#               GenBank_hidden_flag INTEGER,
#               hidden_subtree_root_flag INTEGER,
#               comments TEXT)''')

# with open('nodes.dmp', 'r') as f:
#     reader = csv.reader(filter(lambda row: row[0] != '#', (line.replace('\t|\t', '\t') for line in f)), delimiter='\t')
#     for row in  reader:
#         row = row[:-1]
#         c.execute("INSERT INTO nodes VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", row)
# 查询句式：SELECT parent_tax_id FROM nodes where tax_id = '683';

# 将names.dmp存入names表
# c.execute('''CREATE TABLE names
#              (tax_id INTEGER,
#               name TEXT,
#               unique_name TEXT,
#               name_class TEXT)''')

# with open('names.dmp', 'r') as f:
#     reader = csv.reader(filter(lambda row: row[0] != '#', (line.replace('\t|\t', '\t') for line in f)), delimiter='\t')
#     for row in  reader:
#         row = row[:-1]
#         c.execute("INSERT INTO names VALUES (?, ?, ?, ?)", row)
# 查询句式：SELECT name FROM names where name_class = 'scientific name' and tax_id = '7' limit 1;

# 存入序列文件到seq_info表中
# c.execute('''CREATE TABLE seq_info
#              (gene_id TEXT,
#               genome_id TEXT,
#               gene_head TEXT,
#               compress_seq TEXT)''')

path = './out_ref_down'
gz_files = glob.glob(os.path.join(path, '*.gz'))

for file in gz_files:
    # 抓取基因组名
    pattern = r'^\./out_ref_down/(GCF_\d+\.\d+)_.*$'
    result = re.match(pattern, file)
    if result:
        genome_id = result.group(1)

    with gzip.open(file, 'rt') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            gene_id = record.id
            gene_info = record.description
            seq = Seq(record.seq)
            compressed_seq = zlib.compress(bytes(seq))
            c.execute("INSERT INTO seq_info VALUES (?, ?, ?, ?)", (gene_id,genome_id,gene_info,compressed_seq))

conn.commit()
conn.close()