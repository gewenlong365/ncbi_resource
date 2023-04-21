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

# 创建数据库
conn = sqlite3.connect('ncbi_resource.sqlite')
c = conn.cursor()

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


with open('assembly_summary_refseq.txt', 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        c.execute("INSERT INTO assembly_summary_refseq VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,?)", row)

conn.commit()
conn.close()