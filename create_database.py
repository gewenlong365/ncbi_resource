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
# 启用外键约束
conn.execute('PRAGMA foreign_keys = ON')
# 指定游标对象
c = conn.cursor()

# 定义assembly_summary表结构
def assembly_summary(source):
    create_table_sql = '''CREATE TABLE {} (
                            assembly_accession text PRIMARY KEY  NOT NULL,
                            bioproject text,
                            biosample text,
                            wgs_master text,
                            refseq_category text,
                            taxid int,
                            species_taxid int,
                            organism_name text,
                            infraspecific_name text,
                            isolate text,
                            version_status text,
                            assembly_level text,
                            release_type text,
                            genome_rep text,
                            seq_rel_date text,
                            asm_name text,
                            submitter text,
                            gbrs_paired_asm text,
                            paired_asm_comp text,
                            ftp_path text,
                            excluded_from_refseq text,
                            relation_to_type_material text,
                            asm_not_live_date text
                            );'''.format(source)
    c.execute(create_table_sql)


names = ["assembly_summary_refseq", "assembly_summary_genbank"]
for name in names:
    assembly_summary(name)

# 建立nodes.dmp表
c.execute('''CREATE TABLE nodes
             (tax_id INTEGER PRIMARY KEY,
              parent_tax_id INTEGER,
              rank TEXT,
              embl_code TEXT,
              division_id INTEGER,
              inherited_div_flag INTEGER,
              genetic_code_id INTEGER,
              inherited_GC_flag INTEGER,
              mitochondrial_genetic_code_id INTEGER,
              inherited_MGC_flag INTEGER,
              GenBank_hidden_flag INTEGER,
              hidden_subtree_root_flag INTEGER,
              comments TEXT)''')

# 建立names.dmp表
c.execute('''CREATE TABLE names
             (tax_id INTEGER,
              name TEXT,
              unique_name TEXT,
              name_class TEXT)''')


# 建立基因组资源表
c.execute('''CREATE TABLE genomo2gene (
    assembly_accession TEXT PRIMARY KEY REFERENCES assembly_summary_refseq(assembly_accession) ON DELETE CASCADE,
    genome_body BLOB
    );''')

# 基因组与基因分离
# c.execute('''CREATE TABLE genomo2gene
#              (genome_id TEXT,
#               gene_id TEXT,
#               CONSTRAINT ref基因组对列表 FOREIGN KEY (genome_id)  REFERENCES assembly_summary_refseq(assembly_accession)
#               )''')
# c.execute('''CREATE TABLE genomo2gene (
#     genome_id TEXT REFERENCES assembly_summary_refseq(assembly_accession) ON DELETE CASCADE,
#     gene_id TEXT,
#     PRIMARY KEY (genome_id, gene_id)
#     );''')


# 创建触发器,删除ref表记录时删除基因组
# c.execute('''CREATE TRIGGER delete_genomo
#             AFTER DELETE ON assembly_summary_refseq
#             BEGIN
#                 DELETE FROM genomo2gene WHERE assembly_accession = OLD.assembly_accession;
#             END;''')


# CREATE TRIGGER delete_genomo2gene
# AFTER DELETE ON assembly_summary_refseq
# BEGIN
#     DELETE FROM genomo2gene WHERE genome_id = OLD.assembly_accession;
# END;

# 将assembly_summary_refseq.txt存入assembly_summary_refseq表
# 将assembly_summary_genbank.txt存入assembly_summary_genbank表
for name in names:
    name_file = name + ".txt"
    with open(name_file, 'r') as f:
        reader = csv.reader(filter(lambda row: row[0] != '#', f), delimiter='\t')
        for row in reader:
            c.execute(
                f"INSERT INTO {name} VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,?)",
                row)

# 查询句式：SELECT taxid FROM assembly_summary_refseq where assembly_accession = 'GCF_000002865.3';


with open('nodes.dmp', 'r') as f:
    reader = csv.reader(
        filter(
            lambda row: row[0] != '#',
            (line.replace(
                '\t|\t',
                '\t') for line in f)),
        delimiter='\t')
    for row in reader:
        row = row[:-1]
        c.execute(
            "INSERT INTO nodes VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            row)
# 查询句式：SELECT parent_tax_id FROM nodes where tax_id = '683';


with open('names.dmp', 'r') as f:
    reader = csv.reader(
        filter(
            lambda row: row[0] != '#',
            (line.replace(
                '\t|\t',
                '\t') for line in f)),
        delimiter='\t')
    for row in reader:
        row = row[:-1]
        c.execute("INSERT INTO names VALUES (?, ?, ?, ?)", row)
# 查询句式：SELECT name FROM names where name_class = 'scientific name' and tax_id = '7' limit 1;

conn.commit()
conn.close()

#