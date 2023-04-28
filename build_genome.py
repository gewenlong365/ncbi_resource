#!/usr/bin/env python
'''
      Copyright (C) Hangzhou
      作    者: 葛文龙
      通讯邮件: gwl9505@163.com
      脚本名称: build_genome.pl
      版    本: 1.0
      创建日期: 2023年04月27日
'''
import sqlite3
import gzip
import glob
import os
import io
import re

# from Bio import SeqIO
# from Bio.Seq import Seq



# 进入数据库
conn = sqlite3.connect('ncbi_resource.sqlite')
# # 启用外键约束
conn.execute('PRAGMA foreign_keys = ON')
# 指定游标对象
c = conn.cursor()

# 存入基因组与基因名对应表
# 获取目录下所有已下载序列
path = './out_ref_down'
gz_files = glob.glob(os.path.join(path, '*.gz'))


for file in gz_files:
    # 抓取基因组名
    pattern = r'^\./out_ref_down/(GCF_\d+\.\d+)_.*$'
    result = re.match(pattern, file)
    genome_id = None
    if result:
        genome_id = result.group(1)
    # 检查基因组在ref表中是否存在
    c.execute('SELECT assembly_accession FROM assembly_summary_refseq where assembly_accession = ?',(genome_id,))
    result_genome_id = c.fetchall()
    if not result_genome_id :
        continue

    with gzip.open(file, 'rb') as handle:
        data = handle.read()
        compressed_content = sqlite3.Binary(gzip.compress(data))
        c.execute("INSERT INTO genomo2gene VALUES (?, ?)",(genome_id, compressed_content))


# for file in gz_files:
#
#
#
#
#     with gzip.open(file, 'rt') as handle:
#         for record in SeqIO.parse(handle, 'fasta'):
#             # 抓取基因名，基因名行，基因序列压缩
#             gene_id = record.id
#             gene_info = record.description
#             seq = Seq(record.seq)
#             compressed_seq = zlib.compress(bytes(seq))
#             c.execute("INSERT INTO genomo2gene VALUES (?, ?)",
#                       (genome_id, gene_id))
#             c.execute("INSERT INTO gene2seq VALUES (?, ?, ?)",
#                       (gene_id, gene_info, compressed_seq))

conn.commit()
conn.close()

# # 从sqlite数据库中读取压缩后的内容
# c.execute('''SELECT data FROM compressed_data''')
# compressed_content = c.fetchone()[0]
# content = gzip.decompress(compressed_content)