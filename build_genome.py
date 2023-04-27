#!/usr/bin/env python
'''
      Copyright (C) Hangzhou
      作    者: 葛文龙
      通讯邮件: gwl9505@163.com
      脚本名称: build_genome.pl
      版    本: 1.0
      创建日期: 2023年04月27日
'''

import gzip
from Bio import SeqIO
from Bio.Seq import Seq
import glob
import os
import re


# 进入数据库
conn = sqlite3.connect('ncbi_resource2.sqlite')
# 启用外键约束
conn.execute('PRAGMA foreign_keys = ON')
# 指定游标对象
c = conn.cursor()

# 存入基因组与基因名对应表
# 获取目录下所有已下载序列
# 用数据库自带压缩
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
            # 抓取基因名，基因名行，基因序列压缩
            gene_id = record.id
            gene_info = record.description
            seq = Seq(record.seq)
            compressed_seq = zlib.compress(bytes(seq))
            c.execute("INSERT INTO genomo2gene VALUES (?, ?)",
                      (genome_id, gene_id))
            c.execute("INSERT INTO gene2seq VALUES (?, ?, ?)",
                      (gene_id, gene_info, compressed_seq))

conn.commit()
conn.close()