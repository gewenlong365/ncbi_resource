#!/usr/bin/env python
'''
      Copyright (C) Hangzhou
      作    者: 葛文龙
      通讯邮件: gwl9505@163.com
      脚本名称: ref_find.pl
      版    本: 1.0
      创建日期: 2023年05月04日
'''
import sqlite3
import gzip
import getopt
import sys


in_arg = ''
out_file = ''
set_level = ''
db_set = ''
out_level = False
help_arg = False

try:
    opts, args = getopt.getopt(sys.argv[1:], "hi:o:l:d:s", ["help"])
except getopt.GetoptError as err:
    print(err)
    sys.exit(2)

for opt, arg in opts:
    if opt in ("-h", "--help"):
        help_arg = True
    elif opt == "-i":
        in_arg = arg
    elif opt == "-o":
        out_file = arg
    elif opt == "-l":
        set_level = arg.upper()
    elif opt == "-d":
        db_set = arg
    elif opt == "-s":
        out_level = True

if help_arg:
    print("用法: {} -i (taxid) -l (层级) -o ./".format(sys.argv[0]))
    print("选项")
    print("\t-i : 输入taxid")
    print("\t-o : 基因组输出的文件夹，默认对应层级taxid_对应层级名")
    print("\t-l : 选择输出层级，可选以下参数(大小写兼容，默认属级)")
    print("\t\t   K   : 界")
    print("\t\t   P   : 门")
    print("\t\t   C   : 纲")
    print("\t\t   O   : 目")
    print("\t\t   F   : 科")
    print("\t\t   G   : 属")
    print("\t\t   S   : 种")
    print("\t\t   T   : 株\n")
    print("无参数项")
    print("\t-s : 输出层级关系文件到当前目录")
    sys.exit()


# 进入数据库
conn = sqlite3.connect('ncbi_resource.sqlite')
# # 启用外键约束
conn.execute('PRAGMA foreign_keys = ON')
# 指定游标对象
c = conn.cursor()

# nodes表id查询
ctaxid = "SELECT parent_tax_id FROM nodes WHERE tax_id = ?"
all_cid = "SELECT tax_id FROM nodes WHERE  parent_tax_id = ?"

# 展开前置id
now_cid = "573"
prefix = now_cid

while now_cid != 1:
    input = (now_cid,)
    # 查询数据库
    c.execute(ctaxid, input)
    rows = c.fetchall()
    # 输出查询结果
    for row in rows:
        now_cid = row[0]
        # 组合到前缀id中
        prefix = "{}|{}".format(now_cid, prefix)


print(prefix)
# 展开后缀id



# 关闭连接
c.close()