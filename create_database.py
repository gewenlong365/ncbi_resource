#!/usr/bin/env python
'''
      Copyright (C) Hangzhou
      作    者: 葛文龙
      通讯邮件: gwl9505@163.com
      脚本名称: create_database.py.pl
      版    本: 1.0
      创建日期: 2023年04月20日
'''

import sqlite3

# 创建数据库
conn = sqlite3.connect('ncbi_resource.sqlite')
conn.close()

