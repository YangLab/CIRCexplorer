#!/usr/bin/env python
# -*- coding: utf-8 -*-

from distutils.core import setup

setup(name='CIRCexplorer',
      version='1.1.4',
      description='Circular RNA analysis toolkits',
      author='Xiao-Ou Zhang',
      author_email='zhangxiaoou@picb.ac.cn',
      maintainer='YangLab',
      maintainer_email='rnacomplex@gmail.com',
      url='https://github.com/YangLab/CIRCexplorer',
      license='MIT',
      scripts=['CIRCexplorer.py', 'star_parse.py', 'fetch_ucsc.py'],
      py_modules=['genomic_interval'],
      )
