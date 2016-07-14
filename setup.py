#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

setup(name='CIRCexplorer',
      version='1.1.10',
      description='A combined strategy to identify circular RNAs (circRNAs and ciRNAs)',
      author='Xiao-Ou Zhang',
      author_email='zhangxiaoou@picb.ac.cn',
      maintainer='Xu-Kai Ma',
      maintainer_email='maxukai@picb.ac.cn',
      url='https://github.com/YangLab/CIRCexplorer',
      license='MIT',
      packages=find_packages(),
      install_requires=[
          'pysam>=0.8.2',
          'docopt'
      ],
      entry_points={
          'console_scripts': [
              'CIRCexplorer.py=circ.CIRCexplorer:main',
              'fetch_ucsc.py=circ.fetch_ucsc:main',
              'star_parse.py=circ.star_parse:main'
          ],
      },
      )
