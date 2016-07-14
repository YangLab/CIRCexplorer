## 1.1.10 (2016-7-14)

Bugfixes:

* case-sensitivity issue for fusion fix step

Improvements:

* use smaller test data

## 1.1.9 (2016-6-7)

Improvements:

* offer more deployment methods (PyPI, bioconda, Docker and Galaxy Tool). Thanks
@bgruening and @anatskiy for their helps!

## 1.1.8 (2016-4-21)

Improvements:

* add test codes

## 1.1.7 (2016-1-28)

Improvements:

* add pysam version checking at running time

## 1.1.6 (2016-1-22)

Improvements:

* add exon/intron index information

## 1.1.5 (2016-1-20)

Improvements:

* use setuptools instead of distutils in setup.py, fix #11

## 1.1.4 (2016-1-16)

Improvements:

* check pysam version at running time
* index fa file automatically if it has not been indexed
* integrate interval.py into it

## 1.1.3 (2015-7-31)

Bugfixes:

* bug that removes fusion junction reads aligned to first/last exons

Improvements:

* support pysam 0.8.2

## 1.1.2 (2015-6-18)

Improvements:

* add script fetch_ref.py

## 1.1.1 (2015-3-13)

Improvements:

* adjust some misleading labels

## 1.1.0 (2015-2-3)

Bugfixes:

* bug that causes a crash when parsing ciRNA junctions

Improvements:

* support aligner STAR

## 1.0.6 (2014-12-11)

Bugfixes:

* bug that prints redundant circular RNA information

## 1.0.5 (2014-11-23)

Improvements:

* add option '--no-fix' to not use fix step

## 1.0.4 (2014-10-24)

Improvements:

* add option '--tmp' to keep temporary files
* change option '--output' from output file name to output file prefix

## 1.0.3 (2014-10-11)

Improvements:

* remove redundant source codes
* modify source codes according to PEP8

## 1.0.2 (2014-09-15)

Bugfixes:

* bug that ignores mismatches of fusion junction reads near the first or last exon
* bug that removes fusion junction reads due to align to a random isoform

## 1.0.0 (2014-09-12)

Features:

* first release for CIRCexplorer
* all the source codes were derived from the CELL paper and reorganized for general usage
