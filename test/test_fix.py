'''
test_fix.py: Test fix_fusion
'''

import os
import pysam
from utils import check_file
from circ.CIRCexplorer import fix_fusion


class TestFix(object):

    def setup(self):
        '''
        Run fix_fusion
        '''
        print('#%s: Start testing fix_fusion' % __name__)
        ref = 'data/ref.txt'
        genome = pysam.FastaFile('data/chr21.fa')
        input = 'data/annotated_junction.txt'
        output = 'data/test_circular_RNA.txt'
        fix_fusion(ref, genome, input, output, False)

    def testFix(self):
        '''
        Check file
        '''
        print('#%s: Test fix_fusion' % __name__)
        test_file = 'data/test_circular_RNA.txt'
        result_file = 'data/circular_RNA.txt'
        check_file(test_file, result_file)

    def teardown(self):
        '''
        Delete fix file
        '''
        print('#%s: End testing fix_fusion' % __name__)
        os.remove('data/test_circular_RNA.txt')
