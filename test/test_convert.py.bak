'''
test_annotate.py: Test convert_fusion
'''

import os
import pysam
from utils import check_file
from circ.CIRCexplorer import convert_fusion


class TestConvert(object):

    def setup(self):
        '''
        Run convert_fusion
        '''
        print('#%s: Start testing convert_fusion' % __name__)
        fusion_bam = pysam.AlignmentFile('data/tophat_fusion.bam', 'rb')
        output = 'data/test_fusion_junction.txt'
        convert_fusion(fusion_bam, output)

    def testConvert(self):
        '''
        Check file
        '''
        print('#%s: Test convert_fusion' % __name__)
        test_file = 'data/test_fusion_junction.txt'
        result_file = 'data/fusion_junction.txt'
        check_file(test_file, result_file)

    def teardown(self):
        '''
        Delete convert file
        '''
        print('#%s: End testing convert_fusion' % __name__)
        os.remove('data/test_fusion_junction.txt')
