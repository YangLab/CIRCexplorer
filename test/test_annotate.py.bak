'''
test_annotate.py: Test annotate_fusion
'''

import os
from utils import check_file
from circ.CIRCexplorer import annotate_fusion


class TestAnnotate(object):

    def setup(self):
        '''
        Run annotate_fusion
        '''
        print('#%s: Start testing annotate_fusion' % __name__)
        ref = 'data/ref.txt'
        input = 'data/fusion_junction.txt'
        output = 'data/test_annotated_junction.txt'
        annotate_fusion(ref, input, output)

    def testAnnotate(self):
        '''
        Check file
        '''
        print('#%s: Test annotate_fusion' % __name__)
        test_file = 'data/test_annotated_junction.txt'
        result_file = 'data/annotated_junction.txt'
        check_file(test_file, result_file)

    def teardown(self):
        '''
        Delete annotate file
        '''
        print('#%s: End testing annotate_fusion' % __name__)
        os.remove('data/test_annotated_junction.txt')
