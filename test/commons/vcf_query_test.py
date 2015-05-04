"""
The MIT License (MIT) for varex package

Copyright (c) 2015 Wei-Yi Cheng

Permission is hereby granted, free of charge, to any person obtaining a copy 
of this software and associated documentation files (the "Software"), to deal 
in the Software without restriction, including without limitation the rights 
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
copies of the Software, and to permit persons to whom the Software is 
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
SOFTWARE.
"""

import unittest
import varex.commons.vcf_query as vcfq

class VcfQueryTestCase(unittest.TestCase):
    """
    Base TestCase class for testing varex.commons.vcf_query
    """
    def setUp(self):
        self.ve = vcfq.VcfExtractor(debug=True)
        self.pos_file = "test/resources/test_positions.txt"
        self.vcf_file = "test/resources/test_vcf_list.txt"

class PositionsTestCase(VcfQueryTestCase):
    """
    Test the vcf_query.load_position function
    """
    def runTest(self):
        pos = self.ve.load_positions(self.pos_file)
        self.assertEqual(pos[0], ['1', '10506', '10506'], 
            "pos[0] was not load correctly!")
        self.assertEqual(pos[1], ['1', '15274', '15274'],
            "pos[1] was not load correctly!")


if __name__ == '__main__':
    unittest.main()

