"""
The MIT License (MIT) for varex package

Copyright (c) 2015 Wei-Yi Cheng

Permission is hereby granted, free of charge, to any person obtaining a copy 
of this software and associated documentation files (the "Software"), to deal 
in the Software without restriction, including without limitation the rights 
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
copies of the Software, and to permit persons to whom the Software is furnished 
to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

import sys
import subprocess

def load_positions(position_file):
    """
    load position files and output a list of positions
    
    position_file: input file containing a list of positions
        format: chrom, start, end, ...
    """
    with open(position_file, "rb") as fi:
        return map(lambda l: l.strip().split('\t', 3)[:3], fi)


def query_vcf(position_file, vcf_files):
    """
    given a list of positions [chrom, start, end], output a file of 
        query output
    
    position_file: file of positions to be queried, in the format of 
        [chrom, start, end]
    vcf_files: list of paths to vcf files
    """
    subprocess.call(["module load vcftools/0.1.12b"], shell=True)
    subprocess.call(["module load tabix"], shell=True)
    positions = load_positions(position_file)
    out_file='.'.join(position_file.split('.')[:-1]) + ".genotype.txt"
    with open(out_file, "wb") as fo:
        for vf in open(vcf_files, "rb"):
            vf = vf.strip()
            print >> sys.stderr, "Querying file " + vf + "..."
            for p in positions:
                q = p[0] + ":" + p[1] + "-" + p[2]
                print >> sys.stderr, "    Querying position " + q + "..."
                subprocess.call(["vcf-query", vf, q], stdout=fo)

