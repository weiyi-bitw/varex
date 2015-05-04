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

import sys
import os
import re
import subprocess

class VcfExtractor(object):
  """
  VcfExtractor

  This class extract from VCF files a list of positions and convert the 
  genotypes to a matrix of aac count. 
  """
  def __init__(self, 
               debug=False,
               log=sys.stderr, 
               clear_garbage=False):
      """
      Initializer of VcfExtractor
      
      debug: run using debug mode, which will give detailed output to 
          stderr
      log: output stream of information
      clear_garbage: delete temporary files before destroy the class
      """
      self.debug = debug
      self.log = log
      self.__garbage_list__ = []
      self.clear_garbage=clear_garbage
  
  def load_positions(self, position_file):
      """
      load position files and output a list of positions
      
      position_file: input file containing a list of positions
          format: chrom, start, end, ...
      """
      self.__position_file__ = position_file
      with open(position_file, "rb") as fi:
          return map(lambda l: l.strip().split('\t', 3)[:3], fi)


  def query_vcf(self, position_file, vcf_files):
      """
      IMPORTANT: before calling this function, make sure you have 
          vcftools and tabix working in your shell environment
      given a list of positions [chrom, start, end], output a file of 
          query output. 
      
      position_file: file of positions to be queried, in the format of 
          [chrom, start, end]
      vcf_files: list of paths to vcf files
      """
      #subprocess.call(["module load vcftools/0.1.12b"], shell=True)
      #subprocess.call(["module load tabix"], shell=True)
      positions = self.load_positions(position_file)
      out_file='.'.join(position_file.split('.')[:-1]) + ".qout"
      self.__garbage_list__.append(out_file)
      with open(out_file, "wb") as fo:
          self.__vcf_files_source__ = vcf_files
          self.__vcf_files__ = []
          for vf in open(vcf_files, "rb"):
              vf = vf.strip()
              self.__vcf_files__.append(vf)
              if self.debug:
                  print >> self.log, "Querying file %s ..."%(vf)
              for p in positions:
                  q = p[0] + ":" + p[1] + "-" + p[2]
                  if self.debug:
                      #print >> self.log, "current directory: %s"%(os.getcwd())
                      print >> self.log, "    Querying position %s ..."%(q)
                  subprocess.call(["vcf-query", vf, q], stdout=fo)

  def query_out_to_matrix(self, query_out_file):
      """
      convert the vcf-query output to matrix format

      query_out_file: file contains query output from vcf-query, can be the 
          output file from self.query_vcf
      """
      with open(query_out_file, "rb") as fi:
          self.__samples__ = set()
          self.__variants__ = set()
          self.genotypes = {}
          for l in fi:
              [chrom_pos, ref, gt] = l.strip().split("\t", 2)
              [chrom, pos] = chrom_pos.split(":")
              tmp_list = [ele.split("=") for ele in gt.split("\t")]
              self.__samples__ = \
                  self.__samples__.union(set([s[0] for s in tmp_list]))
              alts = set([a for tt in tmp_list 
                  for a in re.split("/|\|", tt[1])])
              alts.remove(ref)
              for alt in alts:
                  variant = "_".join([chrom, pos, ref, alt])
                  #print >> self.log, variant
                  self.__variants__.add(variant)
                  tmp_dict = dict(
                      [[ele[0], ele[1].count(alt) if "." not in ele else "NA"] 
                        for ele in tmp_list]
                  )
                  if variant not in self.genotypes:
                      self.genotypes[variant] = tmp_dict
                  else:
                      self.genotypes[variant].update(tmp_dict)
      
      out_file='.'.join(query_out_file.split('.')[:-1]) + ".matrix.txt"
      final_samples = list(self.__samples__)
      if self.debug:
          print >> self.log, "Number of samples: %d"%(len(final_samples))
          print >> self.log, "Number of variants: %d"%(len(self.__variants__))
      with open(out_file, "wb") as fo:
          print >> fo, "#VAR\\SAMPLE\t" + "\t".join(final_samples)
          for v in self.__variants__:
              t = [str(self.genotypes[v][s]) \
                   if s in self.genotypes[v] else "NA" \
                   for s in final_samples]
              print >> fo, v + "\t" + "\t".join(t)


