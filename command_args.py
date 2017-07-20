'''
Created on June 2, 2017
This module parses the command line arguments
This change was made on lab system- testing git
@author: ashu
'''
import argparse
import os
#from distutils import spawn
#import subprocess

def parseargs():
    parser = argparse.ArgumentParser(description='Parse command line arguments')
    parser.add_argument('-f',dest='xtc',action='store', help='trajectory file for the analysis')
    parser.add_argument('-s',dest='tpr',action='store', help='structure file for the analysis')
    parser.add_argument('-n',dest='ndx_file',action='store', help='index file')
    parser.add_argument('-skip',dest='skip',action='store', help='frames to skip', default='0', nargs='?', const='0')
    parser.add_argument('-b',dest='beg',action='store', help='starting time to extract in ps', default='0', nargs='?', const='0')
    parser.add_argument('-e',dest='end',action='store', help='Last frame to extract in ps')
    parser.add_argument('-o',dest='outfile_name',action='store', help='output file for the analysis', default='outfile.txt', nargs='?', const='outfile.txt')
    
    args=parser.parse_args()
    if os.path.isfile(args.xtc):
        xtc=os.path.abspath(args.xtc)
    else:
        print(args.xtc + " trajectory file does not exist")
    if os.path.isfile(args.tpr):
        tpr=os.path.abspath(args.tpr)
    else:
        print(args.tpr+ " structure file does not exist" )
    if os.path.isfile(args.ndx_file):
        ndx_file=os.path.abspath(args.ndx_file)
    else:
        print(args.ndx_file+" residue list file does not exist" )
    
    return xtc,tpr,ndx_file,args.skip,args.beg,args.end,args.outfile_name
#def gromacs_exe():
#    if spawn.find_executable("gmx"):
#        return '\'gmx\', \'hbond\''
#        #subprocess.call(['gmx','hbond', '-f', xtc, '-s', tpr, '-num', '-g', '-dist', '-ang', '-hbn', '-hbm', '-life', '-tu', 'ns'],stdin=hbchoice)
#    #subprocess.call(['gmx','hbond', '-f', xtc, '-s', tpr, '-num', '-g', '-dist', '-ang', '-hbn', '-hbm', '-life', '-tu', 'ns', '-b', '0','-e','2'],stdin=hbchoice) ######TEST
#    elif spawn.find_executable("g_hbond"):
#        return 'g_hbond'
#        #subprocess.call(['g_hbond', '-f', xtc, '-s', tpr, '-num', '-g', '-dist', '-ang', '-hbn', '-hbm', '-life', '-tu', 'ns'],stdin=hbchoice)
#    