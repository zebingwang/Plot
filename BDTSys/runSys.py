#!/usr/bin/python
import glob
import sys, os, pwd, commands
import optparse, shlex, re
import time
from time import gmtime, strftime
import math
import subprocess

import datetime
import threading
import multiprocessing

# define function for processing the external os commands
def processCmd(cmd, quite = 0):
    #    print cmd
    status, output = commands.getstatusoutput(cmd)
    if (status !=0 and not quite):
        print 'Error in processing command:\n   ['+cmd+']'
        print 'Output:\n   ['+output+'] \n'
        return "ERROR!!! "+output
    else:
        return output

def execCmd(cmd):
    try:
        print "cmd: %s start running%s" % (cmd,datetime.datetime.now())
        #os.system(cmd)
        subprocess.call(cmd)
        print "cmd: %s end running%s" % (cmd,datetime.datetime.now())
    except Exception, e:
        print '%s\t failed,reason: \r\n%s' % (cmd,e)

def process(*cmd):
    subprocess.call(cmd)

def hadd():
    start = time.time()
    verbose = 0

    sysNames = ['norm']#, 'ShowerShape', 'pho_scale', 'pho_smear', 'lep_scale', 'lep_smear']
    q = multiprocessing.Queue()

    cmds = {}
    sub_p = {}

    for sysName in sysNames:

        cmds[sysName] = ['python', 'ALP_BDTSys.py', '-C', '-s'] + [sysName]
        #cmds[sysName] = 'python ALP_BDTSys.py -C ' + sysName
        print cmds[sysName]
        #sub_p[sysName] = multiprocessing.Process(target=process,args=(cmds[sysName]))

    '''
    for sysName in sysNames:
        sub_p[sysName].start()

    for sysName in sysNames:
        sub_p[sysName].join()

    results = [q.get() for j in sub_p]
    print results
    '''
    end = time.time()
    print str(round(end-start,3))+'s'




# run the submitAnalyzer() as main()
if __name__ == "__main__":
    hadd()
