#!/usr/bin/python
#------------------------------------------------
# Script runs several MEKD producers as threads
# Latest update: 2012.11.04
# More info: predrag.milenovic@cern.ch
#------------------------------------------------
import sys, os, commands, threading
import optparse

# define function for argument/options parsing
def parseOptions():

    usage = ('usage: %prog [options] rootFilesDirectory\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    theRootDir=""
    theRootTree="passedEvents"
    theSelection=""
    theppEnergy="8"
  
    parser.add_option('-d', '--dir',     dest='rootDir',    type='string', default=theRootDir,     help='The name of the directory with input tree.')
    parser.add_option('-t', '--tree',    dest='rootTree',   type='string', default=theRootTree,    help='The name of the input tree.')
    parser.add_option('-e', '--energy',  dest='energy',     type='string', default=theppEnergy,    help='pp collsiion energy in TeV.')
    parser.add_option('-s', '--select',  dest='select',     type='string', default=theSelection,   help='The string that must be in the name of the input file.')

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

    # sanity check for some of the arguments
    if len(args)!=1:
        parser.error('Wrong number of arguments. Exiting...')
    if (not os.path.exists(args[0])):
        parser.error('File directory ['+args[0]+'] does not exist. Exiting...')

    if (not opt.rootDir):
        print 'No root dir specified. Using the default one.'

    if (not opt.rootTree):
        print 'No root tree specified. Using the default one.'

#define function for processing of os command
def processCmd(cmd):
    print cmd
    status, output = commands.getstatusoutput(cmd)
    if status !=0:
        print 'Error in processing command:\n   ['+cmd+'] \nExiting...'
        sys.exit()

#define function for running os command via threading
def runCommand(cmd):
    semaphore.acquire()
    try:
        os.system(cmd)
    finally:
        semaphore.release()

# the main procedure
def AddKDs_root():
    
    # parse the arguments and options
    global opt, args, semaphore
    parseOptions()

    # build root dir/tree
    if (opt.rootDir==""):
        rootDirTree = opt.rootTree
    else:
        rootDirTree = opt.rootDir + '/' + opt.rootTree

    # print main config info
    print '--------------------------------------------------'
    print '--------------------------------------------------'
    print 'Directory to process:  ' + args[0];
    print 'Root dir/tree:         ' + rootDirTree
    print 'File selection mask:   ' + opt.select
    print '--------------------------------------------------'
    print '--------------------------------------------------'
    
    # configure threading
    max_processes = 32
    semaphore = threading.Semaphore(max_processes)
    threads = []
    
    # go through directory and get .root files
    for filename in os.listdir(args[0]):
        if ((opt.select in filename) and ('.root' in filename)):
            # build command
            command = './runKD_MAD -f "' + args[0] + '/' + filename + '" -t "' + rootDirTree + '"' + ' -s ' + opt.energy
            # run as a thread
            t = threading.Thread(target=runCommand, args=(command, ))
            t.start()
            threads.append(t)
        else:
            continue

    # Wait for all of them to finish
    [x.join() for x in threads]

    sys.exit()

# run the AddKDs_root() as main()
if __name__ == "__main__":
    AddKDs_root()
