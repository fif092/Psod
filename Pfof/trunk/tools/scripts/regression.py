#!/usr/bin/env python3
# Regression test for the pfof package
import argparse
import os
import shutil
import fileinput
import subprocess

# colors definitions
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

# paths definitions
homepath=os.environ.get("HOME")
trunkpath=homepath+"/Travail/Devel/Cosmologie/Pfof/trunk"
releasepath=homepath+"/Travail/Devel/Cosmologie/Pfof/tags/Release-1.0"
datapath=homepath+"/Travail/Devel/Cosmologie/test/Data/output_00013"
testpath=homepath+"/Travail/Devel/Cosmologie/test"


# edit a .nml file:
# give the value parametervalue to the parameter parametername 
def editnamelist( filename, parametername, parametervalue ):
    "Edit a pfof namelist .nml file"
    
    for line in fileinput.input(filename,inplace=True):
        info=line.split("=")
        if info[0].strip() == parametername:
            print(line.replace(info[1].strip(), parametervalue),end="")
        else:
            print(line,end="")
    return ;

# list all the files of a certain type (cube/mpicube/halo/halomass)
# produced by the trunk and the release version of pfof
def listfiles(filetype) :
    shellcmd="ls trunk_"+filetype+"*.h5"
    buffershell = subprocess.check_output(shellcmd, shell=True)
    filestrunk = buffershell.decode(encoding='UTF-8').split("\n")
    shellcmd="ls release_"+filetype+"*.h5"
    buffershell = subprocess.check_output(shellcmd, shell=True)
    filesrelease = buffershell.decode(encoding='UTF-8').split("\n")

    return filestrunk, filesrelease;

# print a succes message if HDF5 files match
def h5diffpassed(filetype):
    print(filetype+" files:"+bcolors.OKGREEN+"[PASSED]"+bcolors.ENDC)
    return;

# print a failure message if HDF5 file don't match
def h5difffailed(filetype,message):
    print(message)
    print(filetype+" files:"+bcolors.FAIL+"[FAILED]"+bcolors.ENDC)
    return ;

# parse h5diff output: if the files don't match, return a non empty error string
def parseh5diff(filestrunk,filesrelease) :
    for i in range(0,len(filestrunk)):
        if filestrunk[i]:
            shellcmd="h5diff "+filestrunk[i]+" "+filesrelease[i]
            p = subprocess.Popen(shellcmd, shell=True, stdin=subprocess.PIPE, 
                                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
            buffershell = p.stdout.read()
            lines = buffershell.decode(encoding='UTF-8').split("\n")
            error=""
            for i in range(0,len(lines)):
                if lines[i].find("dataset")!= -1 :
                    error="dataset"
                elif (lines[i].find("attribute")!= -1 and
                (lines[i].find("date") == -1 and 
                 lines[i].find("time") == -1 )):
                    error="attribute"
                elif lines[i].find("Some objects are not comparable") != -1:
                    error="objects"
    return error;


# compare HDF5 files of a certain type written by the trunk and release versions of pfof  
def comparehdf5( filetype ):

    filestrunk,filesrelease = listfiles(filetype)
    if len(filestrunk) != len(filesrelease):
        h5difffailed(filetype,"Trunk and release versions did not write the same number of "+filetype+" files.")
    else:
        error=parseh5diff(filestrunk,filesrelease)
    
    if not error:
        h5diffpassed(filetype)
    else:
        h5difffailed(filetype,"Some "+error+" differ in "+filetype+" files.")

    return;

# Regression test for pfofhdf5
def testpfofhdf5( ):
    "Regression test for pfofhdf5."
    print("Regression test for pfofhdf5.")
    print("This test in not fully implemented.")

    cwd = os.getcwd() # get current directory

    build_dir = trunkpath+"/pfof_hdf5/src"
    try:
        os.chdir(build_dir)
        os.system("make")
        trunkexe=testpath+"/pfofhdf5trunk"
        shutil.copy2("pfofhdf5",trunkexe)
        shutil.copy2("../pfof.nml",testpath+"/pfoftrunk.nml")
    finally:
        os.chdir(cwd)

    build_dir = releasepath+"/pfof_hdf5/src"
    try:
        os.chdir(build_dir)
        os.system("make")
        releaseexe=testpath+"/pfofhdf5release"
        shutil.copy2("pfofhdf5",releaseexe)
        shutil.copy2("../pfof.nml",testpath+"/pfofrelease.nml")
    finally:
        os.chdir(cwd)

    print("1st test: read from Ramses, HDF5 output, b=0.2, Mmin=100")
    editnamelist( "pfoftrunk.nml","output_root","'trunk'" )
    editnamelist( "pfoftrunk.nml","pathinput","'"+datapath+"'")
    shutil.copy2("pfoftrunk.nml","pfof.nml")
    os.system("mpirun -np 8 "+trunkexe+" > trunk1.log")

    editnamelist( "pfofrelease.nml","root","'release'" )
    editnamelist( "pfofrelease.nml","pathinput","'"+datapath+"'")
    shutil.copy2("pfofrelease.nml","pfof.nml")
    os.system("mpirun -np 8 "+releaseexe+" > release1.log")

    comparehdf5("cube")
    comparehdf5("halo")
    comparehdf5("halomass")

    print("2nd test: read from HDF5 cube files, HDF5 output, mpicube output, b=0.2, Mmin=100")
    
    editnamelist("pfoftrunk.nml", "gatherwrite","2")
    editnamelist("pfoftrunk.nml", "readfromcube",".true.")
    shutil.copy2("pfoftrunk.nml", "pfof.nml")
    os.system("mpirun -np 8 "+trunkexe+" > trunk2.log")

    editnamelist("pfofrelease.nml", "gatherwrite","2")
    editnamelist("pfofrelease.nml", "readfromcube",".true.")
    shutil.copy2("pfofrelease.nml", "pfof.nml")
    os.system("mpirun -np 8 "+releaseexe+" > release2.log")

    comparehdf5("mpicube")
    comparehdf5("halo")
    comparehdf5("halomass")

    print("3rd test: read from HDF5 mpicube files, HDF5 output, cube output, b=0.3, Mmin=200")
    
    editnamelist("pfoftrunk.nml", "gatherwrite","1")
    editnamelist("pfoftrunk.nml", "gatherread","2")
    editnamelist("pfoftrunk.nml", "perco","0.30")
    editnamelist("pfoftrunk.nml", "Mmin","200")
    editnamelist("pfoftrunk.nml", "readfromcube",".true.")
    shutil.copy2("pfoftrunk.nml", "pfof.nml")
    os.system("mpirun -np 8 "+trunkexe+" > trunk3.log")

    editnamelist("pfofrelease.nml", "gatherwrite","1")
    editnamelist("pfofrelease.nml", "gatherread","2")
    editnamelist("pfofrelease.nml", "perco","0.30")
    editnamelist("pfofrelease.nml", "Mmin","200")
    editnamelist("pfofrelease.nml", "readfromcube",".true.")
    shutil.copy2("pfofrelease.nml", "pfof.nml")
    os.system("mpirun -np 8 "+releaseexe+" > release3.log")

    comparehdf5("cube")
    comparehdf5("halo")
    comparehdf5("halomass")

    return;


def testpfofcone( ):
    "Regression test for pfofcone."
    print("Regression test for pfofcone.")
    print("This test in not fully implemented.")
    return;

def testconecreator( ):
    "Regression test for conecreator."
    print("Regression test for conecreator.")
    print("This test in not fully implemented.")
    return;

def testconemapper( ):
    "Regression test for conemapper."
    print("Regression test for conemapper.")
    print("This test in not fully implemented.")
    return;

def testhaloanalyzer( ):
    "Regression test for haloanalyzer."
    print("Regression test for haloanalyzer.")
    print("This test in not fully implemented.")
    return;



# options: name of the code 
codelist = ["all","pfofcone","pfofhdf5","conecreator",
                             "conemapper","haloanalyzer"]

parser = argparse.ArgumentParser(description="Regression test for the pfof package",
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument("-c","--code",
                    default="all",
                    choices=codelist,
                    help="Code to test. Default value is %(default)s. Allowed values are "+
                    ", ".join(codelist), metavar='')

args = parser.parse_args()

print ("Regression test for the pfof package.")
print (" ")

if args.code == "pfofhdf5":
    testpfofhdf5();
elif args.code == "pfofcone":
    testpfofcone();
elif args.code == "conecreator":
    testconecreator();
elif args.code == "conemapper":
    testconemapper();
elif args.code == "haloanalyzer":
    testhaloanalyzer();
elif args.code == "all":
    testpfofhdf5();
    testpfofcone();
    testconecreator();
    testconemapper();
    testhaloanalyzer();

