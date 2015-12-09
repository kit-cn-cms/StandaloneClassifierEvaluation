import ROOT
from array import array
from subprocess import call
from QueHelper import QueHelper
import time as timer
import sys
import json
import os
from subprocess import check_output
import stat

currentPath = sys.path[0]+"/"
cmsswPath=os.environ["CMSSW_BASE"]
queHelper=QueHelper("NAFSL6")

infiles=sys.argv[2:]
outdir=sys.argv[1]

maxEventsPerJob=100

listOfInputTrees=[]

if not os.path.exists(outdir):
	os.makedirs(outdir)

if not os.path.exists(currentPath+"/jobScripts"):
	os.makedirs(currentPath+"/jobScripts")
if not os.path.exists(currentPath+"/logs"):
	os.makedirs(currentPath+"/logs")

#print infiles
for f in infiles:
  listOfInputTrees.append(f)

print listOfInputTrees

mainrunline=queHelper.GetRunLines()
runlines=[]
JobIDs=[]

for job in listOfInputTrees:
  inf=ROOT.TFile(job,"READ")
  intree=inf.Get("MVATree")
  inNEvents=intree.GetEntries()
  print inNEvents
  inf.Close()

  beginEvt=0
  for ievt in range(inNEvents):
    if ((ievt+1)%maxEventsPerJob)==0 or ievt==(inNEvents-1):
      print (ievt%maxEventsPerJob), ievt
      endEvt=ievt
   
      outfname=job.rsplit("/",1)[1]
      outfname=outfname.replace(".root","_"+str(beginEvt)+"to"+str(endEvt)+".root")
    
      joblines=[]
      #print joblines
      jj=queHelper.GetExecLines()
      for jjj in jj:
        joblines.append(jjj)
 	  #print joblines
      joblines.append("cd "+currentPath+"/jobScripts\n")
      joblines.append("export FILENAME="+job+"\n")

      jobname=outfname.replace(".root","")

      joblines.append("export OUTFILENAME="+outdir+"/"+outfname+"\n")
      joblines.append("export BEGINEVENT="+str(beginEvt)+"\n")
      joblines.append("export ENDEVENT="+str(endEvt)+"\n")
      joblines.append(cmsswPath+"/test/slc6_amd64_gcc491/standaloneClassifier\n") 
      outfile=open("jobScripts/"+jobname+".sh","w")
      for line in joblines:
        outfile.write(line)
      outfile.close()
      st = os.stat("jobScripts/"+jobname+".sh")
      os.chmod("jobScripts/"+jobname+".sh", st.st_mode | stat.S_IEXEC)
	
      runlines=[]
      thisrl=mainrunline[0]
      runlines.append(thisrl)
      runlines[-1]=runlines[-1].replace("INSERTPATHHERE",currentPath)
      runlines[-1]=runlines[-1].replace("INSERTEXECSCRIPTHERE","jobScripts/"+jobname+".sh")
      #print runlines
      runfile=open("run.sh","w")
      for rl in runlines:
	runfile.write(rl)
      runfile.close()
      
      st = os.stat("run.sh")
      os.chmod("run.sh", st.st_mode | stat.S_IEXEC)
      thisID=queHelper.StartJob("./run.sh")
      print "submitted job ", thisID, job, outfname
      JobIDs.append(thisID)

      beginEvt=ievt+1
      #break

print JobIDs







