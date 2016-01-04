import sys

infile=sys.argv[1]

inf=open(infile,"r")
inlist=list(inf)

jetstrings=[]
lepstrings=[]
metstrings=[]
csvString="{"
leptoncharges=[]


for l in inlist:
  thisline=l.split()
  print thisline
  if "JET"==thisline[0]:
    jstr="p4("
    for w in thisline[1:-1]:
      jstr+=w+","
    jstr+="),"
    csvString+=thisline[-1]+","
    jetstrings.append(jstr)
 
  if "LEP"==thisline[0]:
    lstr="p4("
    for w in thisline[1:-1]:
      lstr+=w+","
    lstr+="),"
    leptoncharges.append(thisline[-1])
    lepstrings.append(lstr)

  if "MET"==thisline[0]:
    mstr="lv_met.SetPtEtaPhiE(("
    for w in thisline[1:-1]:
      mstr+=w+","
    mstr+="),"


csvString+="}"

for i in jetstrings:
  print i
print csvString
for i in lepstrings:
  print i
print leptoncharges

for i in metstrings:
  print i

