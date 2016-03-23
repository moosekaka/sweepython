import re,os,glob
parDir=os.path.dirname(os.getcwd())
files = glob.glob('[0[0-5]*\[AyY]*\*meta*')

Name=[];PowerRFP=[];PowerGFP=[];Date=[];
for i in files:
    lines = open( i, "r" ).readlines()
    for line in lines:
        if re.search(r'ChNames.*BF_Confocal*',line):
            break
    else:
        Name.append(i.split('\\')[1])
        for line in lines:
            if re.search(r".561-V.",line):
                PowerRFP.append(line)     
                break

        for line in lines:    
            if re.search(r".488-V.",line):
                PowerGFP.append(line)     
                break
        
        for line in lines:
            if re.search(r".Date.", line):
                Date.append(line) 
                break
s={}

for h,i in enumerate(Name):
    temp=PowerGFP[h].split(':')[1]
    temp2=PowerRFP[h].split(':')[1]
    temp3=Date[h].split(':')[1]
    s[i]=(temp2.strip()[1:-2],temp3.strip()[1:-2],temp.strip()[1:-2])
    
    
    
    
