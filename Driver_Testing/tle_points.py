#Written by Ayush Rautwar
import re
import time
from datetime import datetime
import math

def checksum(line):
	hey = ''.join(filter(lambda x: x.isdigit() or x == "-", line)).replace(" ","")
	hey2 = list(hey)
	b = 0
	for a in hey2:
		if(a=="-"):
			b=b+1
		elif(int(a)<10):
			b=b+int(a)
	return b%10

def propagate():
	GM = 3.986004418*(10**14)
	file = open("tjreverbtle.txt", "r")
	lines = file.readlines()
	eachline = list()
	for line in lines:
		eachline.append(line.split())
	d = datetime.utcnow()
	poskep = [1.2,0.007129387,1.2,1.2,1.5708534575885496,6.283131851920631]
	print(eachline[1][3])
	edays = str(round((d-datetime(2019, 1, 1, 0)).total_seconds()/(24*60*60)+1, 8))
	parts = edays.split(".")
	parts[0] = parts[0].rjust(3, '0')
	parts[1] = '{:.8f}'.format(float("."+parts[1]))
	print(parts[0]+"."+parts[1].split(".")[1])
	eachline[1][3] = str(d.strftime("%y")) + parts[0]+"."+parts[1].split(".")[1]
	print(eachline[1][3])

	i = str(round(poskep[2],4)).split(".")
	eachline[2][2] = i[0].rjust(3, " ") + "." + i[1].ljust(4, '0')
	print(eachline[2][2])

	raan = str(round(poskep[3],4)).split(".")
	eachline[2][3] = raan[0].rjust(3, " ") + "." + raan[1].ljust(4, '0')
	print(eachline[2][3])

	e = str(round(poskep[1],7)).split(".")
	eachline[2][4] = e[1].rjust(7, '0')
	print(eachline[2][4])

	argp = str(round(poskep[4],4)).split(".")
	eachline[2][5] = argp[0].rjust(3, " ") + "." + argp[1].ljust(4, '0')
	print(eachline[2][5])


	alt = 400 #will be gps altitude
	meanmot = (GM/((alt+6378000)**3))**(1/2)/(2*math.pi)*(24*60*60) # gps alt
	meanmot = str(round(meanmot,8)).split(".")
	meanmot = meanmot[0].rjust(2, " ") + "." + str(round(int(meanmot[1]),8)).ljust(8, '0')
	m1 = meanmot
	meanmot = float(meanmot)


	yamlmeananom = 0        #test values
	yamllastyear = 2019     #test values
	yamllastday = 43        #test values


	if yamllastyear%4==0 and yamllastyear%100 !=0:
		days = 366
	else:
		days = 365

	meananom = yamlmeananom + meanmot*(days*abs(yamllastyear-datetime.utcnow().year)+(float(parts[0]+"."+parts[1].split(".")[1])-yamllastday)) #meananomaly = oldmeananomaly + meanmotion * dayssincelasttime
	eachline[2][6] = str(meananom).lstrip('0').split(".")[0].rjust(3, " ") + "." + str(round(meananom, 4)).split(".")[1].ljust(4, '0')
	#write new mean anomaly, last year, and last day to config yaml file
	print(eachline[2][6])

	yamlmeanmot = 15        #test value

	firstd = (meanmot - yamlmeanmot)/((days*abs(yamllastyear-datetime.utcnow().year)+(float(parts[0]+"."+parts[1].split(".")[1])-yamllastday)))
	firstd = firstd/2 #apparently tle is half of this
	#write meanmot to yaml
	randomsign = ""
	if(str(firstd)=="-"):
		randomsign = "-"
	eachline[1][4] = randomsign+str(firstd).lstrip('-').lstrip('0').split(".")[0].rjust(1, ' ')+"."+str(round(firstd,8)).split(".")[1].ljust(8, '0')
	print(eachline[1][4])

	eachline[1][5] = " 00000-0"
	print (eachline[1][5])

	eachline[1][7] = "0"
	print (eachline[1][7])

	yamlrevnum = "1" #test value

	eachline[2][7]=m1+yamlrevnum.rjust(5, '0')
	#every time mean anomaly is 0 update yamlrevnum+1
	print(eachline[2][7])
	lines = [eachline[1], eachline[2]]

	c1 = " 999"+str(checksum(lines[0][:-1]))
	eachline[1][8]= c1
	print(eachline[1][8])

	c2 = str(checksum(lines[1]))
	eachline[2][7]=eachline[2][7]+c2
	print(eachline[2][7])
	eachline[1][6]=" "+eachline[1][6]
	#print(eachline[2])
	print(eachline)
	eachline[1][2] = eachline[1][2]+"  "
	lines = [eachline[1], eachline[2]]

	for q in lines:
			for p in q:
					print(p, end = " ")
			print()


