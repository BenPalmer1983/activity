# !/usr/bin/python
import os, sys, time


#---------------------------------------------------------------------------------------------
# USEFUL FUNCTIONS 11th December 2014
#---------------------------------------------------------------------------------------------

def subString(inputStr, start, length):
  a = start-1
  output = inputStr[a:a+length]
  return output   

def replace(haystack,needle,replace):  # Replace, case insensitive
  output = ""
  i=0
  while i<len(haystack):
    i = i + 1
    testNeedle = subString(haystack,i,len(needle))
    if needle==testNeedle:
      i = i + len(needle)-1    
      output = output + replace
    else:
      output = output + subString(haystack,i,1)
  return output
  
def stristr(haystack,needle):  # Replace, case insensitive
  output = 0
  i=0
  while i<len(haystack):
    i = i + 1
    testNeedle = subString(haystack,i,len(needle))
    if needle==testNeedle:  
      output = 1
      i = len(haystack) 
  return output

def replaceCI(haystack,needle,replace):  # Replace, case insensitive
  haystackU = haystack.upper()
  needleU = needle.upper()
  output = ""
  i=0
  while i<len(haystackU):
    i = i + 1
    testNeedleU = subString(haystackU,i,len(needle))
    if needleU==testNeedleU:
      i = i + len(needle)-1    
      output = output + replace
    else:
      output = output + subString(haystack,i,1)
  return output
  
def removeSpaces(inputStr):
# Removes spaces
  output = ""
  i=0
  while i<len(inputStr):
    i = i + 1
    testChar = subString(inputStr,i,1)
    if ord(testChar) != 32:
      output = output + testChar
  return output  
  
def removeBlanks(inputStr):
# Removes spaces, tabs, returns
  output = ""
  i=0
  while i<len(inputStr):
    i = i + 1
    testChar = subString(inputStr,i,1)
    if ord(testChar) not in (9,10,13,32):
      output = output + testChar
  return output  
  
  
def removeTabs(inputStr):
# Removes tabs
  output = ""
  i=0
  while i<len(inputStr):
    i = i + 1
    testChar = subString(inputStr,i,1)
    if ord(testChar) != 9:
      output = output + testChar
  return output    
  

def removeCR(inputStr):
# Removes spaces
  output = ""
  i=0
  while i<len(inputStr):
    i = i + 1
    testChar = subString(inputStr,i,1)
    if ord(testChar) not in (10,13):
      output = output + testChar
  return output      
  
def trim(inputStr):
  output = ""
  i = len(inputStr)
  lastSpace = 0
  while i>0:
    charSelected = subString(inputStr,i,1)
    if lastSpace==0:
      if charSelected!=" ":
        lastSpace = 1
        output = charSelected + output
    else:    
      output = charSelected + output
    i = i - 1
  return output
  
def trimAll(inputStr):
  output = ""
  i = len(inputStr)
  lastSpace = 0
  while i>0:
    charSelected = subString(inputStr,i,1)
    if lastSpace==0:
      if ord(charSelected) not in (9,10,13,32):
        lastSpace = 1
        output = charSelected + output
    else:    
      output = charSelected + output
    i = i - 1
  return output  
  
def trimLeading(inputStr):
  output = ""
  i = 0
  iEnd = len(inputStr)
  firstSpace = 1
  while i<iEnd:
    i = i + 1
    charSelected = subString(inputStr,i,1)
    if firstSpace==1:
      if ord(charSelected)!=(32):
        firstSpace = 0
        output = output + charSelected
    else:    
      output = output + charSelected
  return output  

def countLeading(inputStr):
  i = 0
  iEnd = len(inputStr)
  firstSpace = 1
  counter = 0
  while i<iEnd:
    i = i + 1
    charSelected = subString(inputStr,i,1)
    if firstSpace==1:
      if ord(charSelected)!=(32):
        firstSpace = 0
      else:    
        counter = counter + 1
  return counter  
  
def spacesString(length):
  i = 0
  output = ""
  while i<length:  
    i = i + 1
    output = output + " "
  return output      
  
def removeTrailingComma(inputStr):    
  inputStr = trimAll(inputStr)
  if subString(inputStr,len(inputStr),1)==",":
    inputStr = subString(inputStr,1,len(inputStr)-1)
  return inputStr
  
def floor(input):
  output = int(input)
  return output
  
def ceil(input):
  inInt = int(input)
  inIntFloat = float(inInt)  
  if input==inIntFloat:
    output = inInt
  else:  
    output = inInt + 1
  return output  
    
def strLen(input):    
  output = len(input)  
  
def convertString(input):
  inputStripped = removeBlanks(input)
  if(len(inputStripped)==0):
    type = "string"
  else:
    i=1
    type = "int"
    dpCount = 0
    while i<=len(inputStripped):
      if ord(subString(inputStripped,i,1)) not in (46,48,49,50,51,52,53,54,55,56,57):
        type = "string"
        i = len(inputStripped)+1
      if subString(inputStripped,i,1)==".":
        dpCount = dpCount + 1
      i = i + 1  
    if dpCount==1 and type!="string":
      type = "float"
    if dpCount>1:
      type = "string"
  if type=="string":
    output = input
  if type=="float":
    output = float(inputStripped)
  if type=="int":
    output = int(inputStripped)
  return output  
    
def stringType(input):
  inputStripped = removeBlanks(input)
  if(len(inputStripped)==0):
    type = "string"
  else:
    i=1
    type = "int"
    dpCount = 0
    while i<=len(inputStripped):
      if ord(subString(inputStripped,i,1)) not in (43,45,46,48,49,50,51,52,53,54,55,56,57,101):
        type = "string"
        i = len(inputStripped)+1
      if subString(inputStripped,i,1)==".":
        dpCount = dpCount + 1
      i = i + 1  
    if dpCount==1 and type!="string":
      type = "float"
    if dpCount>1:
      type = "string"
  return type    
    
  
#---------------------------------------------------------------------------------------------
# End of useful functions section
#---------------------------------------------------------------------------------------------


def readDataLine(inputStr, format):    
  dataArray=[None]*10  
  dataArrayOut=[None]*10
  dataArray[1] = trimLeading(trim(subString(inputStr,1,11)))
  dataArray[2] = trimLeading(trim(subString(inputStr,12,11)))
  dataArray[3] = trimLeading(trim(subString(inputStr,23,11)))
  dataArray[4] = trimLeading(trim(subString(inputStr,34,11)))
  dataArray[5] = trimLeading(trim(subString(inputStr,45,11)))
  dataArray[6] = trimLeading(trim(subString(inputStr,56,11)))
  if format==1:
    for i in range(1,7):
      if stringType(dataArray[i])=="string":
        dataArrayOut[i] = 0
      else:
        factor = 1
        if subString(dataArray[i],1,1)=="-":
          dataArray[i] = subString(dataArray[i],2,len(dataArray[i]))
          factor=-1
        dataArray[i] = replace(dataArray[i],"-","e-")
        dataArray[i] = replace(dataArray[i],"+","e")
        dataArrayOut[i] = factor*float(dataArray[i])
  else:
    for i in range(1,7):  
      dataArrayOut[i] = dataArray[i]
  return dataArrayOut

dataArray=[None]*10
xVal=[None]*1000
yVal=[None]*1000
xValA=[None]*1000
yValA=[None]*1000
directory = "/data/TENDL/TENDL-n"
#directory = "/data/TENDL/TENDL-p"
#directory = "/data/TENDL/test"
for file in os.listdir(directory):
  if file.endswith(".tendl"):  
    tendlFile = open(directory+"/"+file, 'r')
    print "Processing: "+directory+"/"+file
# set projectile    
    if stristr(file,"-p")==1:
      projZ = 1  #Projectile proton number
      projA = 1  #Projectile nucleon number
    readData = 0   
    if stristr(file,"-n")==1:
      projZ = 0  #Projectile proton number
      projA = 1  #Projectile nucleon number   
    readData = 0     
#-----------------------
# Read PROTON file
    if projZ==1 and projA==1:  
      rowCounter = 0
      for line in tendlFile:      
        rowCounter = rowCounter + 1      
        newLine = line            
        newLine = replace(newLine,chr(10),"")
  # read header
        if readData==0:
  # save line markers
          mf = int(subString(newLine,71,2))
          mt = int(subString(newLine,73,3))
          row = int(subString(newLine,76,5))
  # target details      
          if mf==1 and mt==451 and row==1:
            dataArray = readDataLine(newLine,1)
            #print dataArray[1]
            target = dataArray[1]
            tarZ = floor(target/1000.0)
            tarA = int(target-1000*tarZ)
            tarM = 0
            if stristr(file,"m-")==1:
              tarM = 1
  # output file name        
            outFile = str(projZ)+"_"+str(projA)+"-"+str(tarZ)+"_"+str(tarA)+"_"+str(tarM)+".dat"
            #print "Output file: "+outFile    
            outputFile = open(outFile, 'w')
  # read data header
          if mf==3 and mt==5 and row==3:
            dataArray = readDataLine(newLine,0)
            readData = 1
            xsPoints = int(dataArray[1])
            #print "Data points to read: ",xsPoints
            xsCounter = -1
  # read data header
          if mf==6 and mt==5:
            dataArray = readDataLine(newLine,1)
            if int(dataArray[3])==0 and int(dataArray[4])==1 and int(dataArray[5])==1:
              if dataArray[1]>0.0:
                #print "Product: "+str(dataArray[1])              
                product = int(1*float(dataArray[1]))
                productZ = floor(product/1000.0)
                productA = int(product-1000*productZ)
                productM = 0
                readData = 2
                ypLine = 0
            
  # read in cross section data    
        if readData==1:
          if xsCounter==-1:
            xsCounter = xsCounter + 1  # Skip this line
          else:  
            dataArray = readDataLine(newLine,1)
            for i in range(1,4):    # 1 2 3
              xsCounter = xsCounter + 1 
              if xsCounter<=xsPoints:
                #print i,xsCounter,str(dataArray[2*i-1])+","+str(dataArray[2*i])
                xVal[xsCounter] = float(dataArray[2*i-1])  # Energy
                yVal[xsCounter] = float(dataArray[2*i])    # XS
              if xsCounter==xsPoints: 
                readData = 0  # break out

  # read in number of points for product yield
        if readData==2:
          if ypLine==1:
            dataArray = readDataLine(newLine,1)
            pyPoints = int(dataArray[1])
            pyCounter = 0
            readData = 3
            #print "Yield Points: ",pyPoints
            ypLine = 0
          ypLine = ypLine + 1


  # read in points for product yield
        if readData==3:
          if ypLine==1:
            outputLine = "#Header "
            outputFile.write(outputLine+"\n")
            outputLine = "#Target "+str(tarZ)+" "+str(tarA)+" "+str(tarM)
            outputFile.write(outputLine+"\n")
            outputLine = "#Product "+str(productZ)+" "+str(productA)+" "+str(productM)
            outputFile.write(outputLine+"\n")
            if pyPoints>xsPoints:
              outputLine = "#Datapoints "+str(xsPoints)
              outputFile.write(outputLine+"\n")
            else:
              outputLine = "#Datapoints "+str(pyPoints)
              outputFile.write(outputLine+"\n")
          if ypLine>1:
            dataArray = readDataLine(newLine,1)
            for i in range(1,4):    # 1 2 3
              pyCounter = pyCounter + 1
              if pyCounter<=pyPoints:      
                xValA[pyCounter] = float(dataArray[2*i-1])   #Energy
                yValA[pyCounter] = float(dataArray[2*i])     #XS
              if pyCounter==pyPoints:
                readData = 0 
                for j in range(1,xsPoints+1):
                  energyVal = xVal[j]
                  xsVal = 0.0
                  for k in range(1,pyPoints+1):
                    if energyVal==xValA[k]:
                      xsVal = yVal[j]*yValA[k]     
                  outputLine = str(energyVal)+" "+str(xsVal) 
                  outputFile.write(outputLine+"\n")                     
          ypLine = ypLine + 1

#-----------------------
# Read NEUTRON file
    if projZ==0 and projA==1:  
      rowCounter = 0
      for line in tendlFile:   
        rowCounter = rowCounter + 1 
        newLine = line            
        newLine = replace(newLine,chr(10),"")
  # read header
        if readData==0:
  # save line markers
          mf = int(subString(newLine,71,2))
          mt = int(subString(newLine,73,3))
          row = int(subString(newLine,76,5))
  # target details      
          if mf==1 and mt==451 and row==1:
            dataArray = readDataLine(newLine,1)
            #print dataArray[1]
            target = dataArray[1]
            tarZ = floor(target/1000.0)
            tarA = int(target-1000*tarZ)
            tarM = 0
            if stristr(file,"m-")==1:
              tarM = 1
  # output file name        
            outFile = str(projZ)+"_"+str(projA)+"-"+str(tarZ)+"_"+str(tarA)+"_"+str(tarM)+".dat"
            #print "Output file: "+outFile    
            outputFile = open(outFile, 'w')
          writeToFile = 0  
          if mf==3 and mt==16 and row==1:   # MT16 (n,2n)  
            writeToFile = 1  
            readData = 1
            lineCount = 0
            productZ = tarZ
            productA = tarA-1
            productM = 0
          if mf==3 and mt==17 and row==1:   # MT17 (n,3n)  
            writeToFile = 1  
            readData = 1
            lineCount = 0
            productZ = tarZ
            productA = tarA-2
            productM = 0
          if mf==3 and mt==22 and row==1:   # MT22 (n,na)  
            writeToFile = 1  
            readData = 1
            lineCount = 0
            productZ = tarZ-2
            productA = tarA-4
            productM = 0
          if mf==3 and mt==28 and row==1:   # MT28 (n,np)  
            writeToFile = 1  
            readData = 1
            lineCount = 0
            productZ = tarZ-1
            productA = tarA-1
            productM = 0
          if mf==3 and mt==103 and row==1:   # MT103 (n,p)  
            writeToFile = 1  
            readData = 1
            lineCount = 0
            productZ = tarZ-1
            productA = tarA
            productM = 0
          if mf==3 and mt==104 and row==1:   # MT104 (n,d)  
            writeToFile = 1  
            readData = 1
            lineCount = 0
            productZ = tarZ-1
            productA = tarA-1
            productM = 0
          if mf==3 and mt==105 and row==1:   # MT105 (n,t) 
            writeToFile = 1   
            readData = 1
            lineCount = 0
            productZ = tarZ-1
            productA = tarA-2
            productM = 0
          if mf==3 and mt==106 and row==1:   # MT106 (n,h)  
            writeToFile = 1  
            readData = 1
            lineCount = 0
            productZ = tarZ-2
            productA = tarA-2
            productM = 0
          if mf==3 and mt==107 and row==1:   # MT107 (n,a)  
            writeToFile = 1  
            readData = 1
            lineCount = 0
            productZ = tarZ-2
            productA = tarA-3
            productM = 0
# Save header file         
          if writeToFile==1:
            outputLine = "#Header "
            outputFile.write(outputLine+"\n")
            outputLine = "#Target "+str(tarZ)+" "+str(tarA)+" "+str(tarM)
            outputFile.write(outputLine+"\n")
            outputLine = "#Product "+str(productZ)+" "+str(productA)+" "+str(productM)
            outputFile.write(outputLine+"\n")
  # read header
        if readData==1:
          lineCount = lineCount+1
          if lineCount==3:  # data points count
            dataArray = readDataLine(newLine,1)
            xsPoints = int(dataArray[1])
            xsCounter = 0
            outputLine = "#Datapoints "+str(xsPoints)
            outputFile.write(outputLine+"\n")
          if lineCount>=4:  # data points count
            dataArray = readDataLine(newLine,1)
            for i in range(1,4):    # 1 2 3
              xsCounter = xsCounter + 1
              if xsCounter<=xsPoints:      
                energyVal = float(dataArray[2*i-1])   #Energy
                xsVal = float(dataArray[2*i])     #XS
                outputLine = str(energyVal)+" "+str(xsVal) 
                outputFile.write(outputLine+"\n")  
              if xsCounter==xsPoints:
                readData = 0
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          



