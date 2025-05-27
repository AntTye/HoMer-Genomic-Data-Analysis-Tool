from flask import Flask, render_template, send_file
import pandas as pd
import io

app = Flask(__name__)

localStruct = {}
wantedSegment = ""
classification = {}

#####################################
# Contents of HoMer Leaf to Leafs
def parse_file_leaf(filepath):
   #opens and reads the file and put contents of that file into contents var
   with open(filepath, 'r') as file:
      contents = file.read()
   
   #captures each Donor: block. 
   #Donor: AspAMC34, Recipient: AveroniiBAQ135 ... -> ASPAMC34, RecRecipient: AveroniiBAQ135 ... 
   blocks = contents.strip().split("Donor: ")[1:] 

   global wantedSegment
   summaryBlock = contents.split("\n")
   wantedSegment = summaryBlock[-9:-3]
   wantedSegment.append(summaryBlock[-9])
   wantedSegment = "\n".join(wantedSegment)

   for block in blocks:
      #gets the first line when split with \n
      #splits donor and recipient into respective vars
      blockLines = block.splitlines()
      header = blockLines[0].strip()
      donor, recipient = header.split(", Recipient: ")
      
      donor = donor.strip()
      recipient = recipient.strip()

      if donor not in localStruct: #add to the local dictionary
         localStruct[donor] = {}

      localStruct[donor][recipient] = "".join(["Donor: " + block.strip()])

      if block == blocks[-1]:
         newContext = deleteSummary(localStruct, donor, recipient)
         
         localStruct[donor][recipient] = newContext

def deleteSummary(struct, donor, recipient):
    last_block = struct[donor][recipient]
    rid_summary_block = last_block.split("*")[0].strip()
    return rid_summary_block

parse_file_leaf("AeromonasHMGTs\HoMer_Leaf_to_leaf.txt")
#####################################

#####################################
# Classification of Specific Leaf Transfers

def getMultiContig(content):
   splitthat = content.split("\n")
   contigBlock = [[]]              
   currBlock = 0
   for i in range(3, len(splitthat)):
      if splitthat[i].strip() == "":
         currBlock += 1
         contigBlock.append([])
         continue
      contigBlock[currBlock].append(splitthat[i])
   return contigBlock[:-2]                          #returns amount of HMGTs for this pair

def getFamilyDonor(content):
   content = "".join(content)
   getFamily = content.split("Gene Family: ")                # get the geneFamily by splitting it at the gene family so first one is at index 1
   geneFamilyArr = []                                        #arr to store family numbers
   for i in range(1, len(getFamily), 2):                     # start at 1 and count by 2 until length of getfamily is achieved
      geneFamilyArr.append(str(getFamily[i]).split(" -")[0]) #add that to the gene family   
   return geneFamilyArr                                      # returns a array of donor genes 

def getEachContigDonor(content):
   contigBlock = [] # 4 contigs in content [[1,2,3],[],[],[]]

   for contig in content:
      currBlock = []
      for line in contig:
         currBlock += getFamilyDonor (line)
      contigBlock.append(currBlock)
   return contigBlock # get each donors from each contig block

def getSpecContigDonorRecipient(multipleBlocks): #use if want to index 0 for donor and index 1 for recipient
    specificContigBlocksD = [[]]
    specificContigBlocksR = [[]]
    currIndex = 0
    for block in multipleBlocks:
        for contig in block:
            oneItem = ""
            oneItem = oneItem.join(contig)
            atArrow = oneItem.split("-> ")
            specificContigBlocksD[currIndex].append(atArrow[0].strip())
            specificContigBlocksR[currIndex].append(atArrow[1].strip())
    
        specificContigBlocksD.append([])
        specificContigBlocksR.append([])
        currIndex += 1
    return (specificContigBlocksD[:-1], specificContigBlocksR[:-1]) 

def getSpecContig(multipleBlocks):
    specificContigBlocks = [[]]
    currIndex = 0
    for block in multipleBlocks:
        for contig in block:
            oneItem = ""
            oneItem = oneItem.join(contig)
            atArrow = oneItem.split("-> ")
            specificContigBlocks[currIndex].append(atArrow[1].strip())
    
        specificContigBlocks.append([])
        currIndex += 1
    return (specificContigBlocks[:-1])              #returns the contig line pertaining to each line in block. "[[Contig: 7, Gene Number: 127, Gene Name: 3873, Gene Family: 22661]]"

def getSpecContigDonor(multipleBlocks):
    specificContigBlocks = [[]]
    currIndex = 0
    for block in multipleBlocks:
        for contig in block:
            oneItem = ""
            oneItem = oneItem.join(contig)
            atArrow = oneItem.split("-> ")
            specificContigBlocks[currIndex].append(atArrow[0].strip())
    
        specificContigBlocks.append([])
        currIndex += 1
    return (specificContigBlocks[:-1])              #returns the contig line pertaining to each line in block. "[[Contig: 7, Gene Number: 127, Gene Name: 3873, Gene Family: 22661]]"


def getNameContigFamily(multipleBlocks):
    nameContigFamily = [[]]
    currIndex = 0
    checkIndex = 0
    for block in multipleBlocks:
        for contig in block:
            getList = contig.split(", ")
            contigNumber = getList[0][8:]
            nameNumber = getList[2][11:]
            if nameNumber[checkIndex].isalpha() or nameNumber[checkIndex] == ".":
                for i in range(len(nameNumber)):
                    nameNumber = nameNumber[i:]
                    if nameNumber[checkIndex].isdigit():
                        break
            familyNumber = getList[-1][13:]
            nameContigFamily[currIndex].append([contigNumber,nameNumber,familyNumber])
            checkIndex = 0
        nameContigFamily.append([])
        currIndex +=1

    return nameContigFamily[:-1]                    # Turns "[..[Contig: 7, Gene Number: 127, Gene Name: 3873, Gene Family: 22661]..]" to "[..[[7], [3873], [22661]]..]"

def getContentsFromFile(filepath, specificName, nameContigFamily):
    with open(filepath, "r") as f:
        content = f.read()
    splitIntoSpecificNames = content.split("\t")
    extract = nameContigFamily[0]
    getNameNumber = [[]]
    getNameNumberInt = []
    currIndex = 0
    for block in nameContigFamily:
        toInt = []
        for contigLine in block:
            toInt.append(int(contigLine[1]))
        toIntSorted = sorted(toInt)
        toStr = []
        
        for integer in toIntSorted:
            strInt = str(integer)
            toStr.append(strInt)

        getNameNumberInt.append(toIntSorted)

        getNameNumber[currIndex].append(toStr)
        getNameNumber.append([])
        currIndex += 1

    compareArray = [[]]

    currIndex = 0
    secondIndex = 0
    for block in getNameNumber[:-1]:
        for nNumber in block[0]:
            correctString = [f":{nNumber}:_contig_", f".{nNumber}:_contig_", f":peg.{nNumber}:_contig_"]
            compareArray[currIndex].append(correctString)
            secondIndex += 1
        compareArray.append([])
        currIndex += 1
        secondIndex = 0
   
    transferArray = [[]]
    indexToCompare = [[]]
    currIndex = 0
    for block in compareArray[:-1]:
        for singleTransfer in block:
            for i in range(len(splitIntoSpecificNames)):
                if any(element in splitIntoSpecificNames[i] for element in singleTransfer ):
                    transferArray[currIndex].append(splitIntoSpecificNames[i])
                    indexToCompare[currIndex].append(i)
               
        currIndex += 1
        indexToCompare.append([])
        transferArray.append([])
    startSurroundGene = [[]]
    endSurroundGene = [[]]
    currIndex = 0

    for block in indexToCompare[:-1]:
        if block[0] - 2 >= 0:
            startSurroundGene[currIndex].append(splitIntoSpecificNames[block[0] - 2])
        if block[0] - 1 >= 0:
            startSurroundGene[currIndex].append(splitIntoSpecificNames[block[0] - 1])
        if block[-1] + 1 < len(splitIntoSpecificNames):
            endSurroundGene[currIndex].append(splitIntoSpecificNames[block[-1] + 1])
        if block[-1] + 2 < len(splitIntoSpecificNames):
            endSurroundGene[currIndex].append(splitIntoSpecificNames[block[-1] + 2])
        
        startSurroundGene.append([])
        endSurroundGene.append([])
        currIndex += 1
    
    return(transferArray[:-1],startSurroundGene[:-1],endSurroundGene[:-1]) #results, starting 2 or 1 extra, ending 2 or 1 extra

def parseResults(rawResult):
    transferArray = rawResult[0]
    startSurroundGene = rawResult[1]
    endSurroundGene = rawResult [2]

    finalArr = [[]]
    finalStart = [[]]
    finalEnd = [[]]
    for i in range(len(transferArray)):
        tranItem = ""
        startItem = ""
        endItem = ""
        
        tranItem = tranItem.join(transferArray[i])
        startItem = startItem.join(startSurroundGene[i])
        endItem = endItem.join(endSurroundGene[i])

        finalArr[i].append(tranItem)
        finalStart[i].append(startItem)
        finalEnd[i].append(endItem)

        finalArr.append([])
        finalStart.append([])
        finalEnd.append([])

    numberArr = [[]]
    numberStart = [[]]
    numberEnd = [[]]
    currIndexArr = 0
    currIndexStartEnd = 0

    for i in range(len(transferArray)):
      splitFinal = [s.strip() for s in finalArr[i][0].split(":") if s.strip()]
      splitStart = [s.strip() for s in finalStart[i][0].split(":") if s.strip()]
      splitEnd = [s.strip() for s in finalEnd[i][0].split(":") if s.strip()]

      for j in range(3, len(splitFinal), 3):
         actualString = ""
         for char in splitFinal[j]:
               if char.isalpha():
                  break
               actualString += char
         numberArr[currIndexArr].append(actualString.strip())
      numberArr.append([])
      currIndexArr += 1

      def extract_number(val):
         actualString = ""
         for char in val:
               if char.isalpha():
                  break
               actualString += char
         return actualString.strip()

      collected_start = []
      for j in range(3, len(splitStart), 3):
         collected_start.append(extract_number(splitStart[j]))
         if len(collected_start) == 2:
               break
      while len(collected_start) < 2:
         collected_start.append("NA")
      numberStart[currIndexStartEnd].extend(collected_start[::-1])

      collected_end = []
      for j in range(3, len(splitEnd), 3):
         collected_end.append(extract_number(splitEnd[j]))
         if len(collected_end) == 2:
               break
      while len(collected_end) < 2:
         collected_end.append("NA")
      numberEnd[currIndexStartEnd].extend(collected_end)

      numberStart.append([])
      numberEnd.append([])
      currIndexStartEnd += 1

    return (numberArr[:-1], numberStart[:-1], numberEnd[:-1])

def classification(recipients, donors):
   results = []
   for i in range(len(recipients)):
      result = checkSymmetry(recipients[i], donors[i])
      results.append(result)
   return results     
    

def checkSymmetry(recipientFamilyList, donorFamilyList):
   if recipientFamilyList == donorFamilyList: return "In Order"
   if recipientFamilyList[::-1] == donorFamilyList: return "Reversed"
   if recipientFamilyList != donorFamilyList: return "Different"

#####################################

@app.route('/')
def base():
    total_recipients = sum(len(recipients) for recipients in localStruct.values())
    return render_template('base.html', data = localStruct, total_recipients = total_recipients, summary = wantedSegment)

@app.route('/<donor>/<recipient>')
def look_recipient(donor,recipient):
   recipientPath = f"AeromonasDataset\Genomes\{recipient}.synteny"
   donorPath = f"AeromonasDataset\Genomes\{donor}.synteny"
   contents = localStruct[donor][recipient] #entire transfer ->
   listOfContigs = getMultiContig(contents) #contig blocks ->
   listOfFamily = getEachContigDonor(listOfContigs) #array of donors for each contig

   ########################################################################
   recipientLines = getSpecContig(listOfContigs)
   nameContigFamily = getNameContigFamily(recipientLines)

   Results = getContentsFromFile(recipientPath, recipient, nameContigFamily)
   ParsedResults = parseResults(Results)

   classifications = classification(ParsedResults[0], listOfFamily)

   ########################################################################

   donorLines = getSpecContigDonor(listOfContigs)
   donorNameContigFamily = getNameContigFamily(donorLines)

   dResults = getContentsFromFile(donorPath, donor, donorNameContigFamily)
   dParsedResults = parseResults(dResults)

   zipped_results = list(zip(*ParsedResults))  # [(transfer1, start1, end1), ...]
   combined = list(zip(classifications, listOfFamily, dParsedResults[1], dParsedResults[2], zipped_results))  # [(class1, (transfer1, start1, end1)), ...]

   return render_template(
    'recipient.html',
    donor=donor,
    recipient=recipient,
    contents=contents,
    classifications=classifications,
    results=combined)


@app.route('/viewStruct')
def look_struct():
   return render_template('viewStruct.html', localStruct = localStruct)
