from flask import Flask, render_template

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
   splitthat = content.split("\n")                           # split lines
   contigBlock = [[]]                                        # start with one empty block
   currBlock = 0
   for i in range(3, len(splitthat)):
      if splitthat[i].strip() == "":
         currBlock += 1
         contigBlock.append([])                              # adds the new block
         continue
      contigBlock[currBlock].append(splitthat[i])
   return contigBlock[:-2]                                   #returns the contig blocks, where each contig is a elem

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

def look_specific_contig(filepath, donorFamilyList):
   with open (filepath, 'r') as f:
      content = f.read()
   lines = (content.split("\n"))
   contigLines = []
   
   frontParseFix = []

   for i in range(len(donorFamilyList)):
      currentDonor = ""
      for j in range(len(donorFamilyList[i])):
         currentDonor = ":" + donorFamilyList[i][j] + "\t"
         frontParseFix.append(currentDonor)
   
   for contig in donorFamilyList:
      
      startAcc = False
      acceptedLines = ""

      for line in lines:
         if startAcc == True:
            acceptedLines += line

         if startAcc == False:
            if any(element in line for element in frontParseFix):
               acceptedLines += line
               startAcc = True
         
         if all(element in acceptedLines for element in contig):
            break

      startAcc = False
      contigLines.append(acceptedLines)
   return (contigLines) #contigLines

def parse_file_big (listOfLines, donorFamilyList):   

   contigClassification = [] #end results
   specificDonor = 0

   for contigLine in listOfLines:
      classification = parseContigLine(contigLine,donorFamilyList[specificDonor])
      contigClassification.append(classification)
      specificDonor += 1
   return contigClassification

def parseContigLine(contigLine, donorFamilyList): #singular line, array of numbers for contig
   longestStart = ""
   frontParseFix = []
   for i in range(len(donorFamilyList)):
      currentDonor = ":" + donorFamilyList[i]
      frontParseFix.append(currentDonor)
      
   for i in range(len(donorFamilyList)):                             
      splitted = contigLine.split(donorFamilyList[i])[0] + donorFamilyList[i]       
      if len(splitted) > len(longestStart):
         longestStart = splitted
   shortestStart = longestStart

   catchCurr = []
   for i in range(len(donorFamilyList)):
      catchCurr.append(frontParseFix[i])

      if all(element in splitted for element in catchCurr):
         split_result = longestStart.split(frontParseFix[i])
         if len(split_result) > 1:
               splitted = donorFamilyList[i] + split_result[1]

      if len(shortestStart) > len(splitted):
         if all(element in splitted for element in donorFamilyList):
               shortestStart = splitted

   
   recipientFamilyList = getRecipientDonor(shortestStart) #get array of recipient genes

   result = checkSymmetry(recipientFamilyList,donorFamilyList)
   start = result[1]
   end = result[2] + 1
   fullReturn = [result[0], donorFamilyList, recipientFamilyList[start:end], recipientFamilyList]
   return (fullReturn)

def checkSymmetry(recipientFamilyList, donorFamilyList):
   positions = []
   for donor in donorFamilyList:
      index = 0
      for recipient in recipientFamilyList:
         if donor == recipient:
            positions.append(index)
         index += 1
   sortedPos = sorted(positions)
   if sortedPos == positions:
      return ["Reversed", sortedPos[0], sortedPos[-1]]
   if sortedPos == positions[::-1]:
      return ["In Order", sortedPos[0], sortedPos[-1]]
   return ["Different", sortedPos[0], sortedPos[-1]]

def getRecipientDonor (context):
   listcontext = context.split(":")
   splitAt = []

   for i in range(0, len(listcontext), 3):
      splitAt.append(listcontext[i])
   for i in range(len(splitAt)):
      stripyou = splitAt[i].split("\t")[0]
      splitAt[i] = stripyou
   if splitAt[0][0].isalpha():
      splitAt = splitAt[1:]
   return splitAt
#####################################













@app.route('/')
def base():
    total_recipients = sum(len(recipients) for recipients in localStruct.values())
    return render_template('base.html', data = localStruct, total_recipients = total_recipients, summary = wantedSegment)

@app.route('/<donor>/<recipient>')
def look_recipient(donor,recipient):
   contents = localStruct[donor][recipient] #entire transfer ->
   listOfContigs = getMultiContig(contents) #contig blocks ->
   listOfFamily = getEachContigDonor(listOfContigs) #array of donors for each contig
   
   recipientPath = f"AeromonasDataset\Genomes\{recipient}.synteny"
   listOfLines = look_specific_contig(recipientPath, listOfFamily)
   
   classification = parse_file_big(listOfLines, listOfFamily)
   return render_template('recipient.html', donor = donor, recipient = recipient, contents = contents, classification = classification)

@app.route('/viewStruct')
def look_struct():
   return render_template('viewStruct.html', localStruct = localStruct)


