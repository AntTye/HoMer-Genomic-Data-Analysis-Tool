content = "Donor: AveroniiAK241, Recipient: AveroniiAMC35\nNumber of nodes between AveroniiAK241 and AveroniiAMC35: 6.0\n \nContig: 6, Gene Number: 352, Gene Name: 3206, Gene Family: 15002 -----\u003e Contig: 1, Gene Number: 3133, Gene Name: 1445, Gene Family: 15002\nContig: 6, Gene Number: 353, Gene Name: 3207, Gene Family: 17061 -----\u003e Contig: 1, Gene Number: 2914, Gene Name: 1216, Gene Family: 17061\nContig: 6, Gene Number: 355, Gene Name: 3210, Gene Family: 2630 -----\u003e Contig: 1, Gene Number: 2916, Gene Name: 1219, Gene Family: 2630\n \nNumber of transfers in HMGT 1 = 3\n \nTotal number of HMGTs = 1\nTotal number of genes in HMGTs for this pair = 3\nTotal number of HGTs for this pair = 68\n_______________________________________________________________________"
donorPath = "AeromonasDataset\Genomes\AveroniiAER397.synteny"
recipientPath = "AeromonasDataset\Genomes\AveroniiAMC35.synteny"



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

def getNameContigFamily(multipleBlocks):
    nameContigFamily = [[]]
    currIndex = 0
    for block in multipleBlocks:
        for contig in block:
            getList = contig.split(", ")
            contigNumber = getList[0][8:]
            nameNumber = getList[2][11:]
            familyNumber = getList[-1][13:]
            nameContigFamily[currIndex].append([contigNumber,nameNumber,familyNumber])
        nameContigFamily.append([])
        currIndex +=1
        print(nameContigFamily[:-1])
    return nameContigFamily[:-1]                    # Turns "[..[Contig: 7, Gene Number: 127, Gene Name: 3873, Gene Family: 22661]..]" to "[..[[7], [3873], [22661]]..]"

def getDonorNameContigFamily(multipleBlocks):
    nameContigFamily = [[]]
    currIndex = 0
    for block in multipleBlocks:
        for contig in block:
            getList = contig.split(", ")
            contigNumber = getList[0][8:]
            nameNumber = getList[2][11:]
            familyNumber = getList[-1][13:]
            nameContigFamily[currIndex].append([contigNumber,nameNumber,familyNumber])
        nameContigFamily.append([])
        currIndex +=1
        print(nameContigFamily[:-1])
    return nameContigFamily[:-1]                    # Turns "[..[Contig: 7, Gene Number: 127, Gene Name: 3873, Gene Family: 22661]..]" to "[..[[7], [3873], [22661]]..]"


def getContentsFromFile(filepath, specificName, nameContigFamily):
    with open(filepath, "r") as f:
        content = f.read()
    splitIntoSpecificNames = content.split("\t")

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
    for block in getNameNumber[:-1]:
        for nNumber in block[0]:
            correctString = f"{specificName}:{nNumber}:_contig_"
            compareArray[currIndex].append(correctString)
        compareArray.append([])
        currIndex += 1

    transferArray = [[]]
    indexToCompare = [[]]
    currIndex = 0
    for block in compareArray[:-1]:
        for singleTransfer in block:
            for i in range(len(splitIntoSpecificNames)):
                if len(transferArray[currIndex]) == len(block):
                    break 
                if singleTransfer in splitIntoSpecificNames[i]:
                    transferArray[currIndex].append(splitIntoSpecificNames[i])
                    indexToCompare[currIndex].append(i)
        currIndex += 1
        indexToCompare.append([])
        transferArray.append([])
    
    startSurroundGene = [[]]
    endSurroundGene = [[]]
    currIndex = 0

    for block in indexToCompare[::-1]:  # safer backward iteration
        if 0 <= block[0] - 2 < len(splitIntoSpecificNames):
            startSurroundGene[currIndex].append(splitIntoSpecificNames[block[0] - 2])
        if 0 <= block[0] - 1 < len(splitIntoSpecificNames):
            startSurroundGene[currIndex].append(splitIntoSpecificNames[block[0] - 1])
        if 0 <= block[-1] + 1 < len(splitIntoSpecificNames):
            endSurroundGene[currIndex].append(splitIntoSpecificNames[block[-1] + 1])
        if 0 <= block[-1] + 2 < len(splitIntoSpecificNames):
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

    print(finalStart[:-1], "\n\n\n")
    print(finalEnd[:-1], "\n\n\n")

    numberArr = [[]]
    numberStart = [[]]
    numberEnd = [[]]
    currIndexArr = 0
    currIndexStartEnd = 0

    for i in range(len(transferArray)):
        splitFinal = finalArr[i][0].split(":")
        splitStart = finalStart[i][0].split(":")
        splitEnd = finalEnd[i][0].split(":")
        
        for i in range(3, len(splitFinal), 3):
            actualString = ""
            for char in splitFinal[i]:
                if char.isalpha():
                    break
                actualString += char
            numberArr[currIndexArr].append(actualString)
        numberArr.append([])
        currIndexArr += 1

        for i in range(3, len(splitStart), 3):

            actualString = ""

            for char in splitStart[i]:
                if char.isalpha():
                    break
                actualString += char
            numberStart[currIndexStartEnd].append(actualString)

            actualString = ""

            for char in splitEnd[i]:
                if char.isalpha():
                    break
                actualString += char
            numberEnd[currIndexStartEnd].append(actualString)

        numberStart.append([])
        numberEnd.append([])
        currIndexStartEnd += 1

    print(numberStart[:-1], "\n\n\n")
    print(numberEnd[:-1])

    return (numberArr[:-1], numberStart[:-1], numberEnd[:-1])
        




z = getMultiContig(content)
x = getSpecContig(z)
c = getNameContigFamily(x)
v = getContentsFromFile(recipientPath, "AveroniiAMC35", c)
b = parseResults(v)