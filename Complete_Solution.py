import math
import os
import sys
from collections import Counter

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from PyQt6.QtWidgets import (QApplication, QWidget,
                             QVBoxLayout, QLabel,
                             QPushButton, QFileDialog)
from cycler import cycler
from matplotlib import pyplot

def addHitCounts(originalList: list, newList: list) -> None:
    for itemName, itemQuantity in Counter(originalList).items():
        newList.append(itemQuantity)

def calculateSDI_Value(aList: list) -> int:
    requiredValues = []
    for alignmentName, hitCount in Counter(aList).items():
        requiredValues.append((hitCount / len(aList)) * math.log(hitCount / len(aList)))
    return abs(sum(requiredValues))

def getLowerSum(listOne: list, listTwo: list) -> int:
    return Counter(listOne).most_common()[-1][1] + Counter(listTwo).most_common()[-1][1]

def calculateBrayCurtis(populationOne: list, populationTwo: list) -> float:
    return 1 - (2 * getLowerSum(populationOne, populationTwo) / (len(populationOne) + len(populationTwo)))

def wrapText(aList: list) -> list:
    for i in range(len(aList)):
        aList[i] = aList[i].replace(' ', '\n')
    return aList

def shortenFileNames(aList):
    for i in range(len(aList)):
        aList[i] = aList[i].split('/')[-1]


def centerBarLabel(listX, listY):
    for i in range(len(listX)):
        pyplot.text(i, listY[i] // 2, listY[i], ha ='center')

class MainWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.labelText = "Status: Empty"
        self.tooManyFilesTxt = "You have entered too many files."

        self.inputOneLabel = None
        self.inputOneButton = None
        self.inputTwoLabel = None
        self.inputTwoButton = None

        self.blastSearchButton = None
        self.createBarChartButton = None
        self.SDI_Button = None
        self.brayCurtisButton = None

        self.fileNames = []

        self.alignmentNames01 = []
        self.alignmentNames02 = []

        self.hitScores01 = []
        self.hitScores02 = []

        self.fileOneHitCounts = []
        self.fileTwoHitCounts = []

        self.fileOneSpeciesQuantity = None
        self.fileTwoSpeciesQuantity = None

        self.totalHits = None

        self.barWidth = 0.5

        self.createWindow()

    def createWindow(self) -> None:
        self.setWindowTitle("Solution Attempt 8")
        self.setGeometry(500, 250, 500, 500)

        layout = QVBoxLayout()
        self.setLayout(layout)

        # Creating the Labels and Buttons for file input.
        self.inputOneLabel = QLabel()
        self.inputOneLabel.setText(self.labelText)
        layout.addWidget(self.inputOneLabel)

        self.inputOneButton = QPushButton("Enter a .fasta file.")
        self.inputOneButton.clicked.connect(lambda: self.loadFile(0))
        layout.addWidget(self.inputOneButton)

        self.inputTwoLabel = QLabel()
        self.inputTwoLabel.setText(self.labelText)
        layout.addWidget(self.inputTwoLabel)

        self.inputTwoButton = QPushButton("Enter a .fasta file.")
        self.inputTwoButton.clicked.connect(lambda: self.loadFile(1))
        layout.addWidget(self.inputTwoButton)

        self.blastSearchButton = QPushButton("Blast Search")
        self.blastSearchButton.clicked.connect(self.blastSearch)
        layout.addWidget(self.blastSearchButton)

        # Creating buttons for functionality
        self.createBarChartButton = QPushButton("Create bar chart.")
        self.createBarChartButton.clicked.connect(self.createBarChart)
        layout.addWidget(self.createBarChartButton)

        self.SDI_Button = QPushButton("Compare Shannon Diversity indexes.")
        self.SDI_Button.clicked.connect(self.createSDI_Graph)
        layout.addWidget(self.SDI_Button)

        self.brayCurtisButton = QPushButton("Calculate Bray Curtis Index.")
        self.brayCurtisButton.clicked.connect(self.outputBrayCurtis)
        layout.addWidget(self.brayCurtisButton)

    # Creating the load file functionality.
    def loadFile(self, fileNumber: int) -> None:

        fileFilter = 'Data File (*.fasta)'

        response = QFileDialog.getOpenFileName(
            parent=self,
            caption='Select a file',
            directory=os.getcwd(),
            filter=fileFilter,
            initialFilter='Data File (*.fasta)'
        )

        fileName = str(response[0])
        self.fileNames.append(fileName)

        if fileNumber == 0:
            self.inputOneLabel.setText(fileName)
        elif fileNumber == 1:
            self.inputTwoLabel.setText(fileName)

    # Conduct a BLAST search.
    def blastSearch(self) -> None:
        if len(self.fileNames) == 2: # Making sure there have been two .fasta files inputted.

            for fileIndex, fileName in enumerate(self.fileNames): # Loops through each .fasta file.

                print(f"Running BLAST search on file: {fileName.split('/')[-1]}")

                sequenceLimit = 2
                for blastQueryIndex, blastQuery in enumerate(SeqIO.parse(fileName, 'fasta')): # Loops through each sequence in the .fasta file.
                    if blastQueryIndex == sequenceLimit: break

                    # Perform BLAST search
                    resultHandle = NCBIWWW.qblast("blastn", "nt", blastQuery.seq)
                    blastRecords = NCBIXML.parse(resultHandle)

                    alignmentLimit = 3
                    for blastRecord in blastRecords:
                        for blastRecordIndex, alignment in enumerate(blastRecord.alignments):  # Loop through each alignment in the sequence
                            if blastRecordIndex == alignmentLimit:
                                break

                            print(f"Hit ID: {alignment.hit_id}")
                            alignmentName = " ".join(alignment.hit_def.split()[:2])
                            print(f"Hit Description: {alignmentName}\n")

                            if fileIndex == 0:
                                for hsp in alignment.hsps:
                                    self.alignmentNames01.append(alignmentName)
                                    self.hitScores01.append(hsp.score)
                            elif fileIndex == 1:
                                for hsp in alignment.hsps:
                                    self.alignmentNames02.append(alignmentName)
                                    self.hitScores02.append(hsp.score)

                            print(f"\nFile name: {fileName.split('/')[-1]}, Sequence: {blastQueryIndex}, Alignment: {blastRecordIndex}.\n")

                print(f"\nFile name: {fileName.split('/')[-1]} BLAST completed.")
            print("\nBlast search complete.\n")

            addHitCounts(self.alignmentNames01, self.fileOneHitCounts)
            addHitCounts(self.alignmentNames02, self.fileTwoHitCounts)

        else:
            print(self.tooManyFilesTxt)

    def createBarChart(self) -> None:
        if len(self.alignmentNames01) == len(self.hitScores01) and len(self.alignmentNames02) == len(self.hitScores02):
            pyplot.figure(figsize=(20, 9))

            wrappedAlignmentNames01 = wrapText(self.alignmentNames01)
            wrappedAlignmentNames02 = wrapText(self.alignmentNames02)

            pyplot.barh(wrappedAlignmentNames01, self.hitScores01, label=self.fileNames[0].split("/")[-1], height=self.barWidth)
            pyplot.barh(wrappedAlignmentNames02, self.hitScores02, label=self.fileNames[1].split("/")[-1], height=self.barWidth)

            pyplot.legend()
            pyplot.title("BLAST High Scoring Points")
            pyplot.xlabel("Alignment Hit Score")
            pyplot.ylabel("Alignment Names")
            centerBarLabel(self.alignmentNames01, self.hitScores01)
            pyplot.show()
        else:
            print("Mismatch of alignment names and scores.")

    def createSDI_Graph(self) -> None:
        yValues = [calculateSDI_Value(self.fileOneHitCounts), calculateSDI_Value(self.fileTwoHitCounts)]

        shortFileNames = []
        for fileName in self.fileNames:
            shortFileNames.append(fileName.split('/')[-1])

        graphColours = ['#eb3434', '#3c36ff']

        pyplot.barh(shortFileNames, yValues, label="SDI Comparison", color = graphColours, height = self.barWidth)
        pyplot.title("SDI Values of Samples")
        pyplot.xlabel("SDI Value")
        pyplot.ylabel("Samples")
        pyplot.show()

    def outputBrayCurtis(self) -> None:
        brayCurtisOutputLabel = QLabel(f"Bray-Curtis Value: {round(calculateBrayCurtis(self.fileOneHitCounts, self.fileTwoHitCounts), 2)}")
        self.layout().addWidget(brayCurtisOutputLabel)

def main() -> None:
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()