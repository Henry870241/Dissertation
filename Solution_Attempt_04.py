# Objective: This aims to plot a bar chart using the BLAST results: alignment.hit_def, and hsp.score.
# Status: Success.
# The bar chart does not accept redundant data; any organisms with the same name (even if they have different IDs) will be seen as one.

import os
import sys

from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from PyQt6.QtWidgets import QApplication, QWidget, QVBoxLayout, QLabel, QPushButton, QFileDialog

from matplotlib import pyplot

organismNames = []
hitScores = []

class MainWindow(QWidget):
    def __init__(self):
        super().__init__()

        self.createWindow()

    def createWindow(self):
        self.setWindowTitle("Solution Attempt 04")
        self.setGeometry(250, 250, 500, 500)

        layout = QVBoxLayout()
        self.setLayout(layout)

        # Creating buttons.
        inputFileButton = QPushButton("Enter a .fasta/.fastq file.")
        inputFileButton.clicked.connect(self.loadFile)
        layout.addWidget(inputFileButton)

        self.filePathnameLabel = QLabel()
        layout.addWidget(self.filePathnameLabel)

        blastSearchButton = QPushButton("Run a BLAST search.")
        blastSearchButton.clicked.connect(self.conductBlastSearch)
        layout.addWidget(blastSearchButton)

        createBarChartButton = QPushButton("Create Bar Chart.")
        createBarChartButton.clicked.connect(createBarChart)
        layout.addWidget(createBarChartButton)

    def loadFile(self):
        file_filter = 'Data File (*.fasta *.fastq)'
        response = QFileDialog.getOpenFileName(
            parent=self,
            caption='Select a file',
            directory=os.getcwd(),
            filter=file_filter,
            initialFilter= 'Data File (*.fasta *.fastq)'
        )

        fileName = str(response[0])
        self.fastaFile = fileName

        self.filePathnameLabel.setText("File: " + self.fastaFile + " Uploaded")

    def conductBlastSearch(self):

        print("Running BLAST search...")

        try:
            blastQuery = SeqIO.read(self.fastaFile, "fasta")
        except AttributeError:
            print("Please enter a .fasta file first.")

        resultHandle = NCBIWWW.qblast("blastn", "nt", blastQuery.seq)

        blastRecords = NCBIXML.parse(resultHandle)

        for blastRecord in blastRecords:
            for i, alignment in enumerate(blastRecord.alignments):
                if i == 5: break
                print(f"Hit ID: {alignment.hit_id}")
                print(f"Hit Description: {alignment.hit_def.split()[:2]}")

                hitDefToString = " ".join(alignment.hit_def.split()[:2])

                organismNames.append(hitDefToString)

                for hsp in alignment.hsps:
                    print(f"Hit Score: {hsp.score}")
                    hitScores.append(hsp.score)
                    print("-----")

        print("BLAST search complete.")

def createBarChart():
    pyplot.bar(organismNames, hitScores)
    pyplot.title("BLAST results.")
    pyplot.xlabel("Organism Names")
    pyplot.ylabel("Hit Scores")
    pyplot.show()

def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()