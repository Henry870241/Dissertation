# Objective: To add BLAST search capability.
# Status: In progress.
# Notes: Does not have fastq ability yet, only fasta.

import os
import sys
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from PyQt6.QtWidgets import QApplication, QWidget, QVBoxLayout, QLabel, QPushButton, QFileDialog

class SimpleWindow(QWidget):
    def __init__(self):
        super().__init__()

        self.createWindow()

    def createWindow(self):
        self.setWindowTitle("Solution Attempt 02")
        self.setGeometry(250, 250, 500, 500)

        layout = QVBoxLayout()
        self.setLayout(layout)

        # Creating buttons.
        inputFileButton = QPushButton("Enter a .fasta/.fastq file.")
        inputFileButton.clicked.connect(self.loadFile)
        layout.addWidget(inputFileButton)

        blastSearchButton = QPushButton("Run a BLAST search.")
        blastSearchButton.clicked.connect(self.conductBlastSearch)
        layout.addWidget(blastSearchButton)

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

    def conductBlastSearch(self):

        blastQuery = SeqIO.read(self.fastaFile, "fasta")

        resultHandle = NCBIWWW.qblast("blastn", "nt", blastQuery.seq)

        blastRecords = NCBIXML.parse(resultHandle)
        for blastRecord in blastRecords:
            for alignment in blastRecord.alignments:
                print("Hit:", alignment.hit_id)
                for hsp in alignment.hsps:
                    print(f"Score: {hsp.score}")
                    print(f"Alignments: {hsp.sbjct}")
                    print("-----")

def main():
    app = QApplication(sys.argv)
    window = SimpleWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
