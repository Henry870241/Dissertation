import os
import sys
from PyQt6.QtWidgets import QApplication, QWidget, QVBoxLayout, QLabel, QPushButton, QFileDialog


class SimpleWindow(QWidget):
    def __init__(self):
        super().__init__()

        self.createWindow()

    def createWindow(self):
        self.setWindowTitle("Solution Attempt 01")
        self.setGeometry(250, 250, 500, 500)

        layout = QVBoxLayout()
        self.setLayout(layout)

        # Creating buttons.
        inputFileButton = QPushButton("Enter a .fasta/.fastq file.")
        inputFileButton.clicked.connect(self.loadFile)
        layout.addWidget(inputFileButton)

        blastSearchButton = QPushButton("Run a BLAST search.")
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
        fileName = str(response)
        print(fileName)


def main():
    app = QApplication(sys.argv)
    window = SimpleWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
