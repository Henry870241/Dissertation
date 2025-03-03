1. Solution_Attempt_01.py is a simple PyQt program that can retrieve the pathname of a .fasta/.fastq file.
  
2. Solution_Attempt_02.py shares the GUI of Solution_Attempt_01.py but includes a BLAST search function. So far the program has no error handling and does not support .fastq BLAST queries.

3. Solution_Attempt_03.py included some new but minor stylistic changes not worth uploading.

4. Solution_Attempt_04.py incorporates the matlablib library to create bar graphs using data from a BLAST search output. The data used includes hit_def (the description of a matching sequence from a reference database), and hsp.score (a metric which represents how well the query sequence matches to the reference sequence). 
