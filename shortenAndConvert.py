from Bio import SeqIO

def shortenAndConvert(input_file, output_file, num_sequences):
    with open(output_file, "w") as out_f:
        count = 0
        with open(input_file, "r") as in_f:
            for record in SeqIO.parse(in_f, "fastq"):
                if count < num_sequences:
                    SeqIO.write(record, out_f, "fasta")
                    count += 1
                else:
                    break

    print(f"Saved {count} sequences to {output_file}.")