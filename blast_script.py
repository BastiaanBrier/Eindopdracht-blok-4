# Version:  1.0
# Date:     02-06-2018
# Author:   Project group 10
# Function: This module accepts 2 fastq files, one containing forward reads
#           and one containing reverse reads. It runs BLAST searches with these sequences
#           and passes the results on to the insert_script module for insertion into a mysql
#           database.

import warnings
import insert_script as insert
import time
from Bio.Blast import NCBIWWW
from Bio import Entrez, BiopythonExperimentalWarning, SeqIO
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio import SearchIO


def main():
    """This function calls the functions to read the files containing the forward and reverse reads in fastq format,
    BLAST the sequences in the files, retrieve the BLAST data and to store it in a mysql database by using the
    insert_script module. If something is wrong with the file names or file contents, an error message is printed and
    the program shuts down.
    """

    forward_file = "txt1.txt"
    reverse_file = "txt2.txt"
    try:
        headers1, seqs1, ascii_score1 = read_fastq(forward_file)
        headers2, seqs2, ascii_score2 = read_fastq(reverse_file)
    except IOError:
        print("Couldn't find the specified fastq files, shutting down the program.")
    except ValueError as e:
        print("Something is wrong with the contents of one of the files, the following error occurred: \n" + e.args[0])
        print("Please enter correct data files.")
    else:
        headers = headers1 + headers2
        seqs = seqs1 + seqs2
        ascii_score = ascii_score1 + ascii_score2

        try:
            with open("latest_header.txt", "r") as latest_header_file:
                latest_header = latest_header_file.readline().strip()
        except IOError:
            print("No file was found containing the latest header that was blasted. Starting from the beginning.")
        else:
            while headers[0] != latest_header:
                del headers[0]
                del seqs[0]
                del ascii_score[0]
            del headers[0]
            del seqs[0]
            del ascii_score[0]
            print("Continuing with blasting starting from header " + headers[0])

        for (header, sequence, quality) in zip(headers, seqs, ascii_score):
            blast_query(sequence)
            score, query_cover, identity, positives, evalue, protein_codes = read_xml()
            protein_names, protein_comments, organisms = prot_org_info(protein_codes)
            condata = insert.connect("sql7241825", "sql7.freemysqlhosting.net", "sql7241825", "7cCBxT27sc")
            send_to_database(header, sequence, quality, score, query_cover, identity, positives, evalue,
                             organisms, protein_names, protein_comments, protein_codes, condata[1])
            condata[0].commit()
            with open("latest_header.txt", "w") as header_save:
                header_save.write(header)
            condata[1].close()
            condata[0].close()
            time.sleep(12)


def send_to_database(header, sequence, quality, score, query_cover, identity, positives, evalue,
                     organisms, protein_names, protein_comments, protein_codes, cursor):
    """Accepts a single read header, sequence and quality score along with lists containing BLAST result
    data corresponding to that read and a cursor object and passes it along to the insert_script module.

    This function accepts the data of a single read and lists containing all BLAST results corresponding to that
    read. It then restructures these lists to the format accepted by the insert_read_and_data function of the
    insert_script module, where it is inserted into a mysql database.
    """
    read_list = [header, sequence, quality]
    match_list = []
    for i in range(len(score)):
        match_list.append([score[i], query_cover[i], identity[i], positives[i], evalue[i], organisms[i],
                           protein_names[i], protein_comments[i], protein_codes[i]])
    insert.insert_read_and_data(cursor, read_list, match_list)


def read_fastq(filename):
    """Accepts the name of a fastq file, reads it and returns the contents.

    This function accepts the name of a fastq file containing multiple DNA sequences and reads it.
    For each sequence, the corresponding ascii quality score is calculated. This function returns
    3 lists containing the headers, sequences and ascii quality scores.
    """
    seqs = []
    headers = []
    ascii_score = []
    records = SeqIO.parse(filename, "fastq")
    for record in records:
        seqs.append(record.seq)
        headers.append(record.id)
        subascii_score = 0
        for letter in record.letter_annotations["phred_quality"]:
            subascii_score += int(letter)
        ascii_score.append(subascii_score)
    return headers, seqs, ascii_score


def blast_query(seq):
    """Accepts a DNA sequence, BLASTS it and writes the BLAST results to an XML file.

    This function accepts a DNA sequence as input and runs a BLAST against the non redundant protein database using
    the blastx algorithm. The results of this BLAST search are then written to the 'blast_result' xml file.
    """
    print("Running BLAST search...")
    result_handle = NCBIWWW.qblast("blastx", "nr", seq, word_size=6, expect=0.0001, filter=True, matrix_name='BLOSUM62',
                                   gapcosts='11 1', hitlist_size=10)
    with open("blast_result.xml", "w") as out_handle:
        out_handle.write(result_handle.read())
    result_handle.close()
    print("BLAST search finished.")


def read_xml():
    """Reads the blast_result xml file and returns the contents.

    This function reads the blast_result xml file containing the results of a single BLAST search.
    Of every hit, the score, percentage query coverages, percentage identity, percentage positives, expect value and
    protein accession codes are stored in lists and returned.
    """
    score = []
    query_cover = []
    identity = []
    positives = []
    evalue = []
    protein_codes = []
    blast_qresult = SearchIO.read('blast_result.xml', 'blast-xml')
    for i in range(len(blast_qresult)):
        score.append(blast_qresult[i][0].bitscore)
        query_cover.append(round((float(blast_qresult[i][0].query_span) / float(blast_qresult[i][0].query_end) *100), 1))
        identity.append(round((float(blast_qresult[i][0].ident_num) / float(blast_qresult[i][0].hit_span) * 100), 1))
        positives.append(round((float(blast_qresult[i][0].pos_num) / float(blast_qresult[i][0].hit_span) * 100), 1))
        evalue.append(blast_qresult[i][0].evalue)
        protein_codes.append(blast_qresult[i].accession)

    return score, query_cover, identity, positives, evalue, protein_codes


def prot_org_info(protein_codes):
    """Accepts protein codes and returns corresponding names, comments and taxonomy information.

    This function accepts a list of protein accession codes and uses these to retrieve records from
    the NCBI protein database. From these records, the protein names, comments and source organisms are retrieved
    and returned in lists.
    """
    protein_names = []
    protein_comments = []
    organisms = []
    Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
    for i in range(len(protein_codes)):
        handle = Entrez.efetch(db="protein", id=protein_codes[i], rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        protein_names.append(record.description.split('[')[0])

        try:
            protein_comments.append(record.annotations["comment"])
        except KeyError:
            protein_comments.append("")

        organisms.append(record.annotations["organism"])
    return protein_names, protein_comments, organisms


main()
