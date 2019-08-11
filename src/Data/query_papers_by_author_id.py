#!/usr/bin/python3

import os
import sys
import csv

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage:", sys.argv[0], " <author ID>")
        exit(1)

    author_id = int(sys.argv[1])
    author_name = ""

    # find the name of the author
    author_table = "author_normalization_table_reversed_indexed.tsv"
    with open(author_table, 'r') as infile:
        author_reader = csv.reader(infile, delimiter="\t")
        for row in author_reader:
            if (int(row[0]) == author_id):
                author_name = row[1]
                print("Searching for papers by:", author_name)
                break

    if (author_name == ""):
        print("Error: Invalid author ID:", author_id)
        exit(2)

    count = 0
    paper_table = "all_hpec_records_cleaned_normalized_authors_indexed.tsv"
    with open(paper_table, 'r') as infile:
        paper_reader = csv.reader(infile, delimiter="\t")
        for row in paper_reader:
            author_list = row[2].split("|")
            for author_affiliation in author_list:
                author, affiliation = author_affiliation.split("/")
                if (author == author_name):
                    print("Paper ID:", row[0])
                    print("Year:", row[3])
                    print("Title: \"" + row[1] + "\"")
                    print("Authors:\n" + "\n".join(row[2].split('|')))

                    print(" ")
                    count += 1
                    break

    print("Found", count, "papers.")
