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

    # find the all possible names of the author
    author_table = "author_normalization_table_reversed_indexed.tsv"
    with open(author_table, 'r') as infile:
        author_reader = csv.reader(infile, delimiter="\t")
        for row in author_reader:
            if (author_id == int(row[0])):
                print(row[0] + '\t' + row[1])
                exit(0)


    print("Error: Author ID %d not found" % author_id)
