#!/usr/bin/python3

import os
import sys
import csv

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage:", sys.argv[0], " <portion of name>")
        exit(1)

    count = 0
    author_query = sys.argv[1].lower()
    author_name = ""

    # find the all possible names of the author
    author_table = "author_normalization_table_reversed_indexed.tsv"
    with open(author_table, 'r') as infile:
        author_reader = csv.reader(infile, delimiter="\t")
        for row in author_reader:
            if (author_query in row[1].lower()):
                print(row[0] + '\t' + row[1])
                count += 1

    print("Found", count, "records.")
