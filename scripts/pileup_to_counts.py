import sys
import pandas as pd

# Input pileup file from command line
pileup_file = sys.argv[1]

# We will store results here
rows = []

with open(pileup_file) as f:
    for line in f:
        fields = line.strip().split("\t")

        ref_name = fields[0]
        position = int(fields[1])
        ref_base = fields[2].upper()
        depth = int(fields[3])
        bases = fields[4]

        # Clean the bases string step by step

        i = 0
        cleaned = ""

        while i < len(bases):
            base = bases[i]

            if base == "^":
                # skip next character (mapping quality)
                i += 2
                continue

            if base == "$":
                i += 1
                continue

            if base in "+-":
                # insertion or deletion
                i += 1
                indel_length = ""
                while bases[i].isdigit():
                    indel_length += bases[i]
                    i += 1
                indel_length = int(indel_length)
                i += indel_length
                continue

            cleaned += base
            i += 1

        # Now convert matches to actual base
        cleaned = cleaned.replace(".", ref_base)
        cleaned = cleaned.replace(",", ref_base)

        cleaned = cleaned.upper()

        A = cleaned.count("A")
        C = cleaned.count("C")
        G = cleaned.count("G")
        T = cleaned.count("T")

        rows.append([position, ref_base, A, C, G, T])

df = pd.DataFrame(rows, columns=["pos", "ref", "A", "C", "G", "T"])

df.to_csv(pileup_file + ".counts.tsv", sep="\t", index=False)
