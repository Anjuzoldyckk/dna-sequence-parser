# dna-sequence-parser

### Python script to flex on DNA FASTA files. Not just GC content, bro—we go way deeper.

---

## What does this thing do?

- Reads your FASTA file, doesn’t matter if it’s one long boi or a bunch of sequences
- Gives you all the stats: length, base counts, percentages, GC/AT content, you name it
- Plots GC content vs. length (scatter plot, looks sick)
- Lets you choose bar or pie charts for base composition if you want
- Dumps all results into a CSV file (so you can vibe with it in Excel or Google Sheets)
- Asks you what you want at every step—total control
- Runs on Python, uses pandas and matplotlib (install with `pip install pandas matplotlib` if you don’t have them, duh)

---

## How to run it

1. Drop your FASTA file in the same folder as the script
2. Run:
    ```bash
    python dna_parser.py
    ```
3. Follow the prompts (type in your file name, pick if you want individual or one big sequence, choose your graph, whatever)
4. Open the `sequence_stats.csv` file to flex your results in Excel, Google Sheets, Notepad, whatever

---

## Test Cases

All the test FASTA files you need are chilling in the `/test_cases` folder.  
Run the script on any of those files to see what’s up.  


## Example FASTA to copy-paste

seq1
ATGCGTACGTTAGCTAGCTAGCTAGCTGACTGATCGTAGCTAGCTAGC
seq2
ATATATATATATATATATATATATATATATATATATATATATATATA
seq3
GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG


## Example CSV output

length,a_count,t_count,g_count,c_count,a_perc,t_perc,g_perc,c_perc,at_frac,at_perc,gc_frac,gc_perc
61,15,15,16,15,24.59,24.59,26.23,24.59,0.4918,49.18%,0.5082,50.82%
60,15,15,15,15,25.0,25.0,25.0,25.0,0.5,50.00%,0.5,50.00%

## Requirements

- Python 3.x
- `pip install pandas matplotlib`

---

## License

MIT. Do whatever, just don’t blame me if your computer explodes.

---

Big love. PRs and memes welcome.


