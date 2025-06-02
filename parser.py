import os
import sys
import re
import matplotlib.pyplot as plt
import pandas as pd
def read_fasta(file,x=2): #x is tells if user wants individual or concat seqs
  """
    Read a FASTA file and grab sequences, bro.
    Args:
        file (str): Path to your FASTA file.
        x (int): 1 to get individual seqs, else concat all.
    Returns:
        If x == 1: (headers, [seq1, seq2, ...])
        Else:     (headers, "all_sequences_combined")
    """
  current_seq = ""
  seqs = []
  headers=[]
  with open(file,'r') as r:
    for line in r:
      line=line.strip()
      if line.startswith(">"):
        headers.append(line[1:])
        if current_seq:
          seqs.append(current_seq)
        current_seq = ""
      else:
        cleaned = re.sub(r'[^ATGCatgc]', '', line)
        current_seq += cleaned
    if current_seq: #for adding last seq in current seq
      seqs.append(current_seq)
    final_seq = "".join(seqs)
    if x == 1:
        return headers,seqs
    else:
        return headers,final_seq

def gc_content(seq):
  """
    Compute GC fraction & percent, dude.
    Args:
        seq (str): DNA string (can be mixed case/whitespace).
    Returns:
        (gc_fraction, "xx.xx%")
    """
  seq = seq.upper().strip()
  total_len = len(seq)
  if total_len == 0:
    return 0, "0.00%"
  g = seq.count('G')
  c = seq.count('C')
  gc= g+c
  gc_fraction = gc / total_len
  gc_percent = gc_fraction * 100
  return round(gc_fraction, 4), f"{gc_percent:.2f}%"

def at_content(seq):
  """
    Compute AT fraction & percent, fam.
    Args:
        seq (str): DNA string (can be mixed case/whitespace).
    Returns:
        (at_fraction, "xx.xx%")
    """
  seq = seq.upper().strip()
  total_len = len(seq)
  if total_len == 0:
    return 0, "0.00%"
  a = seq.count('A')
  t = seq.count('T')
  at = a + t
  at_fraction = at / total_len
  at_percent = (at_fraction * 100)
  return round(at_fraction, 4), f"{at_percent:.2f}%"

def analyze_seq(sequence: str)-> dict:
  """
    Print & collect basic stats for one DNA seq, brotha.
    Args:
        sequence (str): Single DNA string (mixed case/whitespace OK).
    Returns:
        dict: {
            "length": int,
            "a_count": int, "t_count": int,
            "g_count": int, "c_count": int,
            "a_perc": float, "t_perc": float,
            "g_perc": float, "c_perc": float,
            "at_frac": float, "at_perc": "xx.xx%",
            "gc_frac": float, "gc_perc": "xx.xx%"
        }
        Prints warning & returns None if empty.
    """
  #it gets called multiple times if the fasta read result is a list,
  #it return none for now but later if u want to use for graphs and all u have to denote that it returns a dict
  sequence = sequence.upper().strip()
  length = len(sequence)
  if length == 0:
        print("empty sequence bro!!")
        return
  a = sequence.count('A')
  t = sequence.count('T')
  g = sequence.count('G')
  c = sequence.count('C')
  #cal each base percentage
  a_perc = (a / length) * 100
  t_perc = (t / length) * 100
  g_perc = (g / length) * 100
  c_perc = (c / length) * 100
  #calling defined funs
  at_frac, at_perc = at_content(sequence)
  gc_frac, gc_perc = gc_content(sequence)
  print(f"Sequence length: {length}")
  print(f"A count: {a},({a_perc:.2f}%)")
  print(f"T count: {t},({t_perc:.2f}%)")
  print(f"G count: {g},({g_perc:.2f}%)")
  print(f"C count: {c},({c_perc:.2f}%)")
  print(f"AT content: {at_frac},({at_perc})")
  print(f"GC content: {gc_frac},({gc_perc})")
  result = {"length" : length,
            "a_count" : a,
            "t_count" : t,
            "g_count" : g,
            "c_count" : c,
            "a_perc" : a_perc,
            "t_perc" : t_perc,
            "g_perc" : g_perc,
            "c_perc" : c_perc,
            "at_frac" : at_frac,
            "at_perc" : at_perc,
            "gc_frac" : gc_frac,
            "gc_perc" : gc_perc}
  return result

def plot_gc_vs_length_scatter(result):
   """
    Scatter: seq length vs GC fraction, y'all.
    Args:
        result (list of dict): output from analyze_seq per seq.
    """
  lengths = [item['length'] for item in result]
  gc_contents = [item['gc_frac'] for item in result]
  plt.scatter(lengths, gc_contents)
  plt.xlabel("Sequence Length")
  plt.ylabel("GC Content (fraction)")
  plt.title("GC Content vs Sequence Length")
  plt.show()

def plot_gc_content_bar(result):
  """
    Bar chart of GC fraction per seq, fam.
    Args:
        result (list of dict): output from analyze_seq per seq.
    """
  labels = [f"seq{i+1}" for i in range(len(result))]
  gc_contents = [item['gc_frac'] for item in result]
  plt.bar(labels, gc_contents)
  plt.xlabel("sequence")
  plt.ylabel("gc content (fraction)")
  plt.title("GC Content per Sequence")
  plt.show()

def plot_base_pie(result, seq_index, headers=None):
   """
    Pie chart of A/T/G/C counts for one seq, bro.
    Args:
        result (list of dict): output from analyze_seq per seq.
        seq_index (int): which seq to plot (0-based).
        headers (list of str, optional): seq names to show in title.
    """
  base_counts = [
      result[seq_index]['a_count'],
      result[seq_index]['t_count'],
      result[seq_index]['g_count'],
      result[seq_index]['c_count'],
      ]
  base_labels = ['A', 'T', 'G', 'C']
  plt.pie(base_counts, labels=base_labels, autopct='%1.1f%%')
  if headers and seq_index < len(headers):
    plt.title(f"Base Composition ({headers[seq_index]})")
  else:
    plt.title(f"Base Composition (Sequence {seq_index+1})")
  plt.show()

def plot_gc_histogram(result,bins = 8):
  """
    Histogram of GC fraction across all seqs, dude.
    Args:
        result (list of dict): output from analyze_seq per seq.
        bins (int): number of bins for histogram.
    """
  gc_list = [r['gc_frac'] for r in result]
  plt.figure()
  plt.hist(gc_list, bins = bins, edgecolor = "black")
  plt.xlabel("GC Content(in frac)")
  plt.ylabel("Frequency/no of seqs")
  plt.title("GC Content distribution in histogram")
  plt.show()


def plot_gc_boxplot(result):
   """
    Boxplot of GC fraction for all seqs, brotha.
    Args:
        result (list of dict): output from analyze_seq per seq.
    """
  gc_list = [r['gc_frac'] for r in result]
  plt.figure()
  plt.boxplot(gc_list, vert = True) # vertical by default, change this to false for horizontal boxplot
  plt.xlabel("all seqs")
  plt.ylabel("GC Content")
  plt.title("GC Content distribution in boxplot")
  plt.show()
def translater(seq, frame=0):
    table = """
TTT F Phe      TCT S Ser      TAT Y Tyr      TGT C Cys  
TTC F Phe      TCC S Ser      TAC Y Tyr      TGC C Cys  
TTA L Leu      TCA S Ser      TAA * Ter      TGA * Ter  
TTG L Leu i    TCG S Ser      TAG * Ter      TGG W Trp  

CTT L Leu      CCT P Pro      CAT H His      CGT R Arg  
CTC L Leu      CCC P Pro      CAC H His      CGC R Arg  
CTA L Leu      CCA P Pro      CAA Q Gln      CGA R Arg  
CTG L Leu i    CCG P Pro      CAG Q Gln      CGG R Arg  

ATT I Ile      ACT T Thr      AAT N Asn      AGT S Ser  
ATC I Ile      ACC T Thr      AAC N Asn      AGC S Ser  
ATA I Ile      ACA T Thr      AAA K Lys      AGA R Arg  
ATG M Met i    ACG T Thr      AAG K Lys      AGG R Arg  

GTT V Val      GCT A Ala      GAT D Asp      GGT G Gly  
GTC V Val      GCC A Ala      GAC D Asp      GGC G Gly  
GTA V Val      GCA A Ala      GAA E Glu      GGA G Gly  
GTG V Val      GCG A Ala      GAG E Glu      GGG G Gly  
"""
    table = table.upper().split()
    table = [t for t in table if t != "I"] #to remove stray tokens
    token = [table[i:i+3] for i in range(0, len(table), 3)]
    codon_table = {codon: (single, three) for codon, single, three in token}
    proteins=[]
    for i in range(frame,len(seq)-2,3):
      codon = seq[i:i+3]
      aa = codon_table.get(codon, ("X","X")) #x being default, and not to crash the prog by crashing
      aa = aa[0] #to only allow the single in the list omitting the tuple 
      proteins.append(aa)
    proteins = "".join(proteins)
    print(f"the translated protein is {proteins}")
    return proteins

if __name__ == "__main__":
  file_fasta = input("Enter the name of the FASTA file (e.g., example.fasta): ")
  if not os.path.isfile(file_fasta):
    print(f" XXXXX File not found: {file_fasta}")
    sys.exit()
  x = int(input("hello :)\nif u want individual seq then enter 1 or if u want all seqs in the file combined then enter 2 : "))
  headers, reading_fasta_result = read_fasta(file_fasta, x)
  if isinstance(reading_fasta_result, list) and len(reading_fasta_result) == 0:
      print("empty sequence bro!!")
      sys.exit()
  if isinstance(reading_fasta_result, str) and reading_fasta_result.strip() == "":
      print("empty sequence bro!!")
      sys.exit()
  result = []
  if isinstance(reading_fasta_result, list):
    for i, seq in enumerate(reading_fasta_result):
      print(f"---Sequence {i + 1}---({headers[i]})")
      result_dict = analyze_seq(seq)   #capture dict into variable
      protein = translater(seq, frame=0)
      result_dict["protein_frame1"] = protein
      result.append(result_dict)
  else:
    print("---Concatenated Sequence---")
    result_dict = analyze_seq(reading_fasta_result)
    protein = translater(reading_fasta_result, frame=0)
    result_dict["protein_frame1"] = protein
    result.append(result_dict)
  print(result)
  graph = input("bro u wanna see cool grapht with len and gc content? (y/n)")
  if graph.lower().startswith("y") or graph == "y":
    graph_type = input("which type graph u want brotha? (1: Scatter, 2: Bar, 3: Pie, 4: histogram, 5: boxplot): ")
    if graph_type == "1":
      plot_gc_vs_length_scatter(result)
    elif graph_type == "2":
      plot_gc_content_bar(result)
    elif graph_type == "3":
      max_to_show = 20 #how many headers can be seen
      num_to_show = min(len(headers), max_to_show)
      print(f"there are {len(headers)} in ur file brotha")
      print(f"showing the headers of {num_to_show} seqs")
      for i in range(num_to_show):
        print(f"{i+1}: {headers[i]}")
      seq_index = None
      while True:
        user_input = input(f"gimme the seq num (1 to {num_to_show} for pie chart, or q for quitting)").strip()
        if user_input == 'q':
          break
        try:
          choice_seq = int(user_input)
          if 1 <= choice_seq <= num_to_show:
            seq_index = choice_seq - 1
            break
        except ValueError:
          print("print a valid num bruh, or q for quit")
      if seq_index is not None:
        plot_base_pie(result, seq_index, headers)

    elif graph_type == '4':
      plot_gc_histogram(result)
    elif graph_type == '5':
      plot_gc_boxplot(result)
    else:
        print("huh?")
  else:
    print("okay, cya <3 ")
  #panda export for csv
  df = pd.DataFrame(result)
  df.to_csv("seq_stats.csv", index=False)
  print("check for seq_stats.csv")

'''
give max gc content as output, annotate the outliers
learn about argparse and change it fro0m input() to argparse
Sequence filtering:Let the user keep only sequences above/below a certain length or GC%.
Write some basic test sequences to test your script automatically.
Refactor plotting: Later, consider putting all plot functions in a plots.py file (when the script grows).
motif search,
Write filtered sequences back to a new FASTA file.
Show a table (or output file) with frequencies of all bases for each sequence.
Implement a function to get the reverse complement of a DNA sequence.
Count and report on non-ATGC bases (e.g., N, R, Y).
HEATMAP
set up vs code
Write a tiny FASTA file (e.g. one sequence “ATGCNN”), load it, and assert counts/percentages with pytest'''
