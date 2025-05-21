import os
import sys
import re
import matplotlib.pyplot as plt
import pandas as pd
def read_fasta(file,x=2): #x is tells if user wants individual or concat seqs
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
        return seqs
    else:
        return final_seq

def gc_content(seq):
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


if __name__ == "__main__":
  file_fasta = input("Enter the name of the FASTA file (e.g., example.fasta): ")
  if not os.path.isfile(file_fasta):
    print(f" XXXXX File not found: {file_fasta}")
    sys.exit()
  x = int(input("hello :)\nif u want individual seq then enter 1 or if u want all seqs in the file combined then enter 2 : "))
  reading_fasta_result = read_fasta(file_fasta, x)
  if isinstance(reading_fasta_result, list) and len(reading_fasta_result) == 0:
      print("empty sequence bro!!")
      sys.exit()
  if isinstance(reading_fasta_result, str) and reading_fasta_result.strip() == "":
      print("empty sequence bro!!")
      sys.exit()
  result = []
  if isinstance(reading_fasta_result, list):
    for i, seq in enumerate(reading_fasta_result):
      print(f"---Sequence {i + 1}---")
      result_dict = analyze_seq(seq)   #capture dict into variable
      result.append(result_dict)
  else:
    print("---Concatenated Sequence---")
    result_dict = analyze_seq(reading_fasta_result)
    result.append(result_dict)
  print(result)
  graph = input("bro u wanna see cool grapht with len and gc content? (y/n)")
  if graph.lower().startswith("y") or graph == "y":
    graph_type = input("which type graph u want brotha? (1: Scatter, 2: Bar, 3: Pie): ")
    if graph_type == "1":
      lengths = [item['length'] for item in result]
      gc_contents = [item['gc_frac'] for item in result]
      plt.scatter(lengths, gc_contents)
      plt.xlabel("Sequence Length")
      plt.ylabel("GC Content (fraction)")
      plt.title("GC Content vs Sequence Length")
      plt.show()
    elif graph_type == "2":
      pass # bar code here
    elif graph_type == "3":
      pass # pie code here
    else:
        print("huh?")
  else:
    print("okay, cya <3 ")
  #panda shit
  df = pd.DataFrame(result)
  df.to_csv("seq_stats.csv", index=False)
  print("check for seq_stats.csv")

'''Sequence filtering:Let the user keep only sequences above/below a certain length or GC%.
histogram
motif search,
dna to protein '''
