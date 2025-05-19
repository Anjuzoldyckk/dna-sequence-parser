import os
import sys
def read_fasta(file,x=2): #x is tells if user wants individual or concat seqs
  current_seq = ""
  seqs = []
  with open(file,'r') as r:
    for line in r:
      line=line.strip()
      if line.startswith(">"):
        if current_seq:
          seqs.append(current_seq)
        current_seq = ""
      else:
        current_seq += line
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

def analyze_seq(sequence: str)->str: # u have to call it multiple times if the fasta read result is a list
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



if __name__ == "__main__": 
  file_fasta = input("Enter the name of the FASTA file (e.g., example.fasta): ")
  if not os.path.isfile(file_fasta):
    print(f" XXXXX File not found: {file_fasta}")
    sys.exit()
  x = int(input("hello :)\nif u want individual seq then enter 1 or if u want all seqs in the file combined then enter 2 : "))
  reading_fasta_result = read_fasta(file_fasta, x)
  if isinstance(reading_fasta_result, list):
    for i, seq in enumerate(reading_fasta_result):
      print(f"---Sequence {i + 1}---")
      analyze_seq(seq)
  else:
    analyze_seq(reading_fasta_result)

# edit dealing with empty files along with using regex to filter wild cases and returning dic so i can plot or output it in csv instead of js printing

    
