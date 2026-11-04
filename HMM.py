from collections import Counter 

#Draft 

# DNA w/ CPG islands from: 

# Chromosome range: chr12:56,741,761-56,761,390 

# Bases input into txt file (CpG island bases are capitalized)

# https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm39&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr12%3A56741761%2D56761390&hgsid=4249935935_g40he2bWdsP5AXPaa4uR8ELa1pbw 

def pull_data(path):
    raw = "".join(line.strip() for line in open(path) if not line.startswith(">"))
    extracted_observations = [b.upper() for b in raw]
    extracted_hidden_states = ["CpG" if b.isupper() else "non-CpG" for b in raw] 
    return extracted_observations, extracted_hidden_states 

extracted_observations, extracted_hidden_states = pull_data("mm39ExtendedDNAOutput.txt")

#print(extracted_observations)
#print(extracted_hidden_states)


hidden_states = extracted_hidden_states 

dNA_bases = extracted_observations 

#Figuring out emission probs: 

#A, T, C, G Counts for CpG and A, T, C, G Counts for non-CpG 

counts = {"CpG": Counter(), "non-CpG": Counter()}

for base, state in zip(dNA_bases, hidden_states):
    counts[state][base] += 1 

#print(counts)

#A, T, C, G Probs for CpG and A, T, C, G Probs for non-CpG  

emission = {}

for state, count in counts.items():
    total = sum(count.values())
    emission[state] = [round(count["A"] / total, 4), 
                        round(count["T"] / total, 4),
                        round(count["C"] / total, 4),
                        round(count["G"] / total, 4),] 

for state in emission:
    print(state, emission[state])









