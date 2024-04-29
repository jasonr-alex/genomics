### Hi there 👋, my name is Jason 
#### I am a passionate researcher exploring the world of genomics and bioinformatics. 

I created this page to document my coding journey and share my love of genomics and bioinformatics. The following details describe the code featured in this repository with examples. 

### Complement and Reverse Complement  
In biology, the central dogma is the main process by which DNA is transcribed to RNA and translated to a protein (not all though). The importance of determining the complement of DNA is that it allows us to carry out the same process by which our body determines what bases to add to a coding strand based on a template strand. 

![The Central Dogma Path](https://github.com/jasonr-alex/genomics/blob/main/Central-dogma.jpeg)

Reverse transcription, however, conflicts against this central dogma; reverse transcription is an important viral replication process that allows for insertion in the host's genome. 

### Complementation and Reverse Complementation Code (Python)
```
def Complement(s):
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = t + complement[base]
    return t

def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

#### Example of Complement Use:
def Complement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = t + complement[base]
    return t

print(Complement("TGTGGTGACACATG")) 
##### Result: ACACCACTGTGTAC

#### Example of Reverse-Complement Use:

def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t  
print(reverseComplement("TGTGGTGACACATG"))

##### Result: CATGTGTCACCACA
```
#### Handling FASTQ  and FASTA Files
The major files in which genetic data is stored are referred to as fastq and fasta files. The primary format mentioned is a file format that stores the reads and the quality score associated with the reads; this is the major source of data created by NGS sequencers. The last format mentioned, fasta files, are found on the National Center of Biotechnology Information (NCBI) and store assembled genomic data. 


##### Parsing FASTA Files
```
def readGenome(filename): #Function for parsing a fasta file. 
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome
```
##### Parsing FASTQ Files
```
def readFastq(filename): #Function for parsing a fastq file with sequence and quality score.
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities
```
Both these functions when executed should return a list of the genomic data, but in the case of the fastq function, it will also have the addition of a quality score associated with each base or read. It is also important to state that sequencing errors are common occurrences, and having a higher quality score guarantees that the observed read is not due to error. 
 
#### Algorithms for Genomics 
Genomic data is vast and being able to extract meaningful data requires the necessary algorithms to help process that data or search for a particular parameter. By implementing indexes and matrices, which work well together, we can analyze data in a reasonable time without taking years. The algorithms detailed below allow for determining the changes for two different strings by edit distance and boyer-moore; while also, being able to find a similar string of interest in a genome of interest by naive exact matching.  

```
def editDistance(x, y): #The amount of edits necessary to make two sequences similar. 
    # Create a distance matrix
    D = []
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    # Initialize first row and column of matrix
    for i in range(len(x)+1):
        D[i][0] = i
    for i in range(len(y)+1):
        D[0][i] = i
    # Fill in the rest of the matrix
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    # Edit distance is the value in the bottom right corner of the matrix
    return D[-1][-1]

def boyer_moore(p, p_bm, t): #Method termed Boyer-Moore that only allows substitutions for getting the sequences to become similar.  
    """ Do Boyer-Moore matching. p=pattern, t=text,
        p_bm=BoyerMoore object for p """
    i = 0
    occurrences = []
    while i < len(t) - len(p) + 1:
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift

def naive(p, t): #Algorithm that checks every read but guarantees accuracy when determining if similar sequences are present when comparing genomes. 
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences
```



