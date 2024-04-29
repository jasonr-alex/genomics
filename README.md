### Hi there 👋, my name is Jason 
#### I am a passionate researcher exploring the world of genomics and bioinformatics. 

I created this page to document my coding journey and share my love of genomics and bioinformatics. The following details describe the code featured in this repository with examples. 

### Complement and Reverse Complement  
In biology, the central dogma is the main process by which DNA is transcribed to RNA and translated to a protein (not all though). The importance of determining the complement of DNA is that it allows us to carry out the same process by which our body determines what bases to add to a coding strand based on a template strand. 

![The Central Dogma Path](https://github.com/jasonr-alex/genomics/blob/main/Central-dogma.jpeg)

Reverse transcription, however, conflicts against this central dogma; reverse transcription is an important viral replication process that allows for insertion in the host's genome. 

### Complementation and Reverse Complementation Code (Python)

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




 




