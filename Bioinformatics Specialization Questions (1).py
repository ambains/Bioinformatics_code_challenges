#!/usr/bin/env python
# coding: utf-8

# The below code challenges are from the Bioinformatics specialization from Coursera.  

# First problem presented is Vibrio cholerae, the pathogenic bacterium that causes cholera. Here is the nucleotide sequence for Vibrio cholerae. How can we find the "hidden message" in the ori region? We understand that the replication of a sequence is done with the help of DnaA and a protein known as DnaA box. DnaA box tells the DnaA protein to bind at the ori. How do we find the ori region in the sequence of Vibrio cholerae?

# Question in Coursera: Hidden Message Problem: Find a “hidden message” in the replication origin.
# 
# Input: A string Text (representing the replication origin of a genome).
# 
# Output: A hidden message in Text.
# 

# To present this as a computational problem, a function can be used to find the frequently occuring substring in the sequence because "if there are more occureses of the string, then it is more likely that binding will successfully occur". 

# In[1]:


def count (Text, Pattern):
        result = Text.count(Pattern)
        print(result)


# In[2]:


#lets test the function
Text = 'ACAACTATGCATACTATCGGGAACTATCCT'
Pattern = 'ACTAT'
count(Text, Pattern)


# ## Course Code Challenge 1

# Implement the count function for the dataset provided in the course. 
# 1. Download the dataset.
# 2. Import the dataset into Jupyter Notebook.
# 3. Read in the dataset. 

# Searching the dataframe did not work, so a for loop is needed to read the lines of the text file. 
# 
# The below code is provided by the Coursera instructors to understand the process of searching for a pattern in a nucleotide sequence. 

# In[3]:


def PatternCount(Text, Pattern):
    count = 0  # Initialize count variable
    for i in range(0, len(Text)- len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count + 1  # Use single '=' for assignment
    return count


# Pattern: CATGCTACA

# In[4]:


with open('dataset_2_6 (3).txt', 'r') as f:
    Text = f.read()


# In[5]:


Text


# In[6]:


Pattern ='CATGCTACA'
print(PatternCount(Text, Pattern))


# ## The Frequent Words Problem
# 
# Finding the most frequent k-mers in a sequence. 
# 
# Using the FrequentWords algorithim below is not the best way to calculate the number of repeated patterns in a sequence because of the runtime.
# 
# source for this for loop is the Bioinformatics Specialization by Coursera: 
# FrequentWords(Text, k)
#     FrequentPatterns ← an empty set
#     for i ← 0 to |Text| − k
#         Pattern ← the k-mer Text(i, k)
#         Count(i) ← PatternCount(Text, Pattern)
#     maxCount ← maximum value in array Count
#     for i ← 0 to |Text| − k
#         if Count(i) = maxCount
#             add Text(i, k) to FrequentPatterns
#     remove duplicates from FrequentPatterns
#     return FrequentPatterns
#     
#     
#     
# Instead using a key value pairs concept and creating a map or a dictionary is a better method to count the number of occurences of a k-mer in a sequence. 
# 
# For instance:
# 
# The following for loop will iterate through the text, assigning the frequency of each k-mer. The resulting map will consist of keys, each representing a specific k-mer sequence, and values corresponding to the number of times that particular k-mer occurs in the entire nucleotide sequence. 
# 
# source for this for loop is the Bioinformatics Specialization by Coursera:
# 
# FrequencyTable(Text, k)
#     freqMap ← empty map
#     n ← |Text|
#     for i ← 0 to n − k
#         Pattern ← Text(i, k)
#         if freqMap[Pattern] doesn't exist
#             freqMap[Pattern]← 1
#         else
#            freqMap[pattern] ←freqMap[pattern]+1 
#     return freqMap

# In[7]:


def FrequencyTable(Text, k):
    freqMAP = {}
    n = len(Text)
    for i in range(n - k + 1): # the i represents the index value. we are adding 1 because the end index is exclusive and it will not include the last k-mer value if the 1 is not added.
        Pattern = Text[i:i+k] # start and stop index 
        if Pattern not in freqMAP:
            freqMAP[Pattern] = 1
    else:
            freqMAP[Pattern] += 1
    return freqMAP


# In[8]:


Text = 'ACGTTTCACGTTTTACGG'
k = 3
result = FrequencyTable(Text, k)
print(result)


# In[9]:


text = "ATCGATCGAT"
k = 2

result = FrequencyTable(text, k)

print(result)


# In[10]:


# finding the k-mer with the highest number of occurences hence the max value.

def MaxMAP(freqMAP):
    max_value = max(freqMAP.values())
    
    return(max_value)


# The MaxMAP function is working on the map generated from the nucleotide sequence in line 36, and it is working on the values assigned to each key pair, which has the length of k-mer set by the researcher. 

# In[11]:


def MaxValue (Text, k):
    FrequentPatterns = {}
    freqMAP = FrequencyTable(Text,k)
    max_value = MaxMAP(freqMAP)
    for Pattern in freqMAP:
        if freqMAP[Pattern] == max_value:
            FrequentPatterns[Pattern] = freqMAP[Pattern]
    return FrequentPatterns
    


# In[12]:


Text = 'ACGTTTCACGTTTTACGG'
k = 3
freqMAP1 = FrequencyTable(Text, k)
print(freqMAP1)


# In[13]:


result = MaxMAP(freqMAP1)
result


# In[14]:


with open('dataset_2_13.txt', 'r') as f:
    Text = f.read()
    Text = Text.replace('\n','')
    print(Text)


# In[15]:


with open('dataset_2_13 (11).txt', 'r') as file:
    Text5 = file.read()
    Text5 = Text5.replace('\n','')


# In[16]:


Text5


# In[17]:


k5 = 12
MaxValue(Text5, k5)


# In[18]:


freqMAP = FrequencyTable(Text5, k5)


# In[19]:


freqMAP


# In[20]:


Example = "CGCCTCTAGTGTCACGCCTCTCGCCTCTCCTCATCCTCCTCATCCTCGCCTCTCGCGTGATGCGCCTCTCGCGTGATGCGCCTCTAGTGTCACGGACTCCAGTGTCACGGACTCCCGCGTGATGCGGACTCCCGGACTCCCGCGTGATGCGCGTGATGCGGACTCCCCTCATCCTAGTGTCACGCCTCTAGTGTCAAGTGTCACGCCTCTCGCGTGATGAGTGTCAAGTGTCACCTCATCCTCGGACTCCCCTCATCCTCCTCATCCTCGCGTGATGCGGACTCCCGCGTGATGCGCGTGATGAGTGTCAAGTGTCACCTCATCCTAGTGTCACGCCTCTCGCGTGATGAGTGTCACGGACTCCAGTGTCACCTCATCCTCGCCTCTCCTCATCCTCGCGTGATGAGTGTCACCTCATCCTCGGACTCCCGCGTGATGCGCGTGATGCGCGTGATGAGTGTCACCTCATCCTCCTCATCCTCGCCTCTCGCGTGATGCGGACTCCCGCCTCTCGCCTCTAGTGTCACGCCTCTAGTGTCACGGACTCCCGGACTCCAGTGTCACGCCTCTCGCCTCTCGCCTCTAGTGTCAAGTGTCACGGACTCCCGCGTGATGAGTGTCACGCCTCTCGGACTCCAGTGTCACCTCATCCTCGGACTCCCGGACTCCCGCCTCTCGCCTCTCCTCATCCTCGGACTCCCCTCATCCTCGGACTCCCGGACTCCCGCGTGATGAGTGTCACGCGTGATGAGTGTCACGCGTGATGCGCGTGATGCGCGTGATGCGCGTGATGCGCCTCTCCTCATCCTCGCCTCTCGGACTCCCGGACTCCCGGACTCCAGTGTCACGCGTGATGCGCGTGATGCGCGTGATGCGCCTCTAGTGTCA"
kmer = 11
FrequencyTable(Example, kmer)
freqMAP = FrequencyTable(Example, kmer)
MaxMAP(freqMAP)
MaxValue(Example,kmer)


# In[21]:


def find_most_frequent_kmers(text, k):
    kmers = {}
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        if kmer not in kmers:
            kmers[kmer] = 1
        else:
            kmers[kmer] += 1
    max_freq = max(kmers.values())
    return [kmer for kmer, freq in kmers.items() if freq == max_freq]

text = "CGCCTCTAGTGTCACGCCTCTCGCCTCTCCTCATCCTCCTCATCCTCGCCTCTCGCGTGATGCGCCTCTCGCGTGATGCGCCTCTAGTGTCACGGACTCCAGTGTCACGGACTCCCGCGTGATGCGGACTCCCGGACTCCCGCGTGATGCGCGTGATGCGGACTCCCCTCATCCTAGTGTCACGCCTCTAGTGTCAAGTGTCACGCCTCTCGCGTGATGAGTGTCAAGTGTCACCTCATCCTCGGACTCCCCTCATCCTCCTCATCCTCGCGTGATGCGGACTCCCGCGTGATGCGCGTGATGAGTGTCAAGTGTCACCTCATCCTAGTGTCACGCCTCTCGCGTGATGAGTGTCACGGACTCCAGTGTCACCTCATCCTCGCCTCTCCTCATCCTCGCGTGATGAGTGTCACCTCATCCTCGGACTCCCGCGTGATGCGCGTGATGCGCGTGATGAGTGTCACCTCATCCTCCTCATCCTCGCCTCTCGCGTGATGCGGACTCCCGCCTCTCGCCTCTAGTGTCACGCCTCTAGTGTCACGGACTCCCGGACTCCAGTGTCACGCCTCTCGCCTCTCGCCTCTAGTGTCAAGTGTCACGGACTCCCGCGTGATGAGTGTCACGCCTCTCGGACTCCAGTGTCACCTCATCCTCGGACTCCCGGACTCCCGCCTCTCGCCTCTCCTCATCCTCGGACTCCCCTCATCCTCGGACTCCCGGACTCCCGCGTGATGAGTGTCACGCGTGATGAGTGTCACGCGTGATGCGCGTGATGCGCGTGATGCGCGTGATGCGCCTCTCCTCATCCTCGCCTCTCGGACTCCCGGACTCCCGGACTCCAGTGTCACGCGTGATGCGCGTGATGCGCGTGATGCGCCTCTAGTGTCA"
k = 11
print(find_most_frequent_kmers(text, k))
freqMAP = FrequencyTable(Example, kmer)
MaxMAP(freqMAP)


# ### Frequent Words in Vibrio cholerae
# 
# This is the sequence for Vibrio cholerae, bacteria. 
# 
# "atcaatgatcaacgtaagcttctaagcatgatcaaggtgctcacacagtttatccacaac
# ctgagtggatgacatcaagataggtcgttgtatctccttcctctcgtactctcatgacca
# cggaaagatgatcaagagaggatgatttcttggccatatcgcaatgaatacttgtgactt
# gtgcttccaattgacatcttcagcgccatattgcgctggccaaggtgacggagcgggatt
# acgaaagcatgatcatggctgttgttctgtttatcttgttttgactgagacttgttagga
# tagacggtttttcatcactgactagccaaagccttactctgcctgacatcgaccgtaaat
# tgataatgaatttacatgcttccgcgacgatttacctcttgatcatcgatccgattgaag
# atcttcaattgttaattctcttgcctcgactcatagccatgatgagctcttgatcatgtt
# tccttaaccctctattttttacggaagaatgatcaagctgctgctcttgatcatcgtttc"
# 
# 
# ![image.png](attachment:image.png)
# 
# Based on the counts given for each k-mer, do any of the counts seem suprisingly large?

# My first observation is that the k-mers 'atgatcaag' and 'cttgatcat' are possibly reverse complementary. Whether the counts for any of these k-mers seems suprisingly large is a bit challneging to tell without knowing their role in regulating gene expression. 

# ### Reading complementary strand and the template strands of DNA.
# 
# ![image.png](attachment:image.png)
# 
# 
# 
# The complementary strand of "ACTATGCGACT" will be synthesized from 5 prime to 3 prime because each DNA strand is alwys synthesized from 5' → 3' direction.  
# 
# Understanding the directionality of DNA sequences will help me to progress in molecular biology and genomics in the following areas:
# 
# Sequence Directionality- the way DNA and RNA is synthesized is crucial to analyzing PCR reactions. 
# Complementarity- Understanding the pairing rules is fundamental to figuring out the opposite pair of the sequence.
# Transcription and Translation: DNA transcibed to RNA, and then the RNA translation to proteins.
# Mutation Analysis: Looking for changes in the DNA sequences such as insertions, deletions, or rearrangements, and understanding their impact on genetic diseases. 

# ### Reverse Complement Problem
# 
# Where the input is a DNA String = Pattern, and the reverse complementary string = Pattern<sub>rc</sub>
# 
# 
# In a DNA sequence, the two strands run in opposite directions, and they are complementary to each other. When searching for a specific sequence in the DNA, it's important to consider both the sequence itself and its reverse complement because the sequence could originate from either DNA strand.
# 
# The process of finding the reverse complement involves reversing the sequence and then replacing each nucleotide with its complement. This allows for comprehensive searching and identification of the sequence in the DNA, regardless of the strand it comes from.
# 
# 
# 

# In[22]:


def reverse_complement(dna_sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complement_sequence = [complement_dict[nucleotide] for nucleotide in dna_sequence]
    return "".join(complement_sequence[::-1])  # Reverses and joins the sequence
                  


# In[23]:


with open('dna_sequences.txt', 'r') as file:
    dna_sequences = file.readlines()

reverse_complements = [reverse_complement(seq.strip()) for seq in dna_sequences]

# If you want to write the reverse complements to a file:
with open('reverse_complements.txt', 'w') as file:
    for rc in reverse_complements:
        file.write(rc + '\n')


# Per the lecture, the 9-mer of ATGATCAAG and CTTGATCAT are reverse complement of each other, and occur six times in total in the below sequence of vibrio cholerae. For a DNA string of length 500 is suprising and leads to the hypothesis that the 9-mers are the DnaA boxes in Vibrio cholerae and are the origin sites for replication. 

# ### Pattern Matching Problem
# 
# 'However, before concluding that we have found the DnaA box of Vibrio cholerae, the careful bioinformatician should check if there are other short regions in the Vibrio cholerae genome exhibiting multiple occurrences of ATGATCAAG (or CTTGATCAT). After all, maybe these strings occur as repeats throughout the entire Vibrio cholerae genome, rather than just in the ori region. To this end, we need to solve the following problem.'- Instructors of Bioinformatics Specialization on Coursera. 
# 
# Pattern Matching Problem: Find all occurrences of a pattern in a string.
# 
# Input: Strings Pattern and Genome.
# 
# Output: All starting positions in Genome where Pattern appears as a substring.

# #### Why is it important to know the origin of replication?
# 
# In case researchers need to clone the DNA, and manipulate the DNA. Another reason is to understand the duplication of genetic material during cell division. 

# In[24]:


# for loop must have the output as an  index position as a list
def pattern_matching(pattern, genome):
    positions = []
    pattern_length = len(pattern)
    genome_length = len(genome)

    for i in range(genome_length - pattern_length + 1):
        if genome[i:i + pattern_length] == pattern:
            positions.append(i)

    return positions


# In[25]:


pattern = 'CCAGAGCCC'
with open ('dataset_3_5.txt', 'r') as file:
    genome = file.read()
    print(genome)


# In[26]:


print(pattern_matching(pattern, genome))


# In[27]:


pattern_example_2 = 'ATTAGTGAT'
with open ('dataset_3_5_example2.txt', 'r') as file:
    genome_example_2 = file.read()
    print(genome_example_2)


# In[28]:


sequence = pattern_matching(pattern_example_2, genome_example_2)
sequence


# In[29]:


sequence_without_commas = " ".join(str(element) for element in sequence)

print(sequence_without_commas)


# In[30]:


with open ('Vibrio_cholerae.txt', 'r') as file:
    vibrio_sequence = file.read()


# In[31]:


pattern_vibrio = 'CTTGATCAT' 
positions_vibrio = pattern_matching(pattern_vibrio, vibrio_sequence)


# In[32]:


sequence_without_commas = " ".join(str(element) for element in positions_vibrio)

print(sequence_without_commas)


# ### Looking for hidden messages in multiple genomes
# 
# The following is the ori region of Thermotoga petrophila, a bacterium that thrives in extremely hot environments. This bacterium does not contain even a single occurence of ATGATCAAG or CTTGATCAT which were present in the ori region of Vibrio bacterium. Therefore, different bacteteria might use different DnaA boxes for replication. 

# ### The Clump Finding Problem
# 
# Aiming to find  a pattern (k) is repeated in a (L) length of a nucleotide sequence that appears several times (t) in a short succession. Tihs scenario applies to attempting to find the ori in a newly sequenced bacterial genome. 
# 
# The close proximity of these patterns indicates a region where the replication process is likely to be initiated.

# In[33]:


def clump_pattern(Text, k, L, t): # finding a pattern that appears in a genome sequence several times in a short time 
    Patterns = [] # generating an array so I can add the frequency of the patterns to the array
    n = len(Text) # length of the genome that I am analyzing
    for i in range(n - L + 1):  # FOR each window of size L
        window = Text[i:i + L]  # Extract the window of length L from the genome.
        freqMAP = FrequencyTable(window, k) # previously defined Frequency table in line 7 will help to generate a key-value pair for the pattern of k length
        for s in freqMAP:  # for every key in the frequency map, if the keys in the frequency map occur more than once, then add them to the list of the array
            if freqMAP[s] >= t:
                Patterns.append(s) 
    Patterns = list(set(Patterns)) # removing duplicates from the list of patterns using set
    return Patterns # returns an array for the patterns that occur several times in the set length of genome that under analysis 


# In[34]:


with open ('dataset_4_5 (1).txt', 'r') as file:
    clumps_analysis = file.read()
    print(clumps_analysis)


# In[35]:


k_clumps = 8
L_clumps = 30
t_clumps = 3
result = clump_pattern(clumps_analysis, k_clumps, L_clumps, t_clumps)


# In[36]:


print(result)


# Result: I am getting an empty array meaning there are no patterns that satisy the condition of a pattern with length of 8 letters repeating 3 times in a sliding window length of the genome with 30 letters. 

# In[37]:


with open('dataset_4_5 (2).txt', 'r') as file_example2:
    clumps_analysis2 = file_example2.read()
    print(clumps_analysis2)


# In[38]:


k_clumps2 = 9
L_clumps2 = 24
t_clumps2 = 4
result_2 = clump_pattern(clumps_analysis2, k_clumps2, L_clumps2, t_clumps2)


# In[39]:


print(result_2) 


# Result is again no patterns match the conditions provided in the Text sequence.

# The below code was obtained from one the users of the Bioinformatics Specialization Course because the clump_pattern code had a long run time and was not efficient in retrieving the clump patterns in the E_coli text file. 
# 
# While studying both codes: I observed the following similarities and differences:
# 
# The main similarity of the codes clump_pattern and find_clumps:
#             
#     Same parameters(genome, k, L, t)
#     
#     
# Many differences in the codes clump_pattern and find_clumps, hence certain implementation details differ:
# 
#     Data structure for the output is array in clump_pattern is an array, but in find_clumps is a set. Using a set 
#     does have an advantage of not containing the duplicates. 
#     
#     In find_clumps, the code is starting to collect k-mers in the first window, and then updating the count of the 
#     k- mers as the (L) sliding window changes. Whereas, in clump_pattern the k-mers and the count of each k-mer is         being collected through the FrequencyTable function, and then the freqMAP is being adjusted per the threshold (t). 

# In[40]:


def find_clumps(genome, k, L, t):
    clumps = set()

    # Count occurrences of k-mers within the first window
    window = genome[:L]
    kmer_counts = {}
    for i in range(L - k + 1):
        kmer = window[i:i+k]
        kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1
        if kmer_counts[kmer] >= t:
            clumps.add(kmer)

    # Slide the window and update k-mer counts
    for i in range(1, len(genome) - L + 1):
        left_kmer = window[:k]
        kmer_counts[left_kmer] -= 1

        new_kmer = genome[i+L-k:i+L]
        kmer_counts[new_kmer] = kmer_counts.get(new_kmer, 0) + 1
        if kmer_counts[new_kmer] >= t:
            clumps.add(new_kmer)

        window = genome[i:i+L]

    return clumps

# Read genome from file
with open('E_coli.txt', 'r') as file_ecoli:
    genome = file_ecoli.read()

k = 9
L = 500
t = 3

clumps = find_clumps(genome, k, L, t)
num_clumps = len(clumps)

print("Number of different 9-mers forming (500, 3)-clumps:", num_clumps)


# ### Peculiar Statistics of the Forward and Reverse Half-Strands
# 
# Scenario: The origin of replication is unknown in a circular genome, so it is better to linearize it. In order to linearize it, we can use a skew diagram. 
# 
# Skew is equal to the number of G's in the first i nucleotides - number of C's in the first i nucleotides. Calculating the Skewi(Genome) can give intersting insights on the sequence, including possible locations of the orign of replication in some organisms' genomes. 
# 
# ### Problem: Give all values of Skewi (GAGCCACCGCGATA) for i ranging from 0 to 14.

# In[41]:


genome = 'GAGCCACCGCGATA'
def skewi(genome):
    skew = 0
    print(skew, end=" ")
    
    for nucleotide in genome:
        if nucleotide == "C":
            skew -= 1
        elif nucleotide == "G":
            skew += 1
            
        print(skew,end =" ")
        
skewi(genome)


# In[42]:


string = [skewi(genome)]


# In[43]:


print(type(string))


# In[44]:


def minimumval(string):
    min_value = min (string)
    index_list = [index for index in range(len(string)) if string[index] == min_value]
    return index_list


# In[45]:


result_example1 = minimumval(string)
print(result_example1)


# In[46]:


with open('dataset_7_10.txt', 'r') as skew_file:
    minimum_skew_problem = skew_file.read()


# In[47]:


string_example2 = [skewi(minimum_skew_problem)]


# In[48]:


minimumval(string_example2)


# In[49]:


with open('dataset_7_10 (2).txt', 'r') as skew_file2:
    minimum_skew_problem2 = skew_file2.read()


# In[50]:


def skewi_test(genome):
    skew = 0
    skew_values = [0]  # list to store skew values

    for nucleotide in genome:
        if nucleotide == "C":
            skew -= 1
        elif nucleotide == "G":
            skew += 1
        skew_values.append(skew)

    min_skew = min(skew_values)  # find minimum skew value

    min_positions = [i for i, value in enumerate(skew_values) if value == min_skew]  # find positions of minimum skew

    return min_positions

genome = 'GAGCCACCGCGATA'
min_skew_positions = skewi_test(genome)
print(min_skew_positions)


# In[51]:


skewi_test(minimum_skew_problem2)


# In[52]:


genome1 = 'CATTCCAGTACTTCGATGATGGCGTGAAGA'
print(skewi_test(genome1))


# In[53]:


genome2 = 'CATTCCAGTACTTCATGATGGCGTGAAGA'
print(skewi_test(genome1))


# In[54]:


with open('dataset_7_10 (3).txt', 'r') as skew_file3:
    minimum_skew_problem3 = skew_file3.read()


# In[55]:


skewi_test(minimum_skew_problem3)


# Solving the minimum skew problem provides us with an estimated location of the origin at 89793-89795. 

# From the previous functions and lessons of this course, I have learned that the origin of replication for a DNA sequence of an organism if 500 nucleotides. If the origin of replication is known in an organism's genome, then the finding the most occuring k-mers can help in determining the initiation site for replication, where the DnaA protein will bind to the DnaA Box. However, if the origin of replication is not know, then the skew measurement can help in finding the estimated location of the origin of replication. 
# 
# Keeping in mind that the DnaA protein can bind not only to "perfect" DnaA Box but to variations of it. 

# Moving to figuring out the HammingDistance which is the number of mismatches in two equal length of strings.

# ### Hamming Distance Problem:
# 
# Input will be two equal lengths of strings.
# Output will be the number of mismatches in the nucleotides between these two strings of sequences. 

# In[56]:


with open('stringone', 'r') as file:
    p = file.read()


# In[57]:


with open('stringtwo', 'r') as file:
    q = file.read()


# In[58]:


dist = 0
for i in range(len(p)):
     if p[i] != q[i]:
        dist += 1
print(dist)


# In[59]:


def HammingDistance(Pattern, Genome, d):
    positions = []
    freqMAP = FrequencyTable(Pattern, Genome)
    index_list = [index for index in freqMAP if count <= d]
    return index_list


# In[60]:


# revised the code below where p and q are assumed to be equal lengths and we are measuring the number of mismatches in these 
# two strings

def HammingDistance(p,q):
    dist = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            dist +=1
    return dist


# In[61]:


p = 'CAGAAAGGAAGGTCCCCATACACCGACGCACCAGTTTA'
q = 'CACGCCGTATGCATAAACGAGCCGCACGAACCAGAGAG'
HammingDistance(p,q)


# The issues with the above function:
# 
# * Count is not defined so it will raise an error.
# * FrequencyTable returns the number of times a pattern occurs in a genome not the position of the pattern and the number of mismatches for each position.
# * HammingDistance is the number of mismatches between equal lengths of strings, but the function above has string one has a pattern and string 2 has the genome, which are clearly of two different lengths.

# In[62]:


p = 'CTACAGCAATACGATCATATGGGATCCGAGTGGCCGTAGACACACGT'
q = 'CTACCCCGCTGCTCAATGACCGGGACTAAAGAGGCGAAGATTATGGTGTG'

print(HammingDistance(p,q))


# In[63]:


def approximate_pattern_matching(pattern, genome, d):
    for i in range (len(genome) - len(pattern) +1): # scan the genome and find the pattern xyz
        if  HammingDistance(pattern, genome[i:i + len(pattern)]) <= int(d):
            positions.append(i) # add the positions of the pattern to the positions list.
    return positions


# In[64]:


with open ('dataset_9_4 (3).txt', 'r') as file:
    genome = file.read()


# In[65]:


pattern = 'GTGATCGATC'
d = 6
positions = []
output = approximate_pattern_matching(pattern, genome, d)


# In[66]:


output_without_commas = str(output).replace(",", "")


# In[67]:


output_without_commas


# In[68]:


with open('dataset_9_4 (4).txt', 'r') as example_file:
    genome2 = example_file.read()


# In[69]:


pattern2 = 'TCACAGTGAGG'
d = 4
positions = []
output2 = approximate_pattern_matching(pattern2, genome2, d)
output_without_commas_example_2 = str(output2).replace(",", "")


# In[70]:


output_without_commas_example_2


# In[71]:


pattern = 'TCACCAACGCA'
with open('dataset_9_4 (9).txt', 'r') as file:
    genome = file.read()
d = 4
positions = []

print(*approximate_pattern_matching(pattern, genome, d))


# ### Compute Count2(AACAAGCTGATAAACATTTAAAGAG, AAAAA)
# 
# Scenario: Find the total occurences of k-mer in the genomic sequence with at most d mismatches.

# In[72]:


# defining a function of approximate_count_mismatches
                                    
def approx_count_mismatches(genome, pattern):
    d = []
    dist = 0
    for i in range(len(genome) - len(pattern) + 1): 
        dist = HammingDistance(pattern, genome[i:i + len(pattern)])
        d.append(dist)

    return d                           


# In[73]:


genome = 'AACAAGCTGATAAACATTTAAAGAG'
pattern = 'AAAAA'

print(max(*approx_count_mismatches(genome, pattern))) # the maximum number of mismatches


# In[74]:


def approx_count_mismatches_variation(genome, pattern, d):
    occurences = []
    count = 0
    for i in range(len(genome) - len(pattern) + 1): 
        count = HammingDistance(pattern, genome[i:i + len(pattern)]) <= d
        count += 1
        occurences.append(count) # the number of mismatches per position 

    return occurences         


# The next step is to find the occurences of the patterns with at most two mismatches against the frequent k-mer (pattern) in the genome sequence. 

# In[75]:


# HammingDistance to find the number of mismatches 

genome = 'AACAAGCTGATAAACATTTAAAGAG'
pattern = 'AAAAA'
d = 2

print(*approx_count_mismatches_variation(genome, pattern, d))


# This is returning the number of mismatches per position, not the total count of the mismatches that meet the criteria of >=2 (d). 

# In[76]:


def approx_count_mismatches_variation2(genome, pattern, d):
    count = 0
    for i in range(len(genome) - len(pattern) + 1): 
        if HammingDistance(pattern, genome[i:i + len(pattern)]) <= d:
            count += 1

    return count # the total number of k-mers with at most two mismatches in the given genomic sequence. 


# In[77]:


genome = 'AACAAGCTGATAAACATTTAAAGAG'
pattern = 'AAAAA'
d = 2

print(approx_count_mismatches_variation2(genome, pattern, d))


# In[78]:


with open('dataset_9_6.txt', 'r') as file:
    genome = file.read()

pattern = 'TTATGC'
d = 2

print(approx_count_mismatches_variation2(genome, pattern, d))


# In[79]:


genome1 = 'TACGCATTACAAAGCACA'
pattern = 'AA'
d = 1
print(approx_count_mismatches_variation2(genome1, pattern, d))


# In[80]:


genome2 = 'CATGCCATTCGCATTGTCCCAGTGA'
pattern = 'CCC'
d = 2
print(approx_count_mismatches_variation2(genome2, pattern, d))


# A most frequent k-mer with at most d matches is simply a string pattern maximizing the approx_count_mismatches_variation2 function. 

# ### Frequent words (k-mers) with d mismatches also known as d-neighborhood
# 
# The goal is to create a d-neighborhood Neighbors(Pattern, d), which is the set of all k-mers who Hamming distance from Pattern is less than or equal to the given d (integer) value. 
# 
# The pseudocodes for the Neighbors and FrequentWordswithmismatches are provided in the course. 

# In[81]:


def Neighbors(Pattern, d):
    if d == 0: # we start with checking if we are even accepting any mismatches
        return {Pattern} # if d = 0 meaning no mismatches accepted, then return the pattern (k-mer)
    if len(Pattern) == 1:
        return {'A', 'C', 'G', 'T'} # returns all possiblilities for as neighbors when the length of the pattern is 1.

    Neighborhood = set() # creating an empty set to store the variations of the pattern with d mismatches.
    SuffixNeighbors = Neighbors(Pattern[1:], d) # checking what are the neighbors of the pattern except for the first character with allowance for d mismatches
    for Text in SuffixNeighbors: # checking that the difference between the identified variation and the length of the k-mer is less than the given d value.
        if HammingDistance(Pattern[1:], Text) < d:
            for x in ['A', 'C', 'G', 'T']:
                Neighborhood.add(x + Text) # if the difference between the identified variation and the length of the k-mer meets 
                                           # the condition of d mismatches then the variation is added to the neighborhood set.
        else:
            Neighborhood.add(Pattern[0] + Text) # if the number of mismatches are greater than d then it adds the first character of 
                                                # the pattern with the text and adds it to the neighborhood set
    return Neighborhood


# In[82]:


pattern = 'CCAGTCAATG'
d = 1
results = Neighbors(pattern, d)


# In[83]:


def convert(set):
    return list(set)


# In[84]:


len(convert(results))


# In[85]:


pattern = 'ACGT'
d = 3
results = Neighbors(pattern, d)


# In[86]:


len(convert(results))


# In[87]:


def FrequentWordswithmismatches(Text, k, d):
    Patterns = []
    freqMAP = {}
    n = len(Text)
    for i in range(n - k + 1):
        Pattern = Text[i:i+k] 
        neighborhood = Neighbors(Pattern, d) # the collection of all k-mers (most frequent occuring pattern) is called d-neighborhood
        for neighbor in neighborhood: 
            if neighbor not in freqMAP:  
                freqMAP[neighbor] = 1
            else:
                freqMAP[neighbor] += 1 
            
    m = max(freqMAP.values())
    for pattern in freqMAP:  
        if freqMAP[pattern] == m:
            Patterns.append(pattern)

    return Patterns


# In[88]:


with open('dataset_9_9.txt', 'r') as file:
    genome = file.read()


# In[89]:


print(FrequentWordswithmismatches(genome, 7, 3))


# In[90]:


with open('dataset_9_9 (1).txt', 'r') as file:
    genome = file.read()


# In[91]:


print(FrequentWordswithmismatches(genome, 7, 3))


# ### Frequent Words with Mismatches and Reverse Complements Problem: Find the most frequent k-mers (with mismatches and reverse complements) in a string.
# 
# Input: A DNA string Text as well as integers k and d.
# 
# Output: All k-mers Pattern maximizing the sum Countd(Text, Pattern)+ Countd(Text, Patternrc) over all possible k-mers.

# In[92]:


with open('dataset_9_10.txt', 'r') as file:
    genome = file.read()


# In[93]:


print(FrequentWordswithmismatches(genome, 6,2))


# $4^k$ in studying genomics is used to calculate the total number of possibilities in patterns of a particular length. $4^k$ can be used in the function to determine the inputs that will maximize the sum of the count of the pattern with k length and d mismatches and the count of the reverse complement pattern with k length and d mismatches. 
# 
# 1. Figure out the possible k-mers.
# 2. Need to define a function to generate the reverse complement of the k-mers in the genomic sequence. 
# 2. Count the number of times each possible kmer and its reverse complement appear within Text.
# 
# 
# Definition of reverse complement (rc) pair: To find the rc of a sequence, first reverse the sequence, then find its complement nucleotide pair. 

# In[94]:


def rev_comp(pattern):
    rpattern=''
    a=[]
    for i in pattern:
        if i=='A':
            a.append('T')
        if i=='T':
            a.append('A')
        if i=='C':
            a.append('G')
        if i=='G':
            a.append('C')
    rpattern =rpattern.join(a)
    return(rpattern[::-1])


# In[95]:


# Count the number of times each possible kmer and its reverse complement 
# appear within Text 
# FrequencyTable function can be used which is previously defined in this notebook.
# FrequencyTable will need to be modified to include the reverse complement 

def FrequencyTable(Text, k):
    freqMAP = {}
    n = len(Text)

    def ReverseComplement(Pattern):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return "".join(complement[base] for base in reversed(Pattern))

    for i in range(n - k + 1): 
        Pattern = Text[i:i+k]
        RevPattern = ReverseComplement(Pattern) # calculate the reverse complement
        for seq in [Pattern, RevPattern]:
            if seq not in freqMAP:
                freqMAP[seq] = 1
            else:
                freqMAP[seq] += 1
    return freqMAP


# 1. Return the k-mers and its associated reverse complement sequence. 

# In[96]:


Text = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
k = 4
result = FrequencyTable(Text, k)
print(result)


# 2. Find the k-mers and its associated reverse complement sequence per d number of mismatches. I can use FrequentWordswithmismatches function that I defined earlier in the notebook.
# 
# ### The below code was generated through the pseudocode given by the instructors in the lessons, and the fellow students in the comments of this course. 

# In[97]:


def FrequentWordswithmismatches2(Text, k, d):
    Patterns = []
    freqMAP = {}
    n = len(Text)

    def ReverseComplement(Pattern):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return "".join(complement[base] for base in reversed(Pattern))

    for i in range(n - k + 1): 
        Pattern = Text[i:i+k]
        RevPattern = ReverseComplement(Pattern) # calculate the reverse complement
        for seq in [Pattern, RevPattern]:
            if seq not in freqMAP:
                freqMAP[seq] = 1
            else:
                freqMAP[seq] += 1

    for i in range(n - k + 1):
        Pattern = Text[i:i+k] 
        RevPattern = ReverseComplement(Pattern)
        for seq in [Pattern, RevPattern]:
            neighborhood = Neighbors(seq, d) # the collection of all k-mers (most frequent occuring pattern) is called d-neighborhood
        for neighbor in neighborhood: 
            if neighbor not in freqMAP:  
                freqMAP[neighbor] = 1
            else:
                freqMAP[neighbor] += 1 
    
    totalFreqMAP = {}
    for kmer, freq in freqMAP.items():
        rev_kmer = ReverseComplement(kmer)
        total_freq = freq
        if rev_kmer in freqMAP:
            total_freq += freqMAP[rev_kmer]
        totalFreqMAP[kmer] = total_freq
            
    m = max(totalFreqMAP.values())
    m_kmers = [kmer for kmer, freq in totalFreqMAP.items() if freq == m]

    return m_kmers


# In[98]:


with open('dataset_9_10 (2).txt', 'r') as file:
    Text = file.read()


# In[99]:


k = 5
d = 3
print(FrequentWordswithmismatches2(Text, k, d))


# ### Motif Finding algorithm 
# pseudocode provided by coursera instructors
# 
# The goal is find the pattern with d mismatches in the several Dna strings. This differs from the frequentwords problem because that is one sequence, whereas now we are analyzing different sequences at once. 

# In[100]:


def MotifEnumeration(dna, k, d):
    patterns=[]
    for i in range (0,len(dna[0])-k+1):
        neighbors= Neighbors(dna[0][i:i+k],d)
        for j in neighbors:
            count=0
            for l in dna:
                for i in range(0,len(l)-k+1):
                    if HammingDistance(j, l[i:i+k])<=d:
                        count+=1
                        break
            if count==len(dna):
                patterns.append(j)
    Patterns = [] 
    [Patterns.append(x) for x in patterns if x not in Patterns] 
    Result = ""
    for item in Patterns:
        Result = Result + " " + str(item)
    return Result
        
        


# In[101]:


dna = ('ATTTGGC', 'TGCCTTA', 'CGGTATC' ,'GAAAATT')
k = 3
d = 1
result = print(MotifEnumeration(dna, k,d))


# In[117]:


with open('dataset_156_8 (1).txt', 'r') as file:
    dna = file.read()
    dna = dna.replace('\n','')


# In[118]:


type(dna)


# In[116]:


dna = ('CTGAGGTTAAAAGCTCTTAACGTGT', 'CATTAACAACCTTGGCATAGTGTCG', 'CTGCGTAAGAACAGGGGCCCCGTGT', 'TCTTGAGAAGTATGGCCTGGGCTTC' ,'TAGTACTTGAGCCGCAAAGTAGAGC', 'CGTGGGCGATCTATCTGCTGCAGGT')
k = 5
d = 2
result = print(MotifEnumeration(dna, k ,d))


# In[121]:


with open ('dataset_156_8 (3).txt', 'r') as file:
    dna = file.read()
    dna = dna.replace("\n", "")


# In[122]:


dna


# In[123]:


type(dna)


# In[125]:


dna = ('TCTACAATGTTTGTCGCCCAAGGAT', 'TTACGACTTGCTCTGTTGTCGTGTC', 'TAATGTTGCACTTAGACCCATTACC', 'GTCTACTCTGTAACTCTTGCTTGTC', 'AGTATTGTAAATTTTCCATATTGCA', 'CCACGGACCTCTGTGTTGCGATCTG')
k = 5
d = 2
result = print(MotifEnumeration(dna, k ,d))

