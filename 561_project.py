#!/usr/bin/env python
# coding: utf-8

# In[2]:


from __future__ import division, print_function
import time
from contextlib import suppress
import re
import sys
import csv   
import numpy as np
import pandas as pd
import operator
from collections import Counter
from math import ceil
from math import floor
from math import log
from math import exp
import random
import timeit
from IPython.display import Audio, display


scaling_factor = 2
score_SM = 5
delta = 5
word_size = int(sys.argv[2])
gap = -1


# In[3]:


pt ={'match': 1, 'mismatch': -1, 'gap': -1}

def mch(alpha, beta):
    if alpha == beta:
        return pt['match']
    elif alpha == '-' or beta == '-':
        return pt['gap']
    else:
        return pt['mismatch']

def NW(s1, s2):
    m, n = len(s1), len(s2)
    score = np.zeros((m+1, n+1))
    
    #Initialization
    for i in range(m+1):
        score[i][0] = pt['gap'] * i
    for j in range(n+1):
        score[0][j] = pt['gap'] * j
    
    #Fill
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diag = score[i-1][j-1] + mch(s1[i-1], s2[j-1])
            delete = score[i-1][j] + pt['gap']
            insert = score[i][j-1] + pt['gap']
            score[i][j] = max(diag, delete, insert)

    #print('score matrix = \n%s\n' % score)
    align1, align2 = '', ''
    i,j = m,n
    
    #Traceback
    while i > 0 and j > 0:
        score_current = score[i][j]
        score_diag = score[i-1][j-1]
        score_left = score[i][j-1]
        score_up = score[i-1][j]
        
        #print('score_current: ',score_current)
        #print('score_diag: ',score_diag)
        #print('score_left: ',score_left)
        #print('score_up: ',score_up)

        if score_current == score_diag + mch(s1[i-1], s2[j-1]):
            #print('diag')
            a1,a2 = s1[i-1],s2[j-1]
            i,j = i-1,j-1
        elif score_current == score_up + pt['gap']:
            #print('up')
            a1,a2 = s1[i-1],'-'
            i -= 1
        elif score_current == score_left + pt['gap']:
            #print('left')
            a1,a2 = '-',s2[j-1]
            j -= 1
        #print('%s ---> a1 = %s\t a2 = %s\n' % ('Add',a1,a2))
        align1 += a1
        align2 += a2
            

    while i > 0:
        a1,a2 = s1[i-1],'-'
        #print('%s ---> a1 = %s\t a2 = %s\n' % ('Add',a1,a2))
        align1 += a1
        align2 += a2
        i -= 1
        
    while j > 0:
        a1,a2 = '-',s2[j-1]
        #print('%s --> a1 = %s\t a2 = %s\n' % ('Add',a1,a2))
        align1 += a1
        align2 += a2
        j -= 1
    
    align1 = align1[::-1]
    align2 = align2[::-1]
    seqN = len(align1)
    sym = ''
    seq_score = 0
    ident = 0
    for i in range(seqN):
        a1 = align1[i]
        a2 = align2[i]
        if a1 == a2:
            sym += '|'
            ident += 1
            seq_score += mch(a1, a2)
    
        else: 
            seq_score += mch(a1, a2)
            sym += ' '
        
    ident = ident/seqN * 100
    
    #print('Identity = %2.1f percent' % ident)
    #print('Score = %d\n'% seq_score)
    #print(align1)
    #print(sym)
    #print(align2)
    return align1, align2


# In[4]:


def read_files():
    global confList
    global sequence
    global seq_arr
    fa_path = 'chr22.maf.ancestors.42000000.complete.boreo.fa'
    conf_path = 'chr22.maf.ancestors.42000000.complete.boreo.conf'


    with open(fa_path) as f:
        sequence = f.read()
    seq_arr = [c for c in sequence]
    with open(conf_path) as c:
        conf = c.read()
    confList = conf.split(" ")
    confList.pop()


# In[5]:


def make_matrix():
    mat = []
    i=0
    for prob in confList:
        #print(i)
        p_pred = float(prob)
        p_other = float((1-p_pred)/3)
        if p_other == 0.0:
            p_other = 10**-20
        if sequence[i] == 'A':
            mat.append({'A':p_pred, 'C':p_other, 'G':p_other,'T':p_other})
            #mat.append([p_pred, p_other, p_other, p_other])
        elif sequence[i] == 'C':
            mat.append({'A':p_other, 'C':p_pred, 'G':p_other,'T':p_other})
            #mat.append([p_other, p_pred, p_other, p_other])
        elif sequence[i] == 'G':
            mat.append({'A':p_other, 'C':p_other, 'G':p_pred,'T':p_other})    
            #mat.append([p_other, p_other, p_pred, p_other])
        else:
            mat.append({'A':p_other, 'C':p_other, 'G':p_other,'T':p_pred})
            #mat.append([p_other, p_other, p_other, p_pred])
        i+=1

    #print(np.matrix(mat))
    return mat



# In[6]:


def find_seeds(query,word_size):
    global seeds
    global query_seq
    query_seq=[]
    gene = []
    words_seq = []
    seeds = []
    for line in query:
        for words in line[0:len(line)]:
            query_seq.append(words)
    for i in range(0,len(query_seq)):
        if i+word_size-1 < len(query_seq):
            words_seq.append(query_seq[i:i+word_size])
        
    score = 0
    for i in range(0,len(words_seq)):
        k = 0
        score_init = 0
        index = 0
        while(len(matrix) >= k+word_size):
            score = 0
            log_score=0
            for j in range(k,k+word_size):
                nuc = words_seq[i][j-k]
                score += log(scaling_factor*matrix[j][nuc])
                #score += matrix[j][nuc]

                #print(score_exp, score)
            #print(score,score_init)
            if score > score_init:
                score_init = score
                index = k
                #print(index, score_init)
            k = k+1
        if score_init > 0.0:
            seeds.append([score_init,i,index])

    return seeds


# In[7]:


def ungapped_ext():
    HSPs = []
    for seed in seeds:
        #seed = seeds[10]
        #print(seed)
        #print(query_seq)
        seed_seq = query_seq[seed[1]:seed[1]+word_size]
        #print(seed_seq)
        left = query_seq[0:seed[1]]
        #print(left)
        right = query_seq[seed[1]+word_size:len(query_seq)]
        #print(right)

        #extend to the left
        score_left = seed[0]
        max_score_left = seed
        i=len(left)-1
        j=1
        #print(i)
        while(score_left >= max_score_left[0]-delta and i>=0 and (seed[2]-j)>=0): 
            score_left += log(scaling_factor*matrix[seed[2]-j][left[i]])
            #print(matrix[seed[2]-j][left[i]],log_score,score_left,max_score_left[0])
            if score_left>max_score_left[0]:
                max_score_left = [score_left,i,seed[2]-j]
                #print("new max")
            i-=1
            j+=1
        #print("max_score_left= ",max_score_left)

        #extend to the right
        score_right = max_score_left[0]
        log_score = 0
        max_score_right = [score_right,-1,-1]
        k=0
        l=0
        #print(i)
        while(score_right >= max_score_right[0]-delta and k<len(right) and (seed[2]+word_size+l)<len(matrix)):
            score_right += log(scaling_factor*matrix[seed[2]+word_size+l][right[k]])
            #print(matrix[seed[2]+word_size+l][right[k]], log_score, score_right,max_score_right[0])
            if score_right>max_score_right[0]:
                max_score_right = [score_right,k,l]
                #print("new max")
            k+=1
            l+=1
        #print("max_score_right= ",max_score_right)
        #hsp = {'Score': max_score_right[0], 
        #       'Query indices': [max_score_left[1],len(left)+word_size+max_score_right[1]], 
        #       'Matrix indices': [max_score_left[2],seed[2]+word_size+max_score_right[2]]}
        hsp = {'Score': max_score_right[0], 
               'Query start': max_score_left[1],
               'Query end':len(left)+word_size+max_score_right[1], 
               'Matrix start': max_score_left[2],
               'Matrix end': seed[2]+word_size+max_score_right[2]}

        HSPs.append(hsp)
        del score_left, score_right, max_score_left, max_score_right, i,j,k,l
    
    scores = [hsp['Score'] for hsp in HSPs]
    best_score = max(scores)
    for hsp in HSPs:
        if hsp['Score']==best_score:
            best_HSP = hsp
            break
    #print(best_HSP)
    #print(HSPs)
    
    HSP_arr = np.array(HSPs)
    seen = set()
    global new_HSPs
    new_HSPs = []
    for hsp in HSPs:
        t = tuple(hsp.items())
        if t not in seen:
            seen.add(t)
            new_HSPs.append(hsp)
            #print(hsp)
            query_str = ''.join(query_seq[hsp['Query start']:hsp['Query end']+1])
            #dna_str = ''.join(seq_arr[hsp['Matrix start']:hsp['Matrix end']+1])
            #print(query_str)
            #print(dna_str)
    HSP_seqs = []
    for hsp in new_HSPs:
        query_str = ''.join(query_seq[hsp['Query start']:hsp['Query end']+1])
        dna_str = ''.join(seq_arr[hsp['Matrix start']:hsp['Matrix end']+1])
        HSP_seqs.append([query_str,dna_str])
    return HSP_seqs


# In[8]:


def gapped_ext(gap):
    extensions = []
    for hsp in new_HSPs:
        with suppress(Exception):
        #print(hsp)        
            query_R = query_seq[hsp['Query end']:len(query_seq)]
            matrix_R = matrix[hsp['Matrix end']:hsp['Matrix end']+2*len(query_R)]
            #matrix_R = matrix[hsp['Matrix end']:len(matrix)]
            #print(query_R)
            seq_arr_R = seq_arr[hsp['Matrix end']:hsp['Matrix end']+2*len(query_R)]
            #print(seq_arr_R)
            #print(np.array(matrix_R))
            #CALL SMITH WATERMAN ALGORITHM:
            query_aln_R, sym_R, dna_aln_R,matrix_pos_R = SM(query_R,matrix_R,gap)
            matrix_end = hsp['Matrix end']+matrix_pos_R

            query_L = query_seq[0:hsp['Query start']+1]
            query_L_rev = query_L.copy()
            query_L_rev.reverse()
            #print(query_L)
            seq_arr_L = seq_arr[hsp['Matrix start']-2*len(query_L):hsp['Matrix start']+1].copy()
            #print(seq_arr_L)
            seq_arr_L.reverse()
            #print(seq_arr_L)
            #matrix_L = matrix[0:hsp['Matrix start']+1]
            matrix_L = matrix[hsp['Matrix start']-2*len(query_L):hsp['Matrix start']+1]
            matrix_L_rev = matrix_L.copy()
            matrix_L_rev.reverse()
            #print(query_L_rev, matrix_L_rev, gap)
            #CALL SMITH WATERMAN ALGORITHM:
            query_aln_L_rev, sym_L_rev, dna_aln_L_rev, matrix_pos_L = SM(query_L_rev, matrix_L_rev, gap)
            query_aln_L = reverse(query_aln_L_rev)
            sym_L = reverse(sym_L_rev)
            dna_aln_L = reverse(dna_aln_L_rev)
            matrix_start = hsp['Matrix start']-matrix_pos_L+1
            #print(np.array(matrix_L_rev))
            #print(query_aln_L)
            #print(sym_L)
            #print(dna_aln_L)

            #print(matrix_start,matrix_end)
            #print(seq_arr[matrix_start:matrix_end])
            result_align = {'Left query':query_aln_L, 
                         'Left dna': dna_aln_L,
                         'Right query':query_aln_R, 
                         'Right dna': dna_aln_R, 
                         'Matrix indices':[matrix_start,matrix_end]}
            extensions.append(result_align)
            #print(extensions)
    return extensions


# In[9]:


def reverse(s): 
    str = "" 
    for i in s: 
        str = i + str
    return str


def get_score(nuc, probs):
    #print(nuc,probs)
    score = log(score_SM*probs[nuc])
    return(score)

def SM(query, dna,gap):
    m, n = len(query), len(dna)
    H = np.zeros((m+1, n+1))    
    T = np.zeros((m+1, n+1))
    max_score = 0
    # Score, Pointer Matrix
    H[1][1]=10**10
    T[1][1]=3
    max_i, max_j = 1,1
    for i in range(2, m + 1):
        for j in range(2, n + 1):
            sc_diag = H[i-1][j-1] + get_score(query[i-1], dna[j-1])
            sc_up = H[i][j-1] + gap
            sc_left = H[i-1][j] + gap
            H[i][j] = max(0,sc_left, sc_up, sc_diag)
            
            if H[i][j] == 0: T[i][j] = 0
            if H[i][j] == sc_left: T[i][j] = 1
            if H[i][j] == sc_up: T[i][j] = 2
            if H[i][j] == sc_diag: T[i][j] = 3
            if H[i][j] >= max_score:
                max_i = i
                max_j = j
                max_score = H[i][j];

    
    #print('H=\n',H,'\n')
    #print('T=\n',T,'\n')
    align1, align2 = '', ''
    i,j = max_i,max_j
    matrix_pos = j
    #Traceback
    while T[i][j] != 0:
        dna_nuc = max(dna[j-1].items(), key=operator.itemgetter(1))[0]
        if T[i][j] == 3:
            a1 = query[i-1]
            a2 = dna_nuc
            i -= 1
            j -= 1
        elif T[i][j] == 2:
            a1 = '-'
            a2 = dna_nuc
            j -= 1
        elif T[i][j] == 1:
            a1 = query[i-1]
            a2 = '-'
            i -= 1
        #print('%s ---> a1 = %s\t a2 = %s\n' % ('Add',a1,a2))
        align1 += a1
        align2 += a2

    align1 = align1[::-1]
    align2 = align2[::-1]
    sym = ''
    iden = 0
    for i in range(len(align1)):
        a1 = align1[i]
        a2 = align2[i]
        if a1 == a2:                
            sym += '|'
            iden += 1
        elif a1 != a2 and a1 != '-' and a2 != '-': 
            sym += ' '
        elif a1 == '-' or a2 == '-':          
            sym += ' '
    #print(matrix_pos)
    identity = iden / len(align1) * 100
    #print('Identity = %f percent' % identity)
    #print('Score =', max_score)
    #print(align1)
    #print(sym)
    #print(align2)
    return align1, sym, align2, matrix_pos


# In[ ]:





# In[63]:


def generate_query(length,matrix):
    query=""
    len_i_d = 10 #max length of inserts and deletions
    start = random.randrange(len(matrix)-int(length)-1)
    end = start+length
    #print(start,end)
    arr_nuc = ['A','C','G','T']
    nb_ins = 0
    nb_del = 0
    nb_sub = 0
    insert = False
    delete = False
    substit = False
    i = start
    while i<end:
        if np.random.uniform() < 0.07:
            insert = True
        elif np.random.uniform() < 0.07:
            delete = True
        elif np.random.uniform() < 0.07:
            substit = True
            
        if insert == True:
            insert_length = random.randrange(1,len_i_d)
            #print("Insertion %d" %insert_length)
            nb_ins+=insert_length
            to_insert = np.random.choice(arr_nuc, insert_length, p=[0.25,0.25,0.25,0.25])
            query = query + ''.join(to_insert)
            insert = False
        if delete == True:
            delete_length = random.randrange(1,len_i_d)
            i += delete_length
            nb_del +=delete_length
            #print("Deletion %d" %delete_length)
            delete = False            
        d = matrix[i]
        new_nuc = np.random.choice(arr_nuc, 1, p=[d['A'], d['C'],d['G'], d['T']])
        if substit == True:
            old = new_nuc
            while old == new_nuc:
                new_nuc = np.random.choice(arr_nuc, 1, p=[d['C'], d['A'],d['T'], d['G']])
            nb_sub += 1
            #print("Substitution (%s to %s)" %(old, new_nuc))
            substit = False
        query = query+new_nuc[0]            
        i+=1
    #query ="CCTCAGAGCCGATAAATTCTACCGAGTT"
    #print(query, len(query))
    return query, start, end, nb_ins, nb_del, nb_sub

    
    
    
#random_query, s, e, nb_ins, nb_del = generate_query(50)


# In[113]:


if __name__ == '__main__':
    if len(sys.argv)!=3:
        print("This algorithm takes 2 arguments: (1) length of query sequence and (2) word size")
    t0 = time.process_time()
    query_length = int(sys.argv[1])
    read_files()
    #matrix = make_matrix()[:][0:5000]
    matrix = make_matrix()
    query, start, end, nb_ins, nb_del, nb_sub  = generate_query(query_length,matrix)
    print("Generated query: %s" %query)
    print("   Original DNA: %s" %sequence[start:end])
    print("# insertions = %d, # deletions = %d, # substitutions = %d" %(nb_ins,nb_del,nb_sub))

    
    
    #NEEDLEMAN WUNSCH ALGORITHM
    matrix_nw = matrix.copy()
    matrix_nw = matrix[start:end]
    nw_aln_query, nw_aln_dna = NW(query,sequence[start:end])
    #print(nw_aln_query)
    #print(nw_aln_dna)
    for i in range(0,len(nw_aln_dna)-1):
        if nw_aln_dna[i]=='-':
            matrix_nw.insert(i, '-')
    score_nw_aln = 0
    for i in range(0,len(nw_aln_query)-1):
        nw_query_nuc = nw_aln_query[i]
        nw_dna_dict = matrix_nw[i]
        if (nw_query_nuc =='-' or nw_dna_dict=='-'):
            score_nw_aln -= 0
        else:
            score_nw_aln += nw_dna_dict[nw_query_nuc]
    score_nw_aln = score_nw_aln/len(nw_aln_query)
    #print("Score NW: %f" %score_nw_aln)
    
    print()
    seeds = find_seeds(query,word_size)
    seed_arr = np.array(seeds,dtype = int)
    print('Number of seeds found: %d' %len(seeds))
    HSPs = ungapped_ext()  
    print('Number of HSPs found: %d \n' %len(HSPs))
    #print(HSPs)
    gap_extensions = gapped_ext(gap)
    #print(np.array(gap_extensions))
    #print(gap_extensions[0][0],HSPs[0][0][1:-1],gap_extensions[0][3])
    m = 1
    matches = []
    for ext in gap_extensions:
        #print("MATCH #%d:" %m)
        matrix_start = ext['Matrix indices'][0]
        matrix_end = ext['Matrix indices'][1]
        #print(matrix_start,matrix_end)
        query_seq = str(ext['Left query'])+str(HSPs[0][0][1:-1])+str(ext['Right query'])
        query_start = query.find(query_seq)
        query_end = query_start+len(query_seq)
        #print(query_seq)
        dna_seq = str(ext['Left dna'])+str(HSPs[0][0][1:-1])+str(ext['Right dna'])
        #print(dna_seq)
        dna = matrix[matrix_start:matrix_end]
        #print(dna)
        for i in range(0,len(dna_seq)-1):
            if dna_seq[i]=='-':
                dna.insert(i, '-')
        if (len(dna)!=len(dna_seq) or len(dna_seq)!=len(query_seq)):
            break
        score_aln = 0
        score_log = 0
        for i in range(0,len(query_seq)-1):
            query_nuc = query_seq[i]
            dna_dict = dna[i]
            if (query_nuc =='-' or dna_dict=='-'):
                score_aln -= 0
            else:
                score_aln += dna_dict[query_nuc]
                #score_log += log(scaling_factor*dna_dict[query_nuc])
        #print(score_aln)
        score_aln = score_aln/len(query_seq)
        break_loop = 0
        for k in matches:
            if k[0] == score_aln:
                break_loop = 1
        if break_loop == 1:
            break
        #print("Score = %f" %score_aln)
        indices = [matrix_start,matrix_end]
        output = "MATCH #%d:" %m +'\n'+ query_seq + '\n'+dna_seq +'\n'+ "Score = " +str(score_aln)+'\n'+ "Matrix indices: " +str(indices)
        to_csv = [query_length,nb_ins,nb_del,nb_sub,word_size,len(seeds),len(HSPs),scaling_factor,delta,gap,start-matrix_start,end-matrix_end,score_aln,score_nw_aln]
        matches.append([score_aln, output,to_csv])
        
        m+=1
    
    sorted(matches, key=operator.itemgetter(0))   
    for match in matches:
        print(match[1])
        print("Start offset: %d" %match[2][10])
        print("End offset: %d" %match[2][11])
        if match == matches[0]:
            #print("WRITING")
            to_csv = match[2]
            with open(r'project_results.csv', 'a') as f:
                writer = csv.writer(f)
                writer.writerow(to_csv)

        print()
        
elapsed_time = time.process_time() - t0
print("Running time: %f seconds" %elapsed_time)


# In[ ]:




