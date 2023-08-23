import subprocess
import pandas as pd
from   pathlib import Path

candidates = [#"1BKSA",
              #"1BVCA",
              #"1DDRA",
              #"1FMKA",
              #"1IO2A",
              #"1IOBA",
              "1LNIA",
              #"1LZ1A",
              #"1O6XA",
              #"1RHGA",
              #"1RN1A",
              #"1ROPA",
              #"1RTBA",
              #"1RX4A",
              #"1SARA",
              #"1STNA",
              #"1TTGA",
              #"1XXNA",
              #"2JIEA",
              #"2NVHA",
              #"2RN2A",
              #"2TRXA",
              #"2WSYA",
              #"3PGKA",
              #"5CROA",
              #"5DFRA"
              ]


Path("./temp").mkdir(exist_ok=True)
d1        = pd.read_csv("./DATASETS/s669.csv")
d2        = pd.read_csv("./DATASETS/s669_r.csv")
d3        = pd.read_csv("./DATASETS/new_ds_for_s669.csv")
### df - full PROSTATA train + test set including reverse mutations 
df        = pd.concat([d1,d2,d3])

### save all sequences to one file and prepare BLAST database
used = set()
with open("./temp/all_sequences.fasta",'w') as fo:
    for d_ in df.iloc():
        if d_['pdb_id'] in used:
            continue
        fo.write(f">{d_['pdb_id']}\n{d_['wt_seq']}\n")
        used.add(d_['pdb_id'])
subprocess.call(f"makeblastdb -in all_sequences.fasta -dbtype prot", shell=True, cwd="./temp/")

### using BLAST to find sequencies with homology to the analyzed protein to put them to the test set
### the rest is kept as a train set
Path("./cross_validation_datasets").mkdir(exist_ok=True)
used = set()
for pdb_id in candidates:
    pdb_id_sequence =  df[df["pdb_id"]==pdb_id].iloc()[0]["wt_seq"]
    with open("./temp/target.fasta",'w') as fo:
        fo.write(f">{pdb_id}\n{pdb_id_sequence}\n")
    
    ### using BLAST to find sequencies with homology
    subprocess.call(f"blastp -db all_sequences.fasta -query target.fasta -outfmt 6 -out hits.tsv -num_threads 4", shell=True, cwd="./temp/")
    hit_data    = pd.read_csv("./temp/hits.tsv", delimiter='\t', header=None)

    ### filter sequences that has >30% identity and p-value < 0.05
    hit_data    = hit_data[(hit_data.iloc[:,2]>30.0) & (hit_data.iloc[:,-2]<0.05)]
    hit_pdb_ids = set(hit_data.iloc[:,1])
    ids         = df["pdb_id"].isin(hit_pdb_ids)
    if len(hit_pdb_ids & used)!=0:
        continue
    used     |= set(hit_pdb_ids)
    train_df  = df[~ids].reset_index(drop=True)
    test_df   = df[ids].reset_index( drop=True)
    
    ### save train and test sets
    train_df.to_csv(f"./cross_validation_datasets/train_{pdb_id}.csv")
    test_df.to_csv( f"./cross_validation_datasets/test_{pdb_id}.csv")
    print(pdb_id, train_df.shape, test_df.shape)

