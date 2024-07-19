#   author: Nishat Anjum Bristy

# In[ ]:


from Bio import Entrez
import numpy as np
import pandas as pd
import sys

# Let NCBI know who you are
Entrez.email = "your email"

def get_full_refseq_ids(ucsc_refseq_annot_file,ref='hg19'):
    if ref=='hg38':
        column_names = [
        "Chromosome",
        "Start_position",
        "End_position",
        "Gene_Identifier",
        "Score",
        "Strand",
        "Thick_start",
        "Thick_end",
        "Item_RGB",
        "Block_count",
        "Block_sizes",
        "Block_starts"
        ]
        df = pd.read_csv(ucsc_refseq_annot_file,sep='\t',header=None, index_col=None,names=column_names)
    elif ref == 'hg19':
        df = pd.read_csv(ucsc_refseq_annot_file,sep='\t',index_col=None)
        
    return df

def get_W_chrom_pos(W,u,v, F_info): 
    ones = abs(W.iloc[u]-W.iloc[v])>=1
    ones_ind = ones[ones].index.tolist()
    
    return F_info.iloc[ones_ind]

def get_gene_name(refseq_id):
    handle = Entrez.efetch(db="nucleotide", id=refseq_id, retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    try:
        for feature in records[0]['GBSeq_feature-table']:
            if feature['GBFeature_key'] == 'gene':
                for qualifier in feature['GBFeature_quals']:
                    if qualifier['GBQualifier_name'] == 'gene':
                        return qualifier['GBQualifier_value']
    except:
        return None

def get_one_refseq(refseq_full,chrom, pos,chr_col,st_col,end_col,var_type):
    
    chrom_str = 'chr'+str(chrom)
    if var_type =='snv_sv':
        matching_rows = refseq_full[
            refseq_full[chr_col].str.contains(chrom_str) &
            (refseq_full[st_col] <= pos) &
            (refseq_full[end_col] >= pos)
        ]
#    elif var_type=='cnv':
#        matching_rows = refseq_full[
#            refseq_full[chr_col].str.contains(chrom_str) &
#            (refseq_full[st_col] >= pos) &
#            (refseq_full[end_col] <= (pos+500000))
#        ]   
    
    
    return matching_rows

def get_one_gene(refseq_full,chrom, pos,chr_col,st_col,end_col):
    if chrom==23:
        chrom=' X'
    chrom_str = 'chr'+str(chrom)
    
        
    matching_rows = refseq_full[
        refseq_full[chr_col].str.contains(chrom_str) &
        (refseq_full[st_col] <= pos) &
        (refseq_full[end_col] >= pos)
    ]
    return matching_rows

def get_your_refseq_ids(refseq_full,chromes,pos,ref_genome,var_type):
    genes = []
    for i in range(len(chromes)):
        
        if ref_genome=='hg38':
            matching_rows = get_one_refseq(refseq_full,chromes[i],pos[i],'Chromosome','Start_position','End_position',var_type)
            refseq_ids = list(matching_rows['Gene_Identifier'])
            if len(refseq_ids) == 0:
                continue
            gene = get_gene_name(refseq_ids)
            
        elif ref_genome=='hg19':
            matching_rows = get_one_refseq(refseq_full,chromes[i],pos[i],'chrom','txStart','txEnd',var_type)
            
            gene = list(dict.fromkeys(matching_rows['name2'])) # set(list(matching_rows['name2']))
            if len(gene)!=0:
                genes.append(','.join(gene))
                
    return genes



def get_genes_in_edge(W,u,v,F_info,var_type):
    print('Finding genes...')
    ref_genome = 'hg19'
    refseq_full = get_full_refseq_ids('hgTables_ucsc_refseq_hg19.txt',ref_genome)
    
    chorms_pos = get_W_chrom_pos(W,u,v,F_info)
    chroms = list(chorms_pos['Chrom'])
    pos = list(chorms_pos['Pos'])
    
    genes_in_edge = get_your_refseq_ids(refseq_full,chroms,pos,ref_genome,var_type)
    
        
    print('Edge',u,v,genes_in_edge)
    return (u,v,genes_in_edge)

def save_genes(W,file,folder,var_type,F_info):
    #var_type='snv_sv'
    f = open(folder+'/'+file,'w')
    # this part should read all the edges from the tree one by one. 
    # Maybe update this part to read directly from the tree instead of hard-coding.. 
    u,v,genes = get_genes_in_edge(W,6,5,F_info,var_type)
    f.write('Edge '+str(u)+', '+str(v)+','.join(genes)+'\n')
    u,v,genes = get_genes_in_edge(W,5,4,F_info,var_type)
    f.write('Edge '+str(u)+', '+str(v)+','.join(genes)+'\n')
    u,v,genes = get_genes_in_edge(W,5,3,F_info,var_type)
    f.write('Edge '+str(u)+', '+str(v)+','.join(genes)+'\n')
    u,v,genes = get_genes_in_edge(W,6,1,F_info,var_type)
    f.write('Edge '+str(u)+', '+str(v)+','.join(genes)+'\n')
    u,v,genes = get_genes_in_edge(W,6,2,F_info,var_type)
    f.write('Edge '+str(u)+', '+str(v)+','.join(genes)+'\n')
    u,v,genes = get_genes_in_edge(W,6,0,F_info,var_type)
    f.write('Edge '+str(u)+', '+str(v)+','.join(genes)+'\n')
    f.close()

folder = sys.argv[1]
W = pd.read_csv(folder+'W.tsv',sep='\t',index_col=None,header=None) 
W_SV = pd.read_csv(folder+'W_SV.tsv',sep='\t',index_col=None,header=None)
F_info = pd.read_csv(folder+'F_info_phasing.csv',sep='\t',index_col=None,header=None).iloc[:W.shape[1]]
F_info.columns = ['Chrom','Pos','variant']
save_genes(W_SV,'genes_SV.txt',folder,'snv_sv',F_info)


# %%
