#! /usr/bin/env python
import pandas as pd
import numpy as np
import scipy
import cooler
import multiprocess
from multiprocess import Pool
from functools import partial
from scipy import stats
import warnings
import  sys
import os
import argparse
import matplotlib.pyplot as plt
import gc
import psutil
import time
##convert coordinates to locus - bins
def preproc_data(locus_file,res): 
    locus_coord=pd.read_csv(locus_file,sep='\t',header=None)    
    locus_coord.columns=['chrom','start']
    locus_coord['start']=(locus_coord['start']/res).astype(int)
    Locus=locus_coord.drop_duplicates(['chrom','start'])
    
    return(Locus)

##write obs/exp maps to txt.files
def norm(cool_file,ch,res):
    ff=cooler.Cooler(cool_file)
    s=ff.matrix(balance=True).fetch(ch)
    D=np.zeros(len(s))
    for i in range (len(s)):
        D[i]=np.nanmean(s.diagonal(i))
     
    A=np.zeros((len(s),len(s)))
    n=len(s)
    for i in range (len(s)):
        row_id=np.arange(0,(n-i))
        col_id=np.arange(i,n)

        if D[i]!=0:
            A[row_id,col_id]=s[row_id,col_id]/D[i]
            A[col_id,row_id]=s[col_id,row_id]/D[i]

    return(A)



##create array of profiles of interactions for bins with locus of interest +- predefined distance
def create_profiles(hic_prepared,n,locus_list,ch):
        
    locus_tmp=locus_list.loc[locus_list['chrom']==ch]
    locus_tmp=locus_tmp.sort_values(['chrom','start'])
    locus_tmp=locus_tmp.reset_index(drop=True)
        
    col_idx = locus_tmp['start']+n
    col_idx_rep = np.repeat(col_idx, 2*n+1)

    row_inc = np.arange(-n, n+1)
    row_inc_tile = np.tile(row_inc, len(col_idx))
    row_idx = col_idx_rep + row_inc_tile
    profiles=hic_prepared[row_idx, col_idx_rep][:, np.newaxis].reshape(len(col_idx), -1)  
    
    return(profiles)
    
##create array of profiles for background   
def create_background_profiles(hic_prepared,n):
    col_idx = np.arange(n,len(hic_prepared)-n)
    col_idx_rep = np.repeat(col_idx, 2*n+1)
    row_inc = np.arange(-n, n+1)
    row_inc_tile = np.tile(row_inc, len(col_idx))
    row_idx = col_idx_rep + row_inc_tile
    profiles=hic_prepared[row_idx, col_idx_rep][:, np.newaxis].reshape(len(col_idx), -1)  
            
    return(profiles)
 
##calculate p-values by comparison with profiles for background SNP
def pvalue_calculation(profiles,my_profiles,n,k):
    my_p=np.zeros(np.shape(my_profiles))
    for i in np.append(np.arange(0,n-k),np.arange(n+k+1,2*n+1)):
        prof=profiles[:,i]
        prof=prof[~np.isnan(prof)]
        if len(prof)>0:
            q95=np.percentile(prof,q=[95])
            prof=prof[(prof>0)&(prof<=q95)]
            if (len(prof))>20:
                s1,loc1,scale1=scipy.stats.exponnorm.fit(prof)
                my_p[:,i]=1-scipy.stats.exponnorm.cdf(my_profiles[:,i],s1,loc1,scale1) 
 
    return(my_p)



def sign_contact(hic_prepared,locus_tested_preproc,res,ch,dist,dist_start): 
    
    n=int(dist/res)
    k=int(dist_start/res)
    m=hic_prepared
    m=np.pad(m,n)
    back_profiles=create_background_profiles(m,n)
    my_profiles=create_profiles(m,n,locus_tested_preproc,ch)
    pval=pvalue_calculation(back_profiles,my_profiles,n,k)
    
    return(pval)    



def postprocess_prepare(locus_tested_file,res):
    my_locus=pd.read_csv(locus_tested_file,sep='\t',header=None)
    my_locus[2]=my_locus[0].astype(str)+':'+my_locus[1].astype(str)
    my_locus[1]=(my_locus[1]/res).astype(int)
    my_locus=my_locus.groupby([0,1])[2].apply(list).reset_index()
    my_locus.sort_values([0,1],inplace=True)
    return(my_locus)


    
## fdr
def qvalue_calculation(pval_table,fdr):
    
    qval_table=np.zeros(np.shape(pval_table))
    for i in range(np.shape(pval_table)[1]):
        pval=pval_table[:,i]
        pval=pval[pval!=0]
        pval=pval[~(np.isnan(pval))]
        L=len(pval)
        rank=np.array([sorted(pval).index(x) for x in pval])+1
        u=pd.DataFrame({'pval':pval,'rank':rank,'qval':rank*fdr/L})
        u=u.loc[u.pval<u.qval]
        if len(u)>0:
            critical_q=max(u.pval)
        else:
            critical_q=-1

        qval_table[:,i]=critical_q 
    
    return(qval_table)


def significant_contacts_average_heatmap(norm_map,sign_contact_table,chrom,res,d=10):

    a=sign_contact_table[['chr','bin_to_test','interacting_bin_coord']]
    a.columns=[0,1,2]
    a=a.loc[abs(a[2]-a[1])>2*res]
    s=norm_map
    s=np.pad(s,d)
    w=a[a[0]==chrom]
    w.reset_index(drop=True,inplace=True)
    D=np.zeros((d*2+1,d*2+1,len(w)))
    k=0
    for i in range(len(w)):
        x=(w.iloc[i][1]/res).astype(int)
        y=(w.iloc[i][2]/res).astype(int)
        z=s[x:x+2*d+1,y:y+2*d+1]
        D[:,:,k]=z
        k+=1
    return(D)

def post_processing_for_sign_contacts(hic_prepared,pval_array,locus_tested_spec_file,res,dist,fdr,ch):# 
    
    n=int(dist/res)
    my_locus=locus_tested_spec_file.copy()
    my_locus_tmp=my_locus.loc[my_locus[0]==ch]
    pval_table=pval_array
    qval_table=qvalue_calculation(pval_table,fdr)
    snp_significant_contacts=pd.DataFrame()    
    for i in range(len(my_locus_tmp)):
        interacting_locus_coordinates=np.arange(my_locus_tmp.iloc[i][1]-n,my_locus_tmp.iloc[i][1]+n+1,1)*res

        df_tmp=pd.DataFrame({'chr':[my_locus_tmp.iloc[i][0]]*(2*n+1),
                             'bin_to_test':[my_locus_tmp.iloc[i][1]*res]*(2*n+1),
                             'list_of_initial_loci':[my_locus_tmp.iloc[i][2]]*(2*n+1),
                             'interacting_bin_coord':interacting_locus_coordinates,
                             'p-val':pval_table[i],
                             'p-value_critical':qval_table[i]})

        snp_significant_contacts=pd.concat([snp_significant_contacts, df_tmp.loc[(df_tmp['p-val']>0)&(df_tmp['p-val']<=df_tmp['p-value_critical'])]])

    score=[]
    for i in range(len(snp_significant_contacts)):
        x1=int(snp_significant_contacts.iloc[i]['bin_to_test']/res)
        x2=int(snp_significant_contacts.iloc[i]['interacting_bin_coord']/res)
        #tmp=hic_prepared[x1-20:x1+21,x2-20:x2+21]
        #diag=np.concatenate((np.diagonal(tmp,0),np.diagonal(tmp,1),np.diagonal(tmp,-1)))
        ranges=10
        
        n=int(dist/res)
        m=hic_prepared
        m=np.pad(m,n)
        
        tmp=m[x1+n-ranges:x1+n+ranges+1,x2+n-ranges:x2+n+ranges+1].copy()
        tmp[ranges][ranges]=np.nan
        diag=tmp.flatten()
        diag=diag[~np.isnan(diag)]
        #print(hic_prepared[x1,x2])
        if len(diag)>0:
        	score.append(np.float(sum(hic_prepared[x1,x2]>diag))/np.float(len(diag)))
        else:
        	score.append(np.nan)
    snp_significant_contacts['score']=score
    return(snp_significant_contacts[['chr','bin_to_test','list_of_initial_loci','interacting_bin_coord','p-val','score']])       
            
        


def create_final_table(cool,tested_locus,res,dist,dist_start,fdr,ch):
    file=norm(cool,ch,res)

    #back_locus=preproc_data(background_locus,res)
    test_locus=preproc_data(tested_locus,res)
    sign_cont_pval=sign_contact(file,
                            #back_locus,
                            test_locus,
                            res,ch,dist,dist_start) 
    print('..%d'%(len(sign_cont_pval)))
    locus_test_spec_file=postprocess_prepare(tested_locus,res)
    sign_cont_qval=post_processing_for_sign_contacts(file,sign_cont_pval,locus_test_spec_file,res,dist,fdr,ch)
    data_for_heatmap=significant_contacts_average_heatmap(file,sign_cont_qval,ch,res,d=10)
    #print(psutil.virtual_memory()[2])
    return(sign_cont_qval,data_for_heatmap)



def bad_distance(cool,distance,chrom):
    c=cooler.Cooler(cool)
    #print(chrom,c.chromsizes[chrom])
    #d=int(np.round(min(c.chromsizes[chrom])/4,-int(np.log10(min(c.chromsizes[chrom]/4)))))
    d=int(min(c.chromsizes[chrom])/4)    
    if  (d<distance):
        return(1)

def cpu_number(ch_sizes,res):
    mem_lim=(psutil.swap_memory()[2]+psutil.virtual_memory()[4])*0.9/1000000000
    ch=(ch_sizes/res/1000).sort_values(ascending=False)
    y=[0.016*np.power(i,2)+0.65 for i in ch]
    cum=np.cumsum(y)
    cpu=sum(cum<mem_lim)
    return(cpu)

def get_contacts_cli(cool,background_loci,tested_loci,res,dist,dist_start,chromosomes,fdr,heatmap_generate):
    #start_time = time.time()    
    c=cooler.Cooler(cool)
    warnings.filterwarnings("ignore", category=RuntimeWarning) 
 
    coord_background=pd.read_csv(background_loci,sep='\t',header=None)    
    coord_to_test=pd.read_csv(tested_loci,sep='\t',header=None)   

    if chromosomes==['all']:   
        chroms=list(set(coord_background[0])&set(coord_to_test[0]))
    else:
        chroms=list(set(coord_background[0])&set(coord_to_test[0])&set(chromosomes))

    if len(chroms)==0:
        raise ValueError('check the consistency of chromosome names across all input files')
    del(coord_background)
    gc.collect()

    bd=bad_distance(cool,dist,chroms)
    if bd:
        raise ValueError('the distance for contacts identification is too large: reduce distance, exclude contigs and chrM from a file with loci to be tested (if they present there)')


    #### define CPUs number from RAM #####
    cpus_ram=cpu_number(c.chromsizes[chroms],res)

    if cpus_ram<1:
        raise ValueError('the RAM capacity is not enough for the current resolution')
    else:
        proc=min(max(1,multiprocess.cpu_count()-2),len(chroms),cpus_ram)
    ##############################    
    print('The number of CPUs allocated for this task: %s'%proc)
    pool_obj = multiprocess.Pool(proc)

    func = partial(create_final_table,cool,background_loci,tested_loci,res,dist,fdr)
    results = [pool_obj.apply_async(func,args=(i,)) for i in chroms]

    try:
        for result in results:
            result.get()[0]
    except:
        pool_obj.terminate()
        print("error")

    pool_obj.close()
    pool_obj.join()

    df=pd.DataFrame()
    plot=np.zeros((np.shape(results[0].get()[1])[0:2]+(1,)))

    for result in results:
        df=pd.concat([df,result.get()[0]])
        df.sort_values(['chr','bin_to_test','interacting_bin_coord'],inplace=True)
        plot=np.dstack([plot,result.get()[1]])
    #print("--- %s seconds ---" % (time.time() - start_time))
    return(df,plot)

def plot_heatmap(data,path_to_image):
    plt.imshow(np.log(np.nanmean(data,axis=2)+0.01),cmap='coolwarm')
    plt.xticks([0,10,20],[-10,0,10])
    plt.yticks([0,10,20],[-10,0,10])
    plt.xlabel('distance, bins')
    plt.colorbar(label='ln(FC, obs/exp)')
    plt.title('average heatmap\n around significant contacts')
    plt.savefig(path_to_image)

