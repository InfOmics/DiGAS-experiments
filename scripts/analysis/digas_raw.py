"""Re-implementation of DiGAS raw scores (Hill per group, FC), without permutations.
Mirrors scripts/nested/digas3-linux_validation.py + PLINK 1.9 behaviour:
 - SNP->region mapping with strict interval, duplicate SNPs across a gene's rows kept (weights)
 - each group is written with `plink --keep --make-bed`, which re-sets A1 to the minor allele *within that group*
 - recodeA counts A1; DiGAS codes 1/2 -> carrier, 0/missing -> non-carrier
 - chr23 (X) in males: PLINK counts one allele per hemizygous call and skips het calls when choosing the minor allele,
   but recodeA still writes male het calls as 1 -> carrier."""
import numpy as np, pandas as pd, scipy.sparse as sp, os
ROOT=os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),"..",".."))   # repo root

def build(region="exon_5k"):
    """Load ADNI1_QC genotypes for the SNPs inside `region/<region>.txt` and the SNP->gene weight matrix."""
    bim=pd.read_csv(f"{ROOT}/ADNI1_QC.bim",sep="\t",header=None,names=["chr","id","cm","bp","a1","a2"])
    fam=pd.read_csv(f"{ROOT}/ADNI1_QC.fam",sep=r"\s+",header=None)
    reg=pd.read_csv(f"{ROOT}/region/{region}.txt",sep="\t",header=None,names=["g","s","e"])
    reg["chr"]=reg.s.str.split(":").str[0].astype(int)
    reg["s"]=reg.s.str.split(":").str[1].astype(int); reg["e"]=reg.e.str.split(":").str[1].astype(int)
    rows=[];cols=[];genes={}
    for c,rc in reg.groupby("chr"):
        b=bim[bim.chr==c]; bp=b.bp.values; idx=b.index.values; o=np.argsort(bp,kind="stable"); bp=bp[o]; idx=idx[o]
        lo=np.searchsorted(bp,rc.s.values,side="right"); hi=np.searchsorted(bp,rc.e.values,side="left")
        for g,l,h in zip(rc.g.values,lo,hi):
            if h>l:
                j=genes.setdefault(g,len(genes)); rows+=list(idx[l:h]); cols+=[j]*(h-l)
    used=np.unique(rows); pos={s:i for i,s in enumerate(used)}
    W=sp.coo_matrix((np.ones(len(rows)),([pos[r] for r in rows],cols)),shape=(len(used),len(genes))).tocsc()
    Wc=sp.csc_matrix(W)                                            # SNP multiplicity per gene (integer)
    W=sp.csc_matrix(W.multiply(1.0/np.asarray(W.sum(0)).ravel()))
    n=len(fam); nb=(n+3)//4
    raw=np.fromfile(f"{ROOT}/ADNI1_QC.bed",dtype=np.uint8); assert raw[:3].tolist()==[0x6c,0x1b,1]
    M=raw[3:].reshape(len(bim),nb)[used]
    codes=np.stack([(M>>(2*k))&3 for k in range(4)],axis=2).reshape(len(used),nb*4)[:,:n]   # SNP x subj
    male=(fam[4].values==1); hap=(bim.chr.values[used]==23)[:,None]&male[None,:]
    hom1=(codes==0); het=(codes==2); hom2=(codes==3)
    gl=[None]*len(genes)
    for g,j in genes.items(): gl[j]=g
    return dict(hom1=hom1,het=het,hom2=hom2,hap=hap,W=W,Wc=Wc,genes=gl,iids=list(fam[1]))

def hills(D,mask,fixed_allele=False):
    """Hill per gene for the subjects in `mask`. fixed_allele=True counts the same allele (cohort A1) in every group."""
    hp=D["hap"][:,mask]
    h1=D["hom1"][:,mask].sum(1); h2=D["hom2"][:,mask].sum(1); he=D["het"][:,mask].sum(1)
    he_dip=(D["het"][:,mask]&~hp).sum(1)                       # hets counted for allele frequency
    a1=2*h1+he_dip-(D["hom1"][:,mask]&hp).sum(1); a2=2*h2+he_dip-(D["hom2"][:,mask]&hp).sum(1)
    swap=np.zeros_like(a1,dtype=bool) if fixed_allele else (a1>a2)
    carriers=np.where(swap,h2+he,h1+he)/mask.sum()
    return np.asarray(D["W"].T@carriers).ravel()

def fc(h1,h2): return np.abs(np.log1p(h1)-np.log1p(h2))
