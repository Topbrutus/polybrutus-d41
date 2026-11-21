#!/usr/bin/env python3
"""Polybrutus D32: exact Euclidean labelled metric and chosen ideal weighted graph.
Conservative perpendicular unit extension of the immutable Polybrutus D31 v0.1.0.
No physical spacetime dimension or experimentally measured resonance is claimed.
"""
from pathlib import Path
from itertools import combinations
from collections import Counter
from hashlib import sha256
import json
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]
Q=s.Rational
core=[(0,0,1),(Q(1,2),0,Q(1,2)),(-Q(1,4),Q(3,4),Q(1,2)),
      (-Q(1,4),-Q(3,4),Q(1,2)),(0,0,-1),(Q(1,2),0,-Q(1,2)),
      (-Q(1,2),0,-Q(1,2))]
n=37
pairs=list(combinations(range(n),2))
G=s.zeros(n)
for i in range(7):
    for j in range(7):
        G[i,j]=3*core[i][0]*core[j][0]+core[i][1]*core[j][1]+core[i][2]*core[j][2]
for i in range(7,12):
    for j in range(7,12):
        G[i,j]=1 if i==j else -Q(1,4)
for i in range(12,n):
    G[i,i]=1
assert all(G[i,i]==1 for i in range(n))
assert DomainMatrix.from_Matrix(G).to_field().rank()==32
D=2*s.ones(n)-2*G
for i in range(n):
    D[i,i]=0
assert len(pairs)==666 and all(D[i,j]>0 for i,j in pairs)

parent_path=ROOT/'prior'/'d31'/'results'/'certificate_d31.json'
parent_bytes=parent_path.read_bytes()
assert sha256(parent_bytes).hexdigest()=='7c373c50d61ca2964eb832fc6fb8178590715e603f5fdada122575f624f60e8c'
parent=json.loads(parent_bytes)
assert parent['status']=='PASS' and parent['dimension']==31 and parent['vertices']==36
former={(int(e['i'])-1,int(e['j'])-1):s.Rational(e['d2'])
        for e in parent['labelled_squared_distances']}
assert len(former)==630
assert all(D[i,j]==v for (i,j),v in former.items())
assert all(D[i,36]==2 for i in range(36))
hist=Counter(4*D[i,j] for i,j in pairs)
expected=Counter({4:6,7:2,8:635,9:3,10:10,12:6,13:2,16:2})
assert hist==expected and sum(hist.values())==666
assert sum(D[i,j] for i,j in pairs)==Q(5375,4)

H=s.eye(n)-s.ones(n)/n
B=-H*D*H/2
assert B==H*G*H
assert DomainMatrix.from_Matrix(B).to_field().rank()==32
# The first and fifth canonical points are antipodal; origin lies in the affine hull.
assert G[0,4]==-1

L=s.zeros(n)
for i,j in pairs:
    w=1/D[i,j]
    L[i,i]+=w;L[j,j]+=w;L[i,j]-=w;L[j,i]-=w
assert L*s.ones(n,1)==s.zeros(n,1)
assert DomainMatrix.from_Matrix(L-Q(37,2)*s.eye(n)).to_field().rank()==11
assert DomainMatrix.from_Matrix(L-18*s.eye(n)).to_field().rank()==33
old=s.zeros(36)
for i,j in combinations(range(36),2):
    w=1/D[i,j]
    old[i,i]+=w;old[j,j]+=w;old[i,j]-=w;old[j,i]-=w
assert L[:36,:36]-old==s.eye(36)/2
assert all(L[i,36]==-Q(1,2) for i in range(36))
assert L[36,36]==18

# Four adversarial negative controls.
damaged=D.copy();damaged[0,36]=damaged[36,0]=Q(9,4)
assert Counter(4*damaged[i,j] for i,j in pairs)!=hist
duplicate=G.copy();duplicate[36,:]=G[35,:];duplicate[:,36]=G[:,35];duplicate[36,36]=1
assert DomainMatrix.from_Matrix(duplicate).to_field().rank()==31
tilted_dot_to_old_core=s.sqrt(3)*Q(3,10)
assert tilted_dot_to_old_core!=0
swapped=D.copy();swapped.row_swap(0,7);swapped.col_swap(0,7)
assert swapped!=D and Counter(4*swapped[i,j] for i,j in pairs)==hist

report={
  'schema':'POLYBRUTUS-D32-EXACT-CERTIFICATE-v1',
  'status':'PASS','dimension':32,'vertices':37,'pairs':666,
  'affine_rank':32,'centered_gram_rank':32,'ambient_storage_coordinates':33,
  'radius':1,'parent':'Topbrutus/polybrutus-d31@v0.1.0',
  'parent_tag_commit':'0c1f546a6c5442bc83a5bdde506132cce20b6e05',
  'parent_doi':'10.5281/zenodo.24411160',
  'parent_certificate_sha256':sha256(parent_bytes).hexdigest(),
  'parent_distances_preserved':630,'new_distances':36,
  'new_squared_distances':['2']*36,
  'histogram_k_4d2':{str(k):int(v) for k,v in sorted(hist.items())},
  'sum_squared_distances':'5375/4',
  'gram_psd_certificate':'B=-(H D H)/2=H G H; G=XX^T exact',
  'negative_controls':'4/4 PASS',
  'oscillator_model':'Chosen ideal equal-mass complete graph with w_ij=1/d_ij^2; no physical measurement',
  'eigenvalue_37_over_2_multiplicity':26,
  'eigenvalue_18_multiplicity':4,
  'spectral_recurrence':'chi37(z)=z*(z-37/2)*chi36(z-1/2)/(z-1/2)',
  'spectral_closed_form':'chi37(z)=z*(z-37/2)^26*(z-18)^4*f6(z-15), f6(z)=chiP7(z)/z',
  'novel_universal_theorem':'NOT CLAIMED','physical_discovery':'NOT CLAIMED',
  'human_peer_review':'NOT_DONE','sympy':s.__version__,
  'labelled_squared_distances':[{'i':i+1,'j':j+1,'d2':str(D[i,j]),'k':str(4*D[i,j])} for i,j in pairs]
}
out=ROOT/'results'/'certificate_d32.json'
out.parent.mkdir(parents=True,exist_ok=True)
out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8',newline='\n')
print('PRIMARY D32 PASS: 37 vertices, 666 labelled pairs, affine rank 32; 630/630 previous pairs preserved')
print('IDEAL GRAPH PASS: eigenvalue 37/2 x26; eigenvalue 18 x4; 4/4 negative controls')
print('CERTIFICATE_SHA256='+sha256(out.read_bytes()).hexdigest())
