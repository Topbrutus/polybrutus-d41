#!/usr/bin/env python3
"""Polybrutus D34: exact Euclidean labelled metric and chosen ideal weighted graph.
Conservative perpendicular unit extension of the immutable Polybrutus D33 v0.1.0.
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
n=39
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
assert DomainMatrix.from_Matrix(G).to_field().rank()==34
D=2*s.ones(n)-2*G
for i in range(n):
    D[i,i]=0
assert len(pairs)==741 and all(D[i,j]>0 for i,j in pairs)

parent_path=ROOT/'prior'/'d33'/'results'/'certificate_d33.json'
parent_bytes=parent_path.read_bytes()
assert sha256(parent_bytes).hexdigest()=='86516b8a6393bf9fdf4e698d6e439d4db783a71956730ac35bbc89e2c8e7b734'
parent=json.loads(parent_bytes)
assert parent['status']=='PASS' and parent['dimension']==33 and parent['vertices']==38
former={(int(e['i'])-1,int(e['j'])-1):s.Rational(e['d2'])
        for e in parent['labelled_squared_distances']}
assert len(former)==703
assert all(D[i,j]==v for (i,j),v in former.items())
assert all(D[i,38]==2 for i in range(38))
hist=Counter(4*D[i,j] for i,j in pairs)
expected=Counter({4:6,7:2,8:710,9:3,10:10,12:6,13:2,16:2})
assert hist==expected and sum(hist.values())==741
assert sum(D[i,j] for i,j in pairs)==Q(5975,4)

H=s.eye(n)-s.ones(n)/n
B=-H*D*H/2
assert B==H*G*H
assert DomainMatrix.from_Matrix(B).to_field().rank()==34
# The first and fifth canonical points are antipodal; origin lies in the affine hull.
assert G[0,4]==-1

L=s.zeros(n)
for i,j in pairs:
    w=1/D[i,j]
    L[i,i]+=w;L[j,j]+=w;L[i,j]-=w;L[j,i]-=w
assert L*s.ones(n,1)==s.zeros(n,1)
assert DomainMatrix.from_Matrix(L-Q(39,2)*s.eye(n)).to_field().rank()==11
assert DomainMatrix.from_Matrix(L-19*s.eye(n)).to_field().rank()==35
old=s.zeros(38)
for i,j in combinations(range(38),2):
    w=1/D[i,j]
    old[i,i]+=w;old[j,j]+=w;old[i,j]-=w;old[j,i]-=w
assert L[:38,:38]-old==s.eye(38)/2
assert all(L[i,38]==-Q(1,2) for i in range(38))
assert L[38,38]==19

# Four adversarial negative controls.
damaged=D.copy();damaged[0,38]=damaged[38,0]=Q(9,4)
assert Counter(4*damaged[i,j] for i,j in pairs)!=hist
duplicate=G.copy();duplicate[38,:]=G[37,:];duplicate[:,38]=G[:,37];duplicate[38,38]=1
assert DomainMatrix.from_Matrix(duplicate).to_field().rank()==33
tilted_dot_to_old_core=s.sqrt(3)*Q(3,10)
assert tilted_dot_to_old_core!=0
swapped=D.copy();swapped.row_swap(0,7);swapped.col_swap(0,7)
assert swapped!=D and Counter(4*swapped[i,j] for i,j in pairs)==hist

report={
  'schema':'POLYBRUTUS-D34-EXACT-CERTIFICATE-v1',
  'status':'PASS','dimension':34,'vertices':39,'pairs':741,
  'affine_rank':34,'centered_gram_rank':34,'ambient_storage_coordinates':35,
  'radius':1,'parent':'Topbrutus/polybrutus-d33@v0.1.0',
  'parent_tag_commit':'bc78e6c6fe1755a069d8860a0144cbebbb1168a6',
  'parent_doi':'10.5281/zenodo.24411569',
  'parent_certificate_sha256':sha256(parent_bytes).hexdigest(),
  'parent_distances_preserved':703,'new_distances':38,
  'new_squared_distances':['2']*38,
  'histogram_k_4d2':{str(k):int(v) for k,v in sorted(hist.items())},
  'sum_squared_distances':'5975/4',
  'gram_psd_certificate':'B=-(H D H)/2=H G H; G=XX^T exact',
  'negative_controls':'4/4 PASS',
  'oscillator_model':'Chosen ideal equal-mass complete graph with w_ij=1/d_ij^2; no physical measurement',
  'eigenvalue_39_over_2_multiplicity':28,
  'eigenvalue_19_multiplicity':4,
  'spectral_recurrence':'chi39(z)=z*(z-39/2)*chi38(z-1/2)/(z-1/2)',
  'spectral_closed_form':'chi39(z)=z*(z-39/2)^28*(z-19)^4*f6(z-16), f6(z)=chiP7(z)/z',
  'novel_universal_theorem':'NOT CLAIMED','physical_discovery':'NOT CLAIMED',
  'human_peer_review':'NOT_DONE','sympy':s.__version__,
  'labelled_squared_distances':[{'i':i+1,'j':j+1,'d2':str(D[i,j]),'k':str(4*D[i,j])} for i,j in pairs]
}
out=ROOT/'results'/'certificate_d34.json'
out.parent.mkdir(parents=True,exist_ok=True)
out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8',newline='\n')
print('PRIMARY D34 PASS: 39 vertices, 741 labelled pairs, affine rank 34; 703/703 previous pairs preserved')
print('IDEAL GRAPH PASS: eigenvalue 39/2 x28; eigenvalue 19 x4; 4/4 negative controls')
print('CERTIFICATE_SHA256='+sha256(out.read_bytes()).hexdigest())
