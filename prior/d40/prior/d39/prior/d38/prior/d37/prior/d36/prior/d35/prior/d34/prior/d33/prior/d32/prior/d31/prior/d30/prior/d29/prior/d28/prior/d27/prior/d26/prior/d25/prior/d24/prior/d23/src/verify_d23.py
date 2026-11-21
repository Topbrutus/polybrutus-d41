#!/usr/bin/env python3
"""Polybrutus D23: exact Euclidean metric, predecessor preservation, ideal graph.
The Euclidean Gram matrix uses rational scalar products: canonical seven-point
core, a five-point centered regular simplex, and sixteen orthogonal unit axes.
No experimental spacetime or measured oscillations are asserted.
"""
from itertools import combinations
from collections import Counter
from pathlib import Path
import hashlib
import json
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]
Q=s.Rational
core=[(0,0,1),(Q(1,2),0,Q(1,2)),(-Q(1,4),Q(3,4),Q(1,2)),
      (-Q(1,4),-Q(3,4),Q(1,2)),(0,0,-1),(Q(1,2),0,-Q(1,2)),
      (-Q(1,2),0,-Q(1,2))]
n=28
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
assert DomainMatrix.from_Matrix(G).to_field().rank()==23
D=2*s.ones(n)-2*G
for i in range(n):
    D[i,i]=0
assert len(pairs)==378 and all(D[i,j]>0 for i,j in pairs)
assert all(D[i,27]==2 for i in range(27))

# Verify all 351 prior labelled distances against the immutable D22 certificate.
parent_path=ROOT/'prior'/'d22'/'results'/'certificate_d22.json'
parent=json.loads(parent_path.read_text(encoding='utf8'))
assert parent['status']=='PASS' and parent['dimension']==22
former={(int(e['i'])-1,int(e['j'])-1):s.Rational(e['d2'])
        for e in parent['labelled_squared_distances']}
assert len(former)==351
assert all(D[i,j]==v for (i,j),v in former.items())
hist=Counter(4*D[i,j] for i,j in pairs)
assert hist==Counter({4:6,7:2,8:347,9:3,10:10,12:6,13:2,16:2})
assert sum(D[i,j] for i,j in pairs)==Q(3071,4)
H=s.eye(n)-s.ones(n)/n
B=-H*D*H/2
assert B==H*G*H
assert DomainMatrix.from_Matrix(B).to_field().rank()==23
# p1 and p5 are antipodal; their midpoint is the origin in the affine hull.
assert G[0,4]==-1

L=s.zeros(n)
for i,j in pairs:
    w=1/D[i,j]
    L[i,i]+=w; L[j,j]+=w; L[i,j]-=w; L[j,i]-=w
assert L*s.ones(n,1)==s.zeros(n,1)
assert DomainMatrix.from_Matrix(L-14*s.eye(n)).to_field().rank()==11
assert DomainMatrix.from_Matrix(L-Q(27,2)*s.eye(n)).to_field().rank()==24
priorpairs=list(combinations(range(27),2))
old=s.zeros(27)
for i,j in priorpairs:
    w=1/D[i,j]
    old[i,i]+=w; old[j,j]+=w; old[i,j]-=w; old[j,i]-=w
assert L[:27,:27]-old==s.eye(27)/2
assert all(L[i,27]==-Q(1,2) for i in range(27))
assert L[27,27]==Q(27,2)
coreL=old[:7,:7]-10*s.eye(7)
z=s.Symbol('z')
f6=s.cancel(coreL.charpoly(z).as_expr()/z)
assert f6.subs(z,Q(7,2))!=0 and f6.subs(z,3)!=0

# Four deliberately invalid configurations are rejected.
damaged=D.copy()
damaged[0,27]=damaged[27,0]=Q(9,4)
assert Counter(4*damaged[i,j] for i,j in pairs)!=hist
duplicate=G.copy()
duplicate[27,:]=G[26,:];duplicate[:,27]=G[:,26];duplicate[27,27]=1
assert DomainMatrix.from_Matrix(duplicate).to_field().rank()==22
tilted_dot_P2=Q(3,5)*s.sqrt(3)*core[1][0]
assert tilted_dot_P2==3*s.sqrt(3)/10 and s.simplify(2-2*tilted_dot_P2)!=2
swapped=D.copy()
swapped.row_swap(0,7);swapped.col_swap(0,7)
assert swapped!=D and Counter(4*swapped[i,j] for i,j in pairs)==hist

report={
 'schema':'POLYBRUTUS-D23-EXACT-CERTIFICATE-v1',
 'status':'PASS','dimension':23,'vertices':28,'pairs':378,
 'affine_rank':23,'centered_gram_rank':23,
 'ambient_storage_coordinates':24,'radius':1,
 'parent':'Topbrutus/polybrutus-d22@v0.1.0',
 'parent_tag_commit':'4d4e80a7156a118f8f2bebf63aca76ac5f0f06eb',
 'parent_doi':'10.5281/zenodo.24409636',
 'parent_certificate_github_blob_sha':'33b1bd18e01b642b9098832f84eeff6679a71466',
 'parent_distances_preserved':351,'new_distances':27,
 'new_squared_distances':['2']*27,
 'histogram_k_4d2':{str(k):int(v) for k,v in sorted(hist.items())},
 'sum_squared_distances':'3071/4',
 'gram_psd_certificate':'B=-(H D H)/2=H G H; G=XX^T exactly',
 'negative_controls':'4/4 PASS',
 'oscillator_model':'chosen ideal unit masses and w_ij=1/d_ij^2; not measured',
 'eigenvalue_14_multiplicity':17,
 'eigenvalue_27_over_2_multiplicity':4,
 'spectral_recurrence':'chi28(z)=z*(z-14)*chi27(z-1/2)/(z-1/2)',
 'spectral_closed_form':'chi28(z)=z*(z-14)^17*(z-27/2)^4*f6(z-21/2)',
 'novel_universal_theorem':'NOT CLAIMED',
 'physical_discovery':'NOT CLAIMED','human_peer_review':'NOT_DONE',
 'sympy':s.__version__,
 'labelled_squared_distances':[{'i':i+1,'j':j+1,'d2':str(D[i,j]),'k':str(4*D[i,j])}
      for i,j in pairs]
}
out=ROOT/'results'/'certificate_d23.json'
out.parent.mkdir(parents=True,exist_ok=True)
out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf8',newline='\n')
print('PRIMARY D23 PASS: 28 vertices, 378 pairs, rank 23, 351/351 parent distances')
print('IDEAL GRAPH: eigenvalue 14 x17; eigenvalue 27/2 x4; negatives 4/4')
print('CERTIFICATE_SHA256='+hashlib.sha256(out.read_bytes()).hexdigest())
