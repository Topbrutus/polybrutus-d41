#!/usr/bin/env python3
"""Independent D23 rational-metric verifier; does not import primary coordinates."""
from fractions import Fraction as F
from collections import Counter
from itertools import combinations
from pathlib import Path
import hashlib,json
import sympy as s
from sympy.polys.matrices import DomainMatrix
ROOT=Path(__file__).resolve().parents[1]
v=[(F(0),F(0),F(1)),(F(1,2),F(0),F(1,2)),
(F(-1,4),F(3,4),F(1,2)),(F(-1,4),F(-3,4),F(1,2)),
(F(0),F(0),F(-1)),(F(1,2),F(0),F(-1,2)),(F(-1,2),F(0),F(-1,2))]
def dist2(i,j):
    if i==j:return F(0)
    if i<7 and j<7:
        dx,dy,dz=(v[i][k]-v[j][k] for k in range(3))
        return 3*dx*dx+dy*dy+dz*dz
    if 7<=i<12 and 7<=j<12:return F(5,2)
    return F(2)
N=28
pairs=list(combinations(range(N),2))
model={(i+1,j+1):dist2(i,j) for i,j in pairs}
p=ROOT/'results'/'certificate_d23.json'
cert=json.loads(p.read_text(encoding='utf8'))
assert cert['schema']=='POLYBRUTUS-D23-EXACT-CERTIFICATE-v1'
assert cert['status']=='PASS'
record={(int(e['i']),int(e['j'])):F(e['d2']) for e in cert['labelled_squared_distances']}
assert len(record)==len(model)==378 and record==model
parent=json.loads((ROOT/'prior'/'d22'/'results'/'certificate_d22.json').read_text(encoding='utf8'))
old={(int(e['i']),int(e['j'])):F(e['d2']) for e in parent['labelled_squared_distances']}
assert len(old)==351 and all(model[k]==old[k] for k in old)
assert all(dist2(i,27)==2 for i in range(27))
assert sum(model.values(),F(0))==F(3071,4)
hist=Counter(4*dist2(i,j) for i,j in pairs)
assert hist==Counter({F(4):6,F(7):2,F(8):347,F(9):3,F(10):10,F(12):6,F(13):2,F(16):2})
D=s.Matrix(N,N,lambda i,j:s.Rational(dist2(i,j).numerator,dist2(i,j).denominator))
H=s.eye(N)-s.ones(N)/N
assert DomainMatrix.from_Matrix(-H*D*H/2).to_field().rank()==23
P=D[:27,:27]; Hp=s.eye(27)-s.ones(27)/27
assert DomainMatrix.from_Matrix(-Hp*P*Hp/2).to_field().rank()==22
L=s.zeros(N)
for i,j in pairs:
    w=1/D[i,j]
    L[i,i]+=w;L[j,j]+=w;L[i,j]-=w;L[j,i]-=w
assert DomainMatrix.from_Matrix(L-14*s.eye(N)).to_field().rank()==11
assert DomainMatrix.from_Matrix(L-s.Rational(27,2)*s.eye(N)).to_field().rank()==24
assert cert['affine_rank']==cert['centered_gram_rank']==23
assert cert['parent_distances_preserved']==351 and cert['new_distances']==27
damaged=model.copy();damaged[(1,28)]=F(9,4)
assert damaged!=record
swap=lambda i: 8 if i==1 else 1 if i==8 else i
relabel={tuple(sorted((swap(i),swap(j)))):value for (i,j),value in model.items()}
assert relabel!=model and Counter(relabel.values())==Counter(model.values())
print('INDEPENDENT D23 PASS: 378/378 pairs; 351/351 prior pairs, Gram rank 23')
print('INDEPENDENT SPECTRAL PASS: eigenvalue 14 x17, 27/2 x4')
print('CERTIFICATE_SHA256='+hashlib.sha256(p.read_bytes()).hexdigest())
