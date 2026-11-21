#!/usr/bin/env python3
"""Independent D25 rational-distance verifier, no import of primary Gram generator."""
from fractions import Fraction as F
from collections import Counter
from itertools import combinations
from pathlib import Path
import json,hashlib
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]
p=ROOT/'results'/'certificate_d25.json'
cert=json.loads(p.read_text(encoding='utf-8'))
parent_file=ROOT/'prior'/'d24'/'results'/'certificate_d24.json'
parent=json.loads(parent_file.read_text(encoding='utf-8'))
assert hashlib.sha256(parent_file.read_bytes()).hexdigest()=='c2b1bbfd2c761871e81748fc5621c700443c7a7e7c9daba479d9ff18112c22eb'
core=[
    (F(0),F(0),F(1)),(F(1,2),F(0),F(1,2)),
    (F(-1,4),F(3,4),F(1,2)),(F(-1,4),F(-3,4),F(1,2)),
    (F(0),F(0),F(-1)),(F(1,2),F(0),F(-1,2)),(F(-1,2),F(0),F(-1,2))
]
def d2(i,j):
    if i==j:return F(0)
    if i<7 and j<7:
        a,b,c=(core[i][k]-core[j][k] for k in range(3))
        return 3*a*a+b*b+c*c
    if 7<=i<12 and 7<=j<12:return F(5,2)
    return F(2)

N=30
pairs=list(combinations(range(N),2))
computed={(i+1,j+1):d2(i,j) for i,j in pairs}
supplied={(int(e['i']),int(e['j'])):F(e['d2']) for e in cert['labelled_squared_distances']}
assert cert['schema']=='POLYBRUTUS-D25-EXACT-CERTIFICATE-v1' and cert['status']=='PASS'
assert len(computed)==len(supplied)==435 and computed==supplied
old={(int(e['i']),int(e['j'])):F(e['d2']) for e in parent['labelled_squared_distances']}
assert len(old)==406 and all(computed[k]==v for k,v in old.items())
assert all(d2(i,29)==F(2) for i in range(29))
expected=Counter({F(4):6,F(7):2,F(8):404,F(9):3,F(10):10,F(12):6,F(13):2,F(16):2})
assert Counter(4*d2(i,j) for i,j in pairs)==expected
assert sum(computed.values(),F(0))==F(3527,4)

D=s.Matrix(N,N,lambda i,j:s.Rational(d2(i,j).numerator,d2(i,j).denominator))
H=s.eye(N)-s.ones(N)/N
assert DomainMatrix.from_Matrix(-H*D*H/2).to_field().rank()==25
oldD=D[:29,:29];Ho=s.eye(29)-s.ones(29)/29
assert DomainMatrix.from_Matrix(-Ho*oldD*Ho/2).to_field().rank()==24

L=s.zeros(N)
for i,j in pairs:
    w=1/D[i,j]
    L[i,i]+=w;L[j,j]+=w;L[i,j]-=w;L[j,i]-=w
assert DomainMatrix.from_Matrix(L-15*s.eye(N)).to_field().rank()==11
assert DomainMatrix.from_Matrix(L-s.Rational(29,2)*s.eye(N)).to_field().rank()==26
assert cert['affine_rank']==cert['centered_gram_rank']==25
assert cert['parent_distances_preserved']==406 and cert['new_distances']==29
# Negative tests on labelled metric and histogram-insensitive relabelling.
wrong=computed.copy();wrong[(1,30)]=F(9,4)
assert wrong!=supplied
swap=lambda i:8 if i==1 else 1 if i==8 else i
permuted={tuple(sorted((swap(i),swap(j)))):value for (i,j),value in computed.items()}
assert permuted!=computed and Counter(permuted.values())==Counter(computed.values())

print('INDEPENDENT D25 PASS: 435/435 rational pairs, 406/406 parent pairs, Gram rank 25')
print('INDEPENDENT IDEAL GRAPH PASS: eigenvalue 15 x19, eigenvalue 29/2 x4')
print('CERTIFICATE_SHA256='+hashlib.sha256(p.read_bytes()).hexdigest())
