#!/usr/bin/env python3
"""Exact Polybrutus 8D certificate; no network or secrets.

Preserves every labelled edge of the publicly deposited 7D point configuration.
The 9-coordinate analytic representation has a single linear relation in the
five Helmert/simplex coordinates, so its intrinsic linear dimension is 8.
"""
import hashlib
import json
from collections import Counter
from itertools import combinations
from pathlib import Path

import sympy as s

ROOT=Path(__file__).resolve().parents[1]
P=s.Rational
rt3=s.sqrt(3)
base=[
    (0,0,1),
    (rt3/2,0,P(1,2)),
    (-rt3/4,P(3,4),P(1,2)),
    (-rt3/4,-P(3,4),P(1,2)),
    (0,0,-1),
    (rt3/2,0,-P(1,2)),
    (-rt3/2,0,-P(1,2)),
]
t=s.sqrt(5)/2
simplex=[tuple(t*(s.KroneckerDelta(i,j)-P(1,5)) for j in range(5)) for i in range(5)]
original=[s.Matrix([*r,*([0]*6)]) for r in base]
added=[s.Matrix([0,0,0,*r,0]) for r in simplex]
new=[s.Matrix([0,0,0,0,0,0,0,0,1])]
points=original+added+new


def gram_distance(pts):
    n=len(pts);d=s.zeros(n)
    for i,j in combinations(range(n),2):
        d[i,j]=d[j,i]=s.simplify((pts[i]-pts[j]).dot(pts[i]-pts[j]))
    return d


def distance_gram(d):
    n=d.rows; h=s.eye(n)-s.ones(n)/n
    return -h*d*h/2


def rank_coordinate(pts):
    return s.Matrix.hstack(*[z-pts[0] for z in pts[1:]]).rank()


def run():
    assert len(points)==13
    assert all(z.shape==(9,1) for z in points)
    assert all(s.simplify(z.dot(z)-1)==0 for z in points)
    d=gram_distance(points)
    d7=d[:12,:12]
    # Source-formula independent 7D metric check
    old_expected=Counter({s.Integer(4):6,s.Integer(7):2,s.Integer(8):35,s.Integer(9):3,s.Integer(10):10,s.Integer(12):6,s.Integer(13):2,s.Integer(16):2})
    hist7=Counter(s.simplify(4*d7[i,j]) for i,j in combinations(range(12),2))
    assert hist7==old_expected
    assert all(d[12,i]==2 for i in range(12))
    assert d[12,12]==0
    hist=Counter(s.simplify(4*d[i,j]) for i,j in combinations(range(13),2))
    expected=Counter({s.Integer(4):6,s.Integer(7):2,s.Integer(8):47,s.Integer(9):3,s.Integer(10):10,s.Integer(12):6,s.Integer(13):2,s.Integer(16):2})
    assert hist==expected and sum(hist.values())==78
    assert rank_coordinate(points)==8
    b=distance_gram(d)
    assert b.rank()==8
    # Explicit PSD factorization, not numeric eigenvalue heuristics.
    X=s.Matrix.hstack(*points).T
    H=s.eye(13)-s.ones(13)/13
    assert all(s.simplify(x)==0 for x in (b - (H*X)*(H*X).T))
    assert all(s.simplify(x)==0 for x in (d7-gram_distance(points[:12])))
    assert s.simplify(sum(d[i,j] for i,j in combinations(range(12),2))-P(575,4))==0
    # sum-square distances: earlier 575/4 plus 12 * 2.
    sumsq=s.simplify(sum(d[i,j] for i,j in combinations(range(13),2)))
    assert sumsq==P(671,4)
    centroid=sum(points,s.zeros(9,1))/13
    assert list(centroid)==[0,0,P(1,26),0,0,0,0,0,P(1,13)]
    # negative controls; all must be rejected by at least one independent invariant
    tamper=d.copy();tamper[0,12]=tamper[12,0]=P(9,4)
    assert Counter(s.simplify(4*tamper[i,j]) for i,j in combinations(range(13),2)) != expected
    duplicate=points[:12]+[points[7]]
    assert rank_coordinate(duplicate)==7
    bent=points[:12]+[s.Matrix([P(3,5),0,0,0,0,0,0,0,P(4,5)])]
    assert all(s.simplify(z.dot(z)-1)==0 for z in bent)
    assert any(gram_distance(bent)[12,i]!=2 for i in range(12))
    # label order countercontrol: histogram alone is not a certificate.
    swap=points.copy();swap[0],swap[7]=swap[7],swap[0]
    assert Counter(s.simplify(4*gram_distance(swap)[i,j]) for i,j in combinations(range(13),2))==expected
    assert gram_distance(swap)[:7,:7]!=d[:7,:7]

    # Weighted oscillator model: equal masses and illustrative coupling w=R^2/d^2.
    def lap(mat):
        n=mat.rows; L=s.zeros(n)
        for i,j in combinations(range(n),2):
            w=1/mat[i,j];L[i,i]+=w;L[j,j]+=w;L[i,j]-=w;L[j,i]-=w
        return L
    L7=lap(d7);L8=lap(d)
    delta=L8[:12,:12]-L7
    assert delta==s.eye(12)/2
    assert all(L8[i,12]==-s.Rational(1,2) for i in range(12))
    assert L8[12,12]==6
    z=s.symbols('z')
    poly7=s.factor(L7.charpoly(z).as_expr());poly8=s.factor(L8.charpoly(z).as_expr())
    expected8=s.factor(z*(z-s.Rational(13,2))*s.cancel(poly7/z).subs(z,z-s.Rational(1,2)))
    assert s.factor(poly8-expected8)==0
    # full exact factorization relative to canonical P7 Laplacian:
    LP=lap(d[:7,:7]);polyP=s.factor(LP.charpoly(z).as_expr());f6=s.cancel(polyP/z)
    spectral=s.factor(z*(z-s.Rational(13,2))**2*(z-6)**4*f6.subs(z,z-3))
    assert s.factor(poly8-spectral)==0
    assert poly8.subs(z,s.Rational(13,2))==0
    # eigen mult 13/2 exactly 2 from charpoly, unless f6 has root 7/2:
    assert s.simplify(f6.subs(z,s.Rational(7,2)))!=0
    assert s.simplify(f6.subs(z,3))!=0

    report={
      'status':'PASS','dimension':8,'vertices':13,'new_vertices_since_7d':1,'pair_count':78,
      'previous_distance_count_preserved':66,'new_pairs':12,'previous_affine_rank':7,
      'affine_rank':8,'centered_gram_rank':8,'centered_gram_psd_proof':'B=(H*X)(H*X)^T, exact symbolic equality',
      'radius':1,'ambient_storage_coordinates':9,'affine_dimension':8,
      'histogram_k_4d2':{str(k):int(v) for k,v in sorted(hist.items())},
      'sum_squared_distances':str(sumsq),
      'centroid':list(map(str,centroid)),
      'oscillator_coupling':'w_ij=1/d_ij^2; hypothetical simulation model, not a real measurement',
      'weighted_laplacian_charpoly_8d':str(poly8),
      'spectral_recursion':'chi13(z)=z(z-13/2) * (chi12(z-1/2)/(z-1/2))',
      'spectral_closed_form':'chi13(z)=z(z-13/2)^2(z-6)^4 f6(z-3), f6(z)=chi_P7(z)/z',
      'lambda_13_over_2_multiplicity':2,
      'negative_controls':'4/4 PASS',
      'novelty':'General orthogonal lifting and graph-cone spectral shift are known; no original universal theorem claimed',
      'physical_discovery':'NOT CLAIMED',
      'sympy':s.__version__,
      'pairwise_distances':[{'i':i+1,'j':j+1,'d2':str(d[i,j]),'k':str(4*d[i,j])} for i,j in combinations(range(13),2)]
    }
    out=ROOT/'results'/'certificate_d08.json'
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf8',newline='\n')
    print('D08 CERTIFICATE PASS; N=13 L=78 RANK=8 GRAM_PSD=EXACT; spectral mode lambda=13/2 x2; negative controls=PASS')
    print('Sympy:',s.__version__,'SHA256:',hashlib.sha256(out.read_bytes()).hexdigest())
    return report

if __name__=='__main__':run()
