#!/usr/bin/env python3
"""Polybrutus D13: exact affine metric and illustrative graph spectrum.

Frozen source: Topbrutus/polybrutus-d12@v0.1.0,
commit e6ea0d1d174bcd9580e4d1d1c5c5dd2355cc935e.
Adds one orthogonal unit vector. Four negative controls and exact
symbolic Gram/Laplacian certificates. No physical-dimension claim.
"""
from pathlib import Path
from collections import Counter
from itertools import combinations
import hashlib,json,sys
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"prior"))
import verify_d12 as parent

Q=s.Rational
points=[p.col_join(s.zeros(1,1)) for p in parent.points]
points.append(s.Matrix([0]*13+[1]))
N=len(points)
pairs=list(combinations(range(N),2))

def distances(P):
    D=s.zeros(len(P))
    for i,j in combinations(range(len(P)),2):
        D[i,j]=D[j,i]=s.simplify((P[i]-P[j]).dot(P[i]-P[j]))
    return D

def affine_rank(P):
    return s.Matrix.hstack(*[p-P[0] for p in P[1:]]).rank()

def centered_gram(D):
    H=s.eye(D.rows)-s.ones(D.rows)/D.rows
    return -H*D*H/2

def laplacian(D):
    L=s.zeros(D.rows)
    for i,j in combinations(range(D.rows),2):
        w=1/D[i,j]
        L[i,i]+=w;L[j,j]+=w
        L[i,j]-=w;L[j,i]-=w
    return L

def run():
    assert len(parent.points)==17 and N==18
    assert all(p.shape==(14,1) for p in points)
    assert all(s.simplify(p.dot(p)-1)==0 for p in points)
    D=distances(points)
    assert all(D[i,j]>0 for i,j in pairs)
    old=json.loads((ROOT/"prior"/"certificate_d12.json").read_text(encoding="utf8"))
    assert old["status"]=="PASS" and old["dimension"]==12
    old_pairs={(int(e["i"])-1,int(e["j"])-1):s.Rational(e["d2"])
               for e in old["labelled_squared_distances"]}
    assert len(old_pairs)==136
    assert all(D[i,j]==v for (i,j),v in old_pairs.items())
    assert all(D[i,17]==2 for i in range(17))
    hist=Counter(4*D[i,j] for i,j in pairs)
    expected=Counter({s.Integer(4):6,s.Integer(7):2,s.Integer(8):122,
                       s.Integer(9):3,s.Integer(10):10,s.Integer(12):6,
                       s.Integer(13):2,s.Integer(16):2})
    assert hist==expected and sum(hist.values())==153
    assert affine_rank(points)==13
    B=centered_gram(D)
    assert B.rank()==13
    X=s.Matrix.hstack(*points).T
    H=s.eye(N)-s.ones(N)/N
    assert all(s.simplify(z)==0 for z in B-(H*X)*(H*X).T)
    dsum=s.simplify(sum(D[i,j] for i,j in pairs))
    assert dsum==Q(1271,4)
    center=sum(points,s.zeros(14,1))/N
    assert list(center)==[0,0,Q(1,36),0,0,0,0,0,Q(1,18),Q(1,18),
                          Q(1,18),Q(1,18),Q(1,18),Q(1,18)]

    # Four bad cases must be rejected by the indicated check.
    corrupt=D.copy();corrupt[0,17]=corrupt[17,0]=Q(9,4)
    assert Counter(4*corrupt[i,j] for i,j in pairs)!=expected
    duplicate=points[:17]+[points[16]]
    assert affine_rank(duplicate)==12
    tilted=points[:17]+[s.Matrix([Q(3,5),0,0,0,0,0,0,0,0,0,0,0,0,Q(4,5)])]
    assert tilted[-1].dot(tilted[-1])==1
    assert any(distances(tilted)[i,17]!=2 for i in range(17))
    permuted=points.copy();permuted[0],permuted[7]=permuted[7],permuted[0]
    assert Counter(4*distances(permuted)[i,j] for i,j in pairs)==expected
    assert distances(permuted)[:7,:7]!=D[:7,:7]

    # Idealized equal-mass complete weighted-graph model, not measurements.
    L12=laplacian(D[:17,:17])
    L13=laplacian(D)
    assert L13[:17,:17]-L12==s.eye(17)/2
    assert all(L13[i,17]==-Q(1,2) for i in range(17))
    assert L13[17,17]==Q(17,2)
    z=s.Symbol("z")
    chi12=s.factor(L12.charpoly(z).as_expr())
    chi13=s.factor(L13.charpoly(z).as_expr())
    recursive=s.factor(z*(z-9)*s.cancel(chi12/z).subs(z,z-Q(1,2)))
    assert s.factor(chi13-recursive)==0
    f6=s.cancel(s.factor(laplacian(D[:7,:7]).charpoly(z).as_expr())/z)
    closed=s.factor(z*(z-9)**7*(z-Q(17,2))**4*f6.subs(z,z-Q(11,2)))
    assert s.factor(chi13-closed)==0
    assert f6.subs(z,Q(7,2))!=0 and f6.subs(z,3)!=0

    report={
      "schema":"POLYBRUTUS-D13-EXACT-CERTIFICATE-v1",
      "status":"PASS","dimension":13,"vertices":18,"pairs":153,
      "affine_rank":13,"centered_gram_rank":13,"radius":1,
      "ambient_storage_coordinates":14,"intrinsic_affine_dimension":13,
      "parent":"Topbrutus/polybrutus-d12 tag v0.1.0",
      "parent_tag_commit":"e6ea0d1d174bcd9580e4d1d1c5c5dd2355cc935e",
      "parent_doi":"10.5281/zenodo.24408720",
      "parent_distances_preserved":136,"new_distances":17,
      "new_squared_distances":["2"]*17,
      "histogram_k_4d2":{str(k):int(v) for k,v in sorted(hist.items())},
      "sum_squared_distances":str(dsum),
      "centroid":[str(v) for v in center],
      "gram_psd_certificate":"B=(HX)(HX)^T exact symbolic identity",
      "negative_controls":"4/4 PASS",
      "oscillator_model":"Ideal equal-mass complete weighted graph with w_ij=1/d_ij²; no physical observation",
      "spectral_recurrence":"chi18(z)=z*(z-9)*chi17(z-1/2)/(z-1/2)",
      "spectral_closed_form":"chi18(z)=z*(z-9)^7*(z-17/2)^4*f6(z-11/2), f6(z)=chiP7(z)/z",
      "eigenvalue_9_multiplicity":7,
      "eigenvalue_17_over_2_multiplicity":4,
      "spectral_charpoly":str(chi13),
      "historical_base_13":"The dimension label 13 does not establish the historical base-13, 13R, or Brutus boson conjectures",
      "novel_universal_theorem":"NOT CLAIMED","physical_discovery":"NOT CLAIMED",
      "human_peer_review":"NOT_DONE","sympy":s.__version__,
      "labelled_squared_distances":[{"i":i+1,"j":j+1,"d2":str(D[i,j]),"k":str(4*D[i,j])}
                                   for i,j in pairs]
    }
    out=ROOT/"results"/"certificate_d13.json"
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8",newline="\n")
    print("D13 EXACT PASS: 18 points, 153 pairs, affine rank=13, PSD; 136/136 old pairs preserved")
    print("IDEAL GRAPH PASS: eigenvalue 9 multiplicity 7, eigenvalue 17/2 multiplicity 4; negative controls 4/4")
    print("CERTIFICATE_SHA256="+hashlib.sha256(out.read_bytes()).hexdigest())
    return report

if __name__=="__main__":
    run()
