#!/usr/bin/env python3
"""Polybrutus D21 exact labelled Euclidean metric and selected graph model.

Frozen predecessor: Topbrutus/polybrutus-d20@v0.1.0,
commit 3e9ee3117f1c3c3e5ffc80236fa2026a5d825ce3.
This case preserves all 300 predecessor labelled squared distances.
"""
from pathlib import Path
from collections import Counter
from itertools import combinations
import hashlib, json, sys
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"prior"))
import verify_d20 as parent
Q=s.Rational
points=[p.col_join(s.zeros(1,1)) for p in parent.points]
points.append(s.Matrix([0]*21+[1]))
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
    assert len(parent.points)==25 and N==26
    assert all(p.shape==(22,1) for p in points)
    assert all(s.simplify(p.dot(p)-1)==0 for p in points)
    D=distances(points)
    assert all(D[i,j]>0 for i,j in pairs)
    old=json.loads((ROOT/"prior"/"certificate_d20.json").read_text(encoding="utf8"))
    assert old["status"]=="PASS" and old["dimension"]==20
    former={(int(e["i"])-1,int(e["j"])-1):s.Rational(e["d2"])
            for e in old["labelled_squared_distances"]}
    assert len(former)==300
    assert all(D[i,j]==value for (i,j),value in former.items())
    assert all(D[i,25]==2 for i in range(25))
    hist=Counter(4*D[i,j] for i,j in pairs)
    expected=Counter({s.Integer(4):6,s.Integer(7):2,s.Integer(8):294,
                       s.Integer(9):3,s.Integer(10):10,s.Integer(12):6,
                       s.Integer(13):2,s.Integer(16):2})
    assert hist==expected and sum(hist.values())==325
    assert affine_rank(points)==21
    B=centered_gram(D)
    assert B.rank()==21
    X=s.Matrix.hstack(*points).T
    H=s.eye(N)-s.ones(N)/N
    assert all(s.simplify(v)==0 for v in B-(H*X)*(H*X).T)
    dsum=s.simplify(sum(D[i,j] for i,j in pairs))
    assert dsum==Q(2647,4)
    center=sum(points,s.zeros(22,1))/N
    assert list(center)==[0,0,Q(1,52),0,0,0,0,0]+[Q(1,26)]*14

    # Four negative controls: tampered metric, duplicate, tilt, re-labelling.
    damaged=D.copy();damaged[0,25]=damaged[25,0]=Q(9,4)
    assert Counter(4*damaged[i,j] for i,j in pairs)!=expected
    duplicate=points[:25]+[points[24]]
    assert affine_rank(duplicate)==20
    tilt=points[:25]+[s.Matrix([Q(3,5)]+[0]*20+[Q(4,5)])]
    assert tilt[-1].dot(tilt[-1])==1
    assert any(distances(tilt)[i,25]!=2 for i in range(25))
    swapped=points.copy();swapped[0],swapped[7]=swapped[7],swapped[0]
    assert Counter(4*distances(swapped)[i,j] for i,j in pairs)==expected
    assert distances(swapped)[:7,:7]!=D[:7,:7]

    # Deliberately selected ideal graph with unit masses and w_ij = 1/d_ij².
    L20=laplacian(D[:25,:25])
    L21=laplacian(D)
    assert L21[:25,:25]-L20==s.eye(25)/2
    assert all(L21[i,25]==-Q(1,2) for i in range(25))
    assert L21[25,25]==Q(25,2)
    z=s.Symbol("z")
    chi20=s.factor(L20.charpoly(z).as_expr())
    chi21=s.factor(L21.charpoly(z).as_expr())
    recurrence=s.factor(z*(z-13)*s.cancel(chi20/z).subs(z,z-Q(1,2)))
    assert s.factor(chi21-recurrence)==0
    f6=s.cancel(s.factor(laplacian(D[:7,:7]).charpoly(z).as_expr())/z)
    closed=s.factor(z*(z-13)**15*(z-Q(25,2))**4*f6.subs(z,z-Q(19,2)))
    assert s.factor(chi21-closed)==0
    assert f6.subs(z,Q(7,2))!=0 and f6.subs(z,3)!=0

    report={
      "schema":"POLYBRUTUS-D21-EXACT-CERTIFICATE-v1",
      "status":"PASS","dimension":21,"vertices":26,"pairs":325,
      "affine_rank":21,"centered_gram_rank":21,"radius":1,
      "ambient_storage_coordinates":22,"intrinsic_affine_dimension":21,
      "parent":"Topbrutus/polybrutus-d20 tag v0.1.0",
      "parent_tag_commit":"3e9ee3117f1c3c3e5ffc80236fa2026a5d825ce3",
      "parent_doi":"10.5281/zenodo.24409424",
      "parent_distances_preserved":300,"new_distances":25,
      "new_squared_distances":["2"]*25,
      "histogram_k_4d2":{str(k):int(v) for k,v in sorted(hist.items())},
      "sum_squared_distances":str(dsum),
      "centroid":[str(v) for v in center],
      "gram_psd_certificate":"B=(HX)(HX)^T exact symbolic identity",
      "negative_controls":"4/4 PASS",
      "oscillator_model":"Selected ideal equal-mass weighted graph w_ij=1/d_ij²; no physical measurement",
      "spectral_recurrence":"chi26(z)=z*(z-13)*chi25(z-1/2)/(z-1/2)",
      "spectral_closed_form":"chi26(z)=z*(z-13)^15*(z-25/2)^4*f6(z-19/2), f6(z)=chiP7(z)/z",
      "eigenvalue_13_multiplicity":15,"eigenvalue_25_over_2_multiplicity":4,
      "spectral_charpoly":str(chi21),
      "novel_universal_theorem":"NOT CLAIMED","physical_discovery":"NOT CLAIMED",
      "human_peer_review":"NOT_DONE","sympy":s.__version__,
      "labelled_squared_distances":[{"i":i+1,"j":j+1,"d2":str(D[i,j]),"k":str(4*D[i,j])}
                                   for i,j in pairs]
    }
    out=ROOT/"results"/"certificate_d21.json"
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8",newline="\n")
    print("D21 EXACT PASS: 26 points, 325 labelled pairs, affine rank 21, 300/300 parent pairs")
    print("IDEAL GRAPH PASS: eigenvalue 13 x15; eigenvalue 25/2 x4; negative controls 4/4")
    print("CERTIFICATE_SHA256="+hashlib.sha256(out.read_bytes()).hexdigest())
    return report

if __name__=="__main__":
    run()
