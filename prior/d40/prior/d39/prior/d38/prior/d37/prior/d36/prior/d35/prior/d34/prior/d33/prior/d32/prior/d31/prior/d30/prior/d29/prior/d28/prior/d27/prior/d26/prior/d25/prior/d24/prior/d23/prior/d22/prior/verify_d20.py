#!/usr/bin/env python3
"""Polybrutus D20: 25-point exact labelled Euclidean metric and selected model graph.

Frozen predecessor Topbrutus/polybrutus-d19@v0.1.0,
source commit 1aa9721edcf6eb374a301863a5c1c4e566e69838.
Preserves all 276 predecessor labelled pairwise squared distances.
"""
from pathlib import Path
from collections import Counter
from itertools import combinations
import hashlib, json, sys
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"prior"))
import verify_d19 as parent
Q=s.Rational
points=[p.col_join(s.zeros(1,1)) for p in parent.points]
points.append(s.Matrix([0]*20+[1]))
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
    assert len(parent.points)==24 and N==25
    assert all(p.shape==(21,1) for p in points)
    assert all(s.simplify(p.dot(p)-1)==0 for p in points)
    D=distances(points)
    assert all(D[i,j]>0 for i,j in pairs)
    old=json.loads((ROOT/"prior"/"certificate_d19.json").read_text(encoding="utf8"))
    assert old["status"]=="PASS" and old["dimension"]==19
    former={(int(e["i"])-1,int(e["j"])-1):s.Rational(e["d2"])
            for e in old["labelled_squared_distances"]}
    assert len(former)==276
    assert all(D[i,j]==value for (i,j),value in former.items())
    assert all(D[i,24]==2 for i in range(24))
    hist=Counter(4*D[i,j] for i,j in pairs)
    expected=Counter({s.Integer(4):6,s.Integer(7):2,s.Integer(8):269,
                       s.Integer(9):3,s.Integer(10):10,s.Integer(12):6,
                       s.Integer(13):2,s.Integer(16):2})
    assert hist==expected and sum(hist.values())==300
    assert affine_rank(points)==20
    B=centered_gram(D)
    assert B.rank()==20
    X=s.Matrix.hstack(*points).T
    H=s.eye(N)-s.ones(N)/N
    assert all(s.simplify(v)==0 for v in B-(H*X)*(H*X).T)
    dsum=s.simplify(sum(D[i,j] for i,j in pairs))
    assert dsum==Q(2447,4)
    center=sum(points,s.zeros(21,1))/N
    assert list(center)==[0,0,Q(1,50),0,0,0,0,0]+[Q(1,25)]*13

    # Four falsification controls: tamper pair, duplicate, tilt, swap labels.
    damaged=D.copy();damaged[0,24]=damaged[24,0]=Q(9,4)
    assert Counter(4*damaged[i,j] for i,j in pairs)!=expected
    duplicate=points[:24]+[points[23]]
    assert affine_rank(duplicate)==19
    tilt=points[:24]+[s.Matrix([Q(3,5)]+[0]*19+[Q(4,5)])]
    assert tilt[-1].dot(tilt[-1])==1
    assert any(distances(tilt)[i,24]!=2 for i in range(24))
    swapped=points.copy();swapped[0],swapped[7]=swapped[7],swapped[0]
    assert Counter(4*distances(swapped)[i,j] for i,j in pairs)==expected
    assert distances(swapped)[:7,:7]!=D[:7,:7]

    # Chosen *theoretical* equal-mass complete graph; weights 1/d_ij^2.
    L19=laplacian(D[:24,:24])
    L20=laplacian(D)
    assert L20[:24,:24]-L19==s.eye(24)/2
    assert all(L20[i,24]==-Q(1,2) for i in range(24))
    assert L20[24,24]==12
    z=s.Symbol("z")
    chi19=s.factor(L19.charpoly(z).as_expr())
    chi20=s.factor(L20.charpoly(z).as_expr())
    recurrence=s.factor(z*(z-Q(25,2))*s.cancel(chi19/z).subs(z,z-Q(1,2)))
    assert s.factor(chi20-recurrence)==0
    f6=s.cancel(s.factor(laplacian(D[:7,:7]).charpoly(z).as_expr())/z)
    closed=s.factor(z*(z-Q(25,2))**14*(z-12)**4*f6.subs(z,z-9))
    assert s.factor(chi20-closed)==0
    assert f6.subs(z,Q(7,2))!=0 and f6.subs(z,3)!=0

    report={
      "schema":"POLYBRUTUS-D20-EXACT-CERTIFICATE-v1",
      "status":"PASS","dimension":20,"vertices":25,"pairs":300,
      "affine_rank":20,"centered_gram_rank":20,"radius":1,
      "ambient_storage_coordinates":21,"intrinsic_affine_dimension":20,
      "parent":"Topbrutus/polybrutus-d19 tag v0.1.0",
      "parent_tag_commit":"1aa9721edcf6eb374a301863a5c1c4e566e69838",
      "parent_doi":"10.5281/zenodo.24409343",
      "parent_distances_preserved":276,"new_distances":24,
      "new_squared_distances":["2"]*24,
      "histogram_k_4d2":{str(k):int(v) for k,v in sorted(hist.items())},
      "sum_squared_distances":str(dsum),
      "centroid":[str(v) for v in center],
      "gram_psd_certificate":"B=(HX)(HX)^T exact symbolic identity",
      "negative_controls":"4/4 PASS",
      "oscillator_model":"Selected ideal equal-mass weighted graph w_ij=1/d_ij²; no physical measurement",
      "spectral_recurrence":"chi25(z)=z*(z-25/2)*chi24(z-1/2)/(z-1/2)",
      "spectral_closed_form":"chi25(z)=z*(z-25/2)^14*(z-12)^4*f6(z-9), f6(z)=chiP7(z)/z",
      "eigenvalue_25_over_2_multiplicity":14,"eigenvalue_12_multiplicity":4,
      "spectral_charpoly":str(chi20),
      "novel_universal_theorem":"NOT CLAIMED","physical_discovery":"NOT CLAIMED",
      "human_peer_review":"NOT_DONE","sympy":s.__version__,
      "labelled_squared_distances":[{"i":i+1,"j":j+1,"d2":str(D[i,j]),"k":str(4*D[i,j])}
                                   for i,j in pairs]
    }
    out=ROOT/"results"/"certificate_d20.json"
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8",newline="\n")
    print("D20 EXACT PASS: 25 points, 300 labelled pairs, affine rank 20, 276/276 parent pairs")
    print("IDEAL GRAPH PASS: eigenvalue 25/2 x14; eigenvalue 12 x4; negative controls 4/4")
    print("CERTIFICATE_SHA256="+hashlib.sha256(out.read_bytes()).hexdigest())
    return report

if __name__=="__main__":
    run()
