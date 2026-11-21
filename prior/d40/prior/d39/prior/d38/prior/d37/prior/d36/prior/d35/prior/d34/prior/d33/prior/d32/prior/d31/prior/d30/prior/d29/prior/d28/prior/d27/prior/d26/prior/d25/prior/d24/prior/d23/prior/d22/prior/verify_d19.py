#!/usr/bin/env python3
"""Polybrutus D19 exact Euclidean labelled metric and chosen ideal graph model.

Frozen predecessor: Topbrutus/polybrutus-d18@v0.1.0
commit b42aae51d92b817008c823a34162e1eb94c13cae.
Conservatively extends all 253 labelled old distances without modification.
"""
from pathlib import Path
from collections import Counter
from itertools import combinations
import hashlib,json,sys
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"prior"))
import verify_d18 as parent
Q=s.Rational
points=[p.col_join(s.zeros(1,1)) for p in parent.points]
points.append(s.Matrix([0]*19+[1]))
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
    assert len(parent.points)==23 and N==24
    assert all(p.shape==(20,1) for p in points)
    assert all(s.simplify(p.dot(p)-1)==0 for p in points)
    D=distances(points)
    assert all(D[i,j]>0 for i,j in pairs)
    old=json.loads((ROOT/"prior"/"certificate_d18.json").read_text(encoding="utf8"))
    assert old["status"]=="PASS" and old["dimension"]==18
    former={(int(e["i"])-1,int(e["j"])-1):s.Rational(e["d2"])
            for e in old["labelled_squared_distances"]}
    assert len(former)==253
    assert all(D[i,j]==value for (i,j),value in former.items())
    assert all(D[i,23]==2 for i in range(23))
    hist=Counter(4*D[i,j] for i,j in pairs)
    expected=Counter({s.Integer(4):6,s.Integer(7):2,s.Integer(8):245,
                       s.Integer(9):3,s.Integer(10):10,s.Integer(12):6,
                       s.Integer(13):2,s.Integer(16):2})
    assert hist==expected and sum(hist.values())==276
    assert affine_rank(points)==19
    B=centered_gram(D)
    assert B.rank()==19
    X=s.Matrix.hstack(*points).T
    H=s.eye(N)-s.ones(N)/N
    assert all(s.simplify(v)==0 for v in B-(H*X)*(H*X).T)
    dsum=s.simplify(sum(D[i,j] for i,j in pairs))
    assert dsum==Q(2255,4)
    center=sum(points,s.zeros(20,1))/N
    assert list(center)==[0,0,Q(1,48),0,0,0,0,0]+[Q(1,24)]*12

    # Four invalid controls: altered distance, duplicate point, tilt, wrong labels.
    damaged=D.copy();damaged[0,23]=damaged[23,0]=Q(9,4)
    assert Counter(4*damaged[i,j] for i,j in pairs)!=expected
    duplicate=points[:23]+[points[22]]
    assert affine_rank(duplicate)==18
    tilt=points[:23]+[s.Matrix([Q(3,5)]+[0]*18+[Q(4,5)])]
    assert tilt[-1].dot(tilt[-1])==1
    assert any(distances(tilt)[i,23]!=2 for i in range(23))
    swapped=points.copy();swapped[0],swapped[7]=swapped[7],swapped[0]
    assert Counter(4*distances(swapped)[i,j] for i,j in pairs)==expected
    assert distances(swapped)[:7,:7]!=D[:7,:7]

    # Theoretical complete graph with unit masses and w_ij=1/d_ij².
    L18=laplacian(D[:23,:23])
    L19=laplacian(D)
    assert L19[:23,:23]-L18==s.eye(23)/2
    assert all(L19[i,23]==-Q(1,2) for i in range(23))
    assert L19[23,23]==Q(23,2)
    z=s.Symbol("z")
    chi18=s.factor(L18.charpoly(z).as_expr())
    chi19=s.factor(L19.charpoly(z).as_expr())
    recurrence=s.factor(z*(z-12)*s.cancel(chi18/z).subs(z,z-Q(1,2)))
    assert s.factor(chi19-recurrence)==0
    f6=s.cancel(s.factor(laplacian(D[:7,:7]).charpoly(z).as_expr())/z)
    closed=s.factor(z*(z-12)**13*(z-Q(23,2))**4*f6.subs(z,z-Q(17,2)))
    assert s.factor(chi19-closed)==0
    assert f6.subs(z,Q(7,2))!=0 and f6.subs(z,3)!=0

    report={
      "schema":"POLYBRUTUS-D19-EXACT-CERTIFICATE-v1",
      "status":"PASS","dimension":19,"vertices":24,"pairs":276,
      "affine_rank":19,"centered_gram_rank":19,"radius":1,
      "ambient_storage_coordinates":20,"intrinsic_affine_dimension":19,
      "parent":"Topbrutus/polybrutus-d18 tag v0.1.0",
      "parent_tag_commit":"b42aae51d92b817008c823a34162e1eb94c13cae",
      "parent_doi":"10.5281/zenodo.24409207",
      "parent_distances_preserved":253,"new_distances":23,
      "new_squared_distances":["2"]*23,
      "histogram_k_4d2":{str(k):int(v) for k,v in sorted(hist.items())},
      "sum_squared_distances":str(dsum),
      "centroid":[str(v) for v in center],
      "gram_psd_certificate":"B=(HX)(HX)^T exact symbolic identity",
      "negative_controls":"4/4 PASS",
      "oscillator_model":"Chosen ideal unit-mass weighted graph w_ij=1/d_ij²; no physical measurement",
      "spectral_recurrence":"chi24(z)=z*(z-12)*chi23(z-1/2)/(z-1/2)",
      "spectral_closed_form":"chi24(z)=z*(z-12)^13*(z-23/2)^4*f6(z-17/2), f6(z)=chiP7(z)/z",
      "eigenvalue_12_multiplicity":13,"eigenvalue_23_over_2_multiplicity":4,
      "spectral_charpoly":str(chi19),
      "novel_universal_theorem":"NOT CLAIMED","physical_discovery":"NOT CLAIMED",
      "human_peer_review":"NOT_DONE","sympy":s.__version__,
      "labelled_squared_distances":[{"i":i+1,"j":j+1,"d2":str(D[i,j]),"k":str(4*D[i,j])}
                                   for i,j in pairs]
    }
    out=ROOT/"results"/"certificate_d19.json"
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8",newline="\n")
    print("D19 EXACT PASS: 24 points, 276 labelled pairs, affine rank 19; 253/253 parent pairs")
    print("IDEAL GRAPH PASS: eigenvalue 12 x13; eigenvalue 23/2 x4; 4/4 negative controls")
    print("CERTIFICATE_SHA256="+hashlib.sha256(out.read_bytes()).hexdigest())
    return report

if __name__=="__main__":
    run()
