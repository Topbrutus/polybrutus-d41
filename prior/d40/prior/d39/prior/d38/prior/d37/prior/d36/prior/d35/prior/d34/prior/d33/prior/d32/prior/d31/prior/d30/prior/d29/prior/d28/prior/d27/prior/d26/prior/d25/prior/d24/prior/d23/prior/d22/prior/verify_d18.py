#!/usr/bin/env python3
"""Polybrutus D18: exact 23-vertex Euclidean metric, affine rank and graph model.

Frozen predecessor: Topbrutus/polybrutus-d17@v0.1.0
commit 1ef37f156292aa80f8c3f5d07bb889dc4b4492ab.
No old vertex, distance, frozen source or previous DOI is changed.
"""
from pathlib import Path
from collections import Counter
from itertools import combinations
import hashlib,json,sys
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"prior"))
import verify_d17 as parent
Q=s.Rational
points=[p.col_join(s.zeros(1,1)) for p in parent.points]
points.append(s.Matrix([0]*18+[1]))
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
    assert len(parent.points)==22 and N==23
    assert all(p.shape==(19,1) for p in points)
    assert all(s.simplify(p.dot(p)-1)==0 for p in points)
    D=distances(points)
    assert all(D[i,j]>0 for i,j in pairs)
    old=json.loads((ROOT/"prior"/"certificate_d17.json").read_text(encoding="utf8"))
    assert old["status"]=="PASS" and old["dimension"]==17
    former={(int(e["i"])-1,int(e["j"])-1):s.Rational(e["d2"])
            for e in old["labelled_squared_distances"]}
    assert len(former)==231
    assert all(D[i,j]==value for (i,j),value in former.items())
    assert all(D[i,22]==2 for i in range(22))
    hist=Counter(4*D[i,j] for i,j in pairs)
    expected=Counter({s.Integer(4):6,s.Integer(7):2,s.Integer(8):222,
                       s.Integer(9):3,s.Integer(10):10,s.Integer(12):6,
                       s.Integer(13):2,s.Integer(16):2})
    assert hist==expected and sum(hist.values())==253
    assert affine_rank(points)==18
    B=centered_gram(D)
    assert B.rank()==18
    X=s.Matrix.hstack(*points).T
    H=s.eye(N)-s.ones(N)/N
    assert all(s.simplify(v)==0 for v in B-(H*X)*(H*X).T)
    dsum=s.simplify(sum(D[i,j] for i,j in pairs))
    assert dsum==Q(2071,4)
    center=sum(points,s.zeros(19,1))/N
    assert list(center)==[0,0,Q(1,46),0,0,0,0,0]+[Q(1,23)]*11

    # Four invalid controls: alter distance, duplicate, tilt, swap labels.
    damaged=D.copy();damaged[0,22]=damaged[22,0]=Q(9,4)
    assert Counter(4*damaged[i,j] for i,j in pairs)!=expected
    duplicate=points[:22]+[points[21]]
    assert affine_rank(duplicate)==17
    tilt=points[:22]+[s.Matrix([Q(3,5)]+[0]*17+[Q(4,5)])]
    assert tilt[-1].dot(tilt[-1])==1
    assert any(distances(tilt)[i,22]!=2 for i in range(22))
    swapped=points.copy();swapped[0],swapped[7]=swapped[7],swapped[0]
    assert Counter(4*distances(swapped)[i,j] for i,j in pairs)==expected
    assert distances(swapped)[:7,:7]!=D[:7,:7]

    # Chosen theoretical graph: unit masses and w_ij = 1/d_ij².
    L17=laplacian(D[:22,:22])
    L18=laplacian(D)
    assert L18[:22,:22]-L17==s.eye(22)/2
    assert all(L18[i,22]==-Q(1,2) for i in range(22))
    assert L18[22,22]==11
    z=s.Symbol("z")
    chi17=s.factor(L17.charpoly(z).as_expr())
    chi18=s.factor(L18.charpoly(z).as_expr())
    recurrence=s.factor(z*(z-Q(23,2))*s.cancel(chi17/z).subs(z,z-Q(1,2)))
    assert s.factor(chi18-recurrence)==0
    f6=s.cancel(s.factor(laplacian(D[:7,:7]).charpoly(z).as_expr())/z)
    closed=s.factor(z*(z-Q(23,2))**12*(z-11)**4*f6.subs(z,z-8))
    assert s.factor(chi18-closed)==0
    assert f6.subs(z,Q(7,2))!=0 and f6.subs(z,3)!=0

    report={
      "schema":"POLYBRUTUS-D18-EXACT-CERTIFICATE-v1",
      "status":"PASS","dimension":18,"vertices":23,"pairs":253,
      "affine_rank":18,"centered_gram_rank":18,"radius":1,
      "ambient_storage_coordinates":19,"intrinsic_affine_dimension":18,
      "parent":"Topbrutus/polybrutus-d17 tag v0.1.0",
      "parent_tag_commit":"1ef37f156292aa80f8c3f5d07bb889dc4b4492ab",
      "parent_doi":"10.5281/zenodo.24409067",
      "parent_distances_preserved":231,"new_distances":22,
      "new_squared_distances":["2"]*22,
      "histogram_k_4d2":{str(k):int(v) for k,v in sorted(hist.items())},
      "sum_squared_distances":str(dsum),
      "centroid":[str(v) for v in center],
      "gram_psd_certificate":"B=(HX)(HX)^T exact symbolic identity",
      "negative_controls":"4/4 PASS",
      "oscillator_model":"Chosen ideal equal-mass weighted graph w_ij=1/d_ij²; no physical measurement",
      "spectral_recurrence":"chi23(z)=z*(z-23/2)*chi22(z-1/2)/(z-1/2)",
      "spectral_closed_form":"chi23(z)=z*(z-23/2)^12*(z-11)^4*f6(z-8), f6(z)=chiP7(z)/z",
      "eigenvalue_23_over_2_multiplicity":12,"eigenvalue_11_multiplicity":4,
      "spectral_charpoly":str(chi18),
      "novel_universal_theorem":"NOT CLAIMED","physical_discovery":"NOT CLAIMED",
      "human_peer_review":"NOT_DONE","sympy":s.__version__,
      "labelled_squared_distances":[{"i":i+1,"j":j+1,"d2":str(D[i,j]),"k":str(4*D[i,j])}
                                   for i,j in pairs]
    }
    out=ROOT/"results"/"certificate_d18.json"
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8",newline="\n")
    print("D18 EXACT PASS: 23 points, 253 labelled pairs, affine rank 18; 231/231 parent pairs")
    print("IDEAL GRAPH PASS: eigenvalue 23/2 x12; eigenvalue 11 x4; 4/4 negative controls")
    print("CERTIFICATE_SHA256="+hashlib.sha256(out.read_bytes()).hexdigest())
    return report

if __name__=="__main__":
    run()
