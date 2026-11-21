#!/usr/bin/env python3
"""Polybrutus D17: exact labelled Euclidean metric and chosen graph spectrum.

Frozen parent: Topbrutus/polybrutus-d16@v0.1.0,
commit c7700a46cb0b2c3d2556f50a207d86da84055ea0.
The new point is an orthogonal unit vector; all old points are unchanged.
"""
from pathlib import Path
from collections import Counter
from itertools import combinations
import hashlib,json,sys
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"prior"))
import verify_d16 as parent
Q=s.Rational
points=[p.col_join(s.zeros(1,1)) for p in parent.points]
points.append(s.Matrix([0]*17+[1]))
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
    assert len(parent.points)==21 and N==22
    assert all(p.shape==(18,1) for p in points)
    assert all(s.simplify(p.dot(p)-1)==0 for p in points)
    D=distances(points)
    assert all(D[i,j]>0 for i,j in pairs)
    old=json.loads((ROOT/"prior"/"certificate_d16.json").read_text(encoding="utf8"))
    assert old["status"]=="PASS" and old["dimension"]==16
    former={(int(e["i"])-1,int(e["j"])-1):s.Rational(e["d2"])
            for e in old["labelled_squared_distances"]}
    assert len(former)==210
    assert all(D[i,j]==value for (i,j),value in former.items())
    assert all(D[i,21]==2 for i in range(21))
    hist=Counter(4*D[i,j] for i,j in pairs)
    expected=Counter({s.Integer(4):6,s.Integer(7):2,s.Integer(8):200,
                       s.Integer(9):3,s.Integer(10):10,s.Integer(12):6,
                       s.Integer(13):2,s.Integer(16):2})
    assert hist==expected and sum(hist.values())==231
    assert affine_rank(points)==17
    B=centered_gram(D)
    assert B.rank()==17
    X=s.Matrix.hstack(*points).T
    H=s.eye(N)-s.ones(N)/N
    assert all(s.simplify(x)==0 for x in B-(H*X)*(H*X).T)
    dsum=s.simplify(sum(D[i,j] for i,j in pairs))
    assert dsum==Q(1895,4)
    center=sum(points,s.zeros(18,1))/N
    assert list(center)==[0,0,Q(1,44),0,0,0,0,0]+[Q(1,22)]*10

    # Four intentionally invalid controls: corrupted pair, duplicate, tilt, wrong labels.
    damaged=D.copy();damaged[0,21]=damaged[21,0]=Q(9,4)
    assert Counter(4*damaged[i,j] for i,j in pairs)!=expected
    duplicate=points[:21]+[points[20]]
    assert affine_rank(duplicate)==16
    tilt=points[:21]+[s.Matrix([Q(3,5)]+[0]*16+[Q(4,5)])]
    assert tilt[-1].dot(tilt[-1])==1
    assert any(distances(tilt)[i,21]!=2 for i in range(21))
    swapped=points.copy();swapped[0],swapped[7]=swapped[7],swapped[0]
    assert Counter(4*distances(swapped)[i,j] for i,j in pairs)==expected
    assert distances(swapped)[:7,:7]!=D[:7,:7]

    # Ideal equal-mass graph; physical experiments are not claimed.
    L16=laplacian(D[:21,:21])
    L17=laplacian(D)
    assert L17[:21,:21]-L16==s.eye(21)/2
    assert all(L17[i,21]==-Q(1,2) for i in range(21))
    assert L17[21,21]==Q(21,2)
    z=s.Symbol("z")
    chi16=s.factor(L16.charpoly(z).as_expr())
    chi17=s.factor(L17.charpoly(z).as_expr())
    recurrence=s.factor(z*(z-11)*s.cancel(chi16/z).subs(z,z-Q(1,2)))
    assert s.factor(chi17-recurrence)==0
    f6=s.cancel(s.factor(laplacian(D[:7,:7]).charpoly(z).as_expr())/z)
    closed=s.factor(z*(z-11)**11*(z-Q(21,2))**4*f6.subs(z,z-Q(15,2)))
    assert s.factor(chi17-closed)==0
    assert f6.subs(z,Q(7,2))!=0 and f6.subs(z,3)!=0

    report={
      "schema":"POLYBRUTUS-D17-EXACT-CERTIFICATE-v1",
      "status":"PASS","dimension":17,"vertices":22,"pairs":231,
      "affine_rank":17,"centered_gram_rank":17,"radius":1,
      "ambient_storage_coordinates":18,"intrinsic_affine_dimension":17,
      "parent":"Topbrutus/polybrutus-d16 tag v0.1.0",
      "parent_tag_commit":"c7700a46cb0b2c3d2556f50a207d86da84055ea0",
      "parent_doi":"10.5281/zenodo.24408979",
      "parent_distances_preserved":210,"new_distances":21,
      "new_squared_distances":["2"]*21,
      "histogram_k_4d2":{str(k):int(v) for k,v in sorted(hist.items())},
      "sum_squared_distances":str(dsum),
      "centroid":[str(x) for x in center],
      "gram_psd_certificate":"B=(HX)(HX)^T exact symbolic identity",
      "negative_controls":"4/4 PASS",
      "oscillator_model":"Hypothetical unit-mass weighted graph w_ij=1/d_ij²; no experimental physics",
      "spectral_recurrence":"chi22(z)=z*(z-11)*chi21(z-1/2)/(z-1/2)",
      "spectral_closed_form":"chi22(z)=z*(z-11)^11*(z-21/2)^4*f6(z-15/2), f6(z)=chiP7(z)/z",
      "eigenvalue_11_multiplicity":11,"eigenvalue_21_over_2_multiplicity":4,
      "spectral_charpoly":str(chi17),
      "novel_universal_theorem":"NOT CLAIMED","physical_discovery":"NOT CLAIMED",
      "human_peer_review":"NOT_DONE","sympy":s.__version__,
      "labelled_squared_distances":[{"i":i+1,"j":j+1,"d2":str(D[i,j]),"k":str(4*D[i,j])}
                                   for i,j in pairs]
    }
    out=ROOT/"results"/"certificate_d17.json"
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8",newline="\n")
    print("D17 EXACT PASS: 22 vertices, 231 pairs, rank 17, 210/210 old pairs preserved")
    print("MODEL GRAPH PASS: eigenvalue 11 x11; 21/2 x4; negative controls 4/4 PASS")
    print("CERTIFICATE_SHA256="+hashlib.sha256(out.read_bytes()).hexdigest())
    return report

if __name__=="__main__":
    run()
