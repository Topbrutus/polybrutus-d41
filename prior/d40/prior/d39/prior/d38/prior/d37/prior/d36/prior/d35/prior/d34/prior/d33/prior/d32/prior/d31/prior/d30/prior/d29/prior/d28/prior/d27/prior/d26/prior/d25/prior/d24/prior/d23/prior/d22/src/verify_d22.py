#!/usr/bin/env python3
"""Polybrutus D22: exact Euclidean labelled distance metric and chosen ideal graph.

Parent immutable source: Topbrutus/polybrutus-d21@v0.1.0
commit f17e4d02f7442d18cdafbcbb0b62380c557e4ef2.
No old vertex, labelled distance, source tag or predecessor DOI is changed.
"""
from pathlib import Path
from collections import Counter
from itertools import combinations
import hashlib,json,sys
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"prior"))
import verify_d21 as parent
Q=s.Rational
points=[p.col_join(s.zeros(1,1)) for p in parent.points]
points.append(s.Matrix([0]*22+[1]))
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
    assert len(parent.points)==26 and N==27
    assert all(p.shape==(23,1) for p in points)
    assert all(s.simplify(p.dot(p)-1)==0 for p in points)
    D=distances(points)
    assert all(D[i,j]>0 for i,j in pairs)
    old=json.loads((ROOT/"prior"/"certificate_d21.json").read_text(encoding="utf8"))
    assert old["status"]=="PASS" and old["dimension"]==21
    former={(int(e["i"])-1,int(e["j"])-1):s.Rational(e["d2"])
            for e in old["labelled_squared_distances"]}
    assert len(former)==325
    assert all(D[i,j]==value for (i,j),value in former.items())
    assert all(D[i,26]==2 for i in range(26))
    hist=Counter(4*D[i,j] for i,j in pairs)
    expected=Counter({s.Integer(4):6,s.Integer(7):2,s.Integer(8):320,
                       s.Integer(9):3,s.Integer(10):10,s.Integer(12):6,
                       s.Integer(13):2,s.Integer(16):2})
    assert hist==expected and sum(hist.values())==351
    assert affine_rank(points)==22
    B=centered_gram(D)
    assert B.rank()==22
    X=s.Matrix.hstack(*points).T
    H=s.eye(N)-s.ones(N)/N
    assert all(s.simplify(v)==0 for v in B-(H*X)*(H*X).T)
    dsum=s.simplify(sum(D[i,j] for i,j in pairs))
    assert dsum==Q(2855,4)
    center=sum(points,s.zeros(23,1))/N
    assert list(center)==[0,0,Q(1,54),0,0,0,0,0]+[Q(1,27)]*15

    # Detect invalidity: changed distance, duplicated point, tilted point, relabelling.
    damaged=D.copy();damaged[0,26]=damaged[26,0]=Q(9,4)
    assert Counter(4*damaged[i,j] for i,j in pairs)!=expected
    duplicate=points[:26]+[points[25]]
    assert affine_rank(duplicate)==21
    tilt=points[:26]+[s.Matrix([Q(3,5)]+[0]*21+[Q(4,5)])]
    assert tilt[-1].dot(tilt[-1])==1
    assert any(distances(tilt)[i,26]!=2 for i in range(26))
    swapped=points.copy();swapped[0],swapped[7]=swapped[7],swapped[0]
    assert Counter(4*distances(swapped)[i,j] for i,j in pairs)==expected
    assert distances(swapped)[:7,:7]!=D[:7,:7]

    # Deliberately stipulated graph: equal unit masses and weights w_ij=1/d_ij².
    L21=laplacian(D[:26,:26])
    L22=laplacian(D)
    assert L22[:26,:26]-L21==s.eye(26)/2
    assert all(L22[i,26]==-Q(1,2) for i in range(26))
    assert L22[26,26]==13
    z=s.Symbol("z")
    chi21=s.factor(L21.charpoly(z).as_expr())
    chi22=s.factor(L22.charpoly(z).as_expr())
    recurrence=s.factor(z*(z-Q(27,2))*s.cancel(chi21/z).subs(z,z-Q(1,2)))
    assert s.factor(chi22-recurrence)==0
    f6=s.cancel(s.factor(laplacian(D[:7,:7]).charpoly(z).as_expr())/z)
    closed=s.factor(z*(z-Q(27,2))**16*(z-13)**4*f6.subs(z,z-10))
    assert s.factor(chi22-closed)==0
    assert f6.subs(z,Q(7,2))!=0 and f6.subs(z,3)!=0

    report={
      "schema":"POLYBRUTUS-D22-EXACT-CERTIFICATE-v1",
      "status":"PASS","dimension":22,"vertices":27,"pairs":351,
      "affine_rank":22,"centered_gram_rank":22,"radius":1,
      "ambient_storage_coordinates":23,"intrinsic_affine_dimension":22,
      "parent":"Topbrutus/polybrutus-d21 tag v0.1.0",
      "parent_tag_commit":"f17e4d02f7442d18cdafbcbb0b62380c557e4ef2",
      "parent_doi":"10.5281/zenodo.24409527",
      "parent_distances_preserved":325,"new_distances":26,
      "new_squared_distances":["2"]*26,
      "histogram_k_4d2":{str(k):int(v) for k,v in sorted(hist.items())},
      "sum_squared_distances":str(dsum),
      "centroid":[str(v) for v in center],
      "gram_psd_certificate":"B=(HX)(HX)^T exact symbolic identity",
      "negative_controls":"4/4 PASS",
      "oscillator_model":"Chosen ideal equal-mass weighted graph w_ij=1/d_ij²; no physical measurement",
      "spectral_recurrence":"chi27(z)=z*(z-27/2)*chi26(z-1/2)/(z-1/2)",
      "spectral_closed_form":"chi27(z)=z*(z-27/2)^16*(z-13)^4*f6(z-10), f6(z)=chiP7(z)/z",
      "eigenvalue_27_over_2_multiplicity":16,"eigenvalue_13_multiplicity":4,
      "spectral_charpoly":str(chi22),
      "novel_universal_theorem":"NOT CLAIMED","physical_discovery":"NOT CLAIMED",
      "human_peer_review":"NOT_DONE","sympy":s.__version__,
      "labelled_squared_distances":[{"i":i+1,"j":j+1,"d2":str(D[i,j]),"k":str(4*D[i,j])}
                                   for i,j in pairs]
    }
    out=ROOT/"results"/"certificate_d22.json"
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8",newline="\n")
    print("D22 EXACT PASS: 27 points, 351 labelled pairs, affine rank 22, 325/325 predecessor pairs")
    print("IDEAL GRAPH PASS: eigenvalue 27/2 x16; eigenvalue 13 x4; negative controls 4/4")
    print("CERTIFICATE_SHA256="+hashlib.sha256(out.read_bytes()).hexdigest())
    return report

if __name__=="__main__":
    run()
