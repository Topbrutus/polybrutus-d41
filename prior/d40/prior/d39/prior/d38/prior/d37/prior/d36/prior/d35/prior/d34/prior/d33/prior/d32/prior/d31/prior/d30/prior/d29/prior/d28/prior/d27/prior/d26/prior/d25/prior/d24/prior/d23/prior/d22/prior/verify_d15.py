#!/usr/bin/env python3
"""Polybrutus D15 exact labelled Euclidean metric and ideal graph certificate.

Frozen parent: Topbrutus/polybrutus-d14@v0.1.0,
commit 853687b71dd32b05c807aa47b1c8d7b48d7e97e4.
All previous labelled points and distances remain unaltered.
"""
from pathlib import Path
from collections import Counter
from itertools import combinations
import hashlib,json,sys
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"prior"))
import verify_d14 as parent
Q=s.Rational
points=[p.col_join(s.zeros(1,1)) for p in parent.points]
points.append(s.Matrix([0]*15+[1]))
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
    assert len(parent.points)==19 and N==20
    assert all(p.shape==(16,1) for p in points)
    assert all(s.simplify(p.dot(p)-1)==0 for p in points)
    D=distances(points)
    assert all(D[i,j]>0 for i,j in pairs)
    old=json.loads((ROOT/"prior"/"certificate_d14.json").read_text(encoding="utf8"))
    assert old["status"]=="PASS" and old["dimension"]==14
    former={(int(e["i"])-1,int(e["j"])-1):s.Rational(e["d2"])
            for e in old["labelled_squared_distances"]}
    assert len(former)==171
    assert all(D[i,j]==value for (i,j),value in former.items())
    assert all(D[i,19]==2 for i in range(19))
    hist=Counter(4*D[i,j] for i,j in pairs)
    expected=Counter({s.Integer(4):6,s.Integer(7):2,s.Integer(8):159,
                       s.Integer(9):3,s.Integer(10):10,s.Integer(12):6,
                       s.Integer(13):2,s.Integer(16):2})
    assert hist==expected and sum(hist.values())==190
    assert affine_rank(points)==15
    B=centered_gram(D)
    assert B.rank()==15
    X=s.Matrix.hstack(*points).T
    H=s.eye(N)-s.ones(N)/N
    assert all(s.simplify(x)==0 for x in B-(H*X)*(H*X).T)
    dsum=s.simplify(sum(D[i,j] for i,j in pairs))
    assert dsum==Q(1567,4)
    center=sum(points,s.zeros(16,1))/N
    assert list(center)==[0,0,Q(1,40),0,0,0,0,0]+[Q(1,20)]*8

    # Four failing examples: mutated distance, duplicate, tilt and wrong labels.
    damaged=D.copy();damaged[0,19]=damaged[19,0]=Q(9,4)
    assert Counter(4*damaged[i,j] for i,j in pairs)!=expected
    duplicate=points[:19]+[points[18]]
    assert affine_rank(duplicate)==14
    tilt=points[:19]+[s.Matrix([Q(3,5)]+[0]*14+[Q(4,5)])]
    assert tilt[-1].dot(tilt[-1])==1
    assert any(distances(tilt)[i,19]!=2 for i in range(19))
    swapped=points.copy();swapped[0],swapped[7]=swapped[7],swapped[0]
    assert Counter(4*distances(swapped)[i,j] for i,j in pairs)==expected
    assert distances(swapped)[:7,:7]!=D[:7,:7]

    # A chosen hypothetical complete graph; no physical frequencies inferred.
    L14=laplacian(D[:19,:19])
    L15=laplacian(D)
    assert L15[:19,:19]-L14==s.eye(19)/2
    assert all(L15[i,19]==-Q(1,2) for i in range(19))
    assert L15[19,19]==Q(19,2)
    z=s.Symbol("z")
    chi14=s.factor(L14.charpoly(z).as_expr())
    chi15=s.factor(L15.charpoly(z).as_expr())
    recursion=s.factor(z*(z-10)*s.cancel(chi14/z).subs(z,z-Q(1,2)))
    assert s.factor(chi15-recursion)==0
    f6=s.cancel(s.factor(laplacian(D[:7,:7]).charpoly(z).as_expr())/z)
    closed=s.factor(z*(z-10)**9*(z-Q(19,2))**4*f6.subs(z,z-Q(13,2)))
    assert s.factor(chi15-closed)==0
    assert f6.subs(z,Q(7,2))!=0 and f6.subs(z,3)!=0

    report={
      "schema":"POLYBRUTUS-D15-EXACT-CERTIFICATE-v1",
      "status":"PASS","dimension":15,"vertices":20,"pairs":190,
      "affine_rank":15,"centered_gram_rank":15,"radius":1,
      "ambient_storage_coordinates":16,"intrinsic_affine_dimension":15,
      "parent":"Topbrutus/polybrutus-d14 tag v0.1.0",
      "parent_tag_commit":"853687b71dd32b05c807aa47b1c8d7b48d7e97e4",
      "parent_doi":"10.5281/zenodo.24408836",
      "parent_distances_preserved":171,"new_distances":19,
      "new_squared_distances":["2"]*19,
      "histogram_k_4d2":{str(k):int(v) for k,v in sorted(hist.items())},
      "sum_squared_distances":str(dsum),
      "centroid":[str(x) for x in center],
      "gram_psd_certificate":"B=(HX)(HX)^T exact symbolic identity",
      "negative_controls":"4/4 PASS",
      "oscillator_model":"Hypothetical equal-mass weighted graph with w_ij=1/d_ij²; not observed physics",
      "spectral_recurrence":"chi20(z)=z*(z-10)*chi19(z-1/2)/(z-1/2)",
      "spectral_closed_form":"chi20(z)=z*(z-10)^9*(z-19/2)^4*f6(z-13/2), f6(z)=chiP7(z)/z",
      "eigenvalue_10_multiplicity":9,"eigenvalue_19_over_2_multiplicity":4,
      "spectral_charpoly":str(chi15),
      "novel_universal_theorem":"NOT CLAIMED","physical_discovery":"NOT CLAIMED",
      "human_peer_review":"NOT_DONE","sympy":s.__version__,
      "labelled_squared_distances":[{"i":i+1,"j":j+1,"d2":str(D[i,j]),"k":str(4*D[i,j])}
                                   for i,j in pairs]
    }
    out=ROOT/"results"/"certificate_d15.json"
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8",newline="\n")
    print("D15 EXACT PASS: 20 vertices, 190 pairs, affine rank 15; 171/171 parent distances preserved")
    print("MODEL SPECTRUM PASS: eigenvalue 10 x9; 19/2 x4; negative controls 4/4 PASS")
    print("CERTIFICATE_SHA256="+hashlib.sha256(out.read_bytes()).hexdigest())
    return report

if __name__=="__main__":
    run()
