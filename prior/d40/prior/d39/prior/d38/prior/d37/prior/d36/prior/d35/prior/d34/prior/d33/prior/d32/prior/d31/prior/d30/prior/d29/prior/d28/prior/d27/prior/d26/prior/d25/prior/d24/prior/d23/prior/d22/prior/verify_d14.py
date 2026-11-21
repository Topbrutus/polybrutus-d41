#!/usr/bin/env python3
"""Polybrutus D14 - exact conservative Euclidean metric certificate.

Frozen predecessor: Topbrutus/polybrutus-d13@v0.1.0,
commit 4d5f10de586e6677ddbe1d90b18218b46602d14e.
No prior labelled point is moved; an orthogonal radius-one point is added.
Calculations are exact SymPy, not physical measurements or novelty claims.
"""
from pathlib import Path
from collections import Counter
from itertools import combinations
import hashlib,json,sys
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"prior"))
import verify_d13 as parent

Q=s.Rational
points=[p.col_join(s.zeros(1,1)) for p in parent.points]
points.append(s.Matrix([0]*14+[1]))
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
    assert len(parent.points)==18 and N==19
    assert all(p.shape==(15,1) for p in points)
    assert all(s.simplify(p.dot(p)-1)==0 for p in points)
    D=distances(points)
    assert all(D[i,j]>0 for i,j in pairs)
    old=json.loads((ROOT/"prior"/"certificate_d13.json").read_text(encoding="utf8"))
    assert old["status"]=="PASS" and old["dimension"]==13
    old_pairs={(int(e["i"])-1,int(e["j"])-1):s.Rational(e["d2"])
               for e in old["labelled_squared_distances"]}
    assert len(old_pairs)==153
    assert all(D[i,j]==v for (i,j),v in old_pairs.items())
    assert all(D[i,18]==2 for i in range(18))
    hist=Counter(4*D[i,j] for i,j in pairs)
    expected=Counter({s.Integer(4):6,s.Integer(7):2,s.Integer(8):140,
                       s.Integer(9):3,s.Integer(10):10,s.Integer(12):6,
                       s.Integer(13):2,s.Integer(16):2})
    assert hist==expected and sum(hist.values())==171
    assert affine_rank(points)==14
    B=centered_gram(D)
    assert B.rank()==14
    X=s.Matrix.hstack(*points).T
    H=s.eye(N)-s.ones(N)/N
    assert all(s.simplify(z)==0 for z in B-(H*X)*(H*X).T)
    dsum=s.simplify(sum(D[i,j] for i,j in pairs))
    assert dsum==Q(1415,4)
    center=sum(points,s.zeros(15,1))/N
    assert list(center)==[0,0,Q(1,38),0,0,0,0,0]+[Q(1,19)]*7

    # Deliberately incorrect variants must fail independent criteria.
    corrupt=D.copy();corrupt[0,18]=corrupt[18,0]=Q(9,4)
    assert Counter(4*corrupt[i,j] for i,j in pairs)!=expected
    repeated=points[:18]+[points[17]]
    assert affine_rank(repeated)==13
    tilted=points[:18]+[s.Matrix([Q(3,5)]+[0]*13+[Q(4,5)])]
    assert tilted[-1].dot(tilted[-1])==1
    assert any(distances(tilted)[i,18]!=2 for i in range(18))
    swapped=points.copy();swapped[0],swapped[7]=swapped[7],swapped[0]
    assert Counter(4*distances(swapped)[i,j] for i,j in pairs)==expected
    assert distances(swapped)[:7,:7]!=D[:7,:7]

    # Ideal chosen graph model, not experimental frequencies.
    L13=laplacian(D[:18,:18])
    L14=laplacian(D)
    assert L14[:18,:18]-L13==s.eye(18)/2
    assert all(L14[i,18]==-Q(1,2) for i in range(18))
    assert L14[18,18]==9
    z=s.Symbol("z")
    chi13=s.factor(L13.charpoly(z).as_expr())
    chi14=s.factor(L14.charpoly(z).as_expr())
    recursive=s.factor(z*(z-Q(19,2))*s.cancel(chi13/z).subs(z,z-Q(1,2)))
    assert s.factor(chi14-recursive)==0
    f6=s.cancel(s.factor(laplacian(D[:7,:7]).charpoly(z).as_expr())/z)
    closed=s.factor(z*(z-Q(19,2))**8*(z-9)**4*f6.subs(z,z-6))
    assert s.factor(chi14-closed)==0
    assert f6.subs(z,Q(7,2))!=0 and f6.subs(z,3)!=0

    report={
      "schema":"POLYBRUTUS-D14-EXACT-CERTIFICATE-v1",
      "status":"PASS","dimension":14,"vertices":19,"pairs":171,
      "affine_rank":14,"centered_gram_rank":14,"radius":1,
      "ambient_storage_coordinates":15,"intrinsic_affine_dimension":14,
      "parent":"Topbrutus/polybrutus-d13 tag v0.1.0",
      "parent_tag_commit":"4d5f10de586e6677ddbe1d90b18218b46602d14e",
      "parent_doi":"10.5281/zenodo.24408769",
      "parent_distances_preserved":153,"new_distances":18,
      "new_squared_distances":["2"]*18,
      "histogram_k_4d2":{str(k):int(v) for k,v in sorted(hist.items())},
      "sum_squared_distances":str(dsum),
      "centroid":[str(v) for v in center],
      "gram_psd_certificate":"B=(HX)(HX)^T exact symbolic identity",
      "negative_controls":"4/4 PASS",
      "oscillator_model":"Ideal equal-mass complete weighted graph with w_ij=1/d_ij²; no physical observation",
      "spectral_recurrence":"chi19(z)=z*(z-19/2)*chi18(z-1/2)/(z-1/2)",
      "spectral_closed_form":"chi19(z)=z*(z-19/2)^8*(z-9)^4*f6(z-6), f6(z)=chiP7(z)/z",
      "eigenvalue_19_over_2_multiplicity":8,
      "eigenvalue_9_multiplicity":4,
      "spectral_charpoly":str(chi14),
      "novel_universal_theorem":"NOT CLAIMED","physical_discovery":"NOT CLAIMED",
      "human_peer_review":"NOT_DONE","sympy":s.__version__,
      "labelled_squared_distances":[{"i":i+1,"j":j+1,"d2":str(D[i,j]),"k":str(4*D[i,j])}
                                   for i,j in pairs]
    }
    out=ROOT/"results"/"certificate_d14.json"
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8",newline="\n")
    print("D14 EXACT PASS: 19 points, 171 pairs, affine rank 14, Gram PSD, old 153/153 preserved")
    print("IDEAL GRAPH PASS: eigenvalue 19/2 multiplicity 8, eigenvalue 9 multiplicity 4; negative controls 4/4")
    print("CERTIFICATE_SHA256="+hashlib.sha256(out.read_bytes()).hexdigest())
    return report

if __name__=="__main__":
    run()
