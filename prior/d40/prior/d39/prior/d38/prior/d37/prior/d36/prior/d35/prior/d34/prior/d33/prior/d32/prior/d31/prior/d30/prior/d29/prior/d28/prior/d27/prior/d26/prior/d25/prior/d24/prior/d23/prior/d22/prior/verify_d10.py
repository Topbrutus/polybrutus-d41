#!/usr/bin/env python3
"""Exact Polybrutus D10 certificate; preserves D09 frozen labelled metric.

The canonical 7-point core + regular 5-simplex-vertex block from D07 are
not changed. Each higher dimensional step is an orthogonal unit vector.
D10 adds e11 to D09's 10-coordinate constrained representation.
"""
from pathlib import Path
from collections import Counter
from itertools import combinations
import json,hashlib,sys
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"prior"))
import verify_d09 as parent
Q=s.Rational
points=[p.col_join(s.zeros(1,1)) for p in parent.points]
points.append(s.Matrix([0]*10+[1]))
N=len(points)
edges=list(combinations(range(N),2))

def dist(pts):
    d=s.zeros(len(pts))
    for i,j in combinations(range(len(pts)),2):
        d[i,j]=d[j,i]=s.simplify((pts[i]-pts[j]).dot(pts[i]-pts[j]))
    return d

def coordinate_rank(pts):
    return s.Matrix.hstack(*[p-pts[0] for p in pts[1:]]).rank()

def centered_gram(D):
    h=s.eye(D.rows)-s.ones(D.rows)/D.rows
    return -h*D*h/2

def lap(D):
    L=s.zeros(D.rows)
    for i,j in combinations(range(D.rows),2):
        w=1/D[i,j]
        L[i,i]+=w;L[j,j]+=w
        L[i,j]-=w;L[j,i]-=w
    return L

def run():
    assert len(parent.points)==14 and N==15
    assert all(p.shape==(11,1) for p in points)
    assert all(s.simplify(p.dot(p)-1)==0 for p in points)
    D=dist(points)
    assert all(D[i,j]!=0 for i,j in edges)
    parentdoc=json.loads((ROOT/"prior"/"certificate_d09.json").read_text(encoding="utf8"))
    assert parentdoc["status"]=="PASS" and parentdoc["dimension"]==9
    labelled={(int(e["i"])-1,int(e["j"])-1):s.Rational(e["d2"]) for e in parentdoc["labelled_squared_distances"]}
    assert len(labelled)==91
    assert all(D[i,j]==d2 for (i,j),d2 in labelled.items())
    assert all(D[i,14]==2 for i in range(14))
    hist=Counter(4*D[i,j] for i,j in edges)
    expected=Counter({s.Integer(4):6,s.Integer(7):2,s.Integer(8):74,s.Integer(9):3,
                      s.Integer(10):10,s.Integer(12):6,s.Integer(13):2,s.Integer(16):2})
    assert hist==expected and sum(hist.values())==105
    assert coordinate_rank(points)==10
    B=centered_gram(D)
    assert B.rank()==10
    X=s.Matrix.hstack(*points).T
    H=s.eye(N)-s.ones(N)/N
    assert all(s.simplify(v)==0 for v in B-(H*X)*(H*X).T)
    squared_sum=s.simplify(sum(D[i,j] for i,j in edges))
    assert squared_sum==Q(887,4)
    cent=sum(points,s.zeros(11,1))/N
    assert list(cent)==[0,0,Q(1,30),0,0,0,0,0,Q(1,15),Q(1,15),Q(1,15)]

    # Four corruption controls: distance, duplicate, nonorthogonal, permuted labels.
    corrupted=D.copy();corrupted[0,14]=corrupted[14,0]=Q(9,4)
    assert Counter(4*corrupted[i,j] for i,j in edges)!=expected
    cloned=points[:14]+[points[13]]
    assert coordinate_rank(cloned)==9
    tilted=points[:14]+[s.Matrix([Q(3,5),0,0,0,0,0,0,0,0,0,Q(4,5)])]
    assert tilted[-1].dot(tilted[-1])==1
    assert any(dist(tilted)[i,14]!=2 for i in range(14))
    swap=points.copy();swap[0],swap[7]=swap[7],swap[0]
    assert Counter(4*dist(swap)[i,j] for i,j in edges)==expected
    assert dist(swap)[:7,:7]!=D[:7,:7]

    L9=lap(D[:14,:14]);L10=lap(D)
    assert L10[:14,:14]-L9==s.eye(14)/2
    assert all(L10[i,14]==-Q(1,2) for i in range(14))
    assert L10[14,14]==7
    z=s.Symbol("z")
    chi9=s.factor(L9.charpoly(z).as_expr())
    chi10=s.factor(L10.charpoly(z).as_expr())
    recurrence=s.factor(z*(z-Q(15,2))*s.cancel(chi9/z).subs(z,z-Q(1,2)))
    assert s.factor(chi10-recurrence)==0

    P7=lap(D[:7,:7])
    chiP7=s.factor(P7.charpoly(z).as_expr())
    f6=s.cancel(chiP7/z)
    closed_form=s.factor(z*(z-Q(15,2))**4*(z-7)**4*f6.subs(z,z-4))
    assert s.factor(chi10-closed_form)==0
    assert f6.subs(z,Q(7,2))!=0
    assert f6.subs(z,3)!=0

    report={
      "schema":"POLYBRUTUS-D10-EXACT-CERTIFICATE-v1",
      "status":"PASS","dimension":10,"vertices":15,"pairs":105,
      "affine_rank":10,"centered_gram_rank":10,"radius":1,
      "ambient_storage_coordinates":11,"intrinsic_affine_dimension":10,
      "parent":"Topbrutus/polybrutus-d09 tag v0.1.0",
      "parent_tag_commit":"3e435ea423e829d1b6cbc43fe8913e7338d4cb2c",
      "parent_doi":"10.5281/zenodo.24408448",
      "parent_distances_preserved":91,"new_distances":14,
      "new_squared_distances":["2"]*14,
      "histogram_k_4d2":{str(k):int(v) for k,v in sorted(hist.items())},
      "sum_squared_distances":str(squared_sum),
      "centroid":[str(a) for a in cent],
      "gram_psd_proof":"B=(HX)(HX)^T exact symbolic identity",
      "negative_controls":"4/4 PASS",
      "oscillator_model":"Hypothetical equal-mass complete graph with w_ij=1/d_ij^2; NOT physical observation",
      "spectral_recurrence":"chi15(z)=z*(z-15/2)*chi14(z-1/2)/(z-1/2)",
      "spectral_closed_form":"chi15(z)=z*(z-15/2)^4*(z-7)^4*f6(z-4); f6(z)=chiP7(z)/z",
      "eigenvalue_15_over_2_multiplicity":4,
      "eigenvalue_7_multiplicity":4,
      "spectral_charpoly":str(chi10),
      "novel_universal_theorem":"NOT CLAIMED",
      "physical_discovery":"NOT CLAIMED",
      "human_peer_review":"NOT_DONE",
      "sympy":s.__version__,
      "labelled_squared_distances":[{"i":i+1,"j":j+1,"d2":str(D[i,j]),"k":str(4*D[i,j])} for i,j in edges]
    }
    p=ROOT/"results"/"certificate_d10.json"
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8",newline="\n")
    print("D10 EXACT PASS: 15 points, 105 pairs, affine rank=10, Gram PSD; parent 91/91 distances preserved")
    print("SPECTRAL PASS: lambda 15/2 mult 4, lambda 7 mult 4; counterchecks 4/4 PASS")
    print("CERTIFICATE_SHA256="+hashlib.sha256(p.read_bytes()).hexdigest())
    return report

if __name__=="__main__":run()
