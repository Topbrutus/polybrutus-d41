#!/usr/bin/env python3
"""Exact Polybrutus D11 certificate; preserves labelled D10 metric.

Polybrutus D11 adds e12 to the immutable 15-point D10 representation.
A single simplex-block relation means the 12 storage slots encode an
intrinsic 11-dimensional Euclidean affine configuration.
"""
from pathlib import Path
from collections import Counter
from itertools import combinations
import json,hashlib,sys
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"prior"))
import verify_d10 as parent
Q=s.Rational
points=[p.col_join(s.zeros(1,1)) for p in parent.points]
points.append(s.Matrix([0]*11+[1]))
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
    H=s.eye(D.rows)-s.ones(D.rows)/D.rows
    return -H*D*H/2

def lap(D):
    L=s.zeros(D.rows)
    for i,j in combinations(range(D.rows),2):
        w=1/D[i,j]
        L[i,i]+=w;L[j,j]+=w;L[i,j]-=w;L[j,i]-=w
    return L

def run():
    assert len(parent.points)==15 and N==16
    assert all(p.shape==(12,1) for p in points)
    assert all(s.simplify(p.dot(p)-1)==0 for p in points)
    D=dist(points)
    assert all(D[i,j]>0 for i,j in edges)
    previous=json.loads((ROOT/"prior"/"certificate_d10.json").read_text(encoding="utf8"))
    assert previous["status"]=="PASS" and previous["dimension"]==10
    labelled={(int(e["i"])-1,int(e["j"])-1):s.Rational(e["d2"]) for e in previous["labelled_squared_distances"]}
    assert len(labelled)==105
    assert all(D[i,j]==d2 for (i,j),d2 in labelled.items())
    assert all(D[i,15]==2 for i in range(15))
    hist=Counter(4*D[i,j] for i,j in edges)
    expected=Counter({s.Integer(4):6,s.Integer(7):2,s.Integer(8):89,s.Integer(9):3,
                      s.Integer(10):10,s.Integer(12):6,s.Integer(13):2,s.Integer(16):2})
    assert hist==expected and sum(hist.values())==120
    assert coordinate_rank(points)==11
    B=centered_gram(D)
    assert B.rank()==11
    X=s.Matrix.hstack(*points).T
    H=s.eye(N)-s.ones(N)/N
    assert all(s.simplify(v)==0 for v in B-(H*X)*(H*X).T)
    squared_sum=s.simplify(sum(D[i,j] for i,j in edges))
    assert squared_sum==Q(1007,4)
    cent=sum(points,s.zeros(12,1))/N
    assert list(cent)==[0,0,Q(1,32),0,0,0,0,0,Q(1,16),Q(1,16),Q(1,16),Q(1,16)]

    # Four counterchecks: intentional corruption must fail the right criterion.
    corrupted=D.copy();corrupted[0,15]=corrupted[15,0]=Q(9,4)
    assert Counter(4*corrupted[i,j] for i,j in edges)!=expected
    cloned=points[:15]+[points[14]]
    assert coordinate_rank(cloned)==10
    tilted=points[:15]+[s.Matrix([Q(3,5),0,0,0,0,0,0,0,0,0,0,Q(4,5)])]
    assert tilted[-1].dot(tilted[-1])==1
    assert any(dist(tilted)[i,15]!=2 for i in range(15))
    swap=points.copy();swap[0],swap[7]=swap[7],swap[0]
    assert Counter(4*dist(swap)[i,j] for i,j in edges)==expected
    assert dist(swap)[:7,:7]!=D[:7,:7]

    L10=lap(D[:15,:15]);L11=lap(D)
    assert L11[:15,:15]-L10==s.eye(15)/2
    assert all(L11[i,15]==-Q(1,2) for i in range(15))
    assert L11[15,15]==Q(15,2)
    z=s.Symbol("z")
    chi10=s.factor(L10.charpoly(z).as_expr())
    chi11=s.factor(L11.charpoly(z).as_expr())
    recurrence=s.factor(z*(z-8)*s.cancel(chi10/z).subs(z,z-Q(1,2)))
    assert s.factor(chi11-recurrence)==0

    chiP7=s.factor(lap(D[:7,:7]).charpoly(z).as_expr())
    f6=s.cancel(chiP7/z)
    closed_form=s.factor(z*(z-8)**5*(z-Q(15,2))**4*f6.subs(z,z-Q(9,2)))
    assert s.factor(chi11-closed_form)==0
    assert f6.subs(z,Q(7,2))!=0
    assert f6.subs(z,3)!=0
    assert chi11.subs(z,8)==0 and chi11.subs(z,Q(15,2))==0

    report={
      "schema":"POLYBRUTUS-D11-EXACT-CERTIFICATE-v1",
      "status":"PASS","dimension":11,"vertices":16,"pairs":120,
      "affine_rank":11,"centered_gram_rank":11,"radius":1,
      "ambient_storage_coordinates":12,"intrinsic_affine_dimension":11,
      "parent":"Topbrutus/polybrutus-d10 tag v0.1.0",
      "parent_tag_commit":"dc2a45522e5e6794742b8a43de6c105a475d1e89",
      "parent_doi":"10.5281/zenodo.24408569",
      "parent_distances_preserved":105,"new_distances":15,
      "new_squared_distances":["2"]*15,
      "histogram_k_4d2":{str(k):int(v) for k,v in sorted(hist.items())},
      "sum_squared_distances":str(squared_sum),
      "centroid":[str(a) for a in cent],
      "gram_psd_proof":"B=(HX)(HX)^T exact symbolic identity",
      "negative_controls":"4/4 PASS",
      "oscillator_model":"Hypothetical equal-mass complete graph with weights w_ij=1/d_ij^2; NOT a physical observation",
      "spectral_recurrence":"chi16(z)=z*(z-8)*chi15(z-1/2)/(z-1/2)",
      "spectral_closed_form":"chi16(z)=z*(z-8)^5*(z-15/2)^4*f6(z-9/2); f6(z)=chiP7(z)/z",
      "eigenvalue_8_multiplicity":5,"eigenvalue_15_over_2_multiplicity":4,
      "spectral_charpoly":str(chi11),
      "novel_universal_theorem":"NOT CLAIMED","physical_discovery":"NOT CLAIMED",
      "human_peer_review":"NOT_DONE","sympy":s.__version__,
      "labelled_squared_distances":[{"i":i+1,"j":j+1,"d2":str(D[i,j]),"k":str(4*D[i,j])} for i,j in edges]
    }
    output=ROOT/"results"/"certificate_d11.json"
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8",newline="\n")
    print("D11 EXACT PASS: 16 points, 120 pairs, affine rank=11, Gram PSD; parent 105/105 distances preserved")
    print("SPECTRAL PASS: lambda 8 multiplicity 5, lambda 15/2 multiplicity 4; negative controls 4/4 PASS")
    print("CERTIFICATE_SHA256="+hashlib.sha256(output.read_bytes()).hexdigest())
    return report

if __name__=="__main__":run()
