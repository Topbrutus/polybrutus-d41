#!/usr/bin/env python3
"""Polybrutus 9D exact labelled-metric and oscillator certificate.

Source of prior 8D coordinates is the FROZEN PUBLIC v0.1.0 tag copied into
prior/verify_d08.py. All 13 old coordinates gain a trailing zero, then a new
orthogonal unit vector is appended. No external IO or secret is required.
"""
from pathlib import Path
from collections import Counter
from itertools import combinations
import hashlib, json, sys
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"prior"))
import verify_d08 as old

Q=s.Rational
points=[p.col_join(s.zeros(1,1)) for p in old.points]
points.append(s.Matrix([0]*9+[1]))
N=len(points)
pairs=list(combinations(range(N),2))

def distances(pts):
    out=s.zeros(len(pts))
    for i,j in combinations(range(len(pts)),2):
        out[i,j]=out[j,i]=s.simplify((pts[i]-pts[j]).dot(pts[i]-pts[j]))
    return out

def affine_rank(pts):
    return s.Matrix.hstack(*[p-pts[0] for p in pts[1:]]).rank()

def gram(D):
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
    assert len(old.points)==13 and len(points)==14
    assert all(p.rows==10 for p in points)
    assert all(s.simplify(p.dot(p)-1)==0 for p in points)
    D=distances(points)
    D8=D[:13,:13]
    previous=json.loads((ROOT/"prior"/"certificate_d08.json").read_text(encoding="utf8"))
    assert previous["status"]=="PASS" and previous["dimension"]==8
    priorpairs={(int(e["i"])-1,int(e["j"])-1):s.Rational(e["d2"]) for e in previous["pairwise_distances"]}
    assert len(priorpairs)==78
    assert all(D8[i,j]==d2 for (i,j),d2 in priorpairs.items())
    assert all(D[13,i]==2 for i in range(13))
    assert all(D[i,j]>0 for i,j in pairs)
    hist=Counter(4*D[i,j] for i,j in pairs)
    expect=Counter({s.Integer(4):6,s.Integer(7):2,s.Integer(8):60,
                    s.Integer(9):3,s.Integer(10):10,s.Integer(12):6,
                    s.Integer(13):2,s.Integer(16):2})
    assert hist==expect and sum(hist.values())==91
    assert affine_rank(points)==9
    B=gram(D)
    assert B.rank()==9
    X=s.Matrix.hstack(*points).T
    H=s.eye(N)-s.ones(N)/N
    assert all(s.simplify(v)==0 for v in B-(H*X)*(H*X).T)
    S=s.simplify(sum(D[i,j] for i,j in pairs))
    assert S==s.Rational(775,4)
    cent=sum(points,s.zeros(10,1))/N
    assert list(cent)==[0,0,Q(1,28),0,0,0,0,0,Q(1,14),Q(1,14)]
    # Four negative controls, deliberately bad cases that must fail.
    tamper=D.copy();tamper[0,13]=tamper[13,0]=Q(9,4)
    assert Counter(4*tamper[i,j] for i,j in pairs)!=expect
    duplicate=points[:13]+[points[12]]
    assert affine_rank(duplicate)==8
    bent=points[:13]+[s.Matrix([Q(3,5),0,0,0,0,0,0,0,0,Q(4,5)])]
    assert bent[-1].dot(bent[-1])==1
    assert any(distances(bent)[13,i]!=2 for i in range(13))
    sw=points.copy();sw[0],sw[7]=sw[7],sw[0]
    assert Counter(4*distances(sw)[i,j] for i,j in pairs)==expect
    assert distances(sw)[:7,:7]!=D[:7,:7]

    L8=laplacian(D8)
    L9=laplacian(D)
    assert L9[:13,:13]-L8==s.eye(13)/2
    assert all(L9[i,13]==-Q(1,2) for i in range(13))
    assert L9[13,13]==Q(13,2)
    z=s.symbols("z")
    chi8=s.factor(L8.charpoly(z).as_expr())
    chi9=s.factor(L9.charpoly(z).as_expr())
    predicted=s.factor(z*(z-7)*s.cancel(chi8/z).subs(z,z-Q(1,2)))
    assert s.factor(chi9-predicted)==0
    # Independent closed form from the canonical 7-point internal Laplacian.
    chiP=s.factor(laplacian(D[:7,:7]).charpoly(z).as_expr())
    f6=s.cancel(chiP/z)
    formula=s.factor(z*(z-7)**3*(z-Q(13,2))**4*f6.subs(z,z-Q(7,2)))
    assert s.factor(chi9-formula)==0
    # Exactly those multiplicities, not merely at least as large.
    assert f6.subs(z,Q(7,2))!=0
    assert f6.subs(z,3)!=0
    assert chi9.subs(z,7)==0 and chi9.subs(z,Q(13,2))==0
    report={
      "schema":"POLYBRUTUS-D09-EXACT-CERTIFICATE-v1",
      "status":"PASS","dimension":9,"vertices":14,"pairs":91,
      "affine_rank":9,"centered_gram_rank":9,"radius":1,
      "ambient_storage_coordinates":10,"unconstrained_affine_dimension":9,
      "parent":"Topbrutus/polybrutus-d08 tag v0.1.0",
      "parent_tag_commit":"356973b71ec14aa5881a8e26edd375efdb84f0af",
      "parent_doi":"10.5281/zenodo.24408344",
      "parent_distances_preserved":78,
      "new_distances":13,
      "new_squared_distances":["2"]*13,
      "histogram_k_4d2":{str(k):int(v) for k,v in sorted(hist.items())},
      "sum_squared_distances":str(S),
      "centroid":[str(x) for x in cent],
      "gram_psd_certificate":"B=(H X)(H X)^T exactly (symbolic)",
      "negative_controls":"4/4 PASS",
      "spectral_model":"Hypothetical equal-mass scalar oscillator network; all-pair coupling 1/d_ij²",
      "spectral_recursion":"chi14(z)=z*(z-7)*chi13(z-1/2)/(z-1/2)",
      "spectral_closed_form":"chi14(z)=z*(z-7)^3*(z-13/2)^4*f6(z-7/2), f6(z)=chiP7(z)/z",
      "lambda_7_multiplicity":3,
      "lambda_13_over_2_multiplicity":4,
      "spectral_charpoly_9d":str(chi9),
      "claims":{"new_universal_theorem":False,"new_physical_dimensions":False,"human_peer_review":False},
      "sympy":s.__version__,
      "labelled_squared_distances":[{"i":i+1,"j":j+1,"d2":str(D[i,j]),"k":str(4*D[i,j])} for i,j in pairs]
    }
    out=ROOT/"results"/"certificate_d09.json"
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf8",newline="\n")
    print("D09 EXACT PASS: N=14, pairs=91, rank=9, Gram PSD; old 78/78 preserved")
    print("SPECTRAL PASS: lambda=7 x3; lambda=13/2 x4; four negative controls PASS")
    print("CERTIFICATE_SHA256="+hashlib.sha256(out.read_bytes()).hexdigest())
    return report

if __name__=="__main__":run()
