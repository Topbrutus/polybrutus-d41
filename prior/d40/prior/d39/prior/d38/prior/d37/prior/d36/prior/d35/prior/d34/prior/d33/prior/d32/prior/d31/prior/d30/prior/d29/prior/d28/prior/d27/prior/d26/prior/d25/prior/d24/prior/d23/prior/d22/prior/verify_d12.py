#!/usr/bin/env python3
"""Polybrutus D12: exact conservative dimension 12, labelled EDM and spectrum.

Source: frozen Polybrutus D11 v0.1.0, commit
b49252a03760ea55301b15a60df3c7317f8e3160.
All 16 original vertices are embedded unaltered with a trailing zero.
One orthogonal unit vector is appended to the 13 storage coordinates
(only 12 intrinsic affine dimensions because the centered simplex block
of five coordinates has one linear constraint).
"""
from pathlib import Path
from collections import Counter
from itertools import combinations
import hashlib,json,sys
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"prior"))
import verify_d11 as old
Q=s.Rational
points=[p.col_join(s.zeros(1,1)) for p in old.points]
points.append(s.Matrix([0]*12+[1]))
N=len(points)
pairs=list(combinations(range(N),2))

def dists(pts):
    M=s.zeros(len(pts))
    for i,j in combinations(range(len(pts)),2):
        M[i,j]=M[j,i]=s.simplify((pts[i]-pts[j]).dot(pts[i]-pts[j]))
    return M

def rank_of_differences(pts):
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
    assert len(old.points)==16 and N==17
    assert all(p.shape==(13,1) for p in points)
    assert all(s.simplify(p.dot(p)-1)==0 for p in points)
    D=dists(points)
    assert all(D[i,j]>0 for i,j in pairs)
    previous=json.loads((ROOT/"prior"/"certificate_d11.json").read_text(encoding="utf8"))
    assert previous["status"]=="PASS" and previous["dimension"]==11
    former={(int(e["i"])-1,int(e["j"])-1):s.Rational(e["d2"])
            for e in previous["labelled_squared_distances"]}
    assert len(former)==120
    assert all(D[i,j]==d2 for (i,j),d2 in former.items())
    assert all(D[i,16]==2 for i in range(16))
    histogram=Counter(4*D[i,j] for i,j in pairs)
    target=Counter({s.Integer(4):6,s.Integer(7):2,s.Integer(8):105,s.Integer(9):3,
                    s.Integer(10):10,s.Integer(12):6,s.Integer(13):2,s.Integer(16):2})
    assert histogram==target and sum(histogram.values())==136
    assert rank_of_differences(points)==12
    B=gram(D)
    assert B.rank()==12
    X=s.Matrix.hstack(*points).T
    H=s.eye(N)-s.ones(N)/N
    assert all(s.simplify(q)==0 for q in B-(H*X)*(H*X).T)
    total=s.simplify(sum(D[i,j] for i,j in pairs))
    assert total==Q(1135,4)
    center=sum(points,s.zeros(13,1))/N
    assert list(center)==[0,0,Q(1,34),0,0,0,0,0,Q(1,17),Q(1,17),Q(1,17),Q(1,17),Q(1,17)]

    # Deliberately incorrect candidates: must be caught by validation.
    corrupt=D.copy();corrupt[0,16]=corrupt[16,0]=Q(9,4)
    assert Counter(4*corrupt[i,j] for i,j in pairs)!=target
    repeat=points[:16]+[points[15]]
    assert rank_of_differences(repeat)==11
    tilt=points[:16]+[s.Matrix([Q(3,5),0,0,0,0,0,0,0,0,0,0,0,Q(4,5)])]
    assert tilt[-1].dot(tilt[-1])==1
    assert any(dists(tilt)[i,16]!=2 for i in range(16))
    swap=points.copy();swap[0],swap[7]=swap[7],swap[0]
    assert Counter(4*dists(swap)[i,j] for i,j in pairs)==target
    assert dists(swap)[:7,:7]!=D[:7,:7]

    # Ideal weighted complete graph, NOT laboratory physics.
    L11=laplacian(D[:16,:16])
    L12=laplacian(D)
    assert L12[:16,:16]-L11==s.eye(16)/2
    assert all(L12[i,16]==-Q(1,2) for i in range(16))
    assert L12[16,16]==8
    z=s.Symbol("z")
    chi11=s.factor(L11.charpoly(z).as_expr())
    chi12=s.factor(L12.charpoly(z).as_expr())
    recursion=s.factor(z*(z-Q(17,2))*s.cancel(chi11/z).subs(z,z-Q(1,2)))
    assert s.factor(chi12-recursion)==0
    f6=s.cancel(s.factor(laplacian(D[:7,:7]).charpoly(z).as_expr())/z)
    closed=s.factor(z*(z-Q(17,2))**6*(z-8)**4*f6.subs(z,z-5))
    assert s.factor(chi12-closed)==0
    assert f6.subs(z,Q(7,2))!=0 and f6.subs(z,3)!=0

    report={
      "schema":"POLYBRUTUS-D12-EXACT-CERTIFICATE-v1",
      "status":"PASS","dimension":12,"vertices":17,"pairs":136,
      "affine_rank":12,"centered_gram_rank":12,"radius":1,
      "ambient_storage_coordinates":13,"intrinsic_affine_dimension":12,
      "parent":"Topbrutus/polybrutus-d11 tag v0.1.0",
      "parent_tag_commit":"b49252a03760ea55301b15a60df3c7317f8e3160",
      "parent_doi":"10.5281/zenodo.24408650",
      "parent_distances_preserved":120,"new_distances":16,
      "new_squared_distances":["2"]*16,
      "histogram_k_4d2":{str(k):int(v) for k,v in sorted(histogram.items())},
      "sum_squared_distances":str(total),
      "centroid":[str(v) for v in center],
      "gram_psd_proof":"B=(HX)(HX)^T exact symbolic identity",
      "negative_controls":"4/4 PASS",
      "oscillator_model":"Hypothetical complete graph, unit masses, w_ij=1/d_ij²; not physical evidence",
      "spectral_recurrence":"chi17(z)=z*(z-17/2)*chi16(z-1/2)/(z-1/2)",
      "spectral_closed_form":"chi17(z)=z*(z-17/2)^6*(z-8)^4*f6(z-5), f6(z)=chiP7(z)/z",
      "eigenvalue_17_over_2_multiplicity":6,
      "eigenvalue_8_multiplicity":4,
      "spectral_charpoly":str(chi12),
      "relation_to_twelve_pythagoras":"NO CLAIM: dimension number 12 is not a new proof of the separate 12-35-37 relation",
      "novel_universal_theorem":"NOT CLAIMED","physical_discovery":"NOT CLAIMED",
      "human_peer_review":"NOT_DONE","sympy":s.__version__,
      "labelled_squared_distances":[{"i":i+1,"j":j+1,"d2":str(D[i,j]),"k":str(4*D[i,j])} for i,j in pairs]
    }
    path=ROOT/"results"/"certificate_d12.json"
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8",newline="\n")
    print("D12 EXACT PASS: 17 points, 136 pairs, rank 12; parent 120/120 distances preserved")
    print("GRAPH SPECTRUM PASS: eigenvalue 17/2 x6; eigenvalue 8 x4; four counterchecks PASS")
    print("CERTIFICATE_SHA256="+hashlib.sha256(path.read_bytes()).hexdigest())
    return report
if __name__=="__main__":
    run()
