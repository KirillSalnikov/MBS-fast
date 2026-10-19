#!/usr/bin/env python3
"""Compare buffered PO preparation statistics and matrices with a reference binary."""
import argparse,json,subprocess,math,re
from pathlib import Path
from test_fullauto_arbitrary import asymmetric_tetrahedron,rotated_box

def matrices(path):
    text=path.read_text().splitlines()
    rows=[[float(v) for v in line.split()] for line in text[1:] if line.strip()]
    assert rows and all(len(x) in (18,19) for x in rows)
    assert all(math.isfinite(v) for row in rows for v in row)
    return rows

def compare(a,b):
    assert len(a)==len(b)
    worst=0
    for x,y in zip(a,b):
        assert len(x)==len(y) and x[:2]==y[:2]
        first=len(x)-16;scale=max(abs(x[first]),abs(y[first]),1e-100)
        worst=max(worst,max(abs(v-w)/scale for v,w in zip(x[first:],y[first:])))
        if abs(x[first])>1e-100 and abs(y[first])>1e-100:
            worst=max(worst,max(abs(v/x[first]-w/y[first]) for v,w in zip(x[first+1:],y[first+1:])))
    return worst

def main():
    p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--reference-binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--backend',default='cpu');p.add_argument('--case',default='all');p.add_argument('--record-counter-differences',action='store_true',help='Record counter differences for repeatability investigations; matrix and integral gates remain mandatory.');p.add_argument('--tolerance',type=float,default=1e-11);a=p.parse_args()
    root=Path(__file__).resolve().parents[1];out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
    tet=out/'tetrahedron.particle';box=out/'box.particle';asymmetric_tetrahedron(tet);rotated_box(box)
    prism=['--particle','1','10','6'];shape=['--particle-file',str(tet)];absorbing=['--particle-file',str(box)]
    cases=[('prism',prism,[],0),('off',prism,['--cutoff-profile','off'],0),('all_thresholds',prism,['--trace-cutoff-jones','.001','--trace-cutoff-area','.002','--trace-cutoff-importance','.0001'],0),('tetrahedron',shape,[],0),('absorbing',absorbing,[],.02),('concave',['--particle-file',str(root/'examples/particles/concave_hexagonal.particle')],[],0),('paired',prism,['--cutoff-profile','off','--paired-depth','6'],0),('incoherent',shape,['--incoherent'],0)]
    results=[]
    keys=['Branches evaluated','Branches rejected','Candidate output beams','Kept output beams','Rejected output beams','Reason counters','Prepare calls:','Orientations still','Configured beam-limit hits','Hard tree-limit hits','Adaptive retry attempts']
    for name,geometry,extra,imag in cases:
        if a.case!='all' and name not in a.case.split(','):continue
        saved=[]
        for label,binary in [('old',a.reference_binary),('new',a.binary)]:
            extras=list(extra)
            d=out/(name+'_'+label);dest=d/'result';d.mkdir(exist_ok=True)
            cmd=[str(binary.resolve()),'--method','po','--backend',a.backend]+geometry+['--ri','1.3116',str(imag),'--wavelength-um','.532','--symmetry','1','1','--max-reflections','8','--sobol-seed','128','37','--haar-alpha','--scattering-grid','0','180','12','4','--threads','2','--orientation-chunk','16','--cutoff-profile','balanced','--close','--output',str(dest)]
            for k in ['--cutoff-profile']:
                if k in extras:cmd[cmd.index(k)+1]=extras[extras.index(k)+1];extras=extras[:extras.index(k)]+extras[extras.index(k)+2:]
            if any(k.startswith('--trace-cutoff-') for k in extras):
                i=cmd.index('--cutoff-profile');del cmd[i:i+2]
            cmd+=extras;(d/'command.args').write_text('\n'.join(cmd)+'\n')
            with (d/'run.log').open('w') as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=180)
            # Final DAT is angle-averaged; paired output can retain its raw phi grid.
            mats={x.name:matrices(x) for x in dest.glob('*.dat')}
            log=(d/'run.log').read_text();stats=[x for x in log.splitlines() if any(x.startswith(k) for k in keys)]
            assert stats
            integral=(dest/'result_integrals.tsv').read_text()
            saved.append((mats,stats,integral))
        old,new=saved;assert old[0].keys()==new[0].keys()
        if not a.record_counter_differences:assert old[1]==new[1],(name,old[1],new[1])
        error=max(compare(old[0][key],new[0][key]) for key in old[0]);assert error<a.tolerance,(name,error)
        # The integrals file has numeric and status fields; check each numeric
        # field within roundoff, and require exact strings for status labels.
        x=old[2].split();y=new[2].split();assert len(x)==len(y)
        for v,w in zip(x,y):
            try:vf,wf=float(v),float(w)
            except ValueError:assert v==w;continue
            if math.isnan(vf) and math.isnan(wf):continue
            assert math.isfinite(vf) and math.isfinite(wf)
            assert abs(vf-wf)<=a.tolerance*max(1,abs(vf),abs(wf)),(name,v,w)
        results.append(dict(case=name,all16_and_ratios_error=error,counters_identical=old[1]==new[1],counter_differences=[dict(old=v,new=w) for v,w in zip(old[1],new[1]) if v!=w],integrals_match=True))
        (out/'comparison.json').write_text(json.dumps(results,indent=2));print('PASS',name,error,flush=True)

if __name__=='__main__':main()
