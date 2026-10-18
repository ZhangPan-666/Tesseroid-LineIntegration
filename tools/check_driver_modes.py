from pathlib import Path
import struct, subprocess, os, math, json

import argparse
parser=argparse.ArgumentParser(description='Compare serial, OpenMP, MPI and component-selective batch drivers.')
parser.add_argument('--serial-dir',type=Path,required=True,help='Directory of built serial/OpenMP programs')
parser.add_argument('--mpi-dir',type=Path,required=True,help='Directory of built MPI programs')
parser.add_argument('--mpiexec',default='mpiexec')
parser.add_argument('--scratch',type=Path,required=True,help='Dedicated output directory; fixed filenames are overwritten')
args=parser.parse_args()
serial=args.serial_dir.resolve()
mpi=args.mpi_dir.resolve()
mpiexec=args.mpiexec
base=args.scratch.resolve()
base.mkdir(parents=True,exist_ok=True)
env=os.environ.copy()
env['OMP_NUM_THREADS']='2'
suffix='.exe' if os.name=='nt' else ''

def doubles(values): return struct.pack('<'+'d'*len(values),*values)
def ordinary(magnetic):
    fields=[34.8,35.2,119.7,120.3,6341.2,6371.2]+([.5,.2,-.1] if magnetic else [2670])
    return doubles([1e-12,1e-9])+struct.pack('<i',1)+doubles(fields)+struct.pack('<i',2)+doubles([120.08,120.21,35.11,35.17,6771.2,6471.2])
def shell(magnetic):
    fields=[6371.2,6341.2]+([.5,.2,-.1] if magnetic else [1])
    data=doubles([1e-12,1e-9,-175,175,-85,85,10,10])+struct.pack('<ii',36,18)
    for v in fields: data+=doubles([v]*648)
    return data+struct.pack('<i',1)+doubles([138.91,24.37,6771.2])
def run(name,folder,parameter,result,content,use_mpi):
    directory=base/name
    directory.mkdir(exist_ok=True)
    (directory/parameter).write_bytes(content)
    command=([str(mpiexec),'-n','2'] if use_mpi else [])+[str(folder/(name+suffix))]
    p=subprocess.run(command,cwd=directory,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=90)
    (directory/'run.log').write_bytes(p.stdout)
    if p.returncode: raise RuntimeError(f'{name} returned {p.returncode}: {p.stdout.decode(errors="replace")}')
    raw=(directory/result).read_bytes()
    n=struct.unpack_from('<i',raw)[0]
    values=struct.unpack_from('<'+'d'*((len(raw)-4)//8),raw,4)
    return n,values
def compare(a,b,indices=None):
    assert a[0]==b[0]
    n=a[0]
    if indices is None: indices=range(len(a[1])//n)
    worst=0
    for c in indices:
        for j in range(n):
            x=a[1][c*n+j];y=b[1][c*n+j]
            assert math.isfinite(x) and math.isfinite(y)
            err=abs(x-y)/max(abs(x),abs(y),1e-25)
            assert err<1e-9,(c,j,x,y,err)
            worst=max(worst,err)
    return worst
results={}
for field,magnetic in [('Gravity',False),('Magnetic',True)]:
    parameter='TFM.ForPar' if magnetic else 'TFG.ForPar'
    result='Result_TFM.dat' if magnetic else 'Result_TFG.dat'
    ref=run('Work1_'+field+'Estimate',serial,parameter,result,ordinary(magnetic),False)
    omp=run('Work1_'+field+'EstimateOpenMP',serial,parameter,result,ordinary(magnetic),False)
    parallel=run('Work1_'+field+'EstimateMPI',mpi,parameter,result,ordinary(magnetic),True)
    results[field]={'serial_vs_openmp':compare(ref,omp),'serial_vs_mpi_2_ranks':compare(ref,parallel)}
    print(field,results[field],flush=True)
for field,magnetic in [('Mass',False),('Magnetization',True)]:
    parameter='MSIM.ForPar' if magnetic else 'MSIG.ForPar'
    result='ResultM.dat' if magnetic else 'ResultG.dat'
    full=run('Work1_'+field+'ShellIntegrationMPI',mpi,parameter,result,shell(magnetic),True)
    selected_name='Work2_MagnetizationShellIntegrationVectorMPI' if magnetic else 'Work2_MassShellIntegrationVzMPI'
    selected=run(selected_name,mpi,parameter,result,shell(magnetic),True)
    indices=[0,1,2,4,5,6] if magnetic else [0,1,2,6]
    results[field]={'full_vs_selected':compare(full,selected,indices)}
    print(field,results[field],flush=True)
(base/'results.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
