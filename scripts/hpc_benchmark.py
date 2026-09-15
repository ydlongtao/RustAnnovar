#!/usr/bin/env python3
"""Serial HPC benchmark; records load and RSS, never clears global caches."""
import argparse, fcntl, hashlib, itertools, json, os, pathlib, statistics, subprocess, time


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1048576), b''): h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--workspace', type=pathlib.Path, default=pathlib.Path('/DATABANK/users/hflt/RustAnnovar'))
    parser.add_argument('--mode', choices=['smoke', 'scale', 'giab'], default='smoke')
    args = parser.parse_args(); root = args.workspace
    lock = open(root/'runs/benchmark.lock', 'w'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    run = root/'runs'/(time.strftime('%Y%m%dT%H%M%S')+'-'+args.mode); run.mkdir(exist_ok=False)
    binaries = {n: root/'builds'/n/'rust-annovar' for n in ['baseline','candidate']}
    topology = subprocess.check_output(['lscpu','-p=CPU,CORE,SOCKET,NODE'], text=True)
    cpus, physical, node = [], set(), None
    for line in topology.splitlines():
        if line.startswith('#'): continue
        cpu, core, socket, n = map(int, line.split(','))
        if cpu not in os.sched_getaffinity(0): continue
        if node is None: node = n
        if n == node and (socket,core) not in physical: cpus.append(cpu); physical.add((socket,core))
    if len(cpus)<32: raise RuntimeError('Need 32 allowed physical cores on one NUMA node')
    envinfo = {'topology':topology,'uname':subprocess.check_output(['uname','-a'],text=True),
               'storage':subprocess.check_output(['df','-T',str(root)],text=True),
               'binaries':{n:{'sha256':digest(b),'path':str(b)} for n,b in binaries.items()},'cpus':cpus[:32],
               'cache':'first-run followed by warm repeats; not cold-cache'}
    (run/'environment.json').write_text(json.dumps(envinfo,indent=2)); rows=[]
    def measure(name,cmd,threads,label,repeat):
        prefix=run/f'{label}-{name}-t{threads}-r{repeat}'
        before={'load':os.getloadavg(),'meminfo':pathlib.Path('/proc/meminfo').read_text(),'vmstat':pathlib.Path('/proc/vmstat').read_text()}
        with open(str(prefix)+'.stdout','w') as out,open(str(prefix)+'.stderr','w') as err:
            code=subprocess.call(['taskset','-c',','.join(map(str,cpus[:threads])),'/usr/bin/time','-f','%e\t%M','-o',str(prefix)+'.time',*map(str,cmd)],env=dict(os.environ,RAYON_NUM_THREADS=str(threads)),stdout=out,stderr=err)
        row={'implementation':name,'threads':threads,'label':label,'repeat':repeat,'exit_code':code,'before':before,'after_load':os.getloadavg(),'after_vmstat':pathlib.Path('/proc/vmstat').read_text(),'command':list(map(str,cmd))}
        if code==0:
            wall,rss=pathlib.Path(str(prefix)+'.time').read_text().strip().split('\t');row.update(wall_seconds=float(wall),peak_rss_kib=int(rss))
        with open(run/'measurements.jsonl','a') as f:f.write(json.dumps(row)+'\n')
        if code:raise RuntimeError(f'Failed benchmark: {prefix}')
        rows.append(row)
    if args.mode in ['smoke','scale']:
        data=root/'datasets/synthetic';data.mkdir(exist_ok=True);db=data/'hg38_synthetic.txt'
        if not db.exists():
            with open(db,'w') as f:
                f.write('#Chr\tStart\tEnd\tRef\tAlt\tvalue\n')
                for pos in range(1,100001):f.write(f'1\t{pos}\t{pos}\tA\tC\tvalue{pos}\n')
        measure('candidate',[binaries['candidate'],'db','index',db,'--kind','filter','--tmp-dir',root/'tmp'],1,'index-build',0)
        for size in ([100000] if args.mode=='smoke' else [1000000,10000000]):
            vcf=data/f'{size}.vcf'
            if not vcf.exists():
                with open(vcf,'w') as f:
                    f.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
                    for i in range(size):f.write(f'1\t{i%100000+1}\t.\tA\tC\t.\tPASS\t.\n')
            for threads in ([1,2,4,8,16,32] if args.mode=='smoke' else [16]):
                for name,binary in binaries.items():
                    for repeat in range(6):
                        out=run/f'{size}-{name}-t{threads}.tsv';cmd=[binary,'table',vcf,data,'--build','hg38','--protocol','synthetic','--operation','f','--vcf-input','--output',out]
                        if name=='candidate':cmd+=['--threads',str(threads),'--memory-budget','16','--report-json',str(out)+'.json']
                        measure(name,cmd,threads,f'synthetic-{size}',repeat)
            with open(run/f'{size}-baseline-t16.tsv') as a,open(run/f'{size}-candidate-t16.tsv') as b:
                for i,(x,y) in enumerate(itertools.zip_longest(a,b)):
                    if x is None or y is None or x.rstrip('\n').split('\t')!=y.rstrip('\n').split('\t')[:-1]:raise RuntimeError(f'Output mismatch at row {i}')
    else:
        model=root/'baseline/annovar/humandb'
        for assembly,build in [('GRCh37','hg19'),('GRCh38','hg38')]:
            vcf=root/f'datasets/giab/HG002_{assembly}_1_22_v4.2.1_benchmark.vcf.gz'
            missing=[str(p) for p in [vcf,model/f'{build}_refGene.txt',model/f'{build}_refGeneMrna.fa'] if not p.is_file()]
            if missing:(run/f'{build}-not-validated.json').write_text(json.dumps({'missing':missing},indent=2));continue
            for name,binary in binaries.items():
                for repeat in range(6):
                    out=run/f'{build}-{name}.tsv';cmd=[binary,'table',vcf,model,'--build',build,'--protocol','refGene','--operation','g','--vcf-input','--output',out]
                    if name=='candidate':cmd+=['--threads','16','--memory-budget','16']
                    measure(name,cmd,16,build,repeat)
            avinput=run/f'{build}.avinput';subprocess.run([binaries['candidate'],'convert',vcf,'--output',avinput],check=True)
            for repeat in range(6):measure('perl',['perl',root/'baseline/annovar/table_annovar.pl',avinput,model,'-buildver',build,'-protocol','refGene','-operation','g','-nastring','.','-out',run/f'{build}-perl','-remove'],16,build,repeat)
    summary=[]
    for name,threads,label in sorted(set((r['implementation'],r['threads'],r['label']) for r in rows if r['repeat']>0)):
        subset=[r for r in rows if (r['implementation'],r['threads'],r['label'])==(name,threads,label) and r['repeat']>0]
        summary.append({'implementation':name,'threads':threads,'label':label,'median_seconds':statistics.median(r['wall_seconds'] for r in subset),'min_seconds':min(r['wall_seconds'] for r in subset),'max_seconds':max(r['wall_seconds'] for r in subset),'max_rss_kib':max(r['peak_rss_kib'] for r in subset)})
    (run/'summary.json').write_text(json.dumps(summary,indent=2));(run/'SUCCESS').write_text('Selected mode completed. GIAB timings alone do not establish annotation equivalence.\n');print(run)
if __name__=='__main__':main()
