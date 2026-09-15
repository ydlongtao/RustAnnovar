#!/usr/bin/env python3
"""Compare ordered annotation tables; no unexplained differences are suppressed."""
import argparse, csv, itertools, json


def compare(actual, expected):
    counts = {}; examples = []; rows = 0
    with open(actual, newline='') as a, open(expected, newline='') as b:
        ar, br = csv.DictReader(a, delimiter='\t'), csv.DictReader(b, delimiter='\t')
        fields = [f for f in br.fieldnames if not f.startswith('Otherinfo')]
        missing = [f for f in fields if f not in ar.fieldnames]
        if missing: return {'passed':False,'missing_columns':missing}
        for rows, (x, y) in enumerate(itertools.zip_longest(ar, br), 1):
            if x is None or y is None:
                counts['row_count'] = counts.get('row_count',0)+1
                continue
            for field in fields:
                av, bv = x[field], y[field]
                if field.startswith(('Gene.','AAChange.')):
                    av, bv = sorted(av.split(',')), sorted(bv.split(','))
                if av != bv:
                    counts[field] = counts.get(field,0)+1
                    if len(examples)<50: examples.append({'row':rows,'coordinate':[x.get(k) for k in ['Chr','Start','End','Ref','Alt']],'field':field,'actual':x[field],'expected':y[field]})
    return {'passed':not counts,'rows':rows,'differences':counts,'examples':examples,'ignored':'Otherinfo columns; order within Gene/AAChange comma lists only'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('actual');p.add_argument('expected');p.add_argument('--output',required=True);args=p.parse_args()
    result=compare(args.actual,args.expected)
    with open(args.output,'w') as f:json.dump(result,f,indent=2);f.write('\n')
    raise SystemExit(0 if result['passed'] else 1)
