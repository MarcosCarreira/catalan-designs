import sys, math, time
exec(open('lai_sum.py').read().split('for n in map')[0])
for n in map(int,sys.argv[1].split(',')):
    T=time.time(); lo=math.isqrt(8*n)+1
    s=sum(max(bound44(n,k),-12*n) for k in range(lo,6*n+1))
    print(f'n={n}: (1/n^2) sum over all integers k in (sqrt(8n),6n] of the Lemma 4.4 bound = {s/n**2:.3f}  [{time.time()-T:.1f}s]',flush=True)
