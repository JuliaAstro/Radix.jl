#!/usr/bin/env python3
"""Build a standalone reference for XSTAR's ucalc type 63 (CollisionProb).

Usage: extract_type63.py /path/to/xstarsub.f

Writes routines63.f (the Fortran units that type 63 depends on, extracted from
xstarsub.f) and prob63.f (the body of the ucalc `63 continue` block wrapped in
a subroutine). Together with drv.f they build with, e.g.,

    gfortran -std=legacy -w -fno-automatic -fdefault-real-8 -fdefault-double-8 \
        -o ref63 drv.f prob63.f routines63.f

drv.f reads lines `ni li nf lf Z T[K] E1 E2 g1 g2 ne` from stdin and prints
`ans1 ans2` (including the factor ne).
"""
import re, sys

src = open(sys.argv[1], errors='replace').read().split('\n')

# --- program units -----------------------------------------------------------
unit_re = re.compile(r'^\s{6}\s*(?:real\*8\s+|double\s+precision\s+|real\s+|integer\s+)?(subroutine|function)\s+(\w+)', re.I)
units, order, cur = {}, [], None
for i, l in enumerate(src):
    if l[:1] in 'cC*!':
        continue
    m = unit_re.match(l)
    if m and l.lstrip().lower().startswith(('subroutine', 'function', 'real function',
            'real*8 function', 'double precision function', 'integer function')):
        cur = m.group(2).lower(); units[cur] = [i, None]; order.append(cur)
    if cur and re.match(r'^\s+end\s*$', l, re.I):
        units[cur][1] = i; cur = None

names = set(units)
def deps(name):
    s, e = units[name]
    body = '\n'.join(l for l in src[s:e+1] if l[:1] not in 'cC*!').lower()
    out = {m.group(1) for m in re.finditer(r'call\s+(\w+)', body)}
    out |= {n for n in names if n != name and re.search(r'\b' + n + r'\s*\(', body)}
    return {d for d in out if d in names}

need, stack = set(), ['anl1', 'amcrs', 'erc']
while stack:
    n = stack.pop()
    if n not in need:
        need.add(n); stack.extend(deps(n))
open('routines63.f', 'w').write('\n'.join('\n'.join(src[units[n][0]:units[n][1]+1]) for n in order if n in need) + '\n')

# --- the ucalc block ---------------------------------------------------------
s0 = next(i for i, l in enumerate(src) if re.match(r'^\s*subroutine ucalc\(', l))
b = next(i for i in range(s0, len(src)) if re.match(r'^ 63\s+continue', src[i]))
e = next(i for i in range(b + 1, len(src)) if re.match(r'^ 64\s+continue', src[i]))
body = re.sub(r'go to 9000', 'return', '\n'.join(src[b+1:e]))
hdr = '''      subroutine prob63(idat,nidat,ilev,rlev,nlev,t,xnx,lpri,lun11,
     $                  ans1,ans2)
      implicit none
      integer idat(10),nidat,ilev(6,3),nlev,lpri,lun11
      real rlev(4,3),t,xnx,ans1,ans2
      integer idest1,idest2,lpril,ni,li,nf,lf,iq,lff,lii,li1,nn,nnz
      integer il,iz
      integer nll,nu
      real eeup,eelo,elin,hij,ekt,delt,sum,alm,alp,ecm,tbig,z1,rm,psi
      real cn,cno,aa1,se,sd,atmp
      ans1=0.
      ans2=0.
'''
open('prob63.f', 'w').write(hdr + body + '\n      end\n')
