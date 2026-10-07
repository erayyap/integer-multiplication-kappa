import sorted_supp as ss
def structure(n):
    st = ss.structure(n)
    allp = frozenset(range(n))
    st['gsupp'] = {g: (A | {min(allp - A)} if len(A) < n else A) for g, A in st['gsupp'].items()}
    return st
