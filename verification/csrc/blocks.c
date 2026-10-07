// Count side wires for a (source grouping, target grouping) of triples.
// Input (stdin): h  nS  then nS groups: size m1..  ; nT then groups.
// Each group: size followed by triple bitmasks (uint64).
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned long long u64;
static int readgroups(u64 ***g, int **sz){
  int n; if(scanf("%d",&n)!=1) exit(1);
  *g = malloc(sizeof(u64*)*n); *sz = malloc(sizeof(int)*n);
  for(int i=0;i<n;i++){ int m; scanf("%d",&m); (*sz)[i]=m; (*g)[i]=malloc(sizeof(u64)*m);
    for(int j=0;j<m;j++) scanf("%llu",&(*g)[i][j]); }
  return n;
}
int main(){
  int h; scanf("%d",&h);
  u64 **S,**T; int *ss,*ts;
  int nS = readgroups(&S,&ss), nT = readgroups(&T,&ts);
  long long total=0, complete=0, rowc=0, colc=0, sing=0, pairs=0, cpairs=0, linewires=0, linecells=0;
  int nb[64][64];
  for(int a=0;a<nS;a++){
    for(int b=0;b<nT;b++){
      int ma=ss[a], mb=ts[b], n=0;
      for(int i=0;i<ma;i++) for(int j=0;j<mb;j++){ int c=__builtin_popcountll(S[a][i]&T[b][j]); nb[i][j]=(c==1); n+=nb[i][j]; }
      if(!n) continue;
      pairs+=n;
      if(n==ma*mb){ total+=1; complete++; cpairs+=n; continue; }
      // rows
      int fr=0, rest=0; for(int i=0;i<ma;i++){ int s=0; for(int j=0;j<mb;j++) s+=nb[i][j]; if(s==mb) fr++; else rest+=s; }
      int cr=fr+rest;
      int fc=0, restc=0; for(int j=0;j<mb;j++){ int s=0; for(int i=0;i<ma;i++) s+=nb[i][j]; if(s==ma) fc++; else restc+=s; }
      int cc=fc+restc;
      int best=n; if(cr<best) best=cr; if(cc<best) best=cc;
      int cm = fr + fc + n - fr*mb - fc*ma + 2*fr*fc; if(cm<best) best=cm;
      total+=best;
      if(best==n) sing+=n; else if(best==cr){linewires+=fr; linecells+=fr*mb; sing+=rest;} else if(best==cc){linewires+=fc; linecells+=fc*ma; sing+=restc;} else {linewires+=fr+fc; sing+=best-fr-fc;}
    }
  }
  printf("%lld %lld %lld %lld %lld %lld %lld\n", total, pairs, complete, cpairs, linewires, linecells, sing);
  return 0;
}
