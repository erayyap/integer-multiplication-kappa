// Laminar (binary-tree) side-wire groups per kernel; count side wires.
// usage: tree h mode   (stdin: h source petal priorities, then h target priorities)
// petal of a triple = point with smallest priority (source / target separately).
// For each kernel K: P_K = sorted petals; balanced binary tree over P_K; nodes = intervals.
// For every (K,K') pair: min partition of neighbour cells of P_K x Q_K' into node x node rectangles (DP).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define MAXH 128
int h, ps[MAXH], pt[MAXH];
int nP[MAXH*MAXH], P[MAXH*MAXH][MAXH], nQ[MAXH*MAXH], Q[MAXH*MAXH][MAXH];
int pre[MAXH+1][MAXH+1];
// tree nodes: intervals; build list per size n: node i -> lo,hi,left,right
typedef struct {int lo,hi,l,r;} Node;
Node T[MAXH+1][2*MAXH]; int nT[MAXH+1];
int build(int n,int lo,int hi){int id=nT[n]++;T[n][id].lo=lo;T[n][id].hi=hi;
 if(hi-lo>1){int mid=(lo+hi)/2; int l=build(n,lo,mid); int r=build(n,mid,hi);T[n][id].l=l;T[n][id].r=r;} else T[n][id].l=T[n][id].r=-1; return id;}
int memo[2*MAXH][2*MAXH]; unsigned stamp[2*MAXH][2*MAXH]; unsigned cur=0;
int na,nb_; Node *TA,*TB;
static inline int cnt(Node*a,Node*b){return pre[a->hi][b->hi]-pre[a->lo][b->hi]-pre[a->hi][b->lo]+pre[a->lo][b->lo];}
long long full_cells=0, fullrect=0, kkw=0,kkc=0,djw=0,djc=0,ppc=0;
int f(int ai,int bi){ if(stamp[ai][bi]==cur) return memo[ai][bi];
 Node*a=&TA[ai],*b=&TB[bi]; int c=cnt(a,b), area=(a->hi-a->lo)*(b->hi-b->lo), r;
 if(c==0) r=0; else if(c==area) r=1; else { r=1<<30;
   if(a->l>=0){int x=f(a->l,bi)+f(a->r,bi); if(x<r) r=x;}
   if(b->l>=0){int x=f(ai,b->l)+f(ai,b->r); if(x<r) r=x;} }
 stamp[ai][bi]=cur; memo[ai][bi]=r; return r;}
int rs=-1,rt=-1;
static int rk(int p,int a,int b){return (p>a)+(p>b);}
int main(int argc,char**argv){ if(argc>2){rs=atoi(argv[1]);rt=atoi(argv[2]);}
 if(scanf("%d",&h)!=1) return 1;
 for(int i=0;i<h;i++) scanf("%d",&ps[i]);
 for(int i=0;i<h;i++) scanf("%d",&pt[i]);
 for(int n=1;n<=h;n++) build(n,0,n);
 int nk=0; int ka[MAXH*MAXH],kb[MAXH*MAXH];
 for(int a=0;a<h;a++) for(int b=a+1;b<h;b++){ka[nk]=a;kb[nk]=b;
   nP[nk]=nQ[nk]=0;
   for(int p=0;p<h;p++){ if(p==a||p==b) continue;
     if(rs>=0? rk(p,a,b)==rs : (ps[p]<ps[a]&&ps[p]<ps[b])) P[nk][nP[nk]++]=p;
     if(rt>=0? rk(p,a,b)==rt : (pt[p]<pt[a]&&pt[p]<pt[b])) Q[nk][nQ[nk]++]=p; }
   nk++;}
 long long total=0, pairs=0;
 for(int i=0;i<nk;i++){ if(!nP[i]) continue; uint64_t Km=(1ULL<<ka[i])|(1ULL<<kb[i]);
  for(int j=0;j<nk;j++){ if(!nQ[j]) continue; uint64_t Lm=(1ULL<<ka[j])|(1ULL<<kb[j]);
   if(__builtin_popcountll(Km&Lm)==2) continue;
   int n1=nP[i], n2=nQ[j];
   for(int x=0;x<=n1;x++) pre[x][0]=0; for(int y=0;y<=n2;y++) pre[0][y]=0;
   for(int x=0;x<n1;x++){ uint64_t X=Km|(1ULL<<P[i][x]); int row=0;
     for(int y=0;y<n2;y++){ uint64_t S=Lm|(1ULL<<Q[j][y]); row+= (__builtin_popcountll(X&S)==1);
       pre[x+1][y+1]=pre[x][y+1]+row; } }
   if(pre[n1][n2]==0) continue;
   pairs+=pre[n1][n2];
   cur++; TA=T[n1]; TB=T[n2];
   {int w=f(0,0); total+=w; if(__builtin_popcountll(Km&Lm)==1){kkw+=w;kkc+=pre[n1][n2];} else {djw+=w; djc+=pre[n1][n2]; for(int x=0;x<n1;x++) for(int y=0;y<n2;y++) if(P[i][x]==Q[j][y]) ppc++;}}
  }}
 printf("%lld %lld kk %lld %lld dj %lld %lld pp %lld\n",total,pairs,kkw,kkc,djw,djc,ppc); return 0;}
