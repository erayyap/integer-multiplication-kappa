// General block hierarchies.  usage: gen h SIDE_S SIDE_T
// SIDE = Fr  : family blocks: kernel pair fixed, varying point = rank r (0 min,1 mid,2 max) in triple; 1D binary tree.
//        Sr<pol>: slice blocks: point of rank r fixed, other two points (a<b) on 2D grid; k-d tree with policy
//                 pol in {A: split a first, B: split b first, K: longer side}.
// Each element lies in exactly one block per side.  Wires = sum over block pairs of the min partition (DP).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef struct {int a0,a1,b0,b1,l,r; int valid;} Node;
typedef struct {int n0,n1; uint64_t *mask; /* n0*n1, 0 = invalid */ } Block;
int h;
Block *mk(char kind,int r,int *nb){
  Block *B=malloc(sizeof(Block)*h*h); int n=0;
  if(kind=='F'){
    for(int x=0;x<h;x++) for(int y=x+1;y<h;y++){
      int cnt=0; uint64_t tmp[128];
      for(int p=0;p<h;p++){ if(p==x||p==y) continue; int rk=(p>x)+(p>y); if(rk==r) tmp[cnt++]=(1ULL<<x)|(1ULL<<y)|(1ULL<<p);} 
      if(!cnt) continue; B[n].n0=cnt; B[n].n1=1; B[n].mask=malloc(8*cnt); memcpy(B[n].mask,tmp,8*cnt); n++; }
  } else {
    for(int q=0;q<h;q++){
      int pts[128],np=0; // points available for the pair
      int lo[128],nl=0,hi[128],nh=0;
      for(int p=0;p<h;p++){ if(p<q) lo[nl++]=p; else if(p>q) hi[nh++]=p; }
      int *A,*Bv,na,nbv;
      if(r==0){A=hi;na=nh;Bv=hi;nbv=nh;} else if(r==2){A=lo;na=nl;Bv=lo;nbv=nl;} else {A=lo;na=nl;Bv=hi;nbv=nh;}
      if(!na||!nbv) continue;
      B[n].n0=na; B[n].n1=nbv; B[n].mask=calloc(na*nbv,8); int any=0;
      for(int i=0;i<na;i++) for(int j=0;j<nbv;j++){ if(A[i]<Bv[j]){ B[n].mask[i*nbv+j]=(1ULL<<q)|(1ULL<<A[i])|(1ULL<<Bv[j]); any=1;} }
      if(any) n++; }
  }
  *nb=n; return B; }
// tree
static int build(Node*T,int*n,int a0,int a1,int b0,int b1,char pol,int *vpre,int n1){
  // count valid in box via vpre (2D prefix of validity, dims (n0+1)x(n1+1))
  #define VP(i,j) vpre[(i)*(n1+1)+(j)]
  int v=VP(a1,b1)-VP(a0,b1)-VP(a1,b0)+VP(a0,b0);
  if(v==0) return -1;
  int id=(*n)++; T[id]=(Node){a0,a1,b0,b1,-1,-1,v};
  if(v==1) return id;
  int da=a1-a0, db=b1-b0, sa;
  if(pol=='A') sa = da>1; else if(pol=='B') sa = !(db>1); else sa = da>=db;
  int L,R;
  if(sa){int m=(a0+a1)/2; L=build(T,n,a0,m,b0,b1,pol,vpre,n1); R=build(T,n,m,a1,b0,b1,pol,vpre,n1);}
  else {int m=(b0+b1)/2; L=build(T,n,a0,a1,b0,m,pol,vpre,n1); R=build(T,n,a0,a1,m,b1,pol,vpre,n1);}
  if(L<0) return R>=0? (T[id]=T[R], id) : id; // collapse: single child -> same set
  if(R<0) { T[id]=T[L]; return id; }
  T[id].l=L; T[id].r=R; return id; }
int *pre; int D1,D2,D3,D4;
#define IDX(a,b,c,d) ((((size_t)(a)*D2+(b))*D3+(c))*D4+(d))
static inline long long box(const Node*x,const Node*y){
  long long s=0;
  for(int m=0;m<16;m++){
    int A=(m&1)?x->a1:x->a0, B=(m&2)?x->b1:x->b0, C=(m&4)?y->a1:y->a0, E=(m&8)?y->b1:y->b0;
    int sg=((m&1)!=0)+((m&2)!=0)+((m&4)!=0)+((m&8)!=0);
    s+= (sg&1? -1:1) * pre[IDX(A,B,C,E)]; }
  return s; }
Node *TS,*TT; int *memo; unsigned *stamp; unsigned cur=0; int MAXN;
static int f(int ai,int bi){
  size_t k=(size_t)ai*MAXN+bi; if(stamp[k]==cur) return memo[k];
  const Node*a=&TS[ai],*b=&TT[bi];
  long long c=box(a,b), area=(long long)a->valid*b->valid; int r;
  if(c==0) r=0; else if(c==area) r=1; else { r=1<<30;
    if(a->l>=0){int x=f(a->l,bi)+f(a->r,bi); if(x<r) r=x;}
    if(b->l>=0){int x=f(ai,b->l)+f(ai,b->r); if(x<r) r=x;} }
  stamp[k]=cur; memo[k]=r; return r; }
int main(int argc,char**argv){
  h=atoi(argv[1]); char ks=argv[2][0], kt=argv[3][0]; int rs=argv[2][1]-'0', rt=argv[3][1]-'0';
  char ps=argv[2][2]?argv[2][2]:'K', pt=argv[3][2]?argv[3][2]:'K';
  int nbS,nbT; Block*BS=mk(ks,rs,&nbS),*BT=mk(kt,rt,&nbT);
  MAXN=2*h*h+8; TS=malloc(sizeof(Node)*MAXN); TT=malloc(sizeof(Node)*MAXN);
  memo=malloc(sizeof(int)*(size_t)MAXN*MAXN); stamp=calloc((size_t)MAXN*MAXN,sizeof(unsigned));
  pre=malloc(sizeof(int)*(size_t)(h+1)*(h+1)*(h+1)*(h+1));
  int *vS=malloc(sizeof(int)*(h+1)*(h+1)), *vT=malloc(sizeof(int)*(h+1)*(h+1));
  long long total=0,pairs=0;
  for(int i=0;i<nbS;i++){
    Block*X=&BS[i]; int n0=X->n0,n1=X->n1;
    for(int a=0;a<=n0;a++) for(int b=0;b<=n1;b++) vS[a*(n1+1)+b]= (a&&b)? ((X->mask[(a-1)*n1+b-1]!=0) + vS[(a-1)*(n1+1)+b]+vS[a*(n1+1)+b-1]-vS[(a-1)*(n1+1)+b-1]) : 0;
    int nS=0; build(TS,&nS,0,n0,0,n1,ps,vS,n1);
    for(int j=0;j<nbT;j++){
      Block*Y=&BT[j]; int m0=Y->n0,m1=Y->n1;
      D1=n0+1; D2=n1+1; D3=m0+1; D4=m1+1;
      for(int A=0;A<D1;A++) for(int B=0;B<D2;B++) for(int C=0;C<D3;C++) for(int E=0;E<D4;E++){
        int v=0;
        if(A&&B&&C&&E){ uint64_t x=X->mask[(A-1)*n1+B-1], y=Y->mask[(C-1)*m1+E-1];
          v=(x&&y&&__builtin_popcountll(x&y)==1);
          for(int m=1;m<16;m++){ int a2=A-((m&1)!=0), b2=B-((m&2)!=0), c2=C-((m&4)!=0), e2=E-((m&8)!=0);
            int sg=((m&1)!=0)+((m&2)!=0)+((m&4)!=0)+((m&8)!=0);
            v+= (sg&1?1:-1)*pre[IDX(a2,b2,c2,e2)]; } }
        pre[IDX(A,B,C,E)]=v; }
      long long np=pre[IDX(n0,n1,m0,m1)]; if(!np) continue; pairs+=np;
      for(int a=0;a<=m0;a++) for(int b=0;b<=m1;b++) vT[a*(m1+1)+b]= (a&&b)? ((Y->mask[(a-1)*m1+b-1]!=0) + vT[(a-1)*(m1+1)+b]+vT[a*(m1+1)+b-1]-vT[(a-1)*(m1+1)+b-1]) : 0;
      int nT=0; build(TT,&nT,0,m0,0,m1,pt,vT,m1);
      cur++; total+=f(0,0);
    } }
  printf("%lld %lld\n",total,pairs); return 0; }
