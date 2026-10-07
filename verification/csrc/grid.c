// Grid hierarchies: element X=(l<c<u). Source grids G_c = [0,c) x (c,h); nodes = interval products.
// Tree policy per side: 'H' split hi-range to singletons first (top levels), then lo; 'L' the reverse;
// 'K' split the longer side.  Wires = sum over (c,c') of min partition of neighbour cells into node x node.
// usage: grid h polS polT [nofam]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef struct {int l0,l1,u0,u1,a,b;} Node;
int h; char polS, polT;
Node *TS, *TT; int nS, nTt;
static int build(Node*T,int*n,int l0,int l1,int u0,int u1,char pol){
  int id=(*n)++; T[id]=(Node){l0,l1,u0,u1,-1,-1};
  int dl=l1-l0, du=u1-u0; if(dl*du<=1) return id;
  int splitU;
  if(pol=='H') splitU = du>1; else if(pol=='L') splitU = !(dl>1); else splitU = du>=dl;
  int a,b;
  if(splitU){int m=(u0+u1)/2; a=build(T,n,l0,l1,u0,m,pol); b=build(T,n,l0,l1,m,u1,pol);}
  else {int m=(l0+l1)/2; a=build(T,n,l0,m,u0,u1,pol); b=build(T,n,m,l1,u0,u1,pol);}
  T[id].a=a; T[id].b=b; return id; }
int *pre; int D1,D2,D3,D4; // dims+1
#define IDX(a,b,c,d) ((((size_t)(a)*D2+(b))*D3+(c))*D4+(d))
static inline long long box(const Node*x,const Node*y){
  long long s=0;
  for(int m=0;m<16;m++){
    int A=(m&1)?x->l1:x->l0, B=(m&2)?x->u1:x->u0, C=(m&4)?y->l1:y->l0, E=(m&8)?y->u1:y->u0;
    int sg=((m&1)!=0)+((m&2)!=0)+((m&4)!=0)+((m&8)!=0);
    s+= (sg&1? -1:1) * pre[IDX(A,B,C,E)];
  }
  return (s<0?-s:s) ; // sign: with 4 "upper" terms positive
}
int *memo; unsigned *stamp; unsigned cur=0; int MAXN;
int uS0,uT0; // offsets of u ranges (c+1)
static int f(int ai,int bi){
  size_t k=(size_t)ai*MAXN+bi; if(stamp[k]==cur) return memo[k];
  const Node*a=&TS[ai],*b=&TT[bi];
  long long c=box(a,b); long long area=(long long)(a->l1-a->l0)*(a->u1-a->u0)*(b->l1-b->l0)*(b->u1-b->u0); int r;
  if(c==0) r=0; else if(c==area) r=1; else { r=1<<30;
    if(a->a>=0){int x=f(a->a,bi)+f(a->b,bi); if(x<r) r=x;}
    if(b->a>=0){int x=f(ai,b->a)+f(ai,b->b); if(x<r) r=x;} }
  stamp[k]=cur; memo[k]=r; return r; }
int main(int argc,char**argv){
  h=atoi(argv[1]); polS=argv[2][0]; polT=argv[3][0];
  MAXN=2*(h*h)+4;
  TS=malloc(sizeof(Node)*MAXN); TT=malloc(sizeof(Node)*MAXN);
  memo=malloc(sizeof(int)*(size_t)MAXN*MAXN); stamp=calloc((size_t)MAXN*MAXN,sizeof(unsigned));
  pre=malloc(sizeof(int)*(size_t)(h+1)*(h+1)*(h+1)*(h+1));
  long long total=0,pairs=0;
  for(int c=1;c<h-1;c++){
    nS=0; build(TS,&nS,0,c,0,h-1-c,polS); // u index: u = c+1+j
    for(int d=1;d<h-1;d++){
      nTt=0; build(TT,&nTt,0,d,0,h-1-d,polT);
      D1=c+1; D2=h-c; D3=d+1; D4=h-d;
      // 4D prefix sums
      for(int A=0;A<D1;A++) for(int B=0;B<D2;B++) for(int C=0;C<D3;C++) for(int E=0;E<D4;E++){
        int v=0;
        if(A&&B&&C&&E){ uint64_t X=(1ULL<<(A-1))|(1ULL<<c)|(1ULL<<(c+B)); uint64_t S=(1ULL<<(C-1))|(1ULL<<d)|(1ULL<<(d+E));
          v=(__builtin_popcountll(X&S)==1);
          // inclusion-exclusion over lower corners
          for(int m=1;m<16;m++){ int a2=A-((m&1)!=0), b2=B-((m&2)!=0), c2=C-((m&4)!=0), e2=E-((m&8)!=0);
            int sg=((m&1)!=0)+((m&2)!=0)+((m&4)!=0)+((m&8)!=0);
            v+= (sg&1?1:-1)*pre[IDX(a2,b2,c2,e2)]; }
        }
        pre[IDX(A,B,C,E)]=v; }
      long long np=pre[IDX(c,h-1-c,d,h-1-d)]; if(!np) continue; pairs+=np;
      cur++; total+=f(0,0);
    }
  }
  printf("%lld %lld\n",total,pairs); return 0; }
