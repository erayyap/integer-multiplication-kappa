// Simulated annealing over per-triple petal choices (source petal, target petal).
// Families = kernel pairs; tree over sorted petals; cost(i,j) = DP min partition.
// usage: sa_roles h iters T0 seed [init: anti|same]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#define MH 64
int h, nk; int kid[MH][MH]; int ka[MH*MH], kb[MH*MH];
int tid[MH][MH][MH];
int nt; int ta[50000],tb[50000],tc[50000]; int ps[50000], pt[50000];
// family membership: inS[k][p] = 1 if triple k∪{p} has source petal p
unsigned char inS[MH*MH/2+64][MH], inT[MH*MH/2+64][MH];
int *cost; // nk*nk
typedef struct {int lo,hi,l,r;} Node;
Node TA[2*MH],TB[2*MH]; int nA,nB;
int build(Node*T,int*n,int lo,int hi){int id=(*n)++;T[id].lo=lo;T[id].hi=hi;if(hi-lo>1){int m=(lo+hi)/2;int a=build(T,n,lo,m),b=build(T,n,m,hi);T[id].l=a;T[id].r=b;}else T[id].l=T[id].r=-1;return id;}
int pre[MH+1][MH+1]; int memo[2*MH][2*MH]; unsigned stamp[2*MH][2*MH]; unsigned cur=0;
int f(int ai,int bi){ if(stamp[ai][bi]==cur) return memo[ai][bi]; Node*a=&TA[ai],*b=&TB[bi];
 int c=pre[a->hi][b->hi]-pre[a->lo][b->hi]-pre[a->hi][b->lo]+pre[a->lo][b->lo]; int area=(a->hi-a->lo)*(b->hi-b->lo),r;
 if(c==0) r=0; else if(c==area) r=1; else { r=1<<30; if(a->l>=0){int x=f(a->l,bi)+f(a->r,bi); if(x<r)r=x;} if(b->l>=0){int x=f(ai,b->l)+f(ai,b->r); if(x<r)r=x;} }
 stamp[ai][bi]=cur; memo[ai][bi]=r; return r; }
int P[MH],Q[MH];
int pair_cost(int i,int j){
  uint64_t Km=(1ULL<<ka[i])|(1ULL<<kb[i]), Lm=(1ULL<<ka[j])|(1ULL<<kb[j]);
  if(__builtin_popcountll(Km&Lm)==2) return 0;
  int n1=0,n2=0; for(int p=0;p<h;p++){ if(inS[i][p]) P[n1++]=p; if(inT[j][p]) Q[n2++]=p; }
  if(!n1||!n2) return 0;
  for(int x=0;x<=n1;x++) pre[x][0]=0; for(int y=0;y<=n2;y++) pre[0][y]=0;
  for(int x=0;x<n1;x++){ uint64_t X=Km|(1ULL<<P[x]); int row=0; for(int y=0;y<n2;y++){ uint64_t S=Lm|(1ULL<<Q[y]); row+=(__builtin_popcountll(X&S)==1); pre[x+1][y+1]=pre[x][y+1]+row; } }
  if(pre[n1][n2]==0) return 0;
  nA=0; build(TA,&nA,0,n1); nB=0; build(TB,&nB,0,n2); cur++; return f(0,0); }
int other(int t,int p,int *a,int *b){ int v[3]={ta[t],tb[t],tc[t]}; int k=0; int o[2]; for(int i=0;i<3;i++) if(v[i]!=p) o[k++]=v[i]; *a=o[0]; *b=o[1]; return kid[o[0]][o[1]]; }
double urand(){ return rand()/(RAND_MAX+1.0); }
int rowbuf[MH*MH], rowbuf2[MH*MH];
int main(int argc,char**argv){
  h=atoi(argv[1]); long iters=atol(argv[2]); double T0=atof(argv[3]); srand(atoi(argv[4])); int same = argc>5 && !strcmp(argv[5],"same");
  nk=0; for(int a=0;a<h;a++) for(int b=a+1;b<h;b++){ kid[a][b]=kid[b][a]=nk; ka[nk]=a; kb[nk]=b; nk++; }
  nt=0; for(int a=0;a<h;a++) for(int b=a+1;b<h;b++) for(int c=b+1;c<h;c++){ ta[nt]=a;tb[nt]=b;tc[nt]=c; ps[nt]=a; pt[nt]= same? a : c; nt++; }
  for(int t=0;t<nt;t++){ int x,y; int k=other(t,ps[t],&x,&y); inS[k][ps[t]]=1; k=other(t,pt[t],&x,&y); inT[k][pt[t]]=1; }
  cost=malloc(sizeof(int)*nk*nk); long long total=0;
  for(int i=0;i<nk;i++) for(int j=0;j<nk;j++){ cost[i*nk+j]=pair_cost(i,j); total+=cost[i*nk+j]; }
  printf("init %lld\n",total); fflush(stdout); long long best=total;
  for(long it=0;it<iters;it++){
    double T=T0*(1.0-(double)it/iters);
    int t=rand()%nt; int side=rand()%2; int v[3]={ta[t],tb[t],tc[t]};
    int old = side? pt[t]:ps[t]; int np; do{ np=v[rand()%3]; }while(np==old);
    int x,y; int k_old=other(t,old,&x,&y), k_new=other(t,np,&x,&y);
    // apply
    if(!side){ inS[k_old][old]=0; inS[k_new][np]=1; ps[t]=np; } else { inT[k_old][old]=0; inT[k_new][np]=1; pt[t]=np; }
    long long delta=0;
    int ks[2]={k_old,k_new};
    for(int u=0;u<2;u++){ int k=ks[u]; int *buf= u? rowbuf2:rowbuf;
      for(int j=0;j<nk;j++){ int c = side? pair_cost(j,k) : pair_cost(k,j); int oldc = side? cost[j*nk+k] : cost[k*nk+j]; buf[j]=c; delta+=c-oldc; } }
    if(delta<=0 || (T>0 && urand()<exp(-delta/T))){
      for(int u=0;u<2;u++){ int k=ks[u]; int *buf= u? rowbuf2:rowbuf; for(int j=0;j<nk;j++){ if(side) cost[j*nk+k]=buf[j]; else cost[k*nk+j]=buf[j]; } }
      total+=delta; if(total<best){best=total;}
    } else {
      if(!side){ inS[k_new][np]=0; inS[k_old][old]=1; ps[t]=old; } else { inT[k_new][np]=0; inT[k_old][old]=1; pt[t]=old; }
    }
    if(it%20000==0){ printf("it %ld total %lld best %lld\n",it,total,best); fflush(stdout); }
  }
  // verify
  long long chk=0; for(int i=0;i<nk;i++) for(int j=0;j<nk;j++) chk+=pair_cost(i,j);
  printf("final %lld check %lld best %lld\n",total,chk,best);
  FILE*fo=fopen(argv[6]?argv[6]:"/dev/null","w"); if(fo){ for(int t=0;t<nt;t++) fprintf(fo,"%d %d %d %d %d\n",ta[t],tb[t],tc[t],ps[t],pt[t]); fclose(fo);} 
  return 0; }
