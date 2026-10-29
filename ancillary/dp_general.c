// Exact Bellman fixed point for Younis chains at a squarefree modulus m (m <= 127) with -1 a square mod m.
// States: subsets of the negation pairs {s,-s} of nonzero squares.  T(F,F') = clique number of the Cayley
// graph on Z/mZ forbidding F u D(F'), D(F') = {d != 0 : d^2 not in F'}.  W(F) = max_{F'} log T(F,F') + W(F')/2.
// Usage: ./dp_general m   -> prints W(S), gamma, and the optimal policy from S.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

typedef unsigned __int128 u128;
static int M;
static u128 adj[128];
static int best;
static inline int ctz128(u128 x){ uint64_t lo=(uint64_t)x; if(lo) return __builtin_ctzll(lo); return 64+__builtin_ctzll((uint64_t)(x>>64)); }

static void expand(u128 P, int size){
    if(P==0){ if(size>best) best=size; return; }
    int order[128], color[128], n=0; u128 Q=P; int col=0;
    while(Q){ col++; u128 R=Q;
        while(R){ int v=ctz128(R); R&=~((u128)1<<v); R&=~adj[v]; Q&=~((u128)1<<v); order[n]=v; color[n]=col; n++; } }
    for(int i=n-1;i>=0;i--){ if(size+color[i]<=best) return; int v=order[i]; expand(P&adj[v],size+1); P&=~((u128)1<<v); }
}
static int omega(const unsigned char *forb){ // forb[d]=1 if difference d forbidden (symmetric), d in 1..M-1
    int verts[128],k=0; for(int a=1;a<M;a++) if(!forb[a]) verts[k++]=a;
    for(int i=0;i<k;i++){ adj[i]=0; for(int j=0;j<k;j++){ if(i==j) continue; int d=((verts[i]-verts[j])%M+M)%M; if(!forb[d]) adj[i]|=((u128)1<<j);} }
    best=0; u128 P=(k==128)?~(u128)0:(((u128)1<<k)-1); expand(P,0); return best+1;
}
// simple open-addressing hash set on u128 keys
static u128 *hkeys; static int *hval; static size_t hcap; static size_t hcount;
static size_t hfind(u128 key){ uint64_t h=(uint64_t)key*0x9E3779B97F4A7C15ULL ^ (uint64_t)(key>>64)*0xC2B2AE3D27D4EB4FULL; size_t i=h&(hcap-1);
    while(hval[i]!=-1 && hkeys[i]!=key) i=(i+1)&(hcap-1); return i; }

int main(int argc,char**argv){
    M=atoi(argv[1]);
    int issq[128]={0}; for(int x=0;x<M;x++) issq[(x*x)%M]=1;
    int pr[64],npairs=0; int used[128]={0};
    for(int s=1;s<M;s++) if(issq[s]&&!used[s]){ int t=(M-s)%M; if(!issq[t]){fprintf(stderr,"-1 not a square\n");return 1;} pr[npairs++]=s; used[s]=used[t]=1; }
    int NS=1<<npairs; fprintf(stderr,"m=%d pairs=%d states=%d\n",M,npairs,NS);
    // D(F') masks and F masks as u128 over residues
    u128 *Fm=malloc(sizeof(u128)*NS), *Dm=malloc(sizeof(u128)*NS);
    for(int F=0;F<NS;F++){ u128 fm=0; for(int i=0;i<npairs;i++) if(F>>i&1){ fm|=(u128)1<<pr[i]; fm|=(u128)1<<((M-pr[i])%M);} Fm[F]=fm; }
    for(int Fp=0;Fp<NS;Fp++){ u128 dm=0; for(int d=1;d<M;d++){ int s=(d*d)%M; if(!((Fm[Fp]>>s)&1)) dm|=(u128)1<<d; } Dm[Fp]=dm; }
    hcap=1; while(hcap<(size_t)NS*(size_t)NS/4+1024) hcap<<=1; if(hcap>(1ULL<<27)) hcap=1ULL<<27;
    hkeys=malloc(sizeof(u128)*hcap); hval=malloc(sizeof(int)*hcap); for(size_t i=0;i<hcap;i++) hval[i]=-1; hcount=0;
    unsigned char *T=malloc((size_t)NS*NS);
    size_t solved=0;
    for(int F=0;F<NS;F++){ for(int Fp=0;Fp<NS;Fp++){ u128 B=Fm[F]|Dm[Fp]; size_t i=hfind(B);
            if(hval[i]==-1){ unsigned char forb[128]; memset(forb,0,sizeof forb); for(int d=1;d<M;d++) forb[d]=(B>>d)&1; hval[i]=omega(forb); hkeys[i]=B; hcount++; solved++;
                if(hcount>hcap*3/4){fprintf(stderr,"hash full\n");return 1;} }
            T[(size_t)F*NS+Fp]=(unsigned char)hval[i]; }
        if((F&255)==0) fprintf(stderr,"F=%d distinct=%zu\n",F,hcount); }
    fprintf(stderr,"distinct forbidden sets solved: %zu\n",hcount);
    // value iteration
    double *W=calloc(NS,sizeof(double)), *Wn=calloc(NS,sizeof(double)); int *pol=malloc(sizeof(int)*NS);
    for(int it=0;it<200;it++){ double diff=0;
        for(int F=0;F<NS;F++){ double bestv=-1e9; int bp=0; const unsigned char *row=T+(size_t)F*NS;
            for(int Fp=0;Fp<NS;Fp++){ double v=log((double)row[Fp])+0.5*W[Fp]; if(v>bestv){bestv=v;bp=Fp;} }
            Wn[F]=bestv; pol[F]=bp; double d=fabs(bestv-W[F]); if(d>diff) diff=d; }
        memcpy(W,Wn,sizeof(double)*NS); if(diff<1e-13){ fprintf(stderr,"converged after %d iterations\n",it+1); break; } }
    int S=NS-1; double gamma=0.5+W[S]/(4*log((double)M));
    printf("m=%d W(S)=%.15f gamma=%.15f\n",M,W[S],gamma);
    int F=S; printf("policy from S (|R|, next-state pair count):");
    for(int step=0;step<12;step++){ int Fp=pol[F]; printf(" (%d,%d)",T[(size_t)F*NS+Fp],__builtin_popcount(Fp)); F=Fp; }
    printf("\n");
    FILE *f=fopen(argc>2?argv[2]:"T_general.bin","wb"); fwrite(T,1,(size_t)NS*NS,f); fclose(f);
    return 0;
}
