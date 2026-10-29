// Maximum cliques in Cayley graphs on Z/mZ with a symmetric forbidden difference set, m <= 1024.
// Vertex 0 is fixed in the clique (vertex-transitivity), so instances have at most m-1 vertices.
// Input (stdin):  m
//                 then per instance: mode k d1 ... dk   (forbidden differences, symmetric; mode: w = omega only,
//                 c = omega + one maximum clique, eK = enumerate up to K maximum cliques containing 0)
// Output: per instance one line "omega" (w), "omega v1 v2 ..." (c), or for e: first line "omega count", then one clique per line.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define MAXV 1024
#define NW 16                       // 16*64 = 1024 bits
typedef struct { uint64_t w[NW]; } bs;
static int M, K, nw;
static bs adj[MAXV];
static int verts[MAXV];
static int best; static int bestset[MAXV]; static int cur[MAXV];
static int enumtarget, enumcap, enumcount; static FILE *out;
static long long nodes, budget;

static inline int bs_empty(const bs *a){ for(int i=0;i<nw;i++) if(a->w[i]) return 0; return 1; }
static inline int bs_ctz(const bs *a){ for(int i=0;i<nw;i++) if(a->w[i]) return 64*i+__builtin_ctzll(a->w[i]); return -1; }
static inline void bs_clr(bs *a,int v){ a->w[v>>6] &= ~(1ULL<<(v&63)); }
static inline void bs_and(bs *r,const bs *a,const bs *b){ for(int i=0;i<nw;i++) r->w[i]=a->w[i]&b->w[i]; }
static inline void bs_andnot(bs *r,const bs *a,const bs *b){ for(int i=0;i<nw;i++) r->w[i]=a->w[i]&~b->w[i]; }

static void expand(bs P, int size){
    nodes++;
    if(bs_empty(&P)){
        if(enumtarget){ if(size==enumtarget && enumcount<enumcap){ enumcount++; fprintf(out,"0"); for(int i=0;i<size;i++) fprintf(out," %d",verts[cur[i]]); fprintf(out,"\n"); } }
        else if(size>best){ best=size; memcpy(bestset,cur,sizeof(int)*size); }
        return;
    }
    int order[MAXV], color[MAXV], n=0; bs Q=P; int col=0;
    while(!bs_empty(&Q)){
        col++; bs R=Q;
        while(!bs_empty(&R)){ int v=bs_ctz(&R); bs_clr(&R,v); bs_andnot(&R,&R,&adj[v]); bs_clr(&Q,v); order[n]=v; color[n]=col; n++; }
    }
    int bound = enumtarget ? enumtarget-1 : best;
    for(int i=n-1;i>=0;i--){
        if(size+color[i]<=bound) return;
        if(enumtarget && enumcount>=enumcap) return;
        if(budget>0 && nodes>budget) return;
        int v=order[i]; cur[size]=v; bs NP; bs_and(&NP,&P,&adj[v]); expand(NP,size+1); bs_clr(&P,v);
    }
}

int main(void){
    out=stdout;
    if(scanf("%d %lld",&M,&budget)!=2) return 1;
    static unsigned char forb[MAXV];
    char mode[8];
    while(scanf("%7s %d",mode,&K)==2){
        memset(forb,0,sizeof forb);
        for(int i=0;i<K;i++){ int d; if(scanf("%d",&d)!=1) return 1; d%=M; if(d<0) d+=M; forb[d]=1; forb[(M-d)%M]=1; }
        int k=0; for(int a=1;a<M;a++) if(!forb[a]) verts[k++]=a;
        nw=(k+63)/64; if(nw==0) nw=1;
        for(int i=0;i<k;i++){ memset(&adj[i],0,sizeof(bs)); for(int j=0;j<k;j++){ if(i==j) continue; int d=((verts[i]-verts[j])%M+M)%M; if(!forb[d]) adj[i].w[j>>6]|=1ULL<<(j&63);} }
        bs P; memset(&P,0,sizeof P); for(int i=0;i<k;i++) P.w[i>>6]|=1ULL<<(i&63);
        best=0; enumtarget=0; nodes=0; expand(P,0);
        int omega=best+1;
        if(mode[0]=='w'){ printf("%d\n",omega); }
        else if(mode[0]=='c'){ printf("%d 0",omega); for(int i=0;i<best;i++) printf(" %d",verts[bestset[i]]); printf("\n"); }
        else if(mode[0]=='e'){ int cap=atoi(mode+1); if(cap<=0) cap=1000; enumtarget=best; enumcap=cap; enumcount=0;
            if(best==0){ printf("%d 1\n0\n",omega); }
            else { // enumerate all cliques of size best (plus vertex 0) into a buffer via temporary file-less approach: print after header
                // simple: count first with a dry run? Instead print header with placeholder by buffering lines in memory.
                static char *buf=NULL; static size_t cap_b=0; FILE *mem=open_memstream(&buf,&cap_b); out=mem; nodes=0; expand(P,0); fclose(mem); out=stdout;
                printf("%d %d\n%s",omega,enumcount,buf); free(buf); buf=NULL; }
        }
        fflush(stdout);
    }
    return 0;
}
