/*
 * swap_compare.c  --  fast in-memory reimplementation of the per-atom
 *                     "swapping and comparing engine" from
 *                     compare_to_best_runme.com (the `foreach AB ( BA AB )`
 *                     loop, lines 154-247).
 *
 * The original spawns ~12,000 processes per run (flip_to_target_runme.com +
 * rmsd + awk, once per atom per arrangement), each re-parsing whole PDB files.
 * This does the identical computation entirely in memory.
 *
 * Usage:   swap_compare good.pdb BA AB
 *   good.pdb        reference structure
 *   BA, AB ...      arrangement names; for each, reads <name>.pdb and writes
 *                   sortme_<name>.txt, min_rmsds_<name>.txt, needed_swaps_<name>.txt
 *                   and APPENDS "<count> <name>" to number_of_needed_swaps.txt.
 *   min_rmsds.txt is left as the LAST arrangement's copy (matches the awk
 *   loop's leftover that the tcsh tail consumes). The no-suffix sortme.txt /
 *   needed_swaps.txt are produced by the tcsh selection step, not here.
 *
 * Faithful to the awk engine (numerically equivalent):
 *   - rmsd pairs atoms by name+altloc(+resname+chain+resnum); every atom here
 *     is dual-conformer, so each stat = sqrt((dA^2+dB^2)/2). Stats printed %.8f.
 *   - flip_to_target makes an INDEPENDENT decision per altloc over the
 *     [CON][DE] side-chain atoms (using distances rounded to %.4f, as awk does),
 *     adopting the name-swapped labelling iff its summed-squared-distance to the
 *     target is <= the current one (ties favour the flip). The reported
 *     bflip/sflip is the OR across altlocs.
 *   - score = (nstat-sstat)^2 + dstat^2 + nstat^2, printed %.6g (awk OFMT).
 *   - needed_swaps keeps rows with chain != 'S' and sstat < nstat, sorted by
 *     score descending; sorts use numeric-then-strcmp (GNU sort fallback).
 *
 * Build:  gcc -O2 -Wall -o swap_compare swap_compare.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#define MAXATOM 8000
#define MAXRES  1500
#define MAXRIDX 64      /* atoms per residue (heavy, both altlocs) */
#define NAMELEN 8

typedef struct {
    char   name[NAMELEN];   /* stripped atom name, e.g. "OD1" */
    char   altloc;          /* col 17 */
    char   resname[4];      /* cols 18-20, stripped */
    char   chain;           /* col 22 */
    int    resnum;          /* int(cols 23-27) */
    double x, y, z;
} Atom;

typedef struct {
    char chain;
    int  resnum;
    char resname[4];
    int  idx[MAXRIDX];
    int  n;
} Res;

/* ---- per-atomcode result record ---- */
typedef struct {
    char   atomcode[24];        /* name_<chain><resnum> */
    char   chain;
    double sstat, nstat, dstat;
    int    sflip, bflip;
    char   ss[24], ns[24], ds[24];   /* %.8f strings, reused everywhere */
} Rec;

/* ===================== small string helpers ===================== */

/* copy a fixed-column slice [start0, start0+len) of a line of length llen */
static void slice(const char *line, int llen, int start0, int len, char *out)
{
    int i;
    for (i = 0; i < len; i++) {
        int p = start0 + i;
        out[i] = (p < llen) ? line[p] : ' ';
    }
    out[len] = '\0';
}

/* bounded copy (manual loop; avoids snprintf truncation diagnostics) */
static void setstr(char *dst, int dstsz, const char *src)
{
    int i = 0;
    for (; i < dstsz - 1 && src[i]; i++) dst[i] = src[i];
    dst[i] = '\0';
}

/* strip leading and trailing spaces in place */
static void strip(char *s)
{
    int i = 0, j;
    while (s[i] == ' ') i++;
    if (i) memmove(s, s + i, strlen(s + i) + 1);
    j = (int)strlen(s);
    while (j > 0 && s[j - 1] == ' ') s[--j] = '\0';
}

/* awk-style "%.4f then back to number" rounding used by flip_to_target SSDs */
static double round4(double d)
{
    char b[32];
    snprintf(b, sizeof b, "%.4f", d);
    return atof(b);
}

/* ===================== PDB parsing ===================== */

/* returns the last whitespace-delimited token of a (NUL-terminated) line */
static void last_field(const char *line, char *out, int outsz)
{
    int len = (int)strlen(line);
    int i = len - 1, s, e;
    out[0] = '\0';
    while (i >= 0 && isspace((unsigned char)line[i])) i--;
    e = i;
    while (i >= 0 && !isspace((unsigned char)line[i])) i--;
    s = i + 1;
    if (e >= s && (e - s + 1) < outsz - 1) {
        memcpy(out, line + s, e - s + 1);
        out[e - s + 1] = '\0';
    }
}

/* parse a PDB file into out[]; skips hydrogens ($NF=="H"); returns atom count */
static int parse_pdb(const char *path, Atom *out)
{
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "swap_compare: cannot open %s\n", path); exit(1); }

    char line[512], buf[32], last[16];
    int n = 0;

    while (fgets(line, sizeof line, f)) {
        if (strncmp(line, "ATOM", 4) != 0 && strncmp(line, "HETATM", 6) != 0)
            continue;
        int len = (int)strlen(line);
        while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r'))
            line[--len] = '\0';

        last_field(line, last, sizeof last);
        if (strcmp(last, "H") == 0) continue;     /* drop hydrogens */

        if (n >= MAXATOM) { fprintf(stderr, "swap_compare: too many atoms\n"); exit(1); }
        Atom *a = &out[n];
        memset(a, 0, sizeof *a);

        slice(line, len, 11, 5, buf); strip(buf);                 /* name  cols 12-16 */
        setstr(a->name, sizeof a->name, buf);
        a->altloc = (len > 16) ? line[16] : ' ';                  /* altloc col 17    */
        slice(line, len, 17, 3, buf); strip(buf);                 /* resname cols18-20*/
        setstr(a->resname, sizeof a->resname, buf);
        a->chain = (len > 21) ? line[21] : ' ';                   /* chain col 22     */
        slice(line, len, 22, 5, buf); a->resnum = atoi(buf);      /* resnum cols23-27 */
        slice(line, len, 30, 8, buf); a->x = atof(buf);           /* x cols 31-38     */
        slice(line, len, 38, 8, buf); a->y = atof(buf);           /* y cols 39-46     */
        slice(line, len, 46, 8, buf); a->z = atof(buf);           /* z cols 47-54     */
        n++;
    }
    fclose(f);
    return n;
}

/* group atoms into residues by (chain,resnum), first-seen order */
static int build_res(const Atom *atoms, int na, Res *res)
{
    int nr = 0, i, j;
    for (i = 0; i < na; i++) {
        int r = -1;
        for (j = 0; j < nr; j++)
            if (res[j].chain == atoms[i].chain && res[j].resnum == atoms[i].resnum) { r = j; break; }
        if (r < 0) {
            if (nr >= MAXRES) { fprintf(stderr, "swap_compare: too many residues\n"); exit(1); }
            r = nr++;
            res[r].chain = atoms[i].chain;
            res[r].resnum = atoms[i].resnum;
            snprintf(res[r].resname, sizeof res[r].resname, "%s", atoms[i].resname);
            res[r].n = 0;
        }
        if (res[r].n >= MAXRIDX) { fprintf(stderr, "swap_compare: residue too large\n"); exit(1); }
        res[r].idx[res[r].n++] = i;
    }
    return nr;
}

/* ===================== flip logic ===================== */

static int is_flippable(const char *r)
{
    return !strcmp(r,"HIS") || !strcmp(r,"ASP") || !strcmp(r,"ASN") ||
           !strcmp(r,"GLU") || !strcmp(r,"GLN") || !strcmp(r,"TYR") ||
           !strcmp(r,"PHE") || !strcmp(r,"ARG");   /* ARG added 2026-10-22 (JH) */
}

/* atom name matches /^[CON][DE]/ on the stripped name, or ARG's NH1/NH2
 * (2026-10-22: the NH1/NH2 pair is X-ray-equivalent too; only reached for
 * flippable residues) */
static int is_conde(const char *n)
{
    return ((n[0]=='C'||n[0]=='O'||n[0]=='N') && (n[1]=='D'||n[1]=='E'))
        || !strcmp(n,"NH1") || !strcmp(n,"NH2");
}

/* return the equivalent-atom partner name for a flippable residue, or the
 * same name if the atom has no symmetric partner */
static const char *swapname(const char *res, const char *n)
{
    if (!strcmp(res,"HIS")) {
        if (!strcmp(n,"ND1")) return "CD2";
        if (!strcmp(n,"CD2")) return "ND1";
        if (!strcmp(n,"CE1")) return "NE2";
        if (!strcmp(n,"NE2")) return "CE1";
    } else if (!strcmp(res,"ASP")) {
        if (!strcmp(n,"OD1")) return "OD2";
        if (!strcmp(n,"OD2")) return "OD1";
    } else if (!strcmp(res,"ASN")) {
        if (!strcmp(n,"OD1")) return "ND2";
        if (!strcmp(n,"ND2")) return "OD1";
    } else if (!strcmp(res,"GLU")) {
        if (!strcmp(n,"OE1")) return "OE2";
        if (!strcmp(n,"OE2")) return "OE1";
    } else if (!strcmp(res,"GLN")) {
        if (!strcmp(n,"OE1")) return "NE2";
        if (!strcmp(n,"NE2")) return "OE1";
    } else if (!strcmp(res,"TYR") || !strcmp(res,"PHE")) {
        if (!strcmp(n,"CD1")) return "CD2";
        if (!strcmp(n,"CD2")) return "CD1";
        if (!strcmp(n,"CE1")) return "CE2";
        if (!strcmp(n,"CE2")) return "CE1";
    } else if (!strcmp(res,"ARG")) {          /* added 2026-10-22 (JH) */
        if (!strcmp(n,"NH1")) return "NH2";
        if (!strcmp(n,"NH2")) return "NH1";
    }
    return n;
}

static double dist(const Atom *a, const Atom *b)
{
    double dx = a->x - b->x, dy = a->y - b->y, dz = a->z - b->z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

/* find atom with given name+altloc in a residue atom list; NULL if absent */
static const Atom *find_atom(const Atom *atoms, const int *idx, int n,
                             const char *name, char altloc)
{
    int i;
    for (i = 0; i < n; i++) {
        const Atom *a = &atoms[idx[i]];
        if (a->altloc == altloc && !strcmp(a->name, name)) return a;
    }
    return NULL;
}

/*
 * flip a residue (working copy `wa`, indices 0..wn-1) toward the target good
 * residue (gatoms/gidx/gn). Decides independently per altloc; renames the
 * symmetric atoms where the swapped labelling fits the target at least as well.
 * Returns 1 if ANY altloc was flipped (the bflip/sflip boolean).
 */
static int flip_residue(Atom *wa, int wn, const char *resname,
                        const Atom *gatoms, const int *gidx, int gn)
{
    if (!is_flippable(resname) || gn == 0) return 0;

    char altlocs[8];
    int  nalt = 0, i, j;
    for (i = 0; i < wn; i++) {                 /* distinct altlocs in this residue */
        char L = wa[i].altloc;
        int seen = 0;
        for (j = 0; j < nalt; j++) if (altlocs[j] == L) { seen = 1; break; }
        if (!seen && nalt < 8) altlocs[nalt++] = L;
    }

    int any = 0;
    for (j = 0; j < nalt; j++) {
        char L = altlocs[j];
        double cur = 0, swp = 0;
        int ccur = 0, cswp = 0;

        for (i = 0; i < wn; i++) {
            if (wa[i].altloc != L || !is_conde(wa[i].name)) continue;
            const Atom *g  = find_atom(gatoms, gidx, gn, wa[i].name, L);
            if (g)  { double d = round4(dist(&wa[i], g));  cur += d*d; ccur++; }
            const char *sn = swapname(resname, wa[i].name);
            const Atom *g2 = find_atom(gatoms, gidx, gn, sn, L);
            if (g2) { double d = round4(dist(&wa[i], g2)); swp += d*d; cswp++; }
        }

        /* adopt the flip iff a flipped pairing exists and is no worse
         * (tie -> flip, matching GNU `sort -k4g` where "flip" < "noflip") */
        if (cswp > 0 && (ccur == 0 || swp <= cur)) {
            char newn[MAXRIDX][NAMELEN];
            for (i = 0; i < wn; i++) {
                const char *nn = (wa[i].altloc == L && is_conde(wa[i].name))
                               ? swapname(resname, wa[i].name) : wa[i].name;
                snprintf(newn[i], NAMELEN, "%s", nn);
            }
            for (i = 0; i < wn; i++) snprintf(wa[i].name, NAMELEN, "%s", newn[i]);
            any = 1;
        }
    }
    return any;
}

/* ===================== rmsd (pair-by-name-within-residue) ===================== */

/* RMS distance between atoms named `name` in two residue atom lists, pairing by
 * altloc (= rmsd's 18-char key, since chain/resnum/resname are equal here).
 * Returns sqrt(sum(d^2)/npairs); 0 if no pair (does not occur on matched data) */
static double pairdist(const Atom *A, const int *ai, int an,
                       const Atom *B, const int *bi, int bn,
                       const char *name)
{
    double sum = 0;
    int cnt = 0, i;
    for (i = 0; i < an; i++) {
        const Atom *a = &A[ai[i]];
        if (strcmp(a->name, name)) continue;
        const Atom *b = find_atom(B, bi, bn, name, a->altloc);
        if (!b) continue;
        double dx = a->x - b->x, dy = a->y - b->y, dz = a->z - b->z;
        sum += dx*dx + dy*dy + dz*dz;
        cnt++;
    }
    return cnt ? sqrt(sum / cnt) : 0.0;
}

/* ===================== sort comparators ===================== */

typedef struct { double val; char line[64]; } MinRow;
typedef struct { double score; char cmp[128]; char out[96]; } NSRow;

static int cmp_min(const void *pa, const void *pb)         /* ascending by val */
{
    const MinRow *a = pa, *b = pb;
    if (a->val < b->val) return -1;
    if (a->val > b->val) return 1;
    return strcmp(a->line, b->line);
}

static int cmp_ns(const void *pa, const void *pb)          /* descending: GNU sort -gr */
{
    const NSRow *a = pa, *b = pb;
    if (a->score < b->score) return 1;
    if (a->score > b->score) return -1;
    return strcmp(b->cmp, a->cmp);
}

/* ===================== globals ===================== */

static Atom good[MAXATOM];   static int good_n;
static Res  gres[MAXRES];    static int gres_n;
static Atom model[MAXATOM];  static int model_n;
static Res  mres[MAXRES];    static int mres_n;
static Rec  recs[MAXATOM];   static int nrec;

/* locate the good residue matching (chain,resnum); -1 if none */
static int find_gres(char chain, int resnum)
{
    int i;
    for (i = 0; i < gres_n; i++)
        if (gres[i].chain == chain && gres[i].resnum == resnum) return i;
    return -1;
}

/* ===================== per-arrangement engine ===================== */

static void process_arrangement(const char *name, int is_last)
{
    char path[256];
    snprintf(path, sizeof path, "%s.pdb", name);
    model_n = parse_pdb(path, model);
    mres_n  = build_res(model, model_n, mres);

    nrec = 0;

    int r, i, k;
    for (r = 0; r < mres_n; r++) {
        Res *mr = &mres[r];

        /* FM = flipped model residue; FS = altloc-swapped then flipped */
        Atom fm[MAXRIDX], fs[MAXRIDX];
        int  fmidx[MAXRIDX], fsidx[MAXRIDX];
        for (i = 0; i < mr->n; i++) {
            fm[i] = model[mr->idx[i]];
            fs[i] = model[mr->idx[i]];
            if      (fs[i].altloc == 'A') fs[i].altloc = 'B';
            else if (fs[i].altloc == 'B') fs[i].altloc = 'A';
            fmidx[i] = i; fsidx[i] = i;
        }

        int gi = find_gres(mr->chain, mr->resnum);
        const Atom *ga = (gi >= 0) ? good : NULL;
        const int  *gidx = (gi >= 0) ? gres[gi].idx : NULL;
        int gn = (gi >= 0) ? gres[gi].n : 0;

        int bflip = flip_residue(fm, mr->n, mr->resname, ga ? ga : good, gidx ? gidx : fmidx, gn);
        int sflip = flip_residue(fs, mr->n, mr->resname, ga ? ga : good, gidx ? gidx : fsidx, gn);

        /* one record per unique atom name (post-flip == original name set) */
        for (i = 0; i < mr->n; i++) {
            const char *nm = fm[i].name;
            int dup = 0;
            for (k = 0; k < i; k++) if (!strcmp(fm[k].name, nm)) { dup = 1; break; }
            if (dup) continue;

            double nstat = pairdist(fm, fmidx, mr->n, ga, gidx, gn, nm);
            double sstat = pairdist(fs, fsidx, mr->n, ga, gidx, gn, nm);
            double dstat = pairdist(fm, fmidx, mr->n, fs, fsidx, mr->n, nm);

            Rec *rc = &recs[nrec++];
            snprintf(rc->atomcode, sizeof rc->atomcode, "%s_%c%d", nm, mr->chain, mr->resnum);
            rc->chain = mr->chain;
            snprintf(rc->ss, sizeof rc->ss, "%.8f", sstat);
            snprintf(rc->ns, sizeof rc->ns, "%.8f", nstat);
            snprintf(rc->ds, sizeof rc->ds, "%.8f", dstat);
            /* downstream awk/sort see only these %.8f strings, so do all later
             * comparisons/scoring on the rounded values (matches the original) */
            rc->sstat = atof(rc->ss);
            rc->nstat = atof(rc->ns);
            rc->dstat = atof(rc->ds);
            rc->sflip = sflip; rc->bflip = bflip;
        }
    }

    /* ---- sortme_<name>.txt ---- */
    snprintf(path, sizeof path, "sortme_%s.txt", name);
    FILE *f = fopen(path, "w");
    for (i = 0; i < nrec; i++)
        fprintf(f, "%s %s %s %s %d %d\n", recs[i].atomcode,
                recs[i].ss, recs[i].ns, recs[i].ds, recs[i].sflip, recs[i].bflip);
    fclose(f);

    /* ---- min_rmsds_<name>.txt (and min_rmsds.txt for the last arrangement) ---- */
    static MinRow mrow[MAXATOM];
    int nm = 0;
    for (i = 0; i < nrec; i++) {
        Rec *rc = &recs[i];
        char  conf; const char *vstr; double v; int flip;
        if (rc->sstat < rc->nstat) { conf='B'; vstr=rc->ss; v=rc->sstat; flip=rc->sflip; }
        else                       { conf='A'; vstr=rc->ns; v=rc->nstat; flip=rc->bflip; }
        mrow[nm].val = v;
        snprintf(mrow[nm].line, sizeof mrow[nm].line, "%s %s %c %d",
                 rc->atomcode, vstr, conf, flip);
        nm++;
    }
    qsort(mrow, nm, sizeof mrow[0], cmp_min);
    snprintf(path, sizeof path, "min_rmsds_%s.txt", name);
    f = fopen(path, "w");
    for (i = 0; i < nm; i++) fprintf(f, "%s\n", mrow[i].line);
    fclose(f);
    if (is_last) {
        f = fopen("min_rmsds.txt", "w");
        for (i = 0; i < nm; i++) fprintf(f, "%s\n", mrow[i].line);
        fclose(f);
    }

    /* ---- needed_swaps_<name>.txt ---- */
    static NSRow ns[MAXATOM];
    int nn = 0;
    for (i = 0; i < nrec; i++) {
        Rec *rc = &recs[i];
        if (rc->chain == 'S') continue;          /* grep -v _S (water chain) */
        if (!(rc->sstat < rc->nstat)) continue;
        double score = (rc->nstat - rc->sstat) * (rc->nstat - rc->sstat)
                     + rc->dstat * rc->dstat
                     + rc->nstat * rc->nstat;
        char sc[32];
        snprintf(sc, sizeof sc, "%.6g", score);   /* awk OFMT */
        ns[nn].score = score;
        /* tie-break key = the pre-reformat intermediate line `sort -gr` saw */
        snprintf(ns[nn].cmp, sizeof ns[nn].cmp, "%s %s %s %s %s %d %d",
                 sc, rc->atomcode, rc->ss, rc->ns, rc->ds, rc->sflip, rc->bflip);
        snprintf(ns[nn].out, sizeof ns[nn].out, "SWAP %s %s %s %s %s",
                 rc->atomcode, sc, rc->ss, rc->ns, rc->ds);
        nn++;
    }
    qsort(ns, nn, sizeof ns[0], cmp_ns);
    snprintf(path, sizeof path, "needed_swaps_%s.txt", name);
    f = fopen(path, "w");
    for (i = 0; i < nn; i++) fprintf(f, "%s\n", ns[i].out);
    fclose(f);

    /* ---- number_of_needed_swaps.txt (append, like the original tee -a) ---- */
    f = fopen("number_of_needed_swaps.txt", "a");
    fprintf(f, "%d %s\n", nn, name);
    fclose(f);
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        fprintf(stderr, "usage: %s good.pdb arrangement [arrangement ...]\n", argv[0]);
        return 9;
    }
    good_n  = parse_pdb(argv[1], good);
    gres_n  = build_res(good, good_n, gres);

    int k;
    for (k = 2; k < argc; k++)
        process_arrangement(argv[k], k == argc - 1);

    return 0;
}
