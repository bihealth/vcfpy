#define HTS_MIN_MARKER_DIST 0x10000

// Finds the special meta bin
//  ((1<<(3 * n_lvls + 3)) - 1) / 7 + 1
#define META_BIN(idx) ((idx)->n_bins + 1)

#define pair64_lt(a,b) ((a).u < (b).u)
#define pair64max_lt(a,b) ((a).u < (b).u || \
                           ((a).u == (b).u && (a).max < (b).max))

KSORT_INIT_STATIC(_off, hts_pair64_t, pair64_lt)
KSORT_INIT_STATIC(_off_max, hts_pair64_max_t, pair64max_lt)

typedef struct {
    int32_t m, n;
    uint64_t loff;
    hts_pair64_t *list;
} bins_t;

KHASH_MAP_INIT_INT(bin, bins_t)
typedef khash_t(bin) bidx_t;

typedef struct {
    hts_pos_t n, m;
    uint64_t *offset;
} lidx_t;

struct hts_idx_t {
    int fmt, min_shift, n_lvls, n_bins;
    uint32_t l_meta;
    int32_t n, m;
    uint64_t n_no_coor;
    bidx_t **bidx;
    lidx_t *lidx;
    uint8_t *meta; // MUST have a terminating NUL on the end
    int tbi_n, last_tbi_tid;
    struct {
        uint32_t last_bin, save_bin;
        hts_pos_t last_coor;
        int last_tid, save_tid, finished;
        uint64_t last_off, save_off;
        uint64_t off_beg, off_end;
        uint64_t n_mapped, n_unmapped;
    } z; // keep internal states
};

static char * idx_format_name(int fmt) {
    switch (fmt) {
        case HTS_FMT_CSI: return "csi";
        case HTS_FMT_BAI: return "bai";
        case HTS_FMT_TBI: return "tbi";
        case HTS_FMT_CRAI: return "crai";
        default: return "unknown";
    }
}

/*  tbx.c -- tabix API functions.

    Copyright (C) 2009, 2010, 2012-2015, 2017-2020, 2022-2023, 2025 Genome Research Ltd.
    Copyright (C) 2010-2012 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#define HTS_BUILDING_LIBRARY // Enables HTSLIB_EXPORT, see htslib/hts_defs.h
#include <config.h>

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <errno.h>
#include "htslib/tbx.h"
#include "htslib/bgzf.h"
#include "htslib/hts_endian.h"
#include "hts_internal.h"

#include "htslib/khash.h"
KHASH_DECLARE(s2i, kh_cstr_t, int64_t)

HTSLIB_EXPORT
const tbx_conf_t tbx_conf_gff = { 0, 1, 4, 5, '#', 0 };

HTSLIB_EXPORT
const tbx_conf_t tbx_conf_bed = { TBX_UCSC, 1, 2, 3, '#', 0 };

HTSLIB_EXPORT
const tbx_conf_t tbx_conf_psltbl = { TBX_UCSC, 15, 17, 18, '#', 0 };

HTSLIB_EXPORT
const tbx_conf_t tbx_conf_sam = { TBX_SAM, 3, 4, 0, '@', 0 };

HTSLIB_EXPORT
const tbx_conf_t tbx_conf_vcf = { TBX_VCF, 1, 2, 0, '#', 0 };
const tbx_conf_t tbx_conf_gaf = { TBX_GAF, 1, 6, 0, '#', 0 };

typedef struct {
    int64_t beg, end;
    char *ss, *se;
    int tid;
} tbx_intv_t;

static inline int get_tid(tbx_t *tbx, const char *ss, int is_add)
{
    khint_t k;
    khash_t(s2i) *d;
    if ((tbx->conf.preset&0xffff) == TBX_GAF) return(0);
    if (tbx->dict == 0) tbx->dict = kh_init(s2i);
    if (!tbx->dict) return -1; // Out of memory
    d = (khash_t(s2i)*)tbx->dict;
    if (is_add) {
        int absent;
        k = kh_put(s2i, d, ss, &absent);
        if (absent < 0) {
            return -1; // Out of memory
        } else if (absent) {
            char *ss_dup = strdup(ss);
            if (ss_dup) {
                kh_key(d, k) = ss_dup;
                kh_val(d, k) = kh_size(d) - 1;
            } else {
                kh_del(s2i, d, k);
                return -1; // Out of memory
            }
        }
    } else k = kh_get(s2i, d, ss);
    return k == kh_end(d)? -1 : kh_val(d, k);
}

int tbx_name2id(tbx_t *tbx, const char *ss)
{
    return get_tid(tbx, ss, 0);
}

int tbx_parse1(const tbx_conf_t *conf, size_t len, char *line, tbx_intv_t *intv)
{
    size_t i, b = 0;
    int id = 1, getlen = 0, alcnt = 0, use_svlen = 0, lenpos = -1;
    char *s, *t;
    uint8_t svlenals[8192];
    int64_t reflen = 0, svlen = 0, fmtlen = 0, tmp = 0;

    intv->ss = intv->se = 0; intv->beg = intv->end = -1;
    for (i = 0; i <= len; ++i) {
        if (line[i] == '\t' || line[i] == 0) {
            if (id == conf->sc) {
                intv->ss = line + b; intv->se = line + i;
            } else if (id == conf->bc) {
                // here ->beg is 0-based.
                if ((conf->preset&0xffff) == TBX_GAF){
                    // if gaf find the smallest and largest node id
                    char *t;
                    int64_t nodeid = -1;
                    for (s = line + b + 1; s < line + i;) {
                        nodeid = strtoll(s, &t, 0);
                        if(intv->beg == -1){
                            intv->beg = intv->end = nodeid;
                        } else {
                            if(nodeid < intv->beg){
                                intv->beg = nodeid;
                            }

                            if(nodeid > intv->end){
                                intv->end = nodeid;
                            }
                        }
                        s = t + 1;
                    }
                } else {
                    intv->beg = strtoll(line + b, &s, 0);

                    if (conf->bc <= conf->ec) // don't overwrite an already set end point
                        intv->end = intv->beg;

                    if ( s==line+b ) return -1; // expected int

                    if (!(conf->preset&TBX_UCSC))
                        --intv->beg;
                    else if (conf->bc <= conf->ec)
                        ++intv->end;

                    if (intv->beg < 0) {
                        hts_log_warning("Coordinate <= 0 detected. "
                                        "Did you forget to use the -0 option?");
                        intv->beg = 0;
                    }
                    if (intv->end < 1) intv->end = 1;
                }
            } else {
                if ((conf->preset&0xffff) == TBX_GENERIC) {
                    if (id == conf->ec)
                    {
                        intv->end = strtoll(line + b, &s, 0);
                        if ( s==line+b ) return -1; // expected int
                    }
                } else if ((conf->preset&0xffff) == TBX_SAM) {
                    if (id == 6) { // CIGAR
                        int l = 0;
                        char *t;
                        for (s = line + b; s < line + i;) {
                            long x = strtol(s, &t, 10);
                            char op = toupper_c(*t);
                            if (op == 'M' || op == 'D' || op == 'N') l += x;
                            s = t + 1;
                        }
                        if (l == 0) l = 1;
                        intv->end = intv->beg + l;
                    }
                } else if ((conf->preset&0xffff) == TBX_VCF) {
                    if (id == 4) { //ref allele
                        if (b < i) intv->end = intv->beg + (i - b);
                        ++alcnt;
                        reflen = i - b;
                    } if (id == 5) {    //alt allele
                        int lastbyte = 0, c = line[i];
                        svlenals[lastbyte] = 0;
                        line[i] = 0;
                        s = line + b;
                        do {
                            t = strchr(s, ',');
                            if (alcnt >> 3 != lastbyte) {   //initialize insals
                                lastbyte = alcnt >> 3;
                                svlenals[lastbyte] = 0;
                            }
                            ++alcnt;
                            if (t) {
                                *t = 0;
                            }
                            if (svlen_on_ref_for_vcf_alt(s, -1)) {
                                // Need to check SVLEN for this ALT
                                svlenals[lastbyte] |= 1 << ((alcnt - 1) & 7);
                                use_svlen = 1;
                            } else if (!strcmp("<*>", s) ||
                                       !strcmp("<NON_REF>", s)) {  //note gvcf
                                getlen = 1;
                            }
                            if (t) {
                                *t = ',';
                                s = t + 1;
                            }
                        } while (t && alcnt < 65536);   //max allcnt is 65535
                        line[i] = c;
                    } else if (id == 8) { //INFO, look for "END=" / "SVLEN"
                        int c = line[i], d = 1;
                        line[i] = 0;
                        s = strstr(line + b, "END=");
                        if (s == line + b) s += 4;
                        else if (s) {
                            s = strstr(line + b, ";END=");
                            if (s) s += 5;
                        }
                        if (s && *s != '.') {
                            long long end = strtoll(s, &s, 0);
                            if (end <= intv->beg) {
                                static int reported = 0;
                                if (!reported) {
                                    int l = intv->ss ? (int) (intv->se - intv->ss) : 0;
                                    hts_log_warning("VCF INFO/END=%lld is smaller than POS at %.*s:%"PRIhts_pos"\n"
                                                    "This tag will be ignored. "
                                                    "Note: only one invalid END tag will be reported.",
                                                    end, l >= 0 ? l : 0,
                                                    intv->ss ? intv->ss : "",
                                                    intv->beg);
                                    reported = 1;
                                }
                            } else {
                                intv->end = end;
                            }
                        }
                        s = strstr(line + b, "SVLEN=");
                        if (s == line + b) s += 6;  //at start of info
                        else if (s) {               //not at the start
                            s = strstr(line + b, ";SVLEN=");
                            if (s) s += 7;
                        }
                        while (s && d < alcnt) {
                            t = strchr(s, ',');
                            if ((use_svlen) && (svlenals[d >> 3] & (1 << (d & 7)))) {
                                // <DEL> symbolic allele
                                tmp = atoll(s);
                                tmp = tmp < 0 ? llabs(tmp) : tmp;
                            } else {
                                tmp = 1;
                            }
                            svlen = svlen < tmp ? tmp : svlen;
                            s = t ? t + 1 : NULL;
                            ++d;
                        }
                        line[i] = c;
                    } else if (getlen && id == 9 ) {    //FORMAT
                        int c = line[i], pos = -1;
                        line[i] = 0;
                        s = line + b;
                        while (s) {
                            ++pos;
                            if (!(t = strchr(s, ':'))) {    //no further fields
                                if (!strcmp(s, "LEN")) {
                                    lenpos = pos;
                                }
                                break;  //not present at all!
                            } else {
                                *t = '\0';
                                if (!strcmp(s, "LEN")) {
                                    lenpos = pos;
                                    *t = ':';
                                    break;
                                }
                                *t = ':';
                                s = t + 1;  //check next one
                            }
                        }
                        line[i] = c;
                        if (lenpos == -1) { //not present
                            break;
                        }
                    } else if (id > 9 && getlen && lenpos != -1) {
                        //get LEN from sample
                        int c = line[i], d = 0;
                        line[i] = 0; tmp = 0;
                        s = line + b;
                        for (d = 0; d <= lenpos; ++d) {
                            if (d == lenpos) {
                                tmp = atoll(s);
                                break;
                            }
                            if ((t = strchr(s, ':'))) {
                                s = t + 1;
                            } else {
                                break;    //not in sycn with fmt def!
                            }
                        }
                        fmtlen = fmtlen < tmp ? tmp : fmtlen;
                        line[i] = c;
                    }
                }
            }
            b = i + 1;  //beginning if current field
            ++id;
        }
    }
    if ((conf->preset&0xffff) == TBX_VCF) {
        tmp = reflen < svlen ?
                svlen < fmtlen ? fmtlen : svlen :
                reflen < fmtlen ? fmtlen : reflen ;
        tmp += intv->beg;
        intv->end = intv->end < tmp ? tmp : intv->end;

        //NOTE: 'end' calculation be in sync with end/rlen in vcf.c:get_rlen
    }
    if (intv->ss == 0 || intv->se == 0 || intv->beg < 0 || intv->end < 0) return -1;
    return 0;
}

static inline int get_intv(tbx_t *tbx, kstring_t *str, tbx_intv_t *intv, int is_add)
{
    if (tbx_parse1(&tbx->conf, str->l, str->s, intv) == 0) {
        int c = *intv->se;
        *intv->se = '\0';
        if ((tbx->conf.preset&0xffff) == TBX_GAF){
            intv->tid = 0;
        } else {
            intv->tid = get_tid(tbx, intv->ss, is_add);
        }
        *intv->se = c;
        if (intv->tid < 0) return -2;  // get_tid out of memory
        return (intv->beg >= 0 && intv->end >= 0)? 0 : -1;
    } else {
        char *type = NULL;
        switch (tbx->conf.preset&0xffff)
        {
            case TBX_SAM: type = "TBX_SAM"; break;
            case TBX_VCF: type = "TBX_VCF"; break;
            case TBX_GAF: type = "TBX_GAF"; break;
            case TBX_UCSC: type = "TBX_UCSC"; break;
            default: type = "TBX_GENERIC"; break;
        }
        if (hts_is_utf16_text(str))
            hts_log_error("Failed to parse %s: offending line appears to be encoded as UTF-16", type);
        else
            hts_log_error("Failed to parse %s: was wrong -p [type] used?\nThe offending line was: \"%s\"",
                type, str->s);
        return -1;
    }
}

/*
 * Called by tabix iterator to read the next record.
 * Returns    >=  0 on success
 *               -1 on EOF
 *            <= -2 on error
 */
int tbx_readrec(BGZF *fp, void *tbxv, void *sv, int *tid, hts_pos_t *beg, hts_pos_t *end)
{
    tbx_t *tbx = (tbx_t *) tbxv;
    kstring_t *s = (kstring_t *) sv;
    int ret;

    // Get a line until either EOF or a non-meta character
    do {
        ret = bgzf_getline(fp, '\n', s);
    } while (ret >= 0 && s->l && *s->s == tbx->conf.meta_char);

    // Parse line
    if (ret >= 0)  {
        tbx_intv_t intv;
        if (get_intv(tbx, s, &intv, 0) < 0)
            return -2;
        *tid = intv.tid; *beg = intv.beg; *end = intv.end;
    }

    return ret;
}

static int tbx_set_meta(tbx_t *tbx)
{
    int i, l = 0, l_nm;
    uint32_t x[7];
    char **name;
    uint8_t *meta;
    khint_t k;
    khash_t(s2i) *d = (khash_t(s2i)*)tbx->dict;

    memcpy(x, &tbx->conf, 24);
    name = (char**)malloc(sizeof(char*) * kh_size(d));
    if (!name) return -1;
    for (k = kh_begin(d), l = 0; k != kh_end(d); ++k) {
        if (!kh_exist(d, k)) continue;
        name[kh_val(d, k)] = (char*)kh_key(d, k);
        l += strlen(kh_key(d, k)) + 1; // +1 to include '\0'
    }
    l_nm = x[6] = l;
    meta = (uint8_t*)malloc(l_nm + 28);
    if (!meta) { free(name); return -1; }
    if (ed_is_big())
        for (i = 0; i < 7; ++i)
            x[i] = ed_swap_4(x[i]);
    memcpy(meta, x, 28);
    for (l = 28, i = 0; i < (int)kh_size(d); ++i) {
        int x = strlen(name[i]) + 1;
        memcpy(meta + l, name[i], x);
        l += x;
    }
    free(name);
    hts_idx_set_meta(tbx->idx, l, meta, 0);
    return 0;
}

// Minimal effort parser to extract reference length out of VCF header line
// This is used only used to adjust the number of levels if necessary,
// so not a major problem if it doesn't always work.
static void adjust_max_ref_len_vcf(const char *str, int64_t *max_ref_len)
{
    const char *ptr;
    int64_t len;
    if (strncmp(str, "##contig", 8) != 0) return;
    ptr = strstr(str + 8, "length");
    if (!ptr) return;
    for (ptr += 6; *ptr == ' ' || *ptr == '='; ptr++) {}
    len = strtoll(ptr, NULL, 10);
    if (*max_ref_len < len) *max_ref_len = len;
}

// Same for sam files
static void adjust_max_ref_len_sam(const char *str, int64_t *max_ref_len)
{
    const char *ptr;
    int64_t len;
    if (strncmp(str, "@SQ", 3) != 0) return;
    ptr = strstr(str + 3, "\tLN:");
    if (!ptr) return;
    ptr += 4;
    len = strtoll(ptr, NULL, 10);
    if (*max_ref_len < len) *max_ref_len = len;
}

// Adjusts number of levels if not big enough.  This can happen for
// files with very large contigs.
static int adjust_n_lvls(int min_shift, int n_lvls, int64_t max_len)
{
    int64_t s = hts_bin_maxpos(min_shift, n_lvls);
    max_len += 256;
    for (; max_len > s; ++n_lvls, s <<= 3) {}
    return n_lvls;
}

tbx_t *tbx_index(BGZF *fp, int min_shift, const tbx_conf_t *conf)
{
    tbx_t *tbx;
    kstring_t str;
    int ret, first = 0, n_lvls, fmt;
    int64_t lineno = 0;
    uint64_t last_off = 0;
    tbx_intv_t intv;
    int64_t max_ref_len = 0;

    str.s = 0; str.l = str.m = 0;
    tbx = (tbx_t*)calloc(1, sizeof(tbx_t));
    if (!tbx) return NULL;
    tbx->conf = *conf;
    if (min_shift > 0) n_lvls = (TBX_MAX_SHIFT - min_shift + 2) / 3, fmt = HTS_FMT_CSI;
    else min_shift = 14, n_lvls = 5, fmt = HTS_FMT_TBI;
    while ((ret = bgzf_getline(fp, '\n', &str)) >= 0) {
        ++lineno;
        if (str.s[0] == tbx->conf.meta_char && fmt == HTS_FMT_CSI) {
            switch (tbx->conf.preset) {
                case TBX_SAM:
                    adjust_max_ref_len_sam(str.s, &max_ref_len); break;
                case TBX_VCF:
                    adjust_max_ref_len_vcf(str.s, &max_ref_len); break;
                default:
                    break;
            }
        }
        if (lineno <= tbx->conf.line_skip || str.s[0] == tbx->conf.meta_char) {
            last_off = bgzf_tell(fp);
            continue;
        }
        if (first == 0) {
            if (fmt == HTS_FMT_CSI) {
                if (!max_ref_len)
                    max_ref_len = (int64_t)100*1024*1024*1024; // 100G default
                n_lvls = adjust_n_lvls(min_shift, n_lvls, max_ref_len);
            }
            tbx->idx = hts_idx_init(0, fmt, last_off, min_shift, n_lvls);
            if (!tbx->idx) goto fail;
            first = 1;
        }
        ret = get_intv(tbx, &str, &intv, 1);
        if (ret < 0) goto fail;  // Out of memory or unparsable lines
        if (hts_idx_push(tbx->idx, intv.tid, intv.beg, intv.end,
                         bgzf_tell(fp), 1) < 0) {
            goto fail;
        }
    }
    if (ret < -1) goto fail;
    if ( !tbx->idx ) tbx->idx = hts_idx_init(0, fmt, last_off, min_shift, n_lvls);   // empty file
    if (!tbx->idx) goto fail;
    if ( !tbx->dict ) tbx->dict = kh_init(s2i);
    if (!tbx->dict) goto fail;
    if (hts_idx_finish(tbx->idx, bgzf_tell(fp)) != 0) goto fail;
    if (tbx_set_meta(tbx) != 0) goto fail;
    free(str.s);
    return tbx;

 fail:
    free(str.s);
    tbx_destroy(tbx);
    return NULL;
}

void tbx_destroy(tbx_t *tbx)
{
    khash_t(s2i) *d = (khash_t(s2i)*)tbx->dict;
    if (d != NULL)
    {
        khint_t k;
        for (k = kh_begin(d); k != kh_end(d); ++k)
            if (kh_exist(d, k)) free((char*)kh_key(d, k));
    }
    hts_idx_destroy(tbx->idx);
    kh_destroy(s2i, d);
    free(tbx);
}

int tbx_index_build3(const char *fn, const char *fnidx, int min_shift, int n_threads, const tbx_conf_t *conf)
{
    tbx_t *tbx;
    BGZF *fp;
    int ret;
    if ((fp = bgzf_open(fn, "r")) == 0) return -1;
    if ( n_threads ) bgzf_mt(fp, n_threads, 256);
    if ( bgzf_compression(fp) != bgzf ) { bgzf_close(fp); return -2; }
    tbx = tbx_index(fp, min_shift, conf);
    bgzf_close(fp);
    if ( !tbx ) return -1;
    ret = hts_idx_save_as(tbx->idx, fn, fnidx, min_shift > 0? HTS_FMT_CSI : HTS_FMT_TBI);
    tbx_destroy(tbx);
    return ret;
}

int tbx_index_build2(const char *fn, const char *fnidx, int min_shift, const tbx_conf_t *conf)
{
    return tbx_index_build3(fn, fnidx, min_shift, 0, conf);
}

int tbx_index_build(const char *fn, int min_shift, const tbx_conf_t *conf)
{
    return tbx_index_build3(fn, NULL, min_shift, 0, conf);
}

static tbx_t *index_load(const char *fn, const char *fnidx, int flags)
{
    tbx_t *tbx;
    uint8_t *meta;
    char *nm, *p;
    uint32_t l_meta, l_nm;
    tbx = (tbx_t*)calloc(1, sizeof(tbx_t));
    if (!tbx)
        return NULL;
    tbx->idx = hts_idx_load3(fn, fnidx, HTS_FMT_TBI, flags);
    if ( !tbx->idx )
    {
        free(tbx);
        return NULL;
    }
    meta = hts_idx_get_meta(tbx->idx, &l_meta);
    if ( !meta || l_meta < 28) goto invalid;

    tbx->conf.preset = le_to_i32(&meta[0]);
    tbx->conf.sc = le_to_i32(&meta[4]);
    tbx->conf.bc = le_to_i32(&meta[8]);
    tbx->conf.ec = le_to_i32(&meta[12]);
    tbx->conf.meta_char = le_to_i32(&meta[16]);
    tbx->conf.line_skip = le_to_i32(&meta[20]);
    l_nm = le_to_u32(&meta[24]);
    if (l_nm > l_meta - 28) goto invalid;

    p = nm = (char*)meta + 28;
    // This assumes meta is NUL-terminated, so we can merrily strlen away.
    // hts_idx_load_local() assures this for us by adding a NUL on the end
    // of whatever it reads.
    for (; p - nm < l_nm; p += strlen(p) + 1) {
        if (get_tid(tbx, p, 1) < 0) {
            hts_log_error("%s", strerror(errno));
            goto fail;
        }
    }
    return tbx;

 invalid:
    hts_log_error("Invalid index header for %s", fnidx ? fnidx : fn);

 fail:
    tbx_destroy(tbx);
    return NULL;
}

tbx_t *tbx_index_load3(const char *fn, const char *fnidx, int flags)
{
    return index_load(fn, fnidx, flags);
}

tbx_t *tbx_index_load2(const char *fn, const char *fnidx)
{
    return index_load(fn, fnidx, 1);
}

tbx_t *tbx_index_load(const char *fn)
{
    return index_load(fn, NULL, 1);
}

const char **tbx_seqnames(tbx_t *tbx, int *n)
{
    khash_t(s2i) *d = (khash_t(s2i)*)tbx->dict;
    if (d == NULL)
    {
        *n = 0;
        return calloc(1, sizeof(char *));
    }
    int tid, m = kh_size(d);
    const char **names = (const char**) calloc(m,sizeof(const char*));
    khint_t k;
    if (!names) {
        *n = 0;
        return NULL;
    }
    for (k=kh_begin(d); k<kh_end(d); k++)
    {
        if ( !kh_exist(d,k) ) continue;
        tid = kh_val(d,k);
        assert( tid<m );
        names[tid] = kh_key(d,k);
    }
    // sanity check: there should be no gaps
    for (tid=0; tid<m; tid++)
        assert(names[tid]);
    *n = m;
    return names;
}



/*
 * A variant of hts_parse_reg which is reference-id aware.  It uses
 * the iterator name2id callbacks to validate the region tokenisation works.
 *
 * This is necessary due to GRCh38 HLA additions which have reference names
 * like "HLA-DRB1*12:17".
 *
 * All parameters are mandatory.
 *
 * To work around ambiguous parsing issues, eg both "chr1" and "chr1:100-200"
 * are reference names, we may quote using curly braces.
 * Thus "{chr1}:100-200" and "{chr1:100-200}" disambiguate the above example.
 *
 * Flags are used to control how parsing works, and can be one of the below.
 *
 * HTS_PARSE_LIST:
 *     If present, the region is assmed to be a comma separated list and
 *     position parsing will not contain commas (this implicitly
 *     clears HTS_PARSE_THOUSANDS_SEP in the call to hts_parse_decimal).
 *     On success the return pointer will be the start of the next region, ie
 *     the character after the comma.  (If *ret != '\0' then the caller can
 *     assume another region is present in the list.)
 *
 *     If not set then positions may contain commas.  In this case the return
 *     value should point to the end of the string, or NULL on failure.
 *
 * HTS_PARSE_ONE_COORD:
 *     If present, X:100 is treated as the single base pair region X:100-100.
 *     In this case X:-100 is shorthand for X:1-100 and X:100- is X:100-<end>.
 *     (This is the standard bcftools region convention.)
 *
 *     When not set X:100 is considered to be X:100-<end> where <end> is
 *     the end of chromosome X (set to HTS_POS_MAX here).  X:100- and X:-100
 *     are invalid.
 *     (This is the standard samtools region convention.)
 *
 * Note the supplied string expects 1 based inclusive coordinates, but the
 * returned coordinates start from 0 and are half open, so pos0 is valid
 * for use in e.g. "for (pos0 = beg; pos0 < end; pos0++) {...}"
 *
 * On success a pointer to the byte after the end of the entire region
 *            specifier is returned (plus any trailing comma), and tid,
 *            beg & end will be set.
 * On failure NULL is returned.
 */
const char *hts_parse_region(const char *s, int *tid, hts_pos_t *beg,
                             hts_pos_t *end, hts_name2id_f getid, void *hdr,
                             int flags)
{
    if (!s || !tid || !beg || !end || !getid)
        return NULL;

    size_t s_len = strlen(s);
    kstring_t ks = { 0, 0, NULL };

    const char *colon = NULL, *comma = NULL;
    int quoted = 0;

    if (flags & HTS_PARSE_LIST)
        flags &= ~HTS_PARSE_THOUSANDS_SEP;
    else
        flags |= HTS_PARSE_THOUSANDS_SEP;

    const char *s_end = s + s_len;

    // Braced quoting of references is permitted to resolve ambiguities.
    if (*s == '{') {
        const char *close = memchr(s, '}', s_len);
        if (!close) {
            hts_log_error("Mismatching braces in \"%s\"", s);
            *tid = -1;
            return NULL;
        }
        s++;
        s_len--;
        if (close[1] == ':')
            colon = close+1;
        quoted = 1; // number of trailing characters to trim

        // Truncate to this item only, if appropriate.
        if (flags & HTS_PARSE_LIST) {
            comma = strchr(close, ',');
            if (comma) {
                s_len = comma-s;
                s_end = comma+1;
            }
        }
    } else {
        // Truncate to this item only, if appropriate.
        if (flags & HTS_PARSE_LIST) {
            comma = strchr(s, ',');
            if (comma) {
                s_len = comma-s;
                s_end = comma+1;
            }
        }

        colon = hts_memrchr(s, ':', s_len);
    }

    // No colon is simplest case; just check and return.
    if (colon == NULL) {
        *beg = 0; *end = HTS_POS_MAX;
        kputsn(s, s_len-quoted, &ks); // convert to nul terminated string
        if (!ks.s) {
            *tid = -2;
            return NULL;
        }

        *tid = getid(hdr, ks.s);
        free(ks.s);

        return *tid >= 0 ? s_end : NULL;
    }

    // Has a colon, but check whole name first.
    if (!quoted) {
        *beg = 0; *end = HTS_POS_MAX;
        kputsn(s, s_len, &ks); // convert to nul terminated string
        if (!ks.s) {
            *tid = -2;
            return NULL;
        }
        if ((*tid = getid(hdr, ks.s)) >= 0) {
            // Entire name matches, but also check this isn't
            // ambiguous.  eg we have ref chr1 and ref chr1:100-200
            // both present.
            ks.l = 0;
            kputsn(s, colon-s, &ks); // convert to nul terminated string
            if (!ks.s) {
                *tid = -2;
                return NULL;
            }
            if (getid(hdr, ks.s) >= 0) {
                free(ks.s);
                *tid = -1;
                hts_log_error("Range is ambiguous. "
                              "Use {%s} or {%.*s}%s instead",
                              s, (int)(colon-s), s, colon);
                return NULL;
            }
            free(ks.s);

            return s_end;
        }
        if (*tid < -1) // Failed to parse header
            return NULL;
    }

    // Quoted, or unquoted and whole string isn't a name.
    // Check the pre-colon part is valid.
    ks.l = 0;
    kputsn(s, colon-s-quoted, &ks); // convert to nul terminated string
    if (!ks.s) {
        *tid = -2;
        return NULL;
    }
    *tid = getid(hdr, ks.s);
    free(ks.s);
    if (*tid < 0)
        return NULL;

    // Finally parse the post-colon coordinates
    char *hyphen;
    *beg = hts_parse_decimal(colon+1, &hyphen, flags) - 1;
    if (*beg < 0) {
        if (*beg != -1 && *hyphen == '-' && colon[1] != '\0') {
            // User specified zero, but we're 1-based.
            hts_log_error("Coordinates must be > 0");
            return NULL;
        }
        if (isdigit_c(*hyphen) || *hyphen == '\0' || *hyphen == ',') {
            // interpret chr:-100 as chr:1-100
            *end = *beg==-1 ? HTS_POS_MAX : -(*beg+1);
            *beg = 0;
            return s_end;
        } else if (*beg < -1) {
            hts_log_error("Unexpected string \"%s\" after region", hyphen);
            return NULL;
        }
    }

    if (*hyphen == '\0' || ((flags & HTS_PARSE_LIST) && *hyphen == ',')) {
        *end = flags & HTS_PARSE_ONE_COORD ? *beg+1 : HTS_POS_MAX;
    } else if (*hyphen == '-') {
        *end = hts_parse_decimal(hyphen+1, &hyphen, flags);
        if (*hyphen != '\0' && *hyphen != ',') {
            hts_log_error("Unexpected string \"%s\" after region", hyphen);
            return NULL;
        }
    } else {
        hts_log_error("Unexpected string \"%s\" after region", hyphen);
        return NULL;
    }

    if (*end == 0)
        *end = HTS_POS_MAX; // interpret chr:100- as chr:100-<end>

    if (*beg >= *end) return NULL;

    return s_end;
}

// Next release we should mark this as deprecated?
// Use hts_parse_region above instead.
const char *hts_parse_reg64(const char *s, hts_pos_t *beg, hts_pos_t *end)
{
    char *hyphen;
    const char *colon = strrchr(s, ':');
    if (colon == NULL) {
        *beg = 0; *end = HTS_POS_MAX;
        return s + strlen(s);
    }

    *beg = hts_parse_decimal(colon+1, &hyphen, HTS_PARSE_THOUSANDS_SEP) - 1;
    if (*beg < 0) *beg = 0;

    if (*hyphen == '\0') *end = HTS_POS_MAX;
    else if (*hyphen == '-') *end = hts_parse_decimal(hyphen+1, NULL, HTS_PARSE_THOUSANDS_SEP);
    else return NULL;

    if (*beg >= *end) return NULL;
    return colon;
}

const char *hts_parse_reg(const char *s, int *beg, int *end)
{
    hts_pos_t beg64 = 0, end64 = 0;
    const char *colon = hts_parse_reg64(s, &beg64, &end64);
    if (beg64 > INT_MAX) {
        hts_log_error("Position %"PRId64" too large", beg64);
        return NULL;
    }
    if (end64 > INT_MAX) {
        if (end64 == HTS_POS_MAX) {
            end64 = INT_MAX;
        } else {
            hts_log_error("Position %"PRId64" too large", end64);
            return NULL;
        }
    }
    *beg = beg64;
    *end = end64;
    return colon;
}

hts_itr_t *hts_itr_querys(const hts_idx_t *idx, const char *reg, hts_name2id_f getid, void *hdr, hts_itr_query_func *itr_query, hts_readrec_func *readrec)
{
    int tid;
    hts_pos_t beg, end;

    if (strcmp(reg, ".") == 0)
        return itr_query(idx, HTS_IDX_START, 0, 0, readrec);
    else if (strcmp(reg, "*") == 0)
        return itr_query(idx, HTS_IDX_NOCOOR, 0, 0, readrec);

    if (!hts_parse_region(reg, &tid, &beg, &end, getid, hdr, HTS_PARSE_THOUSANDS_SEP))
        return NULL;

    return itr_query(idx, tid, beg, end, readrec);
}
