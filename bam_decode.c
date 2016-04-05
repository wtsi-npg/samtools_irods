/*  bam_decode.c -- split subcommand.

    Copyright (C) 2016 Genome Research Ltd.

    Author: Jennifer Liddle <js10@sanger.ac.uk>

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

#include <config.h>

#include <assert.h>
#include <htslib/sam.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <unistd.h>
#include <regex.h>
#include <htslib/khash.h>
#include <cram/sam_header.h>
#include "sam_opts.h"
#include "samtools.h"

#define DEFAULT_MAX_LOW_QUALITY_TO_CONVERT 15
#define DEFAULT_MAX_NO_CALLS 2
#define DEFAULT_MAX_MISMATCHES 1
#define DEFAULT_MIN_MISMATCH_DELTA 1
#define DEFAULT_BARCODE_TAG "BC"
#define DEFAULT_QUALITY_TAG "QT"

/*
 * details read from barcode file
 * Plus metrics information for each barcode
 */
typedef struct {
    char *seq;
    char *name;
    char *lib;
    char *sample;
    char *desc;
    int reads, pf_reads, perfect, pf_perfect, one_mismatch, pf_one_mismatch;
} bc_details_t;

/*
 * structure to hold options
 */
typedef struct {
    char* input_name;
    char* output_name;
    char* barcode_name;
    char *metrics_name;
    char *barcode_tag_name;
    char *quality_tag_name;
    bool verbose;
    int max_low_quality_to_convert;
    bool convert_low_quality;
    int max_no_calls;
    int max_mismatches;
    int min_mismatch_delta;
    bool change_read_name;
    char *argv_list;
    sam_global_args ga;
} opts_t;

/*
 * hold state information
 */
typedef struct {
    samFile* input_file;
    bam_hdr_t* input_header;
    samFile* output_file;
    bam_hdr_t* output_header;
    char * barcode_name;
    char *metrics_name;
    FILE *metricsFileHandle;
    char *barcode_tag_name;
    char *quality_tag_name;
    int tag_length;
    bool convert_low_quality;
    int max_low_quality_to_convert;
    int max_no_calls;
    int max_mismatches;
    int min_mismatch_delta;
    bool change_read_name;
    char *argv_list;
    bc_details_t *nullMetric;
} state_t;

// Create a hash map to hold barcode information
// Key is 'barcode sequence'
// Val is bc_details_t
KHASH_MAP_INIT_STR(bc, bc_details_t *)

static int cleanup_state(state_t* status);
static void cleanup_opts(opts_t* opts);

/*
 * display usage information
 */
static void usage(FILE *write_to)
{
    fprintf(write_to,
"Usage: samtools decode [options] filename\n"
"\n"
"Options:\n"
"  -o   --output                        output file [default: stdout]\n"
"  -v   --verbose                       verbose output\n"
"  -b   --barcode-file                  file containing barcodes\n"
"  -c   --convert-low-quality           Convert low quality bases in barcode read to 'N'\n"
"  -q   --max-low-quality-to-convert    Max low quality phred value to convert bases in barcode read to 'N'\n"
"  -n   --max-no-calls                  Max allowable number of no-calls in a barcode read before it is considered unmatchable\n"
"  -m   --max-mismatches                Maximum mismatches for a barcode to be considered a match\n"
"  -d   --min-mismatch-delta            Minimum difference between number of mismatches in the best and second best barcodes for\n"
"                                       a barcode to be considered a match\n"
"  -r   --change-read-name              Change the read name by adding #<barcode> suffix\n"
"  -t   --metrics-file                  Per-barcode and per-lane metrics written to this file\n"
"       --barcode-tag-name              Barcode tag name [default: " DEFAULT_BARCODE_TAG "]\n"
"       --quality-tag-name              Quality tag name [default: " DEFAULT_QUALITY_TAG "]\n"
);
    sam_global_opt_help(write_to, ".-.--");
}

/*
 * Takes the command line options and turns them into something we can understand
 */
static opts_t* parse_args(int argc, char *argv[])
{
    if (argc == 1) { usage(stdout); return NULL; }

    const char* optstring = "i:o:vqcb:n:m:d:t:z:y:";

    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS(0, '-', 0, '-', '-'),
        { "input",                      1, 0, 'i' },
        { "output",                     1, 0, 'o' },
        { "verbose",                    0, 0, 'v' },
        { "max-low-quality-to-convert", 1, 0, 'q' },
        { "convert-low-quality",        0, 0, 'c' },
        { "barcode-file",               1, 0, 'b' },
        { "max-no-calls",               1, 0, 'n' },
        { "max-mismatches",             1, 0, 'm' },
        { "min-mismatch-delta",         1, 0, 'd' },
        { "change-read-name",           0, 0, 'r' },
        { "metrics-file",               1, 0, 't' },
        { "barcode-tag-name",           1, 0, 'z' },
        { "quality-tag-name",           1, 0, 'y' },
        { NULL, 0, NULL, 0 }
    };

    opts_t* retval = calloc(sizeof(opts_t), 1);
    if (! retval ) { perror("cannot allocate option parsing memory"); return NULL; }
    retval->argv_list = stringify_argv(argc+1, argv-1);
    if (retval->argv_list[strlen(retval->argv_list)-1] == ' ') retval->argv_list[strlen(retval->argv_list)-1] = 0;

    sam_global_args_init(&retval->ga);

    // set defaults
    retval->max_low_quality_to_convert = DEFAULT_MAX_LOW_QUALITY_TO_CONVERT;
    retval->max_no_calls = DEFAULT_MAX_NO_CALLS;
    retval->max_mismatches = DEFAULT_MAX_MISMATCHES;
    retval->min_mismatch_delta = DEFAULT_MIN_MISMATCH_DELTA;
    retval->verbose = false;
    retval->convert_low_quality = false;
    retval->change_read_name = false;
    retval->barcode_tag_name = DEFAULT_BARCODE_TAG;
    retval->quality_tag_name = DEFAULT_QUALITY_TAG;

    int opt;
    while ((opt = getopt_long(argc, argv, optstring, lopts, NULL)) != -1) {
        switch (opt) {
        case 'i':   retval->input_name = strdup(optarg);
                    break;
        case 'o':   retval->output_name = strdup(optarg);
                    break;
        case 'b':   retval->barcode_name = strdup(optarg);
                    break;
        case 't':   retval->metrics_name = strdup(optarg);
                    break;
        case 'v':   retval->verbose = true;
                    break;
        case 'q':   retval->max_low_quality_to_convert = atoi(optarg);
                    break;
        case 'c':   retval->convert_low_quality = true;
                    break;
        case 'n':   retval->max_no_calls = atoi(optarg);
                    break;
        case 'm':   retval->max_mismatches = atoi(optarg);
                    break;
        case 'd':   retval->min_mismatch_delta = atoi(optarg);
                    break;
        case 'r':   retval->change_read_name = true;
                    break;
        case 'z':   retval->barcode_tag_name = strdup(optarg);
                    break;
        case 'y':   retval->quality_tag_name = strdup(optarg);
                    break;
        default:    if (parse_sam_global_opt(opt, optarg, lopts, &retval->ga) == 0) break;
            /* else fall-through */
        case '?':   usage(stdout); free(retval); return NULL;
        }
    }

    argc -= optind;
    argv += optind;
    optind = 0;

    printf("argc = %d\n", argc);
    printf("argv = %s\n", argv[0]);

    return retval;
}

//
// char *checkBarcodeQuality(char *barcode, char *quality);
//
// return a new barcode read string with low quality bases converted to 'N'
//
static char *checkBarcodeQuality(char * barcode, char *quality, int max_low_quality_to_convert)
{
    if (!quality) return strdup(barcode);

    if (!barcode || strlen(barcode) != strlen(quality)) {
        fprintf(stderr, "checkBarcodeQuality(): barcode and quality are different lengths\n");
        return NULL;
    }

    int mlq = max_low_quality_to_convert ? max_low_quality_to_convert : DEFAULT_MAX_LOW_QUALITY_TO_CONVERT;
    char *newBarcode = strdup(barcode);
    int i;
    for (i=0; i < strlen(quality); i++) {
        int qual = quality[i] - 33;

        if (qual <= mlq) {
            newBarcode[i] = 'N';
        } else {
            newBarcode[i] = barcode[i];
        }
    }

    return newBarcode;
}

// Set the initial state
static state_t* init(opts_t* opts)
{
    state_t* retval = calloc(sizeof(state_t), 1);
    if (!retval) {
        fprintf(stderr, "Out of memory");
        return NULL;
    }

    retval->argv_list = opts->argv_list;
    retval->nullMetric = calloc(1,sizeof(bc_details_t));

    retval->input_file = sam_open_format(opts->input_name, "rb", &opts->ga.in);
    if (!retval->input_file) {
        fprintf(stderr, "Could not open input file (%s)\n", opts->input_name);
        free(retval);
        return NULL;
    }
    retval->input_header = sam_hdr_read(retval->input_file);
    if (retval->input_header == NULL) {
        fprintf(stderr, "Could not read header for file '%s'\n",
                opts->input_name);
        cleanup_state(retval);
        return NULL;
    }

    if (opts->output_name) {
        retval->output_header = bam_hdr_dup(retval->input_header);

        retval->output_file = sam_open_format(opts->output_name, "wb", &opts->ga.out);
        if (retval->output_file == NULL) {
            fprintf(stderr, "Could not open output file: %s\n", opts->output_name);
            cleanup_state(retval);
            return NULL;
        }
    }

    char* dirsep = strrchr(opts->input_name, '/');
    char* input_base_name = strdup(dirsep? dirsep+1 : opts->input_name);
    if (!input_base_name) {
        fprintf(stderr, "Out of memory\n");
        cleanup_state(retval);
        return NULL;
    }
    char* extension = strrchr(input_base_name, '.');
    if (extension) *extension = '\0';

    free(input_base_name);

    retval->barcode_name = strdup(opts->barcode_name);
    retval->metrics_name = opts->metrics_name ? strdup(opts->metrics_name) : NULL;
    retval->max_no_calls = opts->max_no_calls;
    retval->max_mismatches = opts->max_mismatches;
    retval->min_mismatch_delta = opts->min_mismatch_delta;
    retval->convert_low_quality = opts->convert_low_quality;
    retval->max_low_quality_to_convert = opts->max_low_quality_to_convert;
    retval->change_read_name = opts->change_read_name;
    retval->barcode_tag_name = strdup(opts->barcode_tag_name);
    retval->quality_tag_name = strdup(opts->quality_tag_name);

    if (retval->metrics_name) {
        retval->metricsFileHandle = fopen(retval->metrics_name,"w");
        if (!retval->metricsFileHandle) {
            perror("Open failed");
        }
    } else {
        retval->metricsFileHandle = NULL;
    }

    return retval;
}

void writeMetricsLine(bc_details_t *bcd, state_t *state, int total_reads, int max_reads, int total_pf_reads, int max_pf_reads, int total_pf_reads_assigned, int nReads)
{
    FILE *f = state->metricsFileHandle;
    char *Nseq = NULL;

    if (!bcd->seq) {
        Nseq = calloc(1,state->tag_length+1);
        memset(Nseq,'N', state->tag_length);
    }

    fprintf(f, "%s\t", bcd->seq ? bcd->seq : Nseq);
    fprintf(f, "%s\t", bcd->name ? bcd->name : "" );
    fprintf(f, "%s\t", bcd->lib ? bcd->lib : "" );
    fprintf(f, "%s\t", bcd->sample ? bcd->sample : "" );
    fprintf(f, "%s\t", bcd->desc ? bcd->desc : "" );
    fprintf(f, "%d\t", bcd->reads);
    fprintf(f, "%d\t", bcd->pf_reads);
    fprintf(f, "%d\t", bcd->perfect);
    fprintf(f, "%d\t", bcd->pf_perfect);
    fprintf(f, "%d\t", bcd->one_mismatch);
    fprintf(f, "%d\t", bcd->pf_one_mismatch);
    fprintf(f, "%f\t", total_reads ? bcd->reads / (double)total_reads : 0 );
    fprintf(f, "%f\t", max_reads ? bcd->reads / (double)max_reads : 0 );
    fprintf(f, "%f\t", total_pf_reads ? bcd->pf_reads / (double)total_pf_reads : 0 );
    fprintf(f, "%f\t", max_pf_reads ? bcd->pf_reads / (double)max_pf_reads : 0 );
    fprintf(f, "%f", total_pf_reads_assigned ? bcd->pf_reads * nReads / (double)total_pf_reads_assigned : 0);
    fprintf(f, "\n");

    if (Nseq) free(Nseq);
}

/*
 *
 */
void writeMetrics(khash_t(bc) *barcodeHash, state_t *state)
{
    khiter_t iter;
    int total_reads = state->nullMetric->reads;
    int total_pf_reads = state->nullMetric->pf_reads;
    int total_pf_reads_assigned = 0;
    int max_reads = state->nullMetric->reads;
    int max_pf_reads = state->nullMetric->pf_reads;
    int nReads = 0;

    // first loop to count things
    for (iter = kh_begin(barcodeHash); iter != kh_end(barcodeHash); ++iter) {
        if (!kh_exist(barcodeHash,iter)) continue;
        bc_details_t *bcd = kh_val(barcodeHash,iter);
        total_reads += bcd->reads;
        total_pf_reads += bcd->pf_reads;
        total_pf_reads_assigned += bcd->pf_reads;
        if (max_reads < bcd->reads) max_reads = bcd->reads;
        if (max_pf_reads < bcd->pf_reads) max_pf_reads = bcd->pf_reads;
        nReads++;
    }

    // print header
    fprintf(state->metricsFileHandle, "BARCODE\t");
    fprintf(state->metricsFileHandle, "BARCODE_NAME\t");
    fprintf(state->metricsFileHandle, "LIBRARY_NAME\t");
    fprintf(state->metricsFileHandle, "SAMPLE_NAME\t");
    fprintf(state->metricsFileHandle, "DESCRIPTION\t");
    fprintf(state->metricsFileHandle, "READS\t");
    fprintf(state->metricsFileHandle, "PF_READS\t");
    fprintf(state->metricsFileHandle, "PERFECT_MATCHES\t");
    fprintf(state->metricsFileHandle, "PF_PERFECT_MATCHES\t");
    fprintf(state->metricsFileHandle, "ONE_MISMATCH_MATCHES\t");
    fprintf(state->metricsFileHandle, "PF_ONE_MISMATCH_MATCHES\t");
    fprintf(state->metricsFileHandle, "PCT_MATCHES\t");
    fprintf(state->metricsFileHandle, "RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT\t");
    fprintf(state->metricsFileHandle, "PF_PCT_MATCHES\t");
    fprintf(state->metricsFileHandle, "PF_RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT\t");
    fprintf(state->metricsFileHandle, "PF_NORMALIZED_MATCHES\n");


    // second loop to print things
    for (iter = kh_begin(barcodeHash); iter != kh_end(barcodeHash); ++iter) {
        if (!kh_exist(barcodeHash,iter)) continue;
        bc_details_t *bcd = kh_val(barcodeHash,iter);
        writeMetricsLine(bcd, state, total_reads, max_reads, total_pf_reads, max_pf_reads, total_pf_reads_assigned, nReads);
    }
    writeMetricsLine(state->nullMetric, state, total_reads, max_reads, total_pf_reads, max_pf_reads, 0, nReads);
}

/*
 * Read the barcode file into a hash
 */
khash_t(bc) *loadBarcodeFile(state_t *state)
{
    khash_t(bc) *barcodeHash = kh_init(bc);
    FILE *fh = fopen(state->barcode_name,"r");
    if (!fh) {
        fprintf(stderr,"ERROR: Can't open barcode file %s\n", state->barcode_name);
        return NULL;
    }
    
    char *buf = NULL;
    int tag_length = 0;
    size_t n;
    getline(&buf,&n,fh);    // burn first line which is a header
    free(buf); buf=NULL;

    while (getline(&buf, &n, fh) > 0) {
        char *s;
        int res;
        if (buf[strlen(buf)-1] == '\n') buf[strlen(buf)-1]=0;   // remove trailing lf
        bc_details_t *bcd = calloc(1,sizeof(bc_details_t));
        s = strtok(buf,"\t");  bcd->seq     = strdup(s);
        s = strtok(NULL,"\t"); bcd->name    = strdup(s);
        s = strtok(NULL,"\t"); bcd->lib     = strdup(s);
        s = strtok(NULL,"\t"); bcd->sample  = strdup(s);
        s = strtok(NULL,"\t"); bcd->desc    = strdup(s);
        int iter = kh_put(bc, barcodeHash, bcd->seq, &res);
        kh_value(barcodeHash,iter) = bcd;
        free(buf); buf=NULL;

        if (tag_length == 0) {
            tag_length = strlen(bcd->seq);
        } else {
            if (tag_length != strlen(bcd->seq)) {
                fprintf(stderr,"ERROR: Tag '%s' is a different length to the previous tag\n", bcd->seq);
                return NULL;
            }
        }
    }

    state->tag_length = tag_length;     // save this - we'll need it later
    fclose(fh);
    return barcodeHash;
}

/*
 * return true if base is a noCall
 */
int isNoCall(char b)
{
    return b=='N' || b=='n' || b=='.';
}

/*
 * Count the number of noCalls in a sequence
 */
static int noCalls(char *s)
{
    int n=0;
    while (*s) {
        if (isNoCall(*s++)) n++;
    }
    return n;
}
        
/*
 * count number of mismatches between two sequences
 * (ignoring noCalls)
 */
static int countMismatches(char *tag, char *barcode)
{
    char *t, *b;;
    int n = 0;
    for (t=tag, b=barcode; *t; t++, b++) {
        if (!isNoCall(*t)) {
            if (!isNoCall(*b)) {
                if (*t != *b) {
                    n++;
                }
            }
        }
    }
    return n;
}

/*
 * find the best match in the barcode (tag) file for a given barcode
 * return the tag, if a match found, else return NULL
 */
static char *findBestMatch(char *barcode, khash_t(bc) *barcodeHash, state_t *state)
{
    int bcLen = state->tag_length;   // size of barcode sequence in barcode file
    char *best_match = NULL;
    int nmBest = bcLen;             // number of mismatches (best)
    int nm2Best = bcLen;            // number of mismatches (second best)
    int nCalls = noCalls(barcode);

    // for each tag in barcodeHash
    khiter_t iter;
    for (iter = kh_begin(barcodeHash); iter != kh_end(barcodeHash); ++iter) {
        if (!kh_exist(barcodeHash,iter)) continue;
        bc_details_t *bcd = kh_val(barcodeHash,iter);

        int nMismatches = countMismatches(bcd->seq, barcode);
        if (nMismatches < nmBest) {
            if (best_match) nm2Best = nmBest;
            nmBest = nMismatches;
            best_match = bcd->seq;
        } else {
            if (nMismatches < nm2Best) nm2Best = nMismatches;
        }
    }

    bool matched = best_match &&
                   nCalls <= state->max_no_calls &&
                   nmBest <= state->max_mismatches &&
                   nm2Best - nmBest >= state->min_mismatch_delta;

    if (matched) return strdup(best_match);
    return NULL;
}

/*
 * Update the metrics information
 */
void updateMetrics(bc_details_t *bcd, char *seq, bool isPf)
{
    int n = 99;
    if (seq) n = countMismatches(bcd->seq, seq);

    bcd->reads++;
    if (isPf) bcd->pf_reads++;

    if (n==0) {     // count perfect matches
        bcd->perfect++;
        if (isPf) bcd->pf_perfect++;
    }

    if (n==1) {     // count out-by-one matches
        bcd->one_mismatch++;
        if (isPf) bcd->pf_one_mismatch++;
    }
        
}

/*
 * find the best match in the barcode (tag) file, and return the corresponding barcode name
 * return NULL if no match found
 */
static char *findBarcodeName(char *barcode, khash_t(bc) *barcodeHash, state_t *state, bool isPf)
{
    char *seq = findBestMatch(barcode, barcodeHash, state);
    if (!seq) { updateMetrics(state->nullMetric, seq, isPf); return NULL; }
    khint_t k = kh_get(bc,barcodeHash,seq);
    assert(k != kh_end(barcodeHash));
    bc_details_t *bcd = kh_val(barcodeHash,k);
    updateMetrics(bcd, barcode, isPf);
    return bcd->name;
}

/*
 * make a new tag by appending #<name> to the old tag
 */
char *makeNewTag(bam1_t *rec, char *tag, char *name)
{
    char *rg = "";
    uint8_t *p = bam_aux_get(rec,tag);
    if (p) rg = bam_aux2Z(p);
    char *newtag = malloc(strlen(rg)+1+strlen(name)+1);
    strcpy(newtag, rg);
    strcat(newtag,"#");
    strcat(newtag, name);
    return newtag;
}

/*
 * Change the read name by adding "#<suffix>"
 */
void add_suffix(bam1_t *rec, char *suffix)
{
    int oldqlen = strlen((char *)rec->data);
    int newlen = rec->l_data + strlen(suffix) + 1;

    if (newlen > rec->m_data) {
        rec->m_data = newlen;
        kroundup32(rec->m_data);
        rec->data = (uint8_t *)realloc(rec->data, rec->m_data);
    }
    memmove(rec->data + oldqlen + strlen(suffix) + 1,
            rec->data + oldqlen,
            rec->l_data - oldqlen);
    rec->data[oldqlen] = '#';
    memmove(rec->data + oldqlen + 1,
            suffix,
            strlen(suffix) + 1);
    rec->l_data = newlen;
    rec->core.l_qname += strlen(suffix) + 1;
}

/*
 * Add a new @RG line to the header
 */
void addNewRG(SAM_hdr *sh, char *entry, char *bcname, char *lib, char *sample, char *desc)
{
    char *saveptr;
    char *p = strtok_r(entry,"\t",&saveptr);
    char *newtag = malloc(strlen(p)+1+strlen(bcname)+1);
    strcpy(newtag, p);
    strcat(newtag,"#");
    strcat(newtag, bcname);
    sam_hdr_add(sh, "RG", "ID", newtag, NULL, NULL);

    SAM_hdr_type *hdr = sam_hdr_find(sh, "RG", "ID", newtag);
    while (1) {
        char *pu = NULL;
        char *t = strtok_r(NULL, ":", &saveptr);
        if (!t) break;
        char *v = strtok_r(NULL, "\t", &saveptr);
        if (!v) break;

        // handle special cases
        if (strcmp(t,"PU") == 0) {
            // add #bcname
            pu = malloc(strlen(v) + 1 + strlen(bcname) + 1);
            strcpy(pu, v); strcat(pu,"#"); strcat(pu,bcname);
            v = pu;
        }
        if (strcmp(t,"LB") == 0) {
            if (lib) v = lib;        // use library name
        }
        if (strcmp(t,"DS") == 0) {
            if (desc) v = desc;       // use desc
        }
        if (strcmp(t,"SM") == 0) {
            if (sample) v = sample;     // use sample name
        }
        sam_hdr_update(sh, hdr, t, v, NULL);
        if (pu) free(pu);
    }
    free(newtag);
}

/*
 * for each "@RG ID:x" in the header, replace with
 * "@RG IDx#barcode" for each barcode
 *
 * And don't forget to add a @PG header
 */ 
void changeHeader(khash_t(bc) *barcodeHash, state_t *state)
{
    bam_hdr_t *h = state->output_header;
    SAM_hdr *sh = sam_hdr_parse_(h->text, h->l_text);
    char **rgArray = malloc(sizeof(char*) * sh->nrg);
    int nrg = sh->nrg;
    int n;

    sam_hdr_add_PG(sh, "samtools", "VN", samtools_version(), "CL", state->argv_list, NULL);

    // store the RG names
    for (n=0; n < sh->nrg; n++) {
        // store the names and tags as a string <name>:<tag>:<val>:<tag>:<val>...
        // eg 1:PL:Illumina:PU:110608_HS19

        // first pass to determine size of string required
        int sz=strlen(sh->rg[n].name)+1;
        SAM_hdr_tag *rgtag = sh->rg[n].tag;
        while (rgtag) {
            if (strncmp(rgtag->str,"ID:",3)) {  // ignore name
                sz += 3 + strlen(rgtag->str) + 1;
            }
            rgtag = rgtag->next;
        }
        char *entry = malloc(sz+1);

        // second pass to create string
        strcpy(entry,sh->rg[n].name);
        rgtag = sh->rg[n].tag;
        while (rgtag) {
            if (strncmp(rgtag->str,"ID:",3)) {  // ignore name
                strcat(entry,"\t");
                strcat(entry,rgtag->str);
            }
            rgtag = rgtag->next;
        }
        rgArray[n] = entry;
    }

    // delete the old RG lines
    sam_hdr_del(sh, "RG", NULL, NULL);

    // add the new ones
    for (n=0; n<nrg; n++) {
        char *entry = strdup(rgArray[n]);
        addNewRG(sh, entry, "0", NULL, NULL, NULL);
        free(entry);

        // for each tag in barcodeHash
        khiter_t iter;
        for (iter = kh_begin(barcodeHash); iter != kh_end(barcodeHash); ++iter) {
            if (!kh_exist(barcodeHash,iter)) continue;
            bc_details_t *bcd = kh_val(barcodeHash,iter);

            char *entry = strdup(rgArray[n]);
            addNewRG(sh, entry, bcd->name, bcd->lib, bcd->sample, bcd->desc);
            free(entry);
        }
    }

    for (n=0; n<nrg; n++) {
        free(rgArray[n]);
    }
    free(rgArray);

    free(h->text);
    sam_hdr_rebuild(sh);
    h->text = strdup(sam_hdr_str(sh));
    h->l_text = sam_hdr_length(sh);
    sam_hdr_free(sh);
}

/*
 * Main code
 */
static bool decode(state_t* state)
{
    bam1_t* file_read = bam_init1();
    bam1_t* paired_read = bam_init1();
    int r, r2;
    khash_t(bc) *barcodeHash;
    char *name;

    barcodeHash = loadBarcodeFile(state);

    changeHeader(barcodeHash, state);

    if (sam_hdr_write(state->output_file, state->output_header) != 0) {
        fprintf(stderr, "Could not write output file header\n");
        return false;
    }

    /*
     * For each BAM record
     */
    while (1) {
        r = sam_read1(state->input_file, state->input_header, file_read);
    if (r < 0) break;

        // look for barcode tag
        uint8_t *p = bam_aux_get(file_read,state->barcode_tag_name);
        if (p) {
            char *seq = bam_aux2Z(p);
            char *newseq = strdup(seq);
            if (state->convert_low_quality) {
                uint8_t *q = bam_aux_get(file_read,state->quality_tag_name);
                if (q) {
                    char *qual = bam_aux2Z(q);
                    free(newseq);
                    newseq = checkBarcodeQuality(seq,qual,state->max_low_quality_to_convert);
                }
            }
            name = findBarcodeName(newseq,barcodeHash,state,!(file_read->core.flag & BAM_FQCFAIL));
            if (!name) name = "0";
            char * newtag = makeNewTag(file_read,"RG",name);
            bam_aux_update_str(file_read,"RG",strlen(newtag)+1,(uint8_t*)newtag);
            free(newtag);
            if (state->change_read_name) add_suffix(file_read, name);
            free(newseq);
        }
        r2 = sam_write1(state->output_file, state->output_header, file_read);
        if (r2 < 0) {
            fprintf(stderr, "Could not write sequence\n");
            return false;
        }
        
        if (file_read->core.flag & BAM_FPAIRED) {
            r = sam_read1(state->input_file, state->input_header, paired_read);
            if (p) {
                char *newtag = makeNewTag(paired_read,"RG",name);
                bam_aux_update_str(paired_read,"RG",strlen(newtag)+1,(uint8_t*)newtag);
                free(newtag);
            }
            if (state->change_read_name) add_suffix(paired_read, name);
            r2 = sam_write1(state->output_file, state->output_header, paired_read);
            if (r2 < 0) {
                fprintf(stderr, "Could not write sequence\n");
                return false;
            }
        }
    }

    if (state->metricsFileHandle) writeMetrics(barcodeHash, state);

    return true;
}

static int cleanup_state(state_t* status)
{
    int ret = 0;

    if (!status) return 0;
    if (status->output_header) bam_hdr_destroy(status->output_header);
    if (status->output_file) ret |= sam_close(status->output_file);
    sam_close(status->input_file);
    bam_hdr_destroy(status->input_header);
    if (status->metricsFileHandle) fclose(status->metricsFileHandle); 
    status->metricsFileHandle=NULL;
    free(status);

    return ret;
}

static void cleanup_opts(opts_t* opts)
{
    if (!opts) return;
    free(opts->input_name);
    free(opts->output_name);
    sam_global_args_free(&opts->ga);
    free(opts);
}

int main_decode(int argc, char *argv[])
{
    int ret = 1;
    opts_t* opts = parse_args(argc, argv);
    if (opts) {
        if (opts->verbose) fprintf(stderr, "options parsed ok\n");
        state_t* status = init(opts);
        if (status) {
            if (opts->verbose) fprintf(stderr, "state initialised ok\n");
            if (decode(status)) ret = 0;
            ret |= (cleanup_state(status) != 0);
        }
    }
    cleanup_opts(opts);

    return ret;
}
