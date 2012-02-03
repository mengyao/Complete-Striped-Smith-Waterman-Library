/*  main.c
 *  Created by Mengyao Zhao on 06/23/11.
 *	Version 0.1.4
 *  Last revision by Mengyao Zhao on 02/03/12.
 *	New features: make weight as options 
 */

#include <stdlib.h>
#include <stdint.h>
#include <emmintrin.h>
#include <zlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include "ssw.h"
#include "kseq.h"
#include "banded_sw.h"

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

KSEQ_INIT(gzFile, gzread)

char* seq_reverse(const char* seq, int32_t end)	/* end is 0-based alignment ending position */	
{									
	char* reverse = (char*)calloc(end + 2, sizeof(char*));	
	int32_t start = 0;
	reverse[end + 1] = '\0';				
	while (LIKELY(start <= end)) {			
		reverse[start] = seq[end];		
		reverse[end] = seq[start];		
		++ start;					
		-- end;						
	}								
	return reverse;					
}									

int main (int argc, char * const argv[]) {
	clock_t start, end;
	float cpu_time;
	gzFile ref_fp;
	kseq_t *ref_seq;
	int32_t l, m, k, match = 2, mismatch = 2, insert_open = 3, insert_extention = 1, delet_open = 3, delet_extention = 1, path = 0;
	int8_t mat[25];

	// Parse command line.
	while ((l = getopt(argc, argv, "m:x:i:e:d:f:p")) >= 0) {
		switch (l) {
			case 'm': match = atoi(optarg); break;
			case 'x': mismatch = atoi(optarg); break;
			case 'i': insert_open = atoi(optarg); break;
			case 'e': insert_extention = atoi(optarg); break;
			case 'd': delet_open = atoi(optarg); break;
			case 'f': delet_extention = atoi(optarg); break;
			case 'p': path = 1; break;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: ssw_test [-m weight match] [-x abs(weight mismatch)] [-i abs(weight insert_open)] [-e abs(weight insert_extention)] [-d abs(weight delet_open)] [-f abs(weight delet_extention)] [-p] <target.fa> <query.fa>\n");
		return 1;
	}

	/* This table is used to transform nucleotide letters into numbers. */
	int8_t nt_table[128] = {
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4 
	};

	// initialize scoring matrix for genome sequences
	for (l = k = 0; LIKELY(l < 4); ++l) {
		for (m = 0; LIKELY(m < 4); ++m)
			mat[k++] = l == m ? match : -mismatch;	/* weight_match : -weight_mismatch */
		mat[k++] = 0; // ambiguous base
	}
	for (m = 0; LIKELY(m < 5); ++m) mat[k++] = 0;

	ref_fp = gzopen(argv[optind], "r");
	ref_seq = kseq_init(ref_fp);

	start = clock();
	while ((l = kseq_read(ref_seq)) >= 0) {
		gzFile read_fp;
		kseq_t *read_seq;
		printf("ref_name: %s\n", ref_seq->name.s);
		read_fp = gzopen(argv[optind + 1], "r");
		read_seq = kseq_init(read_fp);
		int32_t refLen = strlen(ref_seq->seq.s);
		char* ref_reverse = seq_reverse(ref_seq->seq.s, refLen - 1); 
		while ((m = kseq_read(read_seq)) >= 0) {
			char *read_reverse;
			int32_t readLen, word = 0;
			alignment_end *bests, *bests_reverse;
			__m128i *vP;
			printf("read_name: %s\n", read_seq->name.s);
			printf("read_seq: %s\n", read_seq->seq.s); 
			readLen = strlen(read_seq->seq.s);
			vP = qP_byte(read_seq->seq.s, nt_table, mat, 5, 4);
			bests = sw_sse2_byte(ref_seq->seq.s, nt_table, refLen, readLen, insert_open, insert_extention, delet_open, delet_extention, vP, 0, 4);
			if (bests[0].score == 255) {
				vP = qP_word(read_seq->seq.s, nt_table, mat, 5);
				bests = sw_sse2_word(ref_seq->seq.s, nt_table, refLen, readLen, insert_open, insert_extention, delet_open, delet_extention, vP, 0);
				word = 1;
			}
			free(vP);
			
			if (bests[0].score != 0) {
				char* cigar1;
				int32_t begin_ref, begin_read, band_width;
				read_reverse = seq_reverse(read_seq->seq.s, bests[0].read);
				if (word == 0) {
					vP = qP_byte(read_reverse, nt_table, mat, 5, 4);
					bests_reverse = sw_sse2_byte(ref_reverse + refLen - bests[0].ref - 1, nt_table, bests[0].ref + 1, bests[0].read + 1, insert_open, insert_extention, delet_open, delet_extention, vP, bests[0].score, 4);
				} else {
					vP = qP_word(read_reverse, nt_table, mat, 5);
					bests_reverse = sw_sse2_word(ref_reverse + refLen - bests[0].ref - 1, nt_table, bests[0].ref + 1, bests[0].read + 1, insert_open, insert_extention, delet_open, delet_extention, vP, bests[0].score);
				}
				free(vP);
				free(read_reverse);
			
				begin_ref = bests[0].ref - bests_reverse[0].ref; 
				begin_read = bests[0].read - bests_reverse[0].read; 
				band_width = abs(bests_reverse[0].ref - bests_reverse[0].read) + 1;
			
				fprintf(stdout, "max score: %d, 2nd score: %d, begin_ref: %d, begin_read: %d\n", bests[0].score, bests[1].score, begin_ref + 1, begin_read + 1);
				if (path == 1) {
					if (bests[0].score != bests[1].score) {
						cigar1 = banded_sw(ref_seq->seq.s + begin_ref, read_seq->seq.s + begin_read, bests_reverse[0].ref + 1, bests_reverse[0].read + 1, bests[0].score, match, mismatch, insert_open, insert_extention, delet_open, delet_extention, band_width, nt_table, mat, 5);
						if (cigar1 != 0) {
							fprintf(stdout, "cigar: %s\n", cigar1);
						} else fprintf(stdout, "No alignment is available.\n");	
						free(cigar1);		
					} else fprintf(stdout, "Two alignments available.\n");
				}
			}else fprintf(stdout, "No alignment found for this read.\n");
		}
		free(ref_reverse);
		kseq_destroy(read_seq);
		gzclose(read_fp);
	}
	end = clock();
	cpu_time = ((float) (end - start)) / CLOCKS_PER_SEC;
	fprintf(stdout, "CPU time: %f seconds\n", cpu_time);

	kseq_destroy(ref_seq);
	gzclose(ref_fp);
	return 0;
}
