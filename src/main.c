/*  main.c
 *  Created by Mengyao Zhao on 06/23/11.
 *	Version 0.1.4
 *  Last revision by Mengyao Zhao on 01/24/12.
 *	New features: weight matrix is extracted 
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

KSEQ_INIT(gzFile, gzread)

char* seq_reverse(const char* seq, int32_t end)	/* end is 0-based alignment ending position */	
{									
	char* reverse = (char*)calloc(end + 2, sizeof(char*));	
	int32_t start = 0;
	reverse[end + 1] = '\0';				
	while (start <= end) {			
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
	int32_t l, m, k, match, mismatch, insert_open, insert_extention, delet_open, delet_extention;
	int8_t mat[25];
	//char ref[50], read[50];

	// Parse command line.
	while ((l = getopt(argc, argv, "m:x:i:e:d:f:")) >= 0) {
		switch (l) {
		//	case 'd': ref = optarg; break;
		//	case 'q': read = optarg; break;
			case 'm': match = atoi(optarg); break;
			case 'x': mismatch = atoi(optarg); break;
			case 'i': insert_open = atoi(optarg); break;
			case 'e': insert_extention = atoi(optarg); break;
			case 'd': delet_open = atoi(optarg); break;
			case 'f': delet_extention = atoi(optarg); break;
		//	case 'f': forward_only = 1; break;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: ssw_test [-m weight match] [-x abs(weight mismatch)] [-i abs(weight insert_open)] [-e abs(weight insert_extention)] [-d abs(weight delet_open)] [-f abs(weight delet_extention)] <target.fa> <query.fa>\n");
		return 1;
	}
	// validate argument count
/*	if (argc < 2) {
		fprintf (stderr, "USAGE: %s <ref.fa/fq> <reads.fa/fq>\n", argv[0]);
		exit (1);
	}*/


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
	for (l = k = 0; l < 4; ++l) {
		for (m = 0; m < 4; ++m)
			mat[k++] = l == m ? match : -mismatch;	/* weight_match : -weight_mismatch */
	//		mat[k++] = l == m ? 2 : -2;	/* weight_match : -weight_mismatch */
		mat[k++] = 0; // ambiguous base
	}
	for (m = 0; m < 5; ++m) mat[k++] = 0;

	//fprintf(stderr, "argv[1]: %s\n", argv[1]);
	ref_fp = gzopen(argv[optind], "r");
	//ref_fp = gzopen(ref, "r");
	ref_seq = kseq_init(ref_fp);

	start = clock();
	while ((l = kseq_read(ref_seq)) >= 0) {
		gzFile read_fp;
		kseq_t *read_seq;
		printf("ref_name: %s\n", ref_seq->name.s);
		read_fp = gzopen(argv[optind + 1], "r");
	//	read_fp = gzopen(read, "r");
		read_seq = kseq_init(read_fp);
		int32_t refLen = strlen(ref_seq->seq.s);
		char* ref_reverse = seq_reverse(ref_seq->seq.s, refLen - 1); 
		fprintf(stderr, "reverse_ref: %s\n", ref_reverse); 										
		while ((m = kseq_read(read_seq)) >= 0) {
			char *read_reverse;
			alignment_end *bests, *bests_reverse;
			printf("read_name: %s\n", read_seq->name.s);
			printf("read_seq: %s\n", read_seq->seq.s); 
			
			int32_t readLen = strlen(read_seq->seq.s);
			__m128i* vProfile = queryProfile_constructor(read_seq->seq.s, nt_table, mat, 5, 4);
			bests = sw_sse2_byte(ref_seq->seq.s, nt_table, refLen, readLen, insert_open, insert_extention, delet_open, delet_extention, vProfile, 0, 4);
		//	bests = sw_sse2_byte(ref_seq->seq.s, nt_table, refLen, readLen, 3, 1, 3, 1, vProfile, 0, 4);
			free(vProfile);
			
			if (bests[0].score != 0) {
				char* cigar1;
				int32_t begin_ref, begin_read, band_width;
				read_reverse = seq_reverse(read_seq->seq.s, bests[0].read);
				fprintf(stderr, "reverse_read: %s\n", read_reverse); 										
				vProfile = queryProfile_constructor(read_reverse, nt_table, mat, 5, 4);
				bests_reverse = sw_sse2_byte(ref_reverse + refLen - bests[0].ref - 1, nt_table, bests[0].ref + 1, bests[0].read + 1, insert_open, insert_extention, delet_open, delet_extention, vProfile, bests[0].score, 4);
				free(vProfile);
				free(read_reverse);
			
				begin_ref = bests[0].ref - bests_reverse[0].ref, begin_read = bests[0].read - bests_reverse[0].read, band_width = abs(bests_reverse[0].ref - bests_reverse[0].read);
			
				fprintf(stdout, "max score: %d, 2nd score: %d, end_ref: %d, end_read: %d\nbegin_ref: %d, begin_read: %d\n", 
				bests[0].score, bests[1].score, bests[0].ref + 1, bests[0].read + 1, begin_ref + 1, begin_read + 1);
				if (bests[0].score != bests[1].score) {
					cigar1 = banded_sw(ref_seq->seq.s + begin_ref, read_seq->seq.s + begin_read, bests_reverse[0].ref + 1, bests_reverse[0].read + 1, match, mismatch, insert_open, insert_extention, delet_open, delet_extention, band_width, nt_table, mat, 5);
					if (cigar1 != 0) {
						fprintf(stdout, "cigar: %s\n", cigar1);
					} else fprintf(stdout, "No alignment is available.\n");	
					free(cigar1);		
				} else fprintf(stdout, "Two alignments available.\n");
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
