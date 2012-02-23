/*  main.c
 *  Created by Mengyao Zhao on 06/23/11.
 *	Version 0.1.4
 *  Last revision by Mengyao Zhao on 02/21/12.
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

alignment_end* find_bests(char* ref_seq,
		   char* read_seq, 
	  	   int8_t* ref_table, 
		   int8_t* read_table,
		   int8_t* mat,
	  	   int32_t n, 
	 	   int32_t refLen,
		   int32_t readLen,
		   int32_t insert_open, 
	 	   int32_t insert_extension,
	 	   int32_t delet_open,
	 	   int32_t delet_extension,
			int32_t* word){

	__m128i *vP;
	alignment_end* bests;
	if (*word == 0) {
		vP = qP_byte(read_seq, read_table, mat, n, 4);
		bests = sw_sse2_byte(ref_seq, ref_table, refLen, readLen, insert_open, insert_extension, delet_open, delet_extension, vP, 0, 4);
	}
	if (*word == 1 || bests[0].score == 255) {
		vP = qP_word(read_seq, read_table, mat, n);
		bests = sw_sse2_word(ref_seq, ref_table, refLen, readLen, insert_open, insert_extension, delet_open, delet_extension, vP, 0);
		*word = 1;
	}
	free(vP);
	return bests;	
}

void align(char* ref_seq,	 
		   char* read_seq, 
	  	   int8_t* ref_table,
		   int8_t* read_table, 
		   int8_t* mat,
	 	   char* ref_reverse,
	  	   int32_t n, 
	 	   int32_t refLen,
		   int32_t match,
		   int32_t mismatch, 
		   int32_t insert_open, 
	 	   int32_t insert_extension,
	 	   int32_t delet_open,
	 	   int32_t delet_extension,
	 	   int32_t path,
			int32_t word,
			alignment_end* bests){

	char *read_reverse;
	alignment_end *bests_reverse;
	__m128i *vP;
	printf("read_seq: %s\n", read_seq);
 
	if (bests[0].score != 0) {
		char* cigar1;
		int32_t begin_ref, begin_read, band_width;
		read_reverse = seq_reverse(read_seq, bests[0].read);
		if (word == 0) {
			vP = qP_byte(read_reverse, read_table, mat, n, 4);
			bests_reverse = sw_sse2_byte(ref_reverse + refLen - bests[0].ref - 1, ref_table, bests[0].ref + 1, bests[0].read + 1, insert_open, insert_extension, delet_open, delet_extension, vP, bests[0].score, 4);
		} else {
			vP = qP_word(read_reverse, read_table, mat, n);
			bests_reverse = sw_sse2_word(ref_reverse + refLen - bests[0].ref - 1, ref_table, bests[0].ref + 1, bests[0].read + 1, insert_open, insert_extension, delet_open, delet_extension, vP, bests[0].score);
		}
		free(vP);
		free(read_reverse);
	
		begin_ref = bests[0].ref - bests_reverse[0].ref; 
		begin_read = bests[0].read - bests_reverse[0].read; 
		band_width = abs(bests_reverse[0].ref - bests_reverse[0].read) + 1;
	
		fprintf(stdout, "max score: %d, 2nd score: %d, begin_ref: %d, begin_read: %d, end_ref: %d, end_read: %d\n", bests[0].score, bests[1].score, begin_ref + 1, begin_read + 1, bests[0].ref + 1, bests[0].read + 1);
		if (path == 1) {
			cigar1 = banded_sw(ref_seq + begin_ref, read_seq + begin_read, bests_reverse[0].ref + 1, bests_reverse[0].read + 1, bests[0].score, match, mismatch, insert_open, insert_extension, delet_open, delet_extension, band_width, ref_table, read_table, mat, n);
			if (cigar1 != 0) {
				fprintf(stdout, "cigar: %s\n", cigar1);
			} else fprintf(stdout, "No alignment is available.\n");	
			free(cigar1);		
		}
	}else fprintf(stdout, "No alignment is found for this read.\n");
}

int main (int argc, char * const argv[]) {
	clock_t start, end;
	float cpu_time;
	gzFile ref_fp;
	kseq_t *ref_seq;
	int32_t l, m, k, n = 5, match = 2, mismatch = 2, insert_open = 3, insert_extension = 1, delet_open = 3, delet_extension = 1, path = 0, reverse = 0;
	int8_t* mat = (int8_t*)calloc(25, sizeof(int8_t));
	char mat_name[16];
	mat_name[0] = '\0';

	/* This table is used to transform amino acid letters into numbers. */
	int8_t aa_table[128] = {
		23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 
		23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 
		23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
		23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 
		23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23, 
		14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23, 
		23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23, 
		14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23 
	};

	/* This table is used to transform nucleotide letters into numbers. */
	int8_t nt_table[128] = {
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4 
	};
	
	int8_t* table = nt_table;

	// initialize scoring matrix for genome sequences
	for (l = k = 0; LIKELY(l < 4); ++l) {
		for (m = 0; LIKELY(m < 4); ++m)
			mat[k++] = l == m ? match : -mismatch;	/* weight_match : -weight_mismatch */
		mat[k++] = 0; // ambiguous base
	}
	for (m = 0; LIKELY(m < 5); ++m) mat[k++] = 0;

	// Parse command line.
	while ((l = getopt(argc, argv, "m:x:i:e:d:f:a:cr")) >= 0) {
		switch (l) {
			case 'm': match = atoi(optarg); break;
			case 'x': mismatch = atoi(optarg); break;
			case 'i': insert_open = atoi(optarg); break;
			case 'e': insert_extension = atoi(optarg); break;
			case 'd': delet_open = atoi(optarg); break;
			case 'f': delet_extension = atoi(optarg); break;
			case 'a': strcpy(mat_name, optarg); break;
			case 'c': path = 1; break;
			case 'r': reverse = 1; break;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: ssw_test [options] ... <target.fa> <query.fa>\n");	
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "\t-m N\tN is a positive integer for weight match in genome sequence alignment.\n");
		fprintf(stderr, "\t-x N\tN is a positive integer. -N will be used as weight mismatch in genome sequence alignment.\n");
		fprintf(stderr, "\t-i N\tN is a positive integer. -N will be used as the weight for the insertion opening.\n");
		fprintf(stderr, "\t-e N\tN is a positive integer. -N will be used as the weight for the insertion extension.\n");
		fprintf(stderr,	"\t-d N\tN is a positive integer. -N will be used as the weight for the deletion opening.\n");
		fprintf(stderr, "\t-f N\tN is a positive integer. -N will be used as the weight for the deletion extension.\n");
		fprintf(stderr, "\t-a FILE\tFor protein sequence alignment. FILE is either the Blosum or Pam weight matrix. Recommend to use the matrix\n\t\tincluding B Z X * columns. Otherwise, corresponding scores will be signed to 0.\n"); 
		fprintf(stderr, "\t-c\tReturn the alignment in cigar format.\n");
		fprintf(stderr, "\t-r\tThe best alignment will be picked between the original read alignment and the reverse complement read alignment.\n\n");
		return 1;
	}

	// Parse score matrix.
	if (strcmp(mat_name, "\0"))	{
		FILE *f_mat = fopen(mat_name, "r");
		char line[128];
		mat = realloc(mat, 1024 * sizeof(int8_t));
		k = 0;
		while (fgets(line, 128, f_mat)) {
			if (line[0] == '*' || (line[0] >= 'A' && line[0] <= 'Z')) {
				char str[4], *s = str;
				str[0] = '\0';
				l = 1;
				while (line[l]) {
					if ((line[l] >= '0' && line[l] <= '9') || line[l] == '-') *s++ = line[l];	
					else if (str[0] != '\0') {					
						*s = '\0';
						mat[k++] = (int8_t)atoi(str);
						s = str;
						str[0] = '\0';			
					}
					++l;
				}
				if (str[0] != '\0') {
					*s = '\0';
					mat[k++] = (int8_t)atoi(str);
					s = str;
					str[0] = '\0';			
				}
			}
			m = k%24;
			while ((m%23 == 0 || m%22 == 0 || m%21 == 0 || m%20 == 0) && m != 0) {
				k++;
				m = k%24;
			} // If the weight matrix doesn't BZX*, set their values 0.
		}
		if (k == 0) {
			fprintf(stderr, "Problem of reading the weight matrix file.\n");
			return 1;
		} 
		fclose(f_mat);	
		if (mat[525] <= 0) {
			fprintf(stderr, "Improper weight matrix file format. Please use standard Blosum or Pam files.\n");
			return 1;
		}
		table = aa_table;
		n = 24;
	}

	ref_fp = gzopen(argv[optind], "r");
	ref_seq = kseq_init(ref_fp);

	start = clock();

	// alignment
	while ((l = kseq_read(ref_seq)) >= 0) {
		gzFile read_fp;
		kseq_t *read_seq;
		printf("ref_name: %s\n", ref_seq->name.s);
		read_fp = gzopen(argv[optind + 1], "r");
		read_seq = kseq_init(read_fp);
		int32_t refLen = strlen(ref_seq->seq.s);
		char* ref_reverse = seq_reverse(ref_seq->seq.s, refLen - 1); 
		while ((m = kseq_read(read_seq)) >= 0) {
			alignment_end* plus;
			int32_t word = 0;
			printf("read_name: %s\n", read_seq->name.s);
			int32_t readLen = strlen(read_seq->seq.s);
			plus = find_bests(ref_seq->seq.s, read_seq->seq.s, nt_table, nt_table, mat, n, refLen, readLen, insert_open, insert_extension, delet_open, delet_extension, &word);
			if (reverse == 0) align (ref_seq->seq.s, read_seq->seq.s, nt_table, nt_table, mat, ref_reverse, 5, refLen, match, mismatch, insert_open, insert_extension, delet_open, delet_extension, path, word, plus);
			else if (!strcmp(mat_name, "\0")){
				char* reverse;
				alignment_end* minus;
				int8_t complement_nt_table[128] = {
					4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
					4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
					4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
					4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
					4, 3, 4, 2,  4, 4, 4, 1,  4, 4, 4, 4,  4, 4, 4, 4, 
					4, 4, 4, 4,  0, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
					4, 3, 4, 2,  4, 4, 4, 1,  4, 4, 4, 4,  4, 4, 4, 4, 
					4, 4, 4, 4,  0, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4 
				};
				reverse = seq_reverse(read_seq->seq.s, readLen - 1);
				minus = find_bests(ref_seq->seq.s, reverse, nt_table, complement_nt_table, mat, n, refLen, readLen, insert_open, insert_extension, delet_open, delet_extension, &word);
				if (minus[0].score > plus[0].score) align (ref_seq->seq.s, reverse, nt_table, complement_nt_table, mat, ref_reverse, 5, refLen, match, mismatch, insert_open, insert_extension, delet_open, delet_extension, path, word, minus);
				else align (ref_seq->seq.s, read_seq->seq.s, nt_table, nt_table, mat, ref_reverse, 5, refLen, match, mismatch, insert_open, insert_extension, delet_open, delet_extension, path, word, plus);
			} else {
				fprintf(stderr, "The reverse complement alignment is not available for protein sequences.\n");
				return 1;
			}
		}
		free(ref_reverse);
		kseq_destroy(read_seq);
		gzclose(read_fp);
	}
	end = clock();
	cpu_time = ((float) (end - start)) / CLOCKS_PER_SEC;
	fprintf(stdout, "CPU time: %f seconds\n", cpu_time);

	free(mat);
	kseq_destroy(ref_seq);
	gzclose(ref_fp);
	return 0;
}
