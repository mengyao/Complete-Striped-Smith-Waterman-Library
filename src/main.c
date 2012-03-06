/*  main.c
 *  Created by Mengyao Zhao on 06/23/11.
 *	Version 0.1.4
 *  Last revision by Mengyao Zhao on 03/06/12.
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

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

KSEQ_INIT(gzFile, gzread)

int8_t* char2num (char* seq, int8_t* table, int32_t l) {	// input l: 0; output l: length of the sequence
	int32_t i;
	int8_t* num = (int8_t*)calloc(l, sizeof(int8_t));
	for (i = 0; i < l; ++i) num[i] = table[(int)seq[i]];
	return num;
}

char* reverse_comple(const char* seq) {
	int32_t end = strlen(seq), start = 0;
	char* rc = (char*)calloc(end + 1, sizeof(char));
	int8_t rc_table[128] = {
		4, 4,  4, 4,  4,  4,  4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4,  4, 4,  4,  4,  4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4,  4, 4,  4,  4,  4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4,  4, 4,  4,  4,  4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 54, 4, 47, 4,  4,  4, 43, 4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4,  4, 4,  41, 41, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 54, 4, 47, 4,  4,  4, 43, 4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4,  4, 4,  41, 41, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4 
	};
	rc[end] = '\0';
	-- end;				
	while (LIKELY(start < end)) {			
		rc[start] = (char)rc_table[(int8_t)seq[end]];		
		rc[end] = (char)rc_table[(int8_t)seq[start]];		
		++ start;					
		-- end;						
	}					
	if (start == end) rc[start] = (char)rc_table[(int8_t)seq[start]];			
	return rc;					
}							

int main (int argc, char * const argv[]) {
	clock_t start, end;
	float cpu_time;
	gzFile read_fp;
	kseq_t *read_seq;
	int32_t l, m, k, match = 2, mismatch = 2, insert_open = 3, insert_extension = 1, delet_open = 3, delet_extension = 1, path = 0, reverse = 0, n = 5;
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
		n = 24;
		table = aa_table;
	}

	read_fp = gzopen(argv[optind + 1], "r");
	read_seq = kseq_init(read_fp);
	start = clock();

	// alignment
	while ((m = kseq_read(read_seq)) >= 0) {
		gzFile ref_fp;
		kseq_t *ref_seq;
		init_param* init = (init_param*)calloc(1, sizeof(init_param));
		profile* p, *p_rc = 0;
		int32_t readLen = read_seq->seq.l; 
		char* read_rc = 0;
		
		printf("read_name: %s\n", read_seq->name.s);
		init->read = char2num(read_seq->seq.s, table, readLen);
		init->readLen = readLen;
		init->mat = mat;
		init->score_size = 2;
		init->n = n;
		p = ssw_init(init);
		if (reverse == 1 && n == 5) {
			read_rc = reverse_comple(read_seq->seq.s);
			init->read = char2num(read_rc, table, readLen);
			p_rc = ssw_init(init);
		}else if (reverse == 1 && n == 24) {
			fprintf (stderr, "Reverse complement alignment is not available for protein sequences. \n");
			return 1;
		}

		ref_fp = gzopen(argv[optind], "r");
		ref_seq = kseq_init(ref_fp);
		while ((l = kseq_read(ref_seq)) >= 0) {
			align_param* a = (align_param*)calloc(1, sizeof(align_param));
			align* result, *result_rc = 0;
			int32_t refLen = ref_seq->seq.l;
			int8_t strand = 0;

			a->prof = p;
			a->ref = char2num(ref_seq->seq.s, table, refLen);
			a->refLen = refLen;
			a->weight_insertB = insert_open;
			a->weight_insertE = insert_extension;
			a->weight_deletB = delet_open;
			a->weight_deletE = delet_extension;
			if (path == 1) {
				a->begin = 1;
				a->align = 1;
			} else {
				a->begin = 1;
				a->align = 0;
			}
			printf("ref_name: %s\n", ref_seq->name.s);
			result = ssw_align (a);
			if (reverse == 1) {
				a->prof = p_rc;
				result_rc = ssw_align(a);
			}
			
			if (result_rc && result_rc->score1 > result->score1) {
				fprintf(stdout, "%d\t%s\n", strand, read_rc);
				fprintf(stdout, "score1: %d\tscore2: %d\tref_begin1: %d\tref_end1: %d\tread_begin1: %d\tread_end1: %d\tref_end2: %d\n", result_rc->score1, result_rc->score2, result_rc->ref_begin1, result_rc->ref_end1, result_rc->read_begin1, result_rc->read_end1, result_rc->ref_end2);
				if (path == 1) fprintf(stdout, "cigar: %s\n\n", result_rc->cigar);
			} else {
				fprintf(stdout, "%d\t%s\n", strand, read_seq->seq.s);
				fprintf(stdout, "score1: %d\tscore2: %d\tref_begin1: %d\tref_end1: %d\tread_begin1: %d\tread_end1: %d\tref_end2: %d\n", result->score1, result->score2, result->ref_begin1, result->ref_end1, result->read_begin1, result->read_end1, result->ref_end2);
				if (path == 1) fprintf(stdout, "cigar: %s\n\n", result->cigar);
			}
			if (result_rc) align_destroy(result_rc);
			align_destroy(result);
			free(a->ref);
			free(a);
		}
		
		if(p_rc) init_destroy(p_rc);
		init_destroy(p);
		free(init);		
		kseq_destroy(ref_seq);
		gzclose(ref_fp);
	}
	end = clock();
	cpu_time = ((float) (end - start)) / CLOCKS_PER_SEC;
	fprintf(stdout, "CPU time: %f seconds\n", cpu_time);

	free(mat);
	kseq_destroy(read_seq);
	gzclose(read_fp);
	return 0;
}
