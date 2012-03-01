/*  main.c
 *  Created by Mengyao Zhao on 06/23/11.
 *	Version 0.1.4
 *  Last revision by Mengyao Zhao on 03/01/12.
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

int main (int argc, char * const argv[]) {
	clock_t start, end;
	float cpu_time;
	gzFile read_fp;
	kseq_t *read_seq;
	int32_t l, m, k, n = 5, match = 2, mismatch = 2, insert_open = 3, insert_extension = 1, delet_open = 3, delet_extension = 1, path = 0, reverse = 0;
	int8_t* mat = (int8_t*)calloc(25, sizeof(int8_t));
	char mat_name[16];
	mat_name[0] = '\0';

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


	read_fp = gzopen(argv[optind + 1], "r");
	read_seq = kseq_init(read_fp);
	start = clock();

	// alignment
	while ((m = kseq_read(read_seq)) >= 0) {
		gzFile ref_fp;
		kseq_t *ref_seq;
		init_param* init;
		profile* p;

		printf("read_name: %s\n", read_seq->name.s);
		printf("read_seq: %s\n\n", read_seq->seq.s);
		init->read = read_seq->seq.s;
		init->mat = mat;
		init->score_size = 2;
		init->reverse = 1;
		p = ssw_init(init);
		ref_fp = gzopen(argv[optind], "r");
		ref_seq = kseq_init(ref_fp);

		while ((l = kseq_read(ref_seq)) >= 0) {
			align_param* a;
			align* result;

			align_param->prof = p;
			align_param->ref = ref_seq->seq.s;
			align_param->refLen = strlen(ref_seq->seq.s);
			align_param->weight_insertB = insert_open;
			align_param->weight_insertE = insert_extension;
			align_param->weight_deletB = delet_open;
			align_param->weight_deletE = delet_extension;
			align_param->begin = 1;
			align_param->align = 1;
			printf("ref_name: %s\n", ref_seq->name.s);
			result = ssw_align (a);
			fprintf(stdout, "%d\t%s\n", result->strand, result->read);
			fprintf(stdout, "score1: %d\tscore2: %d\tref_begin1: %d\tref_end1: %d\tread_begin1: %d\tread_end1: %d\tref_end2: %d\n", result->score1, result->score2, result->ref_begin1, result->ref_end1, result->read_begin1, result->read_end1, result->ref_end2);
			fprintf(stdout, "cigar: %s\n\n", result->cigar);
			align_destroy(result);
		}
		
		init_destroy(p);
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
