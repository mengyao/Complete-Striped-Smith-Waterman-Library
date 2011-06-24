/*  main.c
 *  Created by Mengyao Zhao on 06/23/11.
 *	Version 0.1.4
 *	New features: pure c fasta parser, reach Farrar's speed, embeded run time testing function 
 */

#include <stdlib.h>
#include <emmintrin.h>
#include <zlib.h>
#include <stdio.h>
#include <time.h>
#include "ssw.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int main (int argc, char * const argv[]) {
	// validate argument count
	if (argc < 2) {
		fprintf (stderr, "USAGE: %s <ref.fa/fq> <reads.fa/fq>\n", argv[0]);
		exit (1);
	}

	clock_t start, end;
	float cpu_time;
	gzFile ref_fp;
	kseq_t *ref_seq;
	int l;
	int m;
	ref_fp = gzopen(argv[1], "r");
	ref_seq = kseq_init(ref_fp);

	start = clock();
	while ((l = kseq_read(ref_seq)) >= 0) {
		printf("ref_name: %s\n", ref_seq->name.s);
		/* printf("ref_seq: %s\n", ref_seq->seq.s); */
		gzFile read_fp = gzopen(argv[2], "r");
		kseq_t*	read_seq = kseq_init(read_fp);
		int32_t refLen = strlen(ref_seq->seq.s); 
		int32_t* ref_num = ref_nt2num (ref_seq->seq.s, refLen);
		while ((m = kseq_read(read_seq)) >= 0) {
			printf("read_name: %s\n", read_seq->name.s);
			printf("read_seq: %s\n", read_seq->seq.s); 
			
			int32_t end_seg = 0;
			int32_t readLen = strlen(read_seq->seq.s);
			__m128i* vProfile = queryProfile_constructor(read_seq->seq.s, 2, 1, 4);
			alignment_end* bests = smith_waterman_sse2(ref_num, refLen, readLen, 2, 1, 2, 1, vProfile, &end_seg, 4);
			free(vProfile);
			
			if (bests[0].score != 0) {
				fprintf(stdout, "max score: %d, end_ref: %d\nmax2 score: %d, end_ref: %d\n", 
						bests[0].score, bests[0].ref, bests[1].score, bests[1].ref);		
				}else {
				fprintf(stdout, "No alignment found for this read.\n");
			}
		}
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
