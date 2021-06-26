#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include "ketopt.h" // command-line argument parser

// #include <cuda.h>
// #include <cuda_fp16.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <fcntl.h>
#include <limits.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include <string.h>

#include "kseq.h" // FASTA/Q parser
KSEQ_INIT(gzFile, gzread)


const unsigned char seq_nt4_table[256] = { // translate ACGT to 0123
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};


static void parse_file(const char *fn)
{
	gzFile fp;
	kseq_t *ks;

	if ((fp = gzopen(fn, "r")) == 0) exit(EXIT_FAILURE);

	ks = kseq_init(fp); // descriptor fichero fastaq


    FILE * fPtr;
    fPtr = fopen("parsed_data2.txt", "w");
    if(fPtr == NULL)
    {
        printf("Unable to create file.\n");
        exit(EXIT_FAILURE);
    }
    
	while (kseq_read(ks) >= 0) {
        fprintf(fPtr, "%s\n", ks->seq.s);
	}

    fclose(fPtr);

	kseq_destroy(ks);
	gzclose(fp);
}


int main(int argc, char *argv[])
{

    uint32_t k = 15;
	ketopt_t o = KETOPT_INIT;
    int c;

	while ((c = ketopt(&o, argc, argv, 1, "k:", 0)) >= 0){
		if (c == 'k') k = atoi(o.arg);
    }
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: parse-data <input.fa.gz>\n");
		return 1;
	}

    parse_file(argv[o.ind]);

	return 0;
}
