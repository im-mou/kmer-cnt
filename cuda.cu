#include <cuda.h>
#include <cuda_fp16.h>

#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <fcntl.h>
#include <limits.h>

//#define __STDC_FORMAT_MACROS
//#include <inttypes.h>

#include "kseq.h" // FASTA/Q parser
#include "ketopt.h" // command-line argument parser
//extern "C"
//{
	//#include "kseq.h" // FASTA/Q parser
	//#include "ketopt.h" // command-line argument parser
//}

KSEQ_INIT(gzFile, gzread)

typedef unsigned char uint8_c;
typedef unsigned int uint32_c;
typedef unsigned long long int uint64_c;

typedef struct __ReadSeqList {
	char* sequence;
	unsigned length;
	struct __ReadSeqList* next;
} ReadSeqList;

typedef struct HashTable {
	uint32_t bits;
	uint32_t count;
	//uint32_t *collition;
	uint64_t *keys;
    uint32_t *values;
} HashTable;



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


// funcion para calcular un hash de 64 bits
__device__ uint32_t hash_uint64(uint64_t key) {

	key = ~key + (key << 21);
	key = key ^ key >> 24;
	key = (key + (key << 3)) + (key << 8);
	key = key ^ key >> 14;
	key = (key + (key << 2)) + (key << 4);
	key = key ^ key >> 28;
	key = key + (key << 31);
	return (uint32_t)key;
}

HashTable* HashTable_init(uint32_t bits){
    HashTable *ht;
    ht = (HashTable*)calloc(1, sizeof(HashTable));

	uint32_t capacity = 1U << bits;
	ht->bits = capacity;
	ht->count = 0;

    return ht;
}


void HashTable_destory(HashTable *ht) {
	if (!ht) return;
	free((void *)ht->keys);
	free((void *)ht->values);
	free(ht);
}


__device__ void hash_insert(HashTable *ht, uint64_t kmer, uint32_t capacity_d) {

	unsigned int iKey, last;
	bool end = false;

	iKey = last = hash_uint64(kmer) * (2654435769U >> (32 - ht->bits));

	while (true)
	{
		uint32_t prev = atomicCAS(ht->keys[iKey], NULL, iKey);

		if (prev == NULL || prev == kmer) {
			ht->keys[iKey] = kmer;
			atomicAdd(&(ht->values[iKey]), 1);

			if(prev == NULL) atomicAdd(&(ht->count), 1);

			return;
		}

		if(end) return;

		// Collition: Open addressing
		iKey = (iKey + 1U) & (capacity_d - 1);

		// loop back
		end = (iKey == last);

	}

}

// insert k-mers in $seq to hash table $ht
__global__ void kernel_count_seq_kmers(HashTable *ht, int k, char **d_reads, uint32_t read_count)
{
	if(threadIdx.x < read_count) {

        int i, l;
		char *seq = d_reads[threadIdx.x];
		int len = strlen(seq);
        uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;

		for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
			int c = seq_nt4_table[(uint8_t)seq[i]];
			if (c < 4) { // not an "N" base
				x[0] = (x[0] << 2 | c) & mask;                  // forward strand
				x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
				if (++l >= k) { // we find a k-mer

					uint64_t kmer = x[0] < x[1]? x[0] : x[1];
					hash_insert(ht, kmer); // only add one strand!

				}
			} else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
		}
	}
}

__global__ void kernel_print_hist(const HashTable *ht, uint32_t *cnt_d)
{
	unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned pos;


	if(threadIdx.x < 256) {
		if (ht->values[tid] != 0) {
			pos = ht->values[tid] < 256 ? ht->values[tid] : 255;
			atomicAdd(&cnt_d[pos], (uint32_t)1);
		}
	}
}

static void count_file(const char *fn, int k, uint32_t p)
{
	gzFile fp;
	kseq_t *ks;
	HashTable *ht;
	uint32_t capacity = 1U << p;
    uint32_t cnt[256];	

    // inicializar hashtable
	ht = HashTable_init(p);

    // variables para cuda
	HashTable *ht_d;
	char **reads_d;
	uint32_t *cnt_d;
	uint32_t *capacity_d;
	uint32_t *read_count_d;



    // leer los reads
	if ((fp = gzopen(fn, "r")) == 0) return 0;
	ks = kseq_init(fp); // descriptor fichero fastaq


	ReadSeqList *current, *head;
	head = current = NULL;

	// leer los datos del fichero de entrada y guardarlos en memoria
    uint32_t read_count = 0;
	while (kseq_read(ks) >= 0) {
        read_count++;

		ReadSeqList *node = (ReadSeqList*)malloc(sizeof(ReadSeqList));
        node->sequence = (char*)malloc(strlen(ks->seq.s) + 1);
        strcpy(node->sequence, ks->seq.s);
        node->length = ks->seq.l;
        node->next =NULL;

        if(head == NULL){
            current = head = node;
        } else {
            current = current->next = node;
        }

	}

	kseq_destroy(ks);
	gzclose(fp);

    // crear un array de tamaño fijo para almacenar las lecturas
    // char **reads = malloc(read_count * sizeof(char*));



	


	// allocate memory in device
	cudaMalloc((void **)&reads_d, read_count * sizeof(char *));
	cudaMalloc((void **)&ht_d, sizeof(HashTable));
	cudaMalloc((void **)ht_d->keys, capacity * sizeof(uint64_t));
	cudaMalloc((void **)ht_d->values, capacity * sizeof(uint32_t));
	//cudaMalloc((void **)ht_d->collition, capacity * sizeof(uint32_t));
	cudaMalloc((void **)&cnt_d, 256 * sizeof(uint32_t));
	cudaMalloc((void **)&capacity_d, sizeof(uint32_t));
	cudaMalloc((void **)&read_count_d, sizeof(uint32_t));

	cudaMemset(ht_d->keys, 0, capacity * sizeof(uint64_t));
	cudaMemset(ht_d->values, 0, capacity * sizeof(uint32_t));
	cudaMemset(cnt_d, 0, 256 * sizeof(uint32_t));

	// copy data to device
	cudaMemcpy(ht_d, ht, sizeof(HashTable), cudaMemcpyHostToDevice);
	cudaMemcpy(capacity_d, capacity, sizeof(uint32_t), cudaMemcpyHostToDevice);
	cudaMemcpy(read_count_d, read_count, sizeof(uint32_t), cudaMemcpyHostToDevice);





	// copiar los read a la memoria de la GPU
	char **temp_reads_d = (char **)malloc(read_count * sizeof(char *));
    unsigned int i;

	for(i=0, current = head; current; current=current->next){
		cudaMalloc((void **)&(temp_reads_d[i]), strlen(current->sequence) * sizeof(char));
		cudaMemcpy(temp_reads_d[i], current->sequence, strlen(current->sequence) * sizeof(char), cudaMemcpyHostToDevice);
		cudaMemcpy(reads_d + i, &(temp_reads_d[i]), sizeof(char *), cudaMemcpyHostToDevice);
		i++;
    }



    //printf("%s\n", reads[0]);

    printf("total reads: %d\n", read_count);


    // invocar kernels

	kernel_count_seq_kmers<<<ceil(read_count/1024), 1024>>>(ht_d, k, reads_d, read_count_d);

	kernel_print_hist<<<ceil(ht_d->count/256), 256>>>(ht_d, cnt_d, capacity_d);



	cudaMemcpy(ht, ht_d, sizeof(HashTable), cudaMemcpyDeviceToHost);
	cudaMemcpy(ht->keys, ht_d->keys, capacity * sizeof(uint64_t), cudaMemcpyDeviceToHost);
	cudaMemcpy(ht->values, ht_d->values, capacity * sizeof(uint32_t), cudaMemcpyDeviceToHost);
	cudaMemcpy(cnt, cnt_d, capacity * sizeof(uint32_t), cudaMemcpyDeviceToHost);


    printf("COUNT: %d\n\n", ht->count);

    // for(i = 0; i< 5000 ; i++){
    //     printf("key: ");
    //     printf("%"PRIu64"\t", ht->keys[i]);
    //     printf("value: %d\t", ht->values[i]);
    //     printf("collitions: %d\n", ht->collition[i]);
    // }

	for (i = 1; i < 256; ++i)
		printf("%d\t%ld\n", i, (long)cnt[i]);


	// limpieza
	for(current = head; current; current=current->next){
        free(current->sequence);
        free(current);
    }

	cudaFree(reads_d);
	cudaFree(ht_d);
	cudaFree(cnt_d);

	HashTable_destory(ht);
}


int main(int argc, char *argv[])
{
	int c, k = 31;
    uint32_t p = 27;
    
	ketopt_t o = KETOPT_INIT;
	while ((c = ketopt(&o, argc, argv, 1, "k:", 0)) >= 0)
		if (c == 'k') k = atoi(o.arg);
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: kc-c1 [-k %d] <in.fa>\n", k);
		return 1;
	}

	count_file(argv[o.ind], k, p);

	return 0;
}
