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

typedef struct __ReadSeqList {
	char* sequence;
	unsigned length;
	struct __ReadSeqList* next;
} ReadSeqList;

typedef struct HashTable {
	unsigned int bits;
	unsigned int count;
	unsigned int read_count;
	unsigned long long int *keys;
    unsigned int *values;
} HashTable;



__device__ const unsigned char seq_nt4_table[256] = { // translate ACGT to 0123
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
__device__ unsigned int hash_uint64(unsigned long long int key) {

	key = ~key + (key << 21);
	key = key ^ key >> 24;
	key = (key + (key << 3)) + (key << 8);
	key = key ^ key >> 14;
	key = (key + (key << 2)) + (key << 4);
	key = key ^ key >> 28;
	key = key + (key << 31);
	return (unsigned int)key;
}

HashTable* HashTable_init(unsigned int bits, unsigned int read_count){
    HashTable *ht;
    ht = (HashTable*)calloc(1, sizeof(HashTable));

	ht->read_count = read_count;
	ht->bits = bits;
	ht->count = 0;

    return ht;
}


void HashTable_destory(HashTable *ht) {
	if (!ht) return;
	free(ht);
}


__device__ unsigned int h2b(unsigned int hash, unsigned long long int bits) {
    return hash * 2654435769U >> (32 - bits);
}

__device__ void hash_insert(HashTable *ht, unsigned long long int kmer) {

	unsigned int iKey, last;
	//bool end = false;


	//iKey = last = hash_uint64(kmer) * (2654435769U >> (32 - ht->bits));
    iKey = last = h2b(hash_uint64(kmer), ht->bits);
    while (ht->values[iKey] > 0 && ht->keys[iKey] != kmer) {
        iKey = (iKey + 1U) & ((1U << ht->bits) - 1);
        if (iKey == last) break;
    }

    // Comprobar si se ha encontrado un slot vacío
    if (ht->values[iKey] == 0) { // no se ha encontrado la llave

        ht->keys[iKey] = kmer;
        ht->values[iKey] = 1;
        ++ht->count;

    } else {
        ht->values[iKey]++;
    } 

    /*
	while (true)
	{
		unsigned int prev = atomicCAS(ht->keys[iKey], NULL, iKey);

		if (prev == NULL || prev == kmer) {
			ht->keys[iKey] = kmer;
			atomicAdd(&(ht->values[iKey]), 1);

			if(prev == NULL) atomicAdd(&(ht->count), 1);

			return;
		}

		if(end) return;

		// Collition: Open addressing
		iKey = (iKey + 1U) & (ht->bits - 1);

		// loop back
		end = (iKey == last);

	}
    */

}

// insert k-mers in $seq to hash table $ht
__global__ void kernel_count_seq_kmers(HashTable *ht, int k, char **d_reads)
{
	unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if(tid < ht->read_count) {

        int i, l;
		char *seq = d_reads[tid];
		//int len = strlen(seq);
		int len = 100;
        unsigned long long int x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;

		for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
			int c = seq_nt4_table[(unsigned char)seq[i]];
			if (c < 4) { // not an "N" base
				x[0] = (x[0] << 2 | c) & mask;                  // forward strand
				x[1] = x[1] >> 2 | (unsigned long long int)(3 - c) << shift;  // reverse strand
				if (++l >= k) { // we find a k-mer

					unsigned long long int kmer = x[0] < x[1]? x[0] : x[1];
					hash_insert(ht, kmer); // only add one strand!
				}
			} else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
		}
	}
}

__global__ void kernel_print_hist(const HashTable *ht, unsigned int *cnt_d)
{
	unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned pos;


	if(threadIdx.x < ht->bits) {
		if (ht->values[tid] != 0) {
			pos = ht->values[tid] < 256 ? ht->values[tid] : 255;
            cnt_d[pos]++;
			//atomicAdd(&cnt_d[pos], (unsigned int)1);
		}
	}
}

static int count_file(const char *fn, int k, unsigned int p)
{
	//gzFile fp;
	//kseq_t *ks;
	HashTable *ht;
    unsigned int capacity = 1U << p;
    unsigned int cnt[256];
    unsigned int read_count = 0;




    // variables para cuda
	HashTable *ht_d;
	char **reads_d;
	unsigned int *cnt_d;

    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;

    fp = fopen(fn, "r");
    if (fp == NULL) exit(EXIT_FAILURE);

   
	ReadSeqList *current, *head;
	head = current = NULL;

    while ((read = getline(&line, &len, fp)) != -1) {

        read_count++;

		ReadSeqList *node = (ReadSeqList*)malloc(sizeof(ReadSeqList));
        node->sequence = (char*)malloc(strlen(line));
        strcpy(node->sequence, line);
        node->length = read;
        node->next =NULL;

        if(head == NULL){
            current = head = node;
        } else {
            current = current->next = node;
        }


    }

    fclose(fp);
    if (line) free(line);

    printf("%d\n", read_count);
    unsigned int i;

    char **reads = (char**)malloc(read_count * sizeof(char*)); 

	for(i=0, current = head; current; current=current->next){
        reads[i] = (char*)malloc(current->length);
        sprintf(reads[i], "%s", current->sequence);
        i++;
    }
   

    // inicializar hashtable
	ht = HashTable_init(p, read_count);

    
    unsigned long long int *keys_d;
    unsigned int *values_d;


	// allocate memory in device
	cudaMalloc((void **)&ht_d, sizeof(HashTable));
   	cudaMalloc((void **)&keys_d, capacity * sizeof(unsigned long long int));
	cudaMalloc((void **)&values_d, capacity * sizeof(unsigned int));
	cudaMalloc((void **)&cnt_d, 256 * sizeof(unsigned int));
   	cudaMemset(keys_d, 0, capacity * sizeof(unsigned long long int));
	cudaMemset(values_d, 0, capacity * sizeof(unsigned int));
	cudaMemset(cnt_d, 0, 256 * sizeof(unsigned int));

    
	// copy data to device

    ht->keys = keys_d;
    ht->values = values_d;

	cudaMemcpy(ht_d, ht, sizeof(HashTable), cudaMemcpyHostToDevice);

	char **tmp = (char**)malloc (read_count * sizeof (char*));
    for (int i = 0; i < read_count; i++) {
        cudaMalloc ((void **)&tmp[i], head->length * sizeof (char));
    }

	cudaMalloc((void **)&reads_d, read_count * sizeof(char*));

    cudaMemcpy(reads_d, tmp, read_count * sizeof (char*), cudaMemcpyHostToDevice);
    for (int i = 0; i < read_count; i++) {
        cudaMemcpy(tmp[i], reads[i], head->length * sizeof (char), cudaMemcpyHostToDevice);
    }
    free(tmp);



    printf("total reads: %d\n", read_count);


    // invocar kernels
    int thr = 1024;

	//kernel_count_seq_kmers<<<1,1>>>(ht_d, k, reads_d);
	kernel_count_seq_kmers<<<ceil(read_count/thr), thr>>>(ht_d, k, reads_d);

    cudaDeviceSynchronize();


	//kernel_print_hist<<<ceil(capacity/thr), thr>>>(ht_d, cnt_d);

    //cudaDeviceSynchronize();

	cudaMemcpy(ht, ht_d, sizeof(HashTable), cudaMemcpyDeviceToHost);
	cudaMemcpy(ht->keys, keys_d, capacity * sizeof(unsigned long long int), cudaMemcpyDeviceToHost);
	cudaMemcpy(ht->values, values_d, capacity * sizeof(unsigned int), cudaMemcpyDeviceToHost);
	cudaMemcpy(cnt, cnt_d, capacity * sizeof(unsigned int), cudaMemcpyDeviceToHost);


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
    cudaFree(reads_d);
	cudaFree(ht_d);
	cudaFree(cnt_d);
	cudaFree(keys_d);
	cudaFree(values_d);
    return 0;
	// limpieza
    i = 0;
	for(current = head; current; current=current->next){
        free(current->sequence);
        free(current);
        free(reads[i]);
        i++;
    }

    free(reads);
	HashTable_destory(ht);
    return 0;
}


int main(int argc, char *argv[])
{
	int k = 31;
    unsigned int p = 27;

/*
	ketopt_t o = KETOPT_INIT;
	while ((c = ketopt(&o, argc, argv, 1, "k:", 0)) >= 0)
		if (c == 'k') k = atoi(o.arg);
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: kc-c1 [-k %d] <in.fa>\n", k);
		return 1;
	}
*/
    k = (int)strtol(argv[1], NULL, 10);
	count_file(argv[2], k, p);

	return 0;
}
