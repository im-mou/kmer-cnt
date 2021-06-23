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

#include "kseq.h" // FASTA/Q parser
KSEQ_INIT(gzFile, gzread)

#include "khashl.h" // hash table
KHASHL_MAP_INIT(, kc_c1_t, kc_c1, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)

#define CALLOC(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))
#define slot_used(flag, i) (flag[i>>5] >> (i&0x1fU) & 1U)
#define slot_set_used(flag, i)   (flag[i>>5] |= 1U<<(i&0x1fU))
#define ht_end(h) ((h)->keys? 1U<<(h)->bits : 0U)

typedef struct __ReadSeqList {
	char* sequence;
	unsigned length;
	struct __ReadSeqList* next;
} ReadSeqList;

typedef struct __KeyValue {
    uint64_t key;
    uint32_t val;
} __attribute__ ((__packed__)) KeyValue;

typedef struct HashTable {
	uint32_t bits;
	uint32_t count;
	uint32_t *collition;
	uint32_t *used;
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
static inline uint32_t hash_uint64(uint64_t key) {

	key = ~key + (key << 21);
	key = key ^ key >> 24;
	key = (key + (key << 3)) + (key << 8);
	key = key ^ key >> 14;
	key = (key + (key << 2)) + (key << 4);
	key = key ^ key >> 28;
	key = key + (key << 31);
	return (uint32_t)key;
}

static inline uint32_t h2b(uint32_t hash, uint32_t bits) {
	return hash * 2654435769U >> (32 - bits);
}

HashTable *HashTable_init(uint32_t bits){
	HashTable *ht;
	CALLOC(ht, 1);
	CALLOC(ht->keys, 1U<<bits);
	CALLOC(ht->values, 1U<<bits);

    //ht->used =(uint32_t*) malloc(((1U<<bits) < 32 ? 1 : (1U<<bits)>>5) * sizeof(uint32_t));
    //memset(ht->used,0,  ((1U<<bits) < 32 ? 1 : (1U<<bits)>>5) * sizeof(uint32_t));
	//CALLOC(ht->used, 1U<<bits);
	ht->bits = bits;
	ht->count = 0;
	//ht->collition = 0;
    //CALLOC(ht->collition, 1U<<bits);

	return ht;
}

void HashTable_destory(HashTable *ht) {
	if (!ht) return;
	free((void *)ht->keys);
	free((void *)ht->values);
	//free((void *)ht->collition);
	free(ht->used);
	free(ht);
}


void hash_insert(HashTable *ht, uint64_t kmer) {

	unsigned int i, last;
	//unsigned int n_buckets, i, last, mask;
	//n_buckets = 1U << ht->bits;

	//mask = n_buckets - 1;I
	i = last = h2b(hash_uint64(kmer), ht->bits);
	
	// Collition: Open addressing
	//while (slot_used(ht->used, i) && ht->keys[i] != kmer) {
	while (ht->values[i] > 0 && ht->keys[i] != kmer) {
        //ht->collition[i]++;
		i = (i + 1U) & ((1U << ht->bits) - 1);
		if (i == last) break;
	}

	// Comprobar si se ha encontrado un slot vacío
	//if (!slot_used(ht->used, i)) { // no se ha encontrado la llave
	if (ht->values[i] == 0) { // no se ha encontrado la llave

		ht->keys[i] = kmer;

        ht->values[i] = 1;
		//slot_set_used(ht->used, i);

		++ht->count;
	} else {
        // retornar el indice si se he encontrado la llave
        ht->values[i]++;
    }  
	//return i;
}

static void count_seq_kmers(HashTable *ht, int k, int len, char *seq) // insert k-mers in $seq to hash table $ht
{
	int i, l;

	uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;

	for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 2 | c) & mask;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
			if (++l >= k) { // we find a k-mer

				//unsigned int itr;
				uint64_t kmer = x[0] < x[1]? x[0] : x[1];

                //printf("%"PIRIu64"\n", kmer);
				//itr = hash_insert(ht, key, &absent);// only add one strand!
				hash_insert(ht, kmer);// only add one strand!
				//if (absent) {
					//ht->values[itr] = 0;
				//}

				//ht->values[itr]++;
			}
		} else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
	}
}

static HashTable *count_file(const char *fn, int k, uint32_t p)
{
	gzFile fp;
	kseq_t *ks;
	HashTable *ht;

	if ((fp = gzopen(fn, "r")) == 0) return 0;
	ks = kseq_init(fp); // descriptor fichero fastaq
	// ht = kc_c1_init(); // inicializar hashtable
	ht = HashTable_init(p); // inicializar hashtable

	ReadSeqList *current, *head;
	head = current = NULL;

	// leer los datos del fichero de entrada y guardarlos en memoria
    uint32_t read_count = 0; 
	while (kseq_read(ks) >= 0) {
        read_count++;

		ReadSeqList *node = malloc(sizeof(ReadSeqList));
        node->sequence = malloc(strlen(ks->seq.s) + 1);
        strcpy(node->sequence, ks->seq.s);
        node->length = ks->seq.l;
        node->next =NULL;

        if(head == NULL){
            current = head = node;
        } else {
            current = current->next = node;
        }

	}
    
    unsigned int i;
    // crear un array de tamaño fijo para almacenar las lecturas
    char** reads = malloc(read_count * sizeof(char*));
    

	for(i=0, current = head; current; current=current->next){
		//count_seq_kmers(ht, k, current->length, current->sequence);
		reads[i] = current->sequence;
        i++;
    }

    printf("%s\n", reads[0]);

    printf("total reads: %d\n", read_count);

    current = head;
    for(i = 0; i<read_count ; i++){
		count_seq_kmers(ht, k, strlen(reads[i]), reads[i]);
    } 

    printf("COUNT: %d\n\n", ht->count);
/*
    for(i = 0; i< 5000 ; i++){
        printf("key: ");
        printf("%"PRIu64"\t", ht->keys[i]);
        printf("value: %d\t", ht->values[i]);
        printf("collitions: %d\n", ht->collition[i]);
    }
*/
		// limpieza
	for(current = head; current; current=current->next){
        free(current->sequence);
        free(current);
    }

	kseq_destroy(ks);
	gzclose(fp);
	return ht;
}

static void print_hist(const HashTable *ht)
{
	uint32_t k;
	uint64_t cnt[256];
	int i;
	for (i = 0; i < 256; ++i) cnt[i] = 0;
		//if (slot_used(ht->used, k))
	for (k = 0; k < (1U<<(ht)->bits); ++k)
		if (ht->values[k] != 0)
			++cnt[ht->values[k] < 256 ? ht->values[k] : 255];
	for (i = 1; i < 256; ++i)
		printf("%d\t%ld\n", i, (long)cnt[i]);
}

int main(int argc, char *argv[])
{
	// kc_c1_t *h;
	HashTable *ht;
	int c, k = 15;
    uint32_t p = 27;
	ketopt_t o = KETOPT_INIT;
	while ((c = ketopt(&o, argc, argv, 1, "k:", 0)) >= 0)
		if (c == 'k') k = atoi(o.arg);
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: kc-c1 [-k %d] <in.fa>\n", k);
		return 1;
	}
	ht = count_file(argv[o.ind], k, p);
	print_hist(ht);
	// kc_c1_destroy(h);
	HashTable_destory(ht);
	return 0;
}
