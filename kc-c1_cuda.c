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
	uint32_t bits, count;
	uint32_t *used;
	KeyValue *keys;
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
static uint32_t hash_uint64(uint64_t key) {
	key = ~key + (key << 21);
	key = key ^ key >> 24;
	key = (key + (key << 3)) + (key << 8);
	key = key ^ key >> 14;
	key = (key + (key << 2)) + (key << 4);
	key = key ^ key >> 28;
	key = key + (key << 31);
	return (uint32_t)key;
}

static uint32_t h2b(uint32_t hash, uint32_t bits) {
	return hash * 2654435769U >> (32 - bits);
}

HashTable *HashTable_init(){
	HashTable *ht;
	CALLOC(ht, 1);
	return ht;
}

void HashTable_destory(HashTable *ht) {
	if (!ht) return;
	free((void *)ht->keys); free(ht->used);
	free(ht);
}

unsigned int hash_insert(HashTable *ht, const uint64_t *key, int *absent) {

	unsigned int n_buckets, i, last, mask;
	n_buckets = ht->keys ? 1U << ht->bits : 0U;
	*absent = -1;

	// Incrementar el tamaño del hash table si no hay buckets disponibles
	// if (ht->count >= (n_buckets>>1) + (n_buckets>>2)) { /* rehashing */
	// 	if (hash_resize(ht, n_buckets + 1U) < 0)
	// 		return n_buckets;
	// 	n_buckets = 1U<<ht->bits;
	// }

	mask = n_buckets - 1;
	i = last = h2b(hash_uint64(*key), ht->bits);

	// Collition: Open addressing
	while (slot_used(ht->used, i) && !(ht->keys[i].key == *key)) {
		i = (i + 1U) & mask;
		if (i == last) break;
	}

	// Comprobar si se ha encontrado un slot vacío
	if (!slot_used(ht->used, i)) { // no se ha encontrado la llave
		ht->keys[i].key = *key;
		slot_set_used(ht->used, i);
		++ht->count;
		*absent = 1;
	} else *absent = 0; // retornar el indice si se he encontrado la llave
	return i;
}

static void count_seq_kmers(HashTable *ht, int k, int len, char *seq) // insert k-mers in $seq to hash table $ht
{
	int i, l;

	uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;

	for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
		int absent, c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 2 | c) & mask;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
			if (++l >= k) { // we find a k-mer

				unsigned int itr;
				uint64_t key = x[0] < x[1]? x[0] : x[1];

				//itr = kc_c1_put(h, key, &absent);
				itr = hash_insert(ht, key, &absent);// only add one strand!
				if (absent) {
					ht->keys[itr].val = 0;
				}

				ht->keys[itr].val++;

			}
		} else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
	}
}

// ignore
static void count_seq(kc_c1_t *h, int k, int len, char *seq) // insert k-mers in $seq to hash table $h
{
	int i, l;
	uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
	for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
		int absent, c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 2 | c) & mask;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
			if (++l >= k) { // we find a k-mer
				khint_t itr;
				uint64_t key = x[0] < x[1]? x[0] : x[1];

				itr = kc_c1_put(h, key, &absent); // only add one strand!
				if (absent) kh_val(h, itr) = 0;
				++kh_val(h, itr);
			}
		} else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
	}
}

static HashTable *count_file(const char *fn, int k)
{
	gzFile fp;
	kseq_t *ks;
	HashTable *ht;

	if ((fp = gzopen(fn, "r")) == 0) return 0;
	ks = kseq_init(fp); // descriptor fichero fastaq
	// ht = kc_c1_init(); // inicializar hashtable
	ht = HashTable_init(); // inicializar hashtable

	ReadSeqList *current, *head;
	head = current = NULL;

	// leer los datos del fichero de entrada y guardarlos en memoria
	while (kseq_read(ks) >= 0) {

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

	// contar
	for(current = head; current; current=current->next){
        //count_seq(h, k, current->length, current->sequence);
		count_seq_kmers(ht, k, current->length, current->sequence);
    }

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
	for (k = 0; k < ht_end(ht); ++k)
		if (slot_used(ht->used, k))
			++cnt[ht->keys[k].val < 256? ht->keys[k].val : 255];
	for (i = 1; i < 256; ++i)
		printf("%d\t%ld\n", i, (long)cnt[i]);
}

int main(int argc, char *argv[])
{
	// kc_c1_t *h;
	HashTable *ht;
	int c, k = 15;
	ketopt_t o = KETOPT_INIT;
	while ((c = ketopt(&o, argc, argv, 1, "k:", 0)) >= 0)
		if (c == 'k') k = atoi(o.arg);
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: kc-c1 [-k %d] <in.fa>\n", k);
		return 1;
	}
	ht = count_file(argv[o.ind], k);
	print_hist(ht);
	// kc_c1_destroy(h);
	HashTable_destory(ht);
	return 0;
}
